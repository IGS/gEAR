import ipaddress
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from urllib.parse import urlparse

import geardb
import gosling as gos
import requests
from flask import request
from flask_restful import Resource
from gear.trackhub import (
    fetch_trackdb_path,
    parse_hub_from_file,
    parse_tracks_from_trackdb,
)

"""
NOTE: The "gos" documentation kind of sucks, and requires some trial and error to figure out where things go.
Generally you can make a data class type (i.e. BedData) and add that as a property to a gos.Track
There are some shorthand arrangement functions (i.e. gos.overlay) to combine tracks into a gos.View
The spec is validated against the JSON schema and will error if validation fails.  Building in JS only gives warnings.

The "gos" classes and functions also have a lot of UndefinedTypes, which the linter does not like.
Hence all the "type: ignore" statements in this API file.

Documentation at https://gosling-lang.github.io/gos/index.html (the search function is helpful to find examples)
"""

# Goals
# - Convert UCSC hub tracks to Gosling specs
# - Create a Gosling spec using files alone (no UCSC hub)
# - Provide a RESTful API for Gosling specs
# - Establish a structure for Gosling specs

# CONSTANTS
TWO_LEVELS_UP = 2  # Number of parent directories to go up
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP]  # web-root dir
VIEW_PADDING = 10  # Padding between views
CONDENSED_WIDTH = 1200 - VIEW_PADDING * 2  # Width of the condensed view tracks
EXPANDED_WIDTH = 600 - VIEW_PADDING / 2  # Width of the left view tracks

CONDENSED_HEIGHT = 25  # Height for condensed tracks (25 is lowest height to still show axes)
EXPANDED_HEIGHT = 50  # Height for expanded tracks

# These files will be based on Ensembl's annotation naming structure, but sorted in chromosome order.
# I am adding a 2nd column of 1's to allow us to use the files in a "genomic" track.
ASSEMBLY_TO_CHROMSIZES_FILE = {
    "danRer10": "danRer10.chromInfo.txt",  # zebrafish
    "galGal6": "galGal6.chromInfo.txt",  # chicken
    "hg19": "hg19.chromInfo.txt",
    "hg38": "hg38.chromInfo.txt",
    "mm10": "mm10.chromInfo.txt",
    # "mm39": "mm39.chromInfo.txt",
    "rn6": "rn6.chromInfo.txt",  # rat
    # "calJac3": "calJac3.chromInfo.txt", # marmoset
}

GENOMES_ROOT = "https://umgear.org/tracks/genomes/"
#GENOMES_ROOT = "http://localhost:8080/tracks/genomes"


def _resolve_track_url(track: dict, use_gosling: bool = True) -> str | None:
    """
    Resolve the appropriate URL for a track based on context and available fields.

    For Gosling visualization, prefers gos_url if available (supports different file types).
    For UCSC export, uses bigDataUrl directly.

    Args:
        track (dict): Track specification dictionary.
        use_gosling (bool): If True, prefer gos_url for Gosling compatibility. If False, use bigDataUrl.

    Returns:
        str or None: The resolved URL, or None if URL cannot be determined.
    """
    # Gosling context: prefer gos_url if available
    if use_gosling and track.get("gos_url"):
        return track["gos_url"]

    # Fall back to bigDataUrl
    return track.get("bigDataUrl")

def _validate_hub_url(hub_url: str) -> None:
    """
    Validate hub URL to prevent SSRF attacks.

    Args:
        hub_url (str): The hub URL to validate.

    Raises:
        ValueError: If URL is invalid or points to a disallowed host.
    """
    try:
        domain_url = geardb._read_domain_url()
        if not domain_url:
            raise ValueError("Domain URL not configured. Cannot process track hub.")

        ALLOWED_HUB_DOMAINS = [
            domain_url,
            "localhost",
        ]

        parsed = urlparse(hub_url)

        # Ensure scheme is http or https
        if parsed.scheme not in ["http", "https"]:
            raise ValueError("Invalid URL scheme. Only http and https are allowed.")

        hostname = parsed.hostname
        if not hostname:
            raise ValueError("Invalid URL: missing hostname.")

        # Reject private/internal IP addresses
        try:
            ip = ipaddress.ip_address(hostname)
            if ip.is_private or ip.is_loopback or ip.is_link_local:
                raise ValueError(f"Access to private/internal IP address {hostname} is not allowed.")
        except ValueError as e:
            # Not an IP address, check domain whitelist
            error_msg = str(e)
            if "does not appear to be" not in error_msg and "is not a valid" not in error_msg:
                raise

            # Check against whitelist (case-insensitive)
            if not any(hostname.lower().endswith(domain.lower()) for domain in ALLOWED_HUB_DOMAINS):
                raise ValueError(f"Hub domain '{hostname}' is not in the allowed list.")

    except ValueError:
        raise


def _fetch_tracks_from_hub(hub_url: str, assembly: str) -> list[dict]:
    """
    Fetch and parse tracks from a track hub, handling both useOneFile and traditional modes.

    Args:
        hub_url (str): URL to the hub.txt file.
        assembly (str): Genome assembly identifier (e.g., 'hg38', 'mm10').

    Returns:
        list[dict]: List of track dictionaries.

    Raises:
        requests.RequestException: If hub files cannot be retrieved.
        ValueError: If assembly not found or parsing fails.
    """

    # Validate hub URL before making any requests
    # Addresses https://github.com/IGS/gEAR/security/code-scanning/344
    try:
        _validate_hub_url(hub_url)
    except ValueError as e:
        raise

    # Fix for Docker
    hub_url = hub_url.replace("http://localhost:8080", "http://web")

    try:
        hub_response = requests.get(hub_url)
        hub_response.raise_for_status()
    except requests.RequestException as e:
        raise ValueError(f"Failed to retrieve hub.txt from {hub_url}: {str(e)}") from e

    hub_txt = hub_response.text
    hub_json, tracks = parse_hub_from_file(hub_txt)

    # Traditional mode: parse genomes.txt and trackDb.txt
    if not hub_json.get("useOneFile", "") == "on":
        base_url = hub_url.rsplit("/", 1)[0]

        genomes_file_name = hub_json.get("genomesFile", "genomes.txt")
        genomes_url = f"{base_url}/{genomes_file_name}"

        try:
            genomes_response = requests.get(genomes_url)
            genomes_response.raise_for_status()
        except requests.RequestException as e:
            raise ValueError(
                f"Failed to retrieve genomes.txt from {genomes_url}: {str(e)}"
            ) from e

        trackdb_path = fetch_trackdb_path(genomes_response.text, assembly)
        trackdb_url = f"{base_url}/{trackdb_path}"

        try:
            trackdb_response = requests.get(trackdb_url)
            trackdb_response.raise_for_status()
        except requests.RequestException as e:
            raise ValueError(
                f"Failed to retrieve trackDb from {trackdb_url}: {str(e)}"
            ) from e

        tracks = parse_tracks_from_trackdb(trackdb_response.text, trackdb_url)

    return tracks


def build_assembly_array(assembly) -> list:
    """
    Builds a 2D array of chromosome names and sizes for a given genome assembly.
    Given an assembly identifier (e.g., 'hg38', 'mm10'), this function retrieves the corresponding chromosome sizes file,
    parses its contents, and returns a list of [chromosome, size] pairs. The chromosome sizes file is expected to be
    tab-delimited with at least three columns, where the first column is the chromosome name and the third column is the size.
    Args:
        assembly (str): The genome assembly identifier (e.g., 'hg38', 'mm10').
    Returns:
        list: A list of lists, where each sublist contains a chromosome name (str) and its size (int).
    Raises:
        ValueError: If the assembly is not supported, does not have a chromosome sizes file, or if the file cannot be retrieved.
    """

    chromosome_sizes_name = ASSEMBLY_TO_CHROMSIZES_FILE.get(assembly, None)
    if not chromosome_sizes_name:
        raise ValueError(
            f"Assembly {assembly} is not supported or does not have a chromosome sizes file."
        )

    # Strip out the port so we can use requests.get on it. Only applies to localhost reads
    genomes_base_url = GENOMES_ROOT
    if GENOMES_ROOT.startswith("http://localhost:"):
        parsed_url = urlparse(GENOMES_ROOT)
        genomes_base_url = f"{parsed_url.scheme}://{parsed_url.hostname}/{parsed_url.path.lstrip('/')}"

    chromosome_sizes_url = genomes_base_url + "/" + assembly + "/" + chromosome_sizes_name
    if not chromosome_sizes_url:
        raise ValueError(f"Chromosome sizes URL for assembly {assembly} is not valid.")

    assembly_array = []

    # Get the data from the chromosome sizes url
    response = requests.get(chromosome_sizes_url)
    if response.status_code != 200:
        raise ValueError(
            f"Failed to retrieve chromosome sizes from {chromosome_sizes_url}"
        )

    for line in response.text.splitlines():
        chromosome, size = line.strip().split(
            "\t"
        )  # my chrom sizes file has 0 in the middle column
        assembly_array.append([chromosome, int(size)])
    return assembly_array


def build_assembly_gos_obj(assembly_array) -> gos.Assembly:
    """
    Builds and returns a Gosling Assembly object for the given genome assembly.

    Args:
        assembly_array: A 2D array where each sub-array contains a chromosome name and its size.

    Returns:
        gos.Assembly: An object containing chromosome size information for the specified assembly.
    """

    chrom_sizes = gos.ChromSizes(assembly_array)
    return gos.Assembly(chrom_sizes)

def build_assembly_json_from_array(assembly_array) -> list:
    """
    Builds and returns a Gosling Assembly JSON object for the given genome assembly.

    Args:
        assembly_array: A 2D array where each sub-array contains a chromosome name and its size.

    Returns:
        list: A list of dictionaries representing chromosome size information for the specified assembly.
    """

    values = []
    for chromosome, size in assembly_array:
        values.append({"chrom": chromosome, "chromStart": 1, "chromEnd": size})
    return values

def build_bed_annotation_tracks(assembly, zoom=False, title="left"):
    """
    Builds a Gosling view for gene and exon annotation tracks based on BED files for a given genome assembly.

    Parameters
    ----------
    assembly : str
        The genome assembly identifier (e.g., "hg38", "mm10").
    zoom : bool, optional
        If True, returns zoomed annotation tracks using `build_bed_zoomed_annotation_tracks`. Default is False.
    title : str, optional
        The title or identifier for the annotation view. Default is "left".

    Returns
    -------
    gosling.View
        A Gosling overlay view containing gene and exon annotation tracks for the specified assembly.

    Raises
    ------
    ValueError
        If the specified assembly is not supported or does not have a corresponding BED file.

    Notes
    -----
    - Uses Ensembl-style BED files for genes and exons, sorted in chromosome order.
    - If exon BED file is missing for the assembly, only gene tracks are returned and a warning is printed.
    - The view includes text, rule, and rect tracks for gene names, gene regions, and exon regions, respectively.
    """

    # These files will be based on Ensembl's annotation naming structure, but sorted in chromosome order.
    # I am adding a 2nd column of 1's to allow us to use the files in a "genomic" track.

    ANNOTATION_BEDDB_UID = {
        "danRer10": "YwOpmCgUSqKdJSGtGXWeJw", # zebrafish
        #"galGal6": "galGal6.annotation.beddb", # chicken - could not get to load
        "hg19": "XXcPeaTRSiy8_yxNwjtzEQ",
        "hg38": "GhiCXRRHTH2u-24jRq0HRQ",
        "mm10": "VNbLgNO3T8uAcp_5vRFqdQ",
        # "mm39": "mm39.annotation.beddb",
        "rn6": "C6Tw-g54Rl602PqE62qTHw", # rat
        # "calJac3": "calJac3.annotation.beddb", # marmoset
    }

    ANNOTATION_CONDENSED_HEIGHT = EXPANDED_HEIGHT # Always taller for annotations
    #ANNOTATION_EXPANDED_HEIGHT = ANNOTATION_CONDENSED_HEIGHT * 2
    ANNOTATION_EXPANDED_HEIGHT = ANNOTATION_CONDENSED_HEIGHT

    beddb_uid = ANNOTATION_BEDDB_UID.get(assembly, None)
    if not beddb_uid:
        raise ValueError(
            f"Assembly {assembly} is not supported or does not have a BEDdb UID."
        )

    beddb_url = f"https://higlass.umgear.org/api/v1/tileset_info/?d={beddb_uid}"

    # chrom, start, end, name, score, strand, exon_start, exon_end
    beddb_data = gos.beddb(
        url=beddb_url,
        genomicFields=[
            {"index": 1, "name": "start"},
            {"index": 2, "name": "end"}
        ],
        valueFields=[
            {"index": 5, "name": "strand", "type": "nominal"},
            {"index": 3, "name": "name", "type": "nominal"}
        ],
        exonIntervalFields=[
            {"index": 6, "name": "exon_start"},
            {"index": 7, "name": "exon_end"}
        ]
    )

    #MINUS_ROW_POSITION = 50 if zoom else 35
    MINUS_ROW_POSITION = 25

    # These are shared amongst the genes and exons track
    condensed_row = gos.Row(field="strand", type="nominal", domain=["+", "-"], range=[0, MINUS_ROW_POSITION])  # type:ignore
    expanded_row = gos.Row(field="displace_row", type="nominal")    # type: ignore
    #row = expanded_row if zoom else condensed_row
    row=condensed_row

    color = gos.Color(field="strand", type="nominal", domain=["+", "-"], range=["darkblue", "darkred"])  # type:ignore
    tooltip=[
                gos.Tooltip(field="start", type="genomic", alt="Start Position"),  # type: ignore
                gos.Tooltip(field="end", type="genomic", alt="End Position"),  # type: ignore
                gos.Tooltip(field="strand", type="nominal", alt="Strand"),  # type: ignore
                gos.Tooltip(field="name", type="nominal", alt="Name"),  # type: ignore
            ]


    base_track =  (
        gos.Track(data=beddb_data) # type: ignore
            .encode(
                x=gos.X(field="start", type="genomic"),  # type:ignore
                xe=gos.X(field="end", type="genomic"),  # type:ignore
                row=row,
                color=color,
                tooltip=tooltip,
            )
            .visibility_lt(
                measure="width", threshold="|xe-x|", transitionPadding=10, target="mark"
            )
    )

    gene_track = base_track.transform_filter(field="type", oneOf=["gene"])

    """
    if zoom:
        gene_track = (
            gene_track
                .transform_displace(
                    method="pile",
                    boundingBox={"startField": "start", "endField": "end"},
                    newField="displace_row",
                    maxRows=4,
                )
            )
    """

    # Three tracks
    # 1) Text track - uses bed_track
    # 2) Rule track (for tooltips) - uses bed_track
    # 3) Rect track (for exons in a separate BED file)

    text_track = gene_track.mark_text().encode(
        text=gos.Text(field="name", type="nominal"),  # type: ignore
        style=gos.Style(dy=-10),  # type: ignore
    )

    gene_rule_track = (
        gene_track.mark_rule()
        .encode(
            strokeWidth=gos.StrokeWidth(value=5),
        )
        .properties(id=f"{title}-annotation")   # used for zooming to specific region
    )

    exon_track = (
        base_track.mark_rect()
            .encode(
                x=gos.X(field="exon_start", type="genomic"),  # type:ignore
                xe=gos.X(field="exon_end", type="genomic"),  # type:ignore
                size=gos.Size(value=10)
                )
            .transform_filter(field="type", oneOf=["exon"])
        )

    """
    if zoom:
        exon_track = exon_track.transform_displace(
                method="pile",
                boundingBox={"startField": "exon_start", "endField": "exon_end"},
                newField="displace_row",
                maxRows=4,
            )
    """


    tracks = [text_track, gene_rule_track, exon_track]

    panel_title = "Annotation"
    if zoom:
        panel_title = "Panel B" if title == "right" else "Panel A"

    annotation_view = gos.overlay(*tracks).properties(
        title=panel_title,
        width=EXPANDED_WIDTH if zoom else CONDENSED_WIDTH,
        height=ANNOTATION_EXPANDED_HEIGHT if zoom else ANNOTATION_CONDENSED_HEIGHT,
        opacity=gos.Opacity(value=0.8),
        id=f"{title}-annotation-view",
    )

    return annotation_view

def build_genome_wide_view(
    assembly, assembly_array, zoom=False, position: str | None = None, gene_symbol: str | None = None
) -> gos.View:
    """
    Build a Gosling genome-wide view track for a given genome assembly.

    This function creates a visualization track that displays all chromosomes for the specified assembly,
    using chromosome size information from a remote TSV file. The view includes a chromosome
    representation and interactive brush marks for zooming and selection.

    Args:
        assembly (str): The genome assembly identifier (e.g., "hg38", "mm10").
        zoom (bool, optional): If True, adds an additional brush for a secondary zoom panel. Defaults to False.
        position (str, optional): If provided, restricts the view to the specified chromosome (e.g., "chr1").

    Returns:
        gos.View: A Gosling view object representing the genome-wide track.
    """

    assembly_json = build_assembly_json_from_array(assembly_array)

    chrom = None
    start = None
    end = None
    if position is not None:
        _, chrom, start, end = parse_position_str(position)

    data = gos.json(
        chromosomeField="chrom",
        genomicFields=["chromStart", "chromEnd"],
        values=assembly_json,
    )

    base = gos.Track(data)  # type: ignore

    # Set track title
    title = f"{assembly} assembly"
    if chrom:
        title = f"Regional view around {gene_symbol}" if gene_symbol else f"Chromosome {chrom}"
    if zoom and chrom:
        title = f"Regional view around {gene_symbol}" if gene_symbol else f"Panel A chromosome {chrom}"

    # Genome-wide brush should be easier to see that chromosome brush
    brushwidth = 2 if chrom else 5

    # build *tracks for gos.overlay
    # The Gos documentation notes it is more idiomatic to set specific properties after initialization
    tracks = [
        base.mark_rect()
        .encode(
            x=gos.X(field="chromStart", type="genomic"),  # type:ignore
            xe=gos.X(field="chromEnd", type="genomic"),  # type:ignore
            color=gos.Color(field="chrom", type="nominal", range=["#666666", "#999999"]),  # type:ignore
            size=gos.Size(value=20),  # type: ignore
        )
        .properties(
            title=title,
        ),
        base.mark_brush().encode(
            x=gos.X(linkingId="zoom-to-panel-a"),  # type:ignore
            color=gos.Color(value="steelblue"),
            stroke=gos.Stroke(value="steelblue"),
            strokeWidth=gos.StrokeWidth(value=brushwidth),  # type: ignore
        ),
    ]
    if zoom:
        tracks.append(
            base.mark_brush().encode(
                x=gos.X(linkingId="zoom-to-panel-b"),  # type:ignore
                color=gos.Color(value="yellow"),
                stroke=gos.Stroke(value="yellow"),
                strokeWidth=gos.StrokeWidth(value=brushwidth),  # type: ignore
            )
        )

    genome_wide_view = gos.overlay(*tracks).properties(
        static=True,
        #static=False if chrom else True,
        id="chromosome-wide" if chrom else "genome-wide",
        width=CONDENSED_WIDTH,
        height=20,
    )

    if chrom:
        # restrict xDomain to the desired interval on chromosome to make the brush span visibly.
        # Fixes to the left-panel chromosome.
        padding = 5e5  # 50Kb padding on each side
        interval = None
        if start and end:
            # Adjust interval based on a padding.
            # Ending should not extend beyond the chromosome
            start = int(max(0, start - padding))

            end = end + padding
            chrom_end = assembly_json[
                next(i for i, v in enumerate(assembly_json) if v["chrom"] == chrom)
            ]["chromEnd"]
            end = int(min(chrom_end, end))

            interval = [start, end]

            genome_wide_view = genome_wide_view.properties(
                xDomain=gos.GenomicDomain(chromosome=chrom, interval=interval)  # type: ignore
            )


    """
    NOTE: There was a thought of linking the genome-wide brush to the chromosome brush.
    However, the 'betweenLink' mark does not work between two tracks of different scale.
    It also requires its own data source, which would have to be from one of the scales.
    """

    return genome_wide_view


def build_gosling_tracks(parent_tracks_dict, tracks, zoom=False, position_str="NA"):
    """
    Builds and configures Gosling tracks based on the provided track specifications.

    Args:
        parent_tracks_dict (dict): A dictionary with keys "left" and optionally "right", each mapping to a list of track objects.
        tracks (list): A list of dictionaries, each representing a track specification. Each track dict should contain at least:
            - "type" (str): The type of the track (e.g., "bam", "bigWig", "bed", "vcf").
            - "bigDataUrl" (str): The URL to the data file for the track.
            - Optional: "color" (str), and other track-specific attributes.
        zoom (bool, optional): If True, builds both left and right tracks for zoomed-in views. Defaults to False.

    Returns:
        tuple:
            - parent_view_left: The Gosling  view for the left panel.
            - parent_view_right: The Gosling  view for the right panel if zoom is True, otherwise None.

    Notes:
        - Unsupported track types or tracks missing "bigDataUrl" are skipped with a warning.
        - The function uses specific builder classes for each supported track type.
        - Tracks are configured with properties such as title, width, height, opacity, and visibility.
    """

    TRACK_TYPE_2_SPEC = {
        "bam": BamSpec,
        "bigWig": BigWigSpec,
        "bigBed": BedSpec,
        "vcfTabix": VcfSpec,
        "hic": HiCSpec,
    }

    kwargs = {}
    hic_found = False

    # Build each individual track based on its type
    for track in tracks:
        track_type = track.get("type", "")
        spec_builder_class = TRACK_TYPE_2_SPEC.get(track_type, None)
        if not spec_builder_class:
            print(
                f"WARNING: Unsupported track type '{track_type}' for track '{track.get('shortLabel', '')}'; skipping.",
                file=sys.stderr,
            )
            continue

        # Resolve URL for Gosling context
        data_url = _resolve_track_url(track, use_gosling=True)
        if not data_url:
            print(
                f"WARNING: Could not resolve URL for track '{track.get('shortLabel', '')}'; skipping.",
                file=sys.stderr,
            )
            continue

        # Get other attributes to pass to the class
        color = track.get("color", "orange")  # Default color if not specified

        # Title should be based on shortLabel, longLabel, bigDataUrl (in that order)
        title = track.get("shortLabel", track.get("longLabel", "bigDataUrl"))

        visibility = track.get("visibility", "dense")  # Dense, full, hide

        try:
            spec_builder = spec_builder_class(
                data_url=data_url, color=color, zoom=zoom, title=title, position_str=position_str, visibility=visibility
            )
            left_track = spec_builder.add_track(**kwargs)

            # Skip if track validation failed (add_track returns None)
            if left_track is None:
                continue

            parent_tracks_dict["left"].append(left_track)

            if spec_builder_class == HiCSpec:
                hic_found = True

            if zoom:
                right_track = spec_builder.add_track(**kwargs)
                if right_track is not None:
                    right_track.id = f"right-track-{Path(data_url).stem}"
                    parent_tracks_dict["right"].append(right_track)
        except Exception as e:
            print(f"ERROR building track '{title}': {e}", file=sys.stderr)
            continue

    parent_view_left = gos.stack(*parent_tracks_dict["left"]).properties(
        id="left-view", linkingId="zoom-to-panel-a", spacing=0,
    )

    parent_view_right = None
    if zoom:
        parent_view_right = gos.stack(*parent_tracks_dict["right"]).properties(
            id="right-view", linkingId="zoom-to-panel-b", spacing=0
        )

    return parent_view_left, parent_view_right, hic_found


def build_region_view(parent_view_left, parent_view_right=None):
    """
    Build a Gosling region view track for a given genome assembly.

    This function creates a visualization track that displays detailed genomic regions
    for the specified assembly, using chromosome size information from a remote TSV file.
    The view includes interactive brush marks for zooming and selection, and can optionally
    include a secondary zoom panel.

    Args:
        assembly (str): The genome assembly identifier (e.g., "hg38", "mm10").
        zoom (bool, optional): If True, adds an additional brush for a secondary zoom panel. Defaults to False.
        parent_view_left (gos.View): The left view containing tracks to be displayed in the region view.
        parent_view_right (gos.View, optional): The right view containing tracks for zoomed-in display. Defaults to None.
    Returns:
        gos.View: A Gosling view object representing the region track.
    """

    args = [parent_view_left]
    if parent_view_right:
        args.append(parent_view_right)

    region_view = gos.horizontal(*args)

    return region_view

def get_gene_info(gene_symbol, dataset_id, assembly):
    """
    Fetches gene coordinate information for a given gene symbol from the geardb database.

    Args:
        gene_symbol (str): The gene symbol to look up.
        dataset_id (str or int): The identifier for the dataset to query.
        assembly (str): The genome assembly name (e.g., 'hg38', 'mm10').

    Returns:
        tuple:
            - position_str (str): The genomic position in the format '{assembly}.chr:start-end', suitable for UCSC Genome Browser.
            - gene_symbol (str): The canonical gene symbol as stored in the database.

    Notes:
        - If the gene symbol is not found in the dataset, returns ("NA", "NA").
        - Chromosome names are normalized to UCSC-style (e.g., 'chr1', 'chrX', 'chrM').
    """
    # Fetch gene coordinates from geardb
    gene_info = geardb.get_gene_by_gene_symbol(gene_symbol, dataset_id)
    if not gene_info:
        print(
            f"WARNING: Gene symbol '{gene_symbol}' not found in dataset {dataset_id}; cannot zoom track.",
            file=sys.stderr,
        )
        return ("NA", "NA")

    chrom = gene_info.molecule or "unknown"

    # if chrom is a number, convert it to a string with "chr" prefix
    if chrom.isdigit() or chrom in ["X", "Y"]:
        chrom = f"chr{chrom}"
    # if chrom is MT, convert to "chrM"
    if chrom == "MT":
        chrom = "chrM"

    start = gene_info.start
    end = gene_info.stop
    gene_symbol = gene_info.gene_symbol # To use correct naming later

    left_position = f"{chrom}:{start}-{end}"
    position_str = f"{assembly}.{left_position}"  # This is the format if we want to export position to UCSC Genome Browser

    return position_str, gene_symbol

def parse_position_str(position_str: str) -> tuple:
    """
    Parses a position string in the format 'assembly.chromosome:start-end' and returns its components.

    Args:
        position_str (str): The position string to parse.

    Returns:
        tuple: A tuple containing (assembly, chromosome, start, end) if parsing is successful,
            otherwise (None, None, None, None) on failure.
    """
    try:
        assembly_chrom, interval = position_str.split(".")
        chromosome, positions = interval.split(":")
        start, end = positions.split("-")
        return assembly_chrom, chromosome, int(start), int(end)
    except Exception as e:
        print(f"ERROR: Failed to parse position string '{position_str}': {e}", file=sys.stderr)
        return None, None, None, None

def zoom_view_to_domain(view, position_str, hic_found=False):
    """
    Zooms a Gosling view to a specified genomic domain with padding.

    Args:
        view: A Gosling view object to be updated.
        position_str (str): A string representing the genomic position in the format expected by `parse_position_str`.

    Returns:
        The updated Gosling view object with its xDomain set to the specified genomic coordinates plus padding.

    Notes:
        Adds a base padding of 1500 base pairs to both sides of the specified genomic interval.
    """

    BASE_PADDING = 1500  # base padding on each side of gene
    padding = BASE_PADDING * 500 if hic_found else BASE_PADDING

    _, chrom, start, end = parse_position_str(position_str)

    # Set the x domain of the track to the gene coordinates
    view = view.properties(
        xDomain=gos.GenomicDomain(
            chromosome=chrom, interval=[start - padding, end + padding]
        )
    )
    return view

class TrackSpec(ABC):
    def __init__(self, data_url, color="steelblue", zoom=False, title="", position_str="NA", visibility="full"):
        self.data_url = data_url
        self.color = color  # Passed as RGB string
        self.zoom = zoom
        self.width = EXPANDED_WIDTH if zoom else CONDENSED_WIDTH
        self.height = EXPANDED_HEIGHT if zoom else CONDENSED_HEIGHT
        self.title = title
        self.position_str = position_str
        self.visibility = visibility
        self.track = None

        if self.visibility == "hide":
            self.height = 0 # effectively hide it
        elif self.visibility == "full":
            self.height *= 1.5

    @abstractmethod
    def add_track(self, **kwargs):
        pass

    def validate_track_url(self, url: str, expected_extensions: list[str]) -> None:
        """
        Validate track URL with scheme, extension, and accessibility checks.

        Performs checks in order of cost (cheapest first) to fail fast:
        1. Scheme validation (no network)
        2. Extension validation (no network)

        Args:
            url (str): URL to validate.
            expected_extensions (list[str]): List of valid extensions (e.g., ['.bam', '.bw']).

        Raises:
            ValueError: If any validation fails.
        """
        # Cheap checks first
        self._validate_url_scheme(url)
        self._validate_url_extension(url, expected_extensions)

    @staticmethod
    def _validate_url_scheme(url: str) -> None:
        """Common URL scheme validation."""
        if not (url.startswith("http://") or url.startswith("https://")):
            raise ValueError("Invalid URL: must start with http:// or https://")

    @staticmethod
    def _validate_url_extension(url: str, valid_extensions: list[str]) -> None:
        """Validate URL ends with one of the valid extensions."""
        if not any(url.endswith(ext) for ext in valid_extensions):
            extensions_str = ", ".join(valid_extensions)
            raise ValueError(f"Invalid URL: must end with {extensions_str}")

class BamSpec(TrackSpec):
    def add_track(self, **kwargs):
        url = self.data_url
        color = self.color

        try:
            self.validate_track_url(url, [".bam"])
            self.validate_track_url(f"{url}.bai", [".bam.bai"])
        except ValueError as e:
            print(f"WARNING: Skipping BAM track '{self.title}': {e}", file=sys.stderr)
            return None

        bamData = gos.bam(url=url, indexUrl=f"{url}.bai")

        track = (
            gos.Track(
                data=bamData,  # pyright: ignore[reportArgumentType]
                width=self.width,
                height=self.height,
                title=self.title,  # Use the file name as the title
                id=f"left-track-{self.title}",  # Use the file name without extension as the ID
            )
            .mark_bar()
            .encode(
                x=gos.X(field="start", type="genomic"),  # pyright: ignore[reportArgumentType]
                xe=gos.X(field="end", type="genomic"),  # pyright: ignore[reportArgumentType]
                y=gos.Y(field="coverage", type="quantitative", axis="right"),  # pyright: ignore[reportArgumentType]
                color=gos.Color(value=color),
            )
        )
        return track


class BedSpec(TrackSpec):
    def add_track(self, **kwargs):
        url = self.data_url
        color = self.color

        try:
            self.validate_track_url(url, [".bed.gz"])
            self.validate_track_url(f"{url}.tbi", [".bed.gz.tbi"])
        except ValueError as e:
            print(f"WARNING: Skipping BED track '{self.title}': {e}", file=sys.stderr)
            return None

        bed_data = gos.bed(url=url, indexUrl=f"{url}.tbi")

        track = (
            gos.Track(
                data=bed_data,  # pyright: ignore[reportArgumentType]
                width=self.width,
                height=self.height,
                title=self.title,  # Use the file name as the title
                id=f"left-track-{self.title}",  # Use the file name without extension as the ID
            )
            .mark_rect()
            .encode(
                x=gos.X(field="start", type="genomic", axis="none"),  # pyright: ignore[reportArgumentType]
                xe=gos.X(field="end", type="genomic"),  # pyright: ignore[reportArgumentType]
                size=gos.Size(value=10),
                color=gos.Color(value=color),
            )
        )
        return track


class BigWigSpec(TrackSpec):

    def add_track(self, **kwargs):
        url = self.data_url
        color = self.color

        try:
            self.validate_track_url(url, [".bw", ".bigwig"])
        except ValueError as e:
            print(f"WARNING: Skipping BigWig track '{self.title}': {e}", file=sys.stderr)
            return None

        # TODO: figure out appropriate binsize when zooming out: Default is 256. It looks blocky briefly
        bigwig_data = gos.bigwig(url=url)

        y_kwargs = {}
        y_field="value"

        track = (
            gos.Track(
                data=bigwig_data,  # pyright: ignore[reportArgumentType]
                width=self.width,
                height=self.height,
                title=self.title,  # Use the file name as the title
                id=f"left-track-{self.title}",  # Use the file name without extension as the ID
            )
            .mark_bar()
            .encode(
                x=gos.X(field="start", type="genomic", axis="none"),  # pyright: ignore[reportArgumentType]
                xe=gos.X(field="end", type="genomic"),  # pyright: ignore[reportArgumentType]
                y=gos.Y(field=y_field, type="quantitative", axis="right", aggregate="count", **y_kwargs),  # pyright: ignore[reportArgumentType]
                color=gos.Color(value=color),
                # Tooltip only works when hovering over the peak
                tooltip=[
                    gos.Tooltip(field="position", type="genomic", alt="Position"),  # type: ignore
                    gos.Tooltip(field="value", type="quantitative", alt="Peak value", format=".2"),  # type: ignore

                ],
            )
        )
        return track

class VcfSpec(TrackSpec):
    def add_track(self, **kwargs):
        url = self.data_url
        color = self.color

        try:
            self.validate_track_url(url, [".vcf.gz"])
            self.validate_track_url(f"{url}.tbi", [".vcf.gz.tbi"])
        except ValueError as e:
            print(f"WARNING: Skipping VCF track '{self.title}': {e}", file=sys.stderr)
            return None

        vcf_data = gos.VcfData(type="vcf", url=url, indexUrl=f"{url}.tbi")  # type: ignore

        track = (
            gos.Track(
                data=vcf_data,  # pyright: ignore[reportArgumentType]
                width=self.width,
                height=self.height,
                title=self.title,  # Use the file name as the title
                id=f"left-track-{self.title}",  # Use the file name without extension as the ID
            )
            .mark_point()
            .encode(
                x=gos.X(field="position", type="genomic", axis="none"),  # pyright: ignore[reportArgumentType]
                y=gos.Y(field="value", type="quantitative", axis="right"),  # pyright: ignore[reportArgumentType]
                color=gos.Color(value=color),
            )
        )
        return track

class HiCSpec(TrackSpec):
    def add_track(self, **kwargs):
        url = self.data_url
        #color = self.color  # colorscale instead of single color
        position_str = self.position_str

        """Accepted HiGlass color ranges (others will not work)
            viridis: interpolateViridis,
            grey: interpolateGreys,
            warm: interpolateWarm,
            spectral: interpolateSpectral,
            cividis: interpolateCividis,
            bupu: interpolateBuPu,
            rdbu: interpolateRdBu,
            hot: interpolateYlOrBr,
            pink: interpolateRdPu
        """

        #try:
        # TODO: This requires the HiGlass mcool path
        #    self.validate_track_url(url)
        #except ValueError:
        #    raise

        hic_data = gos.matrix(url=url)

        hic_track = (
            gos.Track(
                data=hic_data,  # pyright: ignore[reportArgumentType]
                width=self.width,
                height=self.width,  # looks best as square aspect ratio
                title=self.title,  # Use the file name as the title
            )
            .mark_bar()
            .encode(
                x=gos.X(field="xs", type="genomic", axis="bottom"),  # pyright: ignore[reportArgumentType]
                xe=gos.Xe(field="xe", type="genomic", axis="none"),  # pyright: ignore[reportArgumentType]
                y=gos.Y(field="ys", type="genomic", axis="none"),  # pyright: ignore[reportArgumentType]
                ye=gos.Ye(field="ye", type="genomic", axis="none"),  # pyright: ignore[reportArgumentType]
                color=gos.Color(field="value", type="quantitative", range="bupu", legend=True),  # pyright: ignore[reportArgumentType]
                style=gos.Style(matrixExtent="full"), # pyright: ignore[reportArgumentType]
            )
        )
        if position_str == "NA":
            return hic_track

        annotation_track = self.add_annotation_track(position_str)

        view = gos.overlay(hic_track, annotation_track, width=self.width, height=self.width, id=f"left-track-{self.title}",  # Use the file name without extension as the ID
)
        return view


    def add_annotation_track(self, position_str):

        _, chrom, start, end = parse_position_str(position_str)

        json_data = gos.json(
            values=[
                {
                    "c": chrom,
                    "x": start,
                    "xe": end,
                    "y": start,
                    "ye": end,
                }
            ],
            chromosomeField="c",
            genomicFields=["x", "xe", "y", "ye"],
            )

        track = (
            gos.Track(
                data=json_data, # type: ignore
            )
            .mark_bar()
            .encode(
                x=gos.X(field="x", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
                xe=gos.Xe(field="xe", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
                y=gos.Y(field="y", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
                ye=gos.Ye(field="ye", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
                color=gos.Color(value="yellow"),
                opacity=gos.Opacity(value=0.2),
                stroke=gos.Stroke(value="yellow"),
                strokeWidth=gos.StrokeWidth(value=4),
            )
        )

        return track

# Assembly track class
class AssemblySpec:
    pass


#####################


class GoslingSpec(Resource):
    def get(self, dataset_id):
        # session_id = request.cookies.get("gear_session_id", "")
        args = request.args
        gene_symbol = args.get("gene")
        assembly = args.get("assembly")
        hub_url = args.get("hub_url", "")
        zoom = args.get("zoom", "false")
        # Set zoom from string to bool
        if zoom.lower() == "true":
            zoom = True
        else:
            zoom = False

        response = {
            "success": 0,
            "spec": {},
            "position": "",  # will be filled with chr:start-end
            "message": "",
            "hic_found": False  # If the spec has a HiC file, we need to adjust zooming behavior
        }

        if assembly is None:
            response["message"] = "Missing required parameter: assembly"
            return response, 400

        if gene_symbol is None:
            response["message"] = "Missing required parameter: gene"
            return response, 400

        if hub_url == "":
            response["message"] = "Missing required parameter: hub_url"
            return response, 400

        # Sanity check to confirm dataset exists
        dataset = geardb.get_dataset_by_id(d_id=dataset_id, include_shape=False)
        if not dataset:
            response["message"] = f"Dataset with ID {dataset_id} not found."
            return response, 404

        try:
            # Fetch and parse hub files
            tracks = _fetch_tracks_from_hub(hub_url, assembly)

        except (ValueError, requests.RequestException) as e:
            response["message"] = str(e)
            return response, 400

        # Let's get the coordinates of the gene_symbol.
        (position_str, new_gene_symbol) = get_gene_info(gene_symbol, dataset_id, assembly)

        if new_gene_symbol != "NA":
            gene_symbol = new_gene_symbol

        gos_tracks = {"left": [], "right": []}
        # Add BED annotation tracks to left (and right if zoom) gos_tracks in the first index position
        # build left and right track
        # Insert into gos_tracks["left"] at index 0 (and gos_tracks["right"] if zoom)
        gos_tracks["left"].append(build_bed_annotation_tracks(assembly, zoom, "left"))
        if zoom:
            gos_tracks["right"].append(
                build_bed_annotation_tracks(assembly, zoom, "right")
            )

        (parent_view_left, parent_view_right, hic_found) = build_gosling_tracks(
            gos_tracks, tracks, zoom=zoom, position_str=position_str
        )

        # Start building the Gosling spec
        region_view = build_region_view(parent_view_left, parent_view_right)

        base_views = []

        assembly_array = build_assembly_array(assembly)

        genome_wide_view = build_genome_wide_view(assembly, assembly_array, zoom=zoom, gene_symbol=gene_symbol)
        base_views.append(genome_wide_view)

        # Let's do some things that depend on whether we found a valid position_str
        if position_str != "NA":
            chromosome_view = build_genome_wide_view(assembly, assembly_array, zoom=zoom, position=position_str, gene_symbol=gene_symbol)
            base_views.append(chromosome_view)

            region_view = zoom_view_to_domain(region_view, position_str, hic_found)


        base_views.append(region_view)
        base_track = gos.vertical(*base_views)

        # Add assembly track to base track
        assembly_obj = build_assembly_gos_obj(assembly_array)
        base_track = base_track.properties(
            assembly=assembly_obj,
        )

        spec = base_track.to_json(indent=2)

        response["success"] = 1
        response["spec"] = spec
        response["position"] = position_str
        response["hic_found"] = hic_found
        return response, 200

    def post(self, dataset_id):
        # Implement your logic to create or update the Gosling spec for the given dataset_id
        pass
        # return view.save("gosling.json")

    def delete(self, dataset_id):
        # Implement your logic to delete the Gosling spec for the given dataset_id
        pass

    def put(self, dataset_id):
        # Implement your logic to replace the Gosling spec for the given dataset_id
        pass
