import json
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from urllib.parse import urljoin

import geardb
import gosling as gos
import requests
from flask import request
from flask_restful import Resource


"""
NOTE: The documentation kind of sucks, and requires some trial and error to figure out where things go.
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
JSON_DIR = abs_path_www.joinpath("tracks")  # Directory for JSON track files
VIEW_PADDING = 10  # Padding between views
CONDENSED_WIDTH = 1200 - VIEW_PADDING  # Width of the condensed view tracks
EXPANDED_WIDTH = 600 - VIEW_PADDING/2  # Width of the left view tracks

CONDENSED_HEIGHT = 20  # Height for condensed tracks
EXPANDED_HEIGHT = 40  # Height for expanded tracks

# These files will be based on Ensembl's annotation naming structure, but sorted in chromosome order.
# I am adding a 2nd column of 1's to allow us to use the files in a "genomic" track.
ASSEMBLY_TO_CHROMSIZES_FILE = {
    "danRer10": "danRer10.chromInfo.txt", # zebrafish
    "galGal6": "galGal6.chromInfo.txt", # chicken
    #"hg19": "hg19.chromInfo.txt",
    "hg38": "hg38.chromInfo.txt",
    "mm10": "mm10.chromInfo.txt",
    #"mm39": "mm39.chromInfo.txt",
    "r6": "r6.chromInfo.txt", # rat
    #"calJac3": "calJac3.chromInfo.txt", # marmoset
    }


def build_assembly_array(assembly):
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
        raise ValueError(f"Assembly {assembly} is not supported or does not have a chromosome sizes file.")

    tracks_db_root = "https://umgear.org/tracks"  # Replace with actual tracks database
    chromosome_sizes_url = tracks_db_root + "/" + chromosome_sizes_name
    if not chromosome_sizes_url:
        raise ValueError(f"Chromosome sizes URL for assembly {assembly} is not valid.")

    assembly_array = []

    # Get the data from the chromosome sizes url
    response = requests.get(chromosome_sizes_url)
    if response.status_code != 200:
        raise ValueError(f"Failed to retrieve chromosome sizes from {chromosome_sizes_url}")

    for line in response.text.splitlines():
        chromosome, _, size = line.strip().split("\t")  # my chrom sizes file has 0 in the middle column
        assembly_array.append([chromosome, int(size)])
    return assembly_array

def build_assembly_gos_obj(assembly):
    """
    Builds and returns a Gosling Assembly object for the given genome assembly.

    Args:
        assembly: The genome assembly identifier or data used to build the assembly array.

    Returns:
        gos.Assembly: An object containing chromosome size information for the specified assembly.
    """

    assembly_array = build_assembly_array(assembly)

    chrom_sizes = gos.ChromSizes(assembly_array)
    return gos.Assembly(chrom_sizes)

def build_bed_annotation_tracks(assembly, zoom=False, title="left"):
    """
    Builds a set of Gosling tracks for visualizing gene annotations and exons from BED files for a given genome assembly.

    This function generates a composite Gosling view containing:
        1. A text track displaying gene names.
        2. A rule track for tooltips showing gene information.
        3. (If available) A rect track displaying exons from a separate BED file.

    The tracks are color-coded by strand and arranged in rows by strand. The function supports only assemblies for which BED and exon BED files are defined.

    Args:
        assembly (str): The genome assembly identifier (e.g., "mm10").
        zoom (bool, optional): If True, uses expanded width for the view. Defaults to False.
        title (str, optional): Position of the panel title ("left" or "right"). Defaults to "left".

    Returns:
        gos.View: A Gosling overlay view containing the annotation tracks.

    Raises:
        ValueError: If the specified assembly does not have a corresponding BED file.

    Notes:
        - If the exon BED file is not available for the assembly, only the text and rule tracks are returned.
        - The function prints a warning to stderr if the exon BED file is missing.
        - The URLs for BED files are constructed based on a predefined root path.
    """


    # These files will be based on Ensembl's annotation naming structure, but sorted in chromosome order.
    # I am adding a 2nd column of 1's to allow us to use the files in a "genomic" track.

    # TODO: make the other files
    ASSEMBLY_TO_BED_FILE = {
        #"danRer10": "danRer10.bed", # zebrafish
        #"galGal6": "galGal6.bed", # chicken
        #"hg19": "hg19.bed",
        #"hg38": "hg38.bed",
        "mm10": "knownGeneM20.bed.gz",
        #"mm39": "mm39.chromInfo.txt",
        #"r6": "r6.chromInfo.txt", # rat
        #"calJac3": "calJac3.chromInfo.txt", # marmoset
        }

    # TODO: make the other files
    ASSEMBLY_TO_EXON_FILE = {
        #"danRer10": "danRer10_exons.bed", # zebrafish
        #"galGal6": "galGal6_exons.bed", # chicken
        #"hg19": "hg19_exons.bed",
        #"hg38": "hg38_exons.bed",
        "mm10": "Mus_musculus.GRCm38.102.exons.bed.gz",
        #"mm39": "mm39_exons.bed",
        #"r6": "r6_exons.bed", # rat
        #"calJac3": "calJac3_exons.bed", # marmoset
    }

    bed_file_stem = ASSEMBLY_TO_BED_FILE.get(assembly, None)
    if not bed_file_stem:
        raise ValueError(f"Assembly {assembly} is not supported or does not have a BED file.")

    bed_tbi_file_stem = bed_file_stem + ".tbi"
    tracks_db_root = "https://umgear.org/tracks"  # Replace with actual tracks database
    bed_file_url = tracks_db_root + "/" + bed_file_stem
    bed_tbi_file_url = tracks_db_root + "/" + bed_tbi_file_stem

    # base bed track to build overlays on
    bed_data = gos.BedData(
        url=bed_file_url, # type: ignore
        indexUrl=bed_tbi_file_url,  # type: ignore
        type="bed", # type: ignore
        customFields=[
                "name2",
                "cdsStartStat",
                "cdsEndStat",
                "exonFrames",
                "type",
                "geneName"
                ] # type: ignore
    )

    # These are shared amongst the genes and exons track
    x=gos.X(field="chromStart", type="genomic") # type:ignore
    xe=gos.X(field="chromEnd", type="genomic") # type:ignore
    row=gos.Row(field="strand", type="nominal", domain=[1, -1], range=[0,20]) # type:ignore
    color=gos.Color(field="strand", type="nominal", domain=[1, -1], range=["darkblue", "darkred"]) # type:ignore

    gene_track = gos.Track(
        data=bed_data  # type: ignore
    ).encode(
        x=x,
        xe=xe,
        row=row,
        color=color
    ).visibility_lt(
        measure="width",
        threshold="|xe-x|",
        transitionPadding=10,
        target="mark"
    )

    # Three tracks
    # 1) Text track - uses bed_track
    # 2) Rule track (for tooltips) - uses bed_track
    # 3) Rect track (for exons in a separate BED file)

    text_track = gene_track.mark_text().encode(
        text=gos.Text(field="geneName", type="nominal"),    # type: ignore
        style=gos.Style(dy=-10) # type: ignore
    )

    tooltip_track = gene_track.mark_rule().encode(
        tooltip=[
            gos.Tooltip(field="chromStart", type="genomic", alt="Start Position"),  # type: ignore
            gos.Tooltip(field="chromEnd", type="genomic", alt="End Position"),  # type: ignore
            gos.Tooltip(field="strand", type="nominal", alt="Strand"),  # type: ignore
            gos.Tooltip(field="geneName", type="nominal", alt="Name")  # type: ignore
        ],
        strokeWidth=gos.StrokeWidth(value=1)
    ).properties(
        id=f"{title}-annotation"
    ).visibility_lt(
        measure="width",
        threshold="|xe-x|",
        transitionPadding=10,
        target="mark"
    )

    tracks = [text_track, tooltip_track]

    exon_file_name = ASSEMBLY_TO_EXON_FILE.get(assembly, None)
    if not exon_file_name:
        # non-fatal
        print(f"WARNING: Assembly {assembly} does not have an exon BED file.", file=sys.stderr)
        return tracks

    exon_tbi_file_name = exon_file_name + ".tbi"
    exon_file_url = tracks_db_root + "/" + exon_file_name
    exon_tbi_file_url = tracks_db_root + "/" + exon_tbi_file_name

    exon_data = gos.BedData(
        url=exon_file_url, # type: ignore
        indexUrl=exon_tbi_file_url,  # type: ignore
        type="bed"  # type: ignore
    )

    exon_track = gos.Track(
        data=exon_data  # type: ignore
    ).mark_rect().encode(
        color=color,
        row=row,
        x=x,
        xe=xe,
        size=gos.Size(value=10) # type:ignore
    ).visibility_lt(
        measure="width",
        threshold="|xe-x|",
        transitionPadding=10,
        target="mark"
    )

    tracks.append(exon_track)

    panel_title = "Panel B" if title == "right" else "Panel A"
    if not zoom:
        panel_title = "Annotation"

    annotation_view = gos.overlay(*tracks).properties(
            title=panel_title,
            width=EXPANDED_WIDTH if zoom else CONDENSED_WIDTH,
            height=EXPANDED_HEIGHT,    # Always taller for annotations
            opacity=gos.Opacity(value=0.8)
            )

    return annotation_view

def build_genome_wide_view(assembly, zoom=False, chromosome_only: str | None = None) -> gos.View:
    """
    Build a Gosling genome-wide view track for a given genome assembly.

    This function creates a visualization track that displays all chromosomes for the specified assembly,
    using chromosome size information from a remote TSV file. The view includes a chromosome
    representation and interactive brush marks for zooming and selection.

    Args:
        assembly (str): The genome assembly identifier (e.g., "hg38", "mm10").
        zoom (bool, optional): If True, adds an additional brush for a secondary zoom panel. Defaults to False.

    Returns:
        gos.View: A Gosling view object representing the genome-wide track.
    """

    # NOTE: This may be tempporary or permanent.  Ideally each localhost will have their own set of chromosome sizes files.
    chrom_sizes_url = f"https://umgear.org/tracks/{assembly}.chromInfo.txt"

    data = gos.CsvData(
        url=chrom_sizes_url, # type: ignore
        type="csv", # type: ignore
        headerNames=["chrom", "chromStart", "chromEnd"], # type: ignore
        separator="\t", # type: ignore
        chromosomeField="chrom", # type: ignore
        genomicFields=["chromStart", "chromEnd"] # type: ignore
    )

    base = gos.Track(data)  # type: ignore

    title = f"Chromosome {chromosome_only}" if chromosome_only else f"{assembly} assembly"
    if zoom and chromosome_only:
        title = f"Panel A chromosome {chromosome_only}"

    # build *tracks for gos.overlay
    # The Gos documentation notes it is more idiomatic to set specific properties after initialization
    tracks = [
        base.mark_rect().encode(
            x=gos.X(field="chromStart", type="genomic"), # type:ignore
            xe=gos.X(field="chromEnd", type="genomic"), # type:ignore
            color=gos.Color(field="chrom", type="nominal", range=["#666666", "#999999"]) # type:ignore
        ).properties(
            title=title,
        ),
        base.mark_brush().encode(
            x=gos.X(linkingId="zoom-to-panel-a"), # type:ignore
            color=gos.Color(value="steelblue"),
            stroke=gos.Stroke(value="steelblue"),
            strokeWidth=gos.StrokeWidth(value=3)  # type: ignore
        )
    ]
    if zoom:
        tracks.append(
            base.mark_brush().encode(
                x=gos.X(linkingId="zoom-to-panel-b"), # type:ignore
                color=gos.Color(value="yellow"),
                stroke=gos.Stroke(value="yellow"),
                strokeWidth=gos.StrokeWidth(value=3)  # type: ignore
            )
        )

    genome_wide_view = gos.overlay(*tracks).properties(
        static=True,
        id="chromosome-wide" if chromosome_only else "genome-wide",
        width=CONDENSED_WIDTH,
        height=20,
        data=gos.Data(
            url=chrom_sizes_url,
            type="csv",
            headerNames=["chrom", "chromStart", "chromEnd"],
            separator="\t",
            chromosomeField="chrom",
            genomicFields=["chromStart", "chromEnd"]
        ),
    )

    if chromosome_only:
        # restrict xDomain to the desired interval on chromosome to make the brush span visibly.
        # Fixes to the left-panel chromosome.
        genome_wide_view = genome_wide_view.properties(
            xDomain=gos.GenomicDomain(chromosome=chromosome_only)
        )

    return genome_wide_view

def build_gosling_tracks(parent_tracks_dict, tracks, zoom=False):
    """
    Builds and configures Gosling tracks based on the provided track specifications.

    Args:
        parent_tracks_dict (dict): A dictionary with keys "left" and optionally "right", each mapping to a list of track objects.
        tracks (list): A list of dictionaries, each representing a track specification. Each track dict should contain at least:
            - "type" (str): The type of the track (e.g., "bam", "bigWig", "bed", "vcf").
            - "bigDataUrl" (str): The URL to the data file for the track.
            - Optional: "color" (str), "group" (any), and other track-specific attributes.
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
        "bed": BedSpec,
        "vcf": VcfSpec,
    }

    # Build each individual track based on its type
    for track in tracks:
        track_type = track.get("type", "")
        spec_builder_class = TRACK_TYPE_2_SPEC.get(track_type, None)
        if not spec_builder_class:
            print(f"WARNING: Unsupported track type '{track_type}' for track '{track.get('shortLabel', '')}'; skipping.", file=sys.stderr)
            continue

        data_url = track.get("bigDataUrl", None)
        if not data_url:
            print(f"WARNING: No bigDataUrl found for track '{track.get('shortLabel', '')}'; skipping.", file=sys.stderr)
            continue

        # Get other attributes to pass to the class
        color = track.get("color", "steelblue")  # Default color if not specified
        group = track.get("group", None)

        spec_builder = spec_builder_class(data_url=data_url, color=color, group=group, zoom=zoom)
        left_track = spec_builder.addTrack()
        parent_tracks_dict["left"].append(left_track)

        if zoom:
            right_track = spec_builder.addTrack()
            right_track.id = f"right-track-{Path(data_url).stem}"
            parent_tracks_dict["right"].append(right_track)

    parent_view_left = gos.vertical(*parent_tracks_dict["left"]).properties(
        id="left-view",
        linkingId="zoom-to-panel-a",
        spacing=0
    )

    parent_view_right = None
    if zoom:
        parent_view_right = gos.vertical(*parent_tracks_dict["right"]).properties(
            id="right-view",
            linkingId="zoom-to-panel-b",
            spacing=0
        )


    return parent_view_left, parent_view_right

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

def convert_track_rgb_to_hex(rgb):
    """
    Convert the tracksdb.txt RGB "color" values (i.e. 31,119,180) to a Hex value
    """

    r, g, b = map(int, rgb.split(","))
    return f"#{r:02x}{g:02x}{b:02x}"

def fetch_trackdb_and_groups_info(genomes_txt, assembly) -> dict:
    """Extract 'trackDb' and 'groups' URLs for an assembly from genomes_txt.

    Looks for a "genome <assembly>" line, then reads the immediate following
    "trackDb" and optional "groups" lines. Returns a dict with keys
    "trackDb" and "groups" (empty string if not found).

    NOTE: These can be relative paths to the genomes.txt location; caller
    must resolve them if needed.
    """
    urls = {"trackDb": "", "groups": ""}

    for line in genomes_txt.splitlines():
        if line.startswith(f"genome {assembly}"):
            # The next line should contain the trackDb URL
            next_line_index = genomes_txt.splitlines().index(line) + 1
            if next_line_index < len(genomes_txt.splitlines()):
                next_line = genomes_txt.splitlines()[next_line_index]
                if next_line.startswith("trackDb"):
                    urls["trackDb"] = next_line.split(" ")[1]
                    # Now look for groups line (next line)
                    next_next_line_index = next_line_index + 1
                    if next_next_line_index < len(genomes_txt.splitlines()):
                        next_next_line = genomes_txt.splitlines()[next_next_line_index]
                        if next_next_line.startswith("groups"):
                            urls["groups"] = next_next_line.split(" ")[1]
                    break
    return urls

def parse_groups_from_groupsdb(groupsdb_txt) -> dict:
    """Parse groups info from a UCSC groups.txt content.

    Returns a dict with group names as keys and their descriptions as values.

    Example group:

    name ATAC
    label ATAC-seq
    defaultIsClosed 0
    """
    groups = {}
    current_group = None
    for line in groupsdb_txt.splitlines():
        if line.startswith("name"):
            current_group = line.split(" ")[1]
            groups[current_group] = {}
        elif current_group:
            if line.startswith("label"):
                groups[current_group]["label"] = line.split(" ")[1]
            elif line.startswith("defaultIsClosed"):
                groups[current_group]["defaultIsClosed"] = line.split(" ")[1]
    return groups

def parse_tracks_from_trackdb(trackdb_txt, trackdb_url) -> list:
    """Parse track names from a UCSC trackDb.txt content.

    Returns a list of dicts with track information.

    Example track:
    track P1HC_ATAC_1
    bigDataUrl P1HC_ATAC_1.bigwig
    shortLabel ATAC-seq 1st replicate
    longLabel ATAC-seq 1st replicate
    group ATAC
    color 31,119,180
    autoscale on
    visibility dense
    type bigWig
    """

    tracks = [] # List of dicts
    current_track = {}
    for line in trackdb_txt.splitlines():
        if line.startswith("track"):
            if current_track:
                tracks.append(current_track)
            current_track = {"name": line.split(" ")[1]}
        elif current_track:
            if line.startswith("bigDataUrl"):
                current_track["bigDataUrl"] = line.split(" ")[1]
                # If not a URL, make it one by replacing the "trackDb.txt" part of trackdb_url
                if not current_track["bigDataUrl"].startswith("http://") \
                    and not current_track["bigDataUrl"].startswith("https://"):
                    current_track["bigDataUrl"] = urljoin(trackdb_url, current_track['bigDataUrl'])
            elif line.startswith("shortLabel"):
                current_track["shortLabel"] = " ".join(line.split(" ")[1:])
            elif line.startswith("longLabel"):
                current_track["longLabel"] = " ".join(line.split(" ")[1:])
            elif line.startswith("group"):
                current_track["group"] = line.split(" ")[1]
            elif line.startswith("color"):
                color = line.split(" ")[1]
                current_track["color"] = convert_track_rgb_to_hex(color)
            elif line.startswith("type"):
                current_track["type"] = line.split(" ")[1]
    if current_track:
        tracks.append(current_track)
    return tracks

def zoom_view_to_gene(view, gene_symbol, dataset_id, assembly):
    """
    Zooms a Gosling view to the genomic coordinates of a specified gene.

    This function fetches the genomic coordinates for a given gene symbol from the geardb database,
    adjusts the chromosome naming conventions as needed, and sets the x-domain of the provided Gosling
    view to focus on the region surrounding the gene (with additional padding). It also constructs a
    position string suitable for export to the UCSC Genome Browser.

    Args:
        view: The Gosling view object to be modified.
        gene_symbol (str): The gene symbol to zoom to.
        dataset_id (str or int): The identifier for the dataset containing gene annotations.
        assembly (str): The genome assembly (e.g., "hg38", "mm10").

    Returns:
        tuple:
            - view: The modified Gosling view object with updated x-domain.
            - position_str (str or None): The UCSC Genome Browser position string, or None if gene not found.
            - chrom (str or None): The chromosome name, or None if not available.

    Notes:
        - If the gene symbol is not found or does not have valid coordinates, the original view is returned
          along with None for position_str and/or chrom.
        - Chromosome naming conventions are adjusted to match UCSC standards (e.g., "chr1", "chrX", "chrM").
    """


    BASE_PADDING = 1500  # base padding on each side of gene

    # Fetch gene coordinates from geardb
    gene_info = geardb.get_gene_by_gene_symbol(gene_symbol, dataset_id)
    if not gene_info:
        print(f"WARNING: Gene symbol '{gene_symbol}' not found in dataset {dataset_id}; cannot zoom track.", file=sys.stderr)
        return (view, None, None)

    chrom = gene_info.molecule or "unknown"

    # if chrom is a number, convert it to a string with "chr" prefix
    if chrom.isdigit() or chrom in ["X", "Y"]:
        chrom = f"chr{chrom}"
    # if chrom is MT, convert to "chrM"
    if chrom == "MT":
        chrom = "chrM"

    start = gene_info.start
    end = gene_info.stop

    if not start or not end:
        print(f"WARNING: Gene symbol '{gene_symbol}' does not have valid start/end coordinates; cannot zoom track.", file=sys.stderr)
        return (view, None, chrom)

    left_position = f"{chrom}:{start}-{end}"
    position_str = f"{assembly}.{left_position}"  # This is the format if we want to export position to UCSC Genome Browser

    # Set the x domain of the track to the gene coordinates
    view = view.properties(
        xDomain=gos.GenomicDomain(chromosome=chrom, interval=[start-BASE_PADDING, end+BASE_PADDING])
    )
    return (view, position_str, chrom)

class TrackSpec(ABC):
    def __init__(self, data_url,color="steelblue", group=None, zoom=False):
        self.data_url = data_url
        self.color = color  # Passed as RGB string
        self.group = group
        self.zoom = zoom
        self.width = EXPANDED_WIDTH if zoom else CONDENSED_WIDTH
        self.height = CONDENSED_HEIGHT
        self.title = Path(self.data_url).stem

        self.track = None

    @abstractmethod
    def addTrack(self):
        pass

    @abstractmethod
    def validate_url(self, url):
        pass

class BamSpec(TrackSpec):
    def addTrack(self):

        url = self.data_url
        color = self.color

        try:
            self.validate_url(url)
            self.validate_index_url(f"{url}.bai")
        except ValueError:
            raise

        bamData = gos.BamData(type="bam", url=url, indexUrl=f"{url}.bai") # type: ignore

        track = gos.Track(
            data=bamData, # pyright: ignore[reportArgumentType]
            width=self.width,
            height=self.height,
            title=self.title,  # Use the file name as the title
            id=f"left-track-{self.title}"  # Use the file name without extension as the ID
        ).mark_bar().encode(
            x=gos.X(field="start", type="genomic"), # pyright: ignore[reportArgumentType]
            xe=gos.X(field="end", type="genomic"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="coverage", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def validate_index_url(self, url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .bai
        if not url.endswith(".bai"):
            raise ValueError("Invalid URL: must end with .bai")
        return True

    def validate_url(self,url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .bam
        if not url.endswith(".bam"):
            raise ValueError("Invalid URL: must end with .bam")
        return True

class BedSpec(TrackSpec):
    def addTrack(self):

        url = self.data_url
        color = self.color

        try:
            self.validate_url(url)
            self.validate_index_url(f"{url}.tbi")
        except ValueError:
            raise

        bedData = gos.BedData(type="bed", url=url, indexUrl=f"{url}.tbi") # type: ignore

        track = gos.Track(
            data=bedData, # pyright: ignore[reportArgumentType]
            width=self.width,
            height=self.height,
            title=self.title,  # Use the file name as the title
            id=f"left-track-{self.title}"  # Use the file name without extension as the ID
        ).mark_area().encode(
            x=gos.X(field="start", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
            xe=gos.X(field="end", type="genomic"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="value", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def validate_index_url(self, url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .tbi
        if not url.endswith(".bed.gz.tbi"):
            raise ValueError("Invalid URL: must end with .bed.gz.tbi")
        return True

    def validate_url(self,url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .bed.gz
        if not url.endswith(".bed.gz"):
            raise ValueError("Invalid URL: must end with .bed.gz")
        return True

class BigWigSpec(TrackSpec):
    def addTrack(self):

        url = self.data_url
        color = self.color

        try:
            self.validate_url(url)
        except ValueError:
            raise

        #TODO: Figure out aggregation when zooming out.  Default is mean, but the initial zoomout looks boxy until you zoom out further.
        bigWigData = gos.BigWigData(type="bigwig", url=url) # type: ignore

        track = gos.Track(
            data=bigWigData, # pyright: ignore[reportArgumentType]
            width=self.width,
            height=self.height,
            title=self.title,  # Use the file name as the title
            id=f"left-track-{self.title}"  # Use the file name without extension as the ID
        ).mark_area().encode(
            x=gos.X(field="start", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
            xe=gos.X(field="end", type="genomic"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="value", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def validate_url(self, url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .bw or .bigwig
        if not (url.endswith(".bw") or url.endswith(".bigwig")):
            raise ValueError("Invalid URL: must end with .bw or .bigwig")
        return True

class VcfSpec(TrackSpec):
    def addTrack(self):

        url = self.data_url
        color = self.color

        try:
            self.validate_url(url)
            self.validate_index_url(f"{url}.tbi")
        except ValueError:
            raise

        vcfData = gos.VcfData(type="vcf", url=url, indexUrl=f"{url}.tbi") # type: ignore

        track = gos.Track(
            data=vcfData, # pyright: ignore[reportArgumentType]
            width=self.width,
            height=self.height,
            title=self.title,  # Use the file name as the title
            id=f"left-track-{self.title}"  # Use the file name without extension as the ID
        ).mark_point().encode(
            x=gos.X(field="position", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="value", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def validate_index_url(self, url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .tbi
        if not url.endswith(".vcf.gz.tbi"):
            raise ValueError("Invalid URL: must end with .vcf.gz.tbi")
        return True

    def validate_url(self, url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .vcf.gz
        if not url.endswith(".vcf.gz"):
            raise ValueError("Invalid URL: must end with .vcf.gz")
        return True

# Assembly track class
class AssemblySpec:
    pass

#####################

class GoslingSpec(Resource):
    def get(self, dataset_id):
        #session_id = request.cookies.get("gear_session_id", "")
        args = request.args
        gene_symbol = args.get("gene")
        assembly = args.get("assembly")
        hub_url = args.get("hub_url", "")
        zoom = args.get("zoom", "false")
        # Set zoom from string to bool
        if zoom.lower() == 'true':
            zoom = True
        else:
            zoom = False

        response = {
            "success": 0,
            "spec": {},
            "position": "", # will be filled with chr:start-end
            "message": ""
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

        # Cut off name of hub_url (hub.txt). This will be used to build more paths
        base_url = hub_url.rsplit('/', 1)[0]  # Get base URL of hub.txt

        # Look for a genomes.txt file to matching the assembly to get the "trackDb" file
        # This file is the same as the one required by the UCSC Genome Browser
        genomes_url = f"{base_url}/genomes.txt"
        try:
            genomes_response = requests.get(genomes_url)
            genomes_response.raise_for_status()
        except requests.RequestException as e:
            response["message"] = f"Failed to retrieve genomes.txt from {genomes_url}: {str(e)}"
            return response, 400

        urls = fetch_trackdb_and_groups_info(genomes_response.text, assembly)
        if not urls:
            response["message"] = f"trackDb URL not found for assembly {assembly} in genomes.txt"
            return response, 400

        # If trackDb and groups URLs are relative to genomes.txt, resolve them
        trackdb_url = urls["trackDb"]
        groups_url = urls["groups"]

        if not trackdb_url.startswith("http://") and not trackdb_url.startswith("https://"):
            trackdb_url = f"{base_url}/{trackdb_url}"
        if groups_url and not groups_url.startswith("http://") and not groups_url.startswith("https://"):
            groups_url = f"{base_url}/{groups_url}"

        # Fetch tracks
        try:
            trackdb_response = requests.get(trackdb_url)
            trackdb_response.raise_for_status()
        except requests.RequestException as e:
            response["message"] = f"Failed to retrieve trackDb from {trackdb_url}: {str(e)}"
            return response, 400
        tracks = parse_tracks_from_trackdb(trackdb_response.text, trackdb_url)

        # Attempt to fetch groups info (non-fatal)
        groups = {}
        if groups_url:
            try:
                groups_response = requests.get(groups_url)
                groups_response.raise_for_status()
                # Parse groups info as needed; here we just store the raw text
                groups = parse_groups_from_groupsdb(groups_response.text)
            except requests.RequestException:
                print("INFO: Could not retrieve groups info; continuing without it.", file=sys.stderr)
                pass  # Ignore errors in fetching groups
        else:
            print("INFO: No groups URL provided; continuing without it.", file=sys.stderr)

        gos_tracks = {"left": [], "right": []}
        # Add BED annotation tracks to left (and right if zoom) gos_tracks in the first index position
        # build left and right track
        # Insert into gos_tracks["left"] at index 0 (and gos_tracks["right"] if zoom)
        gos_tracks["left"].append(build_bed_annotation_tracks(assembly, zoom, "left"))
        if zoom:
            gos_tracks["right"].append(build_bed_annotation_tracks(assembly, zoom, "right"))

        # If groups are present, gather all tracks for this group,
        # and attempt to aggregate the data by mean for each position point.
        # This is an effort to cut down on the number of tracks to render.
        #
        # Otherwise just process all tracks individually
        if groups:
            for group in groups:
                # TODO: Aggregate by group.  May require extra packages. Can cache (maybe upon upload even).
                # Get tracks associated with this group
                group_tracks = [track for track in tracks if track.get("group", "") == group]
                if not group_tracks:
                    continue
                (parent_view_left, parent_view_right) = build_gosling_tracks(gos_tracks, group_tracks, zoom=zoom)
        else:
            (parent_view_left, parent_view_right) = build_gosling_tracks(gos_tracks, tracks, zoom=zoom)

        # Start building the Gosling spec
        region_view = build_region_view(parent_view_left, parent_view_right)
        # At this point, let's zoom this view to the coordinates of the gene_symbol.
        (region_view, position_str, chrom) = zoom_view_to_gene(region_view, gene_symbol, dataset_id, assembly)

        base_views = []

        genome_wide_view = build_genome_wide_view(assembly, zoom)
        base_views.append(genome_wide_view)
        if chrom:
            chromosome_view = build_genome_wide_view(assembly, zoom=zoom, chromosome_only=chrom)
            base_views.append(chromosome_view)

        base_views.append(region_view)
        base_track = gos.vertical(*base_views)

        # Add assembly track to base track
        assembly_obj = build_assembly_gos_obj(assembly)
        base_track = base_track.properties(
            assembly=assembly_obj,
        )

        spec = base_track.to_json(indent=2)

        response["success"] = 1
        response["spec"] = spec
        response["position"] = position_str
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