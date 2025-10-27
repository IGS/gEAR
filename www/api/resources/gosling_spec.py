import json
import sys
from abc import ABC, abstractmethod
from pathlib import Path

import geardb
import gosling as gos
import requests
from flask import request
from flask_restful import Resource

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

# DUMMY VARIABLES

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

    chromosome_sizes_name = ASSEMBLY_TO_CHROMSIZES_FILE.get(assembly, None)
    if not chromosome_sizes_name:
        raise ValueError(f"Assembly {assembly} is not supported or does not have a chromosome sizes file.")

    tracksDbRoot = "https://umgear.org/tracks"  # Replace with actual tracks database
    chromosome_sizes_url = tracksDbRoot + "/" + chromosome_sizes_name
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

def build_gosling_tracks(parent_tracks_dict,tracks, zoom=False):
    """
    Builds Gosling track specifications for a list of track configurations.
    This function iterates over a list of track dictionaries, determines the appropriate
    specification builder class based on the track type, and constructs Gosling track
    specifications. Tracks are organized into 'left' and optionally 'right' groups,
    depending on the `zoom` parameter.
    Args:
        tracks (list): A list of dictionaries, each representing a track configuration.
            Each dictionary should contain at least the keys 'type' and 'bigDataUrl'.
        zoom (bool, optional): If True, builds additional tracks for the 'right' group
            to support zoomed-in views. Defaults to False.
    Returns:
        dict: A dictionary with keys 'left' and 'right', each containing a list of
            constructed Gosling track specifications. The 'right' list is empty if
            `zoom` is False.
    Notes:
        - Supported track types are 'bam', 'bigWig', 'bed', and 'vcf'.
        - Tracks with unsupported types or missing 'bigDataUrl' are skipped with a warning.
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
            print(f"WARNING: Unsupported track type '{track_type}' for track '{track.get('name', '')}'; skipping.", file=sys.stderr)
            continue

        data_url = track.get("bigDataUrl", None)
        if not data_url:
            print(f"WARNING: No bigDataUrl found for track '{track.get('name', '')}'; skipping.", file=sys.stderr)
            continue

        # Get other attributes to pass to the class
        color = track.get("color", "steelblue")  # Default color if not specified
        group = track.get("group", None)

        spec_builder = spec_builder_class(data_url=data_url, color=color, group=group, zoom=zoom)
        track = spec_builder.addTrack()
        parent_tracks_dict["left"].append(track.track)

        if zoom:
            right_track = spec_builder.addTrack()
            right_track.id = f"right-track-{Path(data_url).stem}"
            parent_tracks_dict["right"].append(right_track)
    return parent_tracks_dict

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

def parse_tracks_from_trackdb(trackdb_txt) -> list:
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
            elif line.startswith("shortLabel"):
                current_track["shortLabel"] = " ".join(line.split(" ")[1:])
            elif line.startswith("longLabel"):
                current_track["longLabel"] = " ".join(line.split(" ")[1:])
            elif line.startswith("group"):
                current_track["group"] = line.split(" ")[1]
            elif line.startswith("color"):
                current_track["color"] = line.split(" ")[1]
            elif line.startswith("type"):
                current_track["type"] = line.split(" ")[1]
    if current_track:
        tracks.append(current_track)
    return tracks

def zoom_track_to_gene(track, gene_symbol, dataset_id, assembly):
    """
    Zooms a genomic track to the coordinates of a specified gene, with optional padding, and returns the updated track and a UCSC Genome Browser position string.

    Args:
        track: The genomic track object to be zoomed.
        gene_symbol (str): The gene symbol to zoom to.
        dataset_id (str or int): The identifier for the dataset containing gene annotations.
        assembly (str): The genome assembly (e.g., 'hg38', 'mm10').

    Returns:
        tuple: A tuple containing the updated track object and a position string formatted for the UCSC Genome Browser.
               If the gene is not found or coordinates are invalid, returns the original track and does not set the position string.

    Notes:
        - Adds a base padding of 1500 bases on each side of the gene.
        - Handles chromosome naming conventions (e.g., adds 'chr' prefix, converts 'MT' to 'chrM').
        - Prints warnings to stderr if the gene is not found or coordinates are invalid.
    """

    BASE_PADDING = 1500  # base padding on each side of gene

    # Fetch gene coordinates from geardb
    gene_info = geardb.get_gene_by_gene_symbol(gene_symbol, dataset_id)
    if not gene_info:
        print(f"WARNING: Gene symbol '{gene_symbol}' not found in dataset {dataset_id}; cannot zoom track.", file=sys.stderr)
        return track

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
        return track

    left_position = f"{chrom}:{start}-{end}"
    position_str = f"{assembly}.{left_position}"  # This is the format if we want to export position to UCSC Genome Browser

    # Set the x domain of the track to the gene coordinates
    track.properties(
        xDomain=gos.GenomicDomain(chromosome=chrom, interval=[start-BASE_PADDING, end+BASE_PADDING])
    )
    return (track, position_str)

class TrackSpec(ABC):
    def __init__(self, data_url,color="steelblue", group=None, zoom=False):
        self.data_url = data_url
        self.color = color  # Passed as RGB string
        self.group = group
        self.zoom = zoom

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

        bamData = gos.BamData(url=url, indexUrl=f"{url}.bai") # type: ignore
        bamData.color = color

        track = gos.Track(
            data=bamData, # pyright: ignore[reportArgumentType]
            width=EXPANDED_WIDTH if self.zoom else CONDENSED_WIDTH,
            height=EXPANDED_HEIGHT if self.zoom else CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
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

        bedData = gos.BedData(url=url, indexUrl=f"{url}.tbi") # type: ignore
        bedData.color = color

        track = gos.Track(
            data=bedData, # pyright: ignore[reportArgumentType]
            width=EXPANDED_WIDTH if self.zoom else CONDENSED_WIDTH,
            height=EXPANDED_HEIGHT if self.zoom else CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
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
        bigWigData = gos.BigWigData(url=url)
        bigWigData.color = color

        track = gos.Track(
            data=bigWigData, # pyright: ignore[reportArgumentType]
            width=EXPANDED_WIDTH if self.zoom else CONDENSED_WIDTH,
            height=EXPANDED_HEIGHT if self.zoom else CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
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

        vcfData = gos.VcfData(url=url, indexUrl=f"{url}.tbi") # type: ignore
        vcfData.color = color

        track = gos.Track(
            data=vcfData, # pyright: ignore[reportArgumentType]
            width=EXPANDED_WIDTH if self.zoom else CONDENSED_WIDTH,
            height=EXPANDED_HEIGHT if self.zoom else CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
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
        spec_dir_url = args.get("spec_dir_url", "") # may not be needed. Still in evaluation mode
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
        tracks = parse_tracks_from_trackdb(trackdb_response.text)

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


        # If groups are present, gather all tracks for this group,
        # and attempt to aggregate the data by mean for each position point.
        # This is an effort to cut down on the number of tracks to render.
        #
        # Otherwise just process all tracks individually
        if groups:
            for group in groups:
                # Get tracks associated with this group
                group_tracks = [track for track in tracks if track.get("group", "") == group]
                if not group_tracks:
                    continue
                gos_tracks = build_gosling_tracks(gos_tracks, group_tracks, zoom=zoom)
        else:
            gos_tracks = build_gosling_tracks(gos_tracks, tracks, zoom=zoom)

        # At this point, let's zoom the left track to the coordinates of the gene_symbol.
        parent_view_left = gos.View(tracks=gos_tracks["left"])
        parent_view_right = gos.View(tracks=gos_tracks["right"]) if zoom else None

        (parent_view_left, position_str) = zoom_track_to_gene(parent_view_left, gene_symbol, dataset_id, assembly)
        # If zoom is True, also zoom the right track
        if zoom:
            (parent_view_right, _) = zoom_track_to_gene(parent_view_right, gene_symbol, dataset_id, assembly)

        # Start building the Gosling spec
        genome_wide_view = gos.View()
        region_view = gos.View()
        base_track = gos.vertical(genome_wide_view, region_view)

        # Add assembly track to base track
        assembly_obj = build_assembly_gos_obj(assembly)
        base_track.properties(assembly=assembly_obj)




        spec = base_track.to_json(indent=2)

        response["success"] = 1
        response["spec"] = spec
        response["position"] = position_str
        return response, 200

        # Implement your logic to retrieve the Gosling spec for the given dataset_id
        # return view.to_json(indent)

        ### Monkeypatched code ###

        spec_path = JSON_DIR.joinpath("Litao_bigwigs/condensed.json")
        if zoom:
            spec_path = JSON_DIR.joinpath("Litao_bigwigs/expanded.json")

        # return the JSON
        if spec_path.exists():
            with open(spec_path, 'r') as file:
                spec = json.load(file)
            response["success"] = 1
            response["spec"] = spec
            return response, 200
        else:
            return {"message": "Spec not found"}, 404


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