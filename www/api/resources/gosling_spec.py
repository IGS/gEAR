import json
from pathlib import Path

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
LEFT_TRACK_WIDTH = 600 - VIEW_PADDING/2  # Width of the left view tracks
RIGHT_TRACK_WIDTH = 600 - VIEW_PADDING/2  # Width of the right view tracks

CONDENSED_HEIGHT = 20  # Height for condensed tracks
EXPANDED_HEIGHT = 40  # Height for expanded tracks

TRACKS_URL = "https://umgear.org/tracks/"

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

# DUMMY VARIABLES
tracksDbRoot = "https://umgear.org/tracks"  # Replace with actual tracks database

def build_assembly_array(assembly):
    # Take the chromosome sizes and create a 2D array of chromosome names and sizes (e.g., [['chr1', 248956422], ['chr2', 242193529], ...])
    chromosome_sizes_url = ASSEMBLY_TO_CHROMSIZES_FILE.get(assembly, None)
    if not chromosome_sizes_url:
        raise ValueError(f"Assembly {assembly} is not supported or does not have a chromosome sizes file.")

    assembly_array = []

    # Get the data from the chromosome sizes url
    response = requests.get(chromosome_sizes_url)
    if response.status_code != 200:
        raise ValueError(f"Failed to retrieve chromosome sizes from {chromosome_sizes_url}")

    for line in response.text.splitlines():
        chromosome, _, size = line.strip().split("\t")  # my chrom sizes file has 0 in the middle column
        assembly_array.append([chromosome, int(size)])
    return assembly_array


def get_db_display_config():
    # This function should return the display configuration for the database
    # For now, we will return a dummy configuration
    config = {}
    config["gosling_json_path"] = "Litao_bigwigs.json"  # Will add the actual path later
    config["ucsc_hub_url"] = "https://umgear.org/tracks/Litao_bigwigs/hub.txt"
    config["assembly"] = "mm10"
    return config

class BamSpec:
    def addTrack(self, url, color):

        try:
            self.validate_url(url)
            self.validate_index_url(f"{url}.bai")
        except ValueError:
            raise

        bamData = gos.BamData(url=url, indexUrl=f"{url}.bai")
        bamData.color = color

        track = gos.Track(
            data=bamData, # pyright: ignore[reportArgumentType]
            width=LEFT_TRACK_WIDTH,
            height=CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
        ).mark_bar().encode(
            x=gos.X(field="start", type="genomic"), # pyright: ignore[reportArgumentType]
            xe=gos.X(field="end", type="genomic"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="coverage", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def addTracks(self, tracks, color, zoom=False, viewIndex = 0):
        if not isinstance(tracks, list) or len(tracks) == 0:
            print("No tracks provided or tracks is not an array.")
            return

        for track in tracks:
            newTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
            #this.trackViews[viewIndex].tracks.push(newTrack)
            if zoom:
                newTrack["height"] = EXPANDED_HEIGHT

                # copy the object for the right view
                rightTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
                rightTrack["id"] = f"right-track-{track}"
                rightTrack["width"] = RIGHT_TRACK_WIDTH
                rightTrack["height"] = EXPANDED_HEIGHT
                #this.trackViews[viewIndex + 1].tracks.push(rightTrack)

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

class BedSpec:
    def addTrack(self, url, color):

        try:
            self.validate_url(url)
            self.validate_index_url(f"{url}.tbi")
        except ValueError:
            raise

        bedData = gos.BedData(url=url, indexUrl=f"{url}.tbi")
        bedData.color = color

        track = gos.Track(
            data=bedData, # pyright: ignore[reportArgumentType]
            width=LEFT_TRACK_WIDTH,
            height=CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
        ).mark_area().encode(
            x=gos.X(field="start", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
            xe=gos.X(field="end", type="genomic"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="value", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def addTracks(self, tracks, color, zoom=False, viewIndex = 0):
        if not isinstance(tracks, list) or len(tracks) == 0:
            print("No tracks provided or tracks is not an array.")
            return

        for track in tracks:
            newTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
            #this.trackViews[viewIndex].tracks.push(newTrack)
            if zoom:
                newTrack["height"] = EXPANDED_HEIGHT

                # copy the object for the right view
                rightTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
                rightTrack["id"] = f"right-track-{track}"
                rightTrack["width"] = RIGHT_TRACK_WIDTH
                rightTrack["height"] = EXPANDED_HEIGHT
                #this.trackViews[viewIndex + 1].tracks.push(rightTrack)

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


class BigWigSpec:
    def addTrack(self, url, color):

        try:
            self.validate_url(url)
        except ValueError:
            raise

        #TODO: Figure out aggregation when zooming out.  Default is mean, but the initial zoomout looks boxy until you zoom out further.
        bigWigData = gos.BigWigData(url=url)
        bigWigData.color = color

        track = gos.Track(
            data=bigWigData, # pyright: ignore[reportArgumentType]
            width=LEFT_TRACK_WIDTH,
            height=CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
        ).mark_area().encode(
            x=gos.X(field="start", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
            xe=gos.X(field="end", type="genomic"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="value", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def addTracks(self, tracks, color, zoom=False, viewIndex = 0):
        # loop through each of the bigWig files and add their spec to the outer track for both the left and right views
        if not isinstance(tracks, list) or len(tracks) == 0:
            print("No tracks provided or tracks is not an array.")
            return

        for track in tracks:
            newTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
            #this.trackViews[viewIndex].tracks.push(newTrack)
            if zoom:
                newTrack["height"] = EXPANDED_HEIGHT

                # copy the object for the right view
                rightTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
                rightTrack["width"] = RIGHT_TRACK_WIDTH
                rightTrack["height"] = EXPANDED_HEIGHT
                #this.trackViews[viewIndex + 1].tracks.push(rightTrack)

    def validate_url(self, url):
        # Basic URL validation
        if not url.startswith("http://") and not url.startswith("https://"):
            raise ValueError("Invalid URL: must start with http:// or https://")
        # URL must end with .bw or .bigwig
        if not (url.endswith(".bw") or url.endswith(".bigwig")):
            raise ValueError("Invalid URL: must end with .bw or .bigwig")
        return True

class VcfSpec:
    def addTrack(self, url, color):

        try:
            self.validate_url(url)
            self.validate_index_url(f"{url}.tbi")
        except ValueError:
            raise

        vcfData = gos.VcfData(url=url, indexUrl=f"{url}.tbi")
        vcfData.color = color

        track = gos.Track(
            data=vcfData, # pyright: ignore[reportArgumentType]
            width=LEFT_TRACK_WIDTH,
            height=CONDENSED_HEIGHT,
            title=Path(url).name,  # Use the file name as the title
            id=f"left-track-{Path(url).stem}"  # Use the file name without extension as the ID
        ).mark_point().encode(
            x=gos.X(field="position", type="genomic", axis="none"), # pyright: ignore[reportArgumentType]
            y=gos.Y(field="value", type="quantitative", axis="right"), # pyright: ignore[reportArgumentType]
            color=gos.Color(value=color)
        )
        return track

    def addTracks(self, tracks, color, zoom=False, viewIndex = 0):
        if not isinstance(tracks, list) or len(tracks) == 0:
            print("No tracks provided or tracks is not an array.")
            return

        for track in tracks:
            newTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
            #this.trackViews[viewIndex].tracks.push(newTrack)
            if zoom:
                newTrack["height"] = EXPANDED_HEIGHT

                # copy the object for the right view
                rightTrack = self.addTrack(f"{tracksDbRoot}/{track}", color)
                rightTrack["id"] = f"right-track-{track}"
                rightTrack["width"] = RIGHT_TRACK_WIDTH
                rightTrack["height"] = EXPANDED_HEIGHT
                #this.trackViews[viewIndex + 1].tracks.push(rightTrack)

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

#####################

class GoslingSpec(Resource):
    def get(self, dataset_id):
        #session_id = request.cookies.get("gear_session_id", "")
        req = request.args
        #args = parser.parse_args()
        #gene_symbol = args.get("gene")
        #genome = args.get("genome")
        zoom = req.get("zoom", "false")
        # Set zoom from string to bool
        if zoom.lower() == 'true':
            zoom = True
        else:
            zoom = False

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
            return spec, 200
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