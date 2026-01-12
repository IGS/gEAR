from pathlib import Path
import shlex
import subprocess
import sys

import requests
from Bio import bgzf
from flask import request
from flask_restful import Resource

gear_root = Path(__file__).resolve().parents[3]  # web-root dir
src_path = gear_root / "src"

VALID_TYPES = ["bigWig", "bigBed"]
BIGBED_EXTENSIONS = [".bb", ".bigbed"]

user_upload_file_base = gear_root / 'www' / 'uploads' / 'files'

def bigbed_to_bed(bigbed_path: Path, outdir_path: Path) -> bool:
    """
    Convert a bigBed file to a bed file using the UCSC tool bigBedToBed.
    Next, bgzip the file and tabix the bgzipped file.
    The bed file will be saved in the output_dir
    """

    bigbed_file = bigbed_path.as_posix()
    bed_path = bigbed_path.with_suffix('.bed')
    bed_path = outdir_path / bed_path.name
    bed_file = bed_path.as_posix()

    exec_file = src_path / "bigBedToBed"

    try:
        subprocess.run([exec_file, bigbed_file, bed_file], check=True)
        print(f"Converted {bigbed_file} to {bed_file}.", file=sys.stderr)

        gz_path = bed_path.with_suffix(bed_path.suffix + '.gz')
        with bed_path.open('rb') as f_in, bgzf.BgzfWriter(gz_path) as f_out:
            # write the entirety of the input
            f_out.write(f_in.read())
        create_tabix_indexed_file(gz_path, file_type="bed")

        return True
    except subprocess.CalledProcessError as e:
        print(f"Error converting {bigbed_file} to bed: {e}", file=sys.stderr)
        return False

def create_tabix_indexed_file(bgzipped_path: Path, file_type : str="bed") -> None:
    """
    Tabix-indexes a genomic file.
    bgzipped_file: Path to the bgzipped input file.
    file_type: Type of the file (e.g., "bed", "vcf", "gff").
    """

    tabix_cmd = f"tabix -s 1 -b 2 -e 3 -p {file_type} {bgzipped_path.as_posix()}"
    subprocess.run(shlex.split(tabix_cmd), check=True)


def fetch_trackdb_and_groups_info(genomes_txt, assembly) -> dict:
    """Extract 'trackDb' and 'groups' URLs for an assembly from genomes_txt.

    Looks for a "genome <assembly>" line, then reads the immediate following
    "trackDb" and optional "groups" lines. Returns a dict with keys
    "trackDb" and "groups" (empty string if not found).

    NOTE: These can be relative paths to the genomes.txt location; caller
    must resolve them if needed.
    """
    urls = {"trackDb": "", "groups": ""}

    for genome_line in genomes_txt.splitlines():
        if genome_line.startswith(f"genome {assembly}"):
            # The next line should contain the trackDb URL
            next_line_index = genomes_txt.splitlines().index(genome_line) + 1
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
        raise ValueError(f"Assembly {assembly} not found in genomes.txt")
    return urls

def validate_track_types_from_db(trackdb_txt, trackdb_url) -> list:
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

    invalid_tracks = []

    for track_line in trackdb_txt.splitlines():
        if track_line.startswith("type"):
            tracktype = track_line.split(" ")[1]
            if tracktype not in VALID_TYPES:
                invalid_tracks.append({"type": tracktype, "trackdb_url": trackdb_url})
    return invalid_tracks

class TrackHubValidate(Resource):
    def post(self, share_uid):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        hub_url = req.get("trackhub_url")
        assembly = req.get("assembly")

        result = {
            "success": False,
            "message": "",
            "num_tracks": 0,
        }

        if not hub_url or not hub_url.startswith("http"):
            result["message"] = "Invalid URL"
            return result, 400

        # use the "hubCheck" utiliy to validae the passed in hub file
        hubcheck_exe = src_path / "hubCheck"
        try:
            completed_process = subprocess.run(
                shlex.split(f"{hubcheck_exe} {hub_url}"),
                check=True
            )
        except subprocess.CalledProcessError as e:
            result["message"] = f"hubCheck failed: {str(e)}"
            return result, 500

        # Cut off name of hub_url (hub.txt). This will be used to build more paths
        base_url = hub_url.rsplit("/", 1)[0]  # Get base URL of hub.txt

        # Look for a genomes.txt file to matching the assembly to get the "trackDb" file
        # This file is the same as the one required by the UCSC Genome Browser
        genomes_url = f"{base_url}/genomes.txt"
        try:
            genomes_response = requests.get(genomes_url)
            genomes_response.raise_for_status()
        except requests.RequestException as e:
            result["message"] = f"Error fetching genomes.txt: {str(e)}"
            return result, 500

        urls = fetch_trackdb_and_groups_info(genomes_response.text, assembly)
        if not urls:
            result["message"] = f"No trackDb found for assembly {assembly}"
            return result, 400
        trackdb_url = urls["trackDb"]
        if not trackdb_url.startswith("http://") and not trackdb_url.startswith("https://"):
            trackdb_url = f"{base_url}/{trackdb_url}"

        # Fetch tracks
        try:
            trackdb_response = requests.get(trackdb_url)
            trackdb_response.raise_for_status()
        except requests.RequestException as e:
            result["message"] = f"Error fetching trackDb: {str(e)}"
            return result, 500
        result["num_tracks"] = trackdb_response.text.count("track ")
        invalid_tracks = validate_track_types_from_db(trackdb_response.text, trackdb_url)
        if len(invalid_tracks):
            result["message"] = f"Invalid track types found. Currently accepted types are [{', '.join(VALID_TYPES)}]."
            return result, 400

        # All good
        result["success"] = True
        result["message"] = "Track hub is valid"
        return result, 200

class TrackHubCopy(Resource):
    def post(self, share_uid):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        hub_url = req.get("trackhub_url")
        assembly = req.get("assembly")

        result = {
            "success": False,
            "message": ""
        }

        if not session_id:
            result["message"] = "Missing session_id"
            return result, 400

        if not share_uid:
            result["message"] = "Missing upload ID"
            return result, 400

        if not hub_url or not hub_url.startswith("http"):
            result["message"] = "Invalid URL"
            return result, 400

        track_upload_dir: Path = user_upload_file_base / session_id / share_uid
        track_upload_dir.mkdir(parents=True, exist_ok=True)

        # use the "hubClone" utility to clone the passed in hub file to our upload directory
        hubclone_exe = src_path / "hubClone"
        try:
            completed_process = subprocess.run(
                shlex.split(f"{hubclone_exe} {hub_url} -download -udcDir={track_upload_dir}"),
                check=True
            )
        except subprocess.CalledProcessError as e:
            result["message"] = f"hubCheck failed: {str(e)}"
            return result, 500

        # Test if the hub.txt file was downloaded
        hub_txt_file = track_upload_dir / "hub.txt"
        if not hub_txt_file.is_file():
            result["message"] = "hubClone failed to download hub.txt"
            return result, 500

        # Find the tracks and convert any BigBed to a tabixed and bgzf compressed bed file
        assembly_dir: Path = track_upload_dir / assembly
        if not assembly_dir.is_dir():
            result["message"] = f"hubClone failed to download assembly directory for {assembly}"
            return result, 500

        # find all bigBed files in the assembly directory
        bigbed_paths = []
        for ext in BIGBED_EXTENSIONS:
            bigbed_paths.extend(assembly_dir.rglob(f"*{ext}"))
        for bigbed_path in bigbed_paths:
            conversion_success = bigbed_to_bed(bigbed_path, assembly_dir)
            if not conversion_success:
                result["message"] = f"Failed to convert {bigbed_path.as_posix()} to bed"
                return result, 500

        return result, 200

class TrackHubStatus(Resource):
    def post(self, share_uid):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        url = req.get("url")
        # ...status logic...
        return {"success": True, "message": "Status update"}