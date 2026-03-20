import configparser
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from uuid import uuid4

import pyBigWig
import requests
from Bio import bgzf
from flask import request
from flask_restful import Resource
from hic2cool import hic2cool_convert

gear_root = Path(__file__).resolve().parents[3]  # web-root dir
src_path = gear_root / "src"

VALID_TYPES = ["bigWig", "bigBed", "hic", "vcfTabix"]
VALID_CONTAINER_TYPES = ["multiWig"]

user_upload_file_base = gear_root / 'www' / 'uploads' / 'files'

def append_higlass_url_to_trackdb(trackdb_file: Path, hic_track_name: str, higlass_url: str) -> Path:
    """
    Append a Higlass URL to a TrackDB file by injecting a gos_higlass_url property.

    This function reads a TrackDB file, finds lines where the "bigDataUrl" property
    ends with the specified HiC track name, and inserts a new "gos_higlass_url" property
    with the provided Higlass URL. The modified content is written to a temporary file,
    which then replaces the original file.

    Args:
        trackdb_file (Path): Path to the TrackDB file to be modified.
        hic_track_name (str): The name of the HiC track to match at the end of bigDataUrl lines.
        higlass_url (str): The Higlass URL to be added as the gos_higlass_url property value.

    Returns:
        Path: The path to the modified TrackDB file (same as input trackdb_file).

    Raises:
        FileNotFoundError: If the trackdb_file does not exist.
        IOError: If the file cannot be read or written.
    """

    with open(new_trackdb_file := trackdb_file.with_suffix('.modified.trackDb.txt'), 'w') as out_f:

        for line in str(trackdb_file).splitlines():
            line = line.strip()
            # if the "bigDataUrl" property ends with the hic track name, add to the higlass URL property to the line under it
            if line.startswith("bigDataUrl") and line.endswith(hic_track_name):
                # add the gos_higlass_url property after the type property
                out_f.write(f"gos_higlass_url {higlass_url}\n")
            else:
                out_f.write(line + "\n")
    # Move new file to replace old file
    new_trackdb_file.rename(trackdb_file)

    return trackdb_file

def bigbed_to_bed(bigbed_path: Path, outdir_path: Path) -> bool:
    """
    Converts a BigBed file to a BED file and compresses it using BGZF.

    This function reads a BigBed file, extracts its intervals, and writes them to a BED file.
    The resulting BED file is then compressed using BGZF, and a Tabix index is created for it.

    Args:
        bigbed_path (Path): The path to the input BigBed file.
        outdir_path (Path): The directory where the output BED file and compressed files will be saved.

    Returns:
        bool: True if the conversion, compression, and indexing are successful; False otherwise.

    Side Effects:
        - Writes the converted BED file to the specified output directory.
        - Compresses the BED file using BGZF and creates a Tabix index.
        - Prints status messages and errors to stderr.

    Raises:
        Exception: If an error occurs during file conversion, compression, or indexing.
    """
    bigbed_file = bigbed_path.as_posix()
    bed_path = bigbed_path.with_suffix('.bed')
    bed_path = outdir_path / bed_path.name
    bed_file = bed_path.as_posix()

    # Open the bigBed file for reading
    bb = pyBigWig.open(bigbed_file)

    if not bb.isBigBed():
        print(f"{bigbed_file} is not a bigBed file.")
        bb.close()
        return False

    try:
        # Open the output BED file for writing
        with open(bed_file, 'w') as bed_out:
            # Iterate over all intervals in the bigBed file
            # get_intervals returns a list of tuples (chrom, start, end, rest_of_line)
            for chrom, start, end, rest in bb.intervals():
                # The 'rest' part may contain additional BED columns as a tab-separated string
                bed_line = f"{chrom}\t{start}\t{end}\t{rest}\n"
                bed_out.write(bed_line)

        # Close the bigBed file
        bb.close()

        print(f"Converted {bigbed_file} to {bed_file}.", file=sys.stderr)
    except Exception as e:
        print(f"Error converting {bigbed_file} to bed.: {e}", file=sys.stderr)
        return False

    try:
        gz_path = bed_path.with_suffix(bed_path.suffix + '.gz')
        with bed_path.open('rb') as f_in, bgzf.BgzfWriter(gz_path) as f_out:
            # write the entirety of the input
            f_out.write(f_in.read())
        create_tabix_indexed_file(gz_path, file_type="bed")

        return True
    except Exception as e:
        print(f"Error compressing and indexing {bed_file}: {e}", file=sys.stderr)
        return False

def download_large_file(url, destination):
    """
    Downloads a large file from a URL using streaming to save memory.
    """
    try:
        # Set stream=True to download content in chunks
        with requests.get(url, stream=True) as response:
            response.raise_for_status() # Check if the download was successful
            with open(destination, 'wb') as file:
                # Iterate over the response in chunks and write to the file
                shutil.copyfileobj(response.raw, file)
        print(f"File downloaded successfully to {destination}!")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")


def hic_to_mcool(hic_path: Path, outdir_path: Path) -> bool:
    """
    Convert a .hic file to a .mcool file using the hic2cool tool.
    The mcool file will be saved in the output_dir
    """

    hic_file = hic_path.as_posix()
    mcool_path = hic_path.with_suffix('.mcool')
    mcool_path = outdir_path / mcool_path.name
    mcool_file = mcool_path.as_posix()

    # command to convert .hic to .mcool is: hic2cool_convert(input.hic, output.mcool, 0)
    # The '0' argument is to use all resolutions when building the output

    try:
        hic2cool_convert(hic_file, mcool_file, 0)
        print(f"Converted {hic_file} to {mcool_file}.", file=sys.stderr)
        return True
    except Exception as e:
        print(f"Error converting {hic_file} to mcool: {e}", file=sys.stderr)
        return False

def ingest_mcool_into_higlass(mcool_path: Path) -> str:
    """
    Ingest a .mcool file into HiGlass using the higlass API. Return the URL of the ingested dataset in HiGlass.
    """

    mcool_file = mcool_path.as_posix()

    config = configparser.ConfigParser()
    config.read('../../../gear.ini')

    user = config["higlass"]["admin_user"]
    pw = config["higlass"]["admin_pass"]
    hostname = config["higlass"]["hostname"]    # includes https://

    file_uuid = uuid4()

    """ example POST for chromatin data
    curl --user user:pass -F 'datafile=@/path/to/file.mcool' \
        -F 'filetype=cooler' -F 'datatype=matrix' -F 'coordSystem=""' \
        hostname/api/v1/tilesets/
    """

    url = f"{hostname}/api/v1/tilesets/"
    auth = (user, pw)
    files = {'datafile': (mcool_path.name, mcool_path.open('rb'))}
    data = {
        'filetype': 'cooler',
        'datatype': 'matrix',
        'coordSystem': '',
        'uid': str(file_uuid)
    }

        # response should be 201 status and JSON should include uid property
    try:
        response = requests.post(url, auth=auth, files=files, data=data)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)
        print(f"Ingested {mcool_file} into HiGlass at {hostname}", file=sys.stderr)

        # Assuming the response contains the dataset URL or UID
        uid = response.json().get('uid', None)
        if not uid:
            raise ValueError("Failed to retrieve UID from HiGlass response")

        resp_url = f"{hostname}/api/v1/tileset_info/?d={uid}"
        return resp_url

    except requests.RequestException as e:
        print(f"Error ingesting {mcool_file} into HiGlass: {e}", file=sys.stderr)
        raise


def create_tabix_indexed_file(bgzipped_path: Path, file_type : str="bed") -> None:
    """
    Tabix-indexes a genomic file.
    bgzipped_file: Path to the bgzipped input file.
    file_type: Type of the file (e.g., "bed", "vcf", "gff").
    """

    tabix_cmd = f"tabix -s 1 -b 2 -e 3 -p {file_type} {bgzipped_path.as_posix()}"
    subprocess.run(shlex.split(tabix_cmd), check=True)

def validate_hub_contents(hub_json: dict, track_stanzas: list) -> bool:
    """
    Validates the contents of a track hub by checking the hub.txt content and the track stanzas.

    This function performs basic validation to ensure that the hub.txt content contains required fields
    and that each track stanza includes necessary properties. It checks for the presence of essential
    properties such as "hub", "shortLabel", "longLabel", "email", and "genomesFile" in the hub.txt content,
    as well as "track", "type", and "bigDataUrl" in each track stanza.

    Args:
        hub_json (str): The content of the hub.txt file as a string.
        track_stanzas (list): A list of dictionaries representing the track stanzas.

    Returns:
        bool: True if the hub contents are valid, False otherwise.
    """

    # Basic validation for hub.txt content. Pre-validated in the client code
    required_hub_fields = ["hub", "shortLabel", "longLabel", "email", "useOneFile", "genome"]
    for field in required_hub_fields:
        if field not in hub_json:
            print(f"Missing required field '{field}' in hub.txt content.", file=sys.stderr)
            return False

    # Basic validation for each track stanza. Pre-validated in the client code
    for track in track_stanzas:
        required_track_fields = ["track", "type", "bigDataUrl", "shortLabel", "longLabel", "visibility"]
        for field in required_track_fields:
            if field not in track:
                print(f"Missing required field '{field}' in track stanza: {track}", file=sys.stderr)
                return False
            if track["type"] not in VALID_TYPES + VALID_CONTAINER_TYPES:
                print(f"Invalid track type '{track['type']}' in track stanza: {track}", file=sys.stderr)
                return False

    return True

class TrackHubCopy(Resource):
    def post(self, share_uid):
        req = request.get_json()
        if req is None:
            return {"success": False, "message": "Invalid JSON body"}, 400

        hub_json = req.get("hub_json", None)
        assembly = req.get("assembly", None)
        track_stanzas = req.get("tracks", None)
        dry_run = req.get("dry_run", False)

        result = {
            "success": False,
            "message": ""
        }

        if not hub_json:
            result["message"] = "Missing hub info in request body"
            return result, 400
        if not assembly:
            result["message"] = "Missing assembly in request body"
            return result, 400
        if not track_stanzas:
            result["message"] = "Missing list of tracks in request body"
            return result, 400

        # Our destination hub.txt will write everything in one file.
        if not hub_json.get("useOneFile", "off") == "on":
            hub_json["useOneFile"] = "on"
            hub_json["genome"] = assembly
            hub_json.pop("genomesFile", None)  # Remove genomesFile if it exists

        is_valid = validate_hub_contents(hub_json, track_stanzas)
        if not is_valid:
            result["message"] = "Hub contents failed validation. Please check the hub_json and track stanzas for required fields."
            return result, 400

        staging_area = user_upload_file_base / share_uid

        dest_hub = staging_area / "hub.txt"
        if not dry_run:
            staging_area.mkdir(parents=True, exist_ok=True)
            with open(dest_hub, 'w') as f:
                json.dump(hub_json, f, indent=4)
        else:
            print(f'Dry run - hub.txt content written to staging area {dest_hub}', file=sys.stderr)

        # For each track, we will need to validate that the bigDataUrl is reachable and attempt to copy to our staging area
        for track in track_stanzas:
            big_data_url = track.get("bigDataUrl", None)
            if not big_data_url:
                result["message"] = f"Track {track.get('name', 'unknown')} is missing bigDataUrl"
                return result, 400

            # Validate that the bigDataUrl is reachable
            try:
                response = requests.head(big_data_url)
                response.raise_for_status()
            except requests.RequestException as e:
                result["message"] = f"Error reaching bigDataUrl {big_data_url}: {e}"
                return result, 400

            # Download to staging area
            dest_path = staging_area / Path(big_data_url).name
            if not dry_run:
                download_large_file(big_data_url, dest_path)
            else:
                print(f'Dry run - File {big_data_url} downloaded to staging area {dest_path}', file=sys.stderr)

            track["bigDataUrl"] = Path(big_data_url).name

            if track["type"] == "bigBed":
                # Convert bigBed to bed, bgzip, and tabix index
                if not dry_run:
                    conversion_success = bigbed_to_bed(dest_path, staging_area)
                    if not conversion_success:
                        result["message"] = f"Error converting bigBed file {dest_path} to bed format."
                        return result, 500
                else:
                    print(f'Dry run - Converted bigBed file {dest_path} to bed format in staging area', file=sys.stderr)
            elif track["type"] == "hic":
                # Convert .hic to .mcool and ingest into HiGlass
                if not dry_run:
                    mcool_success = hic_to_mcool(dest_path, staging_area)
                    if not mcool_success:
                        result["message"] = f"Error converting .hic file {dest_path} to .mcool format."
                        return result, 500

                    mcool_path = staging_area / dest_path.with_suffix('.mcool').name
                    higlass_url = ingest_mcool_into_higlass(mcool_path)
                    track["bigDataUrl"] = higlass_url
                else:
                    print(f'Dry run - Converted .hic file {dest_path} to .mcool format and ingested into HiGlass', file=sys.stderr)

        result["success"] = True
        result["message"] = "Track hub copied and processed successfully."
        return result, 200

class TrackHubStatus(Resource):
    def post(self, share_uid):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        if req is None:
            return {"success": False, "message": "Invalid JSON body"}, 400
        url = req.get("url")

        raise(NotImplementedError("Track hub status checking is not yet implemented. This will require tracking the status of the track hub copying and conversion processes, and returning that status here."))

        # ...status logic...
        return {"success": True, "message": "Status update"}
