import configparser
import re
import shlex
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
BIGBED_EXTENSIONS = [".bb", ".bigbed"]

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


def fetch_trackdb_url(genomes_txt, assembly) -> str:
    """
    Extract the trackDb URL for a specified assembly from a genomes.txt file.

    Args:
        genomes_txt (str): The contents of a genomes.txt file as a string.
        assembly (str): The name of the assembly/genome to search for.

    Returns:
        str: The trackDb URL corresponding to the specified assembly.

    Raises:
        ValueError: If the specified assembly is not found in the genomes.txt content.

    Example:
        >>> genomes_txt = "genome hg38\\ntrackDb http://example.com/hg38/trackDb\\ngenome mm10\\ntrackDb http://example.com/mm10/trackDb"
        >>> fetch_trackdb_url(genomes_txt, "hg38")
        'http://example.com/hg38/trackDb'
    """

    lines = genomes_txt.splitlines()
    for i, genome_line in enumerate(lines):
        if genome_line.startswith(f"genome {assembly}"):
            # The next line should contain the trackDb URL
            if i + 1 < len(lines):
                next_line = lines[i + 1]
                if next_line.startswith("trackDb"):
                    return next_line.split(" ")[1]
    raise ValueError(f"Assembly {assembly} not found in genomes.txt")

def validate_track_types_from_db(trackdb_txt, trackdb_url) -> list:
    """Parse track names from a UCSC trackDb.txt content.

    Returns a list of dicts with track information.

    Example track:
    track P1HC_ATAC_1
    bigDataUrl P1HC_ATAC_1.bigwig
    shortLabel ATAC-seq 1st replicate
    longLabel ATAC-seq 1st replicate
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
        if req is None:
            return {"success": False, "message": "Invalid JSON body"}, 400
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

        trackdb_url = fetch_trackdb_url(genomes_response.text, assembly)
        if not trackdb_url:
            result["message"] = f"No trackDb found for assembly {assembly}"
            return result, 400
        if not trackdb_url.startswith("http://") and not trackdb_url.startswith("https://"):
            trackdb_url = f"{base_url}/{trackdb_url}"

        # Fetch tracks
        try:
            trackdb_response = requests.get(trackdb_url)
            trackdb_response.raise_for_status()
        except requests.RequestException as e:
            result["message"] = f"Error fetching trackDb: {str(e)}"
            return result, 500

        # Does not count sub-tracks
        result["num_tracks"] = len(re.findall(r"^track ", trackdb_response.text, re.MULTILINE))
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
        if req is None:
            return {"success": False, "message": "Invalid JSON body"}, 400

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

        trackdb_file: Path = assembly_dir / "trackDb.txt"

        # find all bigBed files in the assembly directory
        bigbed_paths = []
        for ext in BIGBED_EXTENSIONS:
            bigbed_paths.extend(assembly_dir.rglob(f"*{ext}"))
        for bigbed_path in bigbed_paths:
            conversion_success = bigbed_to_bed(bigbed_path, assembly_dir)
            if not conversion_success:
                result["message"] = f"Failed to convert {bigbed_path.as_posix()} to bed"
                return result, 500

        # convert and ingest .hic files
        hic_paths = list(assembly_dir.rglob("*.hic"))
        for hic_path in hic_paths:
            conversion_success = hic_to_mcool(hic_path, assembly_dir)
            if not conversion_success:
                result["message"] = f"Failed to convert {hic_path.as_posix()} to mcool"
                return result, 500

            mcool_path = assembly_dir / hic_path.with_suffix('.mcool').name
            try:
                higlass_url = ingest_mcool_into_higlass(mcool_path)
                # TODO: add higlass_url as the "gos_higlass_url" property for this trackdb entry.
                print(f"Ingested {mcool_path.as_posix()} into HiGlass at {higlass_url}", file=sys.stderr)
                # find "hic" track and add the "gos_higlass_url" property with the higlass_url value
                trackdb_file = append_higlass_url_to_trackdb(trackdb_file, hic_path.name, higlass_url)

            except Exception as e:
                result["message"] = f"Failed to ingest {mcool_path.as_posix()} into HiGlass: {str(e)}"
                return result, 500

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