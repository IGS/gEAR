"""
Track Hub utilities for copying and processing track data.
Shared between API and RabbitMQ consumers.
"""

import json
import shlex
import shutil
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Optional
from uuid import uuid4

import pyBigWig
import requests
from Bio import bgzf

VALID_TYPES = ["bigWig", "bigBed", "hic", "vcfTabix"]
VALID_CONTAINER_TYPES = ["multiWig"]

def write_status(
    status_file: Path,
    job_id: str,
    status: str,
    progress: int,
    completed_tracks: int,
    total_tracks: int,
    message: str = "",
    track_statuses: Optional[dict] = None,
) -> None:
    """Update status.json file with progress. Overwrites existing file."""
    status_data = {
        "job_id": job_id,
        "status": status,
        "progress": progress,
        "completed_tracks": completed_tracks,
        "total_tracks": total_tracks,
        "message": message,
        "track_statuses": track_statuses or {},
    }
    with open(status_file, 'w') as f:
        json.dump(status_data, f, indent=4)

def process_trackhub_synchronously(
    job_id: str,
    share_uid: str,
    staging_area: Path,
    status_file: Path,
    hub_json: dict,
    assembly: str,
    track_stanzas: list,
    dry_run: bool = False,
    higlass_config: Optional[dict] = None,
) -> dict:
    """
    Process trackhub synchronously (blocking call).
    Used when RabbitMQ is not available.

    Returns:
        dict: Result with 'success' bool and 'message' string
    """
    processor = TrackHubProcessor(
        job_id=job_id,
        share_uid=share_uid,
        staging_area=staging_area,
        status_file=status_file,
        higlass_config=higlass_config or {},
    )

    return processor.process(hub_json, assembly, track_stanzas, dry_run)

def download_large_file(url: str, destination: Path) -> None:
    """Download a large file from a URL using streaming to save memory."""
    try:
        # Get expected file size from server
        head_response = requests.head(url, timeout=10)
        head_response.raise_for_status()
        expected_size = int(head_response.headers.get('content-length', 0))

        # Check if file already exists and is complete
        if destination.is_file():
            actual_size = destination.stat().st_size
            if expected_size > 0 and actual_size == expected_size:
                print(
                    f"{destination} already exists and is complete ({actual_size} bytes). "
                    "Skipping download.",
                    file=sys.stderr,
                )
                return
            elif actual_size > 0 and actual_size < expected_size:
                print(
                    f"{destination} is incomplete ({actual_size}/{expected_size} bytes). "
                    "Redownloading...",
                    file=sys.stderr,
                )
                # Delete file to force redownload
                destination.unlink()

        # Download file
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            with open(destination, 'wb') as file:
                shutil.copyfileobj(response.raw, file)

        # Verify downloaded file size matches expected
        actual_size = destination.stat().st_size
        if expected_size > 0 and actual_size != expected_size:
            destination.unlink()
            raise Exception(
                f"Downloaded file size mismatch: expected {expected_size} bytes, "
                f"got {actual_size} bytes"
            )
    except requests.exceptions.RequestException as e:
        raise Exception(f"Error downloading from {url}: {e}") from e


def bigbed_to_bed(bigbed_path: Path, outdir_path: Path) -> bool:
    """Convert a BigBed file to a BED file and compress it using BGZF."""
    bed_path = outdir_path / bigbed_path.with_suffix('.bed').name
    if bed_path.is_file():
        print(f"{bed_path} already exists. Skipping conversion.", file=sys.stderr)
        return True

    try:
        bb = pyBigWig.open(bigbed_path.as_posix())
        if not bb.isBigBed():
            print(f"{bigbed_path} is not a bigBed file.", file=sys.stderr)
            bb.close()
            return False

        with open(bed_path, 'w') as bed_out:
            for chrom, start, end, rest in bb.intervals():
                bed_line = f"{chrom}\t{start}\t{end}\t{rest}\n"
                bed_out.write(bed_line)

        bb.close()
        print(f"Converted {bigbed_path} to {bed_path}.", file=sys.stderr)
    except Exception as e:
        print(f"Error converting {bigbed_path} to bed: {e}", file=sys.stderr)
        return False

    try:
        gz_path = bed_path.with_suffix(bed_path.suffix + '.gz')
        with open(bed_path, 'rb') as f_in, bgzf.BgzfWriter(gz_path) as f_out:
            f_out.write(f_in.read())

        _create_tabix_indexed_file(gz_path, file_type="bed")
        bed_path.unlink()  # Clean up uncompressed version
        return True
    except Exception as e:
        print(f"Error compressing and indexing {bed_path}: {e}", file=sys.stderr)
        return False


def hic_to_mcool(hic_path: Path, outdir_path: Path) -> bool:
    """Convert a .hic file to a .mcool file using the hic2cool tool."""
    mcool_path = outdir_path / hic_path.with_suffix('.mcool').name
    if mcool_path.is_file():
        print(f"{mcool_path} already exists. Skipping conversion.", file=sys.stderr)
        return True

    try:
        from hic2cool import hic2cool_convert
        hic2cool_convert(hic_path.as_posix(), mcool_path.as_posix(), 0)
        print(f"Converted {hic_path} to {mcool_path}.", file=sys.stderr)
        return True
    except Exception as e:
        print(f"Error converting {hic_path} to mcool: {e}", file=sys.stderr)
        return False

def ingest_mcool_into_higlass(mcool_path: Path, config: dict, assembly: str) -> str:
    """Ingest a .mcool file into HiGlass and return the URL."""
    file_uuid = uuid4()

    url = f"{config['higlass_hostname']}/api/v1/tilesets/"
    auth = (config['higlass_admin_user'], config['higlass_admin_pass'])

    try:
        # Ensure file stays open during the duration of the request.
        with open(mcool_path, 'rb') as f:
            # All parameters must be in 'files' dict for multipart/form-data
            # Non-file fields are sent as tuples (None, value)
            files = {
                'datafile': (mcool_path.name, f),
                'filetype': (None, 'cooler'),
                'datatype': (None, 'matrix'),
                'coordSystem': (None, assembly),
                'uid': (None, str(file_uuid)),
            }
            # print an example curl-equivalent command for testing
            #print(f"Example curl command:\ncurl -X POST {url} -u {auth[0]}:{auth[1]} -F 'datafile=@{mcool_path}' -F 'filetype=cooler' -F 'datatype=matrix' -F 'coordSystem={assembly}' -F 'uid={file_uuid}'")

            response = requests.post(url, auth=auth, files=files, timeout=600)
            print(f"HiGlass response status: {response.status_code}", file=sys.stderr)
            print(f"HiGlass response content: {response.text}", file=sys.stderr)
            response.raise_for_status()

        print(f"Ingested {mcool_path} into HiGlass", file=sys.stderr)

        uid = response.json().get('uid')
        if not uid:
            raise ValueError("Failed to retrieve UID from HiGlass response")

        resp_url = f"{config['higlass_hostname']}/api/v1/tileset_info/?d={uid}"
        return resp_url

    except requests.Timeout:
        print(f"Timeout ingesting {mcool_path} into HiGlass. Checking if file is available...", file=sys.stderr)
        return check_higlass_url(config, file_uuid)
    except requests.RequestException as e:
        print(f"Error ingesting {mcool_path} into HiGlass: {e}", file=sys.stderr)
        raise

def check_higlass_url(config, file_uuid) -> str:
    try:
        auth = (config['higlass_admin_user'], config['higlass_admin_pass'])
        check_url = f"{config['higlass_hostname']}/api/v1/tileset_info/?d={file_uuid}"
        check_response = requests.get(check_url, auth=auth, timeout=10)
        if check_response.status_code == 200:
            print(
                f"File was successfully ingested despite timeout. UID: {file_uuid}",
                file=sys.stderr,
            )
            return check_url
    except Exception as check_error:
        print(f"Failed to verify ingestion: {check_error}", file=sys.stderr)
    raise Exception(f"File ingestion likely failed and timed out. UID: {file_uuid}")



def delete_higlass_tileset(tileset_uid: str, config: dict) -> bool:
    """Attempt to delete a tileset from HiGlass."""
    try:
        url = f"{config['higlass_hostname']}/api/v1/tilesets/{tileset_uid}/"
        auth = (config['higlass_admin_user'], config['higlass_admin_pass'])
        response = requests.delete(url, auth=auth, timeout=30)

        if response.status_code in [204, 200]:
            print(f"Deleted HiGlass tileset {tileset_uid}", file=sys.stderr)
            return True
        else:
            print(
                f"Failed to delete HiGlass tileset {tileset_uid}: {response.status_code}",
                file=sys.stderr,
            )
            return False
    except Exception as e:
        print(f"Error deleting HiGlass tileset {tileset_uid}: {e}", file=sys.stderr)
        return False

def _create_tabix_indexed_file(bgzipped_path: Path, file_type: str = "bed") -> None:
    """Tabix-index a genomic file."""
    tabix_cmd = f"tabix -s 1 -b 2 -e 3 -p {file_type} {bgzipped_path.as_posix()}"
    subprocess.run(shlex.split(tabix_cmd), check=True)


def _append_higlass_url_to_track(hub_file: Path, hic_track_name: str, higlass_url: str) -> None:
    """Inject a gos_higlass_url property into a TrackDB file."""
    modified_hub_file = hub_file.with_stem(hub_file.stem + '.modified')

    with open(modified_hub_file, 'w') as out_f:
        for line in hub_file.read_text().splitlines():
            line = line.strip()
            if line.startswith("bigDataUrl") and line.endswith(hic_track_name):
                out_f.write(f"gos_higlass_url {higlass_url}\n")
            out_f.write(line + "\n")

    modified_hub_file.replace(hub_file)


def validate_hub_contents(hub_json: dict, track_stanzas: list) -> bool:
    """Validate hub and track configurations."""
    required_hub_fields = ["hub", "shortLabel", "longLabel", "email", "useOneFile", "genome"]
    for field in required_hub_fields:
        if field not in hub_json:
            print(f"Missing required field '{field}' in hub.txt content.", file=sys.stderr)
            return False

    for track in track_stanzas:
        required_track_fields = ["track", "type", "bigDataUrl", "shortLabel", "longLabel", "visibility"]
        for field in required_track_fields:
            if field not in track:
                print(f"Missing required field '{field}' in track stanza.", file=sys.stderr)
                return False
        if track["type"] not in VALID_TYPES + VALID_CONTAINER_TYPES:
            print(f"Invalid track type '{track['type']}'.", file=sys.stderr)
            return False

        # if human assembly, disallow VCF and HIC types for privacy reasons
        if hub_json.get("genome", "").lower() in ["hg19", "hg38"] and track["type"] in ["vcfTabix", "hic"]:
            print(f"Track type '{track['type']}' is not allowed for human assemblies due to privacy concerns.", file=sys.stderr)
            return False

    return True


class TrackHubProcessor:
    """Process track hub downloads and conversions."""

    def __init__(
        self,
        job_id: str,
        share_uid: str,
        staging_area: Path,
        status_file: Path,
        higlass_config: Optional[dict] = None,
    ):
        self.job_id = job_id
        self.share_uid = share_uid
        self.staging_area = staging_area
        self.status_file = status_file
        self.higlass_config = higlass_config or {}
        self.staging_area.mkdir(parents=True, exist_ok=True)

    def update_status(
        self,
        status: str,
        progress: int,
        completed_tracks: int,
        total_tracks: int,
        message: str = "",
        track_statuses: Optional[dict] = None,
    ) -> None:
        write_status(
            self.status_file,
            self.job_id,
            status,
            progress,
            completed_tracks,
            total_tracks,
            message,
            track_statuses,
        )

    def process(
        self,
        hub_json: dict,
        assembly: str,
        track_stanzas: list,
        dry_run: bool = False,
    ) -> dict:
        """
        Process trackhub: download files, convert formats, ingest into HiGlass.
        Returns result dict with success status and message.
        """

        total_tracks = len(track_stanzas)
        track_statuses = {}
        ingested_tilesets = []  # Track successfully ingested tilesets for cleanup

        try:

            # Initialize status
            self.update_status(
                "processing",
                0,
                0,
                total_tracks,
                "Starting track hub processing...",
                track_statuses,
            )

            # Ensure oneFile mode
            if not hub_json.get("useOneFile") == "on":
                hub_json["useOneFile"] = "on"
            hub_json["genome"] = assembly
            hub_json.pop("genomesFile", None)

            # Validate
            if not validate_hub_contents(hub_json, track_stanzas):
                raise ValueError("Hub contents failed validation")

            # Write hub.txt as a standard hub file (without trackDb.txt) for processing.
            # Add an extra newline before and after the "genome" tag.
            dest_hub = self.staging_area / "hub.txt"
            if not dry_run:
                with open(dest_hub, 'w') as f:
                    for key, value in hub_json.items():
                        if key == "genome":
                            f.write("\n")
                        f.write(f"{key} {value}\n")
                        # tracks will fill out the extra newline

            # Process each track
            completed = 0
            for idx, track in enumerate(track_stanzas):
                track_name = track.get("shortLabel", f"Track {idx + 1}")
                track_statuses[track_name] = "downloading"

                self.update_status(
                    "processing",
                    int((idx / total_tracks) * 100),
                    completed,
                    total_tracks,
                    f"Downloading {track_name}...",
                    track_statuses,
                )

                big_data_url = track.get("bigDataUrl")
                if not big_data_url:
                    raise ValueError(f"Track {track_name} missing bigDataUrl")

                # Validate URL is reachable
                try:
                    response = requests.head(big_data_url, timeout=10)
                    response.raise_for_status()
                except requests.RequestException as e:
                    raise Exception(f"Cannot reach {big_data_url}: {e}")

                # Download
                dest_path = self.staging_area / Path(big_data_url).name
                if not dry_run:
                    download_large_file(big_data_url, dest_path)

                track["bigDataUrl"] = Path(big_data_url).name
                track_statuses[track_name] = "downloaded"

                # Convert based on type
                if track["type"] == "bigBed":
                    track_statuses[track_name] = "converting"
                    self.update_status(
                        "processing",
                        int((idx / total_tracks) * 100),
                        completed,
                        total_tracks,
                        f"Converting {track_name} (bigBed → BED)...",
                        track_statuses,
                    )
                    if not dry_run:
                        if not bigbed_to_bed(dest_path, self.staging_area):
                            raise Exception(f"Failed to convert {track_name} from bigBed to BED")

                elif track["type"] == "hic":
                    track_statuses[track_name] = "converting"
                    self.update_status(
                        "processing",
                        int((idx / total_tracks) * 100),
                        completed,
                        total_tracks,
                        f"Converting {track_name} (HIC → MCool)...",
                        track_statuses,
                    )
                    if not dry_run:
                        if not hic_to_mcool(dest_path, self.staging_area):
                            raise Exception(f"Failed to convert {track_name} from HIC to MCool")

                        track_statuses[track_name] = "ingesting"
                        self.update_status(
                            "processing",
                            int((idx / total_tracks) * 100),
                            completed,
                            total_tracks,
                            f"Ingesting {track_name} into HiGlass...",
                            track_statuses,
                        )

                        print(f"Ingesting {track_name} into HiGlass...", file=sys.stderr)
                        mcool_path = self.staging_area / dest_path.with_suffix('.mcool').name
                        try:
                            higlass_url = ingest_mcool_into_higlass(
                                mcool_path, self.higlass_config, assembly
                            )
                            tileset_uid = higlass_url.split("d=")[-1]
                            ingested_tilesets.append(tileset_uid)
                            _append_higlass_url_to_track(dest_hub, dest_path.name, higlass_url)
                        except Exception:
                            traceback.print_exc()
                            print(
                                f"Ingestion failed for {track_name}. Cleaning up...",
                                file=sys.stderr,
                            )
                            raise

                track_statuses[track_name] = "completed"
                completed += 1

            # append track info to hub.txt (useOneFile on)
            if not dry_run:
                with open(dest_hub, 'a') as f:
                    for track in track_stanzas:
                        f.write("\n")
                        for key, value in track.items():
                            f.write(f"{key} {value}\n")

            self.update_status(
                "complete",
                100,
                completed,
                total_tracks,
                "Track hub processed successfully",
                track_statuses,
            )

            return {"success": True, "message": "Track hub processed successfully"}

        except Exception as e:
            print(f"Error processing track hub ID {self.job_id}: {e}", file=sys.stderr)

            # Clean up ingested tilesets on failure
            if ingested_tilesets and self.higlass_config:
                print(
                    f"Cleaning up {len(ingested_tilesets)} ingested tileset(s)...",
                    file=sys.stderr,
                )
                for tileset_uid in ingested_tilesets:
                    delete_higlass_tileset(tileset_uid, self.higlass_config)

            self.update_status(
                "error",
                0,
                0,
                total_tracks,
                f"Error: {e}",
                track_statuses,
            )
            return {"success": False, "message": str(e)}