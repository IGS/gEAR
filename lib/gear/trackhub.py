"""
Track Hub utilities for copying and processing track data.
Shared between API and RabbitMQ consumers.
"""

import ipaddress
import json
import shlex
import shutil
import socket
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Optional
from urllib.parse import urljoin, urlparse
from uuid import uuid4

import pyBigWig
import requests
from Bio import bgzf

VALID_TYPES = ["bigWig", "bigBed", "hic", "vcfTabix"]
VALID_CONTAINER_TYPES = ["multiWig"]

HUB_FIELDS = ["hub", "shortLabel", "longLabel", "email", "useOneFile", "genome", "genomesFile"]
TRACK_FIELDS = ["type", "bigDataUrl", "shortLabel", "longLabel", "visibility", "color", "autoScale"]

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
    hub_url: str,
    higlass_config: Optional[dict] = None,
    dry_run: bool = False,
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
        hub_url=hub_url
    )

    return processor.process(hub_json, assembly, track_stanzas, dry_run)

def fetch_trackdb_path(genomes_txt: str, assembly: str) -> str:
    """
    Extract the trackDb path for a specified genome assembly from UCSC genomes.txt format.

    Args:
        genomes_txt (str): The contents of a UCSC genomes.txt file as a string.
        assembly (str): The genome assembly identifier to search for (e.g., 'hg38', 'mm10').

    Returns:
        str: The trackDb path for the specified assembly.

    Raises:
        ValueError: If the assembly is not found or trackDb entry is missing.

    """
    for line in genomes_txt.splitlines():
        if line.startswith(f"genome {assembly}"):
            lines = genomes_txt.splitlines()
            next_line_index = lines.index(line) + 1
            if next_line_index < len(lines):
                next_line = lines[next_line_index]
                if next_line.startswith("trackDb"):
                    return next_line.split(" ", 1)[1]

    raise ValueError(f"Assembly {assembly} not found in genomes.txt or trackDb entry is missing.")

def _is_safe_public_http_url(url: str) -> bool:
    """Return True if URL is HTTP(S) and resolves only to public IP addresses."""
    try:
        parsed = urlparse(url)
        if parsed.scheme not in ("http", "https") or not parsed.hostname:
            return False

        host = parsed.hostname
        addrinfo = socket.getaddrinfo(host, None)

        for info in addrinfo:
            ip_str = info[4][0]
            ip_obj = ipaddress.ip_address(ip_str)
            if (
                ip_obj.is_private
                or ip_obj.is_loopback
                or ip_obj.is_link_local
                or ip_obj.is_multicast
                or ip_obj.is_reserved
                or ip_obj.is_unspecified
            ):
                return False

        return True
    except Exception:
        return False


def _normalize_track_dict(track: dict, resolve_urls: bool = False, trackdb_url: str = "") -> dict:
    """
    Normalize a track dictionary to consistent field names and formats.

    Args:
        track (dict): Raw track dictionary from parsing.
        resolve_urls (bool): Whether to resolve relative URLs (for trackDb parsing).
        trackdb_url (str): Base URL for resolving relative URLs.

    Returns:
        dict: Normalized track dictionary with consistent keys.
    """
    normalized = {}

    # Map both "track" and "name" keys to "track"
    if "track" in track:
        normalized["track"] = track["track"]
    elif "name" in track:
        normalized["track"] = track["name"]

    # Copy standard fields
    for key in ["type", "shortLabel", "longLabel", "visibility", "autoScale", "color"]:
        if key in track:
            normalized[key] = track[key]

    # Handle bigDataUrl with optional URL resolution
    if "bigDataUrl" in track:
        url = track["bigDataUrl"]
        if resolve_urls:
            # Resolve relative URLs
            if not url.startswith("http://") and not url.startswith("https://"):
                url = urljoin(trackdb_url, url)
            # Follow redirects only for safe public HTTP(S) URLs
            try:
                if _is_safe_public_http_url(url):
                    response = requests.head(url, allow_redirects=True, timeout=10)
                    if response.status_code == 200 and _is_safe_public_http_url(response.url):
                        url = response.url
                else:
                    print(
                        f"WARNING: Skipping unsafe URL resolution for track '{normalized.get('track')}': {url}",
                        file=sys.stderr,
                    )
            except Exception as e:
                print(f"WARNING: Could not resolve URL for track '{normalized.get('track')}': {e}", file=sys.stderr)
        normalized["bigDataUrl"] = url

    # Normalize color format to rgb()
    if "color" in normalized and not normalized["color"].startswith("rgb("):
        normalized["color"] = f"rgb({normalized['color']})"

    return normalized

def _parse_track_stanzas(lines, resolve_urls: bool = False, trackdb_url: str = "") -> list:
    """
    Internal helper to parse track stanzas from any line-based format.

    Handles the common pattern of parsing "track" blocks with consistent field extraction.
    Used by both trackDb.txt and hub.txt (useOneFile) parsers.

    Args:
        lines (list): Lines to parse (typically from splitlines()).
        resolve_urls (bool): Whether to resolve relative URLs.
        trackdb_url (str): Base URL for URL resolution.

    Returns:
        list: List of normalized track dictionaries.
    """
    track_list = []
    current_track = {}

    for line in lines:
        line = line.strip()

        # Skip empty lines and finalize current track if present
        if not line:
            if current_track:
                track_list.append(_normalize_track_dict(current_track, resolve_urls, trackdb_url))
                current_track = {}
            continue

        # Start of new track stanza
        if line.startswith("track "):
            if current_track:
                track_list.append(_normalize_track_dict(current_track, resolve_urls, trackdb_url))
            current_track = {"track": line.split(" ", 1)[1]}

        # Parse track fields (only if we're inside a track stanza)
        elif current_track:
            for field in TRACK_FIELDS:
                if line.startswith(f"{field} "):
                    current_track[field] = line.split(" ", 1)[1]
                    break

    # Don't forget the last track
    if current_track:
        track_list.append(_normalize_track_dict(current_track, resolve_urls, trackdb_url))

    return track_list

def parse_tracks_from_trackdb(trackdb_txt: str, trackdb_url: str) -> list[dict]:
    """
    Parse track stanzas from a UCSC trackDb.txt file.

    Args:
        trackdb_txt (str): The contents of a UCSC trackDb.txt file.
        trackdb_url (str): The base URL of the trackDb.txt file (used to resolve relative URLs).

    Returns:
        list: A list of dictionaries, each representing a track with keys like 'track', 'type',
              'bigDataUrl', 'shortLabel', 'longLabel', 'color', 'visibility', etc.

    Example track dict:
        {
            "name": "P1HC_ATAC_1",
            "type": "bigWig",
            "bigDataUrl": "https://example.com/P1HC_ATAC_1.bigwig",
            "shortLabel": "ATAC-seq 1st replicate",
            "longLabel": "ATAC-seq 1st replicate",
            "color": "rgb(31,119,180)",
            "visibility": "dense"
        }
    """
    return _parse_track_stanzas(trackdb_txt.splitlines(), resolve_urls=True, trackdb_url=trackdb_url)

def parse_hub_from_file(hub_txt: str) -> tuple[dict, list[dict]]:
    """
    Parse hub.txt file in useOneFile mode into hub metadata and track stanzas.

    Args:
        hub_txt (str): The contents of a hub.txt file in useOneFile mode.

    Returns:
        tuple: (hub_dict, track_stanzas_list)
            - hub_dict: Dictionary with hub metadata (hub, shortLabel, longLabel, email, useOneFile, genome)
            - track_stanzas_list: List of track dictionaries

    Notes:
        In useOneFile mode, the hub.txt file contains both hub metadata and track definitions.
        Track stanzas begin with 'track <name>' and are separated by blank lines.
    """
    hub_json = {}
    track_list = []
    lines = hub_txt.splitlines()

    # Parse hub metadata and tracks
    for line in lines:
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Parse hub metadata
        for field in HUB_FIELDS:
            if line.startswith(f"{field} "):
                hub_json[field] = line.split(" ", 1)[1]
                break

        # If it's a track stanza, stop parsing hub metadata
        if line.startswith("track "):
            break

    # Parse track stanzas (rest of the file after hub metadata)
    track_start_idx = next(
        (i for i, l in enumerate(lines) if l.strip().startswith("track ")),
        len(lines)
    )
    track_list = _parse_track_stanzas(lines[track_start_idx:], resolve_urls=False, trackdb_url="")

    return hub_json, track_list

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


def bigbed_to_bed(bigbed_path: Path, outdir_path: Path) -> str:
    """Convert a BigBed file to a BED file and compress it using BGZF."""
    bed_path = outdir_path / bigbed_path.with_suffix('.bed.gz').name
    if bed_path.is_file():
        print(f"{bed_path} already exists. Skipping conversion.", file=sys.stderr)
        return bed_path.as_posix()

    try:
        bb = pyBigWig.open(bigbed_path.as_posix())
        if not bb.isBigBed():
            print(f"{bigbed_path} is not a bigBed file.", file=sys.stderr)
            bb.close()
            raise

        with open(bed_path, 'w') as bed_out:
            for chrom, start, end, rest in bb.intervals():
                bed_line = f"{chrom}\t{start}\t{end}\t{rest}\n"
                bed_out.write(bed_line)

        bb.close()
        print(f"Converted {bigbed_path} to {bed_path}.", file=sys.stderr)
    except Exception as e:
        print(f"Error converting {bigbed_path} to bed: {e}", file=sys.stderr)
        raise

    try:
        gz_path = bed_path.with_suffix(bed_path.suffix + '.gz')
        with open(bed_path, 'rb') as f_in, bgzf.BgzfWriter(gz_path) as f_out:
            f_out.write(f_in.read())

        _create_tabix_indexed_file(gz_path, file_type="bed")
        bed_path.unlink()  # Clean up uncompressed version
        return gz_path.as_posix()
    except Exception as e:
        print(f"Error compressing and indexing {bed_path}: {e}", file=sys.stderr)
        raise


def hic_to_mcool(hic_path: Path, outdir_path: Path) -> bool:
    """Convert a .hic file to a .mcool file using the hic2cool tool."""
    mcool_path = outdir_path / hic_path.with_suffix('.mcool').name
    # TODO: verify this conversion is complete.
    if mcool_path.is_file():
        print(f"{mcool_path} already exists. Skipping conversion.", file=sys.stderr)
        return True

    # NOTE: hic2cool writes to STDOUT, which will show in RabbitMQ message logs if using that.
    # `sudo journalctl -u gosling-upload-consumer.service -f`
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
    #TODO: Figure out a way to check this as a previously uploaded ingest
    file_uuid = uuid4()

    url = f"{config['higlass_hostname']}/api/v1/tilesets/"
    auth = (config['higlass_admin_user'], config['higlass_admin_pass'])

    try:
        # Ensure file stays open during the duration of the request.
        with open(mcool_path, 'rb') as f:
            # All parameters must be in 'files' dict for multipart/form-data
            # Non-file fields are sent as tuples (None, value)
            files = {
                'datafile': (mcool_path.as_posix(), f),
                'filetype': (None, 'cooler'),
                'datatype': (None, 'matrix'),
                'coordSystem': (None, assembly),
                'uid': (None, str(file_uuid)),
            }

            response = requests.post(url, auth=auth, files=files, timeout=660)  # timeout is slightly longer than server settings

            # Check for 504 Gateway Timeout specifically
            if response.status_code == 504:
                print(
                    f"HiGlass returned 504 Gateway Timeout. File may still be processing. "
                    f"Checking if file was ingested (UID: {file_uuid})...",
                    file=sys.stderr,
                )
                return check_higlass_url(config, file_uuid)

            response.raise_for_status()

        # Example response keys (JSON):
        # "uuid", "datafile", "filetype", "datatype", "name", "coordSystem", "coordSystem2",
        # "created", "project", "project_name", "description", "private"

        uid = response.json().get('uuid', None)
        if not uid:
            raise ValueError("Failed to retrieve UID from HiGlass response")

        print(f"Ingested {mcool_path} into HiGlass", file=sys.stderr)

        resp_url = f"{config['higlass_hostname']}/api/v1/tileset_info/?d={uid}"
        return resp_url
    except requests.HTTPError as e:
        # Catch other HTTP errors (400, 401, 403, 500, etc.)
        print(f"HTTP error ingesting {mcool_path} into HiGlass: {e}", file=sys.stderr)
        raise
    except requests.RequestException as e:
        print(f"Error ingesting {mcool_path} into HiGlass: {e}", file=sys.stderr)
        raise
    except Exception as e:
        print(f"Unexpected error ingesting {mcool_path} into HiGlass: {e}", file=sys.stderr)
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
        hub_url: str = ""
    ):
        self.job_id = job_id
        self.share_uid = share_uid
        self.staging_area = staging_area
        self.status_file = status_file
        self.higlass_config = higlass_config or {}
        self.hub_url = hub_url
        if not self.hub_url:
            raise ValueError("Hub URL base must be provided for TrackHubProcessor")

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
            if not hub_json.get("useOneFile", "") == "on":
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
                with open(dest_hub, 'w', buffering=1) as f:
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

                big_data_url = track.get("bigDataUrl", "").strip()
                uploaded_file_name = track.get("uploadedFileName")
                print(f"Processing track '{track_name}' with bigDataUrl: '{big_data_url}' and uploadedFileName: '{uploaded_file_name}'", file=sys.stderr)
                if big_data_url:
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
                elif uploaded_file_name:
                    # File was already saved to staging area during request handling
                    dest_path = self.staging_area / uploaded_file_name
                    if not dest_path.is_file() and not dry_run:
                        raise ValueError(
                            f"Uploaded file for track '{track_name}' not found in staging area: {dest_path}"
                        )
                    # File is already in place, no need to save again
                else:
                    raise ValueError(
                        f"No bigDataUrl or uploadedFileName specified for track '{track_name}'"
                    )


                # Update bigDataUrl to point to a remote reference.
                track["bigDataUrl"] = f"{self.hub_url}/{dest_path.name}"
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
                        try:
                            bed_path = bigbed_to_bed(dest_path, self.staging_area)
                            track["gos_url"] = bed_path
                        except Exception as e:
                            raise Exception(f"Failed to convert {track_name} from bigBed to BED: {e}")

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
                            track['gos_url'] = higlass_url

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
                with open(dest_hub, 'a', buffering=1) as f:
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