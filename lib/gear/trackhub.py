"""
Track Hub utilities for copying and processing track data.
Shared between API and RabbitMQ consumers.
"""

import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Optional
from uuid import uuid4

import pyBigWig
import requests
import shutil
from Bio import bgzf

VALID_TYPES = ["bigWig", "bigBed", "hic", "vcfTabix"]
VALID_CONTAINER_TYPES = ["multiWig"]

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
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            with open(destination, 'wb') as file:
                shutil.copyfileobj(response.raw, file)
    except requests.exceptions.RequestException as e:
        raise Exception(f"Error downloading from {url}: {e}") from e


def bigbed_to_bed(bigbed_path: Path, outdir_path: Path) -> bool:
    """Convert a BigBed file to a BED file and compress it using BGZF."""
    bed_path = outdir_path / bigbed_path.with_suffix('.bed').name

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

    try:
        from hic2cool import hic2cool_convert
        hic2cool_convert(hic_path.as_posix(), mcool_path.as_posix(), 0)
        print(f"Converted {hic_path} to {mcool_path}.", file=sys.stderr)
        return True
    except Exception as e:
        print(f"Error converting {hic_path} to mcool: {e}", file=sys.stderr)
        return False


def ingest_mcool_into_higlass(mcool_path: Path, config: dict) -> str:
    """Ingest a .mcool file into HiGlass and return the URL."""
    file_uuid = uuid4()

    url = f"{config['higlass_hostname']}/api/v1/tilesets/"
    auth = (config['higlass_admin_user'], config['higlass_admin_pass'])
    files = {'datafile': (mcool_path.name, open(mcool_path, 'rb'))}
    data = {
        'filetype': 'cooler',
        'datatype': 'matrix',
        'coordSystem': '',
        'uid': str(file_uuid)
    }

    try:
        response = requests.post(url, auth=auth, files=files, data=data)
        response.raise_for_status()
        print(f"Ingested {mcool_path} into HiGlass", file=sys.stderr)

        uid = response.json().get('uid')
        if not uid:
            raise ValueError("Failed to retrieve UID from HiGlass response")

        resp_url = f"{config['higlass_hostname']}/api/v1/tileset_info/?d={uid}"
        return resp_url

    except requests.RequestException as e:
        print(f"Error ingesting {mcool_path} into HiGlass: {e}", file=sys.stderr)
        raise


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
        """Update status.json file with progress."""
        status_data = {
            "job_id": self.job_id,
            "status": status,
            "progress": progress,
            "completed_tracks": completed_tracks,
            "total_tracks": total_tracks,
            "message": message,
            "track_statuses": track_statuses or {},
        }
        with open(self.status_file, 'w') as f:
            json.dump(status_data, f)

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
        try:
            total_tracks = len(track_stanzas)
            track_statuses = {}

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

            # Write hub.txt
            dest_hub = self.staging_area / "hub.txt"
            if not dry_run:
                with open(dest_hub, 'w') as f:
                    json.dump(hub_json, f, indent=4)

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
                            raise Exception(f"Failed to convert {track_name}")

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
                            raise Exception(f"Failed to convert {track_name}")

                        track_statuses[track_name] = "ingesting"
                        self.update_status(
                            "processing",
                            int((idx / total_tracks) * 100),
                            completed,
                            total_tracks,
                            f"Ingesting {track_name} into HiGlass...",
                            track_statuses,
                        )

                        mcool_path = self.staging_area / dest_path.with_suffix('.mcool').name
                        higlass_url = ingest_mcool_into_higlass(
                            mcool_path, self.higlass_config
                        )
                        _append_higlass_url_to_track(dest_hub, dest_path.name, higlass_url)

                track_statuses[track_name] = "completed"
                completed += 1

            # Write trackDb.txt
            if not dry_run:
                with open(self.staging_area / "trackDb.txt", 'w') as f:
                    for track in track_stanzas:
                        f.write(f"track {track['track']}\n")
                        for key, value in track.items():
                            if key != 'track':
                                f.write(f"{key} {value}\n")
                        f.write("\n")

            self.update_status(
                "completed",
                100,
                completed,
                total_tracks,
                "Track hub processed successfully",
                track_statuses,
            )

            return {"success": True, "message": "Track hub processed successfully"}

        except Exception as e:
            print(f"Error processing track hub {self.job_id}: {e}", file=sys.stderr)
            self.update_status(
                "failed",
                0,
                0,
                total_tracks,
                f"Error: {e}",
                track_statuses,
            )
            return {"success": False, "message": str(e)}