"""
SpatialProcessor - Process spatial transcriptomics dataset uploads to Zarr.

Reads a platform-specific spatial data archive (Visium, VisiumHD, Curio, GeoMx,
CosMx, Xenium) using the handler class selected by spatial_format, then writes
the result as a Zarr store.
"""

import json
import shutil
from pathlib import Path

import geardb
from gear.anndata_processor import write_status
from gear.spatialhandler import SPATIALTYPE2CLASS


def process_spatial_synchronously(
    job_id: str,
    share_uid: str,
    staging_area: Path,
    status_file: Path,
    spatial_format: str,
    perform_primary_analysis: bool,
) -> dict:
    """Process a spatial dataset upload. Used by both the queued consumer and the synchronous CGI fallback."""
    status = {
        "job_id": job_id,
        "status": "processing",
        "message": "Initializing dataset processing.",
        "progress": 0,
    }
    write_status(status_file, status)

    metadata_file = staging_area / 'metadata.json'
    if not metadata_file.is_file():
        status["status"] = "error"
        status["message"] = "No metadata JSON file found."
        write_status(status_file, status)
        return {"success": 0, "message": status["message"]}

    with open(metadata_file, 'r') as f:
        metadata = json.load(f)

    sample_taxid = metadata.get("sample_taxid", None)
    organism_id = geardb.get_organism_id_by_taxon_id(sample_taxid)
    filepath = staging_area / f"{share_uid}.tar.gz"

    spatial_obj = SPATIALTYPE2CLASS[spatial_format]()

    total_steps = 3 if perform_primary_analysis else 2
    step_counter = 1

    try:
        spatial_obj.process_file(filepath.as_posix(), extract_dir=staging_area, organism_id=organism_id)
    except Exception as e:
        status["status"] = "error"
        status["message"] = f"Error in uploading spatial file: {e}"
        write_status(status_file, status)
        return {"success": 0, "message": status["message"]}

    status["progress"] = int((step_counter / total_steps) * 100)
    write_status(status_file, status)

    output_path = staging_area / f"{share_uid}.zarr"
    # Remove existing Zarr store if present; a safeguard in case a prior attempt
    # failed after the store was partially written.
    if output_path.exists():
        shutil.rmtree(output_path)

    status["message"] = "Writing Zarr store"
    write_status(status_file, status)

    try:
        spatial_obj.write_to_zarr(filepath=output_path)
    except Exception as e:
        status["status"] = "error"
        status["message"] = f"Error writing Zarr store: {e}"
        write_status(status_file, status)
        return {"success": 0, "message": status["message"]}

    step_counter += 1
    status["progress"] = int((step_counter / total_steps) * 100)
    status["status"] = "complete"
    status["message"] = "Dataset processed successfully."
    write_status(status_file, status)

    return {"success": 1, "message": status["message"]}
