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
from gear.utils import (
    flag_ambiguous_obs_columns,
    standardize_and_sanitize_obs,
)


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
    output_path = staging_area / f"{share_uid}.zarr"

    spatial_obj = SPATIALTYPE2CLASS[spatial_format]()

    def _sanitize_and_flag_obs_columns() -> None:
        """
        Scan the final obs table for numeric columns that look like they may
        actually be categorical (e.g. replicate/slide numbers), and record
        them in metadata.json for the uploader's "review column types" step.

        Opens the file in backed mode -- obs is small regardless of dataset
        size, and X is never loaded or touched.
        """

        adata = spatial_obj.sdata.tables["table"]

        # Standarize some columns before proceeding to the ambiguous ones.
        adata.obs = standardize_and_sanitize_obs(adata.obs)

        questionable = flag_ambiguous_obs_columns(adata.obs)

        metadata_file = staging_area / 'metadata.json'
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        metadata['questionable_obs_columns'] = questionable
        # Nothing flagged -- nothing for the user to review, so the uploader
        # can skip straight past that step.
        metadata['obs_dtype_reviewed'] = not bool(questionable)

        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=4)

    def _write_zarr():
        # Remove existing Zarr store if present; a safeguard in case a prior
        # attempt failed after the store was partially written.
        if output_path.exists():
            shutil.rmtree(output_path)
        spatial_obj.write_to_zarr(filepath=output_path)

    # Each SpatialData-modifying stage of the pipeline gets its own status update,
    # rather than one opaque "processing" step covering everything from archive
    # extraction through embeddings. (message, error_label, action)
    steps = [
        (
            "Reading and parsing spatial data archive...",
            "reading spatial data archive",
            lambda: spatial_obj.process_file(filepath.as_posix(), extract_dir=staging_area, organism_id=organism_id),
        ),
        (
            "Subsetting spatial data...",
            "subsetting spatial data",
            spatial_obj.subset_sdata,
        ),
        (
            "Scaling and aligning spatial coordinates...",
            "scaling/aligning spatial coordinates",
            spatial_obj.scale_and_translate_sdata,
        ),
        (
            "Merging spatial coordinates into observations...",
            "merging spatial coordinates into observations",
            spatial_obj.merge_centroids_with_obs,
        ),
        (
            "Computing QC metrics and embeddings...",
            "computing QC metrics and embeddings",
            spatial_obj.compute_qc_and_embeddings,
        ),
        (
            "Flagging ambiguous observation types...",
            "flagging ambiguous observation types",
            _sanitize_and_flag_obs_columns,
        ),
        (
            "Writing Zarr store...",
            "writing Zarr store",
            _write_zarr,
        ),
    ]

    total_steps = len(steps) + (1 if perform_primary_analysis else 0)

    for step_index, (message, error_label, action) in enumerate(steps, start=1):
        status["message"] = message
        status["progress"] = int(((step_index - 1) / total_steps) * 100)
        write_status(status_file, status)

        try:
            action()
        except MemoryError:
            # A bare MemoryError's str() is typically empty/unhelpful on its own -
            # give a specific, actionable message instead of falling through to the
            # generic branch below.
            status["status"] = "error"
            status["message"] = (
                f"This dataset needed more memory than is available on this server to "
                f"complete '{error_label}'. Please contact the gEAR team to resolve this "
                f"issue (share ID: {share_uid})."
            )
            write_status(status_file, status)
            return {"success": 0, "message": status["message"]}
        except Exception as e:
            status["status"] = "error"
            status["message"] = f"Error {error_label}: {e}"
            write_status(status_file, status)
            return {"success": 0, "message": status["message"]}

    status["status"] = "complete"
    status["progress"] = 100
    status["message"] = "Dataset processed successfully."
    write_status(status_file, status)

    return {"success": 1, "message": status["message"]}

