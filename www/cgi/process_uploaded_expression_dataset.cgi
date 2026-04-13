#!/opt/bin/python3

"""
Process the uploaded expression dataset by publishing to RabbitMQ queue.

This CGI endpoint queues dataset processing jobs for asynchronous handling
by the anndata_upload_consumer worker.
"""

import cgi
import json
import os
import sys
import uuid
from pathlib import Path

# Redirect stdout to suppress debug output from dependencies
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = Path(__file__).resolve().parents[2] / 'lib'
sys.path.append(str(lib_path))

import geardb
from gear.anndata_processor import write_status
from gear.serverconfig import ServerConfig
from gear.spatialhandler import SPATIALTYPE2CLASS

servercfg = ServerConfig().parse()
user_upload_file_base = '../uploads/files'


def main() -> dict:
    """Queue expression dataset processing job."""
    result = {'success': 0, "message": ""}

    form = cgi.FieldStorage()
    share_uid = form.getvalue('share_uid')
    session_id = form.getvalue('session_id')
    dataset_format = form.getvalue('dataset_format')
    spatial_format = form.getvalue('spatial_format')  # May be None

    # Validate required parameters
    if not all([share_uid, session_id, dataset_format]):
        result['message'] = 'Missing one or more required parameters.'
        return result

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        return result

    # Validate dataset format
    supported_adata_formats = ['h5ad', 'mex_3tab', 'excel', 'rdata']
    spatial_formats = list(SPATIALTYPE2CLASS.keys())

    if dataset_format not in supported_adata_formats and dataset_format != 'spatial':
        result['message'] = f'Unsupported dataset format: {dataset_format}'
        return result

    if dataset_format == "spatial":
        if spatial_format not in spatial_formats:
            result['message'] = f'Invalid spatial format: {spatial_format}'
            return result

    # Locate dataset upload directory
    dataset_upload_dir = Path(user_upload_file_base) / session_id / share_uid

    if not dataset_upload_dir.is_dir():
        result['message'] = 'Dataset/directory not found.'
        return result

    # Initialize status file
    status = {
        "process_id": None,
        "status": "queued",
        "message": "Job queued for processing.",
        "progress": 0,
    }
    status_file = dataset_upload_dir / 'status.json'
    write_status(status_file, status)

    # Load and update metadata
    metadata_file = dataset_upload_dir / 'metadata.json'
    if not metadata_file.is_file():
        result['message'] = "No metadata JSON file found."
        return result

    try:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        dataset_uid = metadata.get('dataset_uid', '')
        dataset_type = metadata.get('dataset_type', '')

        # Determine if primary analysis should be performed
        perform_primary_analysis = (
            dataset_type in ['single-cell-rnaseq', 'spatial'] and
            dataset_format != 'spatial'  # Exclude spatial from anndata processor
        )

        metadata["dataset_format"] = dataset_format
        metadata["perform_primary_analysis"] = perform_primary_analysis

        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=4)

    except (json.JSONDecodeError, IOError) as e:
        result['message'] = f"Error reading metadata: {str(e)}"
        return result

    # Queue the job (skips spatial format)
    if dataset_format == 'spatial':
        process_spatial(dataset_upload_dir, share_uid, spatial_format, metadata["perform_primary_analysis"])
        return result

    try:
        import gearqueue

        publisher = gearqueue.AsyncConnection(
            host=servercfg["dataset_uploader"]["queue_host"],
            publisher_or_consumer="publisher",
            queue_name="anndata_upload_jobs",
        )

        message = {
            "job_id": str(uuid.uuid4()),
            "share_uid": share_uid,
            "dataset_uid": dataset_uid,
            "dataset_format": dataset_format,
            "perform_primary_analysis": perform_primary_analysis,
        }

        publisher.publish(json.dumps(message))

        result['success'] = 1
        result['message'] = "Job queued for processing."

    except Exception as e:
        result['message'] = f"Failed to queue job: {str(e)}"

    return result

def process_spatial(upload_dir: Path, share_uid: str, spatial_format: str, perform_primary_analysis: bool) -> None:
    """
    Processes a spatial transcriptomics dataset uploaded to a specified directory.

    This function handles the reading and conversion of spatial data files using a handler
    class determined by the spatial_format. It expects a metadata.json file in the upload
    directory to extract the sample's taxonomic ID, which is then used to retrieve the
    organism ID. The function reads the spatial data archive, processes it, and writes the
    output in Zarr format. Status updates and errors are logged using the write_status function.

    Args:
        upload_dir (str): The directory where the uploaded files are located.
        spatial_format (str): The format of the spatial data, used to select the appropriate handler.

    Raises:
        Writes error status if the metadata file is missing or if reading/converting the spatial file fails.
    """

    status = {
        "process_id": None,
        "status": "processing",
        "message": "Initializing dataset processing.",
        "progress": 0,
    }
    status_file = upload_dir / 'status.json'
    write_status(status_file, status)

    spatial_obj = SPATIALTYPE2CLASS[spatial_format]()   # instantiate the appropriate handler class
    metadata_file = upload_dir / 'metadata.json'
    if not metadata_file.is_file():
        status["status"] = "error"
        status["message"] = "No metadata JSON file found."
        write_status(status_file, status)
        return

    # get organism_id by converting sample_taxid(needed for some but not all spatial handlers)
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)

    sample_taxid = metadata.get("sample_taxid", None)
    organism_id=geardb.get_organism_id_by_taxon_id(sample_taxid)
    filepath = upload_dir / f"{share_uid}.tar.gz"

    try:
        spatial_obj.process_file(filepath.as_posix(), extract_dir=upload_dir, organism_id=organism_id)
    except Exception as e:
        import traceback
        traceback.print_exc()
        status["status"] = "error"
        status["message"] = f"Error in uploading spatial file: {e}"
        write_status(status_file, status)
        return

    output_filename = f"{share_uid}.zarr"
    output_path = upload_dir / output_filename
    # Remove existing Zarr store if it exists
    # This is a safeguard; it shouldn't exist at this point, unless there was a failure post-writing
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)

    total_steps = 3 if perform_primary_analysis else 2
    step_counter = 1
    status["progress"] = int((step_counter / total_steps) * 100)
    write_status(status_file, status)

    spatial_obj.write_to_zarr(filepath=output_path)

    step_counter += 1
    status["progress"] = int((step_counter / total_steps) * 100)
    write_status(status_file, status)


if __name__ == '__main__':
    result = main()
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result), flush=True)