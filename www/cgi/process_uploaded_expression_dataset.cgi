#!/opt/bin/python3

"""
Process the uploaded expression dataset by publishing to RabbitMQ queue.

This CGI endpoint queues dataset processing jobs for asynchronous handling
by the anndata_upload_consumer worker.
"""

import cgi
import configparser
import json
import os
import sys
import uuid
from pathlib import Path

# Redirect stdout to suppress debug output from dependencies
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

gear_root = Path(__file__).resolve().parents[2]
lib_path = gear_root / 'lib'
sys.path.append(str(lib_path))

_config = configparser.ConfigParser()
_config.read(gear_root / 'gear.ini')

import geardb
from gear.anndata_processor import process_anndata_synchronously, write_status
from gear.spatial_processor import process_spatial_synchronously
from gear.spatialhandler import SPATIALTYPE2CLASS

user_upload_file_base = '../uploads/files'

class QueueDisabledError(Exception):
    """Custom exception to indicate that the queue is disabled in configuration."""
    pass

def main() -> tuple:
    """Queue expression dataset processing job."""
    result = {'success': 0, "message": ""}

    form = cgi.FieldStorage()
    share_uid = form.getfirst('share_uid')
    session_id = form.getfirst('session_id')
    dataset_format = form.getfirst('dataset_format')
    spatial_format = form.getfirst('spatial_format')  # May be None

    # Validate required parameters
    if share_uid is None:
        result['message'] = 'Share UID is required.'
        return result, 400

    if session_id is None:
        result['message'] = 'Session ID is required.'
        return result, 400

    if dataset_format is None:
        result['message'] = 'Dataset format is required.'
        return result, 400

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        return result, 401

    # Validate dataset format
    supported_adata_formats = ['h5ad', 'mex_3tab', 'excel', 'rds']
    spatial_formats = list(SPATIALTYPE2CLASS.keys())

    if dataset_format not in supported_adata_formats and dataset_format != 'spatial':
        result['message'] = f'Unsupported dataset format: {dataset_format}'
        return result, 400

    if dataset_format == "spatial":
        if spatial_format not in spatial_formats:
            result['message'] = f'Invalid spatial format: {spatial_format}'
            return result, 400

    # Locate dataset upload directory
    dataset_upload_dir = Path(user_upload_file_base) / session_id / share_uid

    if not dataset_upload_dir.is_dir():
        result['message'] = 'Dataset/directory not found.'
        return result, 404

    job_id = str(uuid.uuid4())

    # Initialize status file
    status = {
        "job_id": job_id,
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
        return result, 400

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
        return result, 500

    if dataset_format == 'spatial':
        if spatial_format is None:
            result['message'] = 'Spatial format is required for spatial datasets.'
            return result, 400

        try:
            queue_spatial_job(job_id, share_uid, spatial_format, metadata["perform_primary_analysis"])
            result["success"] = True
            result["message"] = "Dataset upload processing job queued"
            return result, 202  # Accepted
        except QueueDisabledError:
            result_sync = process_spatial_synchronously(
                job_id=job_id,
                share_uid=share_uid,
                staging_area=dataset_upload_dir,
                status_file=status_file,
                spatial_format=spatial_format,
                perform_primary_analysis=metadata["perform_primary_analysis"],
            )
            result["success"] = result_sync["success"]
            result["message"] = result_sync["message"]
            return result, 200 if result["success"] else 500
        except Exception as e:
            result["message"] = f"Error processing spatial dataset: {str(e)}"
            print(f"SpatialUpload error: {str(e)}", file=sys.stderr)
            return result, 500

    # Queue the job
    try:
        queue_anndata_job(job_id, share_uid, dataset_uid, dataset_format, metadata["perform_primary_analysis"])
        result["success"] = True
        result["message"] = "Dataset upload processing job queued"
        return result, 202  # Accepted
    except QueueDisabledError:

        result_sync = process_anndata_synchronously(
            job_id=job_id,
            share_uid=share_uid,
            staging_area=dataset_upload_dir,
            status_file=status_file,
            dataset_uid=dataset_uid,
            dataset_format=dataset_format,
            perform_primary_analysis=metadata["perform_primary_analysis"],

        )

        result["success"] = result_sync["success"]
        result["message"] = result_sync["message"]
        return result, 200 if result["success"] else 500

    except Exception as e:
        result["message"] = f"Error processing track hub: {str(e)}"
        print(f"AnndataUpload error: {str(e)}", file=sys.stderr)
        return result, 500

def queue_spatial_job(job_id: str, share_uid: str, spatial_format: str, perform_primary_analysis: bool) -> None:
    """Queue spatial dataset processing job to RabbitMQ."""

    if not _config.getboolean('dataset_uploader', 'queue_enabled', fallback=False):
        print("Queue is disabled in configuration. Cannot queue spatial job. Falling back to synchronous processing.", file=sys.stderr)
        raise QueueDisabledError()

    import gearqueue
    host = _config["dataset_uploader"]["queue_host"]

    try:
        connection = gearqueue.Connection(
            host=host, publisher_or_consumer="publisher"
        )
    except Exception as e:
        print(f"Error connecting to RabbitMQ: {e}", file=sys.stderr)
        raise Exception(f"Error connecting to RabbitMQ: {e}")

    with connection:
        connection.open_channel()

        payload = {
            "job_id": job_id,
            "share_uid": share_uid,
            "spatial_format": spatial_format,
            "perform_primary_analysis": perform_primary_analysis,
        }

        try:
            connection.publish(
                queue_name="spatial_upload_jobs",
                message=payload,  # method dumps JSON
            )
        except Exception as e:
            print(f"Error publishing message to RabbitMQ: {e}", file=sys.stderr)
            raise
    return

def queue_anndata_job(job_id: str, share_uid: str, dataset_uid: str, dataset_format: str, perform_primary_analysis: bool) -> None:
    """Queue anndata processing job to RabbitMQ."""

    # If queue is not enabled, return False
    if not _config.getboolean('dataset_uploader', 'queue_enabled', fallback=False):
        print("Queue is disabled in configuration. Cannot queue trackhub job. Falling back to synchronous processing.", file=sys.stderr)
        raise QueueDisabledError()

    import gearqueue
    host = _config["dataset_uploader"]["queue_host"]

    try:
        # Connect as a blocking RabbitMQ publisher
        connection = gearqueue.Connection(
            host=host, publisher_or_consumer="publisher"
        )
    except Exception as e:
        print(f"Error connecting to RabbitMQ: {e}", file=sys.stderr)
        raise Exception(f"Error connecting to RabbitMQ: {e}")

    with connection:
        connection.open_channel()

        payload = {
            "job_id": job_id,
            "share_uid": share_uid,
            "dataset_uid": dataset_uid,
            "dataset_format": dataset_format,
            "perform_primary_analysis": perform_primary_analysis,
        }

        try:
            connection.publish(
                queue_name="anndata_upload_jobs",
                message=payload,  # method dumps JSON
            )
        except Exception as e:
            print(f"Error publishing message to RabbitMQ: {e}", file=sys.stderr)
            raise
    return


if __name__ == '__main__':
    result, code = main()
    sys.stdout = original_stdout
    # Return JSON result and status code
    print(f"Status: {code}", flush=True)
    print('Content-Type: application/json', flush=True)
    print('', flush=True)  # blank line ends headers
    print(json.dumps(result), flush=True)