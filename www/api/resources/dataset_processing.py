import json
import os
import subprocess
import sys
from pathlib import Path

import geardb
from flask import request
from flask_restful import Resource

user_upload_file_base = Path(__file__).resolve().parents[3] / 'www' / 'uploads' / 'files'

def _check_process_running(status_data: dict) -> None:
    """
    Verify that the process is still running; mark as error if not.

    Args:
        status_data: The status dictionary to update if process is not running
    """
    process_id = status_data.get('process_id', -1)
    if process_id == -1:
        raise Exception("No process ID found in status data to check if process is running")

    result = subprocess.run(
        ['ps', '-p', str(process_id)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    if result.returncode != 0:
        raise Exception(f"Process with ID {process_id} is not running")


def _check_job_running(status_data: dict) -> None:
    """
    Verify that the job is still running; mark as error if not.

    Args:
        status_data: The status dictionary to update if job is not running
    """
    job_id = status_data.get('job_id', -1)
    if job_id == -1:
        raise Exception("No job ID found in status data to check if job is running")

def _error_response(message: str, status_code: int) -> tuple:
    """
    Generate a standardized error response.

    Args:
        message: Error message
        status_code: HTTP status code

    Returns:
        Tuple[dict, int]: (error_response, http_status_code)
    """
    return {
        'success': False,
        'message': message,
        'status': 'error',
        'progress': 0
    }, status_code

class DatasetProcessingStatus(Resource):
    """
    Unified endpoint for checking processing status of any dataset type.
    Returns status.json contents for expression datasets, spatial datasets, or trackhubs.
    """
    def post(self, share_uid):
        """
        Get processing status for a dataset upload.

        Args:
            share_uid: The share UID for the upload

        Returns:
            Tuple[dict, int]: (status_data, http_status_code)
        """

        session_id = request.cookies.get('gear_session_id', "")
        req = request.get_json() or {}
        dataset_format = req.get("dataset_format", "")  # e.g. "expression", "spatial", "gosling"

        user = geardb.get_user_from_session_id(session_id)
        if not user:
            return _error_response('Invalid session', 401)

        # Load status file
        staging_area = user_upload_file_base / session_id / share_uid
        status_file = staging_area / "status.json"
        if not status_file.is_file():
            return _error_response("Job status file not found", 404)

        with open(status_file, 'r') as f:
            status_data = json.load(f)

            current_status = status_data.get('status', '')

            if current_status in 'complete':
                status_data['progress'] = 100
            elif current_status == 'processing' and dataset_format not in ["gosling"]:
                # ? Remove the "gosling" check
                try:
                    _check_job_running(status_data)
                except Exception as e:
                    print(f"Error checking job status: {e}", file=sys.stderr)
                    # if it doesn't work, check for a process ID (legacy implementation)
                    try:
                        _check_process_running(status_data)
                    except Exception as e:
                        print(f"Error checking process status: {e}", file=sys.stderr)
                        # If both checks fail, we assume the process has failed
                        status_data['status'] = 'error'
                        status_data['message'] = 'The processing step failed. Please contact the gEAR team.'
                        status_data['progress'] = 0
                        return status_data, 400

        return status_data, 200
