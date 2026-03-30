import json
import os
from pathlib import Path

import geardb
from flask import request
from flask_restful import Resource

user_upload_file_base = Path(__file__).resolve().parents[3] / 'www' / 'uploads' / 'files'


class DatasetProcessingStatus(Resource):
    """
    Unified endpoint for checking processing status of any dataset type.
    Returns status.json contents for expression datasets, spatial datasets, or trackhubs.
    """
    def post(self, share_uid):
        """
        Get processing status for a dataset upload.

        Query params:
            - share_uid: The share UID for the upload

        Returns:
            - status.json contents on success
            - 404 if status file not found
            - 401 if not authenticated
        """
        session_id = request.cookies.get('gear_session_id', "")

        status_data = {
            "success": False,
            "message": "",
            "status": "error",
            "progress": 0
        }

        user = geardb.get_user_from_session_id(session_id)
        if not user:
            status_data["message"] = "Invalid session"
            return status_data, 401

        staging_area = user_upload_file_base / session_id / share_uid
        status_file = staging_area / "status.json"

        if not status_file.is_file():
            status_data["message"] = "Job status file not found"
            return status_data, 404

        with open(status_file, 'r') as f:
            status_data = json.load(f)

            state = status_data.get('status', '')

            if state in 'complete':
                status_data['progress'] = 100
                return status_data

            if state == 'processing':
                # Check if the process is still running
                process_id = status_data.get('process_id', -1)
                job_id = status_data.get('job_id', "")
                if process_id > 0:
                    # If job_id is present, then it's a gosling job, so we are fine
                    if job_id:
                        return status_data, 200

                    # TODO: check that the process is the correct name too
                    if os.system(f'ps -p {process_id} > /dev/null') != 0:
                        status_data['status'] = 'error'
                        status_data['message'] = 'The processing step failed. Please contact the gEAR team.'
                        status_data['progress'] = 0

        return status_data, 200
