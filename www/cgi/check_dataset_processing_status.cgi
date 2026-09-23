#!/opt/bin/python3

"""
For a given session_id, and dataset share ID being uploaded,
this returns the status of the dataset processing.  It does
a bit extra, pulling the process ID from the JSON file and
making sure that process is in fact still running.

Structure returned:

{
    "job_id": 1234,
    "status": "processing",
    "message": "Processing the dataset.  This may take a while.",
    "progress": 0
}

Where status can be 'extracting', 'processing', 'error', or 'complete'.
"""

import cgi
import json
import subprocess
import sys
from pathlib import Path


def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    share_uid = form.getfirst('share_uid')

    if session_id is None:
        status = {
            "job_id": -1,
            "status": "error",
            "message": "No session_id provided.",
            "progress": 0
        }
        print(json.dumps(status))
        return status

    if share_uid is None:
        status = {
            "job_id": -1,
            "status": "error",
            "message": "No share_uid provided.",
            "progress": 0
        }
        print(json.dumps(status))
        return status

    user_upload_file_root = Path(__file__).resolve().parents[1] / 'uploads' / 'files'
    user_upload_file_base = user_upload_file_root / session_id / share_uid
    status_file = user_upload_file_base / 'status.json'

    if not status_file.is_file():
        status = {
            "job_id": -1,
            "status": "error",
            "message": "No status file found.  Please upload a dataset first.",
            "progress": 0
        }
        print(f"ERROR: Failed to find status file: {status_file}", file=sys.stderr)
        print(json.dumps(status))
        return status

    with open(status_file, 'r') as f:
        status = json.load(f)

    state = status.get('status', '')

    if state in 'complete':
        status['progress'] = 100
        return status

    if state == 'processing':
        job_id = status.get("job_id", -1)
        if job_id > 0:
            pass

        else:
            # Check if the process is still running.
            # Fallback for legacy uploads
            try:
                process_id = status.get('process_id', -1)
                if process_id == -1:
                    raise Exception("No process ID found in status data to check if process is running")

                result = subprocess.run(
                    ['ps', '-p', str(process_id)],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                if result.returncode != 0:
                    raise Exception(f"Process with ID {process_id} is not running")

            except Exception:
                status['status'] = 'error'
                status['message'] = 'No process nor job ID found for a job in processing state. Please contact the gEAR team.'
                status['progress'] = 0

    return status

if __name__ == '__main__':
    result = main()
    print(json.dumps(result))
