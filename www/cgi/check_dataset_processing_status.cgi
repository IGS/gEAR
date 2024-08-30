#!/opt/bin/python3

"""
For a given session_id, and dataset share ID being uploaded,
this returns the status of the dataset processing.  It does
a bit extra, pulling the process ID from the JSON file and
making sure that process is in fact still running.

Structure returned:

{
    "process_id": 1234,
    "status": "processing",
    "message": "Processing the dataset.  This may take a while.",
    "progress": 0
}

Where status can be 'extracting', 'processing', 'error', or 'complete'.
"""

import cgi, json
import os, sys

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_uid = form.getvalue('share_uid')
    
    user_upload_file_base = os.path.join("../uploads/files", session_id, share_uid)
    status_file = os.path.join(user_upload_file_base, 'status.json')

    if not os.path.exists(status_file):
        status = {
            "process_id": -1,
            "status": "error",
            "message": "No status file found.  Please upload a dataset first.",
            "progress": 0
        }
        print(f"ERROR: Failed to find status file: {status_file}", file=sys.stderr)
        print(json.dumps(status))
        return
    
    with open(status_file, 'r') as f:
        status = json.load(f)

    state = status.get('status', '')

    if state in 'complete':
        print(json.dumps(status))
        return
    
    if state == 'processing':
        # Check if the process is still running
        process_id = status.get('process_id', -1)
        if process_id > 0:
            # TODO: check that the process is the correct name too
            if os.system(f'ps -p {process_id} > /dev/null') != 0:
                status['status'] = 'error'
                status['message'] = 'The processing step failed. Please contact the gEAR team.'
                status['progress'] = 0

        print(json.dumps(status))
        return

    
    print(json.dumps(status))


if __name__ == '__main__':
    main()
