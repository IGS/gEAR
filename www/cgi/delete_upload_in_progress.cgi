#!/opt/bin/python3

"""
This script deletes one of directories representing a user's upload in progress.
"""

import cgi
import json
import os, sys
import shutil
from pathlib import Path

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from werkzeug.utils import secure_filename

share_uid = None
dataset_id = None
session_id = None
user_upload_file_base = '../uploads/files'

def main():
    print('Content-Type: application/json\n\n', flush=True)
    result = {'success':0, 'message':''}
    global share_uid
    global session_id

    form = cgi.FieldStorage()
    share_uid = form.getfirst('share_uid')
    session_id = form.getfirst('session_id')
    dataset_id = form.getfirst('dataset_id')

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print(json.dumps(result))
        return

    safe_share_uid = secure_filename(share_uid or '')
    safe_session_id = secure_filename(session_id or '')
    if not safe_share_uid or safe_share_uid != share_uid:
        result['message'] = 'Invalid share_uid: ' + str(share_uid)
        print(json.dumps(result))
        return
    if not safe_session_id or safe_session_id != session_id:
        result['message'] = 'Invalid session_id.'
        print(json.dumps(result))
        return

    uploads_base = Path(user_upload_file_base).resolve()
    user_upload_file_path = (uploads_base / safe_session_id / safe_share_uid).resolve()

    # Defense in depth: confirm the resolved path is still contained within the
    # uploads base directory before performing a recursive delete.
    if not user_upload_file_path.is_relative_to(uploads_base):
        result['message'] = 'Invalid upload path.'
        print(json.dumps(result))
        return

    if not user_upload_file_path.exists():
        result['message'] = 'Upload directory not found: ' + str(user_upload_file_path)
        print(json.dumps(result))
        return

    try:
        # recursively delete the directory
        shutil.rmtree(user_upload_file_path)

    except Exception as e:
        result['message'] = 'Error deleting file: ' + str(e)
        print(json.dumps(result))
        return

    result['success'] = 1
    result['message'] = 'File deleted successfully.'
    print(json.dumps(result))


if __name__ == '__main__':
    main()
