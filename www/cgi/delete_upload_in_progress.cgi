#!/opt/bin/python3

"""
This script deletes one of directories representing a user's upload in progress.
"""

import cgi
import json
import os, sys
import shutil
import re

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

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

    # Make sure the final directory looks like a share_uid (8 alphanumeric characters)
    if not re.match(r'^[a-zA-Z0-9]{8}$', share_uid or ''):
        result['message'] = 'Invalid share_uid: ' + str(share_uid)
        print(json.dumps(result))
        return

    # session_id is expected to be an existing, server-issued session identifier;
    # enforce a safe filesystem-component format before it's used in a path.
    if not re.match(r'^[a-zA-Z0-9-]+$', session_id or ''):
        result['message'] = 'Invalid session_id.'
        print(json.dumps(result))
        return

    user_upload_file_path = os.path.abspath(os.path.join(user_upload_file_base, session_id, share_uid))
    uploads_base = os.path.abspath(user_upload_file_base)

    # Defense in depth: confirm the resolved path is still contained within the
    # uploads base directory before performing a recursive delete.
    if os.path.commonpath([user_upload_file_path, uploads_base]) != uploads_base:
        result['message'] = 'Invalid upload path.'
        print(json.dumps(result))
        return

    if not os.path.exists(user_upload_file_path):
        result['message'] = 'Upload directory not found: ' + user_upload_file_path
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
