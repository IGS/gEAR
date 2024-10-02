#!/opt/bin/python3

"""
This script deletes one of directories representing a user's upload in progress.
"""

import cgi
import json
import os, sys
import shutil

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
    share_uid = form.getvalue('share_uid')
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print(json.dumps(result))
        return
    
    user_upload_file_path = os.path.join(user_upload_file_base, session_id, share_uid)

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
