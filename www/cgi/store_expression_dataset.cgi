#!/opt/bin/python3

"""
Used by the expression uploader, this stores the actual dataset from the form
and saves it to a file for processing.

Writes a file at: ../uploads/files/<session_id>_<share_uid>.<ext>
"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_uid = form.getvalue('dataset_uid')
    share_uid = form.getvalue('share_uid')
    dataset_format = form.getvalue('dataset_format')
    user_upload_file_base = '../uploads/files'

    user = geardb.get_user_from_session_id(session_id)
    result = {'success': 0}

    filename = form['dataset_file'].filename
    file_extension = filename.split('.')[-1]

    dataset_filename = os.path.join(user_upload_file_base, session_id + '_' + share_uid + '.' + file_extension)

    print("DEBUG: Saving dataset to: " + dataset_filename, file=sys.stderr)

    try:
        with open(dataset_filename, 'wb') as f:
            f.write(form['dataset_file'].file.read())
        result['success'] = 1
    except Exception as e:
        raise Exception("DEBUG: Error saving dataset file: " + str(e))

    print(json.dumps(result))

if __name__ == '__main__':
    main()
