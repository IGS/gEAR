#!/opt/bin/python3

"""
At this point we should have a directory with at least the following:

- JSON file with metadata
- An H5AD file ready to go
- A status.json file

Then, depending on the upload format, we also should have

- If 3tab:
    - A tarball uploaded by the user
- If Excel:
    - An Excel file uploaded by the user

This script does the following:

- Loads the JSON metadata file and stores it in MySQL
- Migrates the H5AD file to the proper directory
- Migrates the original uploaded file so it can be downloaded by users

Returns a status of these steps as JSON data like this:


result = {
    "success": 1,
    "metadata_loaded": 1,
    "h5ad_migrated": 1,
    "userdata_migrated": 1,    
    "message": "All steps completed successfully."
}
"""

import cgi
import json
import os, sys
import shutil

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from gear.metadata import Metadata

user_upload_file_base = '../uploads/files'
dataset_final_dir = '../datasets'

result = {
    "success": 0,
    "metadata_loaded": 0,
    "h5ad_migrated": 0,
    "userdata_migrated": 0,
    "message": ""
}

def main():
    print('Content-Type: application/json\n\n', flush=True)

    form = cgi.FieldStorage()
    share_uid = form.getvalue('share_uid')
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_uid')
    dataset_format = form.getvalue('dataset_format')
    dataset_visibility = form.getvalue('dataset_visibility')

    if dataset_visibility == 'private':
        is_public = 0
    elif dataset_visibility == 'public':
        is_public = 1
    else:
        result['message'] = 'Invalid dataset visibility.'
        print(json.dumps(result))
        sys.exit(0)

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print(json.dumps(result))
        sys.exit(0)

    dataset_upload_dir = os.path.join(user_upload_file_base, session_id, share_uid)

    # if the upload directory doesn't exist, we can't process the dataset
    if not os.path.exists(dataset_upload_dir):
        result['message'] = 'Dataset/directory not found.'
        print(json.dumps(result))
        sys.exit(0)

    # Load the metadata
    metadata_file = os.path.join(dataset_upload_dir, 'metadata.json')
    if not os.path.exists(metadata_file):
        result['message'] = 'Metadata file not found.'
        print(json.dumps(result))
        sys.exit(0)

    with open(metadata_file, 'r') as f:
        metadata = json.load(f)

    # Load the metadata into the database
    metadata = Metadata(file_path=metadata_file)
    try:
        metadata.save_to_mysql(status='completed', is_public=is_public)
        result['metadata_loaded'] = 1
    except Exception as e:
        result['message'] = 'Error saving metadata to MySQL: {}'.format(str(e))
        print(json.dumps(result))
        sys.exit(0)

    # migrate the H5AD file
    h5ad_file = os.path.join(dataset_upload_dir, f'{share_uid}.h5ad')
    h5ad_dest = os.path.join(dataset_final_dir, f'{dataset_id}.h5ad')
    if not os.path.exists(h5ad_file):
        result['message'] = 'H5AD file not found: {}'.format(h5ad_file)
        print(json.dumps(result))
        sys.exit(0)

    try:
        os.rename(h5ad_file, h5ad_dest)
        result['h5ad_migrated'] = 1
    except Exception as e:
        result['message'] = 'Error migrating H5AD file: {}'.format(str(e))
        print(json.dumps(result))
        sys.exit(0)
    
    if dataset_format == '3tab' or dataset_format == 'mex':
        # migrate the tarball
        tarball_file = os.path.join(dataset_upload_dir, f'{share_uid}.tar.gz')
        tarball_dest = os.path.join(dataset_final_dir, f'{dataset_id}.tar.gz')

        try:
            os.rename(tarball_file, tarball_dest)
            result['userdata_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating tarball file: {}'.format(str(e))
            print(json.dumps(result))
            sys.exit(0)
    elif dataset_format == 'excel':
        # migrate the Excel file
        excel_file = os.path.join(dataset_upload_dir, f'{share_uid}.xlsx')
        excel_dest = os.path.join(dataset_final_dir, f'{dataset_id}.xlsx')

        try:
            os.rename(excel_file, excel_dest)
            result['userdata_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating Excel file: {}'.format(str(e))
            print(json.dumps(result))
            sys.exit(0)

    # if we made it this far, all is well, so return success
    result['success'] = 1
    result['message'] = 'All steps completed successfully.'
    print(json.dumps(result))

    # now delete the entire upload directory
    shutil.rmtree(dataset_upload_dir)

if __name__ == '__main__':
    main()
