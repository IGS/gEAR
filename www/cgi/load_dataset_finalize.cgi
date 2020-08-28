#!/opt/bin/python3

"""
-Receives dataset upload data from upload_dataset.js.
-Reads metadata file as pandas dataframe, and populates data from GEO if GEO GSE
    ID is included.
-Stores metadata and dataset details in the mysql 'dataset' table with
    'load_status' set to 'completed'.
-Writes metadata dataframe as JSON to www/datasets
-Renames and moves original metadata file to www/datasets

"""


import cgi
from datetime import datetime
import json
import numpy as np
import os, sys
import re
import shutil


lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.metadata import Metadata

sys.path.append('../..')

def main():
    print('Content-Type: application/json\n\n', flush=True)

    user_upload_file_base = '../uploads/files'
    user_upload_dest_base = '../datasets'

    result = {'success':0}

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_uid = cgi.escape(form.getvalue('dataset_uid'))
    share_uid = cgi.escape(form.getvalue('share_uid'))
    metadata_filename = form.getvalue('uploaded_metadata_name')
    expression_filename = form.getvalue('uploaded_expression_name')
    default_plot_type = form.getvalue('default_plot_type')

    expression_source_path = user_upload_file_base + '/' + expression_filename.replace('.', '_') + '.h5ad'
    expression_dest_path = user_upload_dest_base + '/' + dataset_uid + '.h5ad'
    
    # Metadata file paths to: original, renamed original dest, output json dest
    metadata_user_path = user_upload_file_base + '/' + metadata_filename
    metadata_user_path_renamed = user_upload_dest_base + '/' + dataset_uid + '.xlsx'
    metadata_dest_path = user_upload_dest_base + '/' + dataset_uid + '.json'

    # Must have a gEAR account to upload datasets
    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        # this shouldn't happen, but let's check here just in case.
        # raise Exception("ERROR: must have a user ID to store a dataset")
        result['message'] = 'User ID not found. Please log in to continue.'

    else:
        # Read metadata file & try to populate from GEO
        metadata = Metadata(file_path=metadata_user_path)
        metadata.populate_from_geo()

        # Add IDs to metadata content
        metadata.add_field_value('dataset_uid', dataset_uid)
        metadata.add_field_value('share_uid', share_uid)
        metadata.add_field_value('owner_id', user.id)
        metadata.add_field_value('default_plot_type', default_plot_type)
        metadata.add_field_value('schematic_image', None)

        # Write H5AD file
        shutil.move(expression_source_path, expression_dest_path)

        # Write metadata to JSON
        metadata.write_json(file_path=metadata_dest_path)

        source_user_file = user_upload_file_base + '/' + expression_filename
        
        if expression_filename.endswith('.xlsx'):
            dest_user_file = user_upload_dest_base + '/' + dataset_uid + '.xlsx'
        elif expression_filename.endswith('.tar'):
            dest_user_file = user_upload_dest_base + '/' + dataset_uid + '.tar'
        elif expression_filename.endswith('.gz'):
            # Can't search for .tar.gz because many duplicate upload steps create files like: GSE11347.tar (13).gz
            dest_user_file = user_upload_dest_base + '/' + dataset_uid + '.tar.gz'

        shutil.move(source_user_file, dest_user_file)

        # Write metadata to gEAR MySQL. Set load_status = 'completed'
        metadata.save_to_mysql(status='completed')
        result['success'] = 1

    print(json.dumps(result))


if __name__ == '__main__':
    main()
