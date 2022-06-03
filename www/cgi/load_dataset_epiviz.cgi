#!/opt/bin/python3

"""
-Receives dataset upload data from upload_epigenetic_dataset.js.
for each file the user uploads, creates a database entry into 
dataset & dataset_epiviz

"""


import cgi
import html
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

    user_upload_file_base = '../datasets_epigenetic/uploads/files'
    user_upload_dest_base = '../datasets_epigenetic'

    result = {'success':0}
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_uid = html.escape(form.getvalue('dataset_uid'))
    share_uid = html.escape(form.getvalue('share_uid'))

    file_name = None
    file_url = None
    is_url = False
    if "file_name" in form:
        file_name = form.getvalue('file_name')
        file_type = file_name.split(".")[-1]
    elif "file_url" in form:
        file_url = form.getvalue('file_url')
        file_type = file_url.split(".")[-1]
        is_url = True

    file_title = form.getvalue('file_title')
    file_annoatation = form.getvalue('file_annoatation')
    file_organism = form.getvalue('file_organism')
    file_access = form.getvalue('file_access')
    default_plot_type = "Epiviz"

    file_location = None

    if not is_url:
        # modify this later so that file urls are handled separately
        file_source_path = user_upload_file_base + '/' + file_name
        file_dest_path = user_upload_dest_base + '/' + dataset_uid + "." + file_type
        file_location = os.path.abspath(file_dest_path)
    else:
        file_location = file_url
    
    # Must have a gEAR account to upload datasets
    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        # this shouldn't happen, but let's check here just in case.
        # raise Exception("ERROR: must have a user ID to store a dataset")
        result['message'] = 'User ID not found. Please log in to continue.'

    else:
        
        if not is_url:
            #  move file to dest
            shutil.move(file_source_path, file_dest_path)

        # add file to both dataset and dataset_epiviz
        cnx = geardb.Connection()
        cursor = cnx.get_cursor()

        dataset_epiviz_sql = """
            INSERT INTO dataset_epiviz (id, owner_id, annotation, type, url, title, is_public, share_id, organism)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
        """

        # Insert dataset_epiviz info to database
        cursor.execute(dataset_epiviz_sql, (dataset_uid, user.id, file_annoatation, file_type, file_location, file_title, file_access, share_uid, file_organism, ))
        cnx.commit()

        result['success'] = 1

    print(json.dumps(result))

if __name__ == '__main__':
    main()
