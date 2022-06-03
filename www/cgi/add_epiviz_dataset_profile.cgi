#!/opt/bin/python3

"""
- Receives dataset upload data from epiviz_panel_designer.js.
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

    result = {'success':0}
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_uid = html.escape(form.getvalue('dataset_uid'))
    share_uid = html.escape(form.getvalue('share_uid'))
    dataset_name = form.getvalue('title')
    dataset_description = form.getvalue('description')
    dataset_config = form.getvalue("epiviz_config")
    dataset_access = form.getvalue("access")

    # Must have a gEAR account to upload datasets
    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        # this shouldn't happen, but let's check here just in case.
        # raise Exception("ERROR: must have a user ID to store a dataset")
        result['message'] = 'User ID not found. Please log in to continue.'

    else:
        # add file to both dataset and dataset_epiviz
        cnx = geardb.Connection()
        cursor = cnx.get_cursor()

        # add to dataset
        dataset_sql = """
            INSERT INTO dataset (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, ldesc, date_added, dtype, schematic_image,
                share_id, math_default, load_status, has_h5ad, platform_id, instrument_model, library_selection, library_source, library_strategy, contact_email, contact_institute, contact_name, annotation_source, annotation_release, plot_default)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, NOW(), %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """

        # Insert dataset info to database
        cursor.execute(dataset_sql, (dataset_uid, user.id, dataset_name, 1, None,
                                         None, dataset_access, dataset_description, "epiviz", "img/schematic_epiviz-chipseq.png",
                                         share_uid, "raw", "completed", 0, None,
                                         None, None, None, None, None,
                                         None, None, None, None, None,))
        cnx.commit()

        #  insert into dataset_display
        dataset_display_sql = """
            INSERT INTO dataset_display (dataset_id, user_id, label, plot_type, plotly_config)
            VALUES (%s, %s, %s, %s, %s)
        """

        # Insert dataset_epiviz info to database
        cursor.execute(dataset_display_sql, (dataset_uid, user.id, "Epiviz", "epiviz", dataset_config,))

        cnx.commit()

        #  set preference
        dataset_preference_sql = """
            INSERT INTO dataset_preference (user_id, dataset_id, display_id)
            VALUES (%s, %s, %s)
        """

        # Insert dataset_epiviz info to database
        cursor.execute(dataset_preference_sql, (user.id, dataset_uid, cursor.lastrowid,))

        cnx.commit()

        result['success'] = 1

    print(json.dumps(result))

if __name__ == '__main__':
    main()
