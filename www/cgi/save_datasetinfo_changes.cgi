#!/opt/bin/python3

"""
Changes the information of a dataset via .is-editable elements in the Dataset Explorer

This script first checks if the user owns the dataset, then proceeds with the info change
if they own it.

If the user does not own the dataset, an error is returned stating that.
Successful access_level change returns the dataset JSON that was changed.

Requires:
1) Session id - which contains user_id
2) Dataset id
3) One of the following:
    - Title
    - Description
    - Pubmed ID
    - GEO ID
    - Public/private visibility

"""

import os
import cgi
import json
import os
import sys
import re
import shutil

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')
    visibility = form.getvalue('visibility')
    is_downloadable = form.getvalue('is_downloadable')
    title = form.getvalue('title')
    pubmed_id = form.getvalue('pubmed_id')
    geo_id = form.getvalue('geo_id')
    ldesc = form.getvalue('ldesc')

    user = geardb.get_user_from_session_id(session_id)
    dataset = geardb.get_dataset_by_id(d_id=dataset_id)
    if dataset is None:
        result = {'error': 'Dataset not found', 'success': 0}
        print(json.dumps(result))
        return

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, user.id, dataset.id)

    if owns_dataset == True:
        # see what has changed and execute updates to the DB
        # ? SAdkins - Why are we checking for differences? Can't we just update regardless, or are we trying to reduce transactions?
        if dataset.is_public != visibility:
            dataset.save_change('is_public', visibility)

        if dataset.is_downloadable != is_downloadable:
            dataset.save_change("is_downloadable", is_downloadable)

        if dataset.title != title:
            dataset.save_change('title', title)

        if dataset.pubmed_id != pubmed_id:
            dataset.save_change('pubmed_id', pubmed_id)

        if dataset.geo_id != geo_id:
            dataset.save_change('geo_id', geo_id)

        if dataset.ldesc != ldesc:
            dataset.save_change('ldesc', ldesc)

        result = { 'dataset': dataset, 'success': 1 }

        print(json.dumps(result))

    else:
        result = { 'error':[] }

        error = "Not able to change dataset's information. You do not own this dataset."
        result['error'] = error
        result['owns_dataset'] = owns_dataset
        result['success'] = 0

        print(json.dumps(result))


def check_dataset_ownership(cursor, current_user_id, dataset_id):
    qry = """
       SELECT d.id, d.owner_id
       FROM dataset d
       WHERE d.id = %s
    """
    cursor.execute(qry, (dataset_id,))

    # default: Assume user does not own dataset
    user_owns_dataset = False

    for row in cursor:
        # Change access if user owns the dataset
        if row[1] == current_user_id:
            user_owns_dataset = True

        # Return a statement that the user does not own the dataset (have permission)
        else:
        	user_owns_dataset = False

    return user_owns_dataset

if __name__ == '__main__':
    main()
