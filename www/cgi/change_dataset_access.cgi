#!/opt/bin/python3

"""
Changes the access_level of a dataset via the 'change access' button in the dataset_manager

This script first checks if the user owns the dataset, then proceeds with the access_level change
if they own it.

If the user does not own the dataset, an error is returned stating that.
Successful access_level change returns the dataset JSON that was changed.

Requires:
1) Session id - which contains user_id
2) dataset id to be removed/deleted
3) Access level the dataset is being changed to

"""

import cgi
import json
import os
import sys

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
    access = form.getvalue('access')

    user = geardb.get_user_from_session_id(session_id)

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, user.id, dataset_id)

    if owns_dataset == True:
        result = { 'dataset':[] }

        change_dataset_access(cursor, user.id, dataset_id, access)
        cnx.commit()

        # Make the acceess_level change...
        result['dataset'].extend(get_dataset(cursor, user.id, dataset_id))

        cursor.close()
        cnx.close()

        print(json.dumps(result))

    else:
        result = { 'error':[] }

        error = "Not able to change dataset's access level. User does not own the dataset."
        result['error'] = error

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


def change_dataset_access(cursor, current_user_id, dataset_id, access):

    if access == 'public':
        access_level = 1 #Public
    elif access == 'private':
        access_level = 0 #Private

    qry = """
        UPDATE dataset d
        JOIN guser g ON d.owner_id = g.id
        SET d.is_public = %s
        WHERE d.id = %s AND g.id = %s
    """

    cursor.execute(qry, (access_level, dataset_id, current_user_id,))


def get_dataset(cursor, current_user_id, dataset_id):
    qry = """
       SELECT d.id, d.title, o.label, d.pubmed_id, d.is_public, d.ldesc, d.dtype, d.owner_id
         FROM dataset d
              JOIN organism o ON d.organism_id=o.id
        WHERE d.id = %s
    """
    cursor.execute(qry, (dataset_id,))
    dataset = list()

    for row in cursor:
        if row[4] == 1:
            access_level = 'Public'
        else:
            access_level = 'Private'

        dataset.append({
            'dataset_id': row[0],
            'grid_position': None,
            'mg_grid_position': None,
            'start_col': None,
            'mg_start_col': None,
            'grid_width': 4,
            'mg_grid_width': 6,
            'start_row': None,
            'mg_start_row': None,
            'grid_height': 1,
            'mg_grid_height': 1,
            'title': row[1],
            'organism': row[2],
            'pubmed_id': row[3],
            'access': access_level,
            'ldesc': row[5],
            'dtype': row[6],
            'owner_id': row[7]
        })

    return dataset


if __name__ == '__main__':
    main()
