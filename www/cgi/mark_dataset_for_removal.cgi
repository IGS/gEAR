#!/opt/bin/python3

"""
Removes a dataset via the 'delete' button in the dataset_manager
This requires the user to own the dataset in order to mark it for removal.

This script first checks if the user owns the dataset. If the user owns it, the
script then marks the dataset for removal. Actual removal to be performed by
cron job.

If the user does not own the dataset, an error is returned stating that.

Requires:
1) Session id - which contains user_id
2) dataset id to be removed/deleted

"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')

    user = geardb.get_user_from_session_id(session_id)

    result = { 'success': 0 }

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, user.id, dataset_id)

    if owns_dataset == True:

        # Marks dataset for removal. Actual removal to be performed by cron job
        mark_for_removal(cursor, dataset_id)
        result['success'] = 1
        
        cnx.commit()
        cursor.close()
        cnx.close()

        print(json.dumps(result))

    else:
        error = "Not able to remove dataset. User does not own the dataset."
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

    return user_owns_dataset

def mark_for_removal(cursor, dataset_id):
    qry = """
        UPDATE dataset
        SET marked_for_removal = '1'
        WHERE id = %s
    """
    cursor.execute(qry, (dataset_id,))

if __name__ == '__main__':
    main()
