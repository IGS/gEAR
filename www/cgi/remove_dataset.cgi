#!/opt/bin/python3

"""
Removes a dataset via the 'delete' button in the dataset_manager
This requires the user to own the dataset in order to remove it from the database.

This script first checks if the user owns the dataset, then proceeds with the removal
of the dataset if they own it.

If the user does not own the dataset, an error is returned stating that.

Requires:
1) Session id - which contains user_id
2) dataset id to be removed/deleted

"""

import cgi, json
import sys, os

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

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, current_user_id, dataset_id)

    if owns_dataset == True:
        result = { 'success': 1, 'dataset':[] }

    # Delete dataset from referenced tables
        remove_from_layout_displays(cursor, dataset_id)
        remove_from_dataset_shares(cursor, dataset_id)

        # Mark dataset deleted
        remove_from_dataset(cursor, dataset_id)
        cnx.commit()

        cursor.close()
        cnx.close()

        print(json.dumps(result))

    else:
        error = "Not able to remove dataset. User does not own the dataset."
        result = { 'success': 0, 'error': error }
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

def remove_from_layout_displays(cursor, dataset_id):
    qry = """
        DELETE FROM layout_displays
        JOIN dataset_display ON layout_displays.display_id = dataset_display.id
        WHERE dataset_display.dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))

def remove_from_dataset_shares(cursor, dataset_id):
    qry = """
        DELETE FROM dataset_shares
        WHERE dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))

def remove_from_dataset(cursor, dataset_id):
    qry = """
        UPDATE dataset
           SET marked_for_removal = 1
        WHERE id = %s
    """
    cursor.execute(qry, (dataset_id,))

def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
