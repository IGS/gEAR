#!/opt/bin/python3

"""
Given a dataset_id, returns a list of users that dataset has been shared with.
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
    user = geardb.get_user_from_session_id(session_id)

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, user.id, dataset_id)

    if owns_dataset == True:
        result = { 'user_list':[] }

    	# get the shared user_list
        result['user_list'].extend(get_shared_dataset_list(cursor, dataset_id))
        result['success'] = 1
        cnx.commit()

        cursor.close()
        cnx.close()

        print(json.dumps(result))

    else:
        result = { 'error':[], 'success': 0 }

        error = "Not able to continue. User does not own the dataset."
        result['error'] = error

        print(json.dumps(result))


def check_dataset_ownership(cursor, current_user_id, dataset_id):
    qry = """
       SELECT d.id, d.owner_id
       FROM dataset d
       WHERE d.id = %s
    """

    cursor.execute(qry, (dataset_id,))
    (_, owner_id) = cursor.fetchone()

    if owner_id == current_user_id:
        user_owns_dataset = True
    else:
        user_owns_dataset = False

    return user_owns_dataset

def get_shared_dataset_list(cursor, dataset_id):
    qry = """
        SELECT s.dataset_id, s.user_id, g.user_name, s.is_allowed
        FROM dataset_shares s
        JOIN guser g ON g.id = s.user_id
        WHERE s.dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))
    shared_with_list = list()

    for row in cursor:
        shared_with_list.append({
            'dataset_id': row[0],
            'user_id': row[1],
            'user_name': row[2],
            'is_allowed': row[3]
        })
    return shared_with_list

if __name__ == '__main__':
    main()
