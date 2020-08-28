#!/opt/bin/python3

"""
Given a dataset_id, returns a list of users that dataset has been shared with.
"""

import cgi, json
import configparser
import sys
import mysql.connector
from mysql.connector import errorcode


def main():
    print('Content-Type: application/json\n\n')

    config = configparser.ConfigParser()
    config.read('../../gear.ini')

    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    cursor = cnx.cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')
    user_id = int(form.getvalue('user_id'))
    to_share = int(form.getvalue('to_share'))

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, current_user_id, dataset_id)

    if owns_dataset == True:
        result = { 'user_list':[] }

        #update dataset_shares
        update_dataset_shares(cursor, dataset_id, user_id, to_share)

        #remove from user profiles if dataset has been unshared
        #if to_share == 0:
        #    remove_dataset_from_layout_members(cursor, dataset_id)

    	# get the shared user_list
        result['user_list'].extend(get_shared_dataset_list(cursor, dataset_id, user_id))
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


def update_dataset_shares(cursor, dataset_id, user_id, to_share):
    qry = """
        UPDATE dataset_shares s
        SET s.is_allowed = %s
        WHERE s.dataset_id = %s AND s.user_id = %s
    """
    cursor.execute(qry, (to_share, dataset_id, user_id,))

def remove_dataset_from_layout_members(cursor, dataset_id):
    qry = """
        DELETE FROM layout_members
        WHERE dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))

def get_shared_dataset_list(cursor, dataset_id, user_id):
    qry = """
        SELECT s.dataset_id, s.user_id, g.user_name, s.is_allowed
        FROM dataset_shares s
        JOIN guser g ON g.id = s.user_id
        WHERE s.dataset_id = %s AND s.user_id = %s
    """
    cursor.execute(qry, (dataset_id, user_id,))
    shared_with_list = list()

    for row in cursor:
        shared_with_list.append({
            'dataset_id': row[0],
            'user_id': row[1],
            'user_name': row[2],
            'is_allowed': row[3]
        })
    return shared_with_list



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

def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
