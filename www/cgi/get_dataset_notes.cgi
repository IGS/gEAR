#!/opt/bin/python3

"""
Returns the notes for a given dataset.

First retrieves any notes the user has for the dataset. Then retrieves any public
notes.

"""

import cgi, json
import configparser
from datetime import datetime
import sys
import mysql.connector
from mysql.connector import errorcode

def main():
    print('Content-Type: application/json\n\n')
    result = {'success': 0, 'notes': [], 'dataset_title': [] }

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
    #print(session_id, file=sys.stderr)
    dataset_id = form.getvalue('dataset_id')
    #print(dataset_id, file=sys.stderr)
    scope = form.getvalue('scope')
    #print(scope, file=sys.stderr)

    current_user_id = get_user_id_from_session_id(cursor, session_id)


    # Is a user even logged in?
    if current_user_id is None:
        result['error'] = 'Please log in to access or create notes for this dataset.'

    else:
        # A user is logged in
        # Get user's notes for the dataset
        result['notes'].extend(get_user_notes(cursor, dataset_id, current_user_id))

        # Add any public notes
        result['notes'].extend(get_public_notes(cursor, dataset_id, current_user_id))

        # No notes found? Then get the dataset's title
        if len(result['notes']) < 1:
            result['dataset_title'].extend(get_dataset_title(cursor, dataset_id))

        result['success'] = 1

    cursor.close()
    cnx.close()

    print(json.dumps(result))

def get_dataset_title(cursor, dataset_id):
    qry = """
        SELECT title
        FROM dataset
        WHERE id = %s
    """
    cursor.execute(qry, (dataset_id,))
    dataset_title = []
    for row in cursor:
        dataset_title.append({
            'dataset_title': row[0]
          })
    return dataset_title

def get_user_notes(cursor, dataset_id, current_user_id):
    qry = """
       SELECT n.id, n.title, n.ldesc, n.user_id, n.dataset_id, n.is_public,
       n.date_added, n.date_last_changed, g.user_name, d.title
         FROM note n
       JOIN guser g ON g.id = n.user_id
       JOIN dataset d ON d.id = n.dataset_id
       WHERE n.dataset_id = %s
         AND n.user_id = %s
    """
    cursor.execute(qry, (dataset_id, current_user_id,))
    notes = []

    for row in cursor:

        if row[5] == 1:
            access_level = 'Public'
        else:
            access_level = 'Private'

        # print(row[6], file=sys.stderr)
        date_added = row[6].isoformat()
        date_last_changed = row[7].isoformat()

        notes.append({
            'id': row[0],
            'title': row[1],
            'ldesc': row[2],
            'user_id': row[3],
            'dataset_id': row[4],
            'access': access_level,
            'date_added': date_added,
            'date_last_changed': date_last_changed,
            'is_owner': 1,
            'owner': row[8],
            'dataset_title': row[9] # will provide heading for note panel
        })

    return notes

def get_public_notes(cursor, dataset_id, current_user_id):
    qry = """
       SELECT n.id, n.title, n.ldesc, n.user_id, n.dataset_id, n.is_public,
       n.date_added, n.date_last_changed, g.user_name, d.title
         FROM note n
       JOIN guser g ON g.id = n.user_id
       JOIN dataset d ON d.id = n.dataset_id
       WHERE n.dataset_id = %s
         AND n.is_public = 1
    """
    cursor.execute(qry, (dataset_id,))
    notes = []

    for row in cursor:

        if row[3] == current_user_id:
            continue
        else:
            date_added = row[6].isoformat()
            date_last_changed = row[7].isoformat()

            notes.append({
                'id': row[0],
                'title': row[1],
                'ldesc': row[2],
                'user_id': row[3],
                'dataset_id': row[4],
                'access': 'Public',
                'date_added': date_added,
                'date_last_changed': date_last_changed,
                'is_owner': 0,
                'owner': row[8],
                'dataset_title': row[9] # will provide heading for note panel
            })

    return notes


def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
