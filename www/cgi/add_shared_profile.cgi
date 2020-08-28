#!/opt/bin/python3

"""
Copy the original layout and give it to the user.
Runs when use clicks "Add this profile"
"""

import cgi, json
import configparser
import sys
import mysql.connector
from mysql.connector import errorcode

def main():
    print('Content-Type: application/json\n\n')
    result = {}

    config = configparser.ConfigParser()
    config.read('../../gear.ini')

    cnx = None
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
    layout_share_id = form.getvalue('layout_share_id')

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    if current_user_id == None:
        result = { 'error': [], 'success': 0 }
        error = "Log in to continue."
        result['error'] = error

        print(json.dumps(result))
    else:
        # get layout_id from layout_share_id
        layout_id = int(layout_share_id.split('-')[1]) - 85

        # get layout name
        layout_name = get_layout_name(cursor, layout_id) #+ " (Shared)"
        cnx.commit()

        # check if user already has the profile, if yes append 'copy' to name
        #already_has = False
        already_has = check_user_layouts(cursor, layout_name, current_user_id)
        if already_has == True:
            layout_name = layout_name + " - copy"

        # Get layout_members
        layout_members = []
        layout_members.extend(get_layout_members(cursor, layout_id))
        cnx.commit()

        # make a copy for the user. set the user as owner of their copy
        # returns new_layout_id to add datasets to it
        new_layout_id = add_shared_profile(cursor, layout_name, current_user_id)
        cnx.commit()


        for member in layout_members:
            # add the datasets of that layout to the user's copy
            add_datasets_to_shared_profile(cursor, new_layout_id, member['dataset_id'], member['grid_position'], member['grid_width'], member['math_preference'])
            cnx.commit()

            # if dataset is private, add it to the user's "Shared with me" area
            if member['access'] == 0:
                #check if dataset has already been shared with user
                dataset_already_shared = check_dataset_shares(cursor, member['dataset_id'], current_user_id)

                #only add to 'shared with me' if NOT already shared
                if dataset_already_shared == False or dataset_already_shared == None:
                    add_dataset_to_dataset_shares(cursor, member['dataset_id'], current_user_id)
                    cnx.commit()

        result['success'] = 1
        #cnx.commit()

        cursor.close()
        cnx.close()

        print(json.dumps(result))


def add_dataset_to_dataset_shares(cursor, dataset_id, current_user_id):
    query = """
        INSERT INTO dataset_shares(dataset_id, user_id, is_allowed)
        VALUES (%s, %s, %s)
    """
    cursor.execute(query, (dataset_id, current_user_id, 1))

def check_dataset_shares(cursor, dataset_id, current_user_id):
    query = ( "SELECT dataset_id FROM dataset_shares WHERE dataset_id = %s AND user_id = %s " )
    cursor.execute(query, (dataset_id, current_user_id,))

    for row in cursor:
        if row[0] != dataset_id: # False or None
            return False #already returns 'None'. I can't force it to False
        else: # dataset exists
            return True

def check_user_layouts(cursor, layout_name, current_user_id):
    # added "LIMIT 1" script errors if multiple entries are found. sql error: "Unread result found"
    query = ( "SELECT label FROM layout WHERE label = %s AND user_id = %s LIMIT 1 " )
    cursor.execute(query, (layout_name, current_user_id,))

    for row in cursor:
        if row[0] == layout_name:
            return True
        else:
            return False

def get_layout_name(cursor, layout_id):
    query = ( "SELECT id, label FROM layout WHERE id = %s " )
    cursor.execute(query, (layout_id,))
    label = ''

    for row in cursor:
        label = row[1]
    return label

def get_layout_members(cursor, layout_id):
    query = """
        SELECT lm.dataset_id, lm.grid_position, lm.grid_width, lm.math_preference, d.is_public
        FROM layout_members lm
        JOIN dataset d ON lm.dataset_id=d.id
        WHERE layout_id = %s
    """
    cursor.execute(query, (layout_id,))
    layout_datasets = []
    for row in cursor:
        layout_datasets.append({
              'dataset_id': row[0],
              'grid_position': row[1],
              'grid_width': row[2],
              'math_preference': row[3],
              'access': row[4]
          })
    return layout_datasets

def add_shared_profile(cursor, layout_name, current_user_id):
    qry = """
        INSERT INTO layout
        (label, user_id)
        VALUES (%s, %s)
    """
    cursor.execute(qry, (layout_name, current_user_id,))
    return cursor.lastrowid

def add_datasets_to_shared_profile(cursor, new_layout_id, dataset_id, grid_position, grid_width, math_preference):
    qry = """
        INSERT INTO layout_members
        (layout_id, dataset_id, grid_position, grid_width, math_preference)
        VALUES (%s, %s, %s, %s, %s)
    """
    cursor.execute(qry, (new_layout_id, dataset_id, grid_position, grid_width, math_preference,))


def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
