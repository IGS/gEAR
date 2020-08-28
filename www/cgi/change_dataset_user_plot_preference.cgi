#!/opt/bin/python3

"""
Change the user's plot preference for a dataset in a layout

Input:
------
dataset_id, plot_preference, session_id

Output:
-------
result{'success': 1}


TODO: Needs classes and methods for
    1) getting a layout_id from user_id
    2) change plot_preference in layout_members table
"""

import cgi, json
import configparser
import os, sys
import mysql.connector
from mysql.connector import errorcode


lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    result = {'success': 0 }

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
    plot_preference = form.getvalue('plot_preference')

    # current_user_id = get_user_id_from_session_id(cursor, session_id)
    user = geardb.get_user_from_session_id(session_id)

    layout_id = None

    # Is a user even logged in?  If no, make no changes
    if user.id is None:
        result['error'] = 'Must be logged in to change plot preference.'
    else:
        # A user is logged in

        #Get the user's current profile
        layout_id = get_current_layout(cursor, user.id)

        #Change user's plot preference
        change_plot_preference(cursor, plot_preference, dataset_id, layout_id)
        cnx.commit()

        result['success'] = 1


    cursor.close()
    cnx.close()

    print(json.dumps(result))

#update the user's plot preference
def change_plot_preference(cursor, plot_preference, dataset_id, layout_id):
    qry = """
        UPDATE layout_members
        SET plot_preference = %s
        WHERE dataset_id = %s
            AND layout_id = %s
    """
    cursor.execute(qry, (plot_preference, dataset_id, layout_id,))

#get user's primary layout profile
def get_current_layout(cursor, current_user_id):
    qry = """
        SELECT id
        FROM layout
        WHERE is_current = 1
            AND user_id = %s
    """
    cursor.execute(qry, (current_user_id,))
    for (layout_id,) in cursor:
        return layout_id


if __name__ == '__main__':
    main()
