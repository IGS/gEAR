#!/opt/bin/python3

"""
For a given session_id, returns data on the user's datasets.  Order of priority:

0.  User logged in with a layout passed to this script
1.  User logged in with a submitted search string
    a.  For searching their own datasets
    b.  For searching the open datasets of others
2.  User logged in with current, saved layout
3.  Default layout + user's private datasets
4.  Default layout (or anonymous user)

Data structure returned:

{
   datasets: [
    {
      dataset_id: "dataset12.corrected",
      title: "Cell-specific RNASeq in the ear",
      ldesc: "Super amazing illustration of cell-specific coloring of the ear cells.",
      access: "Private"
    }
   ]
}

"""

import cgi, json
import configparser
import sys
import mysql.connector
from mysql.connector import errorcode

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
    # print(form, file=sys.stderr)
    session_id = form.getvalue('session_id')
    # session_id = '3e3f5b93-41db-401f-9738-7a059a8e7200'
    dataset_id = form.getvalue('dataset_id')
    math_preference = form.getvalue('math_preference')

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    layout_id = None

    # Is a user even logged in?  If no, make no changes
    if current_user_id is None:
        result['error'] = 'Must be logged in to change math transformation preference.'
    else:
        # A user is logged in

        #Get the user's current profile
        layout_id = get_current_layout(cursor, current_user_id)

        change_math_preference(cursor, math_preference, dataset_id, layout_id)
        cnx.commit()

        result['success'] = 1


    cursor.close()
    cnx.close()

    print(json.dumps(result))

#update the user's math preference
def change_math_preference(cursor, math_preference, dataset_id, layout_id):
    qry = """
        UPDATE layout_members
        SET math_preference = %s
        WHERE dataset_id = %s
            AND layout_id = %s
    """
    cursor.execute(qry, (math_preference, dataset_id, layout_id,))

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

def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
