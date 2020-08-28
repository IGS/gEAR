#!/opt/bin/python3

"""
Returns a list of tags from the database.
Tags are used as an autocomplete for the comment form

"""

import cgi
import configparser
import json
import math
import mysql.connector
from mysql.connector import errorcode
import re
import sys
from xml.dom import minidom

sys.path.append('../..')
import loaderutils

def main():
    print('Content-Type: application/json\n\n')

    config = configparser.ConfigParser()
    config.read('../../gear.ini')
    user_upload_file_base = '../uploads/files'
    user_upload_dest_base = '../datasets_uploaded'
    result = { 'success':0, 'tags': list() }

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

    tag_qry = "SELECT * FROM tag"
    cursor.execute(tag_qry)
    for row in cursor:
        result['tags'].append({'id': row[0], 'label': row[1]})



    cnx.commit()
    cursor.close()
    cnx.close()

    result['success'] = 1
    print(json.dumps(result))



def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id


if __name__ == '__main__':
    main()
