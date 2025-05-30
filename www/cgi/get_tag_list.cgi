#!/opt/bin/python3

"""
Returns a list of tags from the database.
Tags are used as an autocomplete for the comment form

"""

import cgi
import json
import math
import re, sys, os
from xml.dom import minidom

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
import loaderutils

def main():
    print('Content-Type: application/json\n\n')
    user_upload_file_base = '../uploads/files'
    user_upload_dest_base = '../datasets_uploaded'
    result = { 'success':0, 'tags': list() }

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

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
