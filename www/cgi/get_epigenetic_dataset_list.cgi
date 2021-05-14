#!/opt/bin/python3

"""
For a given session_id, returns data on the user's accessible epigenetic datasets.

0. files the user is owner of
1. files that are listed as public

Data structure returned:

{
   datasets: [
    {
      dataset_id: "dataset12.corrected",
      url: "url to file",
      title: "Gabaergic",
      type: "bigwig or bigbed"
    }
   ]
}

"""

import cgi, json
from datetime import datetime
from operator import itemgetter

import os, sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')

    current_user_id = get_user_id_from_session_id(cursor, session_id)
    result = { 'datasets':[] }

    if current_user_id is None:
        raise Exception("ERROR: failed to get user ID from session_id {0}".format(session_id))
    else:
        # A user is logged in
        epiviz_file_qery = "SELECT id, type, url, title FROM dataset_epiviz WHERE owner_id = %s OR is_public = 1"

        cursor.execute(epiviz_file_qery, (current_user_id,))
        for row in cursor:
            result['datasets'].append({'id': row[0], 'type': row[1], 'url': row[2], 'title': row[3]})

        cursor.close()
        cnx.close()

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
