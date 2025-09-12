#!/opt/bin/python3

"""
For a given session_id, returns the user information.

Data structure returned:

{
   email: 'you@whereever.foo',
   user_name: 'Latasha Smithe'
   ... and rest of guser table columns
}

"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')

    result = {'email':None, 'user_name':None, 'success':0}

    session_qry = """
       SELECT u.email, u.user_name, u.is_admin, u.id, u.institution, u.colorblind_mode,
              u.updates_wanted, u.default_org_id, l.share_id
         FROM guser u
              JOIN user_session us ON u.id=us.user_id
              LEFT JOIN layout l ON u.layout_id=l.id
        WHERE us.session_id = %s
    """

    cursor.execute(session_qry, (session_id,))
    row = cursor.fetchone()
    if not row:
        print(json.dumps(result))
        cursor.close()
        return
    result = {
        'email':row[0],
        'user_name':row[1],
        'is_admin':row[2],
        'institution':row[4],
        'colorblind_mode': row[5],
        'updates_wanted':row[6],
        'default_org_id':row[7],
        'layout_share_id':row[8],
        'success':1
    }

    # if we didn't get a layout_share_id then pull it from the gear.ini file
    if result['layout_share_id'] is None:
        default_layout_share_id = geardb.servercfg['content']['default_layout_share_id']
        result['layout_share_id'] = default_layout_share_id

    cursor.close()

    # if we got here, there isn't a match
    print(json.dumps(result))


if __name__ == '__main__':
    main()
