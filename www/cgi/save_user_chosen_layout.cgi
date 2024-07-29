#!/opt/bin/python3

"""
There are UI elements in gEAR where the user selects a layout/dataset collection,
and this script can be called to save this selection in the database.

The database column updated is guser.layout_share_id

Replaces: set_primary_layout.cgi
"""

import cgi, json
import sys, os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')
    
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    layout_share_id = form.getvalue('layout_share_id')

    current_user_id = get_user_id_from_session_id(cursor, session_id)
    layout_id = get_layout_id_from_share_id(cursor, layout_share_id)

    result = {'success': 0, 'message': ''}

    if current_user_id is None:
        result['message'] = 'Invalid session_id'
    else:
        qry = "UPDATE guser SET layout_id = %s WHERE id = %s"
        cursor.execute(qry, (layout_id, current_user_id) )
        cnx.commit()
        result['success'] = 1

    cursor.close()

    print(json.dumps(result))

def get_layout_id_from_share_id(cursor, layout_share_id):
    qry = ("SELECT id FROM layout WHERE share_id = %s")
    cursor.execute(qry, (layout_share_id, ) )
    layout_id = None

    for (lid,) in cursor:
        layout_id = lid

    return layout_id    

def get_user_id_from_session_id(cursor, session_id):
    qry = ("SELECT user_id FROM user_session WHERE session_id = %s")
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id


if __name__ == '__main__':
    main()
