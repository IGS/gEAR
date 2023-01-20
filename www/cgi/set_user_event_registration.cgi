#!/opt/bin/python3

"""
This script can be used to both add and remove a user's registration status
to an event.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    event_id = form.getvalue('event_id')
    registration_status = form.getvalue('registration_status')
    user = geardb.get_user_from_session_id(session_id) if session_id else None

    print('Content-Type: application/json\n\n')

    if not user:
        print(json.dumps({'success': 0, 'msg': 'You must be logged in to register or unregister'}))
        sys.exit(0)

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    if registration_status == 1:
        qry = "INSERT INTO event_registration (event_id, user_id) VALUES (%s, %s)"
        cursor.execute(qry, (event_id, user.id))
        print(json.dumps({'success': 1, 'msg': 'You have successfully registered'}))
    
    elif registration_status == 0:
        qry = "DELETE FROM event_registration WHERE event_id = %s AND user_id = %s"
        cursor.execute(qry, (event_id, user.id))
        print(json.dumps({'success': 1, 'msg': 'You have successfully unregistered'}))
    
    else:
        print(json.dumps({'success': 0, 'msg': 'Unrecognized registration status'}))
        sys.exit(0)

    cursor.close()
    cnx.close()

    
if __name__ == '__main__':
    main()
