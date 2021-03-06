#!/opt/bin/python3

"""
Creates a user account with the submitted form data.

ReCAPTCHA to be added yet.

The returned session ID has some informational content:

-1 : User already exists
?  : Any other (uuid4) value is the session ID

"""

import cgi
import json
import hashlib
import os
import sys
import uuid

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    user_name = form.getvalue('inputName')
    institution = form.getvalue('inputInstitution')
    user_email = form.getvalue('inputEmail')
    user_pass = form.getvalue('inputPassword')
    get_updates = form.getvalue('getUpdates')
    remember_me = form.getvalue('rememberMe')
    result = {'session_id': 0, 'long_session': 0}

    if get_updates == 'yes':
        get_updates = 1
    else:
        get_updates = 0

    if remember_me == 'yes':
        result['long_session'] = 1

    add_user_sql = """
       INSERT INTO guser (user_name, email, institution, pass, updates_wanted, help_id)
       VALUES (%s, %s, %s, %s, %s, %s)
    """

    add_session_sql = """
       INSERT INTO user_session (user_id, session_id)
       VALUES (%s, %s)
    """

    if user_already_exists(user_email, cursor) == True:
       result['session_id'] = -1
    else:
        encoded_pass = hashlib.md5(user_pass.encode('utf-8')).hexdigest()
        help_id = str(uuid.uuid4())
        cursor.execute(add_user_sql, (user_name, user_email, institution, encoded_pass, get_updates, help_id))
        user_id = cursor.lastrowid
        session_id = str(uuid.uuid4())
        result['session_id'] = session_id
        cursor.execute(add_session_sql, (user_id, session_id))

    print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()


def user_already_exists(user_email, curs):
    is_found = False

    qry = "SELECT id FROM guser WHERE email = %s"
    curs.execute(qry, (user_email,))
    for row in curs:
        is_found = True

    return is_found


if __name__ == '__main__':
    main()
