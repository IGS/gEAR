#!/opt/bin/python3

"""
For a given username and password, returns a JSON structure:

   { session_id: '??' }

Where the session_id is one of the following:

   0 - The user name wasn't found at all.
  -1 - User was found, but the password was incorrect
   ? - Any other value is the session ID

"""

import cgi, json
import hashlib
import os, sys
import uuid

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    user_email = form.getvalue('user-email')
    user_pass = form.getvalue('user-password')
    encoded_pass = hashlib.md5(user_pass.encode('utf-8')).hexdigest()

    result = {'session_id': 0, 'is_admin': 0}

    user_qry = "SELECT id, user_name, pass, is_admin FROM guser WHERE email = %s"
    user_email_found = False
    user_id = None
    user_name = None
    is_admin = None
    password_correct = False

    session_qry = "SELECT session_id FROM user_session WHERE user_id = %s"
    profile_qry = "SELECT label FROM layout where user_id = %s and is_current = 1"  # Should either be 1 or 0 results
    add_session_qry = "INSERT INTO user_session (user_id, session_id) VALUES (%s, %s)"

    cursor.execute(user_qry, (user_email,))
    for row in cursor:
        user_email_found = True

        if row[2] == encoded_pass:
            password_correct = True
            user_id = row[0]
            user_name = row[1]
            is_admin = row[3]
            break

    if user_email_found == False:
        print(json.dumps(result))
    else:
        # Check the password for the user
        if password_correct == True:
            # If the password is valid, check for an active session
            session_id = None
            cursor.execute(session_qry, (user_id,))
            for row in cursor:
                session_id = row[0]

            # If there isn't an active session, create one
            if session_id == None:
                session_id = str(uuid.uuid4())
                cursor.execute(add_session_qry, (user_id, session_id))

            # Check for default profile
            cursor.execute(profile_qry, (user_id,))
            profile = "Hearing (site default)"
            for row in cursor:
                profile = row[0]

            result['user_id'] = user_id
            result['session_id'] = session_id
            result['name'] = user_name
            result['is_admin'] = is_admin
            if profile:
                result['gear_default_domain']  = profile
            print(json.dumps(result))
        else:
            # If the password was invalid, return -1
            result['session_id'] = -1
            print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()
