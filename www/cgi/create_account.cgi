#!/opt/bin/python3

"""
Creates a user account with the submitted form data.

ReCAPTCHA to be added yet.

The returned session ID has some informational content:

0  : Something in the process failed (see 'error' field)
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
    user_name = form.getvalue('first-last')
    institution = form.getvalue('institution')
    user_email = form.getvalue('email')
    user_pass = form.getvalue('password1')
    colorblind_mode = form.getvalue('colorblind_mode')  # checkbox
    get_updates = form.getvalue('email_updates')
    remember_me = 'yes'  # Leaving here in case we want to add it back to the form
    verification_code_long = form.getvalue('verification_code_long')
    verification_code_short = form.getvalue('verification_code_short')
    result = {'success': 0, 'session_id': 0, 'long_session': 0, 'error': ""}

    if geardb.get_verification_code_short_form(verification_code_long) != verification_code_short:
        result['error'] = "Verification code mismatch. Please refresh and try again."
        print(json.dumps(result))
        return

    if get_updates == 'yes':
        get_updates = 1
    else:
        get_updates = 0

    if remember_me == 'yes':
        result['long_session'] = 1

    colorblind_mode = 1 if colorblind_mode == 'yes' else 0

    add_user_sql = """
       INSERT INTO guser (user_name, email, institution, pass, colorblind_mode, updates_wanted, help_id)
       VALUES (%s, %s, %s, %s, %s, %s, %s)
    """

    add_session_sql = """
       INSERT INTO user_session (user_id, session_id)
       VALUES (%s, %s)
    """

    if user_already_exists(user_email, cursor) == True:
       result['error'] = "User already exists"
       result['session_id'] = -1
    else:
        print("DEBUG: adding user to database", file=sys.stderr)
        try:
            encoded_pass = hashlib.md5(user_pass.encode('utf-8')).hexdigest()
            help_id = str(uuid.uuid4())
            cursor.execute(add_user_sql, (user_name, user_email, institution, encoded_pass, colorblind_mode, get_updates, help_id))
            user_id = cursor.lastrowid
            session_id = str(uuid.uuid4())
            result['session_id'] = session_id
            cursor.execute(add_session_sql, (user_id, session_id))
        except Exception as e:
            result['error'] = "There was an error adding the user to the database."
            print(json.dumps(result))
            return

    # All good if we got this far
    result['success'] = 1
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
