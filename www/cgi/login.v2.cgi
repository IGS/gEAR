#!/opt/bin/python3

"""
For a given username and password, returns a JSON structure:

   { session_id: '??' }

Where the session_id is one of the following:

  -1 - The e-mail wasn't found, or the password was incorrect. Both get the same
       response so the login form can't be used to find registered e-mails; the
       actual reason is logged to stderr.
   ? - Any other value is the session ID

"""

import cgi
import json
import os
import sys
import uuid

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear import passwords

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    user_email = form.getfirst('user-email')
    user_pass = form.getfirst('user-password')
    if user_pass is None:
        user_pass = ""

    result = {}
    result["session_id"] = 0
    result["is_admin"] = 0

    user_qry = "SELECT id, user_name, pass, is_admin FROM guser WHERE email = %s"
    user_email_found = False
    user_id = None
    user_name = None
    is_admin = None
    password_correct = False
    stored_pass = None

    session_qry = "SELECT session_id FROM user_session WHERE user_id = %s"
    profile_qry = "SELECT label, share_id FROM layout where user_id = %s and is_current = 1"  # Should either be 1 or 0 results
    add_session_qry = "INSERT INTO user_session (user_id, session_id) VALUES (%s, %s)"

    cursor.execute(user_qry, (user_email,))
    for row in cursor:
        user_email_found = True

        # Accepts scrypt hashes, and the unsalted SHA3-256 or MD5 hashes of older accounts
        if passwords.verify_password(user_pass, row[2]):
            password_correct = True
            user_id = row[0]
            user_name = row[1]
            stored_pass = row[2]
            is_admin = row[3]
            break

    if password_correct and passwords.needs_rehash(stored_pass):
        upgrade_stored_password(cnx, cursor, user_id, user_pass, stored_pass)

    if not user_email_found:
        # Take as long as a password check would, so timing doesn't reveal which e-mails are registered
        passwords.verify_unknown_user(user_pass)
        print("login.v2.cgi: login failed for {}: e-mail not found".format(user_email), file=sys.stderr)
        result['session_id'] = -1
        print(json.dumps(result))
    else:
        # Check the password for the user
        if password_correct:
            # If the password is valid, check for an active session
            session_id = None
            cursor.execute(session_qry, (user_id,))
            for row in cursor:
                session_id = row[0]

            # If there isn't an active session, create one
            if not session_id:
                session_id = str(uuid.uuid4())
                cursor.execute(add_session_qry, (user_id, session_id))

            # Check for default profile
            cursor.execute(profile_qry, (user_id,))
            profile = "Hearing (site default)"
            profile_share_id = 0

            for row in cursor:
                profile = row[0]
                profile_share_id = row[1]

            result['session_id'] = session_id
            result['user_name'] = user_name
            result['is_admin'] = is_admin
            result['gear_default_domain']  = profile
            result['default_profile_share_id']  = profile_share_id
            print(json.dumps(result))
        else:
            # If the password was invalid, return -1
            print("login.v2.cgi: login failed for {}: incorrect password".format(user_email), file=sys.stderr)
            result['session_id'] = -1
            print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()


def upgrade_stored_password(cnx, cursor, user_id, user_pass, stored_pass):
    """
    Replace an older hash (unsalted SHA3-256 or MD5, or scrypt with weaker settings) with a current
    one. This can only happen at login, the one time the plaintext password is available.

    The login must not fail because of this, and a hash that doesn't fit in the column must never be
    kept: a truncated hash would lock the user out. That happens if guser.pass is still VARCHAR(64)
    and MySQL isn't in strict mode, so the stored length is checked before committing.
    """
    new_pass = passwords.hash_password(user_pass)

    try:
        cursor.execute("UPDATE guser SET pass = %s WHERE id = %s AND pass = %s", (new_pass, user_id, stored_pass))
        cursor.execute("SELECT CHAR_LENGTH(pass) FROM guser WHERE id = %s", (user_id,))
        row = cursor.fetchone()

        if row is None or row[0] != len(new_pass):
            raise ValueError("guser.pass is too short for the new hash; see docs/developer/database_schema.md")

        cnx.commit()
    except Exception as e:
        cnx.rollback()
        print("login.v2.cgi: could not upgrade the password hash of user {}: {}".format(user_id, e), file=sys.stderr)


if __name__ == '__main__':
    main()
