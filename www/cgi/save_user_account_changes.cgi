#!/opt/bin/python3

"""
Save changes to user's gEAR account

Currenty only saves changes to password
"""

import cgi, json
import sys
import os
import hashlib
import re
import uuid

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# Same pattern as the email check in www/js/create_account.js
#  (anchored, so it can't backtrack from every starting position)
EMAIL_PATTERN = re.compile(r"^[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)+$")

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    help_id = form.getfirst('help_id')
    session_id = form.getfirst('session_id')
    new_password = form.getfirst('new_password')
    email = form.getfirst('email')
    institution = form.getfirst('institution')
    colorblind_mode = form.getfirst('colorblind_mode', False)  # string true/false
    updates_wanted = form.getfirst('want_updates', False)  # string true/false
    scope = form.getfirst('scope') # 'password'
    result = {}

    if not session_id:
        result = { 'error':[], 'success': 0 }

        error = "No session id provided, so cannot find user."
        result['error'] = error

        print(json.dumps(result))
        exit()

    user = geardb.get_user_from_session_id(session_id=session_id)
    if user is None:
        result = { 'error':[], 'success': 0 }

        error = "Invalid session id provided, so cannot find user."
        result['error'] = error

        print(json.dumps(result))
        exit()
    user_id = user.id

    # Password is being changed
    if scope == 'password':
        if new_password is None:
            result = { 'error':[], 'success': 0 }

            error = "No new password provided."
            result['error'] = error

            print(json.dumps(result))
            exit()

        # Generate encoded password
        encoded_pass = hashlib.sha3_256(new_password.encode('utf-8')).hexdigest()

        # Save new password
        save_new_password(cursor, help_id, encoded_pass)
        cnx.commit()

        # Update help_id with new string
        generate_new_help_id(cursor, help_id)
        cnx.commit()

        result['success'] = 1

        cursor.close()
        cnx.close()

        print(json.dumps(result))

    elif scope == 'settings_change':
        if email:
            email = email.strip()
            # The profile page always sends the email, so the user's current one (which might predate
            #  these checks) is accepted as is
            email_changed = email.lower() != (user.email or '').strip().lower()
            email_error = check_new_email(cursor, email, user_id) if email_changed else None
            if email_error:
                result = { 'error': email_error, 'success': 0 }

                cursor.close()
                cnx.close()

                print(json.dumps(result))
                return

        begin_query = "UPDATE guser "
        end_query = " WHERE id = %s"
        mid_query = "SET"
        mid_query_settings = []
        field_values = []

        if email:
            mid_query_settings.append(" email = %s")
            field_values.append(email)
        if institution:
            mid_query_settings.append(" institution = %s")
            field_values.append(institution)

        if colorblind_mode:
            colorblind_mode = 1 if colorblind_mode == 'true' else 0
        mid_query_settings.append(" colorblind_mode = %s")
        field_values.append(colorblind_mode)

        if updates_wanted:
            updates_wanted = 1 if updates_wanted == 'true' else 0
        mid_query_settings.append(" updates_wanted = %s")
        field_values.append(updates_wanted)

        if new_password:
            mid_query_settings.append(" pass = %s")
            # Generate encoded password
            encoded_pass = hashlib.sha3_256(new_password.encode('utf-8')).hexdigest()

            field_values.append(encoded_pass)

            # Replace the help_id too, so a password-reset link sent earlier stops working
            mid_query_settings.append(" help_id = %s")
            field_values.append(str(uuid.uuid4()))

        # Build the "SET" values query
        if len(field_values):

            # Safeguard against situation where fields and changed values are not aligned
            if not len(field_values) == len(mid_query_settings):
                result = { 'error':[], 'success': 0 }

                error = "An error has occurred when updating settings. Contact us if this happens again."
                result['error'] = error

                print(json.dumps(result))

            mid_query = mid_query + ",".join(mid_query_settings)
            final_query = begin_query + mid_query + end_query
            # Unpack values so that
            cursor.execute(final_query, (*field_values, user_id,))
            cnx.commit()

        cursor.close()
        cnx.close()

        result['success'] = 1
        print(json.dumps(result))

    else:
        result = { 'error':[], 'success': 0 }

        error = "An error has occurred. Contact us if this happens again."
        result['error'] = error

        print(json.dumps(result))


def check_new_email(cursor, email, user_id):
    """
    Return an error message if the user can't change their email to this address, else None.

    The address must look like an email address (the same pattern the account creation page uses)
    and must not belong to another account (compared case-insensitively, as create_account.cgi's
    lookup is under MySQL's default collation).
    """
    if not EMAIL_PATTERN.fullmatch(email):
        return "Please enter a valid email address."

    qry = "SELECT id FROM guser WHERE LOWER(email) = LOWER(%s) AND id != %s"
    cursor.execute(qry, (email, user_id,))
    if cursor.fetchone() is not None:
        return "Another account already uses this email address."

    return None

def save_new_password(cursor, help_id, encoded_pass):
    """
    Store an encoded password for the user identified by help_id.
    """
    qry = """
        UPDATE guser
        SET pass = %s
        WHERE help_id = %s
    """
    cursor.execute(qry, (encoded_pass, help_id,))

def generate_new_help_id(cursor, help_id):
      """
      Replace the user's help_id with a new random UUID so the old one cannot be reused.
      """
      new_help_id = str(uuid.uuid4())

      qry = """
          UPDATE guser
          SET help_id = %s
          WHERE help_id = %s
      """
      cursor.execute(qry, (new_help_id, help_id,))


if __name__ == '__main__':
    main()
