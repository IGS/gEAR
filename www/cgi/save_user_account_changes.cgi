#!/opt/bin/python3

"""
Save changes to user's gEAR account

Currenty only saves changes to password
"""

import cgi, json
import sys
import os
import hashlib
import uuid

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    help_id = form.getvalue('help_id')
    user_id = form.getvalue('user_id')
    new_password = form.getvalue('new_password')
    email = form.getvalue('email')
    institution = form.getvalue('institution')
    updates_wanted = form.getvalue('wantUpdates')   # Either "on" or "off" because it's a checkbox
    scope = form.getvalue('scope') # 'password'
    result = {}

    # Password is being changed
    if scope == 'password':
        # Generate encoded password
        encoded_pass = hashlib.md5(new_password.encode('utf-8')).hexdigest()

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
        if updates_wanted:
            if updates_wanted == "on":
                updates_wanted = 1
            else:
                updates_wanted = 0
            mid_query_settings.append(" updates_wanted = %s")
            field_values.append(updates_wanted)
        if new_password:
            mid_query_settings.append(" pass = %s")
            # Generate encoded password
            encoded_pass = hashlib.md5(new_password.encode('utf-8')).hexdigest()
            field_values.append(encoded_pass)

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


def save_new_password(cursor, help_id, encoded_pass):
    qry = """
        UPDATE guser
        SET pass = %s
        WHERE help_id = %s
    """
    cursor.execute(qry, (encoded_pass, help_id,))

def generate_new_help_id(cursor, help_id):
      new_help_id = str(uuid.uuid4())

      qry = """
          UPDATE guser
          SET help_id = %s
          WHERE help_id = %s
      """
      cursor.execute(qry, (new_help_id, help_id,))


if __name__ == '__main__':
    main()
