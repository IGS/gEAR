#!/opt/bin/python3

"""
Updates the password for the user, based on the help_id from their password-reset link.

The help_id is replaced with a new one afterwards, so each reset link only works once.
"""

import cgi
import hashlib
import json
import os
import sys
import uuid

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    result = {'error': '', 'success': 0 }

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    password = form.getfirst('password')
    help_id = form.getfirst('help_id')

    if not password:
        result['error'] = "No password provided."
        print(json.dumps(result))
        return

    if not help_id:
        result['error'] = "No help ID provided."
        print(json.dumps(result))
        return

    qry = "SELECT COUNT(*) FROM guser WHERE help_id = %s"
    cursor.execute(qry, (help_id,))
    count = 0
    for row in cursor:
        count = row[0]

    if count == 0:
        result['error'] = "Help ID not found."
        print(json.dumps(result))
        return

    if count != 1:
        result['error'] = ('Duplicate help ID prefixes')
        print(json.dumps(result))
        return

    encoded_pass = hashlib.sha3_256(password.encode('utf-8')).hexdigest()

    # Replace the help_id too, so the reset link (which contains it) cannot be used again
    new_help_id = str(uuid.uuid4())
    qry = "UPDATE guser SET pass = %s, help_id = %s WHERE help_id = %s"
    cursor.execute(qry, (encoded_pass, new_help_id, help_id))

    result['success'] = 1
    print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()



if __name__ == '__main__':
    main()
