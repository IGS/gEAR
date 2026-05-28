#!/opt/bin/python3

"""
Updates the password for the user, based on their session ID.
"""

import cgi
import hashlib
import json
import os
import sys

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

    if password is None:
        result['error'] = "No password provided."
        print(json.dumps(result))
        return

    qry = "SELECT COUNT(*) FROM guser WHERE help_id = %s"
    cursor.execute(qry, (help_id,))
    count = 0
    for row in cursor:
        count = row[0]

    if count != 1:
        result['error'] = ('Duplicate help ID prefixes')
        print(json.dumps(result))
        return

    encoded_pass = hashlib.sha256(password.encode('utf-8')).hexdigest()

    qry = "UPDATE guser SET pass = %s WHERE help_id = %s"
    cursor.execute(qry, (encoded_pass, help_id))

    result['success'] = 1
    print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()



if __name__ == '__main__':
    main()
