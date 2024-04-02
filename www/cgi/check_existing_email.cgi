#!/opt/bin/python3

"""
Checks to see if the passed e-mail address exists for a user already.

Data structure returned:

{
   email_exists: 1/0
}

"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    email_address = form.getvalue('email')

    # remove all whitespace from the email
    email_address = ''.join(email_address.split())
    

    result = {'email_exists': 0}

    session_qry = """
       SELECT id
         FROM guser
        WHERE email = %s
    """

    cursor.execute(session_qry, (email_address,))
    for row in cursor:
        result = {'email_exists': 1}
        break

    cursor.close()
    cnx.close()

    print(json.dumps(result))


if __name__ == '__main__':
    main()
