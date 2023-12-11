#!/opt/bin/python3

"""
For a given session_id, returns the user information.

Data structure returned:

{
   email: 'you@whereever.foo',
   name: 'Latasha Smithe'
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
    session_id = form.getvalue('session_id')

    result = {'email':None, 'user_name':None, 'success':0}

    session_qry = """
       SELECT u.email, u.user_name, u.is_admin, u.id, u.institution, u.colorblind_mode, u.updates_wanted
         FROM guser u
              JOIN user_session us ON u.id=us.user_id
        WHERE us.session_id = %s
    """

    cursor.execute(session_qry, (session_id,))
    for row in cursor:
        result = {'email':row[0]
        , 'user_name':row[1]
        , 'is_admin':row[2]
        , 'id':row[3]
        , 'institution':row[4]
        , 'colorblind_mode': row[5]
        , 'updates_wanted':row[6]
        , 'success':1
        }
        break

    cursor.close()
    cnx.close()

    # remove user_id
    if 'id' in result:
        del result['id']

    # if we got here, there isn't a match
    print(json.dumps(result))


if __name__ == '__main__':
    main()
