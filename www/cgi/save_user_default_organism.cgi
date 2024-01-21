#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    default_org_id = form.getvalue('default_org_id')

    user = geardb.get_user_from_session_id(session_id=session_id)
    user_id = user.id

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    try:
        query = 'UPDATE guser SET default_org_id = %s WHERE id = %s'
        cursor.execute(query, (default_org_id, user_id))
        result = dict(success=True)
    except mysql.connector.Error as err:
        #print("Something went wrong: {}".format(err), file=sys.stderr)
        result = dict(success=False)

    cnx.commit()
    cursor.close()

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
