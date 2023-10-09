#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    display_id = form.getvalue('id')
    session_id = int(form.getvalue('session_id'))

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    display = geardb.get_display_by_id(display_id)

    user = geardb.get_user_from_session_id(session_id=session_id)

    # If no display found for this ID, treat as if it were already deleted.
    if not display:
        result = dict(success=True)
    elif user.id == display['user_id']:
        query = "DELETE FROM dataset_display where id = %s"
        cursor.execute(query, (display_id,))
        cnx.commit()
        result = dict(success=True)
    else:
        result = dict(success=False)

    cursor.close()
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
