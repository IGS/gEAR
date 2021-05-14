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
    user_id = int(form.getvalue('user_id'))

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    display = geardb.get_display_by_id(display_id)

    if user_id == display['user_id']:
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
