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
    user_id = form.getvalue('user_id')
    dataset_id = form.getvalue('dataset_id')
    display_id = form.getvalue('display_id')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    try:
        query = """
            INSERT INTO dataset_preference (user_id, dataset_id, display_id)
            VALUES (%s, %s, %s)
            ON DUPLICATE KEY UPDATE display_id = VALUES(display_id);
        """
        cursor.execute(query, (user_id, dataset_id, display_id))
        result = dict(success=True)
    except mysql.connector.Error as err:
        result = dict(success=False)

    cnx.commit()
    cursor.close()

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
