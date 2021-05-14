#!/opt/bin/python3

"""
Checks if the share_id of the share link is valid
For shared LAYOUT and DATASET
"""

import cgi, json
import sys
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    help_id = form.getvalue('help_id')
    print(help_id, file=sys.stderr)
    result = {}

    # Get user_name if help_id is valid
    user_name = validate_help_id(cursor, help_id)

    if user_name:
        print(user_name, file=sys.stderr)
        result['user_name'] = user_name
        result['success'] = 1
        cnx.commit()

        cursor.close()
        cnx.close()

        print(json.dumps(result))

    else:
        result = { 'error':[], 'success': 0 }

        error = "This link has either expired or the dataset is no longer available."
        result['error'] = error

        print(json.dumps(result))


def validate_help_id(cursor, help_id):
    name = None

    qry = """
       SELECT user_name
       FROM guser
       WHERE help_id = %s
    """
    cursor.execute(qry, (help_id,))

    for row in cursor:
        name = row[0]

    #returns user's name
    return name


if __name__ == '__main__':
    main()
