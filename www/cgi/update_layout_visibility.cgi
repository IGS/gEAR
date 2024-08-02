#!/opt/bin/python3

"""
Updates the visibility of the layout (dataset collection), set between private and public.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    result = {'error': '', 'success': 0 }

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    share_id = form.getvalue('layout_share_id')
    visibility = form.getvalue('visibility')    # 1 for public, 0 for private

    # convert JS string boolean to Python boolean
    if visibility == 'true':
        int_visibility = 1
    else:
        int_visibility = 0

    # update the visibility of the layout
    query = "UPDATE layout SET is_public = %s WHERE share_id = %s"

    try:
        cursor.execute(query, (int_visibility, share_id))
    except Exception as e:
        result['error'] = 'Error: {}'.format(e)
        print(json.dumps(result))
        cursor.close()
        cnx.close()

    result['success'] = 1
    print(json.dumps(result))

    cnx.commit()
    cursor.close()
    cnx.close()



if __name__ == '__main__':
    main()
