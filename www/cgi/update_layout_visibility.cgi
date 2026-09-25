#!/opt/bin/python3

"""
Updates the visibility of the layout (dataset collection), set between private and public.

Input: session_id (required; must own the layout), layout_share_id, visibility ('true' for public).
Output: JSON {success, error}.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    result = {'error': '', 'success': 0 }

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    share_id = form.getfirst('layout_share_id')
    visibility = form.getfirst('visibility')    # 'true' for public, anything else for private

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['error'] = "You must be logged in to change a dataset collection's visibility."
        print(json.dumps(result))
        return

    layout = geardb.get_layout_by_share_id(share_id)
    if not layout:
        result['error'] = "Dataset collection not found."
        print(json.dumps(result))
        return

    # Only the owner may change a collection's visibility
    if layout.user_id != user.id:
        result['error'] = "You can only change the visibility of dataset collections you own."
        print(json.dumps(result))
        return

    # convert JS string boolean to Python boolean
    if visibility == 'true':
        int_visibility = 1
    else:
        int_visibility = 0

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    # update the visibility of the layout
    query = "UPDATE layout SET is_public = %s WHERE id = %s"

    try:
        cursor.execute(query, (int_visibility, layout.id))
        cnx.commit()
        result['success'] = 1
    except Exception as e:
        result['error'] = 'Error: {}'.format(e)
    finally:
        cursor.close()
        cnx.close()

    print(json.dumps(result))



if __name__ == '__main__':
    main()
