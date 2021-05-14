#!/opt/bin/python3

"""
Removes a layout via the 'delete' button in the dataset_manager

This requires the user to own the layout in order to remove it from the database.

Also checks to make sure that somehow they aren't trying to remove the site
default layout

If the user does not own the layout, an error is returned stating that.

Requires:
1) Session id - which contains user_id
2) Layout id to be removed/deleted

"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    layout_id = int(form.getvalue('layout_id'))

    user = geardb.get_user_from_session_id(session_id)
    layout = geardb.Layout(id=layout_id)

    # Does user own the dataset ...
    owns_layout = check_layout_ownership(cursor, user.id, layout.id)

    result = { 'success': None }

    if owns_layout == True and layout.id != 0:
        layout.remove()
        result['success'] = 1
        cnx.commit()
    else:
        error = "Not able to remove layout. User does not own it."
        result['success'] = 0
        result['error'] = error

    cursor.close()
    cnx.close()

    print(json.dumps(result))


def check_layout_ownership(cursor, current_user_id, layout_id):
    qry = """
       SELECT l.id, l.user_id
       FROM layout l
       WHERE l.id = %s
    """
    cursor.execute(qry, (layout_id,))

    # default: Assume user does not own dataset
    user_owns_layout = False

    for row in cursor:
        if row[1] == current_user_id:
            user_owns_layout = True

    return user_owns_layout


if __name__ == '__main__':
    main()
