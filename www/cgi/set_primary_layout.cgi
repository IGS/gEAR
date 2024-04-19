#!/opt/bin/python3

"""

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
    layout_share_id = form.getvalue('layout_share_id')

    current_user_id = geardb.get_user_id_from_session_id(session_id)

    if current_user_id is None:
        result = { 'error':[] }
        error = "Not able to set primary layout. User must be logged in."
        result['error'] = error
    else:
        set_layout(cursor, current_user_id, layout_share_id)
        cnx.commit()

        # Print the label for the current layout id
        get_label_qry = "SELECT label FROM layout WHERE share_id = %s"
        cursor.execute(get_label_qry, (layout_share_id,))
        layout_label = "Hearing (site default)"    # Set to "Hearing (site default)" as a failsafe
        for row in cursor:
            layout_label = row[0]

        result = { 'success': 1, 'label': layout_label }
        cursor.close()
        cnx.close()

    print(json.dumps(result))


def set_layout(cursor, user_id, layout_share_id):
    # no matter what, we set all the user's layouts current flag to 0
    qry = """
            UPDATE layout
               SET is_current = 0
             WHERE user_id = %s
          """
    cursor.execute(qry, (user_id,))

    qry = """
            UPDATE layout
               SET is_current = 1
             WHERE share_id = %s
    """
    # only let the admin set an admin profile
    if user_id != 0:
        qry += " AND user_id != 0"

    cursor.execute(qry, (layout_share_id,))

    # ! This silently does not apply if the user wants to set a site default layout
    # The client side should handle this case by setting their cookie the preferred layout

if __name__ == '__main__':
    main()
