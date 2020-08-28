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
    layout_id = form.getvalue('layout_id')

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    if current_user_id is None:
        result = { 'error':[] }
        error = "Not able to set primary layout. User must be logged in."
        result['error'] = error
    else:
        set_layout(cursor, current_user_id, layout_id)
        cnx.commit()

        # Print the label for the current layout id
        get_label_qry = "SELECT label FROM layout WHERE id = %s"
        cursor.execute(get_label_qry, (layout_id,))
        layout_label = "Hearing (site default)"    # Set to "Hearing (site default)" as a failsafe
        for row in cursor:
            layout_label = row[0]

        result = { 'success': 1, 'label':layout_label }
        cursor.close()
        cnx.close()

    print(json.dumps(result))


def set_layout(cursor, user_id, layout_id):
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
             WHERE id = %s
    """
    # only let the admin set an admin profile
    if user_id != 0:
        qry += " AND user_id != 0"

    cursor.execute(qry, (layout_id,))


def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
