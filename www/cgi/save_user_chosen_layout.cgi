#!/opt/bin/python3

"""
There are UI elements in gEAR where the user selects a layout/dataset collection,
and this script can be called to save this selection in the database.

The database column updated is guser.layout_share_id

Replaces: set_primary_layout.cgi
"""

import cgi, json
import sys, os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')
    
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    layout_share_id = form.getfirst('layout_share_id')

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    result = {'success': 0, 'message': ''}

    if current_user_id is None:
        result['message'] = 'Invalid session_id'
    else:
        layout = get_layout_from_share_id(cursor, layout_share_id)

        if layout is not None and not user_can_choose_layout(cursor, current_user_id, layout):
            result['message'] = "Not able to choose this collection. It belongs to another user and isn't public."
        else:
            # An unknown share ID clears the choice, so the site default is used
            layout_id = layout['id'] if layout else None
            qry = "UPDATE guser SET layout_id = %s WHERE id = %s"
            cursor.execute(qry, (layout_id, current_user_id) )
            cnx.commit()
            result['success'] = 1

    cursor.close()

    print(json.dumps(result))

def get_layout_from_share_id(cursor, layout_share_id):
    """
    Return {id, user_id, is_public, is_domain} for a layout share ID, or None if not found.
    """
    qry = ("SELECT id, user_id, is_public, is_domain FROM layout WHERE share_id = %s")
    cursor.execute(qry, (layout_share_id, ) )
    layout = None

    for (lid, user_id, is_public, is_domain) in cursor:
        layout = {'id': lid, 'user_id': user_id, 'is_public': is_public, 'is_domain': is_domain}

    return layout

def user_can_choose_layout(cursor, user_id, layout):
    """
    Return True if the user may make this layout their chosen one: it's theirs, public, a site
    (domain) layout, or in one of their groups. These are the layouts the collection selectors list
    (see get_user_layouts.cgi), apart from another user's private layout opened by its share link.
    """
    if layout['user_id'] == user_id or layout['is_public'] == 1 or layout['is_domain'] == 1:
        return True

    # Same rule as LayoutCollection.get_by_users_groups()
    qry = """
        SELECT 1
          FROM layout_group_membership lgm
               JOIN user_group_membership ugm ON ugm.group_id = lgm.group_id
         WHERE lgm.layout_id = %s AND ugm.user_id = %s
    """
    cursor.execute(qry, (layout['id'], user_id))
    in_group = False

    for _ in cursor:
        in_group = True

    return in_group

def get_user_id_from_session_id(cursor, session_id):
    """
    Return the user ID for a session ID, or None if not found.
    """
    qry = ("SELECT user_id FROM user_session WHERE session_id = %s")
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id


if __name__ == '__main__':
    main()
