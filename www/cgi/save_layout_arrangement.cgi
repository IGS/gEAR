#!/opt/bin/python3

"""

"""

import cgi
import json
import os
import sys

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
    layout_arrangement_json = form.getvalue('layout_arrangement')

    user = geardb.get_user_from_session_id(session_id)

    layout = geardb.get_layout_by_share_id(layout_share_id)

    if not layout:
        print(json.dumps({'success': 0, 'error': 'Collection not found'}))
        return

    layout_id = layout.id

    # Does user own the dataset ...
    owns_layout = check_layout_ownership(cursor, user.id, layout_id)

    if owns_layout == True and layout_share_id != 0:

        layout_arrangement = json.loads(layout_arrangement_json)

        # Clear the current layout displays for this layout
        qry = """ DELETE FROM layout_displays
                    WHERE layout_id = %s
                """
        cursor.execute(qry, (layout_id,))

        # Loop through all the layout members and update their grid positions
        for display_type in ["single", "multi"]:
            curr_grid_position = 1
            for layout_display in layout_arrangement[display_type]:
                qry = """ INSERT INTO layout_displays
                            (layout_id, display_id, grid_position, grid_width, grid_height, start_col, start_row)
                            VALUES (%s, %s, %s, %s, %s, %s, %s)
                    """

                cursor.execute(qry, (layout_id, layout_display["display_id"], curr_grid_position,
                                    layout_display["grid_width"], layout_display["grid_height"],
                                    layout_display["start_col"], layout_display["start_row"])
                )
                curr_grid_position += 1

        result = { 'success':1 }
        cnx.commit()
        cursor.close()
        cnx.close()
    else:
        error = "Not able to save layout arrangement. User must be logged in and own the layout."
        result = {'success':0, "error": error}


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
