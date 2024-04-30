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

        curr_grid_position = 1

        layout_arrangement = json.loads(layout_arrangement_json)



        # Loop through all the layout members and update their grid positions
        for dataset_id, vals in layout_arrangement["single"].items():
            mgvals = layout_arrangement["multi"][dataset_id]
            vals["grid_position"] = curr_grid_position
            mgvals["grid_position"] = curr_grid_position

            qry = """ UPDATE layout_members
                        SET grid_position = %s, mg_grid_position = %s, grid_width = %s, mg_grid_width = %s,
                        grid_height = %s, mg_grid_height = %s,
                        start_col = %s, mg_start_col = %s,
                        start_row = %s, mg_start_row = %s
                        WHERE layout_id = %s
                            AND dataset_id = %s
                """
            cursor.execute(qry, (vals["grid_position"], mgvals["grid_position"], vals["grid_width"], mgvals["grid_width"],
                                vals["grid_height"], mgvals["grid_height"], vals["start_col"], mgvals["start_col"],
                                vals["start_row"], mgvals["start_row"], layout_id, dataset_id))

            curr_grid_position += 1

        result = { 'success':1 }
        cnx.commit()
        cursor.close()
        cnx.close()
    else:
        result = {'success':0}
        error = "Not able to save layout arrangement. User must be logged in and own the layout."
        result['error'] = error

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
        print(row, file=sys.stderr)
        print("User ID: " + str(current_user_id), file=sys.stderr)
        if row[1] == current_user_id:
            user_owns_layout = True

    return user_owns_layout

if __name__ == '__main__':
    main()
