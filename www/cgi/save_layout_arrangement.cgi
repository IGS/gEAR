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
    layout_id = form.getvalue('layout_id')
    dataset_ids = form.getvalue('dataset_ids')[:-1]
    dataset_widths = form.getvalue('dataset_widths')[:-1]

    user = geardb.get_user_from_session_id(session_id)

    # Does user own the dataset ...
    owns_layout = check_layout_ownership(cursor, user.id, layout_id)

    if owns_layout == True and layout_id != 0:
        d_ids = dataset_ids.split(',')
        d_widths = dataset_widths.split(',')

        #convert widths from pixel values to grid values
        for i, width in enumerate(d_widths):
            # Need float first in case the math created a non-int
            width = int(float(width))

            if width < 210:
                d_widths[i] = 4
            if width > 380 and width < 400:
                d_widths[i] = 8
            if width > 570:
                d_widths[i] = 12
        print(d_widths, file=sys.stderr)

        grid_position = 1
        for d_id, d_width in zip(d_ids, d_widths):
            qry = """ UPDATE layout_members
                      SET grid_position = %s, mg_grid_position = %s, grid_width = %s, mg_grid_width = 12,
                      grid_height = 1, mg_grid_height = 1,
                      start_col = 1, mg_start_col = 1,
                      start_row = 1, mg_start_row = 1
                      WHERE layout_id = %s
                        AND dataset_id = %s
            """
            cursor.execute(qry, (grid_position, grid_position, d_width, layout_id, d_id))
            grid_position += 1

        result = { 'success':1 }
        cnx.commit()
        cursor.close()
        cnx.close()
    else:
        result = { 'error':[] }
        error = "Not able to save layout arrangement. User must be logged in and own the layout."
        result['error'] = error

    print(json.dumps(result))


def add_layout(cursor, user_id, layout_name):
    qry = """
        INSERT INTO layout (user_id, label, is_current)
        VALUES (%s, %s, 0)
    """
    cursor.execute(qry, (user_id, layout_name))
    return cursor.lastrowid

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
