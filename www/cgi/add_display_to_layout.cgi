#!/opt/bin/python3

"""
Requires:
1) Session id - which contains user_id
2) Layout ID to which the display id added
3) The display ID
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def add_display_to_layout(session_id, share_id, display_id, grid_width, grid_height):
    result = { 'success': 0, 'error': '' }

    user = geardb.get_user_from_session_id(session_id)
    layout = geardb.get_layout_by_share_id(share_id)

    if layout is None:
        result['error'] = "Layout not found"
        return result

    layout.load()

    if user is None:
        error = "Not able to add to the layout. User must be logged in."
        result['error'] = error
    else:
        # Determine if we are in "legacy" mode where every member start_col is 1
        legacy = False
        if len(layout.members) > 0:
            if all([m.start_col == 1 for m in layout.members]):
                print("Legacy mode found... rebuilding layout member grid positions...", file=sys.stderr)
                legacy = True

        # If "legacy" mode, adjust the start_col and start_row, as well as mg_start_col and mg_start_row
        if legacy:
            for lm_type in [layout.singlegene_members, layout.multigene_members]:
                current_col = 1
                current_row = 1
                for m in lm_type:
                    width = m.grid_width
                    if current_col + width > 13:
                        current_col = 1
                        current_row += 1
                    m.start_col = current_col
                    m.start_row = current_row
                    current_col += width

                    # update the member
                    m.save(layout)

        # determine if this is a single or multigene display (make a dummy LayoutDisplay)
        lm = geardb.LayoutDisplay(display_id=display_id)
        lm.get_is_multigene()
        is_multigene = lm.is_multigene

        # determine the next start_row and start_col. If the grid is full, start a new row
        # get the last start_row and start_col in that row
        row_to_insert = 1
        col_to_insert = 1
        if not grid_width:
            grid_width = 6 if is_multigene else 4
        if not grid_height:
            grid_height = 1

        members_to_use = layout.multigene_members if is_multigene else layout.singlegene_members

        gpos = len(members_to_use) + 1

        if len(members_to_use) > 0:
            row_to_insert = max([m.start_row for m in members_to_use])

            # get max start_col and max_mg_start_col on the last row
            col_to_insert = max([m.start_col + m.grid_width for m in members_to_use if m.start_row == row_to_insert])

        # If adding this dataset will make the row exceed the grid width, start a new row
        if len(members_to_use) > 0:
            if col_to_insert > 12:
                # get grid height of the last row, first start column
                # Want to ensure that the new row starts below the span of the previous row
                # However, if there are no start_col = 1 members, then we can start at the "row_to_insert" value

                start_col_grid_heights = [m.grid_height for m in members_to_use if m.start_row == row_to_insert and m.start_col == 1]
                if len(start_col_grid_heights) > 0:
                    row_to_insert += start_col_grid_heights[0]

                col_to_insert = 1

        if user.id == layout.user_id:
            lm = geardb.LayoutDisplay(display_id=display_id, grid_position=gpos,
                                    start_col=col_to_insert, grid_width=grid_width,
                                    start_row=row_to_insert, grid_height=grid_height)


            layout.add_member(lm)
            result['success'] = 1
        else:
            error = "Not able to add to the profile. User doesn't own it"
            result['error'] = error

    return result

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_id = form.getvalue('layout_share_id')
    display_id = form.getvalue('display_id')
    grid_width = form.getvalue('grid_width')
    grid_height = form.getvalue('grid_height')

    result = add_display_to_layout(session_id, share_id, display_id)

    print(json.dumps(result))

if __name__ == '__main__':
    main()
