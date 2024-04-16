#!/opt/bin/python3

"""
Requires:
1) Session id - which contains user_id
2) Layout ID to which the dataset id added
3) The dataset ID
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

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_id = form.getvalue('layout_share_id')
    dataset_id = form.getvalue('dataset_id')
    make_public_str = form.getvalue('make_dataset_public', "false")
    result = { 'success': 0, 'error': '' }

    make_public = False
    if make_public_str == "true":
        make_public = True

    user = geardb.get_user_from_session_id(session_id)

    layout = geardb.get_layout_by_share_id(share_id)
    layout.load()

    if user is None:
        error = "Not able to add to the layout. User must be logged in."
        result['error'] = error
    else:
        # If make_public is true, set the dataset to public
        if make_public:
            dataset = geardb.get_dataset_by_id(dataset_id)
            dataset.save_change(attribute='is_public', value=1)

        # Determine if we are in "legacy" mode where every member start_col is 1
        legacy = False
        if len(layout.members) > 0:
            if all([m.start_col == 1 for m in layout.members]):
                print("Legacy mode found... rebuilding layout member grid positions...", file=sys.stderr)
                legacy = True

        # If "legacy" mode, adjust the start_col and start_row, as well as mg_start_col and mg_start_row
        if legacy:
            current_col = 1
            current_row = 1
            current_mg_col = 1
            current_mg_row = 1
            for m in layout.members:
                width = m.grid_width
                if current_col + width > 13:
                    current_col = 1
                    current_row += 1
                m.start_col = current_col
                m.start_row = current_row
                current_col += width

                mg_width = m.mg_grid_width
                if current_mg_col + mg_width > 13:
                    current_mg_col = 1
                    current_mg_row += 1
                m.mg_start_col = current_mg_col
                m.mg_start_row = current_mg_row
                current_mg_col += mg_width
                # update the member
                m.save(layout)

        # make sure the user owns the layout
        gpos = len(layout.members) + 1

        # determine the next start_row and start_col. If the grid is full, start a new row
        # get the last start_row and start_col in that row
        row_to_insert = 1
        mg_row_to_insert = 1
        col_to_insert = 1
        mg_col_to_insert = 1
        grid_width = 4
        mg_grid_width = 12
        grid_height = 1
        mg_grid_height = 1

        if len(layout.members) > 0:
            row_to_insert = max([m.start_row for m in layout.members])
            mg_row_to_insert = max([m.mg_start_row for m in layout.members])

            # get max start_col and max_mg_start_col on the last row
            col_to_insert = max([m.start_col + m.grid_width for m in layout.members if m.start_row == row_to_insert])
            mg_col_to_insert = max([m.mg_start_col + m.mg_grid_width for m in layout.members if m.mg_start_row == mg_row_to_insert])

        # If adding this dataset will make the row exceed the grid width, start a new row
        if len(layout.members) > 0:
            if col_to_insert > 12:
                # get grid height of the last row, first start column
                # Want to ensure that the new row starts below the span of the previous row
                grid_height = [m.grid_height for m in layout.members if m.start_row == row_to_insert and m.start_col == 1][0]
                row_to_insert += grid_height
                col_to_insert = 1
            if mg_col_to_insert > 12:
                mg_grid_height = [m.mg_grid_height for m in layout.members if m.mg_start_row == mg_row_to_insert and m.mg_start_col == 1][0]
                mg_row_to_insert += mg_grid_height
                mg_col_to_insert = 1

        if user.id == layout.user_id:
            lm = geardb.LayoutMember(dataset_id=dataset_id, grid_position=gpos, mg_grid_position=gpos,
                                    start_col=col_to_insert, mg_start_col=mg_col_to_insert, grid_width=grid_width, mg_grid_width=mg_grid_width,
                                    start_row=row_to_insert, mg_start_row=mg_row_to_insert, grid_height=grid_height, mg_grid_height=mg_grid_height)


            layout.add_member(lm)
            result['success'] = 1
        else:
            error = "Not able to add to the profile. User doesn't own it"
            result['error'] = error

    print(json.dumps(result))

if __name__ == '__main__':
    main()
