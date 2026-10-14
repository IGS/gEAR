#!/opt/bin/python3

"""
Requires:
1) Session id - which contains user_id
2) Layout ID to which the display id added
3) The display ID

The user must be able to edit the layout (see user_can_edit_layout) and be allowed to use the
display (see display_permission_error).
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# Number of columns in the layout grid
GRID_COLUMNS = 12

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    share_id = form.getfirst('layout_share_id')
    display_id = form.getfirst('display_id')
    result = { 'success': 0, 'error': '' }

    user = geardb.get_user_from_session_id(session_id)
    layout = geardb.get_layout_by_share_id(share_id)

    if layout is None:
        result['error'] = "Layout not found"
        print(json.dumps(result))
        return

    layout.load()

    # All permission checks come before the legacy rebuild below, which saves the layout's members
    if user is None:
        result['error'] = "Not able to add to the layout. User must be logged in."
    elif not user_can_edit_layout(user, layout):
        result['error'] = "Not able to add to the profile. User doesn't own it"
    else:
        result['error'] = display_permission_error(user, display_id)

    if result['error']:
        print(json.dumps(result))
        return

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
                if current_col + width > GRID_COLUMNS + 1:
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

    # The new display keeps the default size (not the size of its neighbors)
    grid_width = 6 if is_multigene else 4
    grid_height = 1

    members_to_use = layout.multigene_members if is_multigene else layout.singlegene_members

    gpos = len(members_to_use) + 1
    row_to_insert, col_to_insert = next_position(members_to_use, grid_width, grid_height)

    lm = geardb.LayoutDisplay(display_id=display_id, grid_position=gpos,
                              start_col=col_to_insert, grid_width=grid_width,
                              start_row=row_to_insert, grid_height=grid_height)

    layout.add_member(lm)
    result['success'] = 1

    print(json.dumps(result))

def user_can_edit_layout(user, layout):
    """
    Return True if the user may change the layout's members.

    Only the layout's owner can for now. (Group edit access to layouts is planned; add it here.)
    """
    return user.id == layout.user_id

def display_permission_error(user, display_id):
    """
    Return an error message if the user may not add this display to their layout, else "".

    The display must be the user's own or one saved by the dataset's owner (the owner's displays are
    shared with everyone who can see the dataset), and the user must be able to access the dataset.
    """
    display = geardb.get_display_by_id(display_id) if display_id else None
    if not display:
        return "Display not found"

    dataset = geardb.get_dataset_by_id(display.get('dataset_id'))
    if dataset is None:
        return "Display not found"

    if display.get('user_id') not in (user.id, dataset.owner_id):
        return "Not able to add this display. It belongs to another user."

    if not user_can_access_dataset(user, dataset):
        return "Not able to add this display. You don't have access to its dataset."

    return ""

def user_can_access_dataset(user, dataset):
    """
    Return True if the dataset is public, owned by the user, or shared with them.
    """
    if dataset.is_public == 1 or dataset.owner_id == user.id:
        return True

    # Same rule as DatasetCollection.get_shared_with_user()
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    qry = "SELECT 1 FROM dataset_shares WHERE dataset_id = %s AND user_id = %s AND is_allowed = 1"
    cursor.execute(qry, (dataset.id, user.id))
    is_shared = cursor.fetchone() is not None
    cursor.close()
    cnx.close()
    return is_shared

def next_position(members, grid_width, grid_height):
    """
    Return (start_row, start_col) for a new tile: after the last tile in the last row if it fits
    there, else at the start of a new row below every existing tile.
    """
    if not members:
        return 1, 1

    row = max(m.start_row for m in members)
    col = max(m.start_col + m.grid_width for m in members if m.start_row == row)

    fits_in_row = col + grid_width - 1 <= GRID_COLUMNS
    if fits_in_row and not any(overlaps(m, row, col, grid_width, grid_height) for m in members):
        return row, col

    # A tile from any earlier row can be tall enough to reach this far down, so go below them all
    return max(m.start_row + m.grid_height for m in members), 1

def overlaps(member, row, col, grid_width, grid_height):
    """
    Return True if the member's tile overlaps a tile of this size at (row, col).
    """
    return (member.start_row < row + grid_height and row < member.start_row + member.grid_height
            and member.start_col < col + grid_width and col < member.start_col + member.grid_width)

if __name__ == '__main__':
    main()
