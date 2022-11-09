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
    layout_id = int(form.getvalue('layout_id'))
    dataset_id = form.getvalue('dataset_id')
    result = { 'success': 0, 'error': '' }

    user = geardb.get_user_from_session_id(session_id)
    layout = geardb.Layout(id=layout_id)
    layout.load()

    if user is None:
        error = "Not able to add to the layout. User must be logged in."
        result['error'] = error
    else:
        # make sure the user owns the layout
        gpos = len(layout.members) + 1

        if user.id == layout.user_id:
            lm = geardb.LayoutMember(dataset_id=dataset_id, grid_position=gpos, grid_width=4, mg_grid_width=12)
            layout.add_member(lm)
            result['success'] = 1
        else:
            error = "Not able to add to the profile. User doesn't own it"
            result['error'] = error

    print(json.dumps(result))

if __name__ == '__main__':
    main()
