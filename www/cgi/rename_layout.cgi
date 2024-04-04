#!/opt/bin/python3

"""
Renames a layout.

Requires:
1) Session id - which contains user_id
2) Layout share ID
3) Layout name to be added
"""

import cgi
import json
import os
import sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.userhistory import UserHistory

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    layout_share_id = form.getvalue('layout_share_id')
    layout_name = form.getvalue('layout_name')

    user = geardb.get_user_from_session_id(session_id)

    if user is None:
        result = {'error':[]}
        error = "Not able to add dataset collection. User must be logged in."
        result['error'] = error
        print(json.dumps(result))
        return;

    layout = geardb.get_layout_by_share_id(layout_share_id)
    if not layout:
        result = {'error':[]}
        error = "Dataset Collection not found."
        result['error'] = error
        print(json.dumps(result))
        return;


    layout.label = layout_name

    layout.save()
    result = {'layout_id': layout.id,
                'layout_label': layout.label,
                'layout_share_id': layout.share_id
        }

    print(json.dumps(result))

if __name__ == '__main__':
    main()
