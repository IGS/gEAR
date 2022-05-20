#!/opt/bin/python3

"""
Given a specific layout ID returns a list of datasets.

For security purposes, the owner of the layout must match the session ID of the user.

"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    layout_id = form.getvalue('layout_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    result = { 'layout_members':[] }
    layout = geardb.get_layout_by_id(layout_id)

    if layout and layout.user_id == user.id:
        layout.get_members()
        result['layout_members'] = layout.members
        
    #Alphabetize layouts
    print(json.dumps(result))

if __name__ == '__main__':
    main()
