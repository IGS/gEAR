#!/opt/bin/python3

"""
Removes a dataset from a profile via the 'remove from profile' button in the dataset_manager

If the user does not own the dataset, an error is returned stating that.

Requires:
1) Session id - which contains user_id
2) dataset id to be removed/deleted

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
    dataset_id = form.getvalue('dataset_id')
    layout_id = form.getvalue('layout_id')
    
    user = geardb.get_user_from_session_id(session_id)
    layout = geardb.Layout(id=layout_id)
    
    if user == None:
        result = { 'error':[] }
        result['error'] = "Must be logged in to remove dataset from profile."
        print(json.dumps(result))

    # check that layout is not site default
    if layout_id == 0:
        result = { 'error':[] }
        result['error'] = "Not able to remove dataset from profile. The Site Default profile cannot be changed."
        print(json.dumps(result))

    else:
        result = { 'dataset':[] }
        layout.get_members()
        layout.remove_member_by_dataset_id(dataset_id)
        print(json.dumps(result))
    
if __name__ == '__main__':
    main()
