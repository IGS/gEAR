#!/opt/bin/python3

"""
Given a specific layout ID returns a list of display members.

(NOT CURRENTLY DOING THESE AS WE ARE ASSUMING THAT A SHARED LAYOUT OVERRIDES THE PRIVACY)
A layout is not retrieved if it is private and the user is not the owner
A member is not returned if the dataset is private and the user is not the owner

"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    share_id = form.getvalue('layout_share_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    result = { 'layout_members': {"single":[], "multi":[]}, "message":None }
    layout = geardb.get_layout_by_share_id(share_id)

    #if not layout.is_public and layout.user_id != user.id:
    #    result['message'] = 'Layout is private and user is not the owner.'

    layout.get_singlegene_members()
    # Check each member to see if dataset is public or user is owner
    for member in layout.members:
        result['layout_members']["single"].append(member)
        #dataset = geardb.get_dataset_by_id(member.dataset_id)
        #if dataset.is_public or dataset.user_id == user.id:
        #    result['layout_members']["single"].append(member)

    layout.get_multigene_members()
    for member in layout.members:
        result['layout_members']["multi"].append(member)
        #dataset = geardb.get_dataset_by_id(member.dataset_id)
        #if dataset.is_public or dataset.user_id == user.id:
        #    result['layout_members']["multi"].append(member)

    result["is_owner"] = False
    if user and layout.user_id == user.id:
        result["is_owner"] = True

    #Alphabetize layouts
    print(json.dumps(result))

if __name__ == '__main__':
    main()
