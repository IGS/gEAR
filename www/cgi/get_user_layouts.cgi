#!/opt/bin/python3

"""
Queries a list of data for all layouts the user passed (via session ID) is able to see.

Passing the no_domain=1 option will omit layouts from the domain profiles.  Ignored for
anonymous users, who can ONLY see domain and shared ones ones.

Admin will see all
Anonymous users will see only domain, shared
Logged in users will see their own datasets, domain, shared, group

Returns {'user_layouts': [
                        {'id': 123, 'label': 'Layout 1', ... }
                        {'id': 124, 'label': 'Layout 2', ... }
                       ],
         'domain_layouts': [
                        {'id': 312, 'label': 'Layout 3', ... }
                        {'id': 413, 'label': 'Layout 4', ... }
                       ],
         'group_layouts': [
                        {'id': 3122, 'label': 'Layout 5', ... }
                        {'id': 4113, 'label': 'Layout 6', ... }
                       ],
         'shared_layouts': [
                        {'id': 1212, 'label': 'Layout 7', ... }
                       ],
         # key: id, value = folder name
         'folder_names': {},
         # key: id, value = parent_folder_id
         'folder_parents': {}
}

Folder support

- Highlighted profiles
  - Those where gEAR admins have set layout.is_domain = 1
  - Can have nested folders, admin access only
- Your profiles
  - Those uploaded by the user
  - Can have nested folders, owner access only
- Group profiles
  - Any user can create a group and assign users/layouts to it
  - Can have nested folders, group owner can create
- Profiles shared with you
  - Those with entries in dataset_shares or passed via the URL
  - Folders allowed, but current user only
- Other public profiles
  - All other profiles where the user has set layout.is_public=1
  - Can have folders but only admins can create them

All profiles and folders should be nested within these 5 top-level options.

# show all profiles and their labels within a folder
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from gear.serverconfig import ServerConfig

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    no_domain = form.getvalue('no_domain', 0)
    session_id = form.getvalue('session_id')
    layout_share_id = form.getvalue('layout_share_id')
    include_members = form.getvalue('include_members', 1)
    user = geardb.get_user_from_session_id(session_id)

    if no_domain:
        no_domain = int(no_domain)

    if include_members:
        include_members = int(include_members)

    if include_members == 0:
        include_members = False
    elif include_members == 1:
        include_members = True

    result = { 'user_layouts': [],
               'domain_layouts': [],
               'group_layouts': [],
               'shared_layouts': [],
               'public_layouts': [],
               'folders': [],
               'selected': None }

    # Everyone can see public ones
    result['public_layouts'] = geardb.LayoutCollection(include_datasets=include_members).get_public()

    if not no_domain:
        result['domain_layouts'] = geardb.LayoutCollection(include_datasets=include_members).get_domains()

    if user:
        result['user_layouts'] = geardb.LayoutCollection(include_datasets=include_members).get_by_user(user)
        result['group_layouts'] =  geardb.LayoutCollection(include_datasets=include_members).get_by_users_groups(user)

    if layout_share_id:
        result['shared_layouts'] = geardb.LayoutCollection(include_datasets=include_members).get_by_share_id(layout_share_id)

    ## Selected priority:
    ## - A passed share ID
    ## - User has set a saved profile
    ## - Use the site default
    for ltype in ['user', 'domain', 'group', 'shared', 'public']:
        for l in result[ltype + '_layouts']:
            if l.share_id == layout_share_id:
                result['selected'] = l.share_id
                break

    if not result['selected']:
        for l in result['user_layouts']:
            if l.is_current:
                result['selected'] = l.share_id
                break

    if not result['selected']:
        for l in result['domain_layouts']:
            if l.is_current:
                result['selected'] = l.share_id
                break

    # for each layout, determine if the user is the owner
    for ltype in ['user', 'domain', 'group', 'shared', 'public']:
        for l in result[ltype + '_layouts']:
            l.is_owner = True if user and l.user_id == user.id else False
            # delete user_id and layout_id from the layout object
            del l.user_id
            del l.id

    # Doing this so nested objects don't get stringified: https://stackoverflow.com/a/68935297
    print(json.dumps(result, default=lambda o: o.__dict__))

if __name__ == '__main__':
    main()
