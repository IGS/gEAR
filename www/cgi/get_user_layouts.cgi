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
}
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    no_domain = form.getvalue('no_domain')
    session_id = form.getvalue('session_id')
    layout_share_id = form.getvalue('layout_share_id')
    user = geardb.get_user_from_session_id(session_id)

    layout_ids_found = set()
    
    if no_domain:
        no_domain = int(no_domain)

    result = { 'user_layouts': [],
               'domain_layouts': [],
               'group_layouts':[],
               'shared_layouts':[],
               'selected': None }

    if not no_domain:
        result['domain_layouts'] = filter_any_previous(layout_ids_found,
                                                       geardb.LayoutCollection().get_domains())

    if user:
        result['user_layouts'] = filter_any_previous(layout_ids_found,
                                                     geardb.LayoutCollection().get_by_user(user))
        result['group_layouts'] = filter_any_previous(layout_ids_found,
                                                      geardb.LayoutCollection().get_by_users_groups(user))

    if layout_share_id:
        result['shared_layouts'] = filter_any_previous(layout_ids_found,
                                                       geardb.LayoutCollection().get_by_share_id(layout_share_id))

    ## Selected priority:
    ## - A passed share ID
    ## - User has set a saved profile
    ## - Use the site default
    for ltype in ['user', 'domain', 'group', 'shared']:
        for l in result[ltype + '_layouts']:
            if l.share_id == layout_share_id:
                result['selected'] = l.id
                break
    
    if not result['selected']:
        for l in result['user_layouts']:
            if l.is_current:
                result['selected'] = l.id
                break

    if not result['selected']:
        for l in result['domain_layouts']:
            if l.is_current:
                result['selected'] = l.id
                break

    # Doing this so nested objects don't get stringified: https://stackoverflow.com/a/68935297
    print(json.dumps(result, default=lambda o: o.__dict__))


def filter_any_previous(ids, new_layouts):
    layouts = []

    for layout in new_layouts:
        if layout.id not in ids:
            layouts.append(layout)
            ids.add(layout.id)

    layouts.sort(key=lambda l: l.label.upper())
            
    return layouts

if __name__ == '__main__':
    main()
