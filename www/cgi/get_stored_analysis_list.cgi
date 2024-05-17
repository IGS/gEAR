#!/opt/bin/python3

"""
Used by sc_workbench.html, this script gets a list of the H5AD datasets the user can view.

Returns four keyed sets 'primary', 'public', 'user_saved' (if a user can be pulled from the
session) and 'user_unsaved'
"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')
    user = geardb.get_user_from_session_id(session_id)
    result = {'primary': [], 'public': [], 'user_saved': [], 'user_unsaved': []}

    acollection = geardb.AnalysisCollection()
    acollection.get_all_by_dataset_id(user_id=user.id, session_id=session_id, dataset_id=dataset_id)

    result['primary'] = acollection.primary
    result['public'] = acollection.public
    result['user_saved'] = acollection.user_saved
    result['user_unsaved'] = acollection.user_unsaved

    ## get the vetting for each
    for atype in result:
        for ana in result[atype]:
            ana.discover_vetting(current_user_id=user.id)

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
