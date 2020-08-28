#!/opt/bin/python3

"""
Requires:
1) Session id - which contains user_id
2) Layout ID to which the dataset id added
3) The dataset ID
4) is_shared boolean - '1' means the dataset was shared with the user
                     - '0' means the dataset was not shared. (they already own it or its a public dataset)

Also adds the dataset to shared_dataset table. Users can then search/list these
shared datasets
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    layout_id = int(form.getvalue('layout_id'))
    is_shared = int(form.getvalue('is_shared'))

    # Do now allow additions to layout_id=0 for now.
    if layout_id == 0:
       raise Exception("ERROR: Not allowed to add to the Site default profile")

    user = geardb.get_user_from_session_id(session_id)
    layout = geardb.Layout(id=layout_id)

    if user is None:
        result = { 'error':[] }
        error = "Not able to add to the layout. User must be logged in."
        result['error'] = error
    else:
        # dataset was shared, so add it to user's 'Shared With Me'
        if is_shared == 1:
            #dataset_id is ACTUALLY share_id
            share_id = form.getvalue('dataset_id')

            #print("share_id: ", share_id, file=sys.stderr)
            dataset_id = get_dataset_id_from_share_id(cursor, share_id)
            
            #print("got from dataset_shares: ", dataset_id, file=sys.stderr)
            add_to_dataset_shares(cursor, dataset_id, user.id)

        else:
            dataset_id = form.getvalue('dataset_id')

        lm = geardb.LayoutMember(dataset_id=dataset_id, grid_position=1, grid_width=4)
        layout.add_member(lm)
        result = { 'success': 1 }
        cnx.commit()
        cursor.close()
        cnx.close()

    print(json.dumps(result))

def get_dataset_id_from_share_id(cursor, share_id):
    qry = """
        SELECT d.id
        FROM dataset d
        WHERE d.share_id = %s
    """
    cursor.execute(qry, (share_id,))
    for row in cursor:
        dataset_id = row[0]
    return dataset_id

def add_to_dataset_shares(cursor, dataset_id, user_id):
    qry = """
        INSERT INTO dataset_shares(dataset_id, user_id, is_allowed)
        VALUES (%s, %s, %s)
    """
    cursor.execute(qry, (dataset_id, user_id, 1))


if __name__ == '__main__':
    main()
