#!/opt/bin/python3

"""
Checks if the share_id of the share link is valid
For shared LAYOUT and DATASET
"""

import cgi, json
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
    share_id = form.getvalue('share_id')
    scope = form.getvalue('scope') #'permalink', 'profile', or 'dataset'
    result = {}

    user = geardb.get_user_from_session_id(session_id)

    if scope == 'permalink':
        valid_permalink = validate_share_id(cursor, share_id)

        if valid_permalink == True:
            result['success'] = 1
            cnx.commit()

            cursor.close()
            cnx.close()

            print(json.dumps(result))

        else:
            result = { 'error':[], 'success': 0 }

            error = "This permalink has either expired or the dataset is no longer available."
            result['error'] = error

            print(json.dumps(result))

    elif scope == 'profile':
        layout_id = int(share_id.split('-')[1]) - 85
        valid_share = validate_layout_id(cursor, layout_id)

        if valid_share == True:
            result['success'] = 1
            cnx.commit()

            cursor.close()
            cnx.close()

            print(json.dumps(result))

        else:
            result = { 'error':[], 'success': 0 }

            error = "Share link has either expired or the profile is no longer available."
            result['error'] = error

            print(json.dumps(result))

    else: #scope == 'dataset'
        # Check if share_id exists
        valid_share = validate_share_id(cursor, share_id)

        # share_id is valid
        if valid_share == True:

            #Check if user already has the share_id
            already_has = check_dataset_shares(cursor, share_id, user.id)
            cnx.commit()

            #user does NOT already have this dataset shared with them
            if already_has == False or already_has == None:
                result['success'] = 1

                cursor.close()
                cnx.close()

                print(json.dumps(result))

            else: #already_has == True
                result = { 'error':[], 'success': 0 }

                error = "It looks like you already have this dataset. If the owner has not unshared it, you can find in it in your 'Shared with me' datasets."
                result['error'] = error

                print(json.dumps(result))

        else: #valid_share == False
            result = { 'error':[], 'success': 0 }

            error = "Share link has either expired or the dataset is no longer available."
            result['error'] = error

            print(json.dumps(result))


def validate_layout_id(cursor, layout_id):
    qry = ( "SELECT id FROM layout WHERE id = %s" )
    cursor.execute(qry, (layout_id,))
    for row in cursor:
        if row[0] == layout_id:
            return True
        else:
            return False

def validate_share_id(cursor, share_id):
    qry = ( "SELECT share_id FROM dataset WHERE share_id = %s" )
    cursor.execute(qry, (share_id,))

    for row in cursor:
        if row[0] == share_id:
            return True
        else:
            return False

def check_dataset_shares(cursor, share_id, current_user_id):
    qry = """
        SELECT d.share_id
        FROM dataset d
        JOIN dataset_shares s ON d.id = s.dataset_id
        WHERE d.share_id = %s AND s.user_id = %s
    """
    cursor.execute(qry, (share_id, current_user_id,))

    for row in cursor:
        if row[0] == share_id:
            return True
        else:
            return False

def check_dataset_ownership(cursor, current_user_id, dataset_id):
    qry = """
       SELECT d.id, d.owner_id
       FROM dataset d
       WHERE d.id = %s
    """
    cursor.execute(qry, (dataset_id,))

    # default: Assume user does not own dataset
    user_owns_dataset = False

    for row in cursor:
        # Change access if user owns the dataset
        if row[1] == current_user_id:
            user_owns_dataset = True

    return user_owns_dataset

if __name__ == '__main__':
    main()
