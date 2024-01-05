#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')

    user = geardb.get_user_from_session_id(session_id=session_id)

    displays = {"user":[], "owner":[]}

    if user:
        displays["user"] = geardb.get_displays_by_user_id(user_id=user.id,
            dataset_id=dataset_id)

    # Delete user_id from the displays
    for display in displays["user"]:
        del display["user_id"]

    # Get the owner of the dataset
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    query = "SELECT owner_id FROM dataset WHERE id = %s"
    cursor.execute(query, (dataset_id,))
    (dataset_owner,) = cursor.fetchone()

    # If the user is not the owner of the dataset, then retrieve those displays
    if not (user and user.id == dataset_owner):
        displays["owner"] = geardb.get_displays_by_user_id(user_id=dataset_owner,
            dataset_id=dataset_id)

        # Delete user_id from the displays
        for display in displays["owner"]:
            del display["user_id"]

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(displays))

if __name__ == '__main__':
    main()
