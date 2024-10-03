#!/opt/bin/python3

"""

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

    num_entries = form.getvalue('num_entries', 5)

    ## add the user
    user = geardb.get_user_from_session_id(form.getvalue('session_id'))
    entries = list()

    if user is not None:
        history = UserHistory(user_id=user.id)
        entries = history.get_latest_entries(entry_count=num_entries)

    print(json.dumps(entries))

if __name__ == '__main__':
    main()
