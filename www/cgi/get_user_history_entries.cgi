#!/opt/bin/python3

"""
get_user_history_entries.cgi - Return the user's most recent activity history entries.

Input: session_id, num_entries (default 5).
Output: JSON list of history entries (empty list if no valid session).
"""

import cgi
import json
import os
import sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.userhistory import UserHistory

DEFAULT_NUM_ENTRIES = 5


def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()

    # num_entries goes into the query's LIMIT clause, so only accept a positive number
    try:
        num_entries = int(form.getfirst('num_entries', DEFAULT_NUM_ENTRIES))
    except ValueError:
        num_entries = DEFAULT_NUM_ENTRIES
    if num_entries < 1:
        num_entries = DEFAULT_NUM_ENTRIES

    ## add the user
    user = geardb.get_user_from_session_id(form.getfirst('session_id'))
    entries = list()

    if user is not None:
        history = UserHistory(user_id=user.id)
        entries = history.get_latest_entries(entry_count=num_entries)

    print(json.dumps(entries))

if __name__ == '__main__':
    main()
