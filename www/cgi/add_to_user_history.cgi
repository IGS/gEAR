#!/opt/bin/python3

"""
Requires:

1) 'session_id' -> to look up user info
2) 'entry_category' -> enumerated list from userhistory.py
3) 'label' -> Label to be shown in the user history table

"""

import cgi
import os
import sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.userhistory import UserHistory


def main():
    print('Content-Type: application/json\n\n')
    
    form = cgi.FieldStorage()

    # build an argument dict, passing all forward we don't explicitly handle
    historyargs = dict()
    for key in form:
        historyargs[key] = form[key].value

    ## add the user
    user = geardb.get_user_from_session_id(historyargs['session_id'])
    historyargs['user_id'] = user.id
    
    if user is not None:
        # Log the addition
        history = UserHistory()
        history.add_record(**historyargs)

if __name__ == '__main__':
    main()
