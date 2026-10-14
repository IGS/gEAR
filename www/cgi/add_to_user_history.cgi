#!/opt/bin/python3

"""
Requires:

1) 'session_id' -> to look up user info
2) 'entry_category' -> enumerated list from userhistory.py
3) 'label' -> Label to be shown in the user history table

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

    # build an argument dict, passing all forward we don't explicitly handle
    historyargs = dict()
    for key in form:
        historyargs[key] = form[key].value

    ## add the user
    user = geardb.get_user_from_session_id(historyargs.get('session_id'))
    if user is None:
        print(json.dumps({'success': 0, 'error': 'User must be logged in'}))
        return
    historyargs['user_id'] = user.id

    # Log the addition (this used to print no response body at all)
    try:
        history = UserHistory()
        history.add_record(**historyargs)
    except Exception as e:
        print(json.dumps({'success': 0, 'error': str(e)}))
        return

    print(json.dumps({'success': 1}))

if __name__ == '__main__':
    main()
