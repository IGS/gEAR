#!/opt/bin/python3

"""
Queries a list of data for all layouts the user passed (via session ID) is able to see.

There is also a parameter (layout_share_id) which adds an extra entry outside of the
user's usual list.  This is usually when a permalink is shared with the user which
contains a layout ID.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    layout_share_id = form.getvalue('layout_share_id')
    user = geardb.get_user_from_session_id(session_id)
    result = { 'layouts':[] }

    saved_layout_query = """
                         SELECT id, label, is_current, user_id, share_id 
                           FROM layout 
                          WHERE user_id = 0
    """

    if user:
        saved_layout_query += " OR user_id = %s"
        params = [user.id]
    else:
        params = []

    if layout_share_id:
        saved_layout_query += " OR share_id = %s"
        params.append(layout_share_id)
    
    cursor.execute(saved_layout_query, params)
    for row in cursor:
        result['layouts'].append({'id': row[0], 'label': row[1], 'is_current': row[2],
                                  'share_id': row[4], 'is_domain': 0 if int(row[3]) else 1 })

    cursor.close()
    cnx.close()

    #Alphabetize layouts
    result['layouts'].sort(key=lambda a: a['label'].lower())

    print(json.dumps(result))

if __name__ == '__main__':
    main()
