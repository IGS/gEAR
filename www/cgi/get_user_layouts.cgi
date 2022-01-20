#!/opt/bin/python3

"""
Queries a list of data for all layouts the user passed (via session ID) is able to see.

Passing the no_public=1 option will omit layouts from the public profiles.  Ignored for
anonymous users, who can ONLY see public.

Admin will see their own
Anonymous users will see only public
Logged in users will see their own datasets
 and also public unless no_public is passed

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
    no_public = form.getvalue('no_public')
    session_id = form.getvalue('session_id')
    layout_share_id = form.getvalue('layout_share_id')
    user = geardb.get_user_from_session_id(session_id)
    result = { 'layouts':[], 'selected': None }

    if no_public:
        no_public = int(no_public)

    qry = """
          SELECT l.id, l.label, l.is_current, l.user_id, l.share_id, count(lm.id), l.is_domain
            FROM layout l
                 LEFT JOIN layout_members lm ON lm.layout_id=l.id
                 LEFT JOIN dataset d on lm.dataset_id=d.id
                       AND d.marked_for_removal = 0
    """

    if layout_share_id:
        if user:
            if user.id == 0:
                # admin user
                qry += "WHERE (l.user_id = %s or l.share_id = %s)"
                params = [user.id, layout_share_id]
            else:
                # standard user
                if no_public:
                    qry += "WHERE (l.user_id = %s OR l.share_id = %s)"
                    params = [user.id, layout_share_id]
                else:
                    qry += "WHERE (l.user_id = %s OR l.user_id = 0 OR l.share_id = %s)"
                    params = [user.id, layout_share_id]
        else:
            # anonymous user
            qry += "WHERE (l.user_id = %s OR l.share_id = %s)"
            params = [0, layout_share_id]
    else:
        if user:
            if user.id == 0:
                # admin user
                qry += "WHERE l.user_id = %s"
                params = [user.id]
            else:
                # standard user
                if no_public:
                    qry += "WHERE l.user_id = %s"
                    params = [user.id]
                else:
                    qry += "WHERE (l.user_id = %s OR l.user_id = 0)"
                    params = [user.id]
        else:
            # anonymous user
            qry += "WHERE l.user_id = %s"
            params = [0]

    qry += " GROUP BY l.id, l.label, l.is_current, l.user_id, l.share_id"
    
    cursor.execute(qry, params)
    for row in cursor:
        result['layouts'].append({'id': row[0], 'label': row[1], 'is_current': row[2],
                                  'share_id': row[4], 'is_domain': row[6],
                                  'dataset_count': row[5] })
        if row[2] == 1:
            result['selected'] = row[0]

    cursor.close()
    cnx.close()

    #Alphabetize layouts
    result['layouts'].sort(key=lambda a: a['label'].lower())
    print(json.dumps(result))

if __name__ == '__main__':
    main()
