#!/opt/bin/python3

"""
Gets a user's gene carts. Returns {'id': 123, 'label': 'my_gene_cart'}
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    current_user_id = get_user_id_from_session_id(cursor, session_id)
    result = { 'gene_carts':[] }

    # Does the user have a current, saved layout?
    layout_id = None

    if current_user_id is None:
        raise Exception("ERROR: failed to get user ID from session_id {0}".format(session_id))
    else:
        # A user is logged in
        gene_cart_query = "SELECT id, label FROM gene_cart WHERE user_id = %s"

        cursor.execute(gene_cart_query, (current_user_id,))
        for row in cursor:
            result['gene_carts'].append({'id': row[0], 'label': row[1]})

    cursor.close()
    cnx.close()

    #Alphabetize gene carts
    result['gene_carts'].sort(key=lambda a: a['label'].lower())

    print(json.dumps(result))

def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
