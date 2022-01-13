#!/opt/bin/python3

"""
Gets a user's gene carts and any others the user has access
to, with an additional shared one if a share_id is passed

Returns {'user_carts': [
                        {'id': 123, 'label': 'my_gene_cart'}
                        {'id': 124, 'label': 'my_gene_cart'}
                       ],
         'public_carts': [
                        {'id': 12, 'label': 'my_gene_cart'}
                        {'id': 13, 'label': 'my_gene_cart'}
                       ],
         'domain_carts': [
                        {'id': 312, 'label': 'my_gene_cart'}
                        {'id': 413, 'label': 'my_gene_cart'}
                       ],
         'group_carts': [
                        {'id': 3122, 'label': 'my_gene_cart'}
                        {'id': 4113, 'label': 'my_gene_cart'}
                       ],
         'shared_carts': [
                        {'id': 1212, 'label': 'my_gene_cart'}
                       ],
}

Carts will NOT be duplicated across multiple categories, which would cause UI
issues.  Instead each cart is only included in the first category in which
it fits below, starting with the top.

----------------------------------------------------------------------------

- Domain carts are named gene carts created by the site admins, thought to be
of general interest.

- User carts are those owned by the user.

- Group carts are any shared with a group the user belongs to.

- Shared carts are those explicitly shared with the current user by URL, even
if private.  This list will usually only contain one cart since the system
doesn't currently have a mechanism to store multiple cart shares.

- Public carts are those owned by other users who have toggled them as public.

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
    share_id = form.getvalue('share_id')
    current_user = geardb.get_user_from_session_id(session_id)
    current_user_id = current_user.id
    result = { 'domain_carts':[], 'gene_carts':[], 'public_carts':[],
               'shared_carts':[], 'user_carts':[] }
 
    # Does the user have a current, saved layout?
    layout_id = None

    # Track the cart IDs already stored so we don't duplicate
    carts_found = set()

    if current_user_id is None:
        raise Exception("ERROR: failed to get user ID from session_id {0}".format(session_id))
    else:
        # A user is logged in
        gene_cart_query = "SELECT id, label, share_id FROM gene_cart WHERE user_id = %s"
        query_args = [current_user_id,]

        if share_id:
            gene_cart_query += " OR share_id = %s"
            query_args.append(share_id)

        cursor.execute(gene_cart_query, query_args)
        for row in cursor:
            if row[2] == share_id:
                result['shared_carts'].append({'id': row[0], 'label': row[1], 'share_id': row[2]})
            else:
                result['gene_carts'].append({'id': row[0], 'label': row[1], 'share_id': row[2]})

    cursor.close()
    cnx.close()

    #Alphabetize gene carts
    result['gene_carts'].sort(key=lambda a: a['label'].lower())
    print(json.dumps(result))

if __name__ == '__main__':
    main()
