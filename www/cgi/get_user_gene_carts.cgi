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
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_id = form.getvalue('share_id')
    current_user = geardb.get_user_from_session_id(session_id)
    current_user_id = current_user.id
    result = { 'domain_carts':[], 'group_carts':[], 'public_carts':[],
               'shared_carts':[], 'user_carts':[] }
 
    # Track the cart IDs already stored so we don't duplicate
    cart_ids_found = set()

    if current_user_id is None:
        raise Exception("ERROR: failed to get user ID from session_id {0}".format(session_id))
    else:
        domain_carts = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_domain())
        user_carts   = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_by_user(user=current_user))
        group_carts  = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_by_user_groups())
        shared_carts = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_by_share_ids(share_ids=[share_id]))
        public_carts = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_public())

        result = { 'domain_carts':domain_carts,
                   'group_carts':group_carts,
                   'public_carts':public_carts,
                   'shared_carts':shared_carts,
                   'user_carts':user_carts }

    # Doing this so nested objects don't get stringified: https://stackoverflow.com/a/68935297
    print(json.dumps(result, default=lambda o: o.__dict__))

def filter_any_previous(ids, new_carts):
    carts = []

    for cart in new_carts:
        if cart.id not in ids:
            carts.append(cart)
            ids.add(cart.id)

    return carts
    # fails to return rows when I try to sort.
    #return carts.sort(key=lambda a: a.__dict__['label'].lower())
    
if __name__ == '__main__':
    main()
