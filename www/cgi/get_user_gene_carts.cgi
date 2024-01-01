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
    filter_cart_type = form.getvalue('cart_type', None)
    group_by_type = form.getvalue("group_by_type", False)
    current_user = geardb.get_user_from_session_id(session_id)

    result = { 'domain_carts':[], 'group_carts':[], 'public_carts':[],
               'shared_carts':[], 'user_carts':[] }

    # Track the cart IDs already stored so we don't duplicate
    cart_ids_found = set()

    domain_carts = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_domain())
    user_carts = []
    group_carts = []
    if current_user:
        user_carts   = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_by_user(user=current_user))
        group_carts  = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_by_user_groups(user=current_user))
    shared_carts = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_by_share_ids(share_ids=[share_id]))
    public_carts = filter_any_previous(cart_ids_found, geardb.GeneCartCollection().get_public())

    for carts in [domain_carts, user_carts, group_carts, shared_carts, public_carts]:
        for cart in carts:
            cart.label = f"{cart.label} ({cart.num_genes} genes)"

    if group_by_type and not group_by_type == "false":
        # Group all cart results by their cart type and return
        gctypes = ["unweighted-list", "weighted-list"]
        gctypes_result = dict()
        for cart_type in gctypes:
            subset_domain_carts = filter_by_cart_type(domain_carts, cart_type)
            subset_user_carts = filter_by_cart_type(user_carts, cart_type)
            subset_group_carts = filter_by_cart_type(group_carts, cart_type)
            subset_shared_carts = filter_by_cart_type(shared_carts, cart_type)
            subset_public_carts = filter_by_cart_type(public_carts, cart_type)

            gctypes_result[cart_type] = { 'domain_carts':subset_domain_carts,
                'group_carts':subset_group_carts,
                'public_carts':subset_public_carts,
                'shared_carts':subset_shared_carts,
                'user_carts':subset_user_carts }

        # Doing this so nested objects don't get stringified: https://stackoverflow.com/a/68935297
        print(json.dumps(gctypes_result, default=lambda o: o.__dict__))
        sys.exit(0)

    if len(str(filter_cart_type)) > 0 and filter_cart_type != 'undefined':
        domain_carts = filter_by_cart_type(domain_carts, filter_cart_type)
        user_carts = filter_by_cart_type(user_carts, filter_cart_type)
        group_carts = filter_by_cart_type(group_carts, filter_cart_type)
        shared_carts = filter_by_cart_type(shared_carts, filter_cart_type)
        public_carts = filter_by_cart_type(public_carts, filter_cart_type)

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

def filter_by_cart_type(carts, cart_type):
    """ Filters a list of carts by cart_type

    Args:
        carts (List[GeneCart]): A list of GeneCart objects
        cart_type (str): A string representing the cart type

    Returns:
        A list of GeneCart objects
    """
    if cart_type:
        return [cart for cart in carts if cart.gctype == cart_type]
    else:
        return carts

if __name__ == '__main__':
    main()
