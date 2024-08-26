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
         'recent_carts': [
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

- Recent carts are any the user recently created.

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

GENE_LIST_TYPES = ["unweighted-list", "weighted-list", "labeled-list"]

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_id = form.getvalue('share_id')
    filter_cart_type = form.getvalue('cart_type', None)
    group_by_type = form.getvalue("group_by_type", False)
    include_members = form.getvalue("include_members", 1)
    current_user = geardb.get_user_from_session_id(session_id)

    result = { 'domain_carts':[], 'group_carts':[], 'public_carts':[],
               'recent_carts':[], 'shared_carts':[], 'user_carts':[] }

    # Track the cart IDs already stored so we don't duplicate
    cart_ids_found = set()

    if include_members:
        include_members = int(include_members)

    bool_include_members = False
    if include_members == 1:
        bool_include_members = True

    domain_carts = geardb.GeneCartCollection(include_genes=bool_include_members).get_domain()
    user_carts = []
    group_carts = []
    recent_carts = []
    if current_user:
        user_carts = geardb.GeneCartCollection(include_genes=bool_include_members).get_by_user(user=current_user)
        group_carts = geardb.GeneCartCollection(include_genes=bool_include_members).get_by_user_groups(user=current_user)
        recent_carts = geardb.GeneCartCollection(include_genes=bool_include_members).get_by_user_recent(user=current_user, n=10)
    shared_carts = geardb.GeneCartCollection(include_genes=bool_include_members).get_by_share_ids(share_ids=[share_id])
    public_carts = geardb.GeneCartCollection(include_genes=bool_include_members).get_public()

    if group_by_type and not group_by_type == "false":
        # Group all cart results by their cart type and return
        gctypes_result = dict()
        for cart_type in GENE_LIST_TYPES:
            subset_domain_carts = filter_by_cart_type(domain_carts, cart_type)
            subset_user_carts = filter_by_cart_type(user_carts, cart_type)
            subset_group_carts = filter_by_cart_type(group_carts, cart_type)
            subset_recent_carts = filter_by_cart_type(recent_carts, cart_type)
            subset_shared_carts = filter_by_cart_type(shared_carts, cart_type)
            subset_public_carts = filter_by_cart_type(public_carts, cart_type)

            gctypes_result[cart_type] = { 'domain_carts':subset_domain_carts,
                'group_carts':subset_group_carts,
                'public_carts':subset_public_carts,
                'recent_carts':subset_recent_carts,
                'shared_carts':subset_shared_carts,
                'user_carts':subset_user_carts }

        # Doing this so nested objects don't get stringified: https://stackoverflow.com/a/68935297
        print(json.dumps(gctypes_result, default=lambda o: o.__dict__))
        sys.exit(0)

    if filter_cart_type and filter_cart_type in GENE_LIST_TYPES:
        domain_carts = filter_by_cart_type(domain_carts, filter_cart_type)
        user_carts = filter_by_cart_type(user_carts, filter_cart_type)
        group_carts = filter_by_cart_type(group_carts, filter_cart_type)
        recent_carts = filter_by_cart_type(recent_carts, filter_cart_type)
        shared_carts = filter_by_cart_type(shared_carts, filter_cart_type)
        public_carts = filter_by_cart_type(public_carts, filter_cart_type)
    elif filter_cart_type and filter_cart_type not in GENE_LIST_TYPES:
        # log warning that the cart type is not valid and return all carts
        print("WARNING: Invalid cart type: " + filter_cart_type, file=sys.stderr)

    result = { 'domain_carts':domain_carts,
                'group_carts':group_carts,
                'public_carts':public_carts,
                'recent_carts':recent_carts,
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
