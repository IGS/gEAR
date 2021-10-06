#!/opt/bin/python3

"""
Removes a gene cart via the 'delete' button in the gene cart manager
This requires the user to own the gene cart in order to remove it from the database.

This script first checks if the user owns the cart, then proceeds with the removal
if they own it.

If the user does not own the cart, an error is returned stating that.

Requires:
1) Session id - which contains user_id
2) cart id to be removed/deleted

"""

import cgi, json
import sys, os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    gc_id = form.getvalue('gene_cart_id')

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    # Does user own the cart...
    owns_gc = check_cart_ownership(cursor, current_user_id, gc_id)

    if owns_gc == True:
        result = { 'success': 1, 'gc':[] }

        cart = geardb.GeneCart(id=gc_id)
        cart.remove()

        print(json.dumps(result))

    else:
        error = "Not able to remove dataset. User does not own the dataset."
        result = { 'success': 0, 'error': error }
        print(json.dumps(result))

    cursor.close()
    cnx.close()


def check_cart_ownership(cursor, current_user_id, gc_id):
    qry = """
       SELECT id, user_id
       FROM gene_cart
       WHERE id = %s
    """
    cursor.execute(qry, (gc_id,))

    for row in cursor:
        # Change access if user owns it
        if row[1] == current_user_id:
            return True

    return False

def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
