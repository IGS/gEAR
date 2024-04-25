#!/opt/bin/python3

"""
Can be used to export the SQL for a gene cart and all of its members. Probably
only useful if copying data from one server to another (production/devel, etc.)

Has an option to export with a different user ID, which is useful if changing
users when using SQL on a different server.

You need to also pass the next available cart ID on the destination server
so this can create INSERT sql statements with non-conflicting IDs

If your share IDs includ weighted carts, you'll need to migrate these manually,
but a list is provided at the end of this script.

Output:

Redirect the STDOUT to a file for SQL. The rest is printed on STDERR
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser( description='Export SQL for a gene cart and its members')
    parser.add_argument('-c', '--cart_share_ids', type=str, required=True, help='Comma-separated list of cart share_ids' )
    parser.add_argument('-n', '--next_cart_id', type=int, required=True, help='Next gene cart ID this script should use for insert statements' )
    parser.add_argument('-uid', '--user_id', type=int, required=False, help='Use this user ID in SQL instead of the source one' )
    args = parser.parse_args()

    cnx = geardb.Connection()
    cart_cursor = cnx.get_cursor()
    gcm_cursor = cnx.get_cursor()
    next_cart_id = args.next_cart_id

    cart_share_ids = args.cart_share_ids.split(',')

    cart_qry = """
          SELECT id, user_id, organism_id, gctype, label, ldesc, share_id, 
                 is_public, is_domain, date_added
            FROM gene_cart
           WHERE share_id = %s
    """

    gcm_qry = """
          SELECT id, gene_symbol
            FROM gene_cart_member
           WHERE gene_cart_id = %s
    """

    weighted_cart_files = list()

    for cart_share_id in cart_share_ids:
        cart_cursor.execute(cart_qry, [cart_share_id,])
        
        for (id, user_id, organism_id, gctype, label, ldesc, share_id, is_public, is_domain, date_added) in cart_cursor:
            if args.user_id:
                user_id = args.user_id

            if is_domain == 'None' or is_domain is None:
                is_domain = 'NULL'

            if is_public == 'None' or is_public is None:
                is_public = 'NULL'
            
            print("INSERT INTO gene_cart (id, user_id, organism_id, gctype, label, ldesc, share_id, is_public, is_domain, date_added) ")
            print("VALUES ({0}, {1}, {2}, '{3}', '{4}', '{5}', '{6}', {7}, {8}, '{9}');".format(
                next_cart_id, user_id, organism_id, gctype, cnx.mysql_cnx.converter.escape(label), cnx.mysql_cnx.converter.escape(ldesc), share_id, is_public, is_domain, date_added
            ))

            if gctype == 'unweighted-list':
                gcm_cursor.execute(gcm_qry, [id,])
                
                for (id, gene_symbol) in gcm_cursor:
                    print("INSERT INTO gene_cart_member (id, gene_cart_id, gene_symbol)")
                    print("VALUES ({0}, {1}, '{2}');".format(id, next_cart_id, cnx.mysql_cnx.converter.escape(gene_symbol)))

            elif gctype == 'weighted-list':
                weighted_cart_files.append("cart.{0}.tab cart.{0}.h5ad".format(share_id))

            next_cart_id += 1

    if len(weighted_cart_files) > 0:
        print("\nDon't forget to migrate the weighted cart files!:", file=sys.stderr)
        print(" ".join(weighted_cart_files), file=sys.stderr)

    cart_cursor.close()
    gcm_cursor.close()
    cnx.close()

def get_layout_string(cnx, dataset_id):
    cursor = cnx.get_cursor()
    labels = list()

    qry = """
          SELECT l.share_id, l.label
            FROM layout l 
                 JOIN layout_members lm ON lm.layout_id=l.id
           WHERE lm.dataset_id = %s
    """
    cursor.execute(qry, [dataset_id,])

    for (share_id, label) in cursor:
        labels.append(label)

    cursor.close()

    if len(labels):
        return " / ".join(labels)
    else:
        return "None"

    
    
if __name__ == '__main__':
    main()
    
