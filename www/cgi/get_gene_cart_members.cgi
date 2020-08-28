#!/opt/bin/python3

"""
Gets the gene symbols of a gene cart. Returns gene cart member ID's and Gene symbols 
{'id': 1234, 'label': 'Sox9'}
"""

import cgi
import json
import sys

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
    user = geardb.get_user_from_session_id(session_id)
    gene_cart_id = form.getvalue('gene_cart_id')
    result = { 'gene_symbols':[], 'success': 0 }

    # Does the user have a current, saved layout?
    layout_id = None

    if user is None:
        raise Exception("ERROR: failed to get user ID from session_id {0}".format(session_id))
    else:
        gene_cart_query = """
            SELECT id, gene_symbol
              FROM gene_cart_member
             WHERE gene_cart_id = %s
        """

        cursor.execute(gene_cart_query, (gene_cart_id,))
        for row in cursor:
            result['gene_symbols'].append({'id': row[0], 'label': row[1]})

    result['success'] = 1
    cursor.close()
    cnx.close()

    print(json.dumps(result))

if __name__ == '__main__':
    main()
