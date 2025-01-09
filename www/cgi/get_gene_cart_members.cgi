#!/opt/bin/python3

"""
Gets the gene symbols of a gene cart. Returns gene cart member ID's and Gene symbols
{'id': 1234, 'label': 'Sox9'}
"""

import cgi
import json
import sys
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2].joinpath('lib')
sys.path.insert(0, str(lib_path))

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

import geardb


def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)
    gene_cart_share_id = form.getvalue('share_id')
    result = { 'gene_symbols':[], 'success': 0 }

    if not gene_cart_share_id:
        raise Exception("ERROR: missing gene cart share ID")

    gc = geardb.get_gene_cart_by_share_id(gene_cart_share_id)
    if gc is None:
        raise Exception("ERROR: failed to get gene cart share ID {0}".format(gene_cart_share_id))

    if gc.gctype == "unweighted-list":
        gene_cart_query = """
            SELECT id, gene_symbol
                FROM gene_cart_member
                WHERE gene_cart_id = %s
        """

        cursor.execute(gene_cart_query, (gc.id,))
        for row in cursor:
            result['gene_symbols'].append({'id': row[0], 'label': row[1]})
    elif gc.gctype == "weighted-list":
        share_id = gc.share_id
        if not share_id:
            raise Exception("ERROR: weighted-list gene cart {0} has no share ID".format(gc.id))

        # Get the gene symbols from the shared cart file
        file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + share_id))
        import csv

        with open(file_path, "r") as fh:
            reader = csv.reader(fh, delimiter="\t")
            # Skip header
            next(reader)
            result["gene_symbols"] = [{"id":row[0], "label":row[1]} for row in reader]

    elif gc.gctype == "labeled-list":
        raise NotImplementedError("ERROR: labeled-list gene carts not yet implemented")
    else:
        raise Exception("ERROR: unknown gene cart type {0}".format(gc.gctype))

    result['success'] = 1
    cursor.close()
    cnx.close()

    print(json.dumps(result))

if __name__ == '__main__':
    main()
