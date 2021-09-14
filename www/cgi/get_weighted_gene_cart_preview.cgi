#!/opt/bin/python3

"""

"""

import cgi
import json
import string
import sys

import pandas as pd
import scanpy as sc

import os, sys

gene_cart_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'carts')
ROWS_TO_SHOW = 5

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    share_id = form.getvalue('share_id')
    valid_chars = "%s%s" % (string.ascii_letters, string.digits)
    share_id = ''.join(c for c in share_id if c in valid_chars)
    
    result = { 'preview_json':[], 'success': 0 }

    cart_file_path = os.path.join(gene_cart_path, "cart.{0}.h5ad".format(share_id))

    adata = sc.read_h5ad(cart_file_path, backed="r")
    adata.var_names_make_unique()

    df = pd.DataFrame(adata.X[:,0:ROWS_TO_SHOW], adata.obs.index, adata.var.index[0:ROWS_TO_SHOW]).transpose()
    result['preview_json'] = df.to_html(classes=['weighted-list'])
    result['success'] = 1
    print(json.dumps(result))

if __name__ == '__main__':
    main()
