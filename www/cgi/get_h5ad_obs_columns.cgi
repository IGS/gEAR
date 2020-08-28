#!/opt/bin/python3

import cgi
import json
import os
import sys

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
import scanpy as sc
sc.settings.verbosity = 0

import pandas as pd

def main():
    result = { 'success': 0 }

    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
    h5_path = ds.get_file_path()

    # Let's not fail if the file isn't there
    if not os.path.exists(h5_path):
        result['success'] = 0
        result['error'] = "No h5 file found for this dataset"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    adata = sc.read_h5ad(h5_path)
    columns = adata.obs.columns.tolist()
    if 'replicate' in columns:
        columns.remove('replicate')

    result['obs_columns'] = columns

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
