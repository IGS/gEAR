#!/opt/bin/python3

import cgi, json
import sys
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from geardb import Dataset

import scanpy as sc
sc.settings.verbosity = 0

def main():
    form = cgi.FieldStorage()
    projection_source = form.getvalue('projection_source')
    set_of_patterns = form.getvalue('set_of_patterns')
    source_dataset_id = form.getvalue('source_dataset_id')
    target_dataset_id = form.getvalue('target_dataset_id')

    dataset = Dataset(id=dataset_id, has_h5ad=1)
    h5_path = dataset.get_file_path()


    if not os.path.exists(h5_path):
        result = dict()
        result['success'] = 0
        result['error'] = "No h5 file found for this dataset"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    adata = sc.read(h5_path)

    result = {
               'success': 1,
    }

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
