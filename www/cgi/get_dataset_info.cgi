#!/opt/bin/python3

"""
get_dataset_info.cgi - Return metadata for a single dataset.

Input: dataset_id, include_shape (optional; include matrix dimensions).
Output: JSON-serialized Dataset object (as printed by geardb.Dataset).
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    dataset_id = form.getfirst('dataset_id')
    include_shape = form.getfirst('include_shape')

    ds = geardb.get_dataset_by_id(d_id=dataset_id, include_shape=include_shape)

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(ds)

if __name__ == '__main__':
    main()
