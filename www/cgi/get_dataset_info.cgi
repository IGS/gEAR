#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    include_shape = form.getvalue('include_shape')

    ds = geardb.get_dataset_by_id(id=dataset_id, include_shape=include_shape)

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(ds)

if __name__ == '__main__':
    main()
