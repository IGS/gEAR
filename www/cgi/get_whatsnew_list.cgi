#!/opt/bin/python3

"""
Gets new items we want to show in the "What's New" area of the home page
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

DATASETS_TO_SHOW = 3

def main():

    result = {
        'success': 1,
        'new_items': None,
    }

    cnx = geardb.Connection()
    datasets = geardb.DatasetCollection()
    datasets.get_public(has_h5ad=1, n=DATASETS_TO_SHOW, order_by='date_added')
    result['new_items'] = datasets

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

