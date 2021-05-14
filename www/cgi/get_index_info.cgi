#!/opt/bin/python3

"""

"""

import cgi, json
import os, sys, re

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

NUM_NEWEST_DATASETS_TO_SHOW = 3

def main():

    result = {
        'success': 1,
        'dataset_count': 0,
        'user_count': 0,
        'newest_public_datasets': None
    }

    cnx = geardb.Connection()

    public_collection = geardb.DatasetCollection()
    public_collection.get_public(has_h5ad=1, n=NUM_NEWEST_DATASETS_TO_SHOW, order_by='date_added')
    result['newest_public_datasets'] = public_collection

    result['dataset_count'] = geardb.get_dataset_count()
    result['user_count'] = geardb.get_user_count()
    
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

