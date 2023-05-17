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
IMAGE_ROOT = os.path.abspath(os.path.join('..', 'img', 'dataset_previews'))
WEB_IMAGE_ROOT = './img/dataset_previews'

def main():
    new_items = list()
    
    result = {
        'success': 1,
        'new_items': None,
    }

    cnx = geardb.Connection()
    dc = geardb.DatasetCollection()
    dc.get_public(has_h5ad=1, n=DATASETS_TO_SHOW, order_by='date_added')

    for dataset in dc.datasets:
        if os.path.exists("{0}/{1}.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        elif os.path.exists("{0}/{1}.single.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.single.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        elif os.path.exists("{0}/{1}.multi.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.multi.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        else:
            dataset.preview_image_url = "{0}/missing.png".format(WEB_IMAGE_ROOT, dataset.id)

        new_items.append(dataset)
    
    result['new_items'] = new_items
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

