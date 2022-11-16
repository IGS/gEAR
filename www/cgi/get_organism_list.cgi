#!/opt/bin/python3

"""
Gets a simple list of the organisms in the database

Data structure returned:

{
   organisms: [
      {"id": 5, 
       "label": "Chicken", 
       "genus": "Gallus", 
       "species": "gallus",  
       "strain": null, 
       "taxon_id": 9031,
       "has_icon": 1
      }, ...
   ]
}

"""

import json
import os, sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    result = { 'organisms':[] }

    org_collection = geardb.OrganismCollection()
    result['organisms'] = org_collection.get_all()

    for org in result['organisms']:
        icon_dir = os.path.join(lib_path, '..', 'www', 'img', 'icons')
        icon_path = os.path.join(icon_dir, "org-{0}-dark-64.svg".format(org.id))

        if os.path.exists(icon_path):
            org.has_icon = 1
        else:
            org.has_icon = 0
   
    cursor.close()
    cnx.close()

    # does a recursive dump
    # https://yzhong-cs.medium.com/serialize-and-deserialize-complex-json-in-python-205ecc636caa
    print(json.dumps(result, default=lambda o: o.__dict__))

if __name__ == '__main__':
    main()
