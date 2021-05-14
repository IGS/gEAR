#!/opt/bin/python3

"""
Gets a simple list of the organisms in the database

Data structure returned:

{
   organisms: [
    {
      dataset_id: "dataset12.corrected",
      title: "Cell-specific RNASeq in the ear",
      ldesc: "Super amazing illustration of cell-specific coloring of the ear cells.",
      access: "Private"
    }
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
   
    cursor.close()
    cnx.close()

    # does a recursive dump
    # https://yzhong-cs.medium.com/serialize-and-deserialize-complex-json-in-python-205ecc636caa
    print(json.dumps(result, default=lambda o: o.__dict__))

if __name__ == '__main__':
    main()
