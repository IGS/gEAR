#!/opt/bin/python3

"""
Saves a GeneCart object (from genecart.js) serialized as 
formData and submitted to this script.  This assumes we're creating 
a NEW GeneCart object

This script is appropriate to call when passing JSON directly.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    gc = geardb.GeneCart()
    data = json.load(sys.stdin)
    gc.update_from_json(data)
    gc.save()

    result = { 'id': gc.id }
    print(json.dumps(result))

if __name__ == '__main__':
    main()
