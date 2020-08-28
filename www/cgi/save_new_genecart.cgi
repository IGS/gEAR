#!/opt/bin/python3

"""
Saves a GeneCart object (from genecart.js) serialized as a JSON string
and submitted to this script.  This assumes a NEW GeneCart object.

Input looks like this:

{"label": "Test", "session_id": "e385305c-4387-433e-8b62-4bcf7c30ac52", "genes": [{"id": 2804, "gene_symbol": "Otor"}, {"id": 2957, "gene_symbol": "Tgfbi"}, {"id": 15772, "gene_symbol": "H19"}, {"id": 168014, "gene_symbol": "Kctd12"}]}

"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    gc = geardb.GeneCart()
    gc.update_from_json(json.load(sys.stdin))
    gc.save()

    result = { 'id': gc.id }
    print(json.dumps(result))


if __name__ == '__main__':
    main()
