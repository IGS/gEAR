#!/opt/bin/python3

"""
Saves a GeneCart object (from genecart.js) serialized as 
formData and submitted to this script.  This assumes we're creating 
a NEW GeneCart object.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
from pprint import pprint
import geardb

def main():
    print('Content-Type: application/json\n\n')
    gc = geardb.GeneCart()

    ## else it's from a form submission
    form = cgi.FieldStorage()
    gc.update_from_form_data(form)
    gc.save()

    result = { 'id': gc.id }
    print(json.dumps(result))

if __name__ == '__main__':
    main()
