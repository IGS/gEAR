#!/opt/bin/python3

"""
Changes the information of a gene cart via .is-editable elements in the Gene Cart Manager

This script first checks if the user owns the cart, then proceeds with the info change
if they own it.

If the user does not own the cart, an error is returned stating that.
Successful execution returns the cart in JSON format.

Requires:
1) session_id - used to get user_id
2) gene_cart_id
3) One of the following:
    - Title (title)
    - Organism (organism_id)
    - Long description (ldesc)
    - Public/private visibility (visibility)

"""

import os
import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    gc_id = form.getvalue('gc_id')
    visibility = int(form.getvalue('visibility'))
    organism_id = form.getvalue('organism_id')
    label = form.getvalue('title')
    ldesc = form.getvalue('ldesc')

    user = geardb.get_user_from_session_id(session_id)
    gc = geardb.get_gene_cart_by_id(gc_id)

    print("visibility:{0} gc.is_public:{1}".format(visibility, gc.is_public), file=sys.stderr)

    if user.id == gc.user_id:
        # see what has changed and execute updates to the DB
        if visibility == 0 and gc.is_public == 1:
            gc.save_change('is_public', 0)
        elif visibility == 1 and (gc.is_public == 0 or gc.is_public == None):
            gc.save_change('is_public', 1)

        if gc.label != label:
            gc.save_change('label', label)

        if gc.organism_id != organism_id:
            gc.save_change('organism_id', organism_id)

        if gc.ldesc != ldesc:
            gc.save_change('ldesc', ldesc)

        result = { 'gene_cart': gc, 'success': 1 }
            
        print(json.dumps(result))

    else:
        result = { 'error':[] }

        error = "Not able to change gene cart's information. You do not own this cart."
        result['error'] = error
        result['success'] = 0

        print(json.dumps(result))

if __name__ == '__main__':
    main()
