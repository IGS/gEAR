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
    session_id = form.getfirst('session_id')
    gc_id = form.getfirst('gc_id')
    visibility = form.getfirst('visibility')
    visibility = int(visibility or 0)
    organism_id = form.getfirst('organism_id')
    label = form.getfirst('title')
    ldesc = form.getfirst('ldesc')

    result = {}


    if session_id is None:
        error = "Not able to change gene cart's information. No session id provided."
        result['error'] = error
        result['success'] = 0

        print(json.dumps(result))
        return

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        error = "Not able to change gene cart's information. Invalid session id provided."
        result['error'] = error
        result['success'] = 0

        print(json.dumps(result))
        return

    gc = geardb.get_gene_cart_by_id(gc_id)

    if gc is None:
        error = "Not able to change gene cart's information. Invalid gene cart id provided."
        result['error'] = error
        result['success'] = 0

        print(json.dumps(result))
        return

    #print("visibility:{0} gc.is_public:{1}".format(visibility, gc.is_public), file=sys.stderr)

    if user.id == gc.user_id:
        # see what has changed and execute updates to the DB
        # ? SAdkins - Why are we checking for differences? Can't we just update regardless, or are we trying to reduce transactions?
        if gc.is_public != visibility:
            gc.save_change('is_public', visibility)

        if gc.label != label:
            gc.save_change('label', label)

        if gc.organism_id != organism_id:
            gc.save_change('organism_id', organism_id)

        if gc.ldesc != ldesc:
            gc.save_change('ldesc', ldesc)

        result = { 'gene_cart': gc, 'success': 1 }

        print(json.dumps(result))

    else:
        error = "Not able to change gene cart's information. You do not own this cart."
        result['error'] = error
        result['success'] = 0

        print(json.dumps(result))

if __name__ == '__main__':
    main()
