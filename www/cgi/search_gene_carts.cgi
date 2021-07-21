#!/opt/bin/python3

"""
Used by gene_cart_manager.html, this script focuses on searching gene carts and
returns a list of matches with extended attributes.
"""

import cgi
import json
import os, sys
import re

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# limits the number of matches returned
DEFAULT_MAX_RESULTS = 200;
DEBUG_MODE = False

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    custom_list = form.getvalue('custom_list')
    search_terms = form.getvalue('search_terms').split(' ') if form.getvalue('search_terms') else []
    date_added = form.getvalue('date_added')
    ownership = form.getvalue('ownership')
    sort_by = re.sub("[^[a-z]]", "", form.getvalue('sort_by'))
    user = geardb.get_user_from_session_id(session_id) if session_id else None
    result = {'success': 0, 'problem': '', 'gene_carts': []}

    gene_carts = list()
    qry_params = []

    selects = ["gc.id", "g.user_name", "gc.gctype", "gc.label", "gc.ldesc", "gc.share_id",
               "gc.is_public", "gc.date_added"]
    froms = ["gene_cart gc", "guser g"]
    wheres = [
        "gc.user_id = g.id "
    ]
    orders_by = []

    if not user:
        # user not logged in, they can only see public datasets
        wheres.append("AND gc.is_public = 1")
    else:
        # if any ownership filters are defined, collect those
        if ownership:
            # options are 'yours', 'shared', 'public'
            owners = ownership.split(',')
            ownership_bits = []

            if 'yours' in owners:
                ownership_bits.append("gc.user_id = %s")
                qry_params.append(user.id)

            if 'public' in owners:
                ownership_bits.append("gc.is_public = 1")

            wheres.append("AND ({0})".format(' OR '.join(ownership_bits)))

        # otherwise, give the usual self and public.
        else:
            wheres.append("AND (gc.is_public = 1 OR gc.user_id = %s)")
            qry_params.extend([user.id])
        
        if search_terms:
            selects.append(' MATCH(gc.label, gc.ldesc) AGAINST("%s" IN BOOLEAN MODE) as rscore')
            wheres.append(' AND MATCH(gc.label, gc.ldesc) AGAINST("%s" IN BOOLEAN MODE)')

            # this is the only instance where a placeholder can be in the SELECT statement, so it will
            #  be the first qry param
            qry_params.insert(0, ' '.join(search_terms))
            qry_params.append(' '.join(search_terms))

        if date_added:
            date_added = re.sub("[^a-z]", "", date_added)
            wheres.append("AND gc.date_added BETWEEN date_sub(now(), INTERVAL 1 {0}) AND now()".format(date_added))

        if sort_by == 'relevance':
            # relevance can only be ordered if a search term was used
            if search_terms:
                orders_by.append(" rscore DESC")
            else:
                orders_by.append(" gc.date_added DESC")
        elif sort_by == 'title':
            orders_by.append(" gc.label")
        elif sort_by == 'owner':
            orders_by.append(" g.user_name")
        else:
            orders_by.append(" gc.date_added DESC")
            
        # build query
        qry = """
        SELECT {0}
         FROM {1}
         WHERE {2}
        ORDER BY {3}
        """.format(
            ", ".join(selects),
            ", ".join(froms),
            " ".join(wheres),
            " ".join(orders_by)
        )

    if DEBUG_MODE:
        ofh = open('/tmp/debug', 'wt')
        ofh.write("QRY:\n{0}\n".format(qry))
        ofh.write("QRY_params:\n{0}\n".format(qry_params))
        ofh.close()
        
    cursor.execute(qry, qry_params)

    for row in cursor:
        gc = geardb.GeneCart(id=row[0], gctype=row[2], label=row[3], ldesc=row[4], share_id=row[5],
                             is_public=row[6], date_added=row[7])
        gc.user_name = row[1]
        gc.gene_count = len(gc.genes)
        gene_carts.append(gc)

    result['gene_carts'] = gene_carts
    result['success'] = 1

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
