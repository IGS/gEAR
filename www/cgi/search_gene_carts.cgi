#!/opt/bin/python3

"""
Used by gene_cart_manager.html, this script focuses on searching gene carts and
returns a list of matches with extended attributes.
"""

import cgi
import json
import os, sys
import re
from math import ceil

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# limits the number of matches returned
DEFAULT_MAX_RESULTS = 200
DEBUG_MODE = False

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    custom_list = form.getvalue('custom_list')
    search_terms = form.getvalue('search_terms').split(' ') if form.getvalue('search_terms') else []
    organism_ids = form.getvalue('organism_ids')
    date_added = form.getvalue('date_added')
    ownership = form.getvalue('ownership')
    page = form.getvalue('page')    # page starts at 1
    limit = form.getvalue('limit')
    sort_by = re.sub("[^[a-z]]", "", form.getvalue('sort_by'))
    user = geardb.get_user_from_session_id(session_id) if session_id else None
    result = {'success': 0, 'problem': '', 'gene_carts': []}

    if page and not page.isdigit():
        raise ValueError("Page must be a number")

    if page and int(page) < 1:
        raise ValueError("Page must be greater than 0")

    if limit and not limit.isdigit():
        raise ValueError("Limit must be a number")

    if limit and int(limit) < 1:
        raise ValueError("Limit must be greater than 0")


    gene_carts = list()
    qry_params = []

    selects = ["gc.id", "g.user_name", "gc.gctype", "gc.label", "gc.ldesc", "gc.share_id",
               "gc.is_public", "gc.date_added", "o.genus", "o.species", "gc.organism_id",
               "gc.user_id"]
    froms = ["gene_cart gc", "guser g", "organism o"]
    wheres = [
        "gc.user_id = g.id ",
        "AND gc.organism_id = o.id "
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

        if organism_ids:
            ## only numeric characters and the comma are allowed here
            organism_ids = re.sub("[^,0-9]", "", organism_ids)
            wheres.append("AND gc.organism_id in ({0})".format(organism_ids))

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

        # if a limit is defined, add it to the query
        if int(limit):
            qry += " LIMIT {0}".format(limit)

        # if a page is defined, add it to the query
        if int(page):
            offset = int(page) - 1
            qry += " OFFSET {0}".format(offset * int(limit))

    if DEBUG_MODE:
        ofh = open('/tmp/debug', 'wt')
        ofh.write("QRY:\n{0}\n".format(qry))
        ofh.write("QRY_params:\n{0}\n".format(qry_params))
        ofh.close()

    cursor.execute(qry, qry_params)

    for row in cursor:
        gc = geardb.GeneCart(id=row[0], gctype=row[2], label=row[3], ldesc=row[4], share_id=row[5],
                             is_public=row[6], date_added=row[7], organism_id=row[10], user_id=row[11])
        gc.user_name = row[1]
        gc.gene_count = len(gc.genes)
        gc.organism = "{0} {1}".format(row[8], row[9])
        gene_carts.append(gc)

    # Get count of total results
    qry_count = """
        SELECT COUNT(*)
        FROM {0}
        WHERE {1}
        """.format(
            ", ".join(froms),
            " ".join(wheres)
        )

    # if search terms are defined, remove first qry_param (since it's in the SELECT statement)
    if search_terms:
        qry_params.pop(0)

    cursor.execute(qry_count, qry_params)

    # compile pagination information
    result["pagination"] = {}
    result["pagination"]['total_results'] = cursor.fetchone()[0]
    result["pagination"]['current_page'] = page if page else 1
    result["pagination"]['limit'] = limit if limit else DEFAULT_MAX_RESULTS
    result["pagination"]["total_pages"] = ceil(int(result["pagination"]['total_results']) / int(result["pagination"]['limit']))
    result["pagination"]["next_page"] = int(result["pagination"]['current_page']) + 1 if int(result["pagination"]['current_page']) < int(result["pagination"]['total_pages']) else None
    result["pagination"]["prev_page"] = int(result["pagination"]['current_page']) - 1 if int(result["pagination"]['current_page']) > 1 else None

    result['gene_carts'] = gene_carts
    result['success'] = 1

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
