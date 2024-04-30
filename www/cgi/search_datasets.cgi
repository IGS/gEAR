#!/opt/bin/python3

"""
Used by by dataset_explorer.html, this script focuses on searching datasets and
returns a list of matches with extended attributes.

There are a few main categories of search:

1. A custom-named search, given by the custom_list command
2. A search with no terms, where we're filtering on attributes instead
3. A search with terms and optionally attributes on which to filter
4. A listing of datasets in a single layout
"""

import cgi
import json
import os, sys
import re
from math import ceil

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.userhistory import UserHistory

# limits the number of matches returned
DEFAULT_MAX_RESULTS = 20
IMAGE_ROOT = os.path.abspath(os.path.join('..', 'img', 'dataset_previews'))
WEB_IMAGE_ROOT = './img/dataset_previews'
DEBUG_MODE = False

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    custom_list = form.getvalue('custom_list')
    search_terms = form.getvalue('search_terms').split(' ') if form.getvalue('search_terms') else []
    organism_ids = form.getvalue('organism_ids')
    dtypes = form.getvalue('dtypes')
    date_added = form.getvalue('date_added')
    ownership = form.getvalue('ownership')
    layout_share_id = form.getvalue('layout_share_id')
    page = form.getvalue('page', "1")    # page starts at 1
    limit = form.getvalue('limit', str(DEFAULT_MAX_RESULTS))
    sort_by = re.sub("[^[a-z]]", "", form.getvalue('sort_by'))
    user = geardb.get_user_from_session_id(session_id) if session_id else None
    result = {'success': 0, 'problem': '', 'datasets': []}

    if page and not page.isdigit():
        raise ValueError("Page must be a number")

    if page and int(page) < 1:
        raise ValueError("Page must be greater than 0")

    if limit and not limit.isdigit():
        raise ValueError("Limit must be a number")

    if limit and int(limit) < 1:
        raise ValueError("Limit must be greater than 0")

    datasets_collection = geardb.DatasetCollection()
    shared_dataset_id_str = None

    qry_params = []
    selects = ["d.id"]
    froms = ["dataset d"]
    wheres = ["d.marked_for_removal = 0", "d.load_status = 'completed'"]
    orders_by = []

    if user:
        shared_dataset_id_str = get_shared_dataset_id_string(user, cursor)
        ownership_bits = []

        # if any ownership filters are defined, collect those
        if ownership:
            # options are 'yours', 'shared', 'public'
            owners = ownership.split(',')

            if 'yours' in owners:
                ownership_bits.append("d.owner_id = %s")
                qry_params.append(user.id)

            if 'public' in owners:
                ownership_bits.append("d.is_public = 1")

            if 'shared' in owners and shared_dataset_id_str:
                ownership_bits.append(f"d.id IN ({shared_dataset_id_str})")
            elif 'shared' in owners:
                # Query gets thrown off if there are no shares for this user so just put
                #  something in which won't match any.  This allows other params to still
                #  match.
                ownership_bits.append("d.id IN ('XYZ')")

        else:
            ownership_bits.append("d.is_public = 1")
            ownership_bits.append("d.owner_id = %s")
            qry_params.append(user.id)
            if shared_dataset_id_str:
                ownership_bits.append(f"d.id IN ({shared_dataset_id_str})")

        wheres.append(f"({' OR '.join(ownership_bits)})")   # OR accomodates the "not ownership" case

    else:
        wheres.append("d.is_public = 1")

    if search_terms:
        selects.append('MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST("%s" IN BOOLEAN MODE) as rscore')
        wheres.append('MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST("%s" IN BOOLEAN MODE)')

        # this is the only instance where a placeholder can be in the SELECT statement, so it will
        #  be the first qry param
        qry_params.insert(0, ' '.join(search_terms))
        qry_params.append(' '.join(search_terms))

    if organism_ids:
        ## only numeric characters and the comma are allowed here
        organism_ids = re.sub("[^,0-9]", "", organism_ids)
        wheres.append("d.organism_id in ({0})".format(organism_ids))

        # This way seemed safer but I couldn't get it working in time
        #organism_id_str = "AND d.organism_id IN (%s)" % repr(tuple(map(str,organism_ids.split(','))))
        #wheres.append(organism_id_str)

    if dtypes:
        ## only alphanumeric characters and the dash are allowed here
        dtypes = re.sub("[^,\-A-Za-z0-9]", "", dtypes).split(',')
        dtype_str = (', '.join('"' + item + '"' for item in dtypes))
        wheres.append(f"d.dtype in ({dtype_str})")

    if date_added:
        date_added = re.sub("[^a-z]", "", date_added)
        wheres.append(f"d.date_added BETWEEN date_sub(now(), INTERVAL 1 {date_added}) AND now()")

    if sort_by == 'relevance':
        # relevance can only be ordered if a search term was used
        if search_terms:
            orders_by.append(" rscore DESC")
        else:
            orders_by.append(" d.date_added DESC")
    elif sort_by == 'title':
        orders_by.append(" d.title")
    elif sort_by == 'owner':
        selects.append("g.user_name")
        froms.append("guser g")
        orders_by.append(" g.user_name")
    else:
        orders_by.append(" d.date_added DESC")

    qry = ""

    if custom_list:
        if custom_list == 'most_recent':
            orders_by.append('d.date_added DESC')
        else:
            result['success'] = 0
            result['problem'] = "Didn't recognize the custom list requested"
    elif layout_share_id:
        qry_params.append(layout_share_id)

        selects.extend(["lm.grid_position", "lm.start_col", "lm.grid_width", "lm.start_row", "lm.grid_height",
                        "lm.mg_grid_position", "lm.mg_start_col", "lm.mg_grid_width", "lm.mg_start_row", "lm.mg_grid_height"])
        froms.extend(["layout_members lm", "layout l"])
        wheres.extend([
            "d.id = lm.dataset_id",
            "lm.layout_id = l.id",
            "l.share_id = %s"
        ])

        if not len(orders_by):
            orders_by.append("lm.grid_position ASC")

    # build query
    qry = f"""
    SELECT {', '.join(selects)}
        FROM {', '.join(froms)}
        WHERE {' AND '.join(wheres)}
        ORDER BY {' ,'.join(orders_by) if orders_by else 'd.date_added DESC'}
    """

    # if a limit is defined, add it to the query
    if int(limit):
        qry += f" LIMIT {limit}"

    # if a page is defined, add it to the query
    if int(page):
        offset = int(page) - 1
        qry += " OFFSET {0}".format(offset * int(limit))

    if DEBUG_MODE:
        ofh = open('/tmp/dataset.search.debug', 'wt')
        ofh.write(f"QRY:\n{qry}\n")
        ofh.write(f"QRY_params:\n{qry_params}\n")
        ofh.close()

    cursor.execute(qry, qry_params)

    matching_dataset_ids = list()
    # this index keeps track of the size and position of each dataset if a layout was passed
    layout_idx = dict()
    for row in cursor:
        matching_dataset_ids.append(row[0])

        if layout_share_id:
            layout_idx[row[0]] = {'position': row[1], 'start_row': row[2], 'width': row[3], 'start_col': row[4], 'grid_height': row[5],
                                  'mg_grid_position': row[6], 'mg_start_row': row[7], 'mg_grid_width': row[8], 'mg_start_col': row[9], 'mg_grid_height': row[10]
                                 }

    result['datasets'] = datasets_collection.get_by_dataset_ids(matching_dataset_ids)

    for dataset in result['datasets']:
        if layout_share_id:
            dataset.grid_position = layout_idx[dataset.id]['position']
            dataset.start_col = layout_idx[dataset.id]['start_col']
            dataset.grid_width = layout_idx[dataset.id]['width']
            dataset.start_row = layout_idx[dataset.id]['start_row']
            dataset.grid_height = layout_idx[dataset.id]['grid_height']
            dataset.mg_grid_position = layout_idx[dataset.id]['mg_grid_position']
            dataset.mg_start_col = layout_idx[dataset.id]['mg_start_col']
            dataset.mg_grid_width = layout_idx[dataset.id]['mg_grid_width']
            dataset.mg_start_row = layout_idx[dataset.id]['mg_start_row']
            dataset.mg_grid_height = layout_idx[dataset.id]['mg_grid_height']

        dataset.get_layouts(user=user)
        # delete user_id and layout_id from each dataset layout (security reasons)
        for l in dataset.layouts:
            #l.is_owner = True if user and l.user_id == user.id else False
            del l.user_id
            del l.id

        if os.path.exists("{0}/{1}.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        elif os.path.exists("{0}/{1}.single.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.single.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        else:
            dataset.preview_image_url = "{0}/missing.png".format(WEB_IMAGE_ROOT, dataset.id)

        # Multi-gene preview image
        if os.path.exists("{0}/{1}.multi.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.mg_preview_image_url = "{0}/{1}.multi.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        else:
            dataset.mg_preview_image_url = "{0}/missing.png".format(WEB_IMAGE_ROOT, dataset.id)

        # add if the user is the owner of the dataset
        dataset.is_owner = True if user and dataset.owner_id == user.id else False

    # Get count of total results
    qry_count = f"""
        SELECT COUNT(*)
        FROM {', '.join(froms)}
        WHERE {' AND '.join(wheres)}
        """

    # if search terms are defined, remove first qry_param (since it's in the SELECT statement)
    if search_terms:
        qry_params.pop(0)

    cursor.execute(qry_count, qry_params)

    # compile pagination information
    result["pagination"] = {}
    result["pagination"]['total_results'] = cursor.fetchone()[0]
    result["pagination"]['current_page'] = int(page)
    result["pagination"]['limit'] = int(limit)
    result["pagination"]["total_pages"] = ceil(int(result["pagination"]['total_results']) / int(result["pagination"]['limit']))
    result["pagination"]["next_page"] = int(result["pagination"]['current_page']) + 1 if int(result["pagination"]['current_page']) < int(result["pagination"]['total_pages']) else None
    result["pagination"]["prev_page"] = int(result["pagination"]['current_page']) - 1 if int(result["pagination"]['current_page']) > 1 else None


    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

    # Log the search for the user
    if user and len(search_terms):
        history = UserHistory()
        history.add_record(
            user_id=user.id,
            entry_category='dataset_search',
            label="Datasets matching '{0}'".format(' '.join(search_terms)),
            search_terms=search_terms,
        )

def get_shared_dataset_id_string(user, cursor):
    qry = "SELECT dataset_id FROM dataset_shares WHERE is_allowed = 1 AND user_id = %s"
    cursor.execute(qry, [user.id,])

    dataset_ids = list()
    for row in cursor:
        dataset_ids.append(row[0])

    return (', '.join('"' + id + '"' for id in dataset_ids))


if __name__ == '__main__':
    main()
