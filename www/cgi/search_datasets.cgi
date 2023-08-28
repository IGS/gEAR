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

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from gear.userhistory import UserHistory

# limits the number of matches returned
DEFAULT_MAX_RESULTS = 200;
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
    sort_by = re.sub("[^[a-z]]", "", form.getvalue('sort_by'))
    user = geardb.get_user_from_session_id(session_id) if session_id else None
    result = {'success': 0, 'problem': '', 'datasets': []}

    datasets_collection = geardb.DatasetCollection()
    qry_params = []
    shared_dataset_id_str = None

    if user:
        shared_dataset_id_str = get_shared_dataset_id_string(user, cursor)

    if custom_list:
        if custom_list == 'most_recent':
            if user:
                qry = """
                      SELECT d.id
                        FROM dataset d
                       WHERE d.marked_for_removal = 0
                         AND d.load_status = 'completed'
                         AND (d.owner_id = %s or d.is_public = 1)
                    ORDER BY d.date_added DESC LIMIT 10
                """
                qry_params = [user.id]
            else:
                qry = """
                      SELECT d.id
                        FROM dataset d
                       WHERE d.marked_for_removal = 0
                         AND d.load_status = 'completed'
                         AND d.is_public = 1
                    ORDER BY d.date_added DESC LIMIT 10
                """
                qry_params = []
        else:
            result['success'] = 0
            result['problem'] = "Didn't recognize the custom list requested"
    elif layout_share_id:
        qry = """
              SELECT d.id, lm.grid_position, lm.grid_width
                  FROM dataset d
                       JOIN layout_members lm ON d.id=lm.dataset_id
                       JOIN layout l ON lm.layout_id=l.id
                 WHERE d.marked_for_removal = 0
                   AND d.load_status = 'completed'
                   AND l.share_id = %s
                ORDER BY lm.grid_position
        """
        qry_params = [layout_share_id]
    else:
        selects = ["d.id", "g.user_name"]
        froms = ["dataset d", "guser g"]
        wheres = [
            "d.owner_id = g.id "
            "AND d.marked_for_removal = 0 ",
            "AND d.load_status = 'completed'"
        ]
        orders_by = []

        if not user:
            # user not logged in, they can only see public datasets
            wheres.append("AND d.is_public = 1")
        else:
            # if any ownership filters are defined, collect those
            if ownership:
                # options are 'yours', 'shared', 'public'
                owners = ownership.split(',')
                ownership_bits = []

                if 'yours' in owners:
                    ownership_bits.append("d.owner_id = %s")
                    qry_params.append(user.id)

                if 'public' in owners:
                    ownership_bits.append("d.is_public = 1")

                if 'shared' in owners:
                    if shared_dataset_id_str:
                        ownership_bits.append("d.id IN ({0})".format(shared_dataset_id_str))
                    else:
                        # Query gets thrown off if there are no shares for this user so just put
                        #  something in which won't match any.  This allows other params to still
                        #  match.
                        ownership_bits.append("d.id IN ('XYZ')".format(shared_dataset_id_str))

                # Don't add this if there aren't any
                wheres.append("AND ({0})".format(' OR '.join(ownership_bits)))

            # otherwise, give the usual self, public and shared.
            else:
                if shared_dataset_id_str:
                    wheres.append("AND (d.is_public = 1 OR d.owner_id = %s OR d.id IN ({0}))".format(shared_dataset_id_str))
                else:
                    wheres.append("AND (d.is_public = 1 OR d.owner_id = %s)")

                qry_params.extend([user.id])

        if search_terms:
            selects.append(' MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST("%s" IN BOOLEAN MODE) as rscore')
            wheres.append(' AND MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST("%s" IN BOOLEAN MODE)')

            # this is the only instance where a placeholder can be in the SELECT statement, so it will
            #  be the first qry param
            qry_params.insert(0, ' '.join(search_terms))
            qry_params.append(' '.join(search_terms))

        if organism_ids:
            ## only numeric characters and the comma are allowed here
            organism_ids = re.sub("[^,0-9]", "", organism_ids)
            wheres.append("AND d.organism_id in ({0})".format(organism_ids))

            # This way seemed safer but I couldn't get it working in time
            #organism_id_str = "AND d.organism_id IN (%s)" % repr(tuple(map(str,organism_ids.split(','))))
            #wheres.append(organism_id_str)

        if dtypes:
            ## only alphanumeric characters and the dash are allowed here
            dtypes = re.sub("[^,\-A-Za-z0-9]", "", dtypes).split(',')
            dtype_str = (', '.join('"' + item + '"' for item in dtypes))
            wheres.append("AND d.dtype in ({0})".format(dtype_str))

        if date_added:
            date_added = re.sub("[^a-z]", "", date_added)
            wheres.append("AND d.date_added BETWEEN date_sub(now(), INTERVAL 1 {0}) AND now()".format(date_added))

        if sort_by == 'relevance':
            # relevance can only be ordered if a search term was used
            if search_terms:
                orders_by.append(" rscore DESC")
            else:
                orders_by.append(" d.date_added DESC")
        elif sort_by == 'title':
            orders_by.append(" d.title")
        elif sort_by == 'owner':
            orders_by.append(" g.user_name")
        else:
            orders_by.append(" d.date_added DESC")

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
        ofh = open('/tmp/dataset.search.debug', 'wt')
        ofh.write("QRY:\n{0}\n".format(qry))
        ofh.write("QRY_params:\n{0}\n".format(qry_params))
        ofh.close()

    cursor.execute(qry, qry_params)

    matching_dataset_ids = list()
    # this index keeps track of the size and position of each dataset if a layout was passed
    layout_idx = dict()
    for row in cursor:
        matching_dataset_ids.append(row[0])

        if layout_share_id:
            layout_idx[row[0]] = {'position': row[1], 'width': row[2]}

    result['datasets'].extend(datasets_collection.get_by_dataset_ids(matching_dataset_ids))

    for dataset in result['datasets']:
        if layout_share_id:
            dataset.grid_position = layout_idx[dataset.id]['position']
            dataset.grid_width = layout_idx[dataset.id]['width']

        dataset.get_layouts(user=user)
        if os.path.exists("{0}/{1}.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        elif os.path.exists("{0}/{1}.single.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.single.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        elif os.path.exists("{0}/{1}.multi.default.png".format(IMAGE_ROOT, dataset.id)):
            dataset.preview_image_url = "{0}/{1}.multi.default.png".format(WEB_IMAGE_ROOT, dataset.id)
        else:
            dataset.preview_image_url = "{0}/missing.png".format(WEB_IMAGE_ROOT, dataset.id)

    if False:
       try:
           #### The relevent block above needs to be put into here once tested.
           pass
       except:
            result['success'] = 0
            result['problem'] = "There seems to be a database issue."

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

    # Log the search for the user
    if user:
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
