#!/opt/bin/python3

"""
For a given session_id, returns data on the user's datasets.  Order of priority:

0.  User logged in with a layout passed to this script
1.  User logged in with a submitted search string
    a.  For searching their own datasets
    b.  For searching the open datasets of others
2.  User logged in with current, saved layout
3.  Default layout (domain-specific, if cookie set) + user's private datasets
4.  Default layout (Hearing domain, if cookie not set) + user's private datasets
5.  Default layout (domain-specific, if cookie set) or anonymous user
6.  Default layout (Hearing domain, if cookie not set) or anonymous user

Data structure returned:

{
   datasets: [
    {
      dataset_id: "dataset12.corrected",
      title: "Cell-specific RNASeq in the ear",
      ldesc: "Super amazing illustration of cell-specific coloring of the ear cells.",
      access: "Private"
    }
   ]
}

"""

import cgi, json
from datetime import datetime
from operator import itemgetter

import os, sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb

def main():
    cnx = geardb.Connection()

    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    scope = form.getvalue('scope')
    search_terms = form.getvalue('search_terms')

    # temporarily dealing with https://github.com/jorvis/gEAR/issues/350
    if search_terms is not None:
        search_terms = search_terms.translate(str.maketrans('','','+-/@'))

    permalink_id = form.getvalue('permalink_share_id')
    only_types_str = form.getvalue('only_types')
    sort_order = form.getvalue('order')
    default_domain_label = form.getvalue('default_domain')

    layout_id = None
    only_types = None

    if only_types_str:
        only_types = only_types_str.replace(' ', '').split(',')

    if sort_order is None:
        sort_order = 'default'

    current_user_id = get_user_id_from_session_id(cursor, session_id)
    result = { 'datasets':[] }

    # only used to non-redundify
    dataset_ids = list()

    # Permalinks only. Get dataset info and return it
    if permalink_id is not None:
        result['datasets'] = get_permalink_dataset(cursor, permalink_id)
        cursor.close()
        cnx.close()
        print(json.dumps(result))
        return

    if form.getvalue(("layout_share_id")) is not None:
        layout_share_id = form.getvalue('layout_share_id')
        layouts = geardb.LayoutCollection().get_by_share_id(layout_share_id)
        if len(layouts) > 2:
            raise Exception("More than one layout found with ID {0}".format(layout_id))
        if not len(layouts):
            raise Exception("No layout found with ID {0}".format(layout_id))
        layout = layouts[0]
        layout.load()

        dsc = geardb.DatasetCollection()
        dsc.get_by_dataset_ids(ids=layout.dataset_ids(), get_links=True)
        dsc.apply_layout(layout=layout)
        result['datasets'].extend(dsc.datasets)

    # Was a specific layout ID passed?
    if form.getvalue('layout_id') is not None:
        layout_id = form.getvalue('layout_id')
        layout = geardb.Layout(id=layout_id)
        layout.load()

        dsc = geardb.DatasetCollection()
        dsc.get_by_dataset_ids(ids=layout.dataset_ids(), get_links=True)
        dsc.apply_layout(layout=layout)
        result['datasets'].extend(dsc.datasets)

    # If scope is defined, the user is performing a search
    elif scope is not None:
        # If no search terms were defined, we want the whole list
        search_term_qry = ''
        qry_params = [current_user_id]

        # Search terms defined, so search for matching datasets
        if search_terms is not None:
            # > = Include word, and increase rank if found
            search_term_qry = '''   AND MATCH(d.title, d.ldesc, d.geo_id) AGAINST( %s )
                ORDER BY MATCH(d.title, d.ldesc) AGAINST(%s IN BOOLEAN MODE) DESC
            '''
            qry_params.append(search_terms)
            qry_params.append(search_terms)

        matching_dataset_ids = list()
        if scope == 'others':
            query = """
              SELECT d.id
                FROM dataset d
               WHERE d.is_public = 1
                 AND d.owner_id != %s
            """ + search_term_qry
        elif scope == 'self':
            #include shared datasets in search/list all
            query = """
              SELECT d.id
                FROM dataset d
               WHERE d.owner_id = %s
            """ + search_term_qry
            query_shared = """
              SELECT s.dataset_id as id
              FROM dataset_shares s
              JOIN dataset d ON d.id=s.dataset_id
              WHERE s.user_id = %s
                AND s.is_allowed = 1
            """ + search_term_qry
        elif scope == 'shared':
            query = """
              SELECT s.dataset_id as id
              FROM dataset_shares s
              JOIN dataset d ON d.id=s.dataset_id
              WHERE s.user_id = %s
                AND s.is_allowed = 1
            """ + search_term_qry
        elif scope == 'user_all':
            # Gathers all datasets the user has access to
            # Targeted function: compare tool - show no pending datasets
            query = """
              SELECT d.id
                FROM dataset d
               WHERE d.owner_id = %s
                AND load_status = 'completed'
            """ + search_term_qry
            query_shared = """
              SELECT s.dataset_id as id
              FROM dataset_shares s
              JOIN dataset d ON d.id=s.dataset_id
              WHERE s.user_id = %s
                AND s.is_allowed = 1
                AND load_status = 'completed'
            """ + search_term_qry
            query_public = """
              SELECT d.id
                FROM dataset d
               WHERE d.is_public = 1
                 AND d.owner_id != %s
                AND load_status = 'completed'
            """ + search_term_qry

        else:
            raise Exception("Dataset list requested but scope ({0}) wasn't recognized.".format(scope));

        try:
            cursor.execute(query, qry_params)
        except:
            print("The failed SQL was: {0}".format(cursor._executed), file=sys.stderr)
            raise Exception("The failed SQL was: {0}".format(cursor._executed))

        for row in cursor:
            matching_dataset_ids.append(row[0])

        #Search tags for any tagged datasets
        if search_terms is not None:
            search_term_list = search_terms.split()
            for term in search_term_list:
                qry_dataset_tag = """
                    SELECT d.dataset_id
                    FROM dataset_tag d
                    JOIN tag t ON t.id=d.tag_id
                    WHERE t.label=%s;
                """
                cursor.execute(qry_dataset_tag, (term,))
                for row in cursor:
                    if row[0] not in matching_dataset_ids:
                        matching_dataset_ids.append(row[0])

        if scope == 'self':
            # now gather shared datasets
            cursor.execute(query_shared, qry_params)
            for row in cursor:
                matching_dataset_ids.append(row[0])

        if scope == 'user_all':
            #now gather shared datasets
            cursor.execute(query_shared, qry_params)
            for row in cursor:
                matching_dataset_ids.append(row[0])

            #now gather public datasets
            cursor.execute(query_public, qry_params)
            for row in cursor:
                matching_dataset_ids.append(row[0])

        datasets_coll = geardb.DatasetCollection()
        result['datasets'].extend(
            datasets_coll.get_by_dataset_ids(ids=matching_dataset_ids, get_links=True)
        )

    else:
        # Do they have a current layout saved?
        saved_layout_query = "SELECT id FROM layout WHERE user_id = %s AND is_current = 1"
        cursor.execute(saved_layout_query, (current_user_id,))
        for row in cursor:
            layout_id = row[0]
            break

        # No layout saved, used the default
        if layout_id is None:
            for dataset in get_default_layout(cursor, default_domain_label):
                result['datasets'].append(dataset)
                dataset_ids.append(dataset.id)
        else:
            layout = geardb.Layout(id=layout_id)
            layout.load()

            dsc = geardb.DatasetCollection()
            dsc.get_by_dataset_ids(ids=layout.dataset_ids(), get_links=True)
            dsc.apply_layout(layout=layout)
            result['datasets'].extend(dsc.datasets)

    cursor.close()
    cnx.close()

    # apply any post-processing
    if only_types is not None:
        kept_datasets = []

        for dataset in result['datasets']:
            if dataset['dtype'] in only_types:
                kept_datasets.append(dataset)

        result['datasets'] = kept_datasets

    # does the user have a specific search requirement?
    if sort_order == 'alpha':
        datasets_sorted = sorted(result['datasets'], key=itemgetter('title'))
        result['datasets'] = datasets_sorted

    print(json.dumps(result))

def get_default_layout(cursor, domain_label):
    # this is the hearing one
    layout_id = 0

    # These values need to match what's in the database (check create_schema.sql)
    if domain_label == "Brain development (default)":
        layout_id = 10000
    elif domain_label == "Huntington's disease (default)":
        layout_id = 10001

    layout =  geardb.Layout(id=layout_id)
    layout.load()

    dsc = geardb.DatasetCollection()
    dsc.get_by_dataset_ids(ids=layout.dataset_ids(), get_links=True)
    dsc.apply_layout(layout=layout)

    return dsc.datasets

def get_users_datasets(cursor, user_id):
    qry = """
       SELECT d.id, d.title, o.label, d.pubmed_id, d.geo_id, d.is_public, d.ldesc,
       d.dtype, d.schematic_image, d.share_id, d.math_default,
       d.marked_for_removal, d.date_added, d.load_status, d.plot_default,
       IFNULL(GROUP_CONCAT(t.label), 'NULL') as tags, o.id
         FROM dataset d
              JOIN organism o ON d.organism_id=o.id
              LEFT JOIN dataset_tag dt ON dt.dataset_id = IFNULL(d.id, 'NULL')
              LEFT JOIN tag t ON t.id = IFNULL(dt.tag_id, 'NULL')
        WHERE d.owner_id = %s
        GROUP BY d.id, d.title, o.label, d.pubmed_id, d.geo_id, d.is_public, d.ldesc,
        d.dtype, d.schematic_image, d.share_id, d.math_default,
        d.marked_for_removal, d.date_added, d.load_status, d.plot_default, o.id
    """
    cursor.execute(qry, (user_id,))
    datasets = list()

    for row in cursor:
        # skip datasets marked_for_removal
        if row[11] == 1:
            continue
        else:
            if row[5] == 1:
                access_level = 'Public'
            else:
                access_level = 'Private'

            date_added = row[12].isoformat()

            if row[15] == 'NULL':
                tag_list = None
            else:
                tag_list = row[15].replace(',', ', ')

            datasets.append({
                'dataset_id': row[0],
                'grid_position': None,
                'grid_width': 4,
                'mg_grid_width': 6,
                'title': row[1],
                'organism': row[2],
                'organism_id': row[16],
                'pubmed_id': row[3],
                'geo_id': row[4],
                'access': access_level,
                'ldesc': row[6],
                'dtype': row[7],
                'user_id': user_id,
                'user_name': 'You',
                'schematic_image': row[8],
                'share_id': row[9],
                'math_format': row[10],
                'date_added': date_added,
                'load_status': row[13],
                'plot_format': row[14],
                'tags': tag_list
            })

    return datasets

def get_permalink_dataset(cursor, permalink_id):
    dataset_id = geardb.get_dataset_id_from_share_id(permalink_id)
    dsc = geardb.DatasetCollection()
    dsc.get_by_dataset_ids(ids=[dataset_id])
    datasets = list()

    if len(dsc.datasets):
        ds = dsc.datasets[0]
        ds.grid_position = 100
        ds.grid_width = 6
        ds.mg_grid_width = 6
        ds.is_permalink = 1
        datasets.append(ds)

    return datasets


def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
