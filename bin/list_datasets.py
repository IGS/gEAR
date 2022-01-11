#!/opt/bin/python3

"""
This is used to list the datasets in gEAR.  It excludes any marked for removal.

List is sorted by most recently-uploaded first.
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser( description='Lists datasets in the gEAR')
    parser.add_argument('-uid', '--user_id', type=int, required=False, help='Filter by user ID' )
    args = parser.parse_args()

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    qry = """
          SELECT d.id, d.title, d.dtype, g.email, d.pubmed_id, d.geo_id, d.is_public, d.date_added
            FROM dataset d
                 JOIN guser g ON d.owner_id=g.id
           WHERE d.marked_for_removal = 0
    """
    qry_args = []

    if args.user_id:
        qry += " AND d.owner_id = %s"
        qry_args.append(args.user_id)

    qry += " ORDER BY date_added DESC"

    cursor.execute(qry, qry_args)

    for (dataset_id, title, dtype, email, pubmed_id, geo_id, is_public, date_added) in cursor:
        visibility = 'Public' if is_public else 'Private'
        layout_str = get_layout_string(cnx, dataset_id)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(layout_str, dataset_id, title, dtype, email, pubmed_id, geo_id, visibility, date_added))

    cursor.close()
    cnx.close()

def get_layout_string(cnx, dataset_id):
    cursor = cnx.get_cursor()
    labels = list()

    qry = """
          SELECT l.share_id, l.label
            FROM layout l 
                 JOIN layout_members lm ON lm.layout_id=l.id
           WHERE lm.dataset_id = %s
    """
    cursor.execute(qry, [dataset_id,])

    for (share_id, label) in cursor:
        labels.append(label)

    cursor.close()

    if len(labels):
        return " / ".join(labels)
    else:
        return "None"

    
    
if __name__ == '__main__':
    main()
    
