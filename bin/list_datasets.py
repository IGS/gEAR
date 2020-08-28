#!/opt/bin/python3

"""
This is used to list the datasets in gEAR.  It excludes any marked for removal.

List is sorted by most recently-uploaded first.
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join('..', 'lib'))
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
                 JOIN guser g on d.owner_id=g.id
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
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(dataset_id, title, dtype, email, pubmed_id, geo_id, visibility, date_added))

    cursor.close()
    cnx.close()

if __name__ == '__main__':
    main()
    
