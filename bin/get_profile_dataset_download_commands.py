#!/opt/bin/python3

"""
Doesn't actually DO anything aside from printing a list of gcloud commands to download
the datasets contained within a profile for local testing.

Assumes you have gcloud utilities set up locally
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser( description="Generates commands to download a profile's datasets")
    parser.add_argument('-p', '--profile_share_id', type=str, required=True, help='Pass a profile/layout share id' )
    args = parser.parse_args()

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    qry = """SELECT d.id
               FROM layout_members lm
                    JOIN dataset d on lm.dataset_id=d.id
                    JOIN layout l ON lm.layout_id=l.id
              WHERE l.share_id = %s
    """

    cursor.execute(qry, (args.profile_share_id,))

    for (dataset_id) in cursor:
        print("gcloud compute scp gear-prod-20200713:/home/jorvis/git/gEAR/www/datasets/{0}.h5ad .".format(dataset_id[0]))

    cursor.close()
    cnx.close()

if __name__ == '__main__':
    main()
    
