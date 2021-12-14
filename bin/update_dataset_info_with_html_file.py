#!/opt/bin/python3

"""
This loads an HTML file to be used as the long description field in the database,
allowing for longer, more complex info boxes.
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser( description='Changes ownership of datasets or profiles')
    parser.add_argument('-d', '--dataset_id', type=str, required=True, help='Dataset ID to update' )
    parser.add_argument('-i', '--input_file', type=str, required=True, help='HTML file whose contents to use' )
    args = parser.parse_args()

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    html_contents = open(args.input_file).read()

    qry = "UPDATE dataset SET ldesc = %s WHERE id = %s"
    cursor.execute(qry, (html_contents, args.dataset_id))

    cnx.commit()
    cursor.close()
    cnx.close()

if __name__ == '__main__':
    main()
    
