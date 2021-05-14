#!/opt/Python-3.7.3/bin/python3

"""
We want to test and make sure queries from MySQL are returning
strings and not bytarrays.

Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import os
import sys
import gzip
import mysql.connector

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser( description='Queries the database to print strings so users can check strings vs byte-arrays')

    ## output file to be written
    args = parser.parse_args()

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    qry = "SELECT gene_symbol FROM gene LIMIT 10"
    cursor.execute(qry)

    for (gs,) in cursor:
        print(gs)
    
    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()







