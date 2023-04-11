#!/usr/bin/env python3

"""
This lets analysis/devs look up a user ID without getting into the database directly.

Whatever string is passed is used to search the user_name and email fields

Example uses:

./get_user_id joshua
./get_user_id jorvis@gmail.com

"""

import os
import sys
import argparse

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

import geardb

def main():
    parser = argparse.ArgumentParser(description='Get user info from the DB')
    parser.add_argument('search_str', type=str, help='Search string')
    args = parser.parse_args()

    #connect to gEAR MySQL
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    
    qry = "SELECT id, user_name, email FROM guser WHERE user_name LIKE %s OR email LIKE %s"

    #cursor.execute(qry, ('%' + str(term) + '%', ))
    cursor.execute(qry, ('%' + str(args.search_str) + '%', '%' + str(args.search_str) + '%'))

    rows_found = 0
    
    for row in cursor:
        rows_found += 1

        if rows_found == 1:
            print("ID\tUser name\tE-mail")
        
        print("{0}\t{1}\t{2}".format(row[0], row[1], row[2]))

    if not rows_found:
        print("No matches found")

    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()
