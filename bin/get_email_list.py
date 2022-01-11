#!/usr/bin/env python3

"""
This is used to get a comma-separated list of all the e-mail addresses where 
users asked to be notified of updates.
"""

import gzip
import json
import os
import re
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import gear.db

def main():

    #connect to gEAR MySQL
    mysql_cnx = gear.db.MySQLDB().connect()
    cursor = mysql_cnx.cursor()
    
    qry = "SELECT email FROM guser WHERE updates_wanted = 1"
    cursor.execute(qry)

    email_count = 0

    for row in cursor:
        email_count += 1
        print(row[0], end=', ')

    sys.stdout.flush()
    print("\n\nExported {0} e-mail addresses".format(email_count), file=sys.stderr)


if __name__ == '__main__':
    main()
