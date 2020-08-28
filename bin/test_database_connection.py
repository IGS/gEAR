#!/opt/bin/python3

"""
Verifies that a database connection can be made and cursor obtained using the gear.ini file
"""

import os
import sys

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    # perhaps add a basic SELECT query here?

    cursor.close()
    cnx.close()

if __name__ == '__main__':
    main()
