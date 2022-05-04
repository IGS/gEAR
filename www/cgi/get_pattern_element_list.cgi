#!/opt/bin/python3

"""

"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

PATTERN_BASE_DIR = os.path.abspath(os.path.join('..', 'patterns'))
WEIGHTED_CARTS_BASE_DIR = os.path.abspath(os.path.join('..', 'carts'))

VALID_SCOPES = ["projection_pattern", "genecart"]

def main():
    form = cgi.FieldStorage()
    file_name = form.getvalue('file_name')
    scope = form.getvalue('scope', "projection_pattern")

    if scope not in VALID_SCOPES:
        print("Status: 501 Not Implemented")
        print('Content-Type: application/json\n\n')
        print("Invalid scope: {}".format(scope))
        return

    result = []
    file_path = os.path.join(PATTERN_BASE_DIR, file_name)
    if scope == "genecart":
        # The gene cart "share_id" is passed to the CGI script as the file_name, which we can get the actual cart name from
        file_name = "cart.{}.tab".format(file_name)
        file_path = os.path.join(WEIGHTED_CARTS_BASE_DIR, file_name)

    for line in open(file_path):
        line = line.rstrip()
        # on the first line, the file should be all pattern names after the first column
        cols = line.split("\t")

        for col in cols[1:]:
            result.append({'label': col})

        # we only care about the first line
        break

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
