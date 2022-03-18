#!/opt/bin/python3

"""

"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# this obviously needs to change
PATTERN_BASE_DIR = '/var/www/patterns/HuttCtxDevoLMDhs'

def main():
    form = cgi.FieldStorage()
    file_name = form.getvalue('file_name')
    result = []
    file_path = os.path.join(PATTERN_BASE_DIR, file_name)

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
