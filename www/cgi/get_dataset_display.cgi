#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    display_id = form.getvalue('display_id')

    display = geardb.get_display_by_id(display_id=display_id)

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(display))

if __name__ == '__main__':
    main()
