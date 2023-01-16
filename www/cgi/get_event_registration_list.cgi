#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print("{}")

if __name__ == '__main__':
    main()
