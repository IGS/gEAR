#!/opt/bin/python3

# find_existed_submitted_datasets.cgi - Check database if current datasets were loaded in previous submissions

import cgi
import json
import os,sys

def main():
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')

    result = {}
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()