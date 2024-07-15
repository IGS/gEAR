#!/opt/bin/python3

import cgi
import json
import os, sys

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.metadata import Metadata

def main():
    user_upload_file_base = '/tmp'
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('metadata-dataset-id')
    fileitem = form['metadata-file-input']
    user = geardb.get_user_from_session_id(session_id)
    result = {'success':0 }

    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print_and_go(json.dumps(result))
    
    #filename = os.path.basename(form.getvalue('metadata-file-input'))
    dest_filepath = os.path.join(user_upload_file_base, "{0}.xlsx".format(dataset_id))

    fh = open(dest_filepath, 'wb')
    fh.write(fileitem.file.read())
    fh.close()

    metadata = Metadata(file_path=dest_filepath)
    result['metadata'] = json.loads(metadata.write_json())
    result['success'] = 1
    print_and_go(json.dumps(result))

def print_and_go(content):
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n', flush=True)
    print(content)
    sys.exit(0)
    
if __name__ == '__main__':
    main()
