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

from gear.datasetuploader import FileType, DatasetUploader

# This is another attempt to fix the STDOUT hack above in a cleaner way, but I couldn't
#  get it to actually change the loglevel of another imported module.
#import logging
#logging.getLogger("matplotlib").setLevel(logging.WARNING)

def main():
    user_upload_file_base = '../datasets_epigenetic/uploads/files'

    result = {'success':0 }
    result['filename'] = filename
    
    form = cgi.FieldStorage()
    filename = form.getvalue('filename')

    # This means the PHP upload failed
    if filename is None:
        result['message'] = "File upload failed.  Try again and contact us if this continues."
        print_and_go(json.dumps(result))
    
    filepath = user_upload_file_base + "/" + filename
    session_id = form.getvalue('session_id')

    user = geardb.get_user_from_session_id(session_id)

    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
    else:
        # gEAR User found. Read and Validate expression file
        file_type = filename.rsplit('.', 1)[1]
        
        result['message'] = 'File uploaded - ' + str(filename)
        result['success'] = 1

    print_and_go(json.dumps(result))

def print_and_go(content):
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n', flush=True)
    print(content)
    sys.exit(0)
    
if __name__ == '__main__':
    main()
