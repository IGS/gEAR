#!/opt/bin/python3

"""
This is used by the expression uploader and accepts an metadata file upload,
stores the file in temp space, parses it to JSON, and returns that. Not meant
for permanent storage.
"""

import cgi
import json
import os
import sys

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.metadata import Metadata


def main() -> dict:
    user_upload_file_base = '/tmp'
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('metadata-dataset-id')
    fileitem = form['metadata-file-input']
    user = geardb.get_user_from_session_id(session_id)
    result = {'success':0, 'message':'', 'metadata':{}}

    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        return result

    filename = os.path.basename(fileitem.filename)

    # does the file end in .xls or .xlsx?
    if not (filename.endswith('.xls') or filename.endswith('.xlsx')):
        result['message'] = 'Metadata file must be in Excel format (with .xls or .xlsx extension)'
        return result

    #filename = os.path.basename(form.getvalue('metadata-file-input'))
    dest_filepath = os.path.join(user_upload_file_base, "{0}.xlsx".format(dataset_id))

    fh = open(dest_filepath, 'wb')
    fh.write(fileitem.file.read())
    fh.close()

    metadata = Metadata(file_path=dest_filepath)

    try:
        metadata.populate_from_geo()
    except Exception:
        pass

    try:
        metadata_json = metadata.write_json()
        result['metadata'] = json.loads(metadata_json)
    except Exception as e:
        result['message'] = 'Error processing metadata file: {0}'.format(str(e))
        return result

    # TODO: run geo if that property is defined

    result['success'] = 1
    return result

if __name__ == '__main__':
    result = main()
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))
    sys.exit(0)
