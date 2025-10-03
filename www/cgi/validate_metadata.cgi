#!/opt/bin/python3

import cgi
import json
import os, sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from gear.metadata import Metadata


def main():
    print('Content-Type: application/json\n\n', flush=True)

    user_upload_file_base = '../uploads/files'

    result = {'success':0, 'message':'', 'metadata':{}}

    form = cgi.FieldStorage()
    filename = form.getvalue('filename')
    filepath = user_upload_file_base + "/" + filename
    populate_from_geo_success = False

    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
    else:
        # gEAR User found. Read and Validate metadata file
        # Initize metadatauploader & read file
        metadata = Metadata(file_path=filepath)

        # Populates empty fields from GEO (if GEO GSE ID was given)
        try:
            metadata.populate_from_geo()
            populate_from_geo_success = True
        except KeyError:
            result['message'] = 'Unable to process GEO ID.  Please check it and try again.'

        # Validate inputs
        if populate_from_geo_success:
            metadata.validate()

            # Returning metadata and any error and warning messages to JS
            result['metadata'] = metadata.write_json()
            result['success'] = 1

    print(json.dumps(result))



if __name__ == '__main__':
    main()
