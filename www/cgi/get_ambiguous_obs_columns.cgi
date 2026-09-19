#!/opt/bin/python3

"""
Used by the expression uploader's "Review column types" step.

Returns the obs columns flagged as possibly-mislabeled numeric/categorical,
as already written into metadata.json by the anndata/spatial processor when
the dataset finished converting. No H5AD/Zarr file is opened here at all --
that scan already happened once during processing, so this is just a fast
metadata.json read to populate the review dropdowns.

Reads: ../uploads/files/<session_id>/<share_uid>/metadata.json
"""

import cgi
import json
import os
import sys
from pathlib import Path

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

user_upload_file_base = '../uploads/files'


def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    share_uid = form.getfirst('share_uid')

    result = {'success': 0, 'message': '', 'questionable_columns': {}, 'reviewed': False}

    if not session_id or not share_uid:
        result['message'] = 'session_id and share_uid are required.'
        print(json.dumps(result))
        return

    user = geardb.get_user_from_session_id(session_id)
    if not user:
        result['message'] = 'Invalid session_id.'
        print(json.dumps(result))
        return

    metadata_file = Path(user_upload_file_base) / session_id / share_uid / 'metadata.json'
    if not metadata_file.is_file():
        result['message'] = 'No metadata JSON file found for this dataset.'
        print(json.dumps(result))
        return

    try:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        result['message'] = f'Error reading metadata: {e}'
        print(json.dumps(result))
        return

    result['success'] = 1
    result['questionable_columns'] = metadata.get('questionable_obs_columns', {})
    result['reviewed'] = metadata.get('obs_dtype_reviewed', False)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
