#!/opt/bin/python3

"""
This reads an image file created by a gEAR plotting function and prints it as a binary
stream, suitable for serving as an img['src'] from an html page.
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import get_analysis


def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    result = {"success": 0}

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return
    is_spatial = ds.dtype == "spatial"

    analysis_obj = None
    if analysis_id or analysis_type:
        analysis_obj = {
            'id': analysis_id if analysis_id else None,
            'type': analysis_type if analysis_type else None,
        }

    try:
        ana = get_analysis(analysis_obj, dataset_id, session_id, is_spatial=is_spatial)
    except Exception:
        print("Analysis for this dataset is unavailable.", file=sys.stderr)
        result['success'] = 0
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    ## if the analysis is primary images are getting saved as user_unsaved
    if ana.type == 'primary':
        ana.type = 'user_unsaved'

    # this has to be what the scanpy library calls it, since we don't seem to be able to
    #  control it: scanpy GitHub ticket #73
    analysis_name = form.getvalue('analysis_name')

    data_file_path = ana.dataset_path
    ana_directory = os.path.normpath(os.path.dirname(data_file_path))
    image_path = "{0}/figures/{1}.png".format(ana_directory, analysis_name)
    #print("DEBUG: streaming image at path: {0}".format(image_path), file=sys.stderr)

    # Normalize path to avoid directory traversal attacks (e.g. ../../../etc/passwd) and validate
    image_path = os.path.normpath(image_path)
    if not image_path.startswith(ana_directory):
        raise Exception("Invalid filename: {}".format(image_path))

    try:
        with open(image_path, 'rb') as f:
            print("Content-Type: image/png\n")
            sys.stdout.flush() # <---
            sys.stdout.buffer.write(f.read())
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        # ensure a 404 response
        print("Status: 404 Not Found\n")
        print("Content-Type: text/plain\n")
        print("File not found: {0}".format(image_path))
        print("Error: {0}".format(e))

if __name__ == '__main__':
    main()
