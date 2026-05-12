#!/opt/bin/python3

"""
Used by sc_workbench.html, this script gets the JSON for a single stored analysis
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
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')

    result = {"success": 0, "error": ""}

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        result['error'] = "Dataset not found"
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
    except Exception as e:
        print(f"Error occurred: {e}", file=sys.stderr)
        result['success'] = 0
        result['error'] = "Analysis could not be loaded"
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    # try to read the analysis object and raise exception if FileNotFoundError

    try:
        print('Content-Type: application/json\n\n')
        print(ana)
    except FileNotFoundError:
        print('Content-Type: application/json\n\n')
        print('{"error": "Analysis config file not found"}')
        return

if __name__ == '__main__':
    main()
