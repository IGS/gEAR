#!/opt/bin/python3

"""
Pre-filter graphics should have been generated already by this script:

    gear/bin/generate_initial_composition_plots.py

This CGI used to actually generate it live, but not just checks for the
presence of the images and returns success:1 if there, and success:0 if
not present.

"""

import cgi
import json
import os
import sys

import matplotlib

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import get_analysis

# this is needed so that we don't get TclError failures in the underlying modules
matplotlib.use('Agg')

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')

    result = {'success': 1}

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return
    is_spatial = ds.dtype == "spatial"

    result["is_spatial"] = 1 if is_spatial else 0

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
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    source_datafile_path = ana.dataset_path
    violin_path = str(source_datafile_path).replace('.h5ad', '.prelim_violin.png')

    if not os.path.exists(violin_path):
        result['success'] = 0
    else:
        result['success'] = 1

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
