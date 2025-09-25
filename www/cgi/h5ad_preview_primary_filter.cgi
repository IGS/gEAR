#!/opt/bin/python3

"""
Pre-filter graphics should have been generated already by this script:

    gear/bin/generate_initial_composition_plots.py

This CGI used to actually generate it live, but not just checks for the
presence of the images and returns success:1 if there, and success:0 if
not present.

"""

import cgi, json
import os, sys

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
from gear.analysis import Analysis

# this is needed so that we don't get TclError failures in the underlying modules
import matplotlib
matplotlib.use('Agg')

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')

    result = {'success': 1}

    ana = Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id, session_id=session_id)

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
