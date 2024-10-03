#!/opt/bin/python3

"""
Used by sc_workbench.html, this script gets the JSON for a single stored analysis
"""

import cgi
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')

    user = geardb.get_user_from_session_id(session_id)

    if not user:
        print('Content-Type: application/json\n\n')
        print('{"error": "Invalid session_id"}')
        return

    ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id,
                          user_id=user.id, session_id=session_id,
                          type=analysis_type)

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
