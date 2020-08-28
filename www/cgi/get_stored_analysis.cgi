#!/opt/bin/python3

"""
Used by analyze_dataset.html, this script gets the JSON for a single stored analysis 
"""

import cgi
import json
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

    ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, 
                          user_id=user.id, session_id=session_id,
                          type=analysis_type)

    print('Content-Type: application/json\n\n')
    print(ana)

if __name__ == '__main__':
    main()
