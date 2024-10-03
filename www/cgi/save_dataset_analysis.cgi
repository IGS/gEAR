#!/opt/bin/python3

"""
Saves the JSON representation of an analysis pipeline
"""

import cgi, json
import os, sys, re

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    analysis_vetting = form.getvalue('analysis_vetting')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    state = form.getvalue('state')
    label = form.getvalue('label')
    user = geardb.get_user_from_session_id(session_id)

    if user is None:
        result = {'success': 0, 'error': "User not found"}
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    if analysis_vetting in ["undefined", "null", ""]:
        analysis_vetting = None

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id, label=label,
                          vetting=analysis_vetting)

    pipeline_path = ana.settings_path()
    try:
        ofh = open(pipeline_path, 'wt')
        ofh.write(state)
        result = {'success': 1}
    except:
        result = {'success': 0, 'error': "An error occurred when attempting to write analysis pipeline path file"}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

