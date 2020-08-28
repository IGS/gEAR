#!/opt/bin/python3

"""

"""

import cgi, json
import os, sys, re
import shutil

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    source_ana = geardb.Analysis(id=analysis_id, type=analysis_type,
                                 dataset_id=dataset_id, session_id=session_id, user_id=user.id)
    source_pipeline_base = source_ana.base_path()

    result = {'success': 0}

    try:
        # TODO: need to verify ownership here in the future before deleting
        #print("DEBUG: would rmtree this:{0}".format(source_pipeline_base), file=sys.stderr)
        shutil.rmtree(source_pipeline_base)
        result = {'success': 1}
    except:
        result = {'success': 0, 'error': 'Unable to delete dataset'}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

