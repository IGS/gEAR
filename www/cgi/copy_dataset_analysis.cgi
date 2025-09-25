#!/opt/bin/python3

"""
This copies an entire analysis, with images, into another place in the analysis
directory structure.

Use cases:

- Migrating class 'user_unsaved' to 'user_saved'
    Actions:
       - Move entire directory
       - Maintain open permissions
- Migrating class 'user_saved' to 'public'
       - Move entire directory
       - Close permissions so it's not accidentally modified
- Copying class 'public' to 'user_unsaved'
       - Copy from public directory
       - Make sure permissions are writeable
"""

import cgi
import json
import os
import re
import shutil
import sys

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import Analysis
from werkzeug.utils import secure_filename


def main():
    form = cgi.FieldStorage()
    source_analysis_id = form.getvalue('source_analysis_id')
    dest_analysis_id = form.getvalue('dest_analysis_id')
    source_analysis_type = form.getvalue('source_analysis_type')
    dest_analysis_type = form.getvalue('dest_analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    if not user:
        result = {'success': 0, 'error': 'Invalid session'}
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    source_analysis_id = secure_filename(source_analysis_id)
    dest_analysis_id = secure_filename(dest_analysis_id)

    source_ana = Analysis(id=source_analysis_id, type=source_analysis_type, dataset_id=dataset_id, session_id=session_id, user_id=user.id)
    dest_ana = Analysis(id=dest_analysis_id, type=dest_analysis_type, dataset_id=dataset_id, session_id=session_id, user_id=user.id)

    if not source_ana:
        result = {'success': 0, 'error': 'Source analysis does not exist'}
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    if not dest_ana:
        result = {'success': 0, 'error': 'Destination analysis does not exist'}
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    source_pipeline_base = source_ana.base_path
    dest_pipeline_base = dest_ana.base_path

    if source_analysis_type == 'user_unsaved' and dest_analysis_type == 'user_saved':
        shutil.move(source_pipeline_base, dest_pipeline_base, copy_function=open_perm_changing_copy)
        set_config_analysis_type(dest_ana.settings_path, dest_analysis_type, session_id, dest_analysis_id)
    elif source_analysis_type == 'user_saved' and dest_analysis_type == 'public':
        set_config_analysis_type(source_ana.settings_path, dest_analysis_type, session_id, dest_analysis_id)
        shutil.move(source_pipeline_base, dest_pipeline_base, copy_function=closed_perm_changing_copy)
    elif source_analysis_type == 'public' and dest_analysis_type == 'user_unsaved':
        dest_pipeline_base_parent = os.path.abspath(os.path.join(dest_pipeline_base, '../'))
        os.makedirs(dest_pipeline_base_parent, exist_ok=True, mode=0o777)
        shutil.copytree(source_pipeline_base, dest_pipeline_base, copy_function=open_perm_changing_copy)
        set_config_analysis_type(dest_ana.settings_path, dest_analysis_type, session_id, dest_analysis_id)
    else:
        result = {'success': 0, 'error': 'Unrecognized source and dest analysis types'}

    result = {'success': 1}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

def closed_perm_changing_copy(src, dst):
    """
    The standard copy methods in shutil handle permissions strangely.  There doesn't seem to be a native
    way to copy a read-only source to a destination where umask is respected.
    """
    shutil.copy(src, dst)
    return os.chmod(dst, 0o444)

def open_perm_changing_copy(src, dst):
    """
    The standard copy methods in shutil handle permissions strangely.  There doesn't seem to be a native
    way to copy a read-only source to a destination where umask is respected.
    """
    shutil.copy(src, dst)
    return os.chmod(dst, 0o755)

def set_config_analysis_type(config_path, atype, session_id, analysis_id):
    """
    Opens an existing config file, changes the analysis type, then rewrites the file.
    """
    with open(config_path) as json_file:
        state = json.load(json_file)

    state['id'] = analysis_id
    state['type'] = atype
    state['user_session_id'] = session_id

    with open(config_path, 'w') as json_out_file:
        json.dump(state, json_out_file)

if __name__ == '__main__':
    main()

