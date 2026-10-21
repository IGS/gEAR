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
from gear.analysis import Analysis, SpatialAnalysis

UUID_PATTERN = re.compile(r'[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}')


def main() -> dict:
    form = cgi.FieldStorage()
    analysis_id = form.getfirst('analysis_id')
    analysis_type = form.getfirst('analysis_type')
    analysis_vetting = form.getfirst('analysis_vetting')
    dataset_id = form.getfirst('dataset_id')
    session_id = form.getfirst('session_id')
    state = form.getfirst('state')
    label = form.getfirst('label')
    user = geardb.get_user_from_session_id(session_id)

    result = {'success': 0, 'error': ""}

    if user is None:
        result['success'] = 0
        result['error'] = "User not found"
        return result

    if analysis_vetting in ["undefined", "null", ""]:
        analysis_vetting = None

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        return result

    is_spatial = ds.dtype == "spatial"

    # Only the user's own analyses can be saved (the workbench saves nothing else), and only by a
    #  plain UUID: a 'primary' or 'public' type, or an ID like "../..", would write the pipeline file
    #  over a dataset's shared analysis or into another directory, since base_path joins the ID as given
    if analysis_type not in ["user_saved", "user_unsaved"] or not UUID_PATTERN.fullmatch(analysis_id or ""):
        result['error'] = "Only your own saved or unsaved analyses can be saved"
        return result

    if is_spatial:
        ana_cls = SpatialAnalysis
    else:
        ana_cls = Analysis
    ana = ana_cls(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id, label=label,
                          vetting=analysis_vetting)

    pipeline_path = ana.settings_path
    try:
        ofh = open(pipeline_path, 'wt')
        ofh.write(state)
        result['success'] = 1
    except Exception:
        import traceback
        print(traceback.format_exc(), file=sys.stderr)
        result['success'] = 0
        result['error'] = "An error occurred when attempting to write analysis pipeline path file"

    return result


if __name__ == '__main__':
    result = main()
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

