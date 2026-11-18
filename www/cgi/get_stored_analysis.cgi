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
from gear.analysis import Analysis, get_analysis


def main():
    form = cgi.FieldStorage()
    analysis_id = form.getfirst('analysis_id')
    analysis_type = form.getfirst('analysis_type')
    session_id = form.getfirst('session_id')
    dataset_id = form.getfirst('dataset_id')

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

    # The analysis prints as its pipeline JSON. Read it before printing the header, so a missing
    #  pipeline file gives one JSON error rather than a second header inside the body
    try:
        analysis_data = json.loads(str(ana))
    except FileNotFoundError:
        analysis_data = {"success": 0, "error": "Analysis config file not found"}
    else:
        current_user = geardb.get_user_from_session_id(session_id) if session_id else None
        analysis_data["vetting"] = stored_vetting(ana, current_user.id if current_user else None)

    print('Content-Type: application/json\n\n')
    print(json.dumps(analysis_data, indent=4))

def stored_vetting(ana, current_user_id):
    """
    Vetting derived from the analysis's owner, the same way as get_stored_analysis_list.cgi. A vetting
    saved in the pipeline file is not used, since the browser writes that file.
    """
    if ana.type in ("user_saved", "user_unsaved"):
        # Loaded from the current user's own area, so they own it
        owner_ana = ana
    else:
        # The owner of a public analysis is found from the session that published it
        try:
            with open(ana.settings_path) as fh:
                owner_ana = Analysis.from_json(json.load(fh))
            owner_ana.type = ana.type
        except Exception as e:
            print(f"Could not determine the owner of analysis {ana.id}: {e}", file=sys.stderr)
            return None

    try:
        return owner_ana.discover_vetting(current_user_id)
    except Exception as e:
        print(f"Could not determine vetting for analysis {ana.id}: {e}", file=sys.stderr)
        return None

if __name__ == '__main__':
    main()
