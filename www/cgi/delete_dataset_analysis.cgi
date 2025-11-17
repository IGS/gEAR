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
from gear.analysis import Analysis


def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    result = {"success": 0, "error": ""}

    user = geardb.get_user_from_session_id(session_id)

    if not user:
        result["success"] = 0
        result["error"] = "User not found"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    # Safeguard check
    if analysis_type and analysis_type in ["primary"]:
        result["success"] = 0
        result["error"] = "Cannot delete primary analyses"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    try:
        source_ana = Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id)
    except Exception:
        print("Analysis for this dataset is unavailable.", file=sys.stderr)
        result['success'] = 0
        result['error'] = "Analysis could not be found."
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    source_pipeline_base = source_ana.base_path

    # verify this analysis belongs to the user
    if not os.path.exists(source_pipeline_base):
        result["success"] = 0
        result["error"] = "Analysis not found"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    # ? Is this necessary if we already verified the session_id above?
    if not source_ana.user_id == user.id:
        result["success"] = 0
        result["error"] = "Analysis does not belong to user"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    try:
        # move directory to a backup in case of failures
        #shutil.move(source_pipeline_base, source_pipeline_base + '.deleted')

        # find any dataset displays with this analysis id in the plotly_config
        """
        conn = Connection()
        cursor = conn.get_cursor()
        qry = "SELECT * from dataset_display where plotly_config like '%\"id\":\"{}\"%'".format(analysis_id)
        cursor.execute(qry, (self.id,))
        for (display_id, user_id, label, plot_type, plotly_config) in cursor:
            display = DatasetDisplay(
                id=display_id,
                dataset_id=self.id,
                user_id=user_id,
                label=label,
                plot_type=plot_type,
                plotly_config=plotly_config
            )
            display.remove()
        """

        #shutil.rmtree(source_pipeline_base + '.deleted')
        # Only allow deletion if the final path component is a UUID
        base_name = os.path.basename(os.path.normpath(source_pipeline_base))
        if re.match(r'^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$', base_name):
            shutil.rmtree(source_pipeline_base)
            result = {'success': 1}
        else:
            result = {'success': 0, 'error': 'Refusing to delete: invalid target pattern'}

    except Exception:
        # restore the backup
        #shutil.move(source_pipeline_base + '.deleted', source_pipeline_base)
        result["success"] = 0
        result["error"] = "Unable to delete dataset display"

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

