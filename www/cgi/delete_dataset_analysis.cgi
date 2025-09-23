#!/opt/bin/python3

"""

"""

import cgi
import json
import os
import shutil
import sys

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
    user = geardb.get_user_from_session_id(session_id)

    if not user:
        print('Content-Type: application/json\n\n')
        print(json.dumps({'success': 0, 'error': 'User not found'}))
        return

    source_ana = Analysis(id=analysis_id, type=analysis_type,
                                 dataset_id=dataset_id, session_id=session_id, user_id=user.id)
    source_pipeline_base = source_ana.base_path()

    result = {'success': 0}

    # verify this analysis belongs to the user
    if not os.path.exists(source_pipeline_base):
        result = {'success': 0, 'error': 'Analysis not found'}
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return
    if not source_ana.user_id == user.id:
        result = {'success': 0, 'error': 'Analysis does not belong to user'}
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
        shutil.rmtree(source_pipeline_base)
        result = {'success': 1}
    except:
        # restore the backup
        #shutil.move(source_pipeline_base + '.deleted', source_pipeline_base)
        result = {'success': 0, 'error': 'Unable to delete dataset'}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

