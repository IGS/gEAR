#!/opt/bin/python3

"""
Add a nemoarchive submission as a new layout.

Script is a clone of add_layout.cgi with submission information in consideration.

Requires:
1) Session id - which contains user_id
2) Layout name to be added
"""

import json
import cgi
import os
import sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    layout_name = form.getvalue('layout_name')
    submission_id = form.getvalue("submission_id")

    if not submission_id:
        raise "Submissiom ID not found for this request."

    if not layout_name:
        layout_name = f"Submission {submission_id}"

    # Throw error if submission already has layout associated.
    # Submissions must have one-to-one relationships with layouts.
    submission = geardb.get_submission_by_id(submission_id)
    if not submission.layout_id:
        error =  f"Submission does not have an associated layout id. Aborting."
        result = {'error':error}
        print(json.dumps(result))
        return

    layout = geardb.Layout()
    layout.load(submission.layout_id)

    layout.save_change(attribute="label", value=layout_name)
    result = {'layout_id': layout.id,
            'layout_label': layout.label,
            'layout_share_id': layout.share_id
    }

    # Associate layout with submission.
    submission.save_change(attribute="layout_id", value=layout.id)

    print(json.dumps(result))

if __name__ == '__main__':
    main()
