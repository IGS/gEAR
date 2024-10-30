#!/opt/bin/python3

"""
Replace a display in a nemoarchive submission layout
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
    submission_id = form.getvalue("submission_id")
    old_display_id = form.getvalue("old_display_id")
    new_display_id = form.getvalue("new_display_id")

    if not submission_id:
        raise Exception("Submissiom ID not found for this request.")

    if not old_display_id or not new_display_id:
        raise Exception("Display IDs need to be provided in order to update the layout.")

    # Throw error if submission already has layout associated.
    # Submissions must have one-to-one relationships with layouts.
    submission = geardb.get_submission_by_id(submission_id)
    if not submission.layout_id:
        error =  f"Submission does not have an associated layout id. Aborting."
        result = {'error':error}
        print(json.dumps(result))
        return

    layout = submission.get_layout_info()
    if not layout:
        error = f"Layout not found for submission {submission_id}."
        result = {'error':error}
        print(json.dumps(result))
        return

    layout.load()

    members = layout.members
    display_ids = [m.display_id for m in members]
    if old_display_id not in display_ids:
        error = f"Display ID {old_display_id} not found in layout members."
        result = {'error':error}
        print(json.dumps(result))
        return

    # Find the member with the old display id and replace it with the new display id.
    for member in members:
        if member.display_id == old_display_id:
            member.save_change(attribute="display_id", value=new_display_id)
            break

    result = {'layout_id': layout.id,
            'layout_label': layout.label,
            'layout_share_id': layout.share_id
    }

    print(json.dumps(result))

if __name__ == '__main__':
    main()
