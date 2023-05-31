#!/opt/bin/python3

"""
nemoarchive_enable_submission_emails.cgi - Enable email notifications for user of submission
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

    if not submission_id:
        raise "Submissiom ID not found for this request."

    result = {"success": False, "message":""}

    try:
        submission = geardb.get_submission_by_id(submission_id)
        submission.save_change(attribute="email_updates", value=1)
        result["success"] = True
    except Exception as e:
        result["message"] = str(e)

    print(json.dumps(result))

if __name__ == '__main__':
    main()
