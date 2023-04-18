#!/opt/bin/python3

# initialize_new_submission.cgi -
#
# 1) Check database if current datasets were loaded in previous submissions
# 2) Insert new submission and datasets

import cgi
import json
import os,sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def get_or_save_submission_dataset(params):
    # Check if dataset already exists in a previous submission
    s_dataset = geardb.get_submission_dataset_by_dataset_id(params["dataset_id"])
    if s_dataset:
        return s_dataset
    # Not found, so insert new submission dataset
    # NOTE: We do not insert into "dataset" until metadata is validated and dataset converted to H5AD
    s_dataset = geardb.SubmissionDataset(dataset_id=params["dataset_id"],
                        is_restricted=params["is_restricted"], nemo_identifier=params["identifier"],
                        pulled_to_vm_status=0, convert_metadata_status=0, convert_to_h5ad_status=0
                        )
    s_dataset.save()
    return s_dataset

def save_submission(submission_id, user_id, is_restricted):
    """
    Adds a new Submission to the database
    """
    conn = geardb.Connection()
    cursor = conn.get_cursor()

    submission = None
    qry = """
        INSERT INTO submission (id, user_id, is_finished, is_restricted, date_added)
        VALUES (%s, %s, 0, %s, NOW())
    """

    cursor.execute(qry, (submission_id, user_id, is_restricted))
    cursor.close()
    conn.commit()
    conn.close()

def main():
    form = cgi.FieldStorage(environ="post")
    file_info = form.getvalue('file_entities')
    submission_id = form.getvalue('submission_id')
    session_id = form.getvalue('session_id')

    result = {"success":0, "message":"", "dataset_status":{}}
    print('Content-Type: application/json\n\n')

    # Must have a gEAR account to upload datasets
    user = geardb.get_user_from_session_id(session_id)

    # TODO: Set to be obtained from the dataset metadata but currently not in metadata
    is_restricted = 0

    save_submission(submission_id, user.id, is_restricted)

    for file_entity in file_info:
        attributes = file_entity["attributes"]

        submission_dataset_params = {
            "is_restricted": is_restricted
            , "dataset_id": attributes["id"]
            , "identifier": attributes["identifier"]
            }

        s_dataset = get_or_save_submission_dataset(submission_dataset_params)

        result["dataset_status"][s_dataset.dataset_id] = {
            "pulled_to_vm": s_dataset.pulled_to_vm_status
            , "convert_metadata": s_dataset.convert_metadata_status
            , "convert_to_h5ad": s_dataset.convert_to_h5ad_status
            }

        # Insert submission members and create new submisison datasets if they do not exist
        submission_member = geardb.SubmissionMember(submission_id=submission_id, submission_dataset_id=s_dataset.id)
        submission_member.save()


    print(json.dumps(result))

if __name__ == '__main__':
    main()