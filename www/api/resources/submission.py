#!/opt/bin/python3

from flask import request
from flask_restful import Resource, abort

import cgi
import json
import os,sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

class Submission(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id):
        """Retrieve a particular submission based on ID."""
        result = {"self": request.path, "success":0, "datasets":{}, "message":""}

        session_id = request.cookies.get('session_id')
        # Must have a gEAR account to upload datasets
        user = geardb.get_user_from_session_id(session_id)

        submission = geardb.get_submission_by_id(submission_id)
        if not submission:
            return result

        sd = geardb.SubmissionDatasetCollection()
        submission.datasets = sd.get_by_submission_id(submission_id)

        for s_dataset in submission.datasets:
            dataset_id = s_dataset.dataset_id
            result["datasets"][dataset_id] = {"message":""}
            result["datasets"][dataset_id]["dataset_status"] = {
                "pulled_to_vm": s_dataset.pulled_to_vm_status
                , "convert_metadata": s_dataset.convert_metadata_status
                , "convert_to_h5ad": s_dataset.convert_to_h5ad_status
                }
        result["success"] = 1
        return result

class Submissions(Resource):
    """Requests to deal with multiple submissions, including creating new ones"""

    def post(self):
        """Create a new submission in the database."""
        # Must have a gEAR account to create submissions
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)

        req = request.get_json()
        file_info = req.get("file_entities")
        submission_id = req.get("submission_id")

        result = {"self": request.path, "success":0, "datasets":{}, "message":""}

        issues_found = False

        # TODO: Set to be obtained from the dataset metadata but currently not in metadata
        is_restricted = 0

        submission = geardb.get_submission_by_id(submission_id)
        if submission:
            abort(400, message="Submission already exists for this ID.")

        try:
            save_submission(submission_id, user.id, is_restricted)
        except:
            abort(400, message="Could not save new submission to database")

        # TODO: Think about making separate API call for submission datasets and members
        for file_entity in file_info:
            attributes = file_entity["attributes"]

            dataset_id = attributes["id"]
            result["datasets"][dataset_id] = {"message":""}

            if not attributes["identifier"]:
                result["datasets"][dataset_id]["message"] = "NeMO Identifier empty. Cannot load."
                result["datasets"][dataset_id]["dataset_status"] = {
                    "pulled_to_vm": "canceled"
                    , "convert_metadata": "canceled"
                    , "convert_to_h5ad": "canceled"
                    }
                issues_found = True
                continue

            submission_dataset_params = {
                "is_restricted": is_restricted
                , "dataset_id": dataset_id
                , "identifier": attributes["identifier"]
                }

            # Dataset may already be inserted from a previous submission
            try:
                s_dataset = add_submission_dataset(submission_dataset_params)
            except Exception as e:
                print(str(e), file=sys.stderr)
                result["datasets"][dataset_id]["message"] = f"Could not save dataset {attributes['id']} as a new SubmissionDataset in database"
                result["datasets"][dataset_id]["dataset_status"] = {
                    "pulled_to_vm": "canceled"
                    , "convert_metadata": "canceled"
                    , "convert_to_h5ad": "canceled"
                    }
                issues_found = True
                continue

            result["datasets"][dataset_id]["dataset_status"] = {
                "pulled_to_vm": s_dataset.pulled_to_vm_status
                , "convert_metadata": s_dataset.convert_metadata_status
                , "convert_to_h5ad": s_dataset.convert_to_h5ad_status
                }

            # Insert submission members and create new submisison datasets if they do not exist
            try:
                save_submission_member(submission_id, s_dataset)
            except Exception as e:
                print(str(e), file=sys.stderr)
                result["datasets"][dataset_id]["message"] = f"Could not save dataset {attributes['id']} as a new SubmissionMember in database"
                result["datasets"][dataset_id]["dataset_status"] = {
                    "pulled_to_vm": "canceled"
                    , "convert_metadata": "canceled"
                    , "convert_to_h5ad": "canceled"
                    }
                issues_found = True
                continue

        result["success"] = 1
        result["message"] = "Submission database entries have been populated"
        if issues_found:
            result["message"] += " but there were issues found."

        return result



def add_submission_dataset(params):
    """Add existing or new submission dataset. If existing, reset all non-complete steps."""

    # Dataset is a foreign key in SubmissionDataset so we need to ensure we do not duplicate
    s_dataset = geardb.get_submission_dataset_by_dataset_id(params["dataset_id"])
    if s_dataset:
        # In order to ensure the UI looks correct,reset all non-complete back to pending
        # since they will be attempted again
        s_dataset.reset_incomplete_steps()
        return s_dataset

    # Not found, so insert new submission dataset
    # NOTE: We do not insert into "dataset" until metadata is validated and dataset converted to H5AD
    try:
        return save_submission_dataset(params["dataset_id"], params["identifier"], params["is_restricted"])
    except Exception as e:
        raise

def insert_minimal_dataset(dataset_id, identifier):
    """
    Insert a minimally viable dataset. We will populate it with metadata later.

    This is purely to get the SubmissionDataset row inserted.
    """

    add_dataset_sql = """
    INSERT INTO dataset (id, share_id, owner_id, title, organism_id, date_added, load_status)
    VALUES              (%s, %s,       %s,       "PRESTAGE NEMO IMPORT - %s",    1, NOW(), "pending")
    """

    (before, sep, share_id) = identifier.rpartition("nemo:")

    # Insert dataset info to database
    conn = geardb.Connection()
    cursor = conn.get_cursor()
    cursor.execute(add_dataset_sql, (dataset_id, share_id, geardb.find_importer_id().id, share_id))
    cursor.close()
    conn.commit()
    conn.close()

def save_submission(submission_id, user_id, is_restricted):
    """
    Adds a new Submission to the database
    """
    conn = geardb.Connection()
    cursor = conn.get_cursor()

    qry = """
        INSERT INTO submission (id, user_id, is_finished, is_restricted, date_added)
        VALUES (%s, %s, 0, %s, NOW())
    """

    cursor.execute(qry, (submission_id, user_id, is_restricted))
    cursor.close()
    conn.commit()
    conn.close()

def save_submission_dataset(dataset_id, identifier, is_restricted):
    # Dataset is a foreign key in SubmissionDataset so we need to ensure we do not duplicate
    dataset = geardb.get_dataset_by_id(dataset_id)
    if not dataset:
        insert_minimal_dataset(dataset_id, identifier)

    s_dataset = geardb.SubmissionDataset(dataset_id=dataset_id,
                    is_restricted=is_restricted, nemo_identifier=identifier,
                    pulled_to_vm_status="pending", convert_metadata_status="pending", convert_to_h5ad_status="pending"
                    )
    s_dataset.save()
    return s_dataset

def save_submission_member(submission_id, s_dataset):
    submission_member = geardb.SubmissionMember(submission_id=submission_id, submission_dataset_id=s_dataset.id)
    submission_member.save()