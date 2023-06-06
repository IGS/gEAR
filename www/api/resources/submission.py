#!/opt/bin/python3

from flask import request, abort
from flask_restful import Resource
import os,sys
import requests

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

def send_email(result, submission, user_email):
    """Send user an email about the completion status of their dataset imports."""

    import smtplib
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText
    sender = this.servercfg['email_sender']['address']
    password = this.servercfg['email_sender']['password']

    domain_url = geardb.domain_url
    domain_short_label = geardb.domain_short_label

    submission_id = submission.id
    layout = submission.get_layout_info()
    layout_share_id = layout.share_id

    # https://docs.python.org/3/library/email-examples.html
    # user = email
    msg = MIMEMultipart('alternative')
    msg['Subject'] = f'{domain_short_label} - Submission {submission_id} import complete'.format(domain_short_label)
    msg['From'] = sender
    msg['To'] = user_email

    num_failures = result["num_failures"]
    num_datasets = len(submission.datasets)

    text = f"This message is to inform you that your recent dataset import on {domain_short_label} with submission ID {submission_id} has finished."
    if result["all_success"]:
        text += "It looks like every dataset was successfully imported."
    elif result["partial_success"]:
        text += f"While some datasets were successfully imported, we encountered some issues importing {num_failures} of the datasets."
    else:
        text += f"Unfortunately, it appears that all {num_datasets} datasets for this submission ran into issues during the import process."

    text += f"\n\nYou can view the status of all datasets at {domain_url}/nemoarchive_import/import.html?submission_id={submission_id}. "
    if result["partial_success"]:
        text += "Each successfully imported dataset has a link to view the dataset on the gene expression search page. "
        text += "We performed a basic Seurat analysis on the dataset and created a tSNE display, but we encourage you to curate your own displays to support your own research (see https://umgear.org/manual.html?doc=curation for documentation). "

        text += f"\n\nYou can also view all dataset displays as a collection by going to {domain_url}/p?l=${layout_share_id} (you will need to search for a gene to show the displays)."

    text += f"\n\nIf you have any questions about your import, feel free to reach us at {domain_url}/contact.html"

    text += f"\n\nSigned,\n The {domain_short_label} team"

    part1 = MIMEText(text, 'plain')

    msg.attach(part1)

    # http://stackoverflow.com/a/17596848/2900840
    s = smtplib.SMTP('smtp.gmail.com:587')
    s.ehlo()
    s.starttls()
    s.login(sender, password)

    s.sendmail( sender, user_email, msg.as_string() )
    s.quit()

class Submission(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id):
        """Retrieve a particular submission based on ID."""
        result = {"self": request.path, "success":False, "datasets":{}, "layout_share_id":None}

        session_id = request.cookies.get('session_id')
        # Must have a gEAR account to upload datasets
        user = geardb.get_user_from_session_id(session_id)

        submission = geardb.get_submission_by_id(submission_id)
        if not submission:
            return result

        layout = submission.get_layout_info()
        if layout:
            result["layout_share_id"] = layout.share_id

        sd = geardb.SubmissionDatasetCollection()
        submission.datasets = sd.get_by_submission_id(submission_id)

        for s_dataset in submission.datasets:
            dataset_id = s_dataset.dataset_id
            result["datasets"][dataset_id] = {"href":f"/api/submissions/{submission_id}/datasets/{dataset_id}"}
        result["success"] = True
        return result

    def post(self, submission_id):
        """Start the import process for this submission."""
        session_id = request.cookies.get('session_id')
        # Must have a gEAR account to upload datasets
        user = geardb.get_user_from_session_id(session_id)

        req = request.get_json()
        file_metadata = req.get("file_metadata")  # Metadata from neo4j request
        sample_metadata = req.get("sample_metadata")
        action = req.get("action")

        submission = geardb.get_submission_by_id(submission_id)
        submission_datasets = geardb.SubmissionDatasetCollection()
        submission.datasets = submission_datasets.get_by_submission_id(submission_id=submission_id)

        if action == "import":
            dataset_results = {}

            all_success = True  # Every dataset import was successful
            partial_success = False # At least one dataset import was successful
            num_failures = 0

            for s_dataset in submission_datasets:

                # Get first instance of entity match dataset ID or sample ID
                file_attributes = next(filter(lambda entity: entity["id"] == dataset_id, file_metadata))["attributes"]
                sample_attributes = next(filter(lambda entity: entity["name"] == file_attributes["sample_id"], sample_metadata))["attributes"]

                all_attributes = {
                    "dataset": file_attributes
                    , "sample": sample_attributes
                    }

                dataset_id = s_dataset["dataset_id"]
                identifier = s_dataset["identifier"]

                # Call NeMO Archive assets API using identifier
                #result = requests.post(f"https://nemoarchive.org/asset/derived/{identifier}")
                result = requests.post(f"/api/mock_identifier/{identifier}")
                result.raise_for_status()
                decoded_result = result.json()
                api_metadata = decoded_result["metadata"]

                # TODO: Eventually get everything from API
                all_metadata = {
                    "dataset": {**all_attributes["dataset"], **api_metadata["dataset"]}
                    , "sample": {**all_attributes["sample"], **api_metadata["sample"]}

                }

                # Start the import process for each dataset
                params = {"metadata": all_metadata, "action": "import"}
                result = requests.post(request.path + f"datasets/${dataset_id}", json=params, verify=False)
                result.raise_for_status()
                decoded_result = result.json()

                dataset_results[dataset_id] = decoded_result
                dataset_results[dataset_id]["filetype"] = all_metadata["dataset"]["filetype"]

                if dataset_results[dataset_id]["success"]:
                    partial_success = True
                else:
                    num_failures += 1
                    all_success = False

            result = {"self":request.path
                    , "dataset":dataset_results
                    , "all_success":all_success
                    , "partial_success":partial_success
                    , "num_failures": num_failures
                    }

            if submission.email_updates:
                # send email about submission status
                send_email(result, submission, user.email)

            return result
        abort(400)

class SubmissionEmail(Resource):
    def put(self, submission_id):
        """Update email updates."""
        if not submission_id:
            abort(404)  # "Submissiom ID not found for this request."

        result = {"success": False, "message":""}

        try:
            submission = geardb.get_submission_by_id(submission_id)
            submission.save_change(attribute="email_updates", value=1)
            result["success"] = True
        except Exception as e:
            result["message"] = str(e)

        return result

class Submissions(Resource):
    """Requests to deal with multiple submissions, including creating new ones"""

    def post(self):
        """Create a new empty submission in the database."""
        # Must have a gEAR account to create submissions
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)

        req = request.get_json()
        submission_id = req.get("submission_id")
        is_restricted = req.get("is_restricted")

        result = {"self": request.path, "success":False, "datasets":{}, "layout_share_id":None, "message":""}

        submission = geardb.get_submission_by_id(submission_id)
        if submission:
            abort(400, message="Submission already exists for this ID.")

        try:
            conn = geardb.Connection()
            cursor = conn.get_cursor()

            qry = """
                INSERT INTO submission (id, user_id, is_finished, is_restricted, date_added)
                VALUES (%s, %s, 0, %s, NOW())
            """

            cursor.execute(qry, (submission_id, user.id, is_restricted))
            cursor.close()
            conn.commit()
            conn.close()
        except:
            abort(400, message="Could not save new submission to database")

        result["href"] = result["self"] + f"/${submission_id}" # This can be GET-able

        result["success"] = True
        result["message"] = "Submission database entries have been populated"
        return result

