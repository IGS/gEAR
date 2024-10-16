#!/opt/bin/python3

from flask import request
from flask_restful import Resource, abort
import os,sys
import asyncio, aiohttp

from pathlib import Path

# Add the cgi directory to the path
abs_path_gear = Path(__file__).resolve().parents[2]
cgi_path = str(abs_path_gear.joinpath('cgi'))
sys.path.append(cgi_path)

lib_path = Path('..') / '..' / 'lib'
lib_path = str(lib_path.resolve())
sys.path.append(lib_path)
import geardb

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

# Need extra stuff to import .cgi scripts due to not having the .py extensions
from importlib.util import spec_from_loader, module_from_spec
from importlib.machinery import SourceFileLoader

# Source -> https://gist.github.com/mportesdev/afb2ec26021ccabee0f67d6f7d18be3f
# This is used to import .cgi files
def import_from_file(module_name, file_path):
    loader = SourceFileLoader(module_name, file_path)
    spec = spec_from_loader(module_name, loader)
    module = module_from_spec(spec)

    sys.modules[module_name] = module    # this step is not necessary, but it's
                                         # probably not a bad idea to have a manually
                                         # imported module stored in sys.modules
    spec.loader.exec_module(module)
    return module

add_display_to_layout = import_from_file("add_display_to_layout", f"{cgi_path}/add_display_to_layout.cgi")


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
    layout_share_id = layout.share_id if layout else None

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
        text += " It looks like every dataset was successfully imported."
    elif result["partial_success"]:
        text += f" While some datasets were successfully imported, we encountered some issues importing {num_failures} of the datasets."
    else:
        text += f" Unfortunately, it appears that all {num_datasets} datasets for this submission ran into issues during the import process."

    text += f"\n\nYou can view the status of all datasets at {domain_url}/nemoarchive_import/import.html?submission_id={submission_id}."
    if result["partial_success"]:
        text += " Each successfully imported dataset has a link to view the dataset on the gene expression search page."
        text += f" We performed a basic Seurat analysis on the dataset and created a tSNE display, but we encourage you to curate your own displays to support your own research (see {domain_url}/manual.html?doc=curation for documentation)."

    if layout_share_id:
        text += f"\n\nYou can also visualize all datasets as a collection by going to {domain_url}/p?l=${layout_share_id} (you will first need to search for a gene to show the displays)."

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

def get_submission(submission_id, fail_on_exist=False):
    submission = geardb.get_submission_by_id(submission_id)
    if fail_on_exist and submission
        abort(400, message=f"Submission already exists for this ID {submission_id}.")
    if not submission:
        abort(404, message=f"Submission id {submission_id} does not exist.")
    return submission

class Submission(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id):
        """Retrieve a particular submission based on ID."""
        url_path = request.root_path + request.path
        result = {
            "self": url_path
            , "success":False
            , "datasets":{}
            , "layout_share_id":None
            , "is_finished":False
            , "is_submitter":False
            }

        session_id = request.cookies.get('gear_session_id')
        user_id = geardb.get_user_id_from_session_id(session_id)
        if not user_id:
            # You can view a particular submission without being logged in
            user_id = -1

        submission = get_submission(submission_id)

        if user_id == submission.user_id:
            result["is_submitter"] = True

        layout = submission.get_layout_info()
        if layout:
            result["layout_share_id"] = layout.share_id
            result["collection_name"] = layout.label

        result["is_finished"] = bool(submission.is_finished)

        sd = geardb.SubmissionDatasetCollection()
        submission.datasets = sd.get_by_submission_id(submission_id)

        for s_dataset in submission.datasets:
            dataset_id = s_dataset.dataset_id
            result["datasets"][dataset_id] = {"href":f"/api/submissions/{submission_id}/datasets/{dataset_id}"}
        result["success"] = True
        return result

    def delete(self, submission_id):
        """Delete the existing submission from the database."""
        submission = get_submission(submission_id)
        if not submission:
            abort(404, message=f"Submission id {submission_id} does not exist.")
        submission.remove()

    def post(self, submission_id):
        """Start the import process for this submission."""
        url_path = request.root_path + request.path

        session_id = request.cookies.get('gear_session_id')
        # Must have a gEAR account to upload datasets
        user = geardb.get_user_from_session_id(session_id)

        req = request.get_json()
        action = req.get("action")

        submission = get_submission(submission_id)
        submission.datasets = geardb.SubmissionDatasetCollection().get_by_submission_id(submission_id=submission_id)

        if action == "import":
            async def import_dataset(s_dataset):

                dataset_id = s_dataset.dataset_id
                identifier = s_dataset.nemo_identifier

                result = {"success": False, "dataset_id": dataset_id}

                async with aiohttp.ClientSession() as session:
                    # Start the import process for each dataset
                    params = {"identifier": identifier, "action": "import"}
                    try:
                        url = f"http://localhost/{url_path}" + f"/datasets/{dataset_id}"
                        async with session.post(url
                                , json=params
                                , cookies={"gear_session_id":session_id}
                                , verify_ssl=False
                                , raise_for_status=True) as resp:
                            import_result = await resp.json()
                    except Exception as e:
                        print(str(e), file=sys.stderr)
                        return result

                result = import_result # should have "success" = True in here

                # Let's save the display to the submission layout while we are at it
                result = add_display_to_layout.add_display_to_layout(session_id, result['share_id'], result['display_id'], 12, 1)
                if not result["success"]:
                    raise Exception("Write H5AD step failed")

                result["filetype"] = result["dataset"]["filetype"]
                result["dataset_id"] = dataset_id

                return result


            coros = [import_dataset(s_dataset) for s_dataset in submission.datasets]
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            tasks = asyncio.gather(*coros, return_exceptions=True)  # We don't want exceptions to crash the other tasks
            results = loop.run_until_complete(tasks)

            all_success = True  # Every dataset import was successful
            partial_success = False # At least one dataset import was successful
            num_failures = 0
            dataset_results = {}

            for res in results:
                if isinstance(res, Exception):
                    message = res
                    print(str(message), file=sys.stderr)
                    #res = {"success":False, "message": message}

                dataset_id = res["dataset_id"]
                dataset_results[dataset_id] = res

                if res["success"]:
                    partial_success = True
                else:
                    num_failures += 1
                    all_success = False

            submission.save_change(attribute="is_finished", value=1)

            result = {"self":url_path
                    , "datasets":dataset_results
                    , "all_success":all_success
                    , "partial_success":partial_success
                    , "num_failures": num_failures
                    }

            if submission.email_updates:
                # send email about submission status
                send_email(result, submission, user.email)
                # Disable email updates in case people refresh the "import" link.
                # If they want another email, they must click the "Send Updates" button on page
                submission.save_change(attribute="email_updates", value=0)

            return result
        abort(400, message=f"No action or unknown action specified for submission id {submission_id}.")

class SubmissionEmail(Resource):
    def put(self, submission_id):
        """Update email updates."""
        if not submission_id:
            abort(404, message=f"Submission id {submission_id} does not exist.")

        result = {"success": False, "message":""}

        try:
            submission = get_submission(submission_id)
            submission.save_change(attribute="email_updates", value=1)
            result["success"] = True
        except Exception as e:
            result["message"] = str(e)

        return result

class Submissions(Resource):
    """Requests to deal with multiple submissions, including creating new ones"""

    def post(self):
        """Create a new empty submission in the database."""
        url_path = request.root_path + request.path

        # Must have a gEAR account to create submissions
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)

        result = {"self": url_path, "success":False, "datasets":{}, "layout_share_id":None, "message":""}

        if not user:
            result["message"] = "User not logged in."
            return result

        req = request.get_json()
        submission_id = req.get("submission_id")
        is_restricted = req.get("is_restricted")


        submission = get_submission(submission_id, fail_on_exist=True)

        try:
            conn = geardb.Connection()
            cursor = conn.get_cursor()

            qry = """
                INSERT INTO submission (id, user_id, is_finished, is_restricted, date_added)
                VALUES (%s, %s, 0, %s, NOW())
            """

            cursor.execute(qry, (submission_id, user.id, is_restricted))

            # Save submission as a layout
            layout_name = f"Submission {submission_id}"
            layout = geardb.Layout(user_id=user.id, label=layout_name,
                                is_current=1, members=None)
            layout.save()
            result['layout_share_id'] = layout.share_id
            cursor.close()
            conn.commit()   # So Submission transaction is saved

            cursor = conn.get_cursor()
            submission = get_submission(submission_id)
            if not submission:
                abort(400, message=f"New submission {submission_id} could not be found.")
            submission.save_change(attribute="layout_id", value=layout.id)
            cursor.close()
            conn.commit()
            conn.close()
        except Exception as e:
            print(str(e), file=sys.stderr)
            abort(400, message="Could not save new submission to database")

        result["href"] = result["self"] + f"/{submission_id}" # This can be GET-able

        result["success"] = True
        result["message"] = "Submission database entries have been populated"
        return result

