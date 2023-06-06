#!/opt/bin/python3

from flask import request, abort
from flask_restful import Resource
import requests

import json
import os,sys

import geardb

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

from pathlib import Path
abs_path_gear = Path(__file__).resolve().parents[2]
cgi_path = str(abs_path_gear.joinpath('cgi'))
sys.path.append(cgi_path)

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

pull_from_gcp = import_from_file("pull_from_gcp", f"{cgi_path}/nemoarchive_pull_gcp_files_to_vm.cgi")
validate_mdata = import_from_file("validate_mdata", f"{cgi_path}/nemoarchive_validate_metadata.cgi")
write_h5ad = import_from_file("write_h5ad", f"{cgi_path}/nemoarchive_write_h5ad.cgi")
make_display = import_from_file("make_display", f"{cgi_path}/nemoarchive_make_default_display.cgi")

def should_step_run(s_dataset, db_step):
    # Step should only run if it is pending.
    # GET /submission has an option to reset incomplete steps back to pending if enabled
    return s_dataset[db_step] == "pending"

def get_db_status(s_dataset):
    return {
                "pulled_to_vm": s_dataset.pulled_to_vm_status
                , "convert_metadata": s_dataset.convert_metadata_status
                , "convert_to_h5ad": s_dataset.convert_to_h5ad_status
                , "make_tsne": s_dataset.make_tsne_status
            }

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

def save_submission_dataset(dataset_id, identifier, is_restricted):
    # Dataset is a foreign key in SubmissionDataset so we need to ensure we do not duplicate
    dataset = geardb.get_dataset_by_id(dataset_id)
    if not dataset:
        insert_minimal_dataset(dataset_id, identifier)

    s_dataset = geardb.SubmissionDataset(dataset_id=dataset_id,
                    is_restricted=is_restricted, nemo_identifier=identifier,
                    pulled_to_vm_status="pending", convert_metadata_status="pending", convert_to_h5ad_status="pending",
                    make_tsne_status="pending", log_message=""
                    )
    s_dataset.save()
    s_dataset.get_dataset_info()
    return s_dataset

def save_submission_member(submission_id:geardb.Submission, s_dataset:geardb.SubmissionDataset):
    submission_member = geardb.SubmissionMember(submission_id=submission_id, submission_dataset_id=s_dataset.id)
    if not submission_member.does_submission_member_exist():
        submission_member.save()

def insert_minimal_dataset(dataset_id, identifier):
    """
    Insert a minimally viable dataset. We will populate it with metadata later.

    This is purely to get the SubmissionDataset row inserted.
    """

    add_dataset_sql = """
    INSERT INTO dataset (id, share_id, owner_id, title, organism_id, date_added, load_status)
    VALUES              (%s, %s,       %s,       "PRE-STAGE NEMO IMPORT - %s",    1, NOW(), "pending")
    """

    (before, sep, share_id) = identifier.rpartition("nemo:")

    # Insert dataset info to database
    conn = geardb.Connection()
    cursor = conn.get_cursor()
    cursor.execute(add_dataset_sql, (dataset_id, share_id, geardb.find_importer_id().id, share_id))
    cursor.close()
    conn.commit()
    conn.close()

def make_tsne_display(dataset_id, session_id, category=None, gene=None):
    result = make_display.make_default_display(dataset_id, session_id, category, gene)

    plot_type = "tsne_static"
    label = "nemoanalytics import default plot"
    user = geardb.get_user_id_from_session_id(session_id)

    params = {
        "user_id": user.id
        , "dataset_id": dataset_id
        , "plotly_config": result["plot_config"]
        , "plot_type": plot_type
        , "label": label
    }
    result = requests.post("https://localhost/cgi/save_dataset_display.cgi", json=params, verify=False)
    result.raise_for_status()
    decoded_result = result.json()
    display_id = decoded_result["display_id"]
    params = {
        "user_id": user.id
        , "dataset_id": dataset_id
        , "display_id": display_id
    }
    result = requests.post("https://localhost/cgi/save_default_display.cgi", json=params, verify=False)
    result.raise_for_status()
    result = result.json()
    return result

def submission_dataset_callback(dataset_id, metadata, session_id, action=None, category=None, gene=None):
    """Run all steps to import a single dataset."""

    if action == "make_display":
        result = make_tsne_display(dataset_id, session_id, category, gene)
        result["self"] = request.path
        return result

    if action == "import":
        dataset_mdata = metadata["dataset"]
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        #NOTE: Each of the CGI scripts will control loading/failed/complete status of their process

        ###
        db_step = "pulled_to_vm_status"    # step name in database
        if should_step_run(s_dataset, db_step):
            # Component = file format type
            component_files = dataset_mdata["component_fields"]
            for component in component_files:
                bucket_path = dataset_mdata[component]
                result = pull_from_gcp.pull_gcp_files_to_vm(bucket_path, dataset_id)
                if not result["success"]:
                    result["self"] = request.path
                    return result

        ###
        db_step = "convert_metadata_status"
        if should_step_run(s_dataset, db_step):
            result = validate_mdata.validate_metadata(dataset_id, session_id, metadata)
            result["self"] = request.path
            if not result["success"]:
                result["self"] = request.path
                return result

        ###
        db_step = "convert_to_h5ad_status"
        if should_step_run(s_dataset, db_step):
            filetype = dataset_mdata["filetype"]
            result = write_h5ad.run_write_h5ad(dataset_id, filetype)
            if not result["success"]:
                result["self"] = request.path
                return result

        ###
        db_step = "make_tsne_status"
        if should_step_run(s_dataset, db_step):
            result = make_tsne_display(dataset_id, session_id, category, gene)

        result["self"] = request.path
        return result

    abort(400)

class SubmissionDataset(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id, dataset_id):
        """Retrieve a particular submission dataset based on ID."""
        result = {"self":request.path, "href":request.path}
        try:
            s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
            result["success"] = True
            result["dataset_id"] = s_dataset.dataset_id
            result["identifier"] = s_dataset.nemo_identifier
            result["share_id"] = s_dataset.dataset.share_id
            result["status"] = get_db_status(s_dataset)
            return result
        except:
            abort(404)

    def post(self, submission_id, dataset_id):
        # Must have a gEAR account to create submissions
        session_id = request.cookies.get('gear_session_id')

        req = request.get_json()
        metadata = req.get("metadata")
        action = req.get("action", None)
        category = req.get("category", None)
        gene = req.get("gene", None)

        # Create a messaging queue if necessary. Make it persistent across the lifetime of the Flask server.
        # Channels will be spawned during each task.
        if this.servercfg['nemoarchive_import']['queue_enabled'].startswith("1"):
            import gearqueue
            host = this.servercfg['nemoarchive_import']['queue_host']
            try:
                # Connect as a blocking RabbitMQ publisher
                connection = gearqueue.Connection(host=host, publisher_or_consumer="publisher")
            except Exception as e:
                return {
                    'success': -1
                , 'message': str(e)
                }
            # Connect as a blocking RabbitMQ publisher
            with connection:
                connection.open_channel()
                task_finished = False
                response = {}
                def _on_response(channel, method_frame, properties, body):
                    nonlocal task_finished
                    nonlocal response
                    task_finished = True
                    response = json.loads(body)
                    print("[x] - Received response for submission_dataset {}".format(payload["dataset_id"]), file=sys.stderr)

                # Create a "reply-to" consumer
                # see https://pika.readthedocs.io/en/stable/examples/direct_reply_to.html?highlight=reply_to#direct-reply-to-example
                try:
                    connection.replyto_consume(
                        on_message_callback=_on_response
                    )
                except Exception as e:
                    return {
                        'success': -1
                    , 'message': str(e)
                    }

                # Create the publisher
                payload = {}
                payload["dataset_id"] = dataset_id
                payload["metadata"] = metadata
                payload["session_id"] = session_id
                payload["action"] = action
                payload["category"] = category
                payload["gene"] = gene

                try:
                    connection.publish(
                        queue_name="nemoarchive_import"
                        , reply_to="amq.rabbitmq.reply-to"
                        , message=payload   # method dumps JSON
                    )
                    print("[x] Requesting for submission_dataset {}".format(dataset_id), file=sys.stderr)
                except Exception as e:
                    return {
                        'success': -1
                    , 'message': str(e)
                    }
                # Wait for callback to finish, then return the response
                while not task_finished:
                    pass
                print("[x] sending payload response for submission_dataset {} back to client".format(dataset_id), file=sys.stderr)
                return response
        else:
            return submission_dataset_callback(dataset_id, metadata, session_id, category, gene)

class SubmissionDatasetMember(Resource):
    def put(self, submission_id, dataset_id):
        result = {"success": False, "message":"", "self": request.path}
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        # Insert submission members and create new submisison datasets if they do not exist
        try:
            save_submission_member(submission_id, s_dataset)
        except Exception as e:
            print(str(e), file=sys.stderr)
            result["message"] = f"Could not save dataset {dataset_id} as a new SubmissionMember in database"
            s_dataset.update_downstream_steps_to_cancelled()
            result["status"] = get_db_status(s_dataset)
            return result

class SubmissionDatasetStatus(Resource):
    def put(self, submission_id, dataset_id):
        req = request.get_json()
        action = req.get("action")

        result = {"success": False, "self": request.path}
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        if action == "reset_steps":
            # Before initializing a new submission, we check this route to see if the current dataset exists
            # If it does and is not complete or loading, need to reset the steps so it can be resumed.
            s_dataset.reset_incomplete_steps()
            result["status"] = get_db_status(s_dataset)
            result["success"] = True
            return result
        abort(405)

class SubmissionDatasets(Resource):
    """Requests to deal with multiple submissions, including creating new ones"""

    def post(self, submission_id):
        """Create a new submission in the database."""
        req = request.get_json()
        dataset_id = req.get("dataset_id")
        identifier = req.get("identifier")
        is_restricted = req.get("is_restricted")

        result = {"self":request.path, "message":"", "success": False}

        if not identifier:
            result["message"] = "NeMO Identifier empty. Cannot load."
            result["status"] = {
                "pulled_to_vm": "canceled"
                , "convert_metadata": "canceled"
                , "convert_to_h5ad": "canceled"
                , "make_tsne": "canceled"
                }
            return result

        submission_dataset_params = {
            "is_restricted": is_restricted
            , "dataset_id": dataset_id
            , "identifier": identifier
            }

        # Dataset may already be inserted from a previous submission
        try:
            s_dataset = add_submission_dataset(submission_dataset_params)
        except Exception as e:
            print(str(e), file=sys.stderr)
            result["message"] = f"Could not save dataset {dataset_id} as a new SubmissionDataset in database"
            result["status"] = {
                "pulled_to_vm": "canceled"
                , "convert_metadata": "canceled"
                , "convert_to_h5ad": "canceled"
                , "make_tsne": "canceled"
                }
            return result

        result["dataset_id"] = s_dataset.dataset_id
        result["identifier"] = s_dataset.nemo_identifier
        result["share_id"] = s_dataset.dataset.share_id
        result["status"] = get_db_status(s_dataset)

        result["href"] = result["self"] + f"/${dataset_id}" # This can be GET-able

        result["success"] = True
        return result


