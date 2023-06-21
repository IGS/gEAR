#!/opt/bin/python3

from flask import request
from flask_restful import Resource, abort

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
    return getattr(s_dataset, db_step) == "pending"

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
        s_dataset.save_change(attribute="log_message", value="")
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

def submission_dataset_callback(dataset_id, metadata, session_id, url_path, action=None, category=None, gene=None):
    """Run all steps to import a single dataset."""

    #TODO: How to run as non-root when RabbitMQ is running this (at least to move h5ad/json to final area)

    result = {"success" : False}

    if action == "make_display":
        try:
            result = make_display.make_default_display(dataset_id, session_id, category, gene)
            result["self"] = url_path
            if not result["success"]:
                raise Exception("Make tSNE step failed")
        except Exception as e:
            result["message"] = str(e)
        finally:
            result["self"] = url_path
            return result


    if action == "import":
        dataset_mdata = metadata["dataset"]
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        #NOTE: Each of the CGI scripts will control loading/failed/complete status of their process

        ###
        try:
            db_step = "pulled_to_vm_status"    # step name in database
            if should_step_run(s_dataset, db_step):
                # Component = file format type
                component_files = dataset_mdata["component_fields"]
                for component in component_files:
                    bucket_path = dataset_mdata[component]
                    result = pull_from_gcp.pull_gcp_files_to_vm(bucket_path, dataset_id)
                    if not result["success"]:
                        raise Exception("Pull GCP Files step failed")

            ###
            db_step = "convert_metadata_status"
            if should_step_run(s_dataset, db_step):
                result = validate_mdata.validate_metadata(dataset_id, session_id, metadata)
                if not result["success"]:
                    raise Exception("Validate Metadata step failed")

            ###
            db_step = "convert_to_h5ad_status"
            if should_step_run(s_dataset, db_step):
                filetype = dataset_mdata["filetype"]
                result = write_h5ad.run_write_h5ad(dataset_id, filetype)
                if not result["success"]:
                    raise Exception("Write H5AD step failed")

            ###
            db_step = "make_tsne_status"
            if should_step_run(s_dataset, db_step):
                result = make_display.make_default_display(dataset_id, session_id, category, gene)
                if not result["success"]:
                    raise Exception("Make tSNE step failed")
            else:
                # Nothing needs to be run.
                result["success"] = True
        except Exception as e:
            result["message"] = str(e)

            # If something happened that prevents the above scripts from changing step status to "failed",
            # such as RabbitMQ consumer crashing, resolve that here
            loading_step = s_dataset.find_loading_step()
            if loading_step:
                s_dataset.save_change(attribute=loading_step, value="failed")
                s_dataset.update_downstream_steps_to_canceled(attribute=loading_step)

        finally:
            result["self"] = url_path
            return result

    result["message"] = "No action operation requested"
    return result


class SubmissionDataset(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id, dataset_id):
        """Retrieve a particular submission dataset based on ID."""
        url_path = request.root_path + request.path

        result = {"self":url_path, "href":url_path}

        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
        if not s_dataset:
            abort(404)

        result["success"] = True
        result["dataset_id"] = s_dataset.dataset_id
        result["identifier"] = s_dataset.nemo_identifier
        result["share_id"] = s_dataset.dataset.share_id
        result["status"] = get_db_status(s_dataset)
        result["message"] = s_dataset.log_message
        return result

    def post(self, submission_id, dataset_id):
        # Must have a gEAR account to create submissions
        url_path = request.root_path + request.path
        session_id = request.cookies.get('gear_session_id')

        req = request.get_json()
        metadata = req.get("metadata")
        action = req.get("action", None)
        category = req.get("category", None)
        gene = req.get("gene", None)

        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
        if not s_dataset:
            abort(404)

        # Create a messaging queue if necessary. Make it persistent across the lifetime of the Flask server.
        # Channels will be spawned during each task.
        if this.servercfg['nemoarchive_import']['queue_enabled'].startswith("1"):
            import gearqueue
            host = this.servercfg['nemoarchive_import']['queue_host']
            try:
                # Connect as a blocking RabbitMQ publisher
                connection = gearqueue.Connection(host=host, publisher_or_consumer="publisher")
            except Exception as e:
                s_dataset.save_change(attribute="log_message", value="Could not establish connection as RabbitMQ publisher.")
                print(f"{request.path} - ERROR - {str(e)}", file=sys.stderr)
                abort(500)

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
                    s_dataset.save_change(attribute="log_message", value="Could not receive reply response from RabbitMQ consumer.")
                    print(f"{request.path} - ERROR - {str(e)}", file=sys.stderr)
                    abort(500)

                # Create the publisher
                payload = {}
                payload["dataset_id"] = dataset_id
                payload["metadata"] = metadata
                payload["session_id"] = session_id
                payload["url_path"] = url_path
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
                    s_dataset.save_change(attribute="log_message", value="Could not publish payload to RabbitMQ consumer.")
                    print(f"{request.path} - ERROR - {str(e)}", file=sys.stderr)
                    abort(500)

                # Wait for callback to finish, then return the response
                while not task_finished:
                    pass
                print("[x] sending payload response for submission_dataset {} back to client".format(dataset_id), file=sys.stderr)
                if not result["success"]:
                    print(result.get("message", "Something went wrong."), file=sys.stderr)
                    abort(500)
                return result

        else:
            result =  submission_dataset_callback(dataset_id, metadata, session_id, url_path, action, category, gene)
            if not result["success"]:
                print(result.get("message", "Something went wrong."), file=sys.stderr)
                abort(500)
            return result

class SubmissionDatasetMember(Resource):
    def put(self, submission_id, dataset_id):
        url_path = request.root_path + request.path

        session_id = request.cookies.get('gear_session_id')
        user_id = geardb.get_user_id_from_session_id(session_id)

        result = {"success": False, "message":"", "self": url_path}

        if not user_id:
            result["message"] = "User not logged in."
            return result

        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        if not s_dataset:
            abort(404)

        # Insert submission members and create new submisison datasets if they do not exist
        try:
            save_submission_member(submission_id, s_dataset)

            # ? Should this go in a separate API function, like /layout
            # Let's save the dataset to the submission layout while we are at it
            submission = geardb.get_submission_by_id(submission_id)
            if not submission:
                abort(404)

            layout = submission.get_layout_info()
            # make sure the user owns the layout
            gpos = len(layout.members) + 1

            # This shoud never be an issue here (copied from add_dataset_to_layout.cgi)
            if user_id == layout.user_id:
                lm = geardb.LayoutMember(dataset_id=dataset_id, grid_position=gpos, grid_width=4, mg_grid_width=12)
                layout.add_member(lm)

        except Exception as e:
            print(str(e), file=sys.stderr)
            result["message"] = f"Could not save dataset {dataset_id} as a new SubmissionMember in database"
            s_dataset.update_downstream_steps_to_canceled()
            result["status"] = get_db_status(s_dataset)
            return result

class SubmissionDatasetStatus(Resource):
    def put(self, submission_id, dataset_id):
        req = request.get_json()
        action = req.get("action")
        url_path = request.root_path + request.path

        result = {"success": False, "self": url_path}
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        if not s_dataset:
            abort(404)

        if action == "reset_steps":
            # Before initializing a new submission, we check this route to see if the current dataset exists
            # If it does and is not complete or loading, need to reset the steps so it can be resumed.
            s_dataset.reset_incomplete_steps()
            s_dataset.save_change(attribute="log_message", value="")
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

        # NOTE: Eventually the passed in dataset UUID may be abandoned, so
        # we may have to create one and check for the existance of the nemo identifier

        result = {"self":request.url, "message":"", "success": False}

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

        result["href"] = result["self"] + f"/{dataset_id}" # This can be GET-able

        result["success"] = True
        return result


