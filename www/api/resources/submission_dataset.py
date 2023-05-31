#!/opt/bin/python3

from flask import request, abort
from flask_restful import Resource

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

cgi_path = os.path.abspath(os.path.join("..", "cgi"))
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

pull_from_gcp = import_from_file("pull_from_gcp", "../cgi/nemoarchive_pull_gcp_files_to_vm.cgi")
validate_mdata = import_from_file("validate_mdata", "../cgi/nemoarchive_validate_metadata.cgi")
write_h5ad = import_from_file("write_h5ad", "../cgi/nemoarchive_write_h5ad.cgi")
make_display = import_from_file("make_display", "../cgi/nemoarchive_make_default_display.cgi")

def should_step_run(s_dataset, db_step):
    # Step should only run if it is pending.
    # GET /submission has an option to reset incomplete steps back to pending if enabled
    return s_dataset[db_step] == "pending"

def submission_dataset_callback(dataset_id, metadata, session_id, action=None, category=None, gene=None fh):
    """Run all steps to import a single dataset."""

    if action == "make_display":
        result = make_display.make_default_display(dataset_id, session_id, category, gene)
        result["self"] = request.path
        return result

    dataset_mdata = metadata["dataset"]
    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

    ###
    db_step = "pulled_to_vm_status"    # step name in database
    if should_step_run(s_dataset, db_step):
        s_dataset.save_change(attribute=db_step, value="loading")

        # Component = file format type
        component_files = dataset_mdata["component_fields"]
        for component in component_files:
            bucket_path = dataset_mdata[component]
            result = pull_from_gcp.pull_gcp_files_to_vm(bucket_path, dataset_id)
            if not result["success"]:
                result["self"] = request.path
                s_dataset.update_downstream_steps_to_cancelled(attribute=db_step)
                return result

    ###
    db_step = "convert_metadata_status"
    if should_step_run(s_dataset, db_step):
        s_dataset.save_change(attribute=db_step, value="loading")
        result = validate_mdata.validate_metadata(dataset_id, session_id, metadata)
        result["self"] = request.path
        if not result["success"]:
            result["self"] = request.path
            s_dataset.update_downstream_steps_to_cancelled(attribute=db_step)
            return result

    ###
    db_step = "convert_to_h5ad_status"
    if should_step_run(s_dataset, db_step):
        s_dataset.save_change(attribute=db_step, value="loading")
        filetype = dataset_mdata["filetype"]
        result = write_h5ad.run_write_h5ad(dataset_id, filetype)
        if not result["success"]:
            result["self"] = request.path
            return result

    ###
    result = make_display.make_default_display(dataset_id, session_id, category, gene)
    result["self"] = request.path
    return result

class SubmissionDataset(Resource):
    """Requests to deal with a single submission"""

    def get(self, dataset_id):
        """Retrieve a particular submission dataset based on ID."""
        pass

    def post(self, dataset_id):
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
            return submission_dataset_callback(dataset_id, metadata, session_id, category, gene, sys.stderr)



class SubmissionDatasets(Resource):
    """Requests to deal with multiple submissions, including creating new ones"""

    def post(self):
        """Create a new submission in the database."""
        # Must have a gEAR account to create submissions
        session_id = request.cookies.get('gear_session_id')

        req = request.get_json()
        metadata = req.get("metadata")
        abort(501, "POST /submission_datasets has not been implemented yet.")

