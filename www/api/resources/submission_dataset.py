#!/opt/bin/python3

from flask import request
from flask_restful import Resource, abort

import json
import os,sys
import requests

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
pull_from_http = import_from_file("pull_from_http", f"{cgi_path}/nemoarchive_pull_http_files_to_vm.cgi")
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
                , "make_umap": s_dataset.make_umap_status
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

def get_submission_dataset(dataset_id) -> geardb.SubmissionDataset:
    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    if not s_dataset:
        abort(404, message=f"Submission dataset id {dataset_id} does not exist.")
    return s_dataset

def save_submission_dataset(dataset_id, identifier, is_restricted) -> geardb.SubmissionDataset:
    # Dataset is a foreign key in SubmissionDataset so we need to ensure we do not duplicate
    dataset = geardb.get_dataset_by_id(dataset_id)
    if not dataset:
        insert_minimal_dataset(dataset_id, identifier)

    s_dataset = geardb.SubmissionDataset(dataset_id=dataset_id,
                    is_restricted=is_restricted, nemo_identifier=identifier,
                    pulled_to_vm_status="pending", convert_metadata_status="pending", convert_to_h5ad_status="pending",
                    make_umap_status="pending", log_message=""
                    )
    s_dataset.save()
    s_dataset.get_dataset_info()
    return s_dataset

def save_submission_member(submission_id, s_dataset):
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

    result: dict[str] = {"success" : False}

    if action == "make_display":
        try:
            result.update(make_display.make_default_display(dataset_id, session_id, category, gene))
            result["self"] = url_path
            if not result["success"]:
                raise Exception("Make UMAP step failed")
        except Exception as e:
            result["message"] = str(e)
        finally:
            result["self"] = url_path
            return result

    if action == "import":
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)

        #NOTE: Each of the CGI scripts will control loading/failed/complete status of their process

        try:
            result.update(pull_nemoarchive_metadata(s_dataset, metadata["identifier"]))
            if not result["success"]:
                raise Exception("Could not pull metadata from NeMO Archive API")
            dataset_mdata = result["metadata"].get("dataset")

            db_step = "pulled_to_vm_status"    # step name in database
            if should_step_run(s_dataset, db_step):
                bucket_path = dataset_mdata["gcp_uri"]
                if bucket_path:
                    result.update(pull_from_gcp.pull_gcp_files_to_vm(bucket_path, dataset_id))
                    if not result["success"]:
                        raise Exception("Pull GCP Files step failed")
                else:
                    http_path = dataset_mdata["http_uri"]
                    if http_path:
                        result.update(pull_from_http.pull_http_files_to_vm(http_path, dataset_id))
                        if not result["success"]:
                            raise Exception("Pull HTTP Files step failed")
                    else:
                        raise Exception("No GCP or HTTP URI found in metadata")

            ###
            db_step = "convert_to_h5ad_status"
            if should_step_run(s_dataset, db_step):
                filetype = dataset_mdata["filetype"]
                result.update(write_h5ad.run_write_h5ad(dataset_id, filetype))
                if not result["success"]:
                    raise Exception("Write H5AD step failed")

            ###
            db_step = "make_umap_status"
            if should_step_run(s_dataset, db_step):
                result.update(make_display.make_default_display(dataset_id, session_id, category, gene))
                if not result["success"]:
                    raise Exception("Make UMAP step failed")

            # Add display to layout will occur in submission as one display can be added to multiple layouts

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

def pull_nemoarchive_metadata(s_dataset, nemo_id) -> dict:
    """
    Pull metadata from NeMO Archive API
    """

    result = {"success" : False, "metadata":{}}

    # Call NeMO Archive assets API using identifier
    try:
        url = f"https://assets.nemoarchive.org/file/{nemo_id}"
        response = requests.get(url, verify=False)
        response.raise_for_status()

        #url = f"http://localhost/api/mock_identifier/{identifier}"
        api_file_result = response.json()
    except Exception as e:
        s_dataset.save_change(attribute="log_message", value=str(e))
        return result

    if "error" in api_file_result:
        s_dataset.save_change(attribute="log_message", value=api_file_result["error"])
        return result

    if not "access" in api_file_result or api_file_result["access"] != "open":
        s_dataset.save_change(attribute="log_message", value="File is not open access. Cannot import file at this time.")
        return result

    """
    sample_identifier = api_file_result["sample"]
    if not sample_identifier:
        s_dataset.save_change(attribute="log_message", value="No sample identifier found in NeMO Archive API. Cannot get sample metadata.")
        return result

    # Call NeMO Archive assets API using identifier
    try:
        url = f"https://assets.nemoarchive.org/sample/{sample_identifier}"
        response = requests.get(url, verify=False)
        response.raise_for_status()

        #url = f"http://localhost/api/mock_identifier/{identifier}"
        api_sample_result = response.json()
    except Exception as e:
        s_dataset.save_change(attribute="log_message", value=str(e))
        return result

    if "error" in api_sample_result:
        s_dataset.save_change(attribute="log_message", value=api_sample_result["error"])
        return result
    """
    api_sample_result = {}

    # Dataset metadata is used for the dataset entry in the database
    # Sample metadata is not actively used but could be used for curation purposes if linked to counts
    dataset_metadata = process_nemo_assets_api_file_result(api_file_result)
    sample_metadata = process_nemo_assets_api_sample_result(api_sample_result)

    api_metadata = {"dataset":dataset_metadata, "sample":sample_metadata}
    result["success"] = True
    result["metadata"] = api_metadata
    return result

def process_nemo_assets_api_file_result(api_result):
    """
    {
        "id": "nemo:hak-03bgkrw",
        "file_name": "string",
        "data_type": "string",
        "duls": {
            "dul": "string",
            "dul_modifiers": [
            "string"
            ],
            "specific_limits": "string"
        },
        "aliquot": {
            "nemo_id": "nemo:kuk-h42sdnl",
            "library_name": "string",
            "library_type": "aliquot",
            "modality": "string",
            "technique": "string",
            "assay": "string",
            "modality_cv_term_id": "string",
            "technique_cv_term_id": "string",
            "assay_cv_term_id": "string",
            "specimen_type": "string"
        },
        "program": "biccn",
        "grant_short_name": "string",
        "modality": "string",
        "technique": "string",
        "anatomical_regions": [
            {
            "short_name": "string",
            "region_name": "string",
            "cv_term_id": "string"
            }
        ],
        "lab": "string",
        "contact": {
            "name": "string",
            "organization": "string",
            "orcid": "stringstringstrings"
        },
        "contributors": [
            {
            "name": "string",
            "organization": "string",
            "orcid": "stringstringstrings"
            }
        ],
        "access": "open",
        "collections": [
            "string"
        ],
        "md5": "stringstringstringstringstringst",
        "file_format": "string",
        "file_attributes": {
            "name": "string",
            "value": "string",
            "unit": "string",
            "cv_term_id": "string"
        },
        "size": 0,
        "last_modified": "string",
        "manifest_file_urls": [
            {
            "readme": "string",
            "file_location": "string",
            "protocol": "string",
            "file_count": 0,
            "size": 0,
            "url": "string",
            "type": "all"
            }
        ],
        "analysis": [
            {
            "analysis_name": "string",
            "attribute_name": "string",
            "attribute_value": "string"
            }
        ],
        "parent_files": [
            "string"
        ],
        "child_files": [
            "string"
        ],
        "taxa": [
            {
            "name": "string",
            "cv_term_id": "string"
            }
        ],
        "library_pool": "string",
        "alternate_id": "string",
        "sample": "string"
    }
    """

    # ? Eventually we will need to distinguish between file datasets and sample datasets.  Sample-based ones may have multiple entries for an attribute, such as age.

    dataset_metadata = {}

    dataset_metadata["identifier"] = api_result["id"]
    dataset_metadata["contact_name"] = api_result["contact"].get("name")
    dataset_metadata["contact_orcid"] = f'orcid:{api_result["contact"].get("orcid")}'   # Use orcid identifier as email since email is not provided

    dataset_metadata["contact_institute"] = api_result["contact"].get("organization")

    # For title - two options
    # Collection name + file name
    # File name only (may be indecipherable)
    # ? Verify
    dataset_metadata["title"] = api_result["file_name"]
    """
    TEMP TITLE
    grant = attributes["sample"]["project_grant"]
    tissue = attributes["sample"]["tissue_ontology"]
    sex = attributes["sample"]["sex_assigned_at_birth"]
    age = attributes["sample"]["age_value"]
    age_unit = attributes["sample"]["age_unit"]
    json_attributes["value"].append(f"FAKE NAME {grant} - {tissue} - {sex} - {age}  {age_unit}")
    """

    # ? Is there a better summary we can use?
    url = f"https://assets.nemoarchive.org/{api_result['id']}"
    dataset_metadata["summary"] = f"This dataset was derived from NeMO identifier: {api_result['id']}." \
    f"For more information about the original data, see {url}"


    dataset_metadata["tissue_type"] = api_result["aliquot"].get("specimen_type")
    dataset_metadata["technique"] = api_result["aliquot"].get("technique")
    dataset_metadata["dataset_type"] = api_result["data_type"]  # will change based on other factors in nemoarchive_validate_metadata.cgi

    dataset_metadata["file_format"] = api_result["file_format"]   # should be already got

    # ? How to handle multiple taxa?  For now just take the first one as primary organism
    taxa = api_result["taxa"]
    first_taxon = taxa[0] if taxa else {}
    dataset_metadata["organism"] = first_taxon.get("name")

    dataset_metadata["assay"] = api_result["aliquot"].get("assay")

    # Now for the dataset metadata I cannot quite get from the API
    dataset_metadata["normalization_method"] = None
    dataset_metadata["log_transformation"] = "raw"
    # ? Can we derive this from api_result["analysis"]?
    dataset_metadata["primary_analysis_completed"] = False

    dataset_metadata["reference_annot_id"] = None
    #dataset_metadata["reference_annot_id"] = get_reference_annot_id(connection, nemo_id)

    # get GCP bucket path
    dataset_metadata["gcp_uri"] = None
    dataset_metadata["http_uri"] = None

    manifest_file_urls = api_result["manifest_file_urls"]
    if manifest_file_urls:
        # find the entry where protocol is "gcp"
        gcp_manifest = next((entry for entry in manifest_file_urls if entry["protocol"] == "gcp"), None)
        if gcp_manifest:
            dataset_metadata["gcp_uri"] = gcp_manifest["file_location"]

        # find the entry where protocol is "http"
        http_manifest = next((entry for entry in manifest_file_urls if entry["protocol"] == "http"), None)
        if http_manifest:
            dataset_metadata["http_uri"] = http_manifest["file_location"]

    return dataset_metadata

def process_nemo_assets_api_sample_result(api_result):
    """
    {
        "id": "nemo:oka-w7dgf1q",
        "subjects": [
            "string"
        ],
        "sample_attributes": [
            {
            "name": "string",
            "value": "string",
            "unit": "string",
            "cv_term_id": "string"
            }
        ],
        "source_sample_id": "string",
        "sample_source": "string",
        "sample_name": "string",
        "entity_type": "sample",
        "sample_aggregation": "string",
        "sample_subtype": "string",
        "pool_id": "string",
        "alternate_id": "string",
        "lab": "string",
        "libraries": [
            {
            "nemo_id": "nemo:wgn-nrc4w5s",
            "library_name": "string",
            "library_type": "aliquot",
            "modality": "string",
            "technique": "string",
            "assay": "string",
            "modality_cv_term_id": "string",
            "technique_cv_term_id": "string",
            "assay_cv_term_id": "string",
            "specimen_type": "string"
            }
        ],
        "anatomical_regions": [
            {
            "short_name": "string",
            "region_name": "string",
            "cv_term_id": "string"
            }
        ],
        "parent_samples": [
            "string"
        ],
        "child_samples": [
            "string"
        ],
        "files": [
            "string"
        ]
    }
    """

    def find_attribute(sample_attributes: list = [], attr_name: str = ""):
        for attr in sample_attributes:
            if attr_name in attr["name"].lower():
                return attr
        return {}

    def list_all_region_names(anatomical_regions: list = []):
        return ", ".join([region["region_name"] for region in anatomical_regions])

    # ? Eventually we will need to distinguish between file datasets and sample datasets.  Sample-based ones may have multiple entries for an attribute, such as age.
    sample_metadata = {}
    sample_metadata["tissue_ontology"] = list_all_region_names(api_result["anatomical_regions"])
    sample_metadata["treatment"] = find_attribute(api_result["sample_attributes"], "treatment")["value"]
    sample_metadata["age_value"] = find_attribute(api_result["sample_attributes"], "age")["value"]
    sample_metadata["age_unit"] = find_attribute(api_result["sample_attributes"], "age")["unit"]
    sample_metadata["sex_assigned_at_birth"] = find_attribute(api_result["sample_attributes"], "sex_assigned_at_birth")["value"]
    return sample_metadata


class SubmissionDataset(Resource):
    """Requests to deal with a single submission"""

    def get(self, submission_id, dataset_id):
        """Retrieve a particular submission dataset based on ID."""
        url_path = request.root_path + request.path

        result = {"self":url_path, "href":url_path}

        s_dataset = get_submission_dataset(dataset_id)
        if not s_dataset:
            abort(404, message=f"Submission dataset id {dataset_id} does not exist.")

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
        identifier = req.get("identifier")
        action = req.get("action", None)
        category = req.get("category", None)
        gene = req.get("gene", None)

        s_dataset = get_submission_dataset(dataset_id)
        result = {}

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
                abort(500, message=f"Cannot open connection to RabbitMQ publisher for submission dataset id {dataset_id}.")

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
                    abort(500, message=f"Could not open receive reply response from RabbitMQ consumer for submission dataset id {dataset_id}.")

                # Create the publisher
                payload = {}
                payload["dataset_id"] = dataset_id
                payload["identifier"] = identifier
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
                    abort(500, message=f"Could not publish payload to RabbitMQ consumer for submission dataset id {dataset_id}.")

                # Wait for callback to finish, then return the response
                while not task_finished:
                    pass
                print("[x] sending payload response for submission_dataset {} back to client".format(dataset_id), file=sys.stderr)
                result.update(response)
        else:
            result.update(submission_dataset_callback(dataset_id, metadata, session_id, url_path, action, category, gene))

        if not result["success"]:
            print(result.get("message", "Something went wrong."), file=sys.stderr)
            abort(500, message="Submission dataset id {} - {}".format(dataset_id, result.get("message", "Something went wrong.")))
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

        s_dataset = get_submission_dataset(dataset_id)

        # Insert submission members and create new submisison datasets if they do not exist
        try:
            save_submission_member(submission_id, s_dataset)

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
        s_dataset = get_submission_dataset(dataset_id)

        if action == "reset_steps":
            # Before initializing a new submission, we check this route to see if the current dataset exists
            # If it does and is not complete or loading, need to reset the steps so it can be resumed.
            self.reset_dataset_steps(s_dataset)
            result["status"] = get_db_status(s_dataset)
            result["success"] = True
            return result
        abort(405, message=f"No action or unknown action specified for submission dataset id {dataset_id}.")

    def get_submission_dataset(self, dataset_id):
        s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
        if not s_dataset:
            abort(404, message=f"Submission dataset id {dataset_id} does not exist.")
        return s_dataset

    def reset_dataset_steps(self, s_dataset):
        s_dataset.reset_incomplete_steps()
        s_dataset.save_change(attribute="log_message", value="")


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
                , "make_umap": "canceled"
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
                , "make_umap: "canceled"
                }
            return result

        result["dataset_id"] = s_dataset.dataset_id
        result["identifier"] = s_dataset.nemo_identifier
        result["share_id"] = s_dataset.dataset.share_id
        result["status"] = get_db_status(s_dataset)

        result["href"] = result["self"] + f"/{dataset_id}" # This can be GET-able

        result["success"] = True
        return result


