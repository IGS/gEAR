#!/opt/bin/python3

# write_nemoarchive_metadata_to_json.cgi - Write metadata to JSON file

import cgi
import json
import logging
import sys
from pathlib import Path

gear_root = Path(__file__).resolve().parents[2] # web-root dir
bin_path = gear_root.joinpath("bin")
sys.path.append(bin_path)
lib_path = gear_root.joinpath("lib")
sys.path.append(lib_path)

import geardb

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
UPLOAD_BASE_DIR = abs_path_www.joinpath("uploads/files")
DATASETS_DIR = abs_path_www.joinpath("datasets")

DB_STEP = "convert_metadata_status" # name of step in database

def handle_error(s_dataset, result, err_msg):
    logger.error(err_msg)
    s_dataset.save_change(attribute=DB_STEP, value="failed")
    s_dataset.update_downstream_steps_to_cancelled(attribute=DB_STEP)
    result["message"] = err_msg
    print(json.dumps(result))


def setup_logger():
    """Set up the logger."""
    logging.basicConfig(stream=sys.stderr)
    logger = logging.getLogger('tracking_log')
    logger.setLevel(logging.INFO)
    return logger

logger = setup_logger()

def find_importer_id():
    importer_id = this.servercfg["nemoanalytics_import"]["importer_id"]
    return geardb.get_user_by_id(importer_id)

def organism_to_taxon_id(org):
    # Returns a gear-related mapping, or None if not encountered
    ORG2TAX_ID = {
        "human":"9606"
        , "mouse":"10090"
        , "marmoset":"9483"
    }

    return ORG2TAX_ID.get(org.lower())

def share_dataset_with_submission_user(dataset_id, user_id):
    """Since the dataset is owned by importer, share the dataset with the submitter."""
    #TODO: Deal with restricted datasets in the future to ensure submitter has permisison to view it
    dc = geardb.DatasetCollection()
    datasets = dc.get_shared_with_user(user_id)
    our_datasets = [d for d in datasets if d.id == dataset_id]

    # If dataset already shared with user, no need to create a new dataset_shares entry
    if len(our_datasets):
        return

    conn = geardb.Connection()
    cursor = conn.get_cursor()

    qry = "INSERT INTO dataset_shares (dataset_id, user_id, is_allowed) VALUES (%s, %s, 1)"
    cursor.execute(qry, (dataset_id, user_id))
    cursor.close()
    conn.commit()
    conn.close()

def subspecimen_type_to_dataset_type(ss_type):
    # Returns a gear-related mapping, or None if not encountered
    SS_TYPE2_DATASET_TYPE = {
        "bulk":"bulk"
        , "cells":"single-cell-rnaseq"
        , "nuclei":"single-nucleus-rnaseq"
    }
    return SS_TYPE2_DATASET_TYPE.get(ss_type.lower())

def write_json(file, base_dir, filetype):
    """Use supplied metadata to write a JSON file."""

    # 3tab should already have json present
    if filetype == "3tab":
        return base_dir.glob("*EXPmeta.json")[0]

def main():
    form = cgi.FieldStorage(environ="post")
    attributes = json.loads(form.getvalue("attributes"))
    dataset_id = attributes["id"]
    filetype = attributes["filetype"]
    session_id = form.getvalue('session_id')

    # Must have a gEAR account to upload datasets
    user = geardb.get_user_from_session_id(session_id)

    base_dir = Path(UPLOAD_BASE_DIR).joinpath(dataset_id)
    result = {'success':False, 'message': '', 'json_name': None }

    print('Content-Type: application/json\n\n', flush=True)


    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    s_dataset.save_change(attribute=DB_STEP, value="loading")

    # Add some metadata to attributes if not present
    """
            'annotation_release_number',
            'annotation_source',
    """

    # "Importer" service account will be the contact
    importer = find_importer_id()
    if not importer:
        err_msg = "Could not find gEAR importer account"
        handle_error(s_dataset, result, err_msg)
        return

    # ! Eventually the metadata will be taken from an API call to NeMO Archive based on the identifier,
    # !     not from the neo4j database. Also may need to fail file if identifier cannot be found.

    attributes["contact_name"] = importer.user_name
    attributes["contact_institution"] = importer.institution    # ? should this be attributes["organization"]
    attributes["contact_email"] = importer.email

    # Temporary "title" for dataset will be the dataset name
    attributes["title"] = attributes["name"]

    # Summary information
    url = "https://assets.nemoarchive.org/{}".format(attributes["identifier"])
    attributes["summary"] = f"""
    This dataset was derived from nemo identifier: {attributes["identifier"]}.
    For more information about the original data see {url}
    """

    # Dataset type
    dataset_type = subspecimen_type_to_dataset_type(attributes["sample_subspecimen_type"])
    # ATAC-Seq has it's own metadata datatype
    if "ATAC-seq".lower() in attributes["technique"].lower():
        dataset_type == "atac-seq"
    attributes["dataset_type"] = dataset_type
    if not dataset_type:
        err_msg = "Could not find dataset type for subspecimen type {}".format(attributes["sample_subspecimen_type"])
        handle_error(s_dataset, result, err_msg)
        return

    # Organism
    taxon_id = organism_to_taxon_id(attributes["sample_organism"])
    attributes["sample_taxid"] = taxon_id
    if not taxon_id:
        err_msg = "Could not find taxon ID for organism {}".format(attributes["sample_organism"])
        handle_error(s_dataset, result, err_msg)
        return

    # Update status in dataset
    try:
        metadata_file_path = write_json(attributes, base_dir, filetype)
        # Save dataset to mysql
        #metadata = Metadata(file_path=metadata_file_path)
        #metadata.save_to_mysql(status='loading')
        share_dataset_with_submission_user(dataset_id, user.id)
        dataset = geardb.get_dataset_by_id(dataset_id)
        result["dataset"] = dataset
        s_dataset.save_change(attribute=DB_STEP, value="complete")
        result["success"] = True
        result["json_name"] = Path(metadata_file_path).name
    except Exception as e:
        handle_error(s_dataset, result, err_msg)
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        logger.error(str(e))
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        s_dataset.update_downstream_steps_to_cancelled(attribute=DB_STEP)
        result["message"] = "Submission {} - Dataset - {} -- Could not save status to database.".format("test", dataset_id)

    print(json.dumps(result))

if __name__ == '__main__':
    main()