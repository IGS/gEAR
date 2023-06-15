#!/opt/bin/python3

# nemoarchive_validate_metadata.cgi - Write metadata to JSON file

import json
import logging
import os, subprocess, sys
from pathlib import Path

gear_root = Path(__file__).resolve().parents[2] # web-root dir
bin_path = gear_root.joinpath("bin")
sys.path.append(str(bin_path))
lib_path = gear_root.joinpath("lib")
sys.path.append(str(lib_path))

from find_best_ensembl_release_match import find_best_ensembl_release_match

import geardb
from gear.metadata import Metadata
import pandas as pd

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
    s_dataset.update_downstream_steps_to_canceled(attribute=DB_STEP)
    s_dataset.save_change(attribute="log_message", value=err_msg)
    print(json.dumps(result))

def setup_logger():
    """Set up the logger."""
    logging.basicConfig(stream=sys.stderr)
    logger = logging.getLogger('tracking_log')
    logger.setLevel(logging.INFO)
    return logger

logger = setup_logger()

def get_ensembl_release(gene_file, organism_id):
    """Given the list of genes and the organism, determine best ensemble release to use for this dataset."""

    # Couple of notes:
    # 1) Script accepts either Ensembl ID or gene symbol as column 0 and handles it appropriately
    # 2) Script treats first row as header. The odds of two releases having gene counts off by 1 is very slim
    # so I will ignore checking if the gene file has a header or not

    return find_best_ensembl_release_match(gene_file, organism_id, silent=True)

def get_genes_file_path(base_dir:Path, file_format):
    """Grab list of genes depending on file format. Returns filepath. Accepts Ensembl ID as first column too."""

    # This is working under the assumption that only one of the files in the base_dir matches
    if file_format.lower() == "mex":
        # I'm bad with glob patterns
        mex_features_file = list(base_dir.glob(r"features.tsv*"))
        mex_genes_file =  list(base_dir.glob(r"genes.tsv*"))
        genes_file = str([*mex_features_file, *mex_genes_file][0])  # One of these should match
        if genes_file.endswith(".gz"):
            import shutil, gzip
            gunzip_file = genes_file.replace(".gz", "")
            with gzip.open(genes_file, 'rb') as f_in:
                with open(gunzip_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    genes_file = gunzip_file
        return base_dir.joinpath(genes_file)
    if file_format.lower() == "tabcounts":
        genes_file =  list(base_dir.glob(r"*(genes|ROWmeta).tab"))[0]
        return base_dir.joinpath(genes_file)
    if file_format.lower() == "h5ad":
        h5ad_file = list(base_dir.glob(r"*\.h5ad"))[0]
        import anndata
        adata = anndata.read(base_dir.joinpath(h5ad_file))
        genes_file = base_dir.joinpath("h5ad_genes.tsv")
        adata.var.to_csv(genes_file, sep="\t")
        return genes_file
    raise "File format {} not supported".format(file_format)

def organism_to_taxon_id(org):
    # Returns a gear-related mapping, or None if not encountered
    # ! The taxonID will eventually come from the assets db directly
    ORG2TAX_ID = {
        "human":"9606"
        , "mouse":"10090"
        , "marmoset":"9483"
        , "zebrafish":"7955"
        , "chicken":"9031"
    }

    return ORG2TAX_ID.get(org.lower())

def get_gear_organism_id(taxon_id):
    """
    data_organism_id = {'id' : [1, 2, 3, 5, 8],
                        'label' : ['Mouse', 'Human', 'Zebrafish', 'Chicken', 'Macaque'],
                        'taxon_id' : [10090, 9606, 7955, 9031, 9544]
                        }
    """

    """ OLD WAY
    if sample_attributes.lower() in ["human", "homo sapiens","9606"]:
        return 2
    if sample_attributes.lower() in ["mouse","mus musculus", "10090"]:
        return 1
    if sample_attributes.lower() in ["zebrafish", "danio rerio", "7955"]:
        return 3
    if sample_attributes.lower() in ["chicken", "gallus gallus", "9031"]:
        return 5
    #if sample_attributes.lower() in ["macaque", "macaca mulatta", "9544"]:
    #    return 8
    logger.error("Could not associate organism or taxon id {} with a gEAR organism ID".format(sample_attributes))
    return -1
    """
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    organism_qry = ( "SELECT id FROM organism WHERE taxon_id = %s" )
    cursor.execute(organism_qry, ( str(taxon_id), ))
    for (id, ) in cursor:
        organism_id = id
    cursor.close()
    cnx.close()
    if not organism_id:
        raise
    return organism_id

def share_dataset_with_submission_user(dataset_id, user):
    """Since the dataset is owned by importer, share the dataset with the submitter."""
    #TODO: Deal with restricted datasets in the future to ensure submitter has permisison to view it
    dc = geardb.DatasetCollection()
    dc.get_shared_with_user(user=user)

    # If dataset already shared with user, no need to create a new dataset_shares entry
    if dc.datasets and any(d.id == dataset_id for d in dc.datasets):
        return

    conn = geardb.Connection()
    cursor = conn.get_cursor()

    qry = "INSERT INTO dataset_shares (dataset_id, user_id, is_allowed) VALUES (%s, %s, 1)"
    cursor.execute(qry, (dataset_id, user.id))
    cursor.close()
    conn.commit()
    conn.close()

def subspecimen_type_to_dataset_type(ss_type):
    # Returns a gear-related mapping, or None if not encountered
    # ! The mapping keys will change as the assets db API comes into effect
    SS_TYPE2_DATASET_TYPE = {
        "bulk":"bulk"
        , "cells":"single-cell-rnaseq"
        , "nuclei":"single-nucleus-rnaseq"
    }
    return SS_TYPE2_DATASET_TYPE.get(ss_type.lower())

def write_json(attributes, base_dir:Path):
    """Use supplied metadata to write a JSON file."""

    outpath = base_dir.joinpath("EXPmeta.json")
    with open(outpath, "w") as out_fh:
        json.dump(attributes, out_fh, ensure_ascii=False, indent=4)
    return outpath

def validate_metadata(dataset_id, session_id, attributes):
    user = geardb.get_user_from_session_id(session_id)
    result = {'success':False }

    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    s_dataset.save_change(attribute=DB_STEP, value="loading")

    # ! Eventually the metadata will be taken from an API call to NeMO Archive based on the identifier,
    # !     not from the neo4j database. Also may need to fail file if identifier cannot be found.
    json_attributes = {"field":[], "value":[]}  # see Metadata.read_file()

    json_attributes["field"].append("contact_name")
    json_attributes["field"].append("contact_email")
    json_attributes["field"].append("contact_institution")

    json_attributes["value"].append(attributes["dataset"]["contact_name"])
    json_attributes["value"].append(attributes["dataset"]["contact_email"])
    json_attributes["value"].append(attributes["dataset"]["contact_institute"])

    json_attributes["field"].append("title")
    #json_attributes["title"] = attributes["dataset"]["title"]

    # Temporary "title" for dataset will be the dataset name
    grant = attributes["sample"]["project_grant"]
    tissue = attributes["sample"]["tissue_ontology"]
    sex = attributes["sample"]["sex_assigned_at_birth"]
    age = attributes["sample"]["age_value"]
    age_unit = attributes["sample"]["age_unit"]
    json_attributes["value"].append(f"FAKE NAME {grant} - {tissue} - {sex} - {age}  {age_unit}")

    # Summary information
    json_attributes["field"].append("summary")
    identifier = attributes["dataset"]["identifier"]
    url = "https://assets.nemoarchive.org/{}".format(identifier)
    json_attributes["value"].append(f"This dataset was derived from nemo identifier: {identifier}." \
        f"For more information about the original data see {url}")


    # Dataset type
    json_attributes["field"].append("dataset_type")
    subspecimen_type = attributes["sample"]["sample_subspecimen_type"]
    dataset_type = subspecimen_type_to_dataset_type(subspecimen_type)
    # ATAC-Seq has it's own metadata datatype
    if "ATAC-seq".lower() in attributes["sample"]["sample_technique"].lower():
        dataset_type == "atac-seq"
    if not dataset_type:
        err_msg = "Could not find dataset type for subspecimen type {}".format(subspecimen_type)
        handle_error(s_dataset, result, err_msg)
        return
    json_attributes["value"].append(dataset_type)

    # Organism
    json_attributes["field"].append("sample_taxid")
    organism = attributes["sample"]["sample_organism"]
    taxon_id = organism_to_taxon_id(organism)
    if not taxon_id:
        err_msg = "Could not find taxon ID for organism {}".format(organism)
        handle_error(s_dataset, result, err_msg)
        return
    json_attributes["value"].append(taxon_id)

    try:
        organism_id = get_gear_organism_id(taxon_id)
    except:
        err_msg = "Could not associate taxon id {} with a gEAR organism ID".format(taxon_id)
        handle_error(s_dataset, result, err_msg)
        return

    # Not provided by NeMO for now, so we need to figure out best release ourselves
    json_attributes["field"].append("annotation_release_number")
    json_attributes["field"].append("annotation_source")

    base_dir = Path(UPLOAD_BASE_DIR).joinpath(dataset_id)
    try:
        genes_file = get_genes_file_path(base_dir, attributes["dataset"]["filetype"])
    except Exception as e:
        handle_error(s_dataset, result, str(e))
        return

    release_number = get_ensembl_release(genes_file, organism_id)
    json_attributes["value"].append(release_number)
    json_attributes["value"].append("Ensembl")  # Pretty sure these will come from Ensembl

    # Update status in dataset
    try:
        write_json(json_attributes, base_dir)
        # Save dataset to mysql
        df = pd.read_json(json.dumps(json_attributes), orient="columns")
        df.set_index('field', inplace=True)
        metadata = Metadata(metadata=df)
        if not metadata.validate():
            raise ValueError("Metadata failed validation")

        # Some fields need to change because the "Dataset" columns are different from the validation column fields.
        metadata.metadata.rename(index={"summary":"ldesc"
                        , "annotation_release_number": "annotation_release"
                        , "dataset_type": "dtype"
                        , "contact_institution": "contact_institute"
                        }, inplace=True)
        metadata.metadata["organism_id"] = organism_id
        # Drop sample_taxid, now that organism_id is present
        metadata.metadata.drop(["sample_taxid"], inplace=True)


        # Metadata class does not have an "update" command, so I am updating here (since simple dataset was created during submission)
        metadata.update_dataset_in_db(dataset_id, status="loading")
        share_dataset_with_submission_user(dataset_id, user)
        dataset = geardb.get_dataset_by_id(dataset_id)
        result["dataset"] = dataset
        s_dataset.save_change(attribute=DB_STEP, value="completed")
        result["success"] = True
    except ValueError as ve:
        handle_error(s_dataset, result, str(ve))
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        logger.error(str(ve))
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        s_dataset.update_downstream_steps_to_canceled(attribute=DB_STEP)
    except Exception as e:
        err_msg = "Dataset - {} -- Could not save status to database.".format(dataset_id)
        handle_error(s_dataset, result, err_msg)
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        logger.error(str(e))
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        s_dataset.update_downstream_steps_to_canceled(attribute=DB_STEP)
    finally:
        return result

def main():
    data = json.load(sys.stdin)
    attributes = data["attributes"]
    dataset_id = attributes["dataset"]["id"]

    # Must have a gEAR account to upload datasets
    session_id = data['session_id']


    result = validate_metadata(dataset_id, session_id, attributes)
    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))

if __name__ == '__main__':
    main()