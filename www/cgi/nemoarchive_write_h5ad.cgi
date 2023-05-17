#!/opt/bin/python3

# nemoarchive_write_h5ad.cgi - Convert import data archive files into a writable h5ad. This script also validates Metadata as it writes the final h5ad
import cgi
import datetime
import json
import logging
import os,sys, subprocess, shutil
from pathlib import Path


gear_root = Path(__file__).resolve().parents[2] # web-root dir
bin_path = gear_root.joinpath("bin")
sys.path.append(str(bin_path))
lib_path = gear_root.joinpath("lib")
sys.path.append(str(lib_path))

from gear.dataarchive import DataArchive
import geardb

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
UPLOAD_BASE_DIR = abs_path_www.joinpath("uploads/files")
DATASETS_DIR_PATH = abs_path_www.joinpath("datasets")
DATASETS_DIR = str(DATASETS_DIR_PATH)

DB_STEP = "convert_to_h5ad_status"

def run_command(cmd):
    logger.info("Running command: {0}".format(cmd))
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))

def setup_logger():
    """Set up the logger."""
    logging.basicConfig(stream=sys.stderr)
    logger = logging.getLogger('tracking_log')
    logger.setLevel(logging.INFO)
    return logger

logger = setup_logger()


def ensure_ensembl_index(h5_path, organism_id):
    """
    Input: An H5AD ideally with Ensembl IDs as the index.  If instead they are gene
           symbols, this function should perform the mapping.

    Output: An updated (if necessary) H5AD file indexed on Ensembl IDs after mapping.
           Returns nothing.
    """

    add_ensembl_cmd = "python3 {0}/add_ensembl_id_to_h5ad_missing_release.py -i {1} -o {1}_new.h5ad -org {2}".format(bin_path, h5_path, organism_id)
    run_command(add_ensembl_cmd)
    shutil.move("{0}_new.h5ad".format(h5_path), h5_path)


def main():
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    filetype = form.getvalue("filetype")
    base_dir = Path(UPLOAD_BASE_DIR).joinpath(dataset_id)
    result = {'success':False, 'message': '', 'dataset': None }

    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    s_dataset.save_change(attribute=DB_STEP, value="loading")
    dataset = s_dataset.dataset

    try:
        da = DataArchive(base_dir=str(base_dir), format=filetype.lower(), dataset_id=dataset_id)
        organism_id = dataset.organism_id
        if organism_id == -1:
            raise

        logger.debug("Organism ID is {}".format(organism_id))

        h5ad_path = da.convert_to_h5ad(output_dir=DATASETS_DIR)
        if not da.ens_present:
            # If anndata.var only has genes, add ensembl IDs to index
            ensure_ensembl_index(h5ad_path, organism_id, da.ens_present)

        source_metadata_filepath = base_dir.joinpath("EXPmeta.json")
        dest_metadata_filepath = "{0}/{1}.json".format(DATASETS_DIR, dataset_id)
        shutil.copy(source_metadata_filepath, dest_metadata_filepath)
        result["success"] = True
        result['message'] = 'File successfully converted to h5ad.'
        logger.info("H5AD conversion complete!")

    except Exception as e:
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        logger.error(str(e))
        result["message"] = str(e)

    # Update status in dataset
    try:
        if result["success"]:
            dataset.save_change(attribute="load_status", value="completed")
            dataset.save_change(attribute="has_h5ad", value=1)
            s_dataset.save_change(attribute=DB_STEP, value="completed")
            #result["share_id"] = dataset.share_id
            logger.info("Database updated! Exiting")
    except Exception as e:
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        logger.error(str(e))
        result["message"] = "Submission {} - Dataset - {} -- Could not save status to database.".format("test", dataset_id)
        result["success"] = False

    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))

if __name__ == '__main__':
    main()