#!/opt/bin/python3

# convert_and_write_h5ad.cgi - Convert import data archive files into a writable h5ad. This script also validates Metadata as it writes the final h5ad
import cgi
import json
import logging
import os,sys, subprocess, shutil
from pathlib import Path


gear_root = Path(__file__).resolve().parents[2] # web-root dir
bin_path = gear_root.joinpath("bin")
sys.path.append(bin_path)
lib_path = gear_root.joinpath("lib")
sys.path.append(lib_path)

from gear.dataarchive import DataArchive
from gear.metadata import Metadata
import geardb

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
UPLOAD_BASE_DIR = abs_path_www.joinpath("uploads/files")
DATASETS_DIR = abs_path_www.joinpath("datasets")

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

def get_organism_id(metadata_path):
    with open(metadata_path) as json_file:
        jdata = json.load(json_file)
        if 'sample_taxid' in jdata:
            return jdata['sample_taxid']
        elif 'taxon_id' in jdata:
            return jdata['taxon_id']
        raise Exception("No taxon id provided in file {}".format(metadata_path, datetime.datetime.now()))

def get_gear_organism_id(sample_attributes):
    """
    data_organism_id = {'id' : [1, 2, 3, 5, 8],
                        'label' : ['Mouse', 'Human', 'Zebrafish', 'Chicken', 'Macaque'],
                        'taxon_id' : [10090, 9606, 7955, 9031, 9544]
                        }
    """
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

def ensure_ensembl_index(h5_path, organism_id, is_en):
    """
    Input: An H5AD ideally with Ensembl IDs as the index.  If instead they are gene
           symbols, this function should perform the mapping.

    Output: An updated (if necessary) H5AD file indexed on Ensembl IDs after mapping.
           Returns nothing.
    """

    if is_en == False:
        add_ensembl_cmd = "python3 {0}/add_ensembl_id_to_h5ad_missing_release.py -i {1} -o {1}_new.h5ad -org {2}".format(bin_path, h5_path, organism_id)
        run_command(add_ensembl_cmd)
        shutil.move("{0}_new.h5ad".format(h5_path), h5_path)
    else:
        add_ensembl_cmd = "/usr/local/common/Python-3.7.2/bin/python3 {0}/find_best_ensembl_release_from_h5ad.py -i {1} -org {2}".format(bin_path, h5_path, organism_id)
        run_command(add_ensembl_cmd)


def main():
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    filetype = form.getvalue("filetype")
    base_dir = Path(UPLOAD_BASE_DIR).joinpath(dataset_id)
    result = {'success':False, 'message': '', 'dataset': None }

    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    s_dataset.save_change(attribute="convert_to_h5ad_status", value="loading")

    try:
        da = DataArchive(base_dir=str(base_dir), format=filetype.lower(), dataset_id=dataset_id)
        # Validate metadata against gEAR's validator
        metadata = Metadata(file_path=metadata_file_path)
        if metadata.validate():
            logger.info("Metadata file is valid: {0}".format(metadata_file_path))
            try:
                metadata_json_path = "{0}/{1}.json".format(args.output_base, dataset_id)
                metadata.write_json(file_path=metadata_json_path)
                organism_taxa = get_organism_id(metadata_file_path)
                # Ensure organism_taxa is string in case Int is passed through JSON
                organism_id = get_gear_organism_id(str(organism_taxa))
                if organism_id == -1:
                    raise

                logger.debug("Organism ID is {}".format(organism_id))
                # If dataset directory has h5ad file, skip that step
                file_list = os.listdir(dataset_dir)
                h5_path = None
                is_en = False   # assume ENSEMBL IDs are not present if h5ad was already passed to us
                for f in file_list:
                    if f.endswith(".h5ad"):
                        h5_path = "{}/{}".format(dataset_dir, f)
                if not h5_path:
                    h5_path, is_en = convert_to_h5ad(dataset_dir, dataset_id, args.output_base)

                ensure_ensembl_index(h5_path, organism_id, is_en)
            except:
                logger.error("Failed to process file:{0}. Error is below.".format(file_path))
                exctype, value = sys.exc_info()[:2]
                logger.error("{} - {}".format(exctype, value))
                logger.info(file_path, extra={"dataset_id":dataset_id, "status":"FAILED"})
        else:
            logger.error("Metadata file is NOT valid: {0}".format(metadata_file_path))

        h5ad_path = da.convert_to_h5ad(output_dir=DATASETS_DIR)
        result["success"] = True
        result['message'] = 'File successfully read.'
        logger.info("Complete! Exiting.")

    except Exception as e:
        s_dataset.save_change(attribute="convert_to_h5ad_status", value="failed")
        logger.error(str(e))
        result["message"] = str(e)

    # Update status in dataset
    try:
        if result["success"]:
            dataset = geardb.get_dataset_by_id(dataset_id)
            dataset.save_change(attribute="load_status", value="completed")
            result["dataset"] = dataset
            s_dataset.save_change(attribute="convert_to_h5ad_status", value="complete")
    except:
        s_dataset.save_change(attribute="convert_to_h5ad_status", value="failed")
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        logger.error(str(e))
        result["message"] = "Submission {} - Dataset - {} -- Could not save status to database.".format("test", dataset_id)
        result["success"] = False

    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))

if __name__ == '__main__':
    main()