#!/opt/bin/python3

"""

"""

import os
import shutil
import sys
import re
import json

import mysql.connector

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

import numpy as np
import anndata as ad
import scanpy as sc
sc.settings.verbosity = 0

# These should match in cgi/get_embedded_tsne_display.cgi
gear_bin_path = os.path.dirname(os.path.realpath(__file__))
gear_root_path = os.path.dirname(gear_bin_path)
DATASET_BASE_DIR = '{}/www/datasets'.format(gear_root_path)

CMD_FILE = 'fix_cmds.sh'
LOG_FILE = 'investigate.log'
# Any files larger than this (in MB) will be skipped, unless = None
FILE_SIZE_LIMIT = 5000

def main():
    log_fh = open(LOG_FILE, 'wt')
    cmd_fh = open(CMD_FILE, 'wt')
    
    for h5_file in os.listdir(DATASET_BASE_DIR):
        if not h5_file.endswith('.h5ad'):
            continue

        m = re.match("(.+)\.h5ad", h5_file)
        dataset_id = m.group(1)
        
        print("Processing dataset ID: {0}".format(dataset_id), file=sys.stderr)
        dataset = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5ad_path = dataset.get_file_path()

        if FILE_SIZE_LIMIT:
            if os.path.getsize(h5ad_path) >= (FILE_SIZE_LIMIT * 1024 * 1024):
                print("\tSkipping, File is over the size limit")
                continue
        
        adata = sc.read(h5ad_path)

        indexed_on_ensembl_id = False
        gene_symbols_present = False
        org_by_ds = get_dataset_organism_id_map()

        if dataset_id not in org_by_ds:
            print("\tSkipping, I don't have info for this dataset in the database")
            continue

        if idx_looks_like_ensembl(adata):
            print("\tAppears to be indexed by Ensembl ID", file=sys.stderr)
            indexed_on_ensembl_id = True
        else:
            print("\tDoesn't appear to be indexed by Ensembl ID", file=sys.stderr)

        if has_gene_symbols(adata):
            gene_symbols_present = True

        if indexed_on_ensembl_id and gene_symbols_present:
            print("\tIndexed on Ensembl ID and gene symbols found", file=sys.stderr)
            continue

        if gene_symbols_present and not indexed_on_ensembl_id:
            print("\tHas gene symbols but not indexed on ensembl ID, writing to command log", file=sys.stderr)
            cmd_fh.write("./add_ensembl_id_to_h5ad_missing_release__file_based.py -i {0} -o {0}.reindexed -org {1} -idp FAKEID_\n".format(
                h5ad_path, org_by_ds[dataset_id]
            ))
        else:
            print("\tNo gene symbols found. Writing var to log", file=sys.stderr)
            log_fh.write("\n---------------------------------------------\n")
            log_fh.write("-- Dataset ID: {0}\n".format(dataset_id))
            log_fh.write("---------------------------------------------\n")

            log_fh.write("Var columns:\n")
            for col in adata.var.columns:
                log_fh.write("\t'{0}'\n".format(col))

    print("\nLog file written to: {0}".format(LOG_FILE), file=sys.stderr)
    print("Command file written to: {0}\n".format(CMD_FILE), file=sys.stderr)
        
def get_dataset_organism_id_map():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    qry = "SELECT id, organism_id FROM dataset"
    cursor.execute(qry)
    orgmap = dict()

    for (id, organism_id) in cursor:
        orgmap[id] = organism_id

    return orgmap
        
def has_gene_symbols(adata):
    if 'gene_symbol' in adata.var.columns:
        return True
    else:
        return False

def idx_looks_like_ensembl(adata):
    """
    Just checks the first entry of the index list and returns True if it is formatted like
    an Ensembl ID
    """
    m = re.match("ENS.+\d", adata.var.index[0])
    if m:
        return True
    else:
        return False
    

if __name__ == '__main__':
    main()
