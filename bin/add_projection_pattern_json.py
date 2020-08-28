#!/opt/bin/python3

"""

Many datasets exist in the system already with columnar elements which 
make up a user-generated analysis, such as dimensionality reduction (tSNE)
and clustering.

This script reads these and transforms them into a structure within HDF5 
which mimics what scanpy produces.  This allows them to be used throughout
the rest of the interface.

Autodetection isn't perfect.  Currently it checks for these columns for
clustering identification:

VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label']

And these pairs for tSNE:

VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined']]



"""

import os
import shutil
import sys
import re
import json

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import geardb

import numpy as np
import anndata as ad
import scanpy as sc
sc.settings.verbosity = 0

# These should match in cgi/get_embedded_tsne_display.cgi
DATASET_BASE_DIR = '/home/jorvis/git/gEAR/www/datasets'
VALID_NAME_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']
]

def main():
    dataset_ids = get_dataset_ids()

    for dataset_id in dataset_ids:
        print("Processing dataset ID: {0}".format(dataset_id))
        dataset = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5ad_path = dataset.get_file_path()

        analysis = geardb.Analysis(id=dataset_id, dataset_id=dataset_id, type='primary', vetting='owner')
        analysis_json_path = analysis.settings_path()

        if os.path.exists(analysis_json_path):
            with open(analysis_json_path) as json_in:
                analysis_json = json.load(json_in)
        else:
            with open("{0}/../data/analysis_pipeline_template.json".format(lib_path)) as json_in:
                analysis_json = json.load(json_in)

        # Analysis ID and dataset ID are the same for this type
        analysis_json["id"] = dataset_id
        analysis_json["type"] = "primary"
        analysis_json["label"] = "Primary analysis"
        analysis_json["dataset_id"] = dataset_id
        analysis_json["dataset"]["id"] = dataset_id
            
        ana = geardb.Analysis(dataset_id=dataset_id, type='primary')
        adata = ana.get_adata(backed=True)

        changes_made = False

        tsne_detected = detect_tsne(adata)
        if tsne_detected:
            analysis_json['tsne']['tsne_calculated'] = True
            analysis_json['tsne']['plot_tsne'] = 1
            if not has_tsne(adata):
                print("\tAdding tSNE analysis")
                add_tsne_analysis(adata)
                changes_made = True
        
        clustering_detected = detect_clustering(adata)
        if clustering_detected:
            analysis_json['louvain']['calculated'] = True
            if not has_louvain(adata):
                print("\tAdding clustering analysis")
                add_clustering_analysis(adata)
                changes_made = True

        if changes_made:
            try:
                adata.write()

                # we only want to write if it has both clustering and tSNE
                if tsne_detected and clustering_detected:
                    #if any([tsne_detected, clustering_detected]):
                    with open(analysis_json_path, 'w') as outfile:
                        json.dump(analysis_json, outfile, indent=3)
            except:
                print("ERROR: Unable to save write result files back out", file=sys.stderr)
            
            
def add_clustering_analysis(adata):
    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            user_defined_cluster_names = adata.obs[vname].astype('category')
            adata.obs['louvain'] = user_defined_cluster_names
            return
            
def add_tsne_analysis(adata):
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            adata.obsm['X_tsne'] = adata.obs[[pair[0], pair[1]]].values
            return
        
def get_dataset_ids():
    ids = []
    for filename in os.listdir(DATASET_BASE_DIR):
        if filename.endswith('.h5ad'):
            m = re.match("(.+?)\.h5ad", filename)
            if m:
                ids.append(m.group(1))
            else:
                print("\tmatch fail")
    return ids

def detect_clustering(adata):
    """
    Looks first for 'cluster', then 'cell_type'
    """
    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            return True

    return False

def detect_tsne(adata):
    """
    Looks for the combination of tSNE_1/tSNE_2 or tSNE1/tSNE2
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            return True

    return False

def has_louvain(adata):
    cols = adata.obs.columns.tolist()
    if 'louvain' in cols:
        return True
    else:
        return False

def has_tsne(adata):
    try:
        if type(adata.obsm['X_tsne']):
            return True    
    except:
        pass

    return False

if __name__ == '__main__':
    main()
