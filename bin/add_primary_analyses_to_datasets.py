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
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'], 
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]


"""

import os
import shutil
import sys
import re
import json

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
VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']
VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']]
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'], 
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]


def main():
    if len(sys.argv) > 1:
        dataset_ids = [sys.argv[1]]
    else:
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

        h5ad_changes_made = False
        json_changes_made = False

        tsne_detected = detect_tsne(adata)
        if tsne_detected:
            if 'tsne' not in analysis_json or analysis_json['tsne']['tsne_calculated'] == False:
                json_changes_made = True
            
            analysis_json['tsne']['tsne_calculated'] = True
            analysis_json['tsne']['plot_tsne'] = 1
            if not has_tsne(adata):
                print("\tAdding tSNE analysis")
                add_tsne_analysis(adata)
                h5ad_changes_made = True

        umap_detected = detect_umap(adata)
        if umap_detected:
            if 'tsne' not in analysis_json or 'tsne_calculated' not in analysis_json['tsne'] or analysis_json['tsne']['tsne_calculated'] == False:
                json_changes_made = True
            
            # the parent key 'tsne' here should really be renamed to 'dimred' or something like it
            analysis_json['tsne']['umap_calculated'] = True
            analysis_json['tsne']['plot_umap'] = 1
            if not has_umap(adata):
                print("\tAdding UMAP analysis")
                add_umap_analysis(adata)
                h5ad_changes_made = True

        clustering_detected = detect_clustering(adata)
        if clustering_detected:
            if 'louvain' not in analysis_json or analysis_json['louvain']['calculated'] == False:
                json_changes_made = True
            
            analysis_json['louvain']['calculated'] = True
            if not has_louvain(adata):
                print("\tAdding clustering analysis")
                add_clustering_analysis(adata)
                h5ad_changes_made = True

        if h5ad_changes_made:
            # 482cb81f-4816-e4e3-a3fe-514707b847d8.h5ad
            print("\tWriting new H5AD")
            adata.write()

        # Does the JSON need to be updated?
        if json_changes_made:
            print("\tWriting new JSON")
            with open(analysis_json_path, 'w') as outfile:
                json.dump(analysis_json, outfile, indent=3)

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

def add_umap_analysis(adata):
    cols = adata.obs.columns.tolist()

    for pair in VALID_UMAP_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            adata.obsm['X_umap'] = adata.obs[[pair[0], pair[1]]].values
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
    Looks for the combination of pairs in VALID_TSNE_PAIRS
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            return True

    return False

def detect_umap(adata):
    """
    Looks for the combination of pairs in VALID_UMAP_PAIRS
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_UMAP_PAIRS:
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

def has_umap(adata):
    try:
        if type(adata.obsm['X_umap']):
            return True
    except:
        pass

    return False



if __name__ == '__main__':
    main()
