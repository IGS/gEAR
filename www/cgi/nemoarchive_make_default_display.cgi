#!/opt/bin/python3

# nemoarchive_validate_metadata.cgi - Write metadata to JSON file

import json
import sys
from pathlib import Path
import scanpy as sc

gear_root = Path(__file__).resolve().parents[2] # web-root dir
lib_path = gear_root.joinpath("lib")
sys.path.append(str(lib_path))

import geardb
from gear.metadata import Metadata
import pandas as pd

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

def get_analysis(dataset_id):
    ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
    h5_path = ds.get_file_path()

    # Let's not fail if the file isn't there
    if not Path(h5_path).exists():
        raise "No h5 file found for this dataset"
    return geardb.Analysis(type='public', id=dataset_id, dataset_id=dataset_id, vetting="owner", label="Automated Analysis")

def main():
    data = json.load(sys.stdin)
    dataset_id = data["dataset_id"]
    category = data["category"]
    gene = data["gene"]

    submission_id = data["submission_id"]

    # Must have a gEAR account to upload datasets
    session_id = data['session_id']
    user = geardb.get_user_from_session_id(session_id)
    user_id = user.id

    result = {'success':False, 'message': ''}

    ana = get_analysis(dataset_id)
    dest_datafile_path = ana.dataset_path()
    dest_dir = Path(dest_datafile_path).parent
    dest_dir.mkdir(exist_ok=True, parents=True)

    with open("{0}/../data/analysis_pipeline_template.json".format(lib_path)) as json_in:
        analysis_json = json.load(json_in)

    # Analysis ID and dataset ID are the same for this type
    analysis_json["id"] = dataset_id
    analysis_json["type"] = "public"
    analysis_json["label"] = "Automated Analysis"
    analysis_json["dataset_id"] = dataset_id
    analysis_json["dataset"]["id"] = dataset_id

    adata = ana.get_adata(backed=True)
    if 'gene_symbol' not in adata.var.columns:
        sys.exit(f"AnnData.var missing 'gene_symbol' column for dataset ${dataset_id}")

    # ? should i duplicate file
    adata.filename = dest_datafile_path.replace(".h5ad", ".orig.h5ad")

    # We are going to perform a "mini" Seurat analysis just to get something reasonable
    # Defaults from https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

    # primary filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    # mitochondrial genes qc
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    # highly variable genes
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    analysis_json['pca']['pca_calculated'] = True
    # tSNE
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.tSNE(adata)
    analysis_json['tsne']['tsne_calculated'] = True


    try:
        # If gene was not passed, or does not exist, we find a highly variable gene candidate
        if not gene:
            raise

        gene_symbols = (gene,)
        gene_filter = adata.var.gene_symbol.isin(gene_symbols)
        if not gene_filter.any():
            raise
    except:
        highly_variable_genes = adata.var[adata.var.highly_variable].sort_values('dispersions_norm', ascending=False).gene_symbol
        gene = highly_variable_genes[0]


    # Write out updated analysis AnnData file and JSON configuration
    adata.write(dest_datafile_path)
    analysis_json_path = ana.settings_path()
    with open(analysis_json_path, 'w') as outfile:
        json.dump(analysis_json, outfile, indent=3)

    # Set up tSNE plot
    plot_config = {
        "gene_symbol":gene
        ,"colors":{}
        ,"x_axis":"tSNE_1"  # ? Should this be "X_tsne_1"... this is bypassed when analysis is passed in anyways
        ,"y_axis":"tSNE_2"
        ,"order":{}
        ,"colorize_legend_by":category if category else None
        ,"plot_by_group":None
        ,"max_columns":None
        ,"skip_gene_plot":False
        ,"horizontal_legend":False
        }

    result["plot_config"] = plot_config
    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))

if __name__ == '__main__':
    main()