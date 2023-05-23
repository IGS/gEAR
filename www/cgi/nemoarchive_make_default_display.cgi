#!/opt/bin/python3

# nemoarchive_validate_metadata.cgi - Write metadata to JSON file

import json, cgi
import sys, shutil
from pathlib import Path
import scanpy as sc

gear_root = Path(__file__).resolve().parents[2] # web-root dir
lib_path = gear_root.joinpath("lib")
sys.path.append(str(lib_path))

import geardb

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

# Primary Filter
MIN_GENES = 200
MIN_CELLS = 3
# Mitochondria QC
MITO_PREFIX = "MT-"
FILTER_MITO_COUNT =2500
FILTER_MITO_PERCENT = 5
# Highly Variable Genes
NORM_COUNTS_PER_CELL = 1e4
MIN_MEAN = 0.0125
MAX_MEAN = 3
MIN_DISPERSION = 0.5
# tSNE
N_NEIGHBORS = 10
N_PCS = 40

def get_analysis(dataset_id, user_id):
    return geardb.Analysis(type='public', id=dataset_id, dataset_id=dataset_id, user_id=user_id, vetting="owner", label="Automated Analysis")

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_id')
    category = form.getvalue("category")
    gene = form.getvalue("gene")

    # Must have a gEAR account to upload datasets
    user = geardb.get_user_from_session_id(session_id)
    user_id = user.id

    result = {'success':False, 'message': ''}

    # Verify h5ad is in primary area
    ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
    h5_path = ds.get_file_path()
    if not Path(h5_path).exists():
        raise "No h5 file found for this dataset"

    ana = get_analysis(dataset_id, user_id)
    dest_datafile_path = ana.dataset_path()
    dest_dir = Path(dest_datafile_path).parent
    dest_dir.mkdir(exist_ok=True, parents=True)

    shutil.copy(h5_path, dest_dir)

    with open("{0}/../data/analysis_pipeline_template.json".format(lib_path)) as json_in:
        analysis_json = json.load(json_in)

    # Analysis ID and dataset ID are the same for this type
    analysis_json["id"] = dataset_id
    analysis_json["type"] = "public"
    analysis_json["label"] = "Automated Analysis"
    analysis_json["dataset_id"] = dataset_id
    analysis_json["dataset"]["id"] = dataset_id
    analysis_json["user_session_id"] = session_id

    adata = ana.get_adata()
    if 'gene_symbol' not in adata.var.columns:
        sys.exit(f"AnnData.var missing 'gene_symbol' column for dataset ${dataset_id}")

    # We are going to perform a "mini" Seurat analysis just to get something reasonable
    # Defaults from https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # primary filter
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    analysis_json["primary_filter"]["calculated"] = True
    analysis_json["primary_filter"]["filter_cells_gt_n_genes"] = MIN_GENES
    analysis_json["primary_filter"]["filter_cells_gt_n_cells"] = MIN_CELLS
    # mitochondrial genes qc
    adata.var['mt'] = adata.var_names.str.startswith(MITO_PREFIX)  # annotate the group of mitochondrial genes as 'mt'
    analysis_json["qc_by_mito"]["gene_prefix"] = MITO_PREFIX
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < FILTER_MITO_COUNT, :]
    adata = adata[adata.obs.pct_counts_mt < FILTER_MITO_PERCENT, :]
    analysis_json["qc_by_mito"]["filter_mito_count"] = FILTER_MITO_COUNT
    analysis_json["qc_by_mito"]["filter_mito_perc"] = FILTER_MITO_PERCENT
    # highly variable genes
    sc.pp.normalize_total(adata, target_sum=NORM_COUNTS_PER_CELL)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=MIN_MEAN, max_mean=MAX_MEAN, min_disp=MIN_DISPERSION)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    analysis_json["select_variable_genes"]["calculated"] = True
    analysis_json["select_variable_genes"]["norm_counts_per_cell"] = NORM_COUNTS_PER_CELL
    analysis_json["select_variable_genes"]["min_dispersion"] = MIN_DISPERSION
    analysis_json["select_variable_genes"]["min_mean"] = MIN_MEAN
    analysis_json["select_variable_genes"]["max_mean"] = MAX_MEAN
    analysis_json["select_variable_genes"]["flavor"] = "seurat"
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    analysis_json['pca']['calculated'] = True
    # tSNE
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    sc.tl.tsne(adata)
    analysis_json['tsne']['n_neighbors'] = N_NEIGHBORS
    analysis_json['tsne']['n_pcs'] = N_PCS
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