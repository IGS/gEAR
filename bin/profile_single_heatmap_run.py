#!/opt/bin/python3

"""
This script is to process a single heatmap plot run in order to do some memory profiling.

Shaun Adkins
"""

from filprofiler.api import profile

import sys
import pandas as pd
import os
from pathlib import Path


lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import geardb
import gear.mg_plotting as mg
from gear.mg_plotting import PlotError

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

#dataset_ids = ["1b12dde9-1762-7564-8fbd-1b07b750505f"]
#dataset_ids=["c4f16a12-9e98-47be-4335-b8321282919e"]   # p1 kelley
dataset_ids=["efad163b-25cd-2de3-0507-a99e0f66c45f"]    # 176K cells


session_id = 'ee95e48d-c512-4083-bf05-ca9f65e2c12a'
analysis_owner_id = "662"
gene_symbols = [
    "Smpx",
    "Cd164l2",
    "Cldn14",
    "Calb2",
    "Evl",
    "Cdk5rap2",
    "Pcp4",
    "Pou4f3",
    "Gm30191",
    "Mlf1",
    "Acbd7",
    "Pvalb",
    "Lhfpl5"
]
plot_type = "heatmap"
obs_filters =  {
    "cell_type": [
        "DC1/2",
        "DC3",
        "Hensen",
        "IHC",
        "IPC",
        "IPhC",
        "IS",
        "IdC",
        "LGER1",
        "LGER2",
        "LGER3",
        "MGER",
        "OHC",
        "OPC",
        "OS",
        "Oc90",
        "eIHC",
        "eOHC"
    ],
    "louvain": [
        "DC1/2",
        "DC3",
        "Hensen",
        "IHC",
        "IPC",
        "IPhC",
        "IS",
        "IdC",
        "LGER1",
        "LGER2",
        "LGER3",
        "MGER",
        "OHC",
        "OPC",
        "OS",
        "Oc90",
        "eIHC",
        "eOHC"
    ]
}
sort_order = {}
plot_title=""
colorscale=""
reverse_colorscale=False
clusterbar_fields=[]
center_around_zero=False
cluster_obs=False
cluster_genes=False
matrixplot=True
flip_axes=False
hide_obs_labels=False
hide_gene_labels=False
distance_metric="euclidean"


# Stuff that doesn't need to be changed
user = geardb.get_user_from_session_id(session_id)
analysis = None
colorblind_mode = False
colors = {}

ONE_LEVEL_UP = 1
abs_path_git = Path(__file__).resolve().parents[ONE_LEVEL_UP] # git root
abs_path_www = abs_path_git.joinpath("www")

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def order_by_time_point(obs_df):
    """Order observations by time point column if it exists."""
    # check if time point order is intially provided in h5ad
    time_point_order = obs_df.get('time_point_order')
    if (time_point_order is not None and 'time_point' in obs_df.columns):
        sorted_df = obs_df.drop_duplicates().sort_values(by='time_point_order')
        # Safety check. Make sure time point is categorical before
        # calling .cat
        obs_df['time_point'] = pd.Categorical(obs_df['time_point'])
        col = obs_df['time_point'].cat
        obs_df['time_point'] = col.reorder_categories(
            sorted_df.time_point.drop_duplicates(), ordered=True)
        obs_df = obs_df.drop(['time_point_order'], axis=1)
    return obs_df

def get_analysis(analysis, dataset_id, session_id, analysis_owner_id):
    """Return analysis object based on various factors."""
    # If an analysis is posted we want to read from its h5ad
    if analysis:
        ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=analysis_owner_id)

        try:
            ana.type = analysis['type']
        except:
            user = geardb.get_user_from_session_id(session_id)
            ana.discover_type(current_user_id=user.id)
    else:
        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()

        # Let's not fail if the file isn't there
        if not os.path.exists(h5_path):
            raise PlotError("No h5 file found for this dataset")
        ana = geardb.Analysis(type='primary', dataset_id=dataset_id)
    return ana

def create_heatmap(dataset_id):

    """
    try:
        ana = get_analysis(analysis, dataset_id, session_id, analysis_owner_id)
    except PlotError as pe:
        return {
            'success': -1,
            'message': str(pe),
        }

    # Using adata with "backed" mode does not work with volcano plot
    adata = ana.get_adata(backed=False)
    """
    import anndata
    adata = anndata.read("/var/www/datasets/{}.h5ad".format(dataset_id))

    adata.obs = order_by_time_point(adata.obs)

    # get a map of all levels for each column
    columns = adata.obs.select_dtypes(include="category").columns.tolist()

    if 'replicate' in columns:
        columns.remove('replicate')

    if not columns:
        return {
            "success": -1,
            "message": "There are no categorical-datatype conditions found in this dataset."
        }

    # Ensure datasets are not doubly log-transformed
    # In the case of projection inputs, we don't want to log-transform either
    is_log10 = False

    success = 1
    message = ""

    # Success levels
    # -1 Failure
    # 1 Success
    # 2 Warning - duplicate genes found
    # 3 Warning - One or more genes could not be processed
    # NOTE: The success level in a warning can be overridden by another warning or error

    # TODO: How to deal with a gene mapping to multiple Ensemble IDs
    try:
        if not gene_symbols and plot_type in ["dotplot", "heatmap", "mg_violin"]:
            raise PlotError('Must pass in some genes before creating a plot of type {}'.format(plot_type))

        if len(gene_symbols) == 1 and plot_type == "heatmap":
            raise PlotError('Heatmaps require 2 or more genes as input')

        # Some datasets have multiple ensemble IDs mapped to the same gene.
        # Drop dups to prevent out-of-bounds index errors downstream
        gene_filter, success, message = mg.create_dataframe_gene_mask(adata.var, gene_symbols)
    except PlotError as pe:
        return {
            'success': -1,
            'message': str(pe),
        }

    # ADATA - Observations are rows, genes are columns
    selected = adata

    # These plot types filter to only the specific genes.
    # The other plot types use all genes and rather annotate the specific ones.
    if plot_type in ['dotplot', 'heatmap', 'mg_violin'] and gene_filter is not None:
        selected = selected[:, gene_filter]

        try:
            if plot_type == "heatmap" and len(selected.var) == 1:
                raise PlotError("Only one gene from the searched gene symbols was found in dataset.  The heatmap option require at least 2 genes to plot.")
        except PlotError as pe:
            return {
                "success": -1,
                "message": str(pe)
            }
        # Get a list of sorted ensembld IDs based on the specified gene symbol order
        ensm_to_gene = selected.var.to_dict()["gene_symbol"]
        # Since genes were restricted to one Ensembl ID we can invert keys and vals
        gene_to_ensm = {y:x for x,y in ensm_to_gene.items()}

        # Collect all genes from the unfiltered dataset
        dataset_genes = adata.var['gene_symbol'].unique().tolist()
        # Gene symbols list may have genes not in the dataset.
        normalized_genes_list, _found_genes = mg.normalize_searched_genes(dataset_genes, gene_symbols)
        # Sort ensembl IDs based on the gene symbol order
        sorted_ensm = map(lambda x: gene_to_ensm[x], normalized_genes_list)

    # Make a composite index of all categorical types
    selected.obs['composite_index'] = selected.obs[columns].apply(lambda x: ';'.join(map(str,x)), axis=1)
    selected.obs['composite_index'] = selected.obs['composite_index'].astype('category')
    columns.append("composite_index")

    ### Plot Clustergram

    # Filter genes and slice the adata to get a dataframe
    # with expression and its observation metadata
    df = selected.to_df()

    # Sort Ensembl ID columns by the gene symbol order
    df = df[sorted_ensm]

    groupby_index = "composite_index"
    groupby_fields = columns

    # Remove composite index so that the label is not duplicated.
    for cat in columns:
        df[cat] = selected.obs[cat]

    # Groupby to remove the replicates
    # Ensure the composite index is used as the index for plot labeling
    if matrixplot:
        grouped = df.groupby(groupby_fields, observed=True)
        df = grouped.agg('mean') \
            .dropna() \
            .reset_index() \
            .set_index(groupby_index)

    # Since this is the new index in the matrixplot, it does not exist as a droppable series
    # These two statements ensure that the current columns are the same if that option is set or not
    groupby_fields.remove(groupby_index)
    if not matrixplot:
        df = df.drop(columns=groupby_index)

    # Drop the obs metadata now that the dataframe is sorted
    # They cannot be in there when the clustergram is made
    # But save it to add back in later
    df_cols = pd.concat([df.pop(cat) for cat in groupby_fields], axis=1)

    global colorscale

    # Reverse Cividis so that dark is higher expression
    if colorblind_mode:
        colorscale = "cividis_r"

    # "df" must be obs label for rows and genes for cols only
    fig = mg.create_clustergram(df
        , normalized_genes_list
        , is_log10
        , cluster_obs
        , cluster_genes
        , flip_axes
        , center_around_zero
        , distance_metric
        , colorscale
        , reverse_colorscale
        , hide_obs_labels
        , hide_gene_labels
        )
    return fig

def main():
    for dataset_id in dataset_ids:
        res = create_heatmap(dataset_id)

if __name__ == "__main__":
    #main()
    profile(lambda: main(), path="/opt/gEAR/fil-results_single_heatmap")
