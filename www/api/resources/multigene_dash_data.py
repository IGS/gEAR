from flask import request
from flask_restful import Resource
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 0

import numpy as np
import json
import os, sys
import geardb
from gear.plotting import get_config
from plotly.utils import PlotlyJSONEncoder
from itertools import cycle, product

import plotly.express as px
import plotly.graph_objects as go
import dash_bio as dashbio

# SAdkins - 2/15/21 - This is a list of datasets already log10-transformed where if selected will use log10 as the default dropdown option
# This is meant to be a short-term solution until more people specify their data is transformed via the metadata
LOG10_TRANSFORMED_DATASETS = [
"320ca057-0119-4f32-8397-7761ea084ed1"
, "df726e89-b7ac-d798-83bf-2bd69d7f3b52"
, "bad48d04-db27-26bc-2324-e88506f751fd"
, "dbd715bf-778a-4923-6fe7-c587987cdb00"
, "c8d99d13-394f-a87f-5d3a-395968fdb619"
, "bee735e5-d180-332c-7892-dd751dd76bb8"
, "c4f16a12-9e98-47be-4335-b8321282919e"
, "17a07bf4-b41a-d9c3-9aa7-b4729390f57a"
, "6a0a2bca-0f86-59d0-4e3d-4457be3a71ff"
, "39e01b71-415f-afa7-0c64-f0e996be0fb7"
, "6482c608-a6bd-d8b1-6bc1-5b53c34ed61c"
, "0c5a4c18-c2a9-930c-6e52-ef411f54eb67"
, "3c02d449-61ab-4bcd-f100-5f5937b1794e"
, "23e3797f-3016-8142-cbe8-69b03131ad95"
, "b16eeb8d-d68e-c7c9-9dc9-a3f4821e9192"
, "b96f448a-315d-549d-6e8a-83cdf1ce1b5c"
, "b0420910-a0fa-e920-152d-420b6275d3af"
, "f1ce4e63-3577-8020-8307-e88f1fb98953"
, "2f79f784-f7f7-7dc3-9b3e-4c87a4346d91"
, "c32835d3-cac4-bb0e-a90a-0b41dec6617a"
, "fbe1296e-572c-d388-b9d1-6e2a6bf10b0a"
, "1b12dde9-1762-7564-8fbd-1b07b750505f"
, "a2dd9f06-5223-0779-8dfc-8dce7a3897e1"
, "f7de7db2-b4cb-ebe3-7f1f-b278f46f1a7f"
, "e34fa5c6-1083-cacb-eedf-23f59f2e005f"
, "0c5fb6b0-31ab-6bfc-075d-76756ccd56b4"
, "a183b2e6-ab38-458a-52a6-5eb014d073da"
, "c4f16a12-9e98-47be-4335-b8321282919e"
, "2a25e445-2776-8913-076f-9a147a43e8b4"
, "2786d849-f11c-2de6-b22e-12c940aafe07"
, "2e3423b3-74db-d436-8357-abb3031d47e9"
, "4cb2ac62-c283-86a9-83cb-2c1b381948f2"
];

DARK24_COLORS = px.colors.qualitative.Dark24  # 24 colors.  Could be problematic if more groups are chosen
ALPHABET_COLORS = px.colors.qualitative.Alphabet
LIGHT24_COLORS = px.colors.qualitative.Light24
VIVID_COLORS = px.colors.qualitative.Vivid
PALETTE_CYCLER = [DARK24_COLORS, ALPHABET_COLORS, LIGHT24_COLORS, VIVID_COLORS]

COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"  # If colors are in observation

def add_gene_annotations_to_volcano_plot(fig, gene_symbols_list) -> None:
    """Add annotations to point to each desired gene within the volcano plot. Edits in-place."""
    # Very hacky way to add annotations. need to redo.
    for gene in gene_symbols_list:
        # Skip insignificant genes
        for data_idx in range(1, len(fig.data)):
            # TODO: The endswith is fragile and should probably be done differently, maybe with regex.
            gene_indexes = [idx for idx in range(len(fig.data[data_idx].text))
                if fig.data[data_idx].text[idx].endswith("<br>GENE: {}".format(gene))]
            for idx in gene_indexes:
                fig.add_annotation(
                        arg=dict(
                            font=dict(
                                color="white"
                            )
                        )
                        , arrowcolor="black"
                        , bgcolor=fig.data[data_idx]["marker"]["color"]
                        , borderpad=2
                        , showarrow=True
                        , text=gene
                        , x=fig.data[data_idx].x[idx]
                        , y=fig.data[data_idx].y[idx]
                        , xref="x"
                        , yref="y"
                    )

def build_column_group_markers(traces, filter_indexes):
    """Build dictionaries of group annotations for the clustergram."""

    col_group_markers = {}
    trace_column_ids = list(traces["column_ids"])
    # k = obs_category, elem = single observation, i = index position of elem in all observations
    for k, v in filter_indexes.items():
        col_group_markers.setdefault(k, [0 for i in trace_column_ids])
        for elem in v:
            # Filter indexes in order of clustering. Every index should map to a unique column ID
            for i in v[elem]:
                column_id = trace_column_ids.index(i)
                col_group_markers[k][column_id] = {'group': elem}
    return col_group_markers

def build_obs_group_indexes(df, filters):
    """Build dict of group indexes for filtered groups."""
    filter_indexes = {}
    for k, v in filters.items():
        filter_indexes.setdefault(k, {})
        for elem in v:
            obs_index = df.index[df[k] == elem ]
            # Convert dataframe index to ordinal indexes
            filter_indexes[k][elem] = [df.index.get_loc(i) for i in obs_index]
    return filter_indexes

def create_filtered_composite_index(filters):
    """Create an index based on the 'comparison_composite_index' column."""
    all_vals = [v for k, v in filters.items()]  # List of lists

    # itertools.product returns a combation of every value from every list
    # Essentially  ((x,y) for x in A for y in B)
    value_combinations = product(*all_vals)
    string_value_combinations = [";".join(v) for v in value_combinations]

    # NOTE: This returns combinations of indexes that may not exist.  Those are filtered out later
    return string_value_combinations

def create_clustergram(df, gene_symbols, is_log10=False, cluster_cols=False):
    """Generate a clustergram (heatmap+dendrogram).  Returns Plotly figure and dendrogram trace info."""

    # Clustergram (heatmap) plot
    # https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_clustergram.py
    # https://dash.plotly.com/dash-bio/clustergram
    # gene_symbol input will be only genes shown in heatmap, and clustered

    rows = list(df.index)
    columns = list(df.columns.values)

    values = df.loc[rows].values + 1
    if is_log10:
        values = df.loc[rows].values

    return dashbio.Clustergram(
        data=values
        , column_labels=columns
        , row_labels=gene_symbols
        , cluster="all" if cluster_cols else "row"
        , color_map="balance"               # Heatmap colors
        , display_ratio=0.5                 # Make dendrogram slightly bigger relative to plot
        , line_width=1                      # Make dendrogram lines thicker
        , hidden_labels="col"
        , log_transform=False if is_log10 else True
        #, optimal_leaf_order=True
        , return_computed_traces=True
    )

def create_volcano_plot(df):
    # Volcano plot
    # https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_volcano.py
    # https://dash.plotly.com/dash-bio/volcanoplot
    # gene_symbol input will be annotated in Volcano plot

    # df must have the following columns ["EFFECTSIZE", "GENE", "P", "SNP (optional)"]

    # y = log2 p-value
    # x = raw fold-change (effect size)

    return dashbio.VolcanoPlot(
        dataframe=df
        , col="lightgrey"
        , effect_size="logfoldchanges"
        , gene="gene_symbol"
        , genomewideline_value= -np.log10(0.05)
        , highlight_color="black"
        , p="pvals_adj"
        , xlabel="log2 fold-change"
        , ylabel="-log10(adj-P)"

    )

def intersection(lst1, lst2):
    """Intersection of two lists."""
    return list(set(lst1) & set(lst2))

def modify_clustergram(fig, traces, filter_indexes, is_log10=False) -> None:
    """Add column traces for each filtered group.  Edits figure in-place."""

    # Append all traces but the genes dendrogram
    new_data = []
    for data in fig.data:
        if not (data["name"] and "Row" in data["name"]):
            new_data.append(data)
    fig.data = new_data

    # Clustergram "color_list" does not seem to work. Change dendrogram line color here.
    for i in range(len(fig.data) - 1):
        fig.data[i]["marker"]["color"] = "black"

    # Delete dendrogram axis
    fig.layout.pop("xaxis4", None)
    fig.layout.pop("yaxis4", None)

    # Adjust domain of heatmap, and col clusters
    fig.layout["xaxis2"]["domain"] = [0, 0.95]
    fig.layout["xaxis5"]["domain"] = [0, 0.95]
    fig.layout["xaxis8"]["domain"] = [0, 0.95]

    # Move heatmap colorbar to left of plot
    fig.data[len(fig.data) - 1]["colorbar"]["x"] = -0.1
    fig.data[len(fig.data) - 1]["colorbar"]["xanchor"] = "left"
    fig.data[len(fig.data) - 1]["colorbar"].pop("xpad", None)
    fig.data[len(fig.data) - 1]["colorbar"]["y"] = -0.05    # Align with bottom of heatmap
    fig.data[len(fig.data) - 1]["colorbar"]["yanchor"] = "bottom"
    fig.data[len(fig.data) - 1]["colorbar"]["title"]["text"] = "Log10 Gene Expression" if is_log10 else "Log2 Gene Expression"
    fig.data[len(fig.data) - 1]["colorbar"]["title"]["side"] = "right"


    col_group_markers = build_column_group_markers(traces, filter_indexes)
    groups_and_colors = set_obs_groups_and_colors(filter_indexes)

    # Create a 2D-heatmap.  Convert the discrete groups into integers.
    # One heatmap per observation category
    col_group_labels = []

    # This is derived from the heatmap axis in the figure
    x=list(range(int(fig.layout["xaxis5"]["range"][0]+5), int(fig.layout["xaxis5"]["range"][1]+5), 10))

    # Offset the obs group colorbar to the right of the heatmap
    colorbar_x = 1.02

    # Find top of original heatmap and put "groups" heatmap tracks above.  Makes a small space b/t the genes and groups tracks
    mid_y = max(fig.layout["yaxis5"]["tickvals"]) + 12

    for key, val in col_group_markers.items():
        # TODO: col_group_markers is out of order.  Needs to be in the order of traces.column_id
        z = [[ groups_and_colors[key]["groups"].index(cgm["group"]) for cgm in val ]]

        # In order to make the colorscale a discrete one, we must map the start and stop thresholds for our normalized range
        colorscale= []
        for i in range(len(groups_and_colors[key]["colors"])):
            # Start of color thresholds
            colorscale.append(( (i)/len(groups_and_colors[key]["colors"]), groups_and_colors[key]["colors"][i] ))
            # End of color thresholds
            colorscale.append(( (i+1)/len(groups_and_colors[key]["colors"]), groups_and_colors[key]["colors"][i] ))

        trace = go.Heatmap(
            x=x
            , y=[mid_y-5, mid_y+5]
            , z=z
            , colorbar=dict(
                ticktext=[group for group in groups_and_colors[key]["groups"]]
                , tickmode="array"
                , tickvals=[idx for idx in range(len(groups_and_colors[key]["groups"]))]
                , title=key
                , x=colorbar_x
                , y=-0.05   # Align with bottom of heatmap
                , yanchor="bottom"
                )
            , colorscale=colorscale
        )
        col_group_labels.append(trace)

        # Add group label to axis tuples
        fig.layout["yaxis5"]["ticktext"] = fig.layout["yaxis5"]["ticktext"] + (key, )
        fig.layout["yaxis5"]["tickvals"] = fig.layout["yaxis5"]["tickvals"] + (mid_y, )

        colorbar_x += 0.1
        mid_y += 12 # add enough gap to space the "group" tracks

    for cgl in col_group_labels:
        fig.append_trace(cgl, 2, 2)

def modify_volcano_plot(fig):
    """Adjust figure data to show up- and down-regulated data differently.  Edits figure in-place."""
    new_data = []
    sig_data = []
    # Keep non-significant data
    for data in fig.data:
        if not (data["name"] and data["name"] == "Point(s) of interest"):
            new_data.append(data)
        else:
            sig_data.append(data)
    fig.data = new_data

    fig.data[0]["name"] = "Nonsignificant genes"

    #Split the signifcant data into up- and down-regulated traces
    for data in sig_data:
        if data["name"] and data["name"] == "Point(s) of interest":
            downregulated = {
                "name": "Downregulated Genes"
                , "text":[]
                , "x":[]
                , "y":[]
                , "marker":{"color":"#636EFA"}
            }

            upregulated = {
                "name": "Upregulated Genes"
                , "text":[]
                , "x":[]
                , "y":[]
                , "marker":{"color":"#EF553B"}
            }
            for i in range(len(data['x'])):
                if data['x'][i] > 0:
                    upregulated['x'].append(data['x'][i])
                    upregulated['y'].append(data['y'][i])
                    upregulated['text'].append(data['text'][i])
                else:
                    downregulated['x'].append(data['x'][i])
                    downregulated['y'].append(data['y'][i])
                    downregulated['text'].append(data['text'][i])

            for dataset in [upregulated, downregulated]:
                trace = go.Scattergl(
                    x=dataset['x']
                    , y=dataset['y']
                    , text=dataset['text']
                    , marker=dataset["marker"]
                    , mode="markers"
                    , name=dataset["name"]
                )
                fig.add_trace(trace)

    fig.layout["legend"]["y"] = 1.0
    fig.layout["legend"]["yanchor"] = "top"

def normalize_searched_genes(gene_list, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [g for cg in chosen_genes for g in gene_list if cg.lower() == g.lower()]
    found_genes = [cg for cg in chosen_genes for g in gene_list if cg.lower() == g.lower()]
    return case_insensitive_genes, found_genes

def set_obs_groups_and_colors(filter_indexes):
    """Create mapping of groups and colors per observation category."""
    # TODO: Use observation colors if available instead of Dark24."""
    groups_and_colors = {}
    palette_cycler = cycle(PALETTE_CYCLER)
    for k, v in filter_indexes.items():
        groups_and_colors.setdefault(k, {"groups":[], "colors":[]})
        palette = next(palette_cycler)
        color_cycler = cycle(palette)
        groups_and_colors[k]["groups"] = [elem for elem in v]
        groups_and_colors[k]["colors"] = [next(color_cycler) for elem in v]
    return groups_and_colors


class MultigeneDashData(Resource):
    """Resource for retrieving data from h5ad to be used to draw charts on UI.
    Parameters
    ----------
    gene_symbols: str
        Genes to search in adata.
    plot_type: str
        Plot type (bar, violin, scatter or line)
    Returns
    -------
    dict
        Plot data
    """

    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)
        req = request.get_json()
        analysis = req.get('analysis', None)
        analysis_owner_id = req.get('analysis_owner_id', None)
        plot_type = req.get('plot_type')
        gene_symbols = req.get('gene_symbols', [])
        filters = req.get('obs_filters', {})    # Dict of lists
        cluster_cols = req.get('cluster_cols', False)
        sort_filter = req.get('sort_filter', None)
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        # If an analysis is posted we want to read from its h5ad
        if analysis:
            ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                  session_id=session_id, user_id=analysis_owner_id)

            if 'type' in analysis:
                ana.type = analysis['type']
            else:
                ana.discover_type(current_user_id=user.id)
        else:
            ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
            h5_path = ds.get_file_path()

            # Let's not fail if the file isn't there
            if not os.path.exists(h5_path):
                return {
                    'success': -1,
                    'message': "No h5 file found for this dataset"
                }
            ana = geardb.Analysis(type='primary', dataset_id=dataset_id)

        # Using adata with "backed" mode does not work with volcano plot
        adata = ana.get_adata(backed=False)


        # check if time point order is intially provided in h5ad
        time_point_order = adata.obs.get('time_point_order')
        if (time_point_order is not None and 'time_point' in adata.obs.columns):
            sorted_df = adata.obs.drop_duplicates().sort_values(by='time_point_order')
            # Safety check. Make sure time point is categorical before
            # calling .cat
            adata.obs['time_point'] = pd.Categorical(adata.obs['time_point'])
            col = adata.obs['time_point'].cat
            adata.obs['time_point'] = col.reorder_categories(
                sorted_df.time_point.drop_duplicates(), ordered=True)
            adata.obs = adata.obs.drop(['time_point_order'], axis=1)

        # get a map of all levels for each column
        columns = adata.obs.columns.tolist()

        if 'replicate' in columns:
            columns.remove('replicate')

        #gene_symbols_list = gene_symbols.split(',')
        gene_symbols_list = gene_symbols
        success = 1
        message = ""

        if 'gene_symbol' in adata.var.columns:
            dataset_genes = adata.var.gene_symbol.unique().tolist()
            normalized_genes_list, found_genes = normalize_searched_genes(dataset_genes, gene_symbols_list)
            gene_filter = adata.var.gene_symbol.isin(normalized_genes_list)
            if not gene_filter.all():
                # Use message to show a warning
                genes_not_present = [gene for gene in gene_symbols_list if gene not in found_genes]
                success = 3,
                message = 'One or more genes were not found in the dataset: {}'.format(', '.join(genes_not_present)),
        else:
            return {
                'success': -1,
                'message': 'Missing gene_symbol column in adata.var'
            }

        # TODO: How to deal with a gene mapping to multiple Ensemble IDs
        """
        if len(df.columns) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)
            df = df.iloc[:,[0]] # Note, put the '0' in a list to return a DataFrame.  Not having in list returns DataSeries instead
        """

        # If no filters, just have a 1-to-1 mapping with index for later "groupby"
        adata.obs['comparison_composite_index'] = adata.obs.index.tolist()
        groups = "all"
        if filters:
            # Add new column to combine various groups for an eventual "groupby" argument
            composite_index = adata.obs[filters.keys()].apply(lambda x: ';'.join(map(str,x)), axis=1)
            adata.obs['comparison_composite_index'] = composite_index.tolist()
            groups = intersection(create_filtered_composite_index(filters), composite_index.unique().tolist())

        adata.obs['comparison_composite_index'] = adata.obs['comparison_composite_index'].astype('category')

        # Rank the genes in order to get p_values and log-fold changes (effort size)
        # For volcano-plots
        if plot_type == "volcano":

            # Filter composite members for one observation only, since it causes by div-by-zero errors
            # Source; https://github.com/theislab/scanpy/pull/1490
            groups = list(
                adata.obs["comparison_composite_index"]
                .value_counts()
                .loc[lambda x: x > 1]
                .index
                )
            #sc.pp.filter_cells(adata, min_genes=10)
            #sc.pp.filter_genes(adata, min_cells=1)
            sc.tl.rank_genes_groups(adata, 'comparison_composite_index', use_raw=False, groups=groups, reference="rest", n_genes=0, method="wilcoxon", copy=False, corr_method='benjamini-hochberg')#, log_transformed=False)

        # ADATA - Observations are rows, genes are columns
        # adata.uns.rank_genes_groups.<column> will be a 2D array.  Outer dimension is # genes (or n_genes). Inner dimension is per 'groupby' group
        selected = adata

        # Filter the AnnData object based on our criteria
        filtered_composite_index = adata.obs["comparison_composite_index"].unique()
        if filters:
            filtered_composite_index = groups
            condition_filter = adata.obs["comparison_composite_index"].isin(filtered_composite_index)
            selected = adata[condition_filter, :]

        # Ensure datasets are not doubly log-transformed
        is_log10 = False
        if dataset_id in LOG10_TRANSFORMED_DATASETS:
            is_log10 = True

        traces = None
        if plot_type == "volcano":

            # Stolen from https://github.com/theislab/scanpy/blob/8fe1cf9cb6309fa0e91aa5cfd9ed7580e9d5b2ad/scanpy/get/get.py#L17-L93
            # rank_genes_groups colnames
            colnames = ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']

            df = [pd.DataFrame(selected.uns['rank_genes_groups'][c])[filtered_composite_index] for c in colnames]
            df = pd.concat(df, axis=1, names=[None, 'group'], keys=colnames)
            df = df.stack(level=1).reset_index()
            df['group'] = pd.Categorical(df['group'], categories=filtered_composite_index)
            df = df.sort_values(['group', 'level_0']).drop(columns='level_0')
            df = df.join(selected.var.gene_symbol, on="names")
            df["SNP"] = df['group'].astype(str) + '-' + df["gene_symbol"].astype(str)
            df = df.reset_index(drop=True)

            # Volcano plot expects specific parameter names (unless we wish to change the options)
            fig = create_volcano_plot(df)
            modify_volcano_plot(fig)
            add_gene_annotations_to_volcano_plot(fig, gene_symbols_list)

        elif plot_type == "heatmap":
            # Filter genes and slice the adata to get a dataframe
            # with expression and its observation metadata
            df = selected.to_df()
            filter_indexes = build_obs_group_indexes(selected.obs, filters)

            # If sorting by a observation column, adjust group indexes after sorting
            if sort_filter and not cluster_cols:
                sorted_df = selected.obs.sort_values(by=[sort_filter])
                df = df.reindex(sorted_df.index.tolist())
                filter_indexes = build_obs_group_indexes(sorted_df, filters)

            df = df.transpose()
            df = df[gene_filter]

            # Sort gene_symbols_list
            # Wanted to sort dataframe using gene list index as key, but it requires Pandas >1.10.0 (we use 1.0.3)
            gene_index_filter = adata.var.index.isin(df.index)
            filtered_genes = adata.var[gene_index_filter]
            sorted_genes_list = filtered_genes.gene_symbol.tolist()
            plot_stuff = create_clustergram(df, sorted_genes_list, is_log10, cluster_cols)
            fig = plot_stuff[0] # Is tossed in favor of a modified one
            traces = plot_stuff[1]
            modify_clustergram(fig, traces, filter_indexes, is_log10)
            traces = json.dumps(traces, cls=PlotlyJSONEncoder)

        else:
            return {
                'success': -1,
                'message': "Plot type {} is not a valid multi-gene plot option".format(plot_type)
            }

        # If figure is actualy a JSON error message, send that instead
        if "success" in fig and fig["success"] == -1:
            return fig

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # NOTE: With volcano plots, the Chrome "devtools" cannot load the JSON response
        return {
            "success": success
            , "message": message
            , 'gene_symbols': gene_symbols
            , 'plot_json': json.loads(plot_json)
            , 'computed_traces':json.loads(traces) if traces else traces
            , "plot_config": get_config()
        }
