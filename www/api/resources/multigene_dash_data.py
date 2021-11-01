from flask import request
from flask_restful import Resource
import pandas as pd
import diffxpy.api as de

import numpy as np
import json
import os, sys
import geardb
from gear.plotting import get_config
from plotly.utils import PlotlyJSONEncoder
from itertools import cycle, groupby, product

import plotly.express as px
import plotly.graph_objects as go
import dash_bio as dashbio
from plotly.subplots import make_subplots

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
, "d0659d69-1a33-8b84-252c-f7ded46aa3d6"
, "cee5325d-434f-fefe-d2e6-e0be39421951"    # Unsure if log2 or log10 but necessary for correct output
]

DARK24_COLORS = px.colors.qualitative.Dark24  # 24 colors.  Could be problematic if more groups are chosen
ALPHABET_COLORS = px.colors.qualitative.Alphabet
LIGHT24_COLORS = px.colors.qualitative.Light24
VIVID_COLORS = px.colors.qualitative.Vivid
PALETTE_CYCLER = [DARK24_COLORS, ALPHABET_COLORS, LIGHT24_COLORS, VIVID_COLORS]

### Dotplot fxns

def create_dot_plot(df):
    """Creates a dot plot.  Returns the figure."""
    # x = group
    # y = gene
    # color = mean expression
    # size = percent of cells with gene

    # Taking a lot of influence from
    # https://github.com/interactivereport/CellDepot/blob/ec067978dc456d9262c3c59d212d90547547e61c/bin/src/plotH5ad.py#L113

    fig = make_subplots(rows=1, cols=5
        , specs = [[{"colspan":4}, None, None, None, {}]]
        , subplot_titles=(None, "Fraction of cells<br>in group (%)")
    )

    legend_col=5

    fig.add_trace(go.Scatter(
        x=df["cluster"]
        , y=df["gene_symbol"]
        , text = df["value", "count"]
        , hovertemplate="N: %{text}<br>" +
            "Percent: %{marker.size:.2f}<br>" +
            "Mean: %{marker.color:.2f}"
        , mode="markers"
        , marker=dict(
            color=df["value", "mean"]
            , colorscale=['rgb(0,0,255)','rgb(150,0,90)','rgb(255,0,0)']
            , size=df["value", "percent"]
            , sizemode="area"
            , colorbar=dict(
                title="mean"
                )
            )
        , showlegend=False
    ), row=1, col=1)

    # Create a dot size legend
    max_size = df["value", "percent"].max()
    steps = 5
    maxDotSize=20*round(max_size)/100
    dot_legend=list()
    for i in range(steps):
        dot_legend += [["0", i, i*20+20, "{}%".format(i*20+20)]]
    dot_legend = pd.DataFrame(dot_legend,columns=['x','y','percent','text'])

    fig.add_trace(go.Scatter(
        x=dot_legend["x"]
        , y=dot_legend["y"]
        , text=dot_legend["text"]
        , mode="markers+text"
        , marker=dict(
            color="#888"
            , size=dot_legend["percent"]
            , sizemode="area"
            )
        , showlegend=False
        , textposition="top center"
    ), row=1, col=legend_col)

    fig.update_layout(
        template="simple_white"    # change background to pure white
        , title_font_size=4
    )
    # Hide dot size legend axes
    fig.update_xaxes(
        visible=False
        , col=legend_col
    )
    fig.update_yaxes(
        range=[-1, 5]  # Give extra clearance so top dot does not overlap with title
        ,visible=False
        , col=legend_col
    )

    return fig


### Heatmap fxns

def create_clustergram(df, gene_symbols, is_log10=False, cluster_cols=False, flip_axes=False, groupby_filter=None, distance_metric="euclidean"):
    """Generate a clustergram (heatmap+dendrogram).  Returns Plotly figure and dendrogram trace info."""

    # Clustergram (heatmap) plot
    # https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_clustergram.py
    # https://dash.plotly.com/dash-bio/clustergram
    # gene_symbol input will be only genes shown in heatmap, and clustered

    # Gene symbols fail with all genes (dataset too large)

    df = df if flip_axes else df.transpose()
    rows = list(df.index)
    columns = list(df.columns.values)

    values = df.loc[rows].values + 1
    if is_log10:
        values = df.loc[rows].values

    hidden_labels = None
    if not groupby_filter:
        hidden_labels = "row" if flip_axes else "col"

    # Configuring which axes are clustered or not.
    cluster="all"
    col_dist = distance_metric
    row_dist = distance_metric
    if cluster_cols:
        cluster = "col" if flip_axes else "row"
        if flip_axes:
            row_dist = None
        else:
            col_dist = None

    # If just one gene, use go.heatmap instead
    if len(gene_symbols) == 1:
        return go.Figure(data=go.Heatmap(
            z=values
            , x=gene_symbols if flip_axes else columns
            , y=rows if flip_axes else gene_symbols
            , colorscale="balance"
        ))

    return dashbio.Clustergram(
        data=values
        , column_labels=gene_symbols if flip_axes else columns
        , row_labels=rows if flip_axes else gene_symbols
        , hidden_labels=hidden_labels
        , cluster=cluster
        , col_dist=col_dist
        , row_dist=row_dist
        , color_map="balance"               # Heatmap colors
        , display_ratio=0.5                 # Make dendrogram slightly bigger relative to plot
        , line_width=1                      # Make dendrogram lines thicker
        , log_transform=False if is_log10 else True
        , height=700
        , width=700
    )

def modify_clustergram(fig, flip_axes=False, gene_sym_len=1):
    """Curate the clustergram. Edits 'fig' inplace."""

    if gene_sym_len > 1:
        hyperlink_genes = []
        axis = "xaxis5" if flip_axes else "yaxis5"
        for gene in fig.layout[axis]['ticktext']:
            # Inherit from parent tag. Plotly's default style is to make the "a" tag blue which also messes with hovertext
            hyperlink_genes.append("<a style='fill:inherit;>{}</a>".format(gene))
        fig.layout[axis]['ticktext'] = hyperlink_genes

    else:
        axis = "xaxis" if flip_axes else "yaxis"

        # Make heatmap boxes square. Reference anchor should be the observations, else plot leaves empty space on the sides
        fig.layout[axis]['scaleanchor'] = 'y' if flip_axes else 'x'

        fig.layout['yaxis']['constrain'] = "domain"
        fig.layout['yaxis']['constraintoward'] = "top"
        fig.layout['xaxis']['constrain'] = "domain"
        fig.layout['xaxis']['constraintoward'] = "right"

### Violin fxns

def create_violin_plot(df, gene_map, groupby_filter):
    """Creates a violin plot.  Returns the figure."""
    fig = go.Figure()
    grouped = df.groupby([groupby_filter])
    names_in_legend = {}
    color_cycler = cycle(VIVID_COLORS)
    offsetgroup = 0

    # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
    # Group is the 'groupby' dataframe
    for gene in gene_map:
        showlegend = True
        fillcolor = next(color_cycler)
        for name, group in grouped:
            # If facets are present, a legend group trace can appear multiple times.
            # Ensure it only shows once.
            if gene in names_in_legend:
                showlegend = False
            names_in_legend[gene] = True

            fig.add_violin(
                x=group[groupby_filter]
                , y=group[gene_map[gene]]
                , name=gene
                , legendgroup=gene
                , scalegroup="{}_{}".format(gene, name)
                , showlegend=showlegend
                , fillcolor=fillcolor
                , offsetgroup=offsetgroup   # Cleans up some weird grouping stuff, making plots thicker
                , line=dict(color=fillcolor)
                , points=False
                , box=dict(
                    visible=True
                    , fillcolor='slategrey'
                    , line=dict(width=0)
                    )
                , meanline=dict(
                    color="white"
                    )
                , spanmode="hard"   # Do not extend violin tails beyond the min/max values
            )
        offsetgroup += 1

    fig.update_layout(
        # Since each gene/groupby filter is on its own trace,
        # plots are overlayed by group.  So change to "group" mode to stagger each group
        violinmode='group'
    )
    return fig

### Volcano fxns

def add_gene_annotations_to_volcano_plot(fig, gene_symbols_list, annot_nonsig=False) -> None:
    """Add annotations to point to each desired gene within the volcano plot. Edits in-place."""
    for gene in gene_symbols_list:
        # Insignificant genes are at index 0.  If you want to skip annotating them, start at index 1
        for data_idx in range(0 if annot_nonsig else 1, len(fig.data)):
            # TODO: The endswith is fragile and should probably be done differently, maybe with regex.
            gene_indexes = [idx for idx in range(len(fig.data[data_idx].text))
                if fig.data[data_idx].text[idx] == gene]

            for idx in gene_indexes:
                """
                # Determine if the arrow tail to head goes left-to-right (negative value) or the other way
                ax_offset = 0
                if fig.data[data_idx].x[idx] < 0:
                    ax_offset = -2
                elif fig.data[data_idx].x[idx] > 0:
                    ax_offset = 2
                """

                fig.add_annotation(
                        arg=dict(
                            font=dict(
                                color="white"
                            )
                        )
                        , arrowcolor="black"
                        , ax=fig.data[data_idx].x[idx] * 1.05
                        , ay=fig.data[data_idx].y[idx] + 2
                        , axref="x"
                        , ayref="y"
                        , bgcolor=fig.data[data_idx]["marker"]["color"] if data_idx > 0 else "slategrey"
                        , borderpad=2
                        , showarrow=True
                        , text=gene
                        , x=fig.data[data_idx].x[idx]
                        , y=fig.data[data_idx].y[idx]
                        , xref="x"
                        , yref="y"
                    )

def create_volcano_plot(df, query, ref, use_adj_pvals=False):
    """Generate a volcano plot.  Returns Plotly figure."""
    # Volcano plot
    # https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_volcano.py
    # https://dash.plotly.com/dash-bio/volcanoplot
    # gene_symbol input will be annotated in Volcano plot

    # df must have the following columns ["EFFECTSIZE", "GENE", "P", "SNP (optional)"]

    # y = log2 p-value
    # x = raw fold-change (effect size)

    return dashbio.VolcanoPlot(
        dataframe=df
        , title="Differences in {} vs {}".format(query, ref)
        , col="lightgrey"
        , effect_size="logfoldchanges"
        , effect_size_line_color="black"
        , gene="gene_symbol"
        , genomewideline_value= -np.log10(0.05)
        , genomewideline_color="black"
        , highlight_color="black"
        , p="pvals_adj" if use_adj_pvals else "pvals"
        , xlabel="log2FC"
        , ylabel="-log10(adjusted-P)" if use_adj_pvals else "-log10(P)"

    )

def modify_volcano_plot(fig):
    """Adjust figure data to show up- and down-regulated data differently.  Edits figure in-place."""
    new_data = []
    sig_data = []
    # Keep non-significant data
    for data in fig.data:

        # Get rid of unused SNP info, and "GENE: " label.
        data['customdata'] = [text.split(' ')[1].split('<br>')[0] for text in data['text']] # ensembl ID (as passed to SNP property)
        data['text'] = [text.split(' ')[-1] for text in data['text']]   # gene symbol

        if not (data["name"] and data["name"] == "Point(s) of interest"):
            new_data.append(data)
        else:
            sig_data.append(data)
    fig.data = new_data

    fig.data[0]["name"] = "Nonsignificant Genes"

    #Split the signifcant data into up- and down-regulated traces
    for data in sig_data:
        if data["name"] and data["name"] == "Point(s) of interest":
            downregulated = {
                "name": "Upregulated in Ref"
                , "text":[]
                , 'customdata':[]
                , "x":[]
                , "y":[]
                , "marker":{"color":"#636EFA"}
            }

            upregulated = {
                "name": "Upregulated in Query"
                , "text":[]
                , 'customdata':[]
                , "x":[]
                , "y":[]
                , "marker":{"color":"#EF553B"}
            }
            for i in range(len(data['x'])):
                if data['x'][i] > 0:
                    upregulated['x'].append(data['x'][i])
                    upregulated['y'].append(data['y'][i])
                    upregulated['text'].append(data['text'][i])
                    upregulated['customdata'].append(data['customdata'][i])

                else:
                    downregulated['x'].append(data['x'][i])
                    downregulated['y'].append(data['y'][i])
                    downregulated['text'].append(data['text'][i])
                    upregulated['customdata'].append(data['customdata'][i])


            for dataset in [upregulated, downregulated]:
                trace = go.Scattergl(
                    x=dataset['x']
                    , y=dataset['y']
                    , text=dataset['text']
                    , customdata=dataset['customdata']
                    , marker=dataset["marker"]
                    , mode="markers"
                    , name=dataset["name"]
                )
                fig.add_trace(trace)

    fig.update_layout(
        legend={
            "x":1
            ,"xanchor":"left"
            ,"y":0.5
            ,"yanchor":"middle"
            ,"bgcolor":"rgba(0,0,0,0)",   # transparent background
        }
        , title={
            "x":0.5
            ,"xref":"paper"
            ,"y":0.9
        }
        , template="simple_white"    # change background to pure white
    )

### Misc fxns

def create_dataframe_gene_mask(df, gene_symbols, plot_type):
    """Create a gene mask to filter a dataframe."""
    gene_filter = None
    success = 1
    message_list = []
    if 'gene_symbol' in df.columns:
        if gene_symbols:

            # Get list of duplicated genes for the dataset
            gene_counts_df = df['gene_symbol'].value_counts().to_frame()
            dup_genes = gene_counts_df.index[gene_counts_df['gene_symbol'] > 1].tolist()

            # Some genes may map to multiple Ensembl IDs, which can cause issues.  Create a 1-to-1 mapping by dropping dups
            uniq_df = df.drop_duplicates(subset=['gene_symbol'])
            dataset_genes = df['gene_symbol'].unique().tolist()
            normalized_genes_list, found_genes = normalize_searched_genes(dataset_genes, gene_symbols)

            # Use our list of genes to get a single Ensembl ID for each gene
            uniq_gene_filter = uniq_df['gene_symbol'].isin(normalized_genes_list)
            genes_df = uniq_df['gene_symbol'][uniq_gene_filter]

            if genes_df.empty:
                raise PlotError("Genes not found")

            # Now that our mapping is finished, create the gene filter
            gene_filter = df.index.isin(genes_df.index)

            # Note to user which genes were duplicated.
            dup_genes_intersection = intersection(dup_genes, normalized_genes_list)

            if dup_genes_intersection:
                success = 2
                message_list.append('The following genes were mapped to 2 or more Ensembl IDs in this dataset, so one was chosen at random for the plot: {}'.format(', '.join(dup_genes_intersection)))

            # Note to user which genes were not found in the dataset
            genes_not_present = [gene for gene in gene_symbols if gene not in found_genes]
            if genes_not_present:
                success = 3,
                message_list.append('One or more genes were not found in the dataset: {}'.format(', '.join(genes_not_present)))
        else:
            if plot_type in ["heatmap", "mg_violin"]:
                raise PlotError('Must pass in some genes before creating a plot of type {}'.format(plot_type))
    else:
        raise PlotError('Missing gene_symbol column in adata.var')
    message = "\n".join(message_list) if message_list else ""
    return gene_filter, success, message

def create_filtered_composite_indexes(filters, composite_indexes):
    """Create an index based on the 'comparison_composite_index' column."""
    all_vals = [v for k, v in filters.items()]  # List of lists

    # itertools.product returns a combation of every value from every list
    # Essentially  ((x,y) for x in A for y in B)
    filter_combinations = product(*all_vals)
    string_filter_combinations = [";".join(v) for v in filter_combinations]

    # This contains combinations of indexes that may not exist in the dataframe.
    # Use composite indexes from dataframe to return valid filtered indexes
    return intersection(string_filter_combinations, composite_indexes)

def get_analysis(analysis, dataset_id, session_id, analysis_owner_id):
    """Return analysis object based on various factors."""
    # If an analysis is posted we want to read from its h5ad
    if analysis:
        ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=analysis_owner_id)

        if 'type' in analysis:
            ana.type = analysis['type']
        else:
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

def intersection(lst1, lst2):
    """Intersection of two lists."""
    return list(set(lst1) & set(lst2))

def normalize_searched_genes(gene_list, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [g for cg in chosen_genes for g in gene_list if cg.lower() == g.lower()]
    found_genes = [cg for cg in chosen_genes for g in gene_list if cg.lower() == g.lower()]
    return case_insensitive_genes, found_genes

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

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)


class MultigeneDashData(Resource):
    """Resource for retrieving data from h5ad to be used to draw charts on UI.
    Parameters
    ----------
    gene_symbols: str
        Genes to search in adata.
    plot_type: str
        Plot type (heatmap, mg_violin, volcano)
    Returns
    -------
    dict
        Plot data
    """

    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        analysis = req.get('analysis', None)
        analysis_owner_id = req.get('analysis_owner_id', None)
        plot_type = req.get('plot_type')
        gene_symbols = req.get('gene_symbols', [])
        filters = req.get('obs_filters', {})    # Dict of lists
        groupby_filter = req.get('groupby_filter', None)
        # Heatmap opts
        cluster_cols = req.get('cluster_cols', False)
        flip_axes = req.get('flip_axes', False)
        distance_metric = req.get('distance_metric', "euclidean")
        # Volcano plot options
        query_condition = req.get('query_condition', None)
        ref_condition = req.get('ref_condition', None)
        de_test_algo = req.get("de_test_algo")
        use_adj_pvals = req.get('adj_pvals', False)
        annotate_nonsignificant = req.get('annotate_nonsignificant', False)
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        try:
            ana = get_analysis(analysis, dataset_id, session_id, analysis_owner_id)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # Using adata with "backed" mode does not work with volcano plot
        adata = ana.get_adata(backed=False)

        adata.obs = order_by_time_point(adata.obs)

        # get a map of all levels for each column
        columns = adata.obs.columns.tolist()

        if 'replicate' in columns:
            columns.remove('replicate')

        # Ensure datasets are not doubly log-transformed
        is_log10 = False
        if dataset_id in LOG10_TRANSFORMED_DATASETS:
            is_log10 = True

        success = 1
        message = ""

        # TODO: How to deal with a gene mapping to multiple Ensemble IDs
        try:
            # Some datasets have multiple ensemble IDs mapped to the same gene.
            # Drop dups to prevent out-of-bounds index errors downstream
            #var = adata.var.drop_duplicates(subset=['gene_symbol'])
            gene_filter, success, message = create_dataframe_gene_mask(adata.var, gene_symbols, plot_type)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # ADATA - Observations are rows, genes are columns
        selected = adata
        if plot_type in ['dotplot', 'heatmap', 'mg_violin'] and gene_filter is not None:
            selected = selected[:, gene_filter]

        # Filter dataframe on the chosen observation filters
        if filters:
            # Add new column to combine various groups into a single index
            selected.obs['comparison_composite_index'] = selected.obs[filters.keys()].apply(lambda x: ';'.join(map(str,x)), axis=1)
            selected.obs['comparison_composite_index'] = selected.obs['comparison_composite_index'].astype('category')
            unique_composite_indexes = selected.obs["comparison_composite_index"].unique()

            # Only want to keep indexes that match chosen filters
            # However if no filters were chosen, just use everything
            filtered_composite_indexes = create_filtered_composite_indexes(filters, unique_composite_indexes.tolist())
            if filtered_composite_indexes:
                condition_filter = selected.obs["comparison_composite_index"].isin(filtered_composite_indexes)
                selected = selected[condition_filter, :]

        if plot_type == "volcano":

            if not (query_condition and ref_condition):
                return {
                    'success': -1,
                    'message': 'Must pass two conditions in order to generate a volcano plot.'
                }

            (query_key, query_val) = query_condition.split(';-;')
            (ref_key, ref_val) = ref_condition.split(';-;')

            if query_key != ref_key:
                return {
                    'success': -1,
                    'message': "Both comparable conditions must came from same observation group."
                }

            de_filter1 = selected.obs[query_key].isin([query_val])
            selected1 = selected[de_filter1, :]
            de_filter2 = selected.obs[query_key].isin([ref_val])
            selected2 = selected[de_filter2, :]
            # Query needs to be appended onto ref to ensure the test results are not flipped
            de_selected = selected2.concatenate(selected1)

            # Wanted to use de.test.two_sample(test=<>) but you cannot pass is_logged=True
            # which makes the ensuing plot inaccurate
            # TODO: Figure out how to get wald test to work
            if de_test_algo == "rank":
                de_results = de.test.rank_test(
                    de_selected
                    , grouping=query_key
                    , gene_names=de_selected.var["gene_symbol"]
                    , is_logged=is_log10
                )
            else:
                de_results = de.test.t_test(
                    de_selected
                    , grouping=query_key
                    , gene_names=de_selected.var["gene_symbol"]
                    , is_logged=is_log10
                )

            # Cols - ['gene', 'pval', 'qval', 'log2fc', 'mean', 'zero_mean', 'zero_variance']
            df = de_results.summary()
            df["ensm_id"] = de_selected.var.index
            df["pvals"] = df["pval"].fillna(1)      # Unexpressed genes show up as NaN
            df["pvals_adj"] = df["qval"].fillna(1)
            df["logfoldchanges"] = df["log2fc"]
            df["SNP"] = df["ensm_id"]   # To easily get the ensembl id later
            df["gene_symbol"] = df["gene"]

            # Volcano plot expects specific parameter names (unless we wish to change the options)
            fig = create_volcano_plot(df
                , query_val
                , ref_val
                , use_adj_pvals
                )
            modify_volcano_plot(fig)
            if gene_symbols:
                add_gene_annotations_to_volcano_plot(fig, gene_symbols, annotate_nonsignificant)

        elif plot_type == "dotplot":
            df = selected.to_df()

            if not groupby_filter:
                return {
                    'success': -1,
                    'message': "'Groupby filter' option required for dot plots."
                }

            df[groupby_filter] = selected.obs[groupby_filter]

            ### Transform the dataframe to prep for the dot plot
            df = df.melt(id_vars=["cluster"])
            ensm_to_gene = selected.var.to_dict()["gene_symbol"]
            df["gene_symbol"] = df["index"].map(ensm_to_gene)

            # Percent of all cells in this group where the gene has expression
            percent = lambda row: round(len([num for num in row if num > 0]) / len(row) * 100, 2)
            grouped = df.groupby(["gene_symbol", groupby_filter])
            df = grouped.agg(['mean', 'count', ('percent', percent)]) \
                .reset_index()

            fig = create_dot_plot(df)

        elif plot_type == "quadrant":
            pass
        elif plot_type == "heatmap":
            # Filter genes and slice the adata to get a dataframe
            # with expression and its observation metadata
            df = selected.to_df()

            if groupby_filter:
                df[groupby_filter] = selected.obs[groupby_filter]
                grouped = df.groupby([groupby_filter])
                df = grouped.agg('mean') \
                    .dropna()

            fig = create_clustergram(df
                , gene_symbols
                , is_log10
                , cluster_cols
                , flip_axes
                , groupby_filter
                , distance_metric
                )

            modify_clustergram(fig, flip_axes, len(gene_symbols))

        elif plot_type == "mg_violin":
            df = selected.to_df()

            if not groupby_filter:
                return {
                    'success': -1,
                    'message': "'Groupby filter' option required for violin plots."
                }

            df[groupby_filter] = selected.obs[groupby_filter]

            # Naive approach of mapping gene to ensembl ID, in cases of one-to-many mappings
            gene_map = {}
            for gene in gene_symbols:
                gene_map[gene] = selected.var[selected.var.gene_symbol == gene].index.tolist()[0]

            fig = create_violin_plot(df
                , gene_map
                , groupby_filter
                )
        else:
            return {
                'success': -1,
                'message': "Plot type {} is not a valid multi-gene plot option".format(plot_type)
            }

        # If figure is actualy a JSON error message, send that instead
        if "success" in fig and fig["success"] == -1:
            return fig

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # NOTE: With volcano plots, the Chrome "devtools" cannot load the JSON response occasionally
        return {
            "success": success
            , "message": message
            , 'gene_symbols': gene_symbols
            , 'plot_json': json.loads(plot_json)
            , "plot_config": get_config()
        }
