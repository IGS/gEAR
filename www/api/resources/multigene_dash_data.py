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

### Dotplot fxns

def create_dot_legend(fig, legend_col):
    """Creates the dot plot size legend. Edits figure in-place."""

    # Create a dot size legend
    steps = 5
    dot_legend=list()
    for i in range(steps):
        dot_legend += [["0", i, i*20+20, "{}%".format(i*20+20)]]
    dot_legend = pd.DataFrame(dot_legend,columns=['x','y','percent','text'])

    fig.add_scatter(
        x=dot_legend["x"]
        , y=dot_legend["y"]
        , text=dot_legend["text"]
        , hoverinfo="none"  # Do not have hover text
        , mode="markers+text"
        , marker=dict(
            color="#888"
            , size=dot_legend["percent"]
            , sizemode="area"
            )
        , showlegend=False
        , textposition="top center"
        , row=1
        , col=legend_col)

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

def create_dot_plot(df):
    """Creates a dot plot.  Returns the figure."""
    # x = group
    # y = gene
    # color = mean expression
    # size = percent of cells with gene

    # Taking a lot of influence from
    # https://github.com/interactivereport/CellDepot/blob/ec067978dc456d9262c3c59d212d90547547e61c/bin/src/plotH5ad.py#L113

    legend_col=5

    fig = make_subplots(rows=1, cols=legend_col
        , specs = [[{"colspan":legend_col-1}, None, None, None, {}]]   # "None" repeat much be legend_col - 2
        , subplot_titles=("", "Fraction of cells<br>in group (%)")
    )

    fig.add_scatter(
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
        , row=1
        , col=1)

    create_dot_legend(fig, legend_col)

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
    if not cluster_cols:
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

### Quadrant fxns

def add_gene_annotations_to_quadrant_plot(fig, gene_symbols_list) -> None:
    """Add annotations to point to each desired gene within the quadrant plot. Edits in-place."""
    genes_not_found = set()
    genes_none_none = set()
    for gene in gene_symbols_list:
        # Iterate through all the quadrant traces
        for data_idx in range(len(fig.data)):
            gene_indexes = [idx for idx in range(len(fig.data[data_idx].text))
                if fig.data[data_idx].text[idx] == gene]

            for idx in gene_indexes:

                # Do not add annotations at the zero-point of the plot, since they will overlap
                if "NONE/NONE" in fig.data[data_idx].name:
                    genes_none_none.add(gene)
                    break

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
            else:
                # gene wasn't found in filtered plot data
                genes_not_found.add(gene)
    return genes_not_found, genes_none_none


def create_quadrant_plot(df, control_val, compare1_val, compare2_val):
    """Generate a quadrant (fourway) plot.  Returns Plotly figure."""

    # Break the data up into different logfoldchange categories
    traces = [
        {
            "df":df[(df["s1_c_log2FC"] > 0) & (df["s2_c_log2FC"] > 0)]
            , "name":"UP/UP"
            , "color":"red"
        },
        {
            "df":df[(df["s1_c_log2FC"] < 0) & (df["s2_c_log2FC"] < 0)]
            , "name":"DOWN/DOWN"
            , "color":"black"
        },
        {
            "df":df[(df["s1_c_log2FC"] > 0) & (df["s2_c_log2FC"] < 0)]
            , "name":"UP/DOWN"
            , "color":"lightgreen"
        },
        {
            "df":df[(df["s1_c_log2FC"] < 0) & (df["s2_c_log2FC"] > 0)]
            , "name":"DOWN/UP"
            , "color":"orange"
        },
        {
            "df":df[(df["s1_c_log2FC"] > 0) & (df["s2_c_log2FC"] == 0)]
            , "name":"UP/NONE"
            , "color":"brown"
        },
        {
            "df":df[(df["s1_c_log2FC"] == 0) & (df["s2_c_log2FC"] > 0)]
            , "name":"NONE/UP"
            , "color":"cyan"
        },
        {
            "df":df[(df["s1_c_log2FC"] < 0) & (df["s2_c_log2FC"] == 0)]
            , "name":"DOWN/NONE"
            , "color":"green"
        },
        {
            "df":df[(df["s1_c_log2FC"] == 0) & (df["s2_c_log2FC"] < 0)]
            , "name":"NONE/DOWN"
            , "color":"purple"
        },
        {
            "df":df[(df["s1_c_log2FC"] == 0) & (df["s2_c_log2FC"] == 0)]
            , "name":"NONE/NONE"
            , "color":"grey"
        }
    ]


    # Scatter plot
    # x - query condition expression
    # y - ref condition expression
    # Log2-transform the values
    fig = go.Figure()

    for trace in traces:
        if trace["df"].empty:
            continue
        fig.add_scatter(
                x=trace["df"]["s1_c_log2FC"]
                , y=trace["df"]["s2_c_log2FC"]
                , name="{}:{}".format(trace["name"], str(len(trace["df"].index)))
                , text=trace["df"]["gene_symbol"]
                , mode="markers"
                , marker=dict(
                    color=trace["color"]
                    )
            )

    fig.update_xaxes(title="{} vs {} log2FC".format(compare1_val, control_val))
    fig.update_yaxes(title="{} vs {} log2FC".format(compare2_val, control_val))
    fig.update_layout(
        legend_title_text="log2 foldchange and number of genes in group"
        )
    return fig

def prep_quadrant_dataframe(adata, key, control_val, compare1_val, compare2_val, de_test_algo="t_test", fc_threshold=2, fdr_threshold=0.05, include_zero_fc=True, is_log10=False):
    """Prep the AnnData object to be a viable dataframe to use for making volcano plots."""

    # Create some filtered AnnData objects based on each individual comparision group
    de_filter1 = adata.obs[key].isin([compare1_val])
    selected1 = adata[de_filter1, :]
    de_filter2 = adata.obs[key].isin([compare2_val])
    selected2 = adata[de_filter2, :]
    de_filter3 = adata.obs[key].isin([control_val])
    selected3 = adata[de_filter3, :]
    # Query needs to be appended onto ref to ensure the test results are not flipped
    de_selected1 = selected3.concatenate(selected1)
    de_selected2 = selected3.concatenate(selected2)

    # Use diffxpy to compute DE statistics for each comparison
    if de_test_algo == "rank":
        de_results1 = de.test.rank_test(
            de_selected1
            , grouping=key
            , gene_names=de_selected1.var["gene_symbol"]
            , is_logged=is_log10
        )
        de_results2 = de.test.rank_test(
            de_selected2
            , grouping=key
            , gene_names=de_selected2.var["gene_symbol"]
            , is_logged=is_log10
        )
    else:
        de_results1 = de.test.t_test(
            de_selected1
            , grouping=key
            , gene_names=de_selected1.var["gene_symbol"]
            , is_logged=is_log10
        )
        de_results2 = de.test.t_test(
            de_selected2
            , grouping=key
            , gene_names=de_selected2.var["gene_symbol"]
            , is_logged=is_log10
        )

    # Cols - ['gene', 'pval', 'qval', 'log2fc', 'mean', 'zero_mean', 'zero_variance']
    df1 = de_results1.summary()
    df2 = de_results2.summary()

    # Build the data for the final dataframe
    df_data = {
        "gene_symbol" : df1["gene"].tolist()
        , "s1_c_log2FC" : df1["log2fc"]
        , "s2_c_log2FC" : df2["log2fc"]
        , "s1_c_qval" : df1["qval"]
        , "s2_c_qval" : df2["qval"]
        }

    df =  pd.DataFrame.from_dict(df_data)

    # filter dataframe by a specified log2FC threshold and qval (FDR) threshold
    # Also keep LFCs of 0, since we track those
    df["s1_c_log2FC_abs"] = df["s1_c_log2FC"].abs()
    df["s2_c_log2FC_abs"] = df["s2_c_log2FC"].abs()
    log_threshold = np.log2(fc_threshold)
    query = '((s1_c_log2FC_abs > @log_threshold) and (s1_c_qval > @fdr_threshold))' \
        'and ((s2_c_log2FC_abs > @log_threshold) and (s2_c_qval > @fdr_threshold))'
    if include_zero_fc:
        query = '((s1_c_log2FC == 0) or ((s1_c_log2FC_abs > @log_threshold) and (s1_c_qval > @fdr_threshold)))' \
        'and ((s2_c_log2FC == 0) or ((s2_c_log2FC_abs > @log_threshold) and (s2_c_qval > @fdr_threshold)))'
    df = df.query(query)

    return df

def validate_quadrant_conditions(control_condition, compare_group1, compare_group2):
    """Ensure quadrant conditions make sense."""
    if not (control_condition and compare_group1 and compare_group2):
        raise PlotError('Must pass three conditions in order to generate a volcano plot.')

    (control_key, control_val) = control_condition.split(';-;')
    (compare1_key, compare1_val) = compare_group1.split(';-;')
    (compare2_key, compare2_val) = compare_group2.split(';-;')

    if control_key != compare1_key and control_key != compare2_key:
        raise PlotError("All comparable conditions must came from same observation group.")

    return control_key, control_val, compare1_val, compare2_val

### Violin fxns

def create_stacked_violin_plot(df, gene_map, groupby_filter):
    """Create a stacked violin plot.  Returns the figure."""

    # Melt the datafram to make it easier to retrieve the contents for each axis
    df = df.melt(id_vars=[groupby_filter])
    # Create series of gene symbols by reverse-lookup of dict created for violin plot
    df["gene_symbol"] = df["index"].apply(lambda x: next((k for k, v in gene_map.items() if v == x), None))

    grouped = df.groupby([groupby_filter, "gene_symbol"])
    # Add all groupby_filter groups to a list to preserve order
    groupby_groups = []
    for group in grouped.groups.keys():
        if group[0] not in groupby_groups:
            groupby_groups.append(group[0])

    color_cycler = cycle(VIVID_COLORS)
    color_map = {cat: next(color_cycler) for cat in groupby_groups}

    # Map indexes for subplot ordering.  Indexes start at 1 since plotting rows/cols start at 1
    facet_row_indexes = {group: idx for idx, group in enumerate(groupby_groups, start=1)}

    fig = make_subplots(
        rows=len(groupby_groups)
        , cols=1
        , row_titles=groupby_groups
        , shared_xaxes=True
        , shared_yaxes="all"    # to keep the scale the same
        )

    # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
    # Group is the 'groupby' dataframe
    for name, group in grouped:
        row_idx = facet_row_indexes[name[0]]

        fig.add_violin(
            x=group["gene_symbol"]
            , y=group["value"]
            , scalegroup="_".join(name) # Name will be a tuple
            , showlegend=False
            , fillcolor=color_map[name[0]]
            , line=dict(color="slategrey")
            , points=False
            , box=dict(
                visible=False
                )
            , spanmode="hard"   # Do not extend violin tails beyond the min/max values
            , row=row_idx
            , col=1
        )

        fig.update_yaxes(
            side="right"
            , row=row_idx
            , col=1
        )

    # Color the annotations with the fill color
    fig.for_each_annotation(lambda a: a.update(font=dict(color=color_map[a.text])))

    # Row title annotations are on the right currently.  Reposition them to the left side
    fig.update_annotations(
        patch=dict(
            textangle=0
            , x=-0
            , xanchor="right"
        )
    )

    return fig

def create_violin_plot(df, gene_map, groupby_filter):
    """Creates a violin plot.  Returns the figure."""
    grouped = df.groupby([groupby_filter])
    names_in_legend = {}
    color_cycler = cycle(VIVID_COLORS)
    offsetgroup = 0

    fig = go.Figure()

    # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
    # Group is the 'groupby' dataframe
    for gene in gene_map:
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
                , scalegroup="{}_{}".format(gene, name)
                , showlegend=False
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
            gene_indexes = [idx for idx in range(len(fig.data[data_idx].text))
                if fig.data[data_idx].text[idx] == gene]

            for idx in gene_indexes:

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
        , snp=None
        , xlabel="log2FC"
        , ylabel="-log10(adjusted-P)" if use_adj_pvals else "-log10(P)"

    )

def modify_volcano_plot(fig, query, ref):
    """Adjust figure data to show up- and down-regulated data differently.  Edits figure in-place."""
    new_data = []
    sig_data = []
    # Keep non-significant data
    for data in fig.data:
        # Get rid of hover "GENE: " label.
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
                "name": "Upregulated in {}".format(ref)
                , "text":[]
                , "x":[]
                , "y":[]
                , "marker":{"color":"#636EFA"}
            }

            upregulated = {
                "name": "Upregulated in {}".format(query)
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
    )

def prep_volcano_dataframe(adata, key, query_val, ref_val, de_test_algo="ttest", is_log10=False):
    """Prep the AnnData object to be a viable dataframe to use for making volcano plots."""
    de_filter1 = adata.obs[key].isin([query_val])
    selected1 = adata[de_filter1, :]
    de_filter2 = adata.obs[key].isin([ref_val])
    selected2 = adata[de_filter2, :]
    # Query needs to be appended onto ref to ensure the test results are not flipped
    de_selected = selected2.concatenate(selected1)

    # Wanted to use de.test.two_sample(test=<>) but you cannot pass is_logged=True
    # which makes the ensuing plot inaccurate
    # TODO: Figure out how to get wald test to work (and work faster)
    if de_test_algo == "rank":
        de_results = de.test.rank_test(
            de_selected
            , grouping=key
            , gene_names=de_selected.var["gene_symbol"]
            , is_logged=is_log10
        )
    else:
        de_results = de.test.t_test(
            de_selected
            , grouping=key
            , gene_names=de_selected.var["gene_symbol"]
            , is_logged=is_log10
        )

    # Cols - ['gene', 'pval', 'qval', 'log2fc', 'mean', 'zero_mean', 'zero_variance']
    df = de_results.summary()
    df["ensm_id"] = de_selected.var.index
    df["pvals"] = df["pval"].fillna(1)      # Unexpressed genes show up as NaN
    df["pvals_adj"] = df["qval"].fillna(1)
    df["logfoldchanges"] = df["log2fc"]
    df["gene_symbol"] = df["gene"]

    return df

def validate_volcano_conditions(query_condition, ref_condition):
    """Ensure volcano conditions make sense."""
    if not (query_condition and ref_condition):
        raise PlotError('Must pass two conditions in order to generate a volcano plot.')

    (query_key, query_val) = query_condition.split(';-;')
    (ref_key, ref_val) = ref_condition.split(';-;')

    if query_key != ref_key:
        raise PlotError("Both comparable conditions must came from same observation group.")

    return query_key, query_val, ref_val

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
                message_list.append('<li>The following genes were mapped to 2 or more Ensembl IDs in this dataset, so one was chosen at random for the plot: {}</li>'.format(', '.join(dup_genes_intersection)))

            # Note to user which genes were not found in the dataset
            genes_not_present = [gene for gene in gene_symbols if gene not in found_genes]
            if genes_not_present:
                success = 3,
                message_list.append('<li>One or more genes were not found in the dataset: {}</li>'.format(', '.join(genes_not_present)))
        else:
            if plot_type in ["dotplot", "heatmap", "mg_violin"]:
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
        # Quadrant plot options
        compare_group1 = req.get("compare1_condition", None)
        compare_group2 = req.get("compare2_condition", None)
        fc_threshold = float(req.get("fold_change_cutoff", 2))
        fdr_threshold = float(req.get("fdr_cutoff", 0.05))
        include_zero_fc = req.get("include_zero_fc", True)
        # Volcano plot options
        query_condition = req.get('query_condition', None)
        ref_condition = req.get('ref_condition', None)
        de_test_algo = req.get("de_test_algo", "t-test")
        use_adj_pvals = req.get('adj_pvals', True)
        annotate_nonsignificant = req.get('annotate_nonsignificant', True)
        # Violin plot options
        stacked_violin = req.get('stacked_violin', False)
        violin_add_points = req.get('violin_add_points', False)
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

        # Success levels
        # -1 Failure
        # 1 Success
        # 2 Warning - duplicate genes found
        # 3 Warning - One or more genes could not be processed
        # NOTE: The success level in a warning can be overridden by another warning or error

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

        # These plot types filter to only the specific genes.
        # The other plot types use all genes and rather annotate the specific ones.
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
            try:
                key, query_val, ref_val = validate_volcano_conditions(query_condition, ref_condition)
                df = prep_volcano_dataframe(selected
                    , key
                    , query_val
                    , ref_val
                    , de_test_algo
                    , is_log10
                    )
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

            # Volcano plot expects specific parameter names (unless we wish to change the options)
            fig = create_volcano_plot(df
                , query_val
                , ref_val
                , use_adj_pvals
                )
            modify_volcano_plot(fig, query_val, ref_val)
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

            # 1) Flatten to long-form
            # 2) Create a gene symbol column by mapping to the Ensembl IDs
            df = df.melt(id_vars=[groupby_filter])
            ensm_to_gene = selected.var.to_dict()["gene_symbol"]
            df["gene_symbol"] = df["index"].map(ensm_to_gene)

            # Percent of all cells in this group where the gene has expression
            percent = lambda row: round(len([num for num in row if num > 0]) / len(row) * 100, 2)
            grouped = df.groupby(["gene_symbol", groupby_filter])
            df = grouped.agg(['mean', 'count', ('percent', percent)]) \
                .reset_index()

            fig = create_dot_plot(df)

        elif plot_type == "quadrant":
            try:
                key, control_val, compare1_val, compare2_val = validate_quadrant_conditions(ref_condition, compare_group1, compare_group2)
                df = prep_quadrant_dataframe(selected
                        , key
                        , control_val
                        , compare1_val
                        , compare2_val
                        , de_test_algo
                        , fc_threshold
                        , fdr_threshold
                        , include_zero_fc
                        , is_log10
                        )
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

            fig = create_quadrant_plot(df, control_val, compare1_val, compare2_val)
            # Annotate selected genes
            if gene_symbols:
                genes_not_found, genes_none_none = add_gene_annotations_to_quadrant_plot(fig, gene_symbols)
                if genes_not_found:
                    success = 3
                    message += "<li>One or more genes were did not pass cutoff filters to be in the plot: {}</li>".format(', '.join(genes_not_found))
                if genes_none_none:
                    success = 3
                    message += "<li>One or more genes had no fold change in both comparisons and will not be annotated: {}</li>".format(', '.join(genes_none_none))


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

            if stacked_violin:
                fig = create_stacked_violin_plot(df
                    , gene_map
                    , groupby_filter
                    )
            else:
                fig = create_violin_plot(df
                    , gene_map
                    , groupby_filter
                    )

            # Add jitter-based args (to make beeswarm plot)
            if violin_add_points:
                fig.update_traces(
                    jitter=0.25
                    , points="all"
                    , pointpos=0
                    , marker=dict(color="#000000")
                )
        else:
            return {
                'success': -1,
                'message': "Plot type {} is not a valid multi-gene plot option".format(plot_type)
            }

        # If figure is actualy a JSON error message, send that instead
        if "success" in fig and fig["success"] == -1:
            return fig

        # change background to pure white
        fig.update_layout(
            template="simple_white"
        )

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # NOTE: With volcano plots, the Chrome "devtools" cannot load the JSON response occasionally
        return {
            "success": success
            , "message": message
            , 'gene_symbols': gene_symbols
            , 'plot_json': json.loads(plot_json)
            , "plot_config": get_config()
        }
