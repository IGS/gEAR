
from itertools import cycle, product

import dash_bio as dashbio
import diffxpy.api as de
import numpy as np
import pandas as pd
import anndata as ad
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ALPHABET_COLORS = px.colors.qualitative.Alphabet
BOLD_COLORS = px.colors.qualitative.Bold
D3_COLORS = px.colors.qualitative.D3
DARK24_COLORS = px.colors.qualitative.Dark24  # 24 colors.  Could be problematic if more groups are chosen
LIGHT24_COLORS = px.colors.qualitative.Light24
SAFE_COLORS = px.colors.qualitative.Safe
VIVID_COLORS = px.colors.qualitative.Vivid

color_swatch_map = {
    "alphabet": ALPHABET_COLORS
    , "bold": BOLD_COLORS
    , "d3": D3_COLORS
    , "dark24": DARK24_COLORS
    , "light24": LIGHT24_COLORS
    , "safe": SAFE_COLORS
    , "vivid": VIVID_COLORS
}

PALETTE_CYCLER = [DARK24_COLORS, ALPHABET_COLORS, LIGHT24_COLORS, VIVID_COLORS]

# Fractional count to add to all values so log can be computed on non-expressed (0) values
LOG_COUNT_ADJUSTER = 1

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
        range=[-0.5, 5]  # Give extra clearance so top dot does not overlap with title
        , visible=False
        , col=legend_col
    )

def create_dot_plot(df, groupby_filters, is_log10=False, plot_title=None, colorscale="Bluered", reverse_colorscale=False):
    """Creates a dot plot.  Returns the figure."""
    # x = group
    # y = gene
    # color = mean expression
    # size = percent of cells with gene

    # Taking a lot of influence from
    # https://github.com/interactivereport/CellDepot/blob/ec067978dc456d9262c3c59d212d90547547e61c/bin/src/plotH5ad.py#L113

    # Specify the subplot grid
    legend_col=5
    spec_row = [{"colspan":legend_col-1}]
    spec_row.extend([None for i in range(legend_col-2)])
    spec_row.append({})

    fig = make_subplots(rows=1, cols=legend_col
        , specs = [spec_row]   # "None" repeat much be legend_col - 2
        , subplot_titles=(plot_title, "Fraction of cells<br>in group (%)")
    )

    multicategory = create_multicategory_axis_labels(groupby_filters, df)

    # log-transform dataset if it came in raw
    mean = np.log2(df['value', 'mean'] + LOG_COUNT_ADJUSTER)
    if is_log10:
        mean = df['value', 'mean']

    if not colorscale:
        colorscale="Bluered"

    fig.add_scatter(
        x=multicategory
        , y=df["gene_symbol"]
        , text = df["value", "count"]
        , hovertemplate="N: %{text}<br>" +
            "Percent: %{marker.size:.2f}<br>" +
            "Mean: %{marker.color:.2f}"
        , mode="markers"
        , marker=dict(
            color=mean
            , colorscale=colorscale
            , reversescale=reverse_colorscale
            , size=df["value", "percent"]
            , sizemode="area"
            , colorbar=dict(
                title="Log10 Mean Expression" if is_log10 else "Log2 Mean Expression"
                )
            )
        , showlegend=False
        , row=1
        , col=1)

    x_title = groupby_filters[0]
    if len(groupby_filters) > 1:
        x_title += " and {}".format(groupby_filters[1])
    fig.update_xaxes(
        title=x_title.capitalize()
    )
    fig.update_yaxes(
        title="Genes"
    )

    create_dot_legend(fig, legend_col)

    return fig

### Heatmap fxns
#from memory_profiler import profile
#fp=open('/tmp/memory_profiler.log','w+')

def add_clustergram_cluster_bars(fig, clusterbar_indexes, obs_labels=None, is_log10=False, flip_axes=False) -> None:
    """Add column traces for each filtered group.  Edits figure in-place."""

    # Heatmap is located on xaxis5 and yaxis5 in dash-0.6.1
    # Moved to xaxis11 and yaxis11 in dash-1.0.2
    # Left-side dendrogram is xaxis4 and yaxis4 in dash=0.6.1
    # Moved to xaxis9 and yaxis9 in dash=1.0.2
    # Top-sie dendrogram is xaxis2 and yaxis2 in dash=0.6.1
    # Moved to xaxis3 and yaxis3 in dash=1.0.2"

    obs_axis = "xaxis11"
    gene_axis = "yaxis11"
    obs_dendro_axis = "xaxis9"
    gene_dendro_axis = "yaxis9"
    if flip_axes:
        obs_axis = "yaxis11"
        gene_axis = "xaxis11"
        obs_dendro_axis = "yaxis3"
        gene_dendro_axis = "xaxis3"

    title_text = "Log2 Gene Expression"
    colorbar_title = "Log2 Expr."
    if is_log10:
        title_text = "Log10 Gene Expression"
        colorbar_title = "Log10 Expr."

    fig.update_layout(
        title={
            "text":title_text
            ,"x":0.5
            ,"xref":"paper"
            ,"y":0.9
        }
    )

    # Get list of observations in the order they appear on the heatmap
    obs_order = fig.layout[obs_axis]["ticktext"]

    # Update the text with labels based on included/excluded filters
    if obs_labels:
        fig.layout[obs_axis]["ticktext"] = obs_labels

    # Assign observations to their categorical groups and assign colors to the groups
    col_group_markers = build_column_group_markers(clusterbar_indexes, obs_order)
    groups_and_colors = set_obs_groups_and_colors(clusterbar_indexes)

    # Create a 2D-heatmap.  Convert the discrete groups into integers.
    # One heatmap per observation category
    col_group_traces = []

    # Get the position of observations. Individual tickvals are the midpoint for the heatmap square
    obs_positions=fig.layout[obs_axis]["tickvals"]

    # Offset the expresison colorbar to the right of the heatmap
    # At this point it's the last trace added
    # NOTE: The bar may need to be adjusted on the gene search display page
    fig.data[-1]["colorbar"]["x"] = 1
    fig.data[-1]["colorbar"]["xpad"] = 30
    fig.data[-1]["colorbar"]["len"] = 0.5   # Set expression colorbar to be half as tall
    fig.data[-1]["colorbar"]["y"] = 0.5
    fig.data[-1]["colorbar"]["yanchor"] = "middle"
    fig.data[-1]["name"] = "expression" # Name colorbar for easier retrieval
    fig.data[-1]["colorbar"]["title"] = colorbar_title

    # Put "groups" heatmap tracks either above or to the right of the genes in heatmap
    # Makes a small space b/t the genes and groups tracks
    next_bar_position = max(fig.layout[gene_axis]["tickvals"]) + 7

    # Information that will help posittion the "groups" colorbars
    num_colorbars = len(col_group_markers.keys())
    max_heatmap_domain = max(fig.layout["yaxis11"]["domain"])
    curr_colorbar_y = 0

    for key, val in col_group_markers.items():
        z = create_clusterbar_z_value(flip_axes, groups_and_colors, key, val)

        # In order to make the colorscale a discrete one, we must map the start and stop thresholds for our normalized range
        colorscale = []
        # Normally tickvals will stretch to the min and max of the colorbar range. We also need to center the groups in the middle of the color
        tickvals = []
        for i in range(len(groups_and_colors[key]["colors"])):
            # Start of color thresholds
            colorscale.append(( (i)/len(groups_and_colors[key]["colors"]), groups_and_colors[key]["colors"][i] ))
            # End of color thresholds
            colorscale.append(( (i+1)/len(groups_and_colors[key]["colors"]), groups_and_colors[key]["colors"][i] ))
            # Center the group name in its color
            tickvals.append((colorscale[-1][0] + colorscale[-2][0]) / 2 * (len(groups_and_colors[key]["groups"])-1))

        next_x = obs_positions
        next_y = [next_bar_position-2, next_bar_position+2]
        if flip_axes:
            next_x = [next_bar_position-2, next_bar_position+2]
            next_y = obs_positions

        trace = go.Heatmap(
            x=next_x
            , y=next_y
            , z=z
            , colorbar=dict(
                ticktext=[group for group in groups_and_colors[key]["groups"]]
                #, tickvals=[idx for idx in range(len(groups_and_colors[key]["groups"]))]
                , tickvals=tickvals
                , title=key
                , x=1
                , xpad=100  # spaced far enough from the expression bar.  Needs to be adjusted for gene display panels.
                , y=curr_colorbar_y   # Align with bottom of heatmap
                , yanchor="bottom"
                , len=0.9/num_colorbars
                )
            , colorscale=colorscale
            , name="clusterbar"
        )
        col_group_traces.append(trace)

        # Add group label to axis tuples
        fig.layout[gene_axis]["ticktext"] = fig.layout[gene_axis]["ticktext"] + (key, )
        fig.layout[gene_axis]["tickvals"] = fig.layout[gene_axis]["tickvals"] + (next_bar_position, )

        next_bar_position += 4 # add enough gap to space the "group" tracks
        curr_colorbar_y += (max_heatmap_domain * 1/num_colorbars)

    # Shift genes dendropgram to account for new cluster cols
    fig.layout[gene_dendro_axis]["range"] = (min(fig.layout[gene_axis]["tickvals"]), max(fig.layout[gene_axis]["tickvals"]))

    # Discovered "xaxis" range won't autoupdate based on tickvals, so pop it off if it is there
    fig.layout[gene_axis].pop("range", None)

    for cgl in col_group_traces:
        # This was 2, 2 in Dash 0.6.1 but is now 3, 3 in Dash 1.0.2
        fig.append_trace(cgl, 3, 3)

def create_clusterbar_z_value(flip_axes, groups_and_colors, key, val):
    """Create a 2D-heatmap to represent the clusterbar.  Convert the discrete groups into integers."""
    # number of elements in z array needs to equal number of observations in x-axis
    # If axes are flipped, we need one-element arrays equal to number of observations in y-axis
    if flip_axes:
        return [[ groups_and_colors[key]["groups"].index(cgm["group"])] for cgm in val ]
    return [[ groups_and_colors[key]["groups"].index(cgm["group"]) for cgm in val ]]

#@profile(stream=fp)
def create_clustergram(df, gene_symbols, is_log10=False, cluster_obs=False, cluster_genes=False, flip_axes=False, center_around_zero=False
                       , distance_metric="euclidean", colorscale=None, reverse_colorscale=False, hide_obs_labels=False, hide_gene_labels=False):
    """Generate a clustergram (heatmap+dendrogram).  Returns Plotly figure and dendrogram trace info."""

    # Clustergram (heatmap) plot
    # https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_clustergram.py
    # https://dash.plotly.com/dash-bio/clustergram
    # gene_symbol input will be only genes shown in heatmap, and clustered

    # Gene symbols fail with all genes (dataset too large)

    df = df if flip_axes else df.transpose()
    columns = list(df.columns.values)
    rows = list(df.index)
    col_labels = gene_symbols if flip_axes else columns
    row_labels = rows if flip_axes else gene_symbols

    values = df.loc[rows].values + LOG_COUNT_ADJUSTER
    if is_log10:
        values = df.loc[rows].values

    hidden_labels = set_hidden_labels(hide_obs_labels, hide_gene_labels, flip_axes)

    cluster, col_dist, row_dist = determine_axes_cluster_information(cluster_obs, cluster_genes, flip_axes, distance_metric)

    if not colorscale:
        colorscale = "RdYlBu" if center_around_zero else "Reds"
        # In the default colorscheme, reversing the scheme depends on if the plot centers around zero.
        reverse_colorscale = center_around_zero

    fig = dashbio.Clustergram(
        data=values
        , column_labels=col_labels
        , row_labels=row_labels
        , hidden_labels=hidden_labels
        , cluster=cluster
        , col_dist=col_dist
        , row_dist=row_dist
        , center_values=center_around_zero
        , color_map=colorscale
        , display_ratio=0.3                 # Make dendrogram slightly bigger relative to plot
        , line_width=1                      # Make dendrogram lines thicker
        , log_transform=False if is_log10 else True
    )

    fig.data[-1]["reversescale"] = reverse_colorscale

    if center_around_zero:
        fig.data[-1]["zmid"] = 0
    else:
        # both zmin and zmax are required
        fig.data[-1]["zmin"] = 0
        fig.data[-1]["zmax"] = max(map(max, fig.data[-1]["z"])) # Highest z-value in 2D array

    return fig

def determine_axes_cluster_information(cluster_obs=False, cluster_genes=False, flip_axes=False, distance_metric="euclidean"):
    """Configuring which axes are clustered or not.
    Args:
        cluster_obs (bool, optional): Should observations have a dendrogram. Defaults to False.
        cluster_genes (bool, optional): Should genes have a dendrogram. Defaults to False.
        flip_axes (bool, optional): Should axes be transposed. Defaults to False.
        distance_metric (str, optional): Which distance metric to use. Defaults to "euclidean".

    Returns:
        tuple: (cluster, col_dist, row_dist)
        cluster: (str) If row, column or both are clustered.
        col_dist: (str) Distance metric, or None if not clustered.
        row_dist: (str) Distance metric, or None if not clustered.
    """
    cluster = None
    col_dist = None
    row_dist = None
    if cluster_obs and cluster_genes:
        cluster = "all"
        col_dist = distance_metric
        row_dist = distance_metric
    elif cluster_obs:
        cluster = "row" if flip_axes else "col"
        if flip_axes:
            row_dist = distance_metric
        else:
            col_dist = distance_metric
    elif cluster_genes:
        cluster = "col" if flip_axes else "row"
        if flip_axes:
            col_dist = distance_metric
        else:
            row_dist = distance_metric
    return cluster, col_dist, row_dist

def create_clustergram_observation_labels(df, fig, colname="composite_index", flip_axes=False):
    """Create a set of labels to replace the current ticktext in the clustergram."""
    obs_axis = "yaxis11" if flip_axes else "xaxis11"
    obs_order = fig.layout[obs_axis]["ticktext"]    # Gets order of composite index observations
    return df.reindex(obs_order)[colname].tolist()  # reindex based on the observation order, and return the labels

def set_hidden_labels(hide_obs, hide_genes, flip_axes):
    hidden_labels = []
    if hide_obs:
        hidden_labels.append("row" if flip_axes else "col")
    if hide_genes:
        hidden_labels.append("col" if flip_axes else "row")
    if not len(hidden_labels):
        hidden_labels = None
    return hidden_labels

### Quadrant fxns

def add_gene_annotations_to_quadrant_plot(fig, gene_symbols_list) -> None:
    """Add annotations to point to each desired gene within the quadrant plot. Edits in-place."""
    genes_not_found = set()
    genes_none_none = set()

    for gene in gene_symbols_list:
        gene_found = False
        # Iterate through all the quadrant traces
        for data_idx in range(len(fig.data)):
            gene_indexes = [idx for idx in range(len(fig.data[data_idx].text))
                if fig.data[data_idx].text[idx] == gene]

            for idx in gene_indexes:
                gene_found = True

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

        if not gene_found:
            # gene wasn't found in filtered plot data
            genes_not_found.add(gene)
    return genes_not_found, genes_none_none


def create_quadrant_plot(df, control_val, compare1_val, compare2_val, colorscale=None):
    """Generate a quadrant (fourway) plot.  Returns Plotly figure."""

    # Default colors
    colors = ["red", "black", "lightgreen", "orange", "brown", "cyan", "green", "purple"]

    # If scale is sequential, split into equal colors equal to the number of categories
    # If scale is discrete, use the colorscale
    if colorscale:
        if colorscale.lower() in px.colors.named_colorscales():
            px.colors.sample_colorscale(px.colors.get_colorscale(colorscale), len(colors))
        elif colorscale not in color_swatch_map:
            # Not all the quantitivate colorscales available are in the color_swatch_map
            raise Exception("Colorscale {} not a valid colorscale to choose from".format(colorscale))
        else:
            colors = color_swatch_map[colorscale]

    # Break the data up into different logfoldchange categories
    traces = [
        {
            "df":df[(df["s1_c_log2FC"] > 0) & (df["s2_c_log2FC"] > 0)]
            , "name":"UP/UP"
            , "color":colors[0]
        },
        {
            "df":df[(df["s1_c_log2FC"] < 0) & (df["s2_c_log2FC"] < 0)]
            , "name":"DOWN/DOWN"
            , "color":colors[1]
        },
        {
            "df":df[(df["s1_c_log2FC"] > 0) & (df["s2_c_log2FC"] < 0)]
            , "name":"UP/DOWN"
            , "color":colors[2]
        },
        {
            "df":df[(df["s1_c_log2FC"] < 0) & (df["s2_c_log2FC"] > 0)]
            , "name":"DOWN/UP"
            , "color":colors[3]
        },
        {
            "df":df[(df["s1_c_log2FC"] > 0) & (df["s2_c_log2FC"] == 0)]
            , "name":"UP/NONE"
            , "color":colors[4]
        },
        {
            "df":df[(df["s1_c_log2FC"] == 0) & (df["s2_c_log2FC"] > 0)]
            , "name":"NONE/UP"
            , "color":colors[5]
        },
        {
            "df":df[(df["s1_c_log2FC"] < 0) & (df["s2_c_log2FC"] == 0)]
            , "name":"DOWN/NONE"
            , "color":colors[6]
        },
        {
            "df":df[(df["s1_c_log2FC"] == 0) & (df["s2_c_log2FC"] < 0)]
            , "name":"NONE/DOWN"
            , "color":colors[7]
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
                , customdata=trace["df"]["ensm_id"] # Add ensembl ids
                , text=trace["df"]["gene_symbol"]
                , mode="markers"
                , marker=dict(
                    color=trace["color"]
                    )
            )

    fig.update_xaxes(title="{} vs {} log2FC".format(compare1_val, control_val))
    fig.update_yaxes(title="{} vs {} log2FC".format(compare2_val, control_val))
    fig.update_layout(
        legend_title_text="Log2FC: Num Genes in Group"
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
    de_selected1 = ad.concat([selected3, selected1], merge="same")
    de_selected2 = ad.concat([selected3, selected2], merge="same")

    if not is_log10:
        de_selected1.X = de_selected1.X + LOG_COUNT_ADJUSTER
        de_selected2.X = de_selected2.X + LOG_COUNT_ADJUSTER

    # Use diffxpy to compute DE statistics for each comparison
    de_test_func = de.test.t_test
    if de_test_algo == "rank":
        de_test_func = de.test.rank_test

    de_results1 = de_test_func(
        de_selected1
        , grouping=key
        , gene_names=de_selected1.var["gene_symbol"]
        , is_logged=is_log10
    )
    de_results2 = de_test_func(
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
        ,"ensm_id" : de_selected1.var.index
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

def validate_quadrant_conditions(obs_df, control_condition, compare_group1, compare_group2):
    """Ensure quadrant conditions make sense."""
    if not (control_condition and compare_group1 and compare_group2):
        raise PlotError('Must pass three conditions in order to generate a volcano plot.')

    (control_key, control_val) = control_condition.split(';-;')
    (compare1_key, compare1_val) = compare_group1.split(';-;')
    (compare2_key, compare2_val) = compare_group2.split(';-;')

    if control_key != compare1_key and control_key != compare2_key:
        raise PlotError("All comparable conditions must came from same observation group.")

    if control_key not in obs_df.columns:
        raise PlotError(f"Control condition {control_key} not found in observation dataframe. Please update curation.")

    return control_key, control_val, compare1_val, compare2_val

### Violin fxns

def build_violin_x_title(groupby_filters):
    """Create a title for the x axis of the violin plot."""
    x_title = "Genes"
    if groupby_filters[0]:
        x_title += " grouped by {}".format(groupby_filters[0])
    if len(groupby_filters) > 1:
        x_title += " and {}".format(groupby_filters[1])
    return x_title

def create_stacked_violin_plot(df, groupby_filters, is_log10=False, colorscale=None, reverse_colorscale=False):
    """Create a stacked violin plot.  Returns the figure."""

    # Preserve sort order passed to plot, and assign colors to primary category groups
    primary_groups = df[groupby_filters[0]].unique().tolist()
    # Map indexes for subplot ordering.  Indexes start at 1 since plotting rows/cols start at 1
    facet_row_indexes = create_facet_indexes(primary_groups)
    num_rows = len(facet_row_indexes)
    row_titles = primary_groups

    secondary_groups = []
    facet_col_indexes = {}
    num_cols = 1
    col_titles = None
    if len(groupby_filters) > 1:
        secondary_groups = df[groupby_filters[1]].unique().tolist()
        facet_col_indexes = create_facet_indexes(secondary_groups)
        num_cols = len(facet_col_indexes)
        col_titles = secondary_groups

    try:
        colors = get_discrete_colors(primary_groups, colorscale, reverse_colorscale)
    except Exception as e:
        raise PlotError("Error creating colors for violin plot: {}".format(e))
    color_cycler = cycle(colors)

    color_map = {cat: next(color_cycler) for cat in primary_groups}

    fig = make_subplots(
        rows=num_rows
        , cols=num_cols
        , row_titles=row_titles
        , column_titles=col_titles
        , shared_yaxes="all"    # to keep the scale the same for all row facets
        , x_title="Genes"
        , y_title="Log10 Expression" if is_log10 else "Log2 Expression"
        )

    groupby = ["gene_symbol"]
    groupby.extend(groupby_filters)
    # Create groupings for traces
    grouped = df.groupby(groupby)

    # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
    # Group is the 'groupby' dataframe
    for name, group in grouped:
        # name[0] is gene_sym, name[1] is primary category, name[2] is secondary category
        row_idx = facet_row_indexes[name[1]]
        col_idx = facet_col_indexes[name[2]] if len(name) > 2 else 1

        # log-transform dataset if it came in raw
        if not is_log10:
            group['value'] = np.log2(group['value'] + LOG_COUNT_ADJUSTER)

        fig.add_violin(
            x=group["gene_symbol"]
            , y=group["value"]
            , scalegroup=",".join(name) # Name will be a tuple
            , showlegend=False
            , fillcolor=color_map[name[1]]
            , line=dict(color="slategrey")
            , points=False
            , box=dict(
                visible=False
                )
            , spanmode="hard"   # Do not extend violin tails beyond the min/max values
            , row=row_idx
            , col=col_idx
        )

        # Want annotation text to left but axis tick labels on right
        fig.update_yaxes(
            side="right"
            , showticklabels=True
            , row=row_idx
            , col=col_idx
        )

    update_stacked_violin_annotations(fig, primary_groups, color_map)

    plot_title = groupby_filters[0]
    if len(groupby_filters) > 1:
        plot_title += " and {}".format(groupby_filters[1])
    # Thin out the gap between violins. Default is 0.3 for both values.
    fig.update_layout(
        violingap=0.3
        , violingroupgap=0
        , margin={
            "l":200 if len(max(primary_groups, key=len)) > 5 else 80
        }
        , title={
            "text":plot_title.capitalize()
            ,"x":0.5
            ,"xref":"paper"
            ,"y":0.95
        }
    )

    return fig

def create_violin_plot(df, groupby_filters, is_log10=False, colorscale=None, reverse_colorscale=False):
    """Creates a violin plot.  Returns the figure."""

    try:
        colors = get_discrete_colors(df["gene_symbol"].unique().tolist(), colorscale, reverse_colorscale)
    except Exception as e:
        raise PlotError("Error creating colors for violin plot: {}".format(e))
    color_cycler = cycle(colors)

    fig = go.Figure()

    genes_to_color = { gene: next(color_cycler) for gene in df["gene_symbol"].unique().tolist() }
    names_in_legend = {}

    groupby = ["gene_symbol"]
    groupby.extend(groupby_filters)
    # Create groupings for traces
    grouped = df.groupby(groupby, observed=True)

    # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
    # Group is the 'groupby' dataframe
    for name, group in grouped:
        gene_sym = name[0]
        fillcolor = genes_to_color[gene_sym]

        # If facets are present, a legend group trace can appear multiple times.
        # Ensure it only shows once.
        showlegend = False if gene_sym in names_in_legend else True
        names_in_legend[gene_sym] = True

        multicategory = create_multicategory_axis_labels(groupby_filters, group)

        # log-transform dataset if it came in raw
        if not is_log10:
            group['value'] = np.log2(group['value'] + LOG_COUNT_ADJUSTER)

        fig.add_violin(
            x=multicategory
            , y=group["value"]
            , name=gene_sym
            # Scalegroup must be unique for all violins to have max equal width.
            # Hence why we have to do a 3-dimensional groupby to make unique traces
            , scalegroup="{}".format(','.join(name))
            , fillcolor=fillcolor
            , offsetgroup=gene_sym   # Cleans up some weird grouping stuff, making plots thicker
            , showlegend=showlegend
            , line=dict(color="slategrey")
            , points=False
            , box=dict(
                visible=False
                )
            , spanmode="hard"   # Do not extend violin tails beyond the min/max values
        )

    title_text = "Log2 Gene Expression"
    y_title = "Log2 Expression"
    if is_log10:
        title_text = "Log10 Gene Expression"
        y_title = "Log10 Expression"

    fig.update_layout(
        # Since each gene/groupby filter is on its own trace,
        # plots are overlayed by group.  So change to "group" mode to stagger each group
        violinmode='group'
        , title={
            "text":title_text
            ,"x":0.5
            ,"xref":"paper"
            ,"y":0.9
        }
    )
    x_title = build_violin_x_title(groupby_filters)
    fig.update_xaxes(
        title=x_title
    )
    fig.update_yaxes(
        title=y_title
    )
    return fig

def update_stacked_violin_annotations(fig, primary_groups, color_map):
    """Adjust the annotations on the stacked violin plot. Edits Plotly figure.layout in-place"""

    fig.for_each_annotation(
        # Color the row annotations with the fill color
        # Also, row title annotations are on the right currently.  Reposition them to the left side
        # Am attempting to do this based on the assumption that row facet titles will never have yanchor of bottom
        # (or y-pos of 1) or have certain text shared with the axes titles
        lambda a: a.update(
            font=dict(color=color_map[a.text])
            , textangle=0
            , x=0
            , xanchor="right"
            , borderpad=5   # Unsure if this does anything but it should ensure the row titles don't come too close to the edge
        )
        , selector=lambda a: not (a.yanchor == "bottom" or a.text == "Genes" or a.text.endswith("Expression"))
    )

    fig.for_each_annotation(
        # Edit y-axis title
        lambda a: a.update(
            xshift=-170 if len(max(primary_groups, key=len)) > 5 else -40    # Varies based on len of row_facet group names
        )
        , selector=lambda a: a.text.endswith("Expression")
    )

    fig.for_each_annotation(
        # Edit x-axis title
        lambda a: a.update(
            yshift=-50
        )
        , selector=lambda a: a.text == "Genes"
    )

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

def categorize_volcano_datapoint(nonsig_data, sig_data, data):
    """Categorize volcano datapoints based on whether they are significant or not."""
    if not (data["name"] and data["name"] == "Point(s) of interest"):
        nonsig_data.append(data)
    else:
        sig_data.append(data)

def create_volcano_plot(df, query, ref, pval_threshold, logfc_bounds, use_adj_pvals=False):
    """Generate a volcano plot.  Returns Plotly figure."""
    # Volcano plot
    # https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_volcano.py
    # https://dash.plotly.com/dash-bio/volcanoplot
    # gene_symbol input will be annotated in Volcano plot

    # df must have the following columns ["EFFECTSIZE", "GENE", "P", "SNP (optional)"]

    # y = log2 p-value
    # x = raw fold-change (effect size)

    # NOTE: We cannot pass "customdata" to the function, so we must modify the trace afterwards.

    return dashbio.VolcanoPlot(
        dataframe=df
        , title="Differences in {} vs {}".format(query, ref)
        , col="lightgrey"
        , effect_size="logfoldchanges"
        , effect_size_line=logfc_bounds
        , effect_size_line_color="black"
        , gene="ensm_id"    # Will change to 'gene_symbol' in modification step later
        , genomewideline_value= -np.log10(pval_threshold)
        , genomewideline_color="black"
        , highlight_color="black"
        , p="pvals_adj" if use_adj_pvals else "pvals"
        , snp=None
        , xlabel="log2FC"
        , ylabel="-log10(adjusted-P)" if use_adj_pvals else "-log10(P)"
    )

def curate_volcano_datapoint_text(data):
    """Format the text of volcano plot datapoints.  Edits in-place."""
    # Get rid of hover "GENE: " label.
    data['text'] = [text.split(' ')[-1] for text in data['text']]   # gene symbol

def modify_volcano_plot(fig, query, ref, ensm2genesymbol, downcolor=None, upcolor=None):
    """Adjust figure data to show up- and down-regulated data differently.  Edits figure in-place."""
    nonsig_data = []
    sig_data = []
    for data in fig.data:
        curate_volcano_datapoint_text(data)
        categorize_volcano_datapoint(nonsig_data, sig_data, data)

    # Non-significant data does not need to be modified. It is one trace.
    fig.data = nonsig_data

    fig.data[0]["name"] = "Nonsignificant Genes"
    fig.data[0]["customdata"] = fig.data[0]["text"] # Ensembl ID to "customdata" and gene symbols to "text" properties
    gene_symbol_list = [ensm2genesymbol[t] for t in fig.data[0]["text"]]
    fig.data[0]["text"] = gene_symbol_list

    downcolor = downcolor or "#636EFA"
    upcolor = upcolor or "#EF553B"

    #Split the signifcant data into up- and down-regulated traces
    for data in sig_data:
        if data["name"] and data["name"] == "Point(s) of interest":
            downregulated = {
                "name": "Upregulated in {}".format(ref)
                , "text":[] # gene_symbol
                , "customdata":[]   # ensembl id
                , "x":[]
                , "y":[]
                , "marker":{"color":downcolor}
            }

            upregulated = {
                "name": "Upregulated in {}".format(query)
                , "text":[]
                , "customdata":[]
                , "x":[]
                , "y":[]
                , "marker":{"color":upcolor}
            }
            for i in range(len(data['x'])):
                if data['x'][i] > 0:
                    upregulated['x'].append(data['x'][i])
                    upregulated['y'].append(data['y'][i])
                    upregulated['text'].append(ensm2genesymbol[data['text'][i]])
                    upregulated['customdata'].append(data['text'][i])

                else:
                    downregulated['x'].append(data['x'][i])
                    downregulated['y'].append(data['y'][i])
                    downregulated['text'].append(ensm2genesymbol[data['text'][i]])
                    downregulated['customdata'].append(data['text'][i])

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
            ,"y":1
        }
    )

def prep_volcano_dataframe(adata, key, query_val, ref_val, de_test_algo="ttest", is_log10=False):
    """Prep the AnnData object to be a viable dataframe to use for making volcano plots."""
    de_filter1 = adata.obs[key].isin([query_val])
    selected1 = adata[de_filter1, :]
    if ref_val == "rest":
        # Rest is union of rest of conditions that are not the query condition
        selected2 = adata[~de_filter1, :]
    else:
        de_filter2 = adata.obs[key].isin([ref_val])
        selected2 = adata[de_filter2, :]
    # Query needs to be appended onto ref to ensure the test results are not flipped
    de_selected = ad.concat([selected2, selected1], merge="same")

    if not is_log10:
        de_selected.X = de_selected.X + LOG_COUNT_ADJUSTER

    # Wanted to use de.test.two_sample(test=<>) but you cannot pass is_logged=True
    # which makes the ensuing plot inaccurate
    # TODO: Figure out how to get wald test to work (and work faster)
    de_test_func = de.test.t_test
    if de_test_algo == "rank":
        de_test_func = de.test.rank_test

    de_results = de_test_func(
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

def validate_volcano_conditions(obs_df, query_condition, ref_condition):
    """Ensure volcano conditions make sense."""
    if not (query_condition and ref_condition):
        raise PlotError('Must pass two conditions in order to generate a volcano plot.')

    (query_key, query_val) = query_condition.split(';-;')
    (ref_key, ref_val) = ref_condition.split(';-;')

    # Shortening the name for ease
    if ref_val == "Union of the rest of the groups":
        ref_val = "rest"

    if query_key != ref_key:
        raise PlotError("Both comparable conditions must came from same observation group.")

    if query_key not in obs_df.columns:
        raise PlotError("Observation key {} not found in dataset. Please update curation.".format(query_key))

    return query_key, query_val, ref_val

### Misc fxns

def build_column_group_markers(filter_indexes, obs_order):
    """Build dictionaries of group annotations for the clustergram."""
    col_group_markers = {}
    # k = obs_category, elem = single observation, i = index position of indiv. observation in observation order list
    for k, v in filter_indexes.items():
        col_group_markers.setdefault(k, [0 for obs in obs_order])
        for elem in v:
            # Filter indexes in order of clustering
            for obs in set(v[elem]):
                for column_id, value in enumerate(obs_order):
                    if obs == value:
                        col_group_markers[k][column_id] = {'group': elem}
    return col_group_markers

def build_obs_group_indexes(df, filters, clusterbar_fields):
    """Build dict of group indexes for filtered groups."""
    filter_indexes = {}
    for k in clusterbar_fields:
        filter_indexes.setdefault(k, {})
        groups = df[k].unique().tolist()
        if k in filters.keys():
            groups = filters[k]
        for elem in groups:
            obs_index = df.index[df[k] == elem]
            filter_indexes[k][elem] = obs_index.tolist()
    return filter_indexes

def create_dataframe_gene_mask(df, gene_symbols):
    """Create a gene mask to filter a dataframe."""
    if not "gene_symbol" in df:
        raise PlotError('Missing gene_symbol column in adata.var')

    gene_filter = None
    success = 1
    message = ""

    if not gene_symbols:
        return gene_filter, success, message

    try:
        # Some genes may map to multiple Ensembl IDs, which can cause issues.  Create a 1-to-1 mapping by dropping dups
        uniq_df = df.drop_duplicates(subset=['gene_symbol'])
        dataset_genes = df['gene_symbol'].unique().tolist()
        normalized_genes_list, found_genes = normalize_searched_genes(dataset_genes, gene_symbols)

        # Use our list of genes to get a single Ensembl ID for each gene
        uniq_gene_filter = uniq_df['gene_symbol'].isin(normalized_genes_list)
        genes_df = uniq_df['gene_symbol'][uniq_gene_filter]

        # NOTE: While volcanoes and quadrants can be searched without genes, this bit
        # should only execute if the user was searching volcanoes from the main page
        # Most likely the user is not interested in a annotation-less volcano plot from there.
        if genes_df.empty:
            raise PlotError("None of the searched gene symbols were found in this dataset.")

        # Now that our mapping is finished, create the gene filter
        gene_filter = df.index.isin(genes_df.index)

        # Get list of duplicated genes for the dataset
        gene_counts_df = df['gene_symbol'].value_counts().to_frame("count") # adding name to count ensures compatibility between pandas v1.5 and v2.0+
        dup_genes = gene_counts_df.index[gene_counts_df["count"] > 1].tolist()

        # Note to user which genes were duplicated.
        dup_genes_intersection = intersection(dup_genes, normalized_genes_list)

        message_list = []
        if dup_genes_intersection:
            success = 2
            message_list.append('<li>The following genes were mapped to 2 or more Ensembl IDs in this dataset, so one was chosen at random for the plot: {}</li>'.format(', '.join(dup_genes_intersection)))

        # Note to user which genes were not found in the dataset
        genes_not_present = [gene for gene in gene_symbols if gene not in found_genes]
        if genes_not_present:
            success = 2,
            message_list.append('<li>One or more genes were not found in the dataset nor could be mapped: {}</li>'.format(', '.join(genes_not_present)))
        message = "\n".join(message_list)
        return gene_filter, success, message
    except PlotError as pe:
        raise PlotError(str(pe))
    except Exception as e:
        # print stack trace
        import sys
        import traceback
        traceback.print_exc(file=sys.stderr)

        # Catch non-PlotError stuff
        raise PlotError("There was an issue searching genes in this dataset.")

def create_facet_indexes(groups):
    """Create facet indexes for subplots.  Returns a dict of group names to subplot index number."""
    return {group: idx for idx, group in enumerate(groups, start=1)}

def create_multicategory_axis_labels(groupby_filters, df):
    """ Creates the multicategory axis labels for a plot."""
    # If only one groupby column, we must flatten the list or else the x-axis will not plot correctly
    if len(groupby_filters) < 2:
        return df[groupby_filters[0]].tolist()
    multicategory = []
    for col in groupby_filters:
        multicategory.append(df[col].tolist())
    return multicategory

def intersection(lst1, lst2):
    """Intersection of two lists."""
    return list(set(lst1) & set(lst2))

def union(lst1, lst2):
    """Union of two lists."""
    return list(set(lst1) | set(lst2))

def normalize_searched_genes(gene_list, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [str(g) for cg in chosen_genes for g in gene_list if cg.lower() == str(g).lower()]
    found_genes = [cg for cg in chosen_genes for g in gene_list if cg.lower() == str(g).lower()]
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

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def get_discrete_colors(fields, colorscale="vivid", reverse_colorscale=False, ):
    """Get a list of discrete colors equal to the number of fields.

    Args:
        fields (list): List of fields to get colors for.
        colorscale (string, optional): Use the colorscale provided. Can be discrete or continuous. Defaults to "vivid".
        reverse_colorscale (bool, optional): If true, reverse the list of colors. Defaults to False.

    Raises:
        Exception: Passed in colorscale is not continuous but is not in the color swatch dictionary provided.

    Returns:
        list: List of colors.
    """

    if not colorscale:
        colorscale = "vivid"

    # If scale is sequential, split into equal colors equal to the number of categories
    # If scale is discrete, use the colorscale
    colors = None
    if colorscale.lower() in px.colors.named_colorscales():
        num_colors = len(fields)
        px.colors.sample_colorscale(px.colors.get_colorscale(colorscale), num_colors)
    elif colorscale not in color_swatch_map:
        # Not all the quantitivate colorscales available are in the color_swatch_map
        raise Exception("Colorscale {} not a valid colorscale to choose from".format(colorscale))
    else:
        colors = color_swatch_map[colorscale][::-1] if reverse_colorscale else color_swatch_map[colorscale]
    return colors

# NOTE: Currently not used as I refactored the code to use the existing retrieved colorscales
# But keeping this here in case I want to use it again in the future
def get_colorscale(colorscale):
    """Return colorscale 2D list for the selected premade colorscale."""
    if colorscale.lower() in px.colors.named_colorscales():
        return px.colors.get_colorscale(colorscale)
    # Must be a qualitative colorscale
    if colorscale not in color_swatch_map:
        raise Exception("Colorscale {} not a valid colorscale to choose from".format(colorscale))

    # Create a 2d list of colors for the selected colorscale at equal distances
    # to keep consistent with continuous colors
    colorscale_list = []
    length = len(color_swatch_map[colorscale]) - 1
    for i, color in enumerate(color_swatch_map[colorscale]):
        colorscale_list.append([i/length, color])
    return colorscale_list
