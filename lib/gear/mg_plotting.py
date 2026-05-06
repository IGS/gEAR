
import sys
from functools import partial
from itertools import cycle

import anndata as ad
import colorcet as cc
import dash_bio as dashbio
import diffxpy.api as de
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist

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

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

### Dotplot fxns

def create_dot_plot(df:pd.DataFrame, groupby_filters:list, is_log10:bool=False, plot_title:str|None=None, colorscale:str|None="Magma", reverse_colorscale:bool=False, non_interactive:bool=False):
    """Creates a dot plot.  Returns the figure."""
    # x = group
    # y = gene
    # color = mean expression
    # size = percent of cells with gene

    # Taking a lot of influence from
    # https://github.com/interactivereport/CellDepot/blob/ec067978dc456d9262c3c59d212d90547547e61c/bin/src/plotH5ad.py#L113

    fig = go.Figure()

    multicategory = create_multicategory_axis_labels(groupby_filters, df)

    # log-transform dataset if it came in raw
    mean = np.log2(df['mean'] + LOG_COUNT_ADJUSTER)
    if is_log10:
        mean = df['mean']

    # For some reason, the default argument is ignored.
    if not colorscale:
        colorscale="Magma"

    hover_template = "N: %{text}<br>Percent: %{marker.size:.2f}<br>Mean: %{marker.color:.2f}"
    if non_interactive:
        hover_template = None

    fig.add_scatter(
        x=multicategory
        , y=df["gene_symbol"]
        , text = df["count"]
        , hovertemplate=hover_template
        , mode="markers"
        , marker=dict(
            color=mean
            , colorscale=colorscale
            , reversescale=reverse_colorscale
            , size=df["percent"]
            , sizemode="area"
            , colorbar=dict(
                title=dict(
                    text="Log10 Mean Expression" if is_log10 else "Log2 Mean Expression"
                    , font=dict(family="Roboto", size=12)
                    , side="right"
                )
                , tickfont=dict(size=12)
                , len=0.75
                , thickness=15
                , x=1.05
                )
            )
        , showlegend=False
        )

    x_title = groupby_filters[0]
    if len(groupby_filters) > 1:
        x_title += " and {}".format(groupby_filters[1])
    fig.update_xaxes(
        title=x_title.capitalize(),
        domain=[0, 0.9]    # Take most of the plot width, but leave space on the right for the size legend
    )
    # If two x-axes labels are present, constrain the range
    if len(groupby_filters) < 2:
        uniq_cats = set(multicategory)
    else:
        # flatten the multicategory labels to get the unique categories across both axes
        uniq_cats = set([cat for sublist in multicategory for cat in sublist])
    if len(uniq_cats) == 2:
        fig.update_xaxes(range=[-0.5, 1.5])

    fig.update_yaxes(
        title="Genes",
        tickmode="linear",
        dtick=1,    # force tick for each category (gene)
        automargin=True
    )

    create_floating_dot_legend(fig)

    # Truncate faceted column axis labels so annotation can fit
    axis_label_mapping = {}  # Aggregated mapping of truncated -> full label names
    if not non_interactive and len(groupby_filters) == 1:
        # TODO: This does not work with 2D multicategory x-axis
        # For multi-gene plots, categoryarray is not always populated by Plotly automatically.
        # Derive unique ordered categories directly from the dataframe instead.
        x_categories = df[groupby_filters[0]].unique().tolist()  # preserves observed order

        def truncate_and_collect(a):
            # Fall back to dataframe-derived categories if axis hasn't populated categoryarray
            categories = list(a.categoryarray) if a.categoryarray is not None else x_categories
            ticktext, mapping = _truncate_ticktext(categories)
            axis_label_mapping.update(mapping)
            a.update(
                ticktext=ticktext,
                tickvals=categories,
            )
        fig.for_each_xaxis(truncate_and_collect)

        # Store the axis label mapping in the figure metadata for use on the JS side
        if axis_label_mapping:
            existing_meta = fig.layout.meta or {}
            if isinstance(existing_meta, dict):
                existing_meta["axis_label_mapping"] = axis_label_mapping
            else:
                existing_meta = {"axis_label_mapping": axis_label_mapping}
            fig.update_layout(meta=existing_meta)
    return fig

def create_floating_dot_legend(fig: go.Figure):
    """Adds a size legend using paper coordinates to avoid subplot whitespace. Edits in-place."""

    # Create a dot size legend
    steps = 5
    dot_legend=list()
    for i in range(steps):
        dot_legend += [["0", i, i*20+20, "{}%".format(i*20+20)]]
    dot_legend = pd.DataFrame(dot_legend,columns=['x','y','percent','text'])

    # Add dots as a scatter trace on 'paper' coordinates
    fig.add_scatter(
        x=[1 for i in dot_legend["x"].tolist()], # Positioned on the right side
        y=[0.1 + i*0.2 for i in dot_legend["y"].tolist()], # Stacked vertically
        xaxis="x2", yaxis="y2", # Use main axes for positioning
        mode="markers+text",
        text=dot_legend["text"],
        textposition="middle right",
        marker=dict(
            color="#888",
            size=dot_legend["percent"],
            sizemode="area"
        ),
        showlegend=False,
        hoverinfo="skip",
        cliponaxis=False
    )

# Configure the 'Legend' axes to be invisible and fixed
    fig.update_layout(
        xaxis2=dict(
            range=[0, 1],
            visible=False,
            overlaying="x",
            anchor="free",
            position=1,
            automargin=True
            ),
        yaxis2=dict(
            range=[0, 1],
            position=1,
            # setting visible=False hides the title
            showgrid=False,         # Hide the grid lines
            showline=False,         # Hide the vertical axis line
            showticklabels=False,   # Hide the numbers
            zeroline=False,         # Hide the baseline
            ticks="",               # no tick marks
            overlaying="y",
            anchor="free",
            side="right",
            # Use the axis title instead of an annotation
            title=dict(
                text="Percent of cells expressing gene",
                font=dict(size=12, family="Roboto"),
                standoff=0,
            ),
            ),
    )

### Heatmap fxns

def add_clusterbars(fig, all_categories):
    for i, categories in enumerate(all_categories):
        add_clusterbar(fig, categories, i)

def add_clusterbar(fig, metadata_series, bar_index):
    """Layers a clusterbar using a secondary axis. Edits in-place."""

    palette = get_categorical_palette(bar_index)

    unique_vals = metadata_series.unique()
    color_map = {val: palette[i] for i, val in enumerate(unique_vals)}

    # Calculate bar Y-position based on heatmap top
    # Assuming heatmap top is 0.8
    y_pos = 0.82 + (bar_index * 0.04)

    fig.add_scatter(
        x=metadata_series.index,
        y=[y_pos] * len(metadata_series),
        mode='markers',
        marker=dict(
            color=[color_map[v] for v in metadata_series],
            symbol='square',
            size=12
        ),
        name=metadata_series.name, # This becomes the Legend Header
        legendgroup=metadata_series.name,
        showlegend=True,
        xaxis="x", yaxis="y2", # Use an overlaying axis for bar positioning
        cliponaxis=False
    )

def add_left_dendrogram(fig, data, distance_metric:str="euclidean"):
    try:
        distfun = partial(pdist, metric=distance_metric)
    except Exception:
        raise ValueError(f"Invalid distance metric '{distance_metric}' for dendrogram. Please choose a valid metric from scipy.spatial.distance.pdist.")
    dendro = ff.create_dendrogram(data, orientation='right', distfun=distfun)

    # Add dendrogram traces to our existing figure
    for trace in dendro.data:
        trace.showlegend = False    # Removes "trace 1, trace 2..."
        trace.hoverinfo = 'skip'     # Modern "Calm" interaction
        trace.line.color = '#444444' # Single neutral color (dark gray)
        trace.xaxis = "x3" # Move to a dedicated dendrogram axis
        trace.yaxis = "y3"
        fig.add_trace(trace)

    # Align the Dendrogram domain ABOVE the heatmap domain
    fig.update_layout(
        xaxis3=dict(
            domain=[0.02, 0.18], # Sits in the left 16% of the card
            visible=False
        ),
        yaxis3=dict(
            domain=[0.1, 0.9], # Matches heatmap height
            visible=False
        ),
        # Shrink heatmap domain to avoid overlap
        xaxis=dict(domain=[0.2, 0.9]),
        # clear background
        plot_bgcolor='#FFFFFF',
    )

def add_top_dendrogram(fig, data, distance_metric:str="euclidean"):
    try:
        distfun = partial(pdist, metric=distance_metric)
    except Exception:
        raise ValueError(f"Invalid distance metric '{distance_metric}' for dendrogram. Please choose a valid metric from scipy.spatial.distance.pdist.")
    dendro = ff.create_dendrogram(data, orientation='bottom', distfun=distfun)

    # Add dendrogram traces to our existing figure
    for trace in dendro.data:
        trace.showlegend = False    # Removes "trace 1, trace 2..."
        trace.hoverinfo = 'skip'     # Modern "Calm" interaction
        trace.line.color = '#444444' # Single neutral color (dark gray)
        trace.xaxis = "x4" # Move to a dedicated dendrogram axis
        trace.yaxis = "y4"
        fig.add_trace(trace)

    # Align the Dendrogram domain ABOVE the heatmap domain
    fig.update_layout(
        xaxis4=dict(
            domain=[0.1, 0.9], # Matches heatmap width
            visible=False
        ),
        yaxis4=dict(
            domain=[0.82, 0.98], # Sits in the top 16% of the card
            visible=False
        ),
        # Shrink heatmap domain to avoid overlap
        yaxis=dict(domain=[0.1, 0.8]),
        # clear background
        plot_bgcolor='#FFFFFF',
    )

def calculate_domains(num_clusterbars, show_dendrogram=True):
    # Fixed heights in normalized paper units (0 to 1)
    bar_height = 0.03
    gap = 0.01
    dendro_height = 0.15 if show_dendrogram else 0

    # Calculate the top boundary of the heatmap
    heatmap_top = 1.0 - dendro_height - (num_clusterbars * (bar_height + gap))

    return {
        "heatmap": [0.1, heatmap_top],
        "bars_start": heatmap_top + gap,
        "dendro": [1.0 - dendro_height, 1.0]
    }

def create_heatmap(df:pd.DataFrame, groupby_filters: list, clusterbar_fields, is_log10:bool=False, cluster_obs:bool=False,
                    cluster_genes:bool=False, flip_axes:bool=False, center_around_zero:bool=False,
                    distance_metric:str="euclidean", colorscale:str|None=None, reverse_colorscale:bool=False,
                    title:str|None=None, hide_obs_labels:bool=False, hide_gene_labels:bool=False,
                    ) -> go.Figure:

    # df is long form
    # columns include:
    # - gene_symbol
    # - value (expression)
    # - any other groupby or clusterbar series value

    obs_groups = df.index.tolist()
    if groupby_filters:
        id_vars = set(groupby_filters + clusterbar_fields)

        # if id_vars is 3+, we need to create a composite index for the heatmap rows/columns to group by
        # Basically mutlicategory only goes 2 levels deep, so we need to take clusterbar fields
        # and do an index for each groupby filter so there is just 2 categories
        groupby_to_use = groupby_filters
        if len(id_vars) > 2:
            # Find clusterbar fields not in groupby_filters
            groupby_set = set(groupby_filters)
            clusterbar_set = set(clusterbar_fields)
            only_in_clusterbar = list(clusterbar_set - groupby_set)
            groupby_to_use = []
            if len(groupby_filters) > 2:
                # Sanity check.  UI handles this
                raise PlotError("Heatmap with more than 2 groupby filters is not supported. Please reduce the number of groupby filters to 2 or fewer.")

            # We only need the composite on the innermost multicategory
            index = 0
            if len(groupby_filters) == 2:
                index = 1
                groupby_to_use.append(groupby_filters[0])   # Append the outermost category
            composite_name = groupby_filters[index] + "_composite"
            composite_filters = [groupby_filters[index]] + only_in_clusterbar
            df[composite_name] = create_composite_index_column(df, composite_filters)
            groupby_to_use.append(composite_name)

        multicategory = create_multicategory_axis_labels(groupby_to_use, df)
        obs_groups = multicategory

    rows = list(df.index)
    col_labels = df["gene_symbol"] if flip_axes else obs_groups
    row_labels = obs_groups if flip_axes else df["gene_symbol"]

    x_visible = True
    y_visible = True
    if hide_obs_labels:
        if flip_axes:
            y_visible = False
        else:
            x_visible = False
    if hide_gene_labels:
        if flip_axes:
            x_visible = False
        else:
            y_visible = False

    values = df.loc[rows, "value"].to_numpy() + LOG_COUNT_ADJUSTER
    if is_log10:
        # Already log-10 transformed
        values = df.loc[rows, "value"].to_numpy()
    else:
        # log-transform to base-2
        values = np.log2(values)

    if not colorscale:
        colorscale = "RdYlBu" if center_around_zero else "Reds"
        # In the default colorscheme, reversing the scheme depends on if the plot centers around zero.
        reverse_colorscale = center_around_zero

    title_text = "Log2 Gene Expression"
    colorbar_title = "Log2 Expr."
    if is_log10:
        title_text = "Log10 Gene Expression"
        colorbar_title = "Log10 Expr."
    if title:
        title_text = title

    # Initialize a clean figure, not subplots
    fig = go.Figure()

    # --- STEP 1: Main Heatmap ---
    fig.add_heatmap(
        x=col_labels,
        y=row_labels,
        z=values,
        # Using your side-aligned colorbar title trick
        colorbar=dict(
            title=dict(
                font=dict(
                    size=12,
                    family="Roboto"
                ),
                text=colorbar_title,
                side="right"
                ),
            tickfont=dict(size=12),
            x=1.1,         # Sits in the right margin domain
            y=0.5,         # Centered vertically
            xanchor="center",
            yanchor="middle",
            thickness=15
        ),
        reversescale=reverse_colorscale,
        colorscale=colorscale,
        zmin=0 if not center_around_zero else None,
        zmax=float(values.max()) if not center_around_zero else None,
        zmid=0 if center_around_zero else None,
        name="expression"
    )

    fig.update_yaxes(
        tickmode="linear",
        dtick=1,    # force tick for each category
        automargin=True,
        side="right"    # So dendrogram doesn't overlap with labels
    )
    fig.update_xaxes(
        tickmode="linear",
        dtick=1,    # force tick for each category
        automargin=True
    )

    # --- STEP 2: The Domain Fix (The "Sidebar" Area) ---
    # This is how we eliminate the top/left whitespace when dendrograms are off.
    fig.update_layout(
        xaxis=dict(
            domain=[0.1, 0.9] # Heatmap takes 80% width, 10% left margin is safe.
            , visible=x_visible
        ),
        yaxis=dict(
            domain=[0.1, 0.9] # Heatmap takes 80% height, 10% top margin is safe.
            , visible=y_visible
        ),

        # Center the title
        title={"text": title_text, "x": 0.5, "xref": "paper", "y": 0.9}
    )

    show_top_dendrogram = False
    if cluster_obs and not flip_axes:
        show_top_dendrogram = True
    if cluster_genes and flip_axes:
        show_top_dendrogram = True
    calculate_domains(num_clusterbars=0, show_dendrogram=show_top_dendrogram)

    return fig

    try:
        # dendrogram needs 2d values
        values2d = values.reshape(len(row_labels), len(col_labels))

        if cluster_obs:
            if flip_axes:
                add_left_dendrogram(fig, values2d, distance_metric)
            else:
                add_top_dendrogram(fig, values2d, distance_metric)

        if cluster_genes:
            if flip_axes:
                add_top_dendrogram(fig, values2d, distance_metric)
            else:
                add_left_dendrogram(fig, values2d, distance_metric)
    except ValueError as e:
        # If an invalid distance metric was provided, raise an error with a clear message
        raise ValueError(str(e))

    # TODO: Add clusterbars

    return fig

def add_clustergram_cluster_bars(fig, clusterbar_indexes, obs_labels=None, is_log10=False, flip_axes=False) -> None:
    """Add column traces for each filtered group.  Edits figure in-place."""

    # Heatmap is located on xaxis11 and yaxis11 in dash-1.0.2
    # Left-side dendrogram is xaxis9 and yaxis9 in dash=1.0.2
    # Top-sie dendrogram is xaxis3 and yaxis3 in dash=1.0.2

    obs_axis = "xaxis11"
    gene_axis = "yaxis11"
    obs_dendro_axis = "xaxis9"
    gene_dendro_axis = "yaxis9"
    if flip_axes:
        obs_axis = "yaxis11"
        gene_axis = "xaxis11"
        obs_dendro_axis = "yaxis3"
        gene_dendro_axis = "xaxis3"

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

    # Put "groups" heatmap tracks either above or to the right of the genes in heatmap
    # Makes a small space b/t the genes and groups tracks
    next_bar_position = max(fig.layout[gene_axis]["tickvals"]) + 7

    # Information that will help posittion the "groups" colorbars
    curr_colorbar_x = 1.25

    for key, val in col_group_markers.items():

        z = create_clusterbar_z_value(flip_axes, groups_and_colors, key, val)

        # In order to make the colorscale a discrete one, we must map the start and stop thresholds for our normalized range
        colorscale = []
        # Normally tickvals will stretch to the min and max of the colorbar range. We also need to center the groups in the middle of the color
        tickvals = []

        colorlen = len(groups_and_colors[key]["colors"])
        for i in range(colorlen):
            # Start of color thresholds
            colorscale.append(( (i)/colorlen, groups_and_colors[key]["colors"][i] ))
            # End of color thresholds
            colorscale.append(( (i+1)/colorlen, groups_and_colors[key]["colors"][i] ))
            # Center the group name in its color
            tickvals.append((colorscale[-1][0] + colorscale[-2][0]) / 2 * (colorlen-1))

        next_x = obs_positions
        next_y = [next_bar_position-2, next_bar_position+2]
        if flip_axes:
            next_x = [next_bar_position-2, next_bar_position+2]
            next_y = obs_positions

        # Adjust the colorbar length based on the number of groups, but cap at 1
        colorbarlen = min(1, max(0.75, colorlen * 0.15))

        trace = go.Heatmap(
            x=next_x
            , y=next_y
            , z=z
            , colorbar=dict(
                len=colorbarlen
                , tickfont=dict(size=12)
                , ticktext=[group for group in groups_and_colors[key]["truncate"]]
                , tickvals=tickvals
                , thickness = 15
                , title=dict(font=dict(size=12, family="Roboto"), side="right", text=key)
                , x=curr_colorbar_x
                , y=1   # Align with bottom of heatmap
                , yanchor="top"
                )
            , colorscale=colorscale
            , name="clusterbar"
        )
        col_group_traces.append(trace)

        # Add group label to axis tuples
        fig.layout[gene_axis]["ticktext"] = fig.layout[gene_axis]["ticktext"] + (key, )
        fig.layout[gene_axis]["tickvals"] = fig.layout[gene_axis]["tickvals"] + (next_bar_position, )

        next_bar_position += 4 # add enough gap to space the "group" tracks
        curr_colorbar_x += 0.25

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

def create_clustergram_observation_labels(df, fig, colname="composite_index", flip_axes=False):
    """Create a set of labels to replace the current ticktext in the clustergram."""
    obs_axis = "yaxis11" if flip_axes else "xaxis11"
    obs_order = fig.layout[obs_axis]["ticktext"]    # Gets order of composite index observations
    return df.reindex(obs_order)[colname].tolist()  # reindex based on the observation order, and return the labels

### Quadrant fxns

def add_gene_annotations_to_quadrant_plot(fig, gene_symbols_list) -> tuple:
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

def prep_quadrant_dataframe(adata, key, control_val, compare1_val, compare2_val, de_test_algo="t_test", fc_threshold: float=2, fdr_threshold=0.05, include_zero_fc=True, is_log10=False):
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

def create_stacked_violin_plot(df: pd.DataFrame, groupby_filters:list, is_log10: bool=False, colorscale: str | None=None, reverse_colorscale: bool=False, non_interactive: bool=False):
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

    if not colors:
        colors = cc.glasbey

    color_cycler = cycle(colors)

    color_map = {cat: next(color_cycler) for cat in primary_groups}

    fig = make_subplots(
        rows=num_rows
        , cols=num_cols
        , row_titles=row_titles
        , column_titles=col_titles
        , shared_yaxes="all"    # to keep the scale the same for all row facets
        , x_title=None    # This title can sometimes run into gene names + it's implied these are genes
        , y_title="Log10 Expression" if is_log10 else "Log2 Expression"
        )

    groupby = ["gene_symbol"]
    groupby.extend(groupby_filters)
    # Create groupings for traces
    grouped = df.groupby(groupby, observed=True)

    # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
    # Group is the 'groupby' dataframe
    for name, group in grouped:
        # name[0] is gene_sym, name[1] is primary category, name[2] is secondary category
        # Normalize name to a tuple of strings and extract primary/secondary values; if missing, fall back to values from the group dataframe
        if isinstance(name, tuple):
            tuple_name = tuple(str(n) for n in name)
            #gene_sym = tuple_name[0] if len(tuple_name) > 0 else None
            primary_val = tuple_name[1] if len(tuple_name) > 1 else (group[groupby_filters[0]].iloc[0] if groupby_filters else None)
            secondary_val = tuple_name[2] if len(tuple_name) > 2 else None
        else:
            tuple_name = (str(name),)
            #gene_sym = tuple_name[0]
            primary_val = group[groupby_filters[0]].iloc[0] if groupby_filters else None
            secondary_val = group[groupby_filters[1]].iloc[0] if len(groupby_filters) > 1 else None

        row_idx = facet_row_indexes[primary_val]
        col_idx = facet_col_indexes[secondary_val] if secondary_val is not None else 1

        # log-transform dataset if it came in raw
        if not is_log10:
            group['value'] = np.log2(group['value'] + LOG_COUNT_ADJUSTER)

        fig.add_violin(
            x=group["gene_symbol"]
            , y=group["value"]
            , scalegroup=",".join(tuple_name) # Name will be a tuple
            , showlegend=False
            , fillcolor=color_map[primary_val]
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

        # Only show x-axis labels on the bottom row
        fig.update_xaxes(
            showticklabels=False if row_idx < num_rows else True
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

def create_violin_plot(df: pd.DataFrame, groupby_filters:list, is_log10: bool=False, colorscale: str | None=None, reverse_colorscale:bool=False, non_interactive:bool=False):
    """Creates a violin plot.  Returns the figure."""

    try:
        colors = get_discrete_colors(df["gene_symbol"].unique().tolist(), colorscale, reverse_colorscale)
    except Exception as e:
        raise PlotError("Error creating colors for violin plot: {}".format(e))

    if not colors:
        colors = cc.glasbey

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

        # If name is a tuple, create a string name for the scalegroup
        if isinstance(name, tuple):
            name = ','.join(str(n) for n in name)

        fig.add_violin(
            x=multicategory
            , y=group["value"]
            , name=gene_sym
            # Scalegroup must be unique for all violins to have max equal width.
            # Hence why we have to do a 3-dimensional groupby to make unique traces
            , scalegroup=name
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

    # Truncate faceted column axis labels so annotation can fit
    axis_label_mapping = {}  # Aggregated mapping of truncated -> full label names
    if not non_interactive and len(groupby_filters) == 1:
        # TODO: This does not work with 2D multicategory x-axis
        # For multi-gene plots, categoryarray is not always populated by Plotly automatically.
        # Derive unique ordered categories directly from the dataframe instead.
        x_categories = df[groupby_filters[0]].unique().tolist()  # preserves observed order

        def truncate_and_collect(a):
            # Fall back to dataframe-derived categories if axis hasn't populated categoryarray
            categories = list(a.categoryarray) if a.categoryarray is not None else x_categories
            ticktext, mapping = _truncate_ticktext(categories)
            axis_label_mapping.update(mapping)
            a.update(
                ticktext=ticktext,
                tickvals=categories,
            )
        fig.for_each_xaxis(truncate_and_collect)

        # Store the axis label mapping in the figure metadata for use on the JS side
        if axis_label_mapping:
            existing_meta = fig.layout.meta or {}
            if isinstance(existing_meta, dict):
                existing_meta["axis_label_mapping"] = axis_label_mapping
            else:
                existing_meta = {"axis_label_mapping": axis_label_mapping}
            fig.update_layout(meta=existing_meta)

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
        , selector=lambda a: not (a.yanchor == "bottom" or a.text.endswith("Expression"))
    )

    fig.for_each_annotation(
        # Edit y-axis title
        lambda a: a.update(
            xshift=-170 if len(max(primary_groups, key=len)) > 5 else -50    # Varies based on len of row_facet group names
        )
        , selector=lambda a: a.text.endswith("Expression")
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

def build_column_group_markers(filter_indexes, obs_order=None):
    """Build dictionaries of group annotations for the clustergram."""
    col_group_markers = {}
    # k = obs_category, elem = single observation, i = index position of indiv. observation in observation order list
    for k, v in filter_indexes.items():
        # If no obs order is provided, just assign unordered group annotations.
        if obs_order is None:
            col_group_markers[k] = {elem: v[elem] for elem in v}
            continue

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
        if filters and k in filters.keys():
            groups = filters[k]
        for elem in groups:
            obs_index = df.index[df[k] == elem]
            filter_indexes[k][elem] = obs_index.tolist()
    return filter_indexes

def create_dataframe_gene_mask(df, gene_symbols):
    """Create a gene mask to filter a dataframe."""
    if "gene_symbol" not in df:
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

def truncate(group):
    """Truncate group names to 10 characters."""
    return group[:10] + "..." if len(group) > 10 else group

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
        groups_and_colors[k]["truncate"] = [truncate(elem) for elem in v]
    return groups_and_colors

def create_composite_index_column(df, columns):
    """
    Create a composite index column by joining values from multiple columns.

    Args:
        df (pandas.DataFrame): The DataFrame containing the columns.
        columns (list): A list of column names to be joined.

    Returns:
        pandas.Series: A Series containing the composite index values.

    Example:
        >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        >>> create_composite_index_column(df, ['A', 'B'])
        0    1;4
        1    2;5
        2    3;6
        dtype: object
    """
    return df[columns].apply(lambda x: ';'.join(map(str, x)), axis=1)

def get_categorical_palette(index):
    # Cycle through Colorcet categorical palettes
    palettes = [cc.glasbey, cc.glasbey_dark, cc.glasbey_light, cc.glasbey_warm]
    return palettes[index % len(palettes)]

def get_discrete_colors(fields: list, colorscale: str | None="vivid", reverse_colorscale: bool=False):
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

def truncate_xaxis(fig, non_interactive=False, flip_axes=False):
    """ Truncate faceted column axis labels so annotation can fit"""

    obs_axis = "yaxis11" if flip_axes else "xaxis11"

    axis_label_mapping = {}  # Aggregated mapping of truncated -> full label names
    if not non_interactive and not flip_axes:
        # For multi-gene plots, categoryarray is not always populated by Plotly automatically.
        # Derive unique ordered categories directly from the dataframe instead.
        x_categories =  fig.layout[obs_axis]["ticktext"]

        def truncate_and_collect(a):
            # Fall back to dataframe-derived categories if axis hasn't populated categoryarray
            categories = list(a.categoryarray) if a.categoryarray is not None else x_categories
            ticktext, mapping = _truncate_ticktext(categories)
            axis_label_mapping.update(mapping)
            a.update(
                ticktext=ticktext,
            )
        fig.for_each_xaxis(truncate_and_collect)

        # Store the axis label mapping in the figure metadata for use on the JS side
        if axis_label_mapping:
            existing_meta = fig.layout.meta or {}
            if isinstance(existing_meta, dict):
                existing_meta["axis_label_mapping"] = axis_label_mapping
            else:
                existing_meta = {"axis_label_mapping": axis_label_mapping}
            fig.update_layout(meta=existing_meta)

def _truncate_ticktext(group_list: list[str]) -> tuple[list[str] | None, dict[str, str]]:
    """Truncate a group of axis ticks to a specified length."""
    TRUNCATION_LEN = 7  # How much of the original text to use (followed by ellipses)
    MAX_LEN_ALLOWED = 10  # Any text over this limit will be truncated

    # If only 0 or 1 datapoints in group, categoryarray was not present
    if not group_list:
        return None, {}

    new_ticktext = []
    full_name_mapping = {}
    truncated_counts: dict[str, int] = {}  # Track how many times a truncated label has been seen

    for val in group_list:
        if len(val) > MAX_LEN_ALLOWED:
            base_truncated = "{}...".format(val[0:TRUNCATION_LEN])

            if base_truncated in truncated_counts:
                # Collision: append a counter to disambiguate
                truncated_counts[base_truncated] += 1
                truncated = "{}~{}".format(base_truncated, truncated_counts[base_truncated])
            else:
                truncated_counts[base_truncated] = 0
                truncated = base_truncated

            new_ticktext.append(truncated)
            full_name_mapping[truncated] = val
        else:
            new_ticktext.append(val)
            #full_name_mapping[val] = val

    return new_ticktext, full_name_mapping
