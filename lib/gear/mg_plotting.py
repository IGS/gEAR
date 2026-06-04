
import sys
from functools import partial
from itertools import cycle

import anndata as ad
import colorcet as cc
import diffxpy.api as de
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list

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

    # In case colorscale is an empty string
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
            , cmin=0    # Anchor colorscale min at 0
            , cmax=float(mean.max())    # Required if cmin is set.
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
    x_title = underscore_to_space(x_title)
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

def add_clusterbars(fig: go.Figure, obs_columns, all_categories: list, bar_start_pos: float, flip_axes: bool=False, pivot_cols=None, cluster_obs=False):
    curr_bar_pos = bar_start_pos
    curr_legend_pos = 1.2

    # Replicate the same obs
    if cluster_obs:
        obs_groups = obs_columns.map(lambda x: ";".join(str(i) for i in x)).tolist()
    else:
        obs_groups = build_multicategory_obs_labels(obs_columns)
    x = obs_groups

    for i, field in enumerate(all_categories):

        if field not in obs_columns.names:
            raise PlotError(f"Clusterbar field '{field}' not found in heatmap columns. Please ensure all clusterbar fields are included in the heatmap columns.")

        categories = obs_columns.get_level_values(field)

        y = [curr_bar_pos]
        text = [categories]

        palette = get_categorical_palette(i)
        unique_vals = list(dict.fromkeys(categories))
        # map colors to categories
        color_map = {val: i for i, val in enumerate(unique_vals)}

        # In order to make the colorscale a discrete one, we must map the start and stop thresholds for our normalized range
        # You need the stop threshold, otherwise Plotly will try to interpolate some colors and be wrong
        colorscale = []
        tickvals = []
        colorlen = len(unique_vals)
        for i in range(colorlen):
            # Start of color thresholds
            colorscale.append(( (i)/colorlen, palette[i % len(palette)] ))
            # End of color thresholds
            colorscale.append(( (i+1)/colorlen, palette[i % len(palette)] ))
            # Center the group name in its color
            # SAdkins - I honestly don't remember how I arrived at this formula, but it works
            tickvals.append((colorscale[-1][0] + colorscale[-2][0]) / 2 * (colorlen-1))

        # Map groups to their integer index in the palette
        z = [[color_map[v] for v in categories]]
        if flip_axes:
            z = np.array(z).T.tolist()  # Transpose to make it a column vector again

        # Unfortunately, since we want unique colorscales per heatmap, we need a unique 1-row heatmap per field.
        fig.add_heatmap(
            x=y if flip_axes else x,
            y=x if flip_axes else y,
            z=z,
            text=text,
            # Define the discrete colorscale for this specific metadata field
            colorscale=colorscale,
            showscale=True,  # Kills the extra colorbar
            xaxis="x",
            yaxis="y",       # Still use yaxis2 for the domain sandwich
            hoverongaps=False,
            hovertemplate=f"{underscore_to_space(field)}: %{{text}}<extra></extra>",
            name="clusterbar",
            colorbar=dict(
                len=0.75,
                x=curr_legend_pos,         # Sits in the right margin domain
                y=0.5,         # Centered vertically
                xanchor="center",
                yanchor="middle",
                thickness=15,
                title=dict(
                    font=dict(
                        size=12,
                        family="Roboto"
                    ),
                    text=underscore_to_space(field),
                    side="right"
                ),
                tickfont=dict(size=10),
                ticktext=map_underscore_to_space(unique_vals),
                tickvals=tickvals
            )
        )

        curr_bar_pos += 1  # Move to the next bar position   (bar_height + gap)
        curr_legend_pos += 0.25  # Move to the next legend position

def add_left_dendrogram(fig, data, distfun:callable, dendro_domain:dict|None=None):
    dendro = ff.create_dendrogram(data, orientation='right', distfun=distfun)

    # Calculate the internal coordinate limits
    # ff.create_dendrogram leaves are always at 5, 15, 25...
    num_leaves = len(data)
    max_coord = num_leaves * 10

    # Add dendrogram traces to our existing figure
    for trace in dendro.data:
        trace.showlegend = False    # Removes "trace 1, trace 2..."
        trace.hoverinfo = 'skip'     # Modern "Calm" interaction
        trace.line.color = '#444444' # Single neutral color (dark gray)
        trace.xaxis = "x2" # Move to a dedicated dendrogram axis
        trace.yaxis = "y2"
        fig.add_trace(trace)

    dendro_x_domain = (dendro_domain or {}).get("x", {}).get("left", [0.02, 0.2])
    dendro_y_domain = (dendro_domain or {}).get("y", {}).get("left", fig.layout.yaxis.domain)  # Match heatmap height by default

    # Align the Dendrogram domain ABOVE the heatmap domain
    fig.update_layout(
        xaxis2=dict(
            domain=dendro_x_domain, # Sits in the left 16% of the card
            visible=False,
        ),
        yaxis2=dict(
            domain=dendro_y_domain, # Matches heatmap height
            range=[0, max_coord],
            visible=False
        ),
    )

def add_top_dendrogram(fig, data, distfun:callable, dendro_domain:dict|None=None):
    dendro = ff.create_dendrogram(data, orientation='bottom', distfun=distfun)

    num_leaves = len(data)
    max_coord = num_leaves * 10

    # Add dendrogram traces to our existing figure
    for trace in dendro.data:
        trace.showlegend = False    # Removes "trace 1, trace 2..."
        trace.hoverinfo = 'skip'     # Modern "Calm" interaction
        trace.line.color = '#444444' # Single neutral color (dark gray)
        trace.xaxis = "x3" # Move to a dedicated dendrogram axis
        trace.yaxis = "y3"
        fig.add_trace(trace)

    dendro_x_domain = (dendro_domain or {}).get("x", {}).get("top", fig.layout.xaxis.domain)  # Match heatmap width by default
    dendro_y_domain = (dendro_domain or {}).get("y", {}).get("top", [0.77, 0.95])

    # Align the Dendrogram domain ABOVE the heatmap domain
    fig.update_layout(
        xaxis3=dict(
            domain=dendro_x_domain, # Matches heatmap width
            range=[0, max_coord],
            visible=False
        ),
        yaxis3=dict(
            domain=dendro_y_domain, # Sits in the top 16% of the card
            visible=False
        ),
    )

def build_multicategory_obs_labels(obs_columns):
    obs_groups = obs_columns
    if isinstance(obs_columns, pd.MultiIndex):
        # Attempt to make this a multicategory axis of 2 levels (which is the plotly max)
        if obs_columns.nlevels > 2:
            # Collapse level 1 onwards into a composite string, keeping level 0 as outermost
            level_0 = list(obs_columns.get_level_values(0))
            inner_values = [
                ";".join(str(obs_columns[i][j]) for j in range(1, obs_columns.nlevels))
                for i in range(len(obs_columns))
            ]
            obs_groups = [level_0, inner_values]
        else:
            # Convert to multicategory format: list of lists, one per level
            obs_groups = [list(obs_columns.get_level_values(i)) for i in range(obs_columns.nlevels)]
    else:
        obs_groups = obs_columns.tolist()
    return obs_groups

def calculate_x_domains(show_dendrogram=True, gene_domain_ratio=1.0):
    # Fixed heights in normalized paper units (0 to 1)
    dendro_width = 0.18 if show_dendrogram else 0

    # Calculate the start boundary of the heatmap
    heatmap_left = 0.02 + dendro_width if show_dendrogram else 0.05
    heatmap_right = 0.9

    # This should be equal to the heatmap with only the genes
    # And also each leave should sit halfway in the corresponding heatmap cell
    dendro_top_width = gene_domain_ratio * (heatmap_right - heatmap_left)
    dendro_top_left_pos = heatmap_left
    dendro_top_right_pos = heatmap_left + dendro_top_width

    return {
        "heatmap": [heatmap_left, heatmap_right],
        "dendro": {
                "left": [0.02, heatmap_left],
                "top": [dendro_top_left_pos, dendro_top_right_pos]
        }
    }

def calculate_y_domains(show_dendrogram=True, gene_domain_ratio=1.0):
    # Fixed heights in normalized paper units (0 to 1)
    dendro_height = 0.18 if show_dendrogram else 0

    # Calculate the top boundary of the heatmap
    heatmap_bottom = 0.1
    heatmap_top = 0.95 - dendro_height

    # This should be equal to the heatmap with only the genes
    # And also each leave should sit halfway in the corresponding heatmap cell
    dendro_left_height = gene_domain_ratio * (heatmap_top - heatmap_bottom)
    dendro_left_bottom_pos = heatmap_bottom
    dendro_left_top_pos = heatmap_bottom + dendro_left_height

    return {
        "heatmap": [heatmap_bottom, heatmap_top],
        "dendro": {
            "left": [dendro_left_bottom_pos, dendro_left_top_pos],
            "top": [0.95 - dendro_height, 0.95]
        }
    }

def create_heatmap(df:pd.DataFrame, groupby_filters:list=[], clusterbar_fields:list=[], is_log10:bool=False, cluster_obs:bool=False,
                    cluster_genes:bool=False, flip_axes:bool=False, center_around_zero:bool=False,
                    distance_metric:str="euclidean", colorscale:str|None="cividis", reverse_colorscale:bool=False,
                    title:str|None=None, hide_obs_labels:bool=False, hide_gene_labels:bool=False,
                    ) -> go.Figure:

    # df is long form
    # columns include:
    # - gene_symbol
    # - value (expression)
    # - any other groupby or clusterbar series value

    rows = list(df.index)

    values = df.loc[rows, "value"].to_numpy() + LOG_COUNT_ADJUSTER
    if is_log10:
        # Already log-10 transformed
        values = df.loc[rows, "value"].to_numpy()
    else:
        # log-transform to base-2
        values = np.log2(values)

    df["value"] = values

    pivot_cols = df.index.name
    id_vars = set(groupby_filters + clusterbar_fields)
    if len(id_vars):
        # Sets destroy order, so this preserves it.
        pivot_cols = list(dict.fromkeys(groupby_filters + clusterbar_fields))

    # df is now wide form
    # index is gene_symbol
    # multi-index column of groupby cats and clusterbar fields, or the original index by default.
    pivot_df = df.pivot_table(index="gene_symbol", columns=pivot_cols, values="value", aggfunc="mean", observed=False)
    values_2d = pivot_df.to_numpy()

    # Initialize a clean figure
    fig = go.Figure()

    show_left_dendrogram = False
    show_top_dendrogram = False
    if cluster_obs:
        if flip_axes:
            show_left_dendrogram = True
        else:
            show_top_dendrogram = True
    if cluster_genes:
        if flip_axes:
            show_top_dendrogram = True
        else:
            show_left_dendrogram = True

    num_genes = len(pivot_df.index)
    # This is the ratio of genes to clusterbars in the final heatmap.
    gene_domain_ratio = num_genes / (num_genes + len(clusterbar_fields))
    x_gene_domain_ratio = gene_domain_ratio if flip_axes else 1.0
    y_gene_domain_ratio = gene_domain_ratio if not flip_axes else 1.0

    domain_x_dict = calculate_x_domains(show_dendrogram=show_left_dendrogram, gene_domain_ratio=x_gene_domain_ratio)
    domain_y_dict = calculate_y_domains(show_dendrogram=show_top_dendrogram, gene_domain_ratio=y_gene_domain_ratio)
    dendro_dict = {
        "x": domain_x_dict["dendro"],
        "y": domain_y_dict["dendro"]
    }

    # Dendrograms, if present, should be the source of truth for the axis order.
    # So plot them first.
    if cluster_obs or cluster_genes:
        try:
            distfun = partial(pdist, metric=distance_metric)    # type: ignore
        except Exception:
            raise ValueError(f"Invalid distance metric '{distance_metric}' for dendrogram. Please choose a valid metric from scipy.spatial.distance.pdist.")

        if cluster_genes:
            # use same default linkage method as dash-bio clustergrame
            # "complete" linkage considers the furthest distnace between groups for clustering.
            linkage_matrix = linkage(distfun(values_2d), method="complete")
            # 'leaves' gives you the order of indices [5, 2, 0, 11...]
            gene_order = leaves_list(linkage_matrix)
            # Reorder dataframe
            pivot_df = pivot_df.iloc[gene_order]

            if flip_axes:
                add_top_dendrogram(fig, values_2d, distfun, dendro_dict)
            else:
                add_left_dendrogram(fig, values_2d, distfun, dendro_dict)

        if cluster_obs:
            # use same default linkage method as dash-bio clustergrame
            linkage_matrix = linkage(distfun(values_2d.T), method="complete")
            obs_order = leaves_list(linkage_matrix)

            pivot_df = pivot_df.iloc[:, obs_order]

            if flip_axes:
                add_left_dendrogram(fig, values_2d.T, distfun, dendro_dict)
            else:
                add_top_dendrogram(fig, values_2d.T, distfun, dendro_dict)

    # Get a int code for the pivot_df index, so we can use them for the tickvals
    gene_symbol_codes = [i for i, gene in enumerate(pivot_df.index)]

    # Get the observation groups from the pivot table columns.  This will be used for the axis labels and hover info.
    if cluster_obs:
        # Convert the multiindex column tuples into a string
        obs_groups = pivot_df.columns.map(lambda x: ";".join(str(i) for i in x)).tolist()
    else:
        obs_groups = build_multicategory_obs_labels(pivot_df.columns)

    # Now that pivot_df has been sorted, we need to grab values again.
    values_2d = pivot_df.to_numpy()

    x = gene_symbol_codes if flip_axes else obs_groups
    y = obs_groups if flip_axes else gene_symbol_codes
    z = values_2d.T if flip_axes else values_2d

    # If colorscale is empty string
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

    # Generate the heatmap
    fig.add_heatmap(
        x=x,
        y=y,
        z=z,
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
            len=0.75,
            x=1.1 if flip_axes else 1.05,         # Sits in the right margin domain
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

    # We need to use ticks for the axis the clusterbar entries share. That way there is a unified coordinate system
    # and we have full control over the position of things.
    col_ticktext = None
    col_tickvals = None
    row_ticktext = None
    row_tickvals = None
    bars_start = None
    gene_domain_ratio = 1.0
    if flip_axes:
        col_ticktext = pivot_df.index.tolist()
        col_tickvals = gene_symbol_codes
        bars_start = len(col_ticktext)
        bar_tick = bars_start

        for field in clusterbar_fields:
            col_ticktext.append(field)
            col_tickvals.append(bars_start)
            bar_tick += 1

    else:
        row_ticktext = pivot_df.index.tolist()
        row_tickvals = gene_symbol_codes
        bars_start = len(row_ticktext)
        bar_tick = bars_start

        for field in clusterbar_fields:
            row_ticktext.append(field)
            row_tickvals.append(bars_start)
            bar_tick += 1

    # Add clusterbars if they are needed.
    if clusterbar_fields:
        add_clusterbars(fig, pivot_df.columns, clusterbar_fields, bars_start, flip_axes, pivot_cols, cluster_obs)


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

    # This is how we eliminate the top/left whitespace when dendrograms are off.
    fig.update_layout(
        xaxis=dict(
            domain=domain_x_dict["heatmap"] # Heatmap takes 80% width, 10% left margin is safe.
            , visible=x_visible
        ),
        yaxis=dict(
            domain=domain_y_dict["heatmap"] # Heatmap takes 80% height, 10% top margin is safe.
            , visible=y_visible
        ),

        # Center the title
        title={"text": title_text, "x": 0.5, "xref": "paper", "y": 0.9}
    )

    # Set up layout
    if flip_axes:
        fig.update_layout(
            xaxis=dict(
                tickmode="array",
                tickvals=col_tickvals,
                ticktext=map_underscore_to_space(col_ticktext),
                automargin=True
            ),
            yaxis=dict(
                tickmode="linear",
                dtick=1,    # force tick for each category
                automargin=True,
                side="right"    # So dendrogram doesn't overlap with labels
            )
        )
    else:
        fig.update_layout(
            xaxis=dict(
                tickmode="linear",
                dtick=1,    # force tick for each category
                automargin=True
            ),
            yaxis=dict(
                tickmode="array",
                tickvals=row_tickvals,
                ticktext=map_underscore_to_space(row_ticktext),
                automargin=True,
                side="right"    # So dendrogram doesn't overlap with labels
            )
        )

    # This is a band-aid fix to tighten the margins if dendrograms are not present.
    # The dendrograms add extra margins that we need to account for, but if they are not there, we can reclaim that space.
    r_margin=80 + (40 * len(clusterbar_fields))
    fig.update_layout(
        margin=dict(l=10, t=70, b=30, r=r_margin),
        # clear background (mostly to override dendro defaults)
        plot_bgcolor='#FFFFFF',
    )

    return fig


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

    pretty_print_compare1_val = underscore_to_space(compare1_val)
    pretty_print_compare2_val = underscore_to_space(compare2_val)
    pretty_print_control_val = underscore_to_space(control_val)

    fig.update_xaxes(title="{} vs {} log2FC".format(pretty_print_compare1_val, pretty_print_control_val))
    fig.update_yaxes(title="{} vs {} log2FC".format(pretty_print_compare2_val, pretty_print_control_val))
    fig.update_layout(
        legend_title_text="Log2FC: Num Genes in Group"
        )
    return fig

def prep_quadrant_dataframe(adata, key, control_val, compare1_val, compare2_val, de_test_algo="t-test", fc_threshold: float=2, fdr_threshold=0.05, include_zero_fc=True, is_log10=False):
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
    # Wanted to use Wald test, which is generally better. But it requires raw counts
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
    x_title = underscore_to_space(x_title)
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
    plot_title = underscore_to_space(plot_title)
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
        # Insignificant genes are at index 2.  If you want to skip annotating them, end at index 1
        range_end = len(fig.data) if annot_nonsig else 2
        for data_idx in range(0, range_end):
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
                        , bgcolor=fig.data[data_idx]["marker"]["color"] if data_idx < 2 else "slategrey"
                        , borderpad=2
                        , showarrow=True
                        , text=gene
                        , x=fig.data[data_idx].x[idx]
                        , y=fig.data[data_idx].y[idx]
                        , xref="x"
                        , yref="y"
                    )

def create_volcano_plot(df, query, ref, pval_threshold, logfc_bounds, use_adj_pvals=False, downcolor="#636EFA", upcolor="#EF553B"):
    """Generate a volcano plot.  Returns Plotly figure."""

    # In case None is passed in.
    downcolor = downcolor or "#636EFA"
    upcolor = upcolor or "#EF553B"

    color_map = {
        'Upregulated': upcolor,
        'Downregulated': downcolor,
        'Not Significant': '#D3D3D3'  # Light Gray
    }

    status_to_name = {
        'Upregulated': f"Upregulated in {underscore_to_space(query)}",
        'Downregulated': f"Upregulated in {underscore_to_space(ref)}",
        'Not Significant': "Nonsignificant Genes"
    }

    fig = go.Figure()

    for status, color in color_map.items():
        subset = df[df['Status'] == status]
        fig.add_trace(go.Scatter(
            x=subset['logfoldchanges'],
            y=subset['nlog10'],
            mode='markers',
            name=status_to_name[status],
            customdata=subset['ensm_id'], # Add ensembl ids to customdata for potential use in JS hover/click interactions
            text=subset['gene_symbol'],
            marker=dict(color=color, size=6, opacity=0.7),
            # Customizing the hover data for gene exploration
            hovertemplate=(
                "<b>%{text}</b><br>"
                "Ensembl ID: %{customdata}<br>"
                "Log2FC: %{x:.2f}<br>"
                "-Log10(p): %{y:.2f}"
            )
        ))

    # Add Threshold Lines
    line_styles = dict(line_dash="dash", line_color="black", opacity=0.4, line_width=1)
    fig.add_hline(y=-np.log10(pval_threshold), **line_styles)
    fig.add_vline(x=logfc_bounds[0], **line_styles)
    fig.add_vline(x=logfc_bounds[1], **line_styles)

    # Final Layout Adjustments
    yaxis_title = "-log10(adjusted-P)" if use_adj_pvals else "-log10(P)"
    pretty_print_query = underscore_to_space(query)
    pretty_print_ref = underscore_to_space(ref)

    fig.update_layout(
        title=f"Differences in {pretty_print_query} vs {pretty_print_ref}",
        xaxis_title='Log<sub>2</sub> Fold Change',
        yaxis_title=yaxis_title,
        template='plotly_white',
        hovermode='closest',
        font=dict(family="Roboto", size=12)
    )

    return fig

def prep_volcano_dataframe(adata, key, query_val, ref_val, de_test_algo="t-test", is_log10=False, use_adj_pvals=False, pval_threshold=0.05, logfc_bounds=[-1, 1]):
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

    # Use diffxpy to compute DE statistics for each comparison
    # Wanted to use Wald test, which is generally better. But it requires raw counts
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

    key = "pvals_adj" if use_adj_pvals else "pvals"
    df['nlog10'] = -np.log10(df[key])

    # logfc_bounds = [lower, upper]
    nplog10_pval_threshold = -np.log10(pval_threshold)
    conditions = [
        (df['logfoldchanges'] >= logfc_bounds[1]) & (df['nlog10'] >= nplog10_pval_threshold),
        (df['logfoldchanges'] <= logfc_bounds[0]) & (df['nlog10'] >= nplog10_pval_threshold)
    ]

    # Split data into 3 categories: upregulated, downregulated, and not significant
    choices = ['Upregulated', 'Downregulated']
    df['Status'] = np.select(conditions, choices, default='Not Significant')

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

def create_multicategory_axis_labels(groupby_filters, df, return_unique=False):
    """ Creates the multicategory axis labels for a plot."""
    if return_unique:
        # Get the unique values for each groupby filter to create the category arrays for the x-axis.  This also preserves the order of the categories as they appear in the dataframe.
        df = df[groupby_filters].drop_duplicates(keep='first')
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

def normalize_searched_genes(gene_list, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [str(g) for cg in chosen_genes for g in gene_list if cg.lower() == str(g).lower()]
    found_genes = [cg for cg in chosen_genes for g in gene_list if cg.lower() == str(g).lower()]
    return case_insensitive_genes, found_genes


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
    palettes = [cc.glasbey_dark, cc.glasbey_cool, cc.glasbey_warm, cc.glasbey_hv]
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

def map_underscore_to_space(items):
    """Map underscores in a list of strings to spaces."""
    return [underscore_to_space(item) for item in items]

def underscore_to_space(string):
    """Convert underscores in a string to spaces."""
    return string.replace("_", " ")

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
        # Ensure val is a string for length checking and mapping.
        val = str(val)
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
