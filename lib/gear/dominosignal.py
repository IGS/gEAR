"""
Loading and plotting of DominoSignal cell-cell communication output.

Ported from notebooks/domino_plotting.ipynb.  Each plot function returns a plotly Figure,
or raises PlotSkipped with a user-facing reason if the plot can't be drawn for the given
selections.

DominoSignal output files are .parquet files stored in www/cache/domino_signal.  Each may have
a .json file of the same name alongside it, describing the dataset (see load_dataset_info).
"""

import json
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

DOMINO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'www', 'cache', 'domino_signal'))

# DominoSignal output columns are read by position (as in the notebook) and given canonical names
POSITIONAL_COLUMNS = {
    0: "ligand",
    1: "receptor",
    2: "tf",
    3: "ligand_exp",
    4: "receptor_exp",
    6: "sending_cl",
    7: "receiving_cl",
}

# Column headers shown above each river plot stage
STAGE_NAMES = {
    "sending_cl": "Sending Cluster",
    "ligand": "Ligand",
    "receptor": "Receptor",
    "tf": "TF",
    "receiving_cl": "Receiving Cluster",
}

DEFAULT_EXP_RANGE_THRESHOLD = 0.5
DEFAULT_EXP_MAX_THRESHOLD = 1.5

### Color palettes ###
# River links are aggregated (one link per unique node pair, sized by interaction count), so
#  link colors use a fixed opacity rather than the notebook's density-based per-row opacity

# Palette 1 - default colors
PALETTE_1 = {
    "dot_colorscale": "Inferno_r",
    "river_node_color": "blue",
    "river_link_color": "rgba(255,127,0,0.4)",
    "river_link_hover_color": "blue",
}

# Palette 2 - color-blind safe colors, matching gEAR's colorblind mode elsewhere
#  (reversed cividis for continuous scales, dark blue/yellow for two-color contrasts).
#  Links use a higher opacity so the yellow stays visible on a white background.
PALETTE_2 = {
    "dot_colorscale": "cividis_r",
    "river_node_color": "rgb(0,34,78)",
    "river_link_color": "rgba(254,232,56,0.7)",
    "river_link_hover_color": "rgb(0,34,78)",
}

RIVER_MIN_HEIGHT = 500
RIVER_MAX_HEIGHT = 4000
RIVER_HEIGHT_PER_NODE = 15

# Space for the title, axis titles and (rotated) tick labels on top of the per-row dot plot height
DOT_PLOT_PADDING = 300
# Plots with few rows still need room for the rotated tick labels, else Plotly fails to scale the axes
DOT_PLOT_MIN_HEIGHT = 500

# Dot sizes at a dot size scale of 1.  The TF dot plot uses Plotly's default marker size.
DOT_TF_MARKER_SIZE = 6
DOT_CLUSTERS_MARKER_SIZE = 10
DOT_SIZE_SCALE_MIN = 0.25
DOT_SIZE_SCALE_MAX = 3.0

# Thin outline around each dot on the white background, so the palest (lowest expression) dots stay visible
DOT_PLOT_WHITE_BACKGROUND_MARKER_LINE = dict(color="#888888", width=0.5)


class PlotSkipped(Exception):
    """Raised when a plot can't be drawn with the given selections."""


class DominoFileError(Exception):
    """Raised when a requested DominoSignal file is invalid or missing."""


### LOAD ###

def list_domino_files():
    if not os.path.isdir(DOMINO_DIR):
        return []

    return sorted(f for f in os.listdir(DOMINO_DIR) if f.endswith('.parquet'))


def resolve_domino_path(filename):
    """Returns the full path for a bare filename within DOMINO_DIR, refusing anything outside it."""
    if not filename:
        raise DominoFileError("No DominoSignal file provided")

    filename = os.path.basename(filename)
    if not filename.endswith('.parquet'):
        raise DominoFileError("Only .parquet files are supported")

    domino_dir = os.path.realpath(DOMINO_DIR)
    path = os.path.realpath(os.path.join(domino_dir, filename))
    if not path.startswith(domino_dir + os.sep) or not os.path.isfile(path):
        raise DominoFileError("DominoSignal file not found: {0}".format(filename))

    return path


def load_domino_df(filename):
    domino_df = pd.read_parquet(resolve_domino_path(filename), engine='pyarrow')
    domino_df = domino_df.iloc[:, list(POSITIONAL_COLUMNS.keys())]
    domino_df.columns = list(POSITIONAL_COLUMNS.values())
    return domino_df


# Keys in a DominoSignal output's dataset info JSON, mapped to the names returned to the page
DATASET_INFO_FIELDS = {
    "Title": "title",
    "Author": "author",
    "Year": "year",
    "DOI": "doi",
    "PMID": "pmid",
    "gEAR Dataset URL": "gear_url",
    "Description": "description",
}


def load_dataset_info(filename):
    """
    Details of the paper and dataset behind a DominoSignal output, from the .json file with the same
    name as its .parquet file (e.g. domino_network_p21.reformatted.json).  Keys are matched
    case-insensitively, and blank values are returned as None.

    Returns None if there is no such file, or it can't be read as a JSON object.
    """
    json_path = os.path.splitext(resolve_domino_path(filename))[0] + ".json"
    if not os.path.isfile(json_path):
        return None

    try:
        with open(json_path) as json_file:
            raw_info = json.load(json_file)
    except (OSError, ValueError):
        return None

    if not isinstance(raw_info, dict):
        return None

    raw_info = {str(key).strip().lower(): value for key, value in raw_info.items()}
    dataset_info = {}
    for json_key, name in DATASET_INFO_FIELDS.items():
        value = raw_info.get(json_key.lower())
        value = str(value).strip() if value is not None else ""
        dataset_info[name] = value or None

    return dataset_info


def get_domino_options(domino_df):
    """Unique values for each selectable parameter."""
    return {
        "tfs": sorted(domino_df["tf"].unique().tolist()),
        "ligands": sorted(domino_df["ligand"].unique().tolist()),
        "receptors": sorted(domino_df["receptor"].unique().tolist()),
        "sending_clusters": sorted(domino_df["sending_cl"].unique().tolist()),
        "receiving_clusters": sorted(domino_df["receiving_cl"].unique().tolist()),
    }


### Dot plots ###

def build_dot_chart_data(domino_df):
    # Rows differing only by TF share the same cluster pair, ligand/receptor pair and
    #  expression, so they are collapsed into a single dot
    return pd.DataFrame({
        "clusters": domino_df["sending_cl"] + '->' + domino_df["receiving_cl"],
        "lig_rec": domino_df["ligand"] + '->' + domino_df["receptor"],
        # Geometric mean of ligand and receptor expression
        "mean_exp": np.sqrt(domino_df["ligand_exp"].astype(float) * domino_df["receptor_exp"].astype(float)),
    }).drop_duplicates()


def get_lig_rec_pairs_passing_thresholds(dot_chart_data, exp_range_threshold, exp_max_threshold):
    """
    Ligand/receptor pairs where either the range of geometric mean expression across cluster
    pairs exceeds exp_range_threshold, or any single value exceeds exp_max_threshold.
    """
    lig_rec_stats = dot_chart_data.groupby("lig_rec")["mean_exp"].agg(["min", "max"])
    return lig_rec_stats[(lig_rec_stats["max"] - lig_rec_stats["min"] > exp_range_threshold) | (lig_rec_stats["max"] > exp_max_threshold)].index


def thresholds_title_suffix(exp_range_threshold, exp_max_threshold):
    # On a second line of the title, so long titles fit on narrower screens
    return "<br>Expression Range Threshold = " + str(exp_range_threshold) + ", High Expression Threshold = " + str(exp_max_threshold)


def style_dot_plot(dot_fig, dot_chart_data, title, tickangle, row_height, marker_size, white_background=False):
    n_y = dot_chart_data["lig_rec"].nunique()
    n_x = dot_chart_data["clusters"].nunique()

    dot_fig.update_traces(marker=dict(size=marker_size))

    if white_background:
        dot_fig.update_traces(marker=dict(line=DOT_PLOT_WHITE_BACKGROUND_MARKER_LINE))
        dot_fig.update_layout(plot_bgcolor="white")
        # A clean, plain white background: no grid lines within the plot, but keep axis lines and ticks
        dot_fig.update_xaxes(showgrid=False, showline=True, linecolor="black", ticks="outside")
        dot_fig.update_yaxes(showgrid=False, showline=True, linecolor="black", ticks="outside")

    # Make it look like a matrix
    dot_fig.update_layout(
        title=title,
        xaxis_title="Sending -> Receiving Cluster",
        yaxis_title="Ligand -> Receptor Pair",
        # Sorting and adjusting the unnecessary space above and below
        xaxis=dict(categoryorder="category ascending", range=[-0.5, n_x - 0.5]),
        yaxis=dict(categoryorder="category descending", range=[-0.5, n_y - 0.5]),
        xaxis_tickangle=tickangle,
        height=max(DOT_PLOT_MIN_HEIGHT, row_height * n_y + DOT_PLOT_PADDING),
        coloraxis_colorbar=dict(
            len=0.5,
            y=0.5
        )
    )


# Dot plot showing receptor/ligand and sending/receiving pairs, for a selected TF
def domino_dot_plot_tf(domino_df, selected_tf, palette=PALETTE_1, dot_size_scale=1.0, white_background=False):
    if not selected_tf:
        raise PlotSkipped("Select a TF to generate this plot")

    dot_chart_data = build_dot_chart_data(domino_df[domino_df["tf"] == selected_tf])
    if dot_chart_data.empty:
        raise PlotSkipped("No interactions found for TF: {0}".format(selected_tf))

    dot_fig = px.scatter(
        dot_chart_data,
        x="clusters",
        y="lig_rec",
        color="mean_exp",                    # color based on relationship score
        color_continuous_scale=palette["dot_colorscale"],  # heatmap-like gradient
    )

    style_dot_plot(dot_fig, dot_chart_data, selected_tf, tickangle=-90, row_height=30,
                   marker_size=DOT_TF_MARKER_SIZE * dot_size_scale, white_background=white_background)
    return dot_fig


# Dot plot for selected sending/receiving clusters, optionally thresholded to reduce plot size
def domino_dot_plot(domino_df, sending_clusters, receiving_clusters, filter=False,
                    exp_range_threshold=DEFAULT_EXP_RANGE_THRESHOLD, exp_max_threshold=DEFAULT_EXP_MAX_THRESHOLD,
                    palette=PALETTE_1, dot_size_scale=1.0, white_background=False):

    # Select only the rows that involve the sending and/or receiving clusters
    if sending_clusters and receiving_clusters:
        domino_df_selected = domino_df[(domino_df["sending_cl"].isin(sending_clusters)) & (domino_df["receiving_cl"].isin(receiving_clusters))]
    elif sending_clusters:
        domino_df_selected = domino_df[domino_df["sending_cl"].isin(sending_clusters)]
    elif receiving_clusters:
        domino_df_selected = domino_df[domino_df["receiving_cl"].isin(receiving_clusters)]
    else:
        raise PlotSkipped("Select sending and/or receiving clusters to generate this plot")

    dot_chart_data = build_dot_chart_data(domino_df_selected)

    if filter:
        keep = get_lig_rec_pairs_passing_thresholds(dot_chart_data, exp_range_threshold, exp_max_threshold)
        dot_chart_data = dot_chart_data[dot_chart_data["lig_rec"].isin(keep)]

    if dot_chart_data.empty:
        raise PlotSkipped("No interactions found for the selected clusters")

    dot_fig = px.scatter(
        dot_chart_data,
        x="clusters",
        y="lig_rec",
        color="mean_exp",                    # color based on relationship score
        color_continuous_scale=palette["dot_colorscale"],  # heatmap-like gradient
    )

    if filter:
        title = "Geometric Mean of Ligand and Receptor Expression" + thresholds_title_suffix(exp_range_threshold, exp_max_threshold)
        row_height = 30
    else:
        title = "Geometric Mean of Ligand and Receptor Expression"
        row_height = 20

    style_dot_plot(dot_fig, dot_chart_data, title, tickangle=-60, row_height=row_height,
                   marker_size=DOT_CLUSTERS_MARKER_SIZE * dot_size_scale, white_background=white_background)
    return dot_fig


### River plots ###

def domino_river_plot(domino_df, stages, title, palette=PALETTE_1):
    """
    Sankey diagram flowing through the given stages (columns) in order, for all rows of domino_df.
    Each unique pair of nodes in consecutive stages gets one link, sized by the number of interactions.
    """
    if domino_df.empty:
        raise PlotSkipped("No interactions found for the current selection")

    # Nodes for each stage are sorted, and offset by the nodes of the preceding stages
    node_label = []
    node_index = {}
    for stage in stages:
        stage_nodes = sorted(domino_df[stage].unique())
        node_index[stage] = {node: len(node_label) + i for i, node in enumerate(stage_nodes)}
        node_label.extend(stage_nodes)

    source_node = []
    target_node = []
    weight = []
    for from_stage, to_stage in zip(stages, stages[1:]):
        for (from_node, to_node), count in domino_df.groupby([from_stage, to_stage]).size().items():
            source_node.append(node_index[from_stage][from_node])
            target_node.append(node_index[to_stage][to_node])
            weight.append(int(count))

    river_fig = go.Figure(data=[go.Sankey(
        node=dict(
            line=dict(color="black", width=0.1),
            label=node_label,
            color=palette["river_node_color"],
            hovertemplate="%{label}<br>%{value} interactions<extra></extra>",
        ),
        link=dict(
            source=source_node,  # indices correspond to labels
            target=target_node,
            value=weight,
            color=palette["river_link_color"],
            hovercolor=palette["river_link_hover_color"],
            hovertemplate="%{source.label} -> %{target.label}<br>%{value} interactions<extra></extra>",
        ))])

    # Column headers, aligned with the sankey node columns which span the full plot width
    last = len(stages) - 1
    for x_coordinate, stage in enumerate(stages):
        river_fig.add_annotation(
            # Just above the plot area, a fixed distance regardless of plot height (so tall plots don't overlap the title)
            x=x_coordinate / last,
            y=1,
            yanchor="bottom",
            yshift=5,
            xref="paper",
            yref="paper",
            xanchor="left" if x_coordinate == 0 else "right" if x_coordinate == last else "center",
            text=STAGE_NAMES[stage],
            showarrow=False,
            font=dict(
                size=16,
                color="black"
            ),
        )

    max_stage_nodes = max(len(nodes) for nodes in node_index.values())
    river_fig.update_layout(
        title_text=title,
        font_size=10,
        height=min(RIVER_MAX_HEIGHT, max(RIVER_MIN_HEIGHT, RIVER_HEIGHT_PER_NODE * max_stage_nodes)),
        margin=dict(t=100),
        plot_bgcolor='rgba(0,0,0,0)'
    )

    return river_fig


# Links from sending cluster to ligands to receptors to receiving clusters, for a selected TF
def domino_river_plot_tf(domino_df, selected_tf, palette=PALETTE_1):
    if not selected_tf:
        raise PlotSkipped("Select a TF to generate this plot")

    return domino_river_plot(domino_df[domino_df["tf"] == selected_tf],
                             ["sending_cl", "ligand", "receptor", "receiving_cl"], "Selected TF: " + selected_tf, palette)


# Links from sending cluster to receptors to TFs to receiving clusters, for a selected ligand
def domino_river_plot_ligand(domino_df, selected_ligand, palette=PALETTE_1):
    if not selected_ligand:
        raise PlotSkipped("Select a ligand to generate this plot")

    return domino_river_plot(domino_df[domino_df["ligand"] == selected_ligand],
                             ["sending_cl", "receptor", "tf", "receiving_cl"], "Selected ligand: " + selected_ligand, palette)


# Links from sending cluster to ligands to TFs to receiving clusters, for a selected receptor
def domino_river_plot_receptor(domino_df, selected_receptor, palette=PALETTE_1):
    if not selected_receptor:
        raise PlotSkipped("Select a receptor to generate this plot")

    return domino_river_plot(domino_df[domino_df["receptor"] == selected_receptor],
                             ["sending_cl", "ligand", "tf", "receiving_cl"], "Selected receptor: " + selected_receptor, palette)


def filter_by_selections(domino_df, selections):
    """Keep rows where, for every non-empty selection list, the row's value in that column is in the list."""
    mask = pd.Series(True, index=domino_df.index)
    for column, selected in selections.items():
        if selected:
            mask &= domino_df[column].isin(selected)
    return domino_df[mask]


def filter_by_expression_thresholds(domino_df, exp_range_threshold, exp_max_threshold):
    """
    Keep rows whose ligand/receptor pair passes the dot plot expression thresholds.  As in the
    thresholded dot plot, a pair's geometric mean expression range and maximum are computed
    across the sending/receiving cluster pairs within domino_df.
    """
    keep = get_lig_rec_pairs_passing_thresholds(build_dot_chart_data(domino_df), exp_range_threshold, exp_max_threshold)
    return domino_df[(domino_df["ligand"] + '->' + domino_df["receptor"]).isin(keep)]


def domino_river_plot_custom(domino_df, selections, stages, apply_thresholds, exp_range_threshold, exp_max_threshold, palette):
    """
    River plot of the rows matching all non-empty selections, optionally keeping only ligand/receptor
    pairs that pass the dot plot expression thresholds (computed over the selected rows).
    """
    selected_df = filter_by_selections(domino_df, selections)
    title = "Cell-cell communication"

    if apply_thresholds and not selected_df.empty:
        selected_df = filter_by_expression_thresholds(selected_df, exp_range_threshold, exp_max_threshold)
        title += thresholds_title_suffix(exp_range_threshold, exp_max_threshold)

        if selected_df.empty:
            raise PlotSkipped("No ligand/receptor pairs in the current selection pass the dot plot expression thresholds")

    return domino_river_plot(selected_df, stages, title, palette)


# Full pathway, restricted to rows matching all non-empty selections
def domino_river_plot_selected(domino_df, selected_ligands=None, selected_receptors=None, selected_tfs=None,
                               selected_sending=None, selected_receiving=None, apply_thresholds=False,
                               exp_range_threshold=DEFAULT_EXP_RANGE_THRESHOLD, exp_max_threshold=DEFAULT_EXP_MAX_THRESHOLD,
                               palette=PALETTE_1):
    selections = {
        "ligand": selected_ligands,
        "receptor": selected_receptors,
        "tf": selected_tfs,
        "sending_cl": selected_sending,
        "receiving_cl": selected_receiving,
    }

    return domino_river_plot_custom(domino_df, selections, ["sending_cl", "ligand", "receptor", "tf", "receiving_cl"],
                                    apply_thresholds, exp_range_threshold, exp_max_threshold, palette)


# Same as above but with TFs dropped, to simplify plots spanning a wide range of parameters
def domino_river_plot_selected_notf(domino_df, selected_ligands=None, selected_receptors=None,
                                    selected_sending=None, selected_receiving=None, apply_thresholds=False,
                                    exp_range_threshold=DEFAULT_EXP_RANGE_THRESHOLD, exp_max_threshold=DEFAULT_EXP_MAX_THRESHOLD,
                                    palette=PALETTE_1):
    selections = {
        "ligand": selected_ligands,
        "receptor": selected_receptors,
        "sending_cl": selected_sending,
        "receiving_cl": selected_receiving,
    }

    return domino_river_plot_custom(domino_df, selections, ["sending_cl", "ligand", "receptor", "receiving_cl"],
                                    apply_thresholds, exp_range_threshold, exp_max_threshold, palette)
