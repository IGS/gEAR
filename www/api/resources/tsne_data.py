
import base64
import io
import os
import re
from math import ceil

import geardb
import matplotlib as mpl

mpl.use("Agg")  # Prevents the need for a display when plotting, also thread-safe
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from flask import request
from flask_restful import Resource, reqparse
from gear.plotting import PlotError

from .common import create_projection_adata

sc.settings.verbosity = 0

PLOT_TYPE_TO_BASIS = {
    "tsne_static": "tsne",
    "tsne": "tsne",             # legacy
    "umap_static": "umap",
    "pca_static": "pca",
    "mg_tsne_static": "tsne",
    "mg_umap_static": "umap",
    "mg_pca_static": "pca",
}
COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"

NUM_LEGENDS_PER_COL = 12    # Max number of legend items per column allowed in vertical legend

parser = reqparse.RequestParser(bundle_errors=True)

# Common for both MGTSNEData and TSNEData
# NOTE: By default the parser assumes the location is coming from "values", which breaks lists, so we set the location to "json" explicitly for that.
parser.add_argument('plot_type', type=str, default="tsne_static")
parser.add_argument('analysis', type=str, default=None)
parser.add_argument('colorize_legend_by', type=str, default=None)
parser.add_argument('max_columns', type=int, default=None)   # Max number of columns before plotting to a new row
parser.add_argument('expression_palette', type=str, default="YlOrRd")
parser.add_argument('reverse_palette', type=bool, default=False)
parser.add_argument('colors', type=dict, default={})
parser.add_argument('order', type=dict, default={})
parser.add_argument('x_axis', type=str, default='tSNE_1')   # Add here in case old tSNE plotly configs are missing axes data
parser.add_argument('y_axis', type=str, default='tSNE_2')
parser.add_argument('flip_x', type=bool, default=False)
parser.add_argument('flip_y', type=bool, default=False)
parser.add_argument('horizontal_legend', type=bool, default=False)
parser.add_argument('marker_size', type=int, default=None)
parser.add_argument('center_around_median', type=bool, default=False)
parser.add_argument('obs_filters', type=dict, default={})    # dict of lists
parser.add_argument('projection_id', type=str, default=None)    # projection id of csv output
parser.add_argument('colorblind_mode', type=bool, default=False)
parser.add_argument('high_dpi', type=bool, default=False)
parser.add_argument('grid_spec', type=str, default="1/1/2/2")  # start_row/start_col/end_row/end_col (end not inclusive)

single_gene_parser = parser.copy()
single_gene_parser.add_argument('gene_symbol', type=str, default=None)
single_gene_parser.add_argument("plot_by_group", type=str, default=None)  # If true, plot by group
single_gene_parser.add_argument('skip_gene_plot', type=bool, default=False)  # If true, skip the gene expression plot
single_gene_parser.add_argument('two_way_palette', type=bool, default=False)  # If true, data extremes are in the forefront

multi_gene_parser = parser.copy()
multi_gene_parser.add_argument('gene_symbols', type=list, default=[], location="json")

def calculate_figure_height(num_plots, span=1):
    """Determine height of tsne plot based on number of group elements."""
    return ((num_plots * 4) * span) + (num_plots - 1)

def calculate_figure_width(num_plots, span=1):
    """Determine width of tsne plot based on number of group elements."""
    # The + (num_plots - 1) is to account for the space between plots
    return ((num_plots * 2) * span) + (num_plots - 1)

def calculate_num_legend_cols(group_len):
    """Determine number of columns legend should have in tSNE plot."""
    return ceil(group_len / NUM_LEGENDS_PER_COL)

def create_bluered_colorscale():
    """Create the continuous plotly bluered colorscale."""
    bluered = mcolors.LinearSegmentedColormap.from_list("bluered", [(0, "blue"), (1, "red")])
    bluered_r = bluered.reversed()

    # register the color map with matplotlib
    register_colormap("bluered", bluered)
    register_colormap("bluered_r", bluered_r)

def create_colorscale_with_zero_gray(colorscale):
    """Take a predefined colorscale, and change the 0-value color to gray, and return."""

    # Create custom colorscale with gray at the 0.0 level
    # Src: https://matplotlib.org/tutorials/colors/colormap-manipulation.html
    cmap = plt.get_cmap(colorscale, 256)
    newcolors = cmap(np.linspace(0, 1, 256))  # split colormap into 256 parts over 0:1 range
    gray = np.array([192/256, 192/256, 192/256, 1])
    newcolors[0, :] = gray
    return mcolors.ListedColormap(newcolors)

def create_bublrd_colorscale():
    """Create a diverging colorscale but with black in the middle range."""
    nodes = [0.0, 0.25, 0.4, 0.5, 0.6, 0.75, 1.0]
    colors = ["lightblue", "blue", "darkblue", "black", "darkred", "red", "lightcoral"]
    register_colormap("bublrd", mcolors.LinearSegmentedColormap.from_list("bublrd", list(zip(nodes, colors))))
    register_colormap("bublrd_r", mcolors.LinearSegmentedColormap.from_list("bublrd_r", list(zip(nodes, colors[::-1]))))

def create_projection_colorscale():
    """Create a diverging colorscale but with black in the middle range."""
    # Src: https://matplotlib.org/stable/tutorials/colors/colormap-manipulation.html#directly-creating-a-segmented-colormap-from-a-list
    nodes = [0.0, 0.12, 0.25, 0.38, 0.5, 0.62, 0.75, 0.88, 1.0]
    colors = ["violet", "blue", "indigo", "darkblue", "black", "darkred", "red", "orange", "yellow"]
    #create a colormap with the name "multicolor_diverging"
    register_colormap("multicolor_diverging", mcolors.LinearSegmentedColormap.from_list("multicolor_diverging", list(zip(nodes, colors))))
    #create a colormap with the name "multicolor_diverging_r"
    register_colormap("multicolor_diverging_r", mcolors.LinearSegmentedColormap.from_list("multicolor_diverging_r", list(zip(nodes, colors[::-1]))))

def create_two_way_sorting(adata, gene_symbol):
    """
    Sorts the data in `adata` based on the absolute difference between the median value of `gene_symbol` and the values in `adata`.

    Parameters:
        adata (AnnData): Annotated data matrix.
        gene_symbol (str): Symbol of the gene to sort by.

    Returns:
        AnnData: Sorted data matrix.
    """

    median = np.median(adata[:, gene_symbol].X.squeeze())
    sort_order = np.argsort(np.abs(median - adata[:, gene_symbol].X.squeeze()))
    ordered_obs = adata.obs.iloc[sort_order].index
    return adata[ordered_obs, :]

def get_colorblind_scale(n_colors):
    """Get a colorblind friendly colorscale (Viridis). Return n colors spaced equidistantly."""
    cividis = plt.get_cmap("viridis", n_colors)
    colors = cividis.colors
    # convert to hex since I ran into some issues using rpg colors
    return [mcolors.rgb2hex(color) for color in colors]

def is_categorical(series):
    """Return True if Dataframe series is categorical."""
    return series.dtype.name == 'category'

def sort_legend(figure, sort_order, num_cols, horizontal_legend=False):
    """Sort legend of plot."""
    handles, labels = figure.get_legend_handles_labels()
    new_handles = [handles[idx] for idx, name in enumerate(sort_order)]
    new_labels = [labels[idx] for idx, name in enumerate(sort_order)]

    # If horizontal legend, we need to sort in a way to have labels read from left to right
    if horizontal_legend:
        leftover_cards = len(new_handles) % num_cols
        num_chunks = int(len(new_handles) / num_cols)

        # If number of groups is less than num_cols, they can just be put on a single line
        if num_chunks == 0:
            return (new_handles, new_labels)

        # Split into relatively equal chumks.
        handles_sublists = np.array_split(np.array(new_handles), num_chunks)
        labels_sublists = np.array_split(np.array(new_labels), num_chunks)

        # Zipping gets weird if there's a remainder so remove those leftover items and add back later
        if leftover_cards:
            handles_leftover = new_handles[-leftover_cards:]
            labels_leftover = new_labels[-leftover_cards:]

            handles_sublists = np.array_split(np.array(new_handles[:-leftover_cards]), num_chunks)
            labels_sublists = np.array_split(np.array(new_labels[:-leftover_cards]), num_chunks)

        # Zip numpy arrays, then flatten into a 1D list.  Add leftover cards as well to end
        new_handles = np.column_stack(handles_sublists).flatten().tolist()
        new_labels = np.column_stack(labels_sublists).flatten().tolist()

        # Insert leftover cards back into the list. Start from back to front so indexes are not manipulated.
        for i in reversed(range(leftover_cards)):
            new_handles.insert(num_chunks * (i+1), handles_leftover[i])
            new_labels.insert(num_chunks * (i+1), labels_leftover[i])

    return (new_handles, new_labels)

def register_colormap(name, cmap):
    """Handle changes to matplotlib colormap interface in 3.6."""
    try:
        if name not in mpl.colormaps:
            mpl.colormaps.register(cmap, name=name)
    except AttributeError:
        mpl.cm.register_cmap(name, cmap)

# Rename axes labels to be whatever x and y fields were passed in
# If axis labels are from an analysis, just use default embedding labels
def rename_axes_labels(ax, x_axis, y_axis):
    if not x_axis.startswith("X_"):
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)

def validate_args(dataset_id, gene_symbol, analysis, session_id, projection_id):
    """
    Validates and processes arguments required for t-SNE data analysis.

    Parameters:
        max_columns (int or str): The maximum number of columns to use; will be converted to int if provided.
        dataset_id (str): Identifier for the dataset. Required.
        gene_symbol (str or list): Gene symbol or symbols to analyze. Required.
        analysis (str): Analysis identifier or object.
        session_id (str): Session identifier.
        projection_id (str or None): Optional projection identifier.

    Returns:
        dict: A dictionary containing:
            - "success" (int): 1 if validation is successful, -1 otherwise.
            - "message" (str): Description of the validation result or error.
            - "adata" (AnnData, optional): The AnnData object if validation is successful.
            - "ana" (object, optional): The analysis object if validation is successful.

    Raises:
        Does not raise exceptions; all errors are caught and returned in the result dictionary.
    """

    if not dataset_id:
        return {
            "success": -1,
            "message": "Request needs a dataset id."
        }

    if not gene_symbol:
        return {
            "success": -1,
            "message": "Request needs a gene symbol or symbols."
        }

    try:
        ana = geardb.get_analysis(analysis, dataset_id, session_id)
    except Exception:
        import traceback
        traceback.print_exc()
        return {
            "success": -1,
            "message": "Could not retrieve analysis."
        }

    try:
        adata = ana.get_adata(backed=True)
    except Exception:
        import traceback
        traceback.print_exc()
        return {
            "success": -1,
            "message": "Could not retrieve AnnData object."
        }

    if projection_id:
        try:
            adata = create_projection_adata(adata, dataset_id, projection_id)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

    if 'gene_symbol' not in adata.var.columns:
        return {
            "success": -1,
            "message": "The h5ad is missing the gene_symbol column."
        }

    return {
        "success": 1,
        "message": "Validation successful.",
        "adata": adata,
        "ana": ana
    }

def filter_adata_genes(adata, gene_symbols: list):
    """
    Filters the input AnnData object to retain only the specified gene symbols and returns a new AnnData object
    containing the filtered genes and associated observation metadata.

    Parameters:
        adata (AnnData): The input AnnData object containing gene expression data.
        gene_symbols (list): List of gene symbols to filter for.

    Returns:
        AnnData: A new AnnData object containing only the specified genes and their expression data.

    Raises:
        ValueError: If none of the specified gene symbols are found in the dataset.

    Notes:
        - If the input AnnData object is backed, its file handle will be closed after filtering.
        - The returned AnnData object's .X matrix is converted to a dense matrix if it was originally sparse.
    """
    # Filter genes and slice the adata to get a dataframe
    # with expression and its observation metadata
    gene_filter = adata.var.gene_symbol.isin(gene_symbols)
    if not gene_filter.any():
        message = 'The searched gene symbols could not be found in the dataset.'
        if len(gene_symbols) == 1:
            message = 'The searched gene symbol "{}" could not be found in the dataset.'.format(gene_symbols[0])
        raise ValueError(message)

    # Filter genes and slice the adata to get a dataframe
    # with expression and its observation metadata
    try:
        selected = adata[:, gene_filter].to_memory()
    except ValueError:
        # The "try" may fail for projections as it is already in memory
        selected = adata[:, gene_filter]

    # Close adata so that we do not have a stale opened object
    if adata.isbacked:
        adata.file.close()

    # convert adata.X to a dense matrix if it is sparse
    # This prevents potential downstream issues
    try:
        selected.X = selected.X.todense()
    except AttributeError:
        pass

    return selected

def validate_analysis_inputs(adata, analysis, dataset_id, x_axis, y_axis):
    """
    Validates and prepares analysis input data for dimensionality reduction visualizations.

    This function checks that the specified x_axis and y_axis columns are present in the provided AnnData object (`adata`)
    and that the required dimensionality reduction results (tSNE, UMAP, PCA) are available in `adata.obsm` as needed.
    If user-provided columns are selected instead of standard analysis columns, it verifies their presence in `adata.obs`
    and populates the corresponding entries in `adata.obsm`.

    Parameters:
        adata (AnnData): The annotated data matrix.
        analysis (str or None): The analysis identifier or None if not specified.
        dataset_id (str): The identifier for the primary dataset.
        x_axis (str): The name of the column to use for the x-axis.
        y_axis (str): The name of the column to use for the y-axis.

    Returns:
        AnnData: The validated and possibly updated AnnData object.

    Raises:
        ValueError: If the required analysis columns are selected but not present in `adata.obsm`,
                    or if user-provided columns are not present in `adata.obs`.
    """
    # Primary dataset - find tSNE_1 and tSNE_2 in obs and build X_tsne
    if analysis is None or analysis in ["null", "undefined", dataset_id]:
        # Valid analysis column names from api/resources/h5ad.py
        analysis_tsne_columns = ['X_tsne_1', 'X_tsne_2']
        analysis_umap_columns = ['X_umap_1', 'X_umap_2']
        analysis_pca_columns = ['X_pca_1', 'X_pca_2']
        # Safety check to ensure analysis was populated in adata.obsm if user selects those data series
        if x_axis in analysis_tsne_columns and y_axis in analysis_tsne_columns:
            if 'X_tsne' not in adata.obsm:
                raise ValueError(
                    'Analysis tSNE columns were selected but values not present in adata.obsm'
                )
        elif x_axis in analysis_umap_columns and y_axis in analysis_umap_columns:
            if 'X_umap' not in adata.obsm:
                raise ValueError(
                    'Analysis UMAP was selected but values not present in adata.obsm'
                )
        elif x_axis in analysis_pca_columns and y_axis in analysis_pca_columns:
            if 'X_pca' not in adata.obsm:
                raise ValueError(
                    'Analysis PCA was selected but values not present in adata.obsm'
                )
        else:
            # If using user-provided columns, perform safety check and create adata.obsm entry
            for ds in [x_axis, y_axis]:
                if ds not in adata.obs:
                    raise ValueError(
                        'Dataseries {} was selected but not present in adata.obs'.format(ds)
                    )
            adata.obsm['X_tsne'] = adata.obs[[x_axis, y_axis]].values
            adata.obsm['X_umap'] = adata.obs[[x_axis, y_axis]].values
            adata.obsm['X_pca'] = adata.obs[[x_axis, y_axis]].values

    return adata

def prepare_adata(adata, flip_x=False, flip_y=False, order=None, filters=None):
    """
    Prepares an AnnData object by optionally flipping embedding axes, reordering categorical observation columns, and filtering observations.
    Args:
        adata (AnnData): The AnnData object to be prepared.
        flip_x (bool, optional): If True, flip the x-axis (first dimension) of embeddings in `obsm` (e.g., 'X_tsne', 'X_umap', 'X_pca'). Defaults to False.
        flip_y (bool, optional): If True, flip the y-axis (second dimension) of embeddings in `obsm`. Defaults to False.
        order (dict, optional): Dictionary mapping observation column names to a list specifying the desired order of categories. Only applies to categorical columns. Defaults to None.
        filters (dict, optional): Dictionary mapping observation column names to a list of values to retain. Only observations matching these values are kept. Defaults to None.
    Returns:
        AnnData: The modified AnnData object with updated embeddings, reordered categories, and/or filtered observations.
    """
    # Flip x or y axis if requested
    for key in ['X_tsne', 'X_umap', 'X_pca']:
        if key in adata.obsm:
            if flip_x:
                adata.obsm[key][:,0] = -1 * adata.obsm[key][:,0]
            if flip_y:
                adata.obsm[key][:,1] = -1 * adata.obsm[key][:,1]

    # Reorder the categorical values in the observation dataframe
    # Currently in UI only "plot_by_group" has reordering capabilities
    if order:
        order = order
        obs_keys = order.keys()
        for key in obs_keys:
            col = adata.obs[key]
            try:
                # Some columns might be numeric, therefore
                # we don't want to reorder these
                reordered_col = col.cat.reorder_categories(
                    order[key], ordered=True)
                adata.obs[key] = reordered_col
            except AttributeError:
                pass

    # Filter by obs filters
    if filters:
        for col, values in filters.items():
            selected_filter = adata.obs[col].isin(values)
            adata = adata[selected_filter, :]

    return adata

def modify_genes_found_as_obs_columns(adata, gene_symbols):
    """
    Modifies the AnnData object to ensure that gene symbols found in observation columns are renamed with an '_orig' suffix.

    Parameters:
        adata (AnnData): The annotated data matrix.
        gene_symbols (list): List of gene symbols to check against observation columns.

    Returns:
        AnnData: The modified AnnData object with '_orig' suffix added to observation columns that match gene symbols.
    """
    # If selected name in adata.var is also an observation column append _orig to the column name
    for selected_gene in gene_symbols:
        if selected_gene in adata.obs.columns:
            adata.obs["{}_orig".format(selected_gene)] = adata.obs[selected_gene]
            # delete the original column
            adata = adata.obs.drop(selected_gene, axis=1)
    return adata

def dedup_genes(adata, ana):
    # Drop duplicate gene symbols so that only 1 ensemble ID is used in scanpy
    adata.var = adata.var.reset_index().set_index('gene_symbol')
    # Currently the ensembl_id column is still called 'index', which could be confusing when looking at the new .index
    # Rename to end the confusion
    adata.var = adata.var.rename(columns={adata.var.columns[0]: "ensembl_id"})
    # Modify the AnnData object to not include any duplicated gene symbols (keep only first entry)
    dedup_copy = ana.dataset_path().replace('.h5ad', '.dups_removed.h5ad')

    success = 1
    message = ""

    if (adata.var.index.duplicated(keep="first")).any():
        success = 2
        message = "WARNING: The dataset contains gene symbols that are duplicated across multiple Ensembl IDs. " \
                  "The first Ensembl ID will be used for each gene symbol."

        if os.path.exists(dedup_copy):
            os.remove(dedup_copy)
        adata = adata[:, ~adata.var.index.duplicated()].copy(filename=dedup_copy)
    return {
            'success': success,
            'message': message,
            'adata': adata,
            "dedup_copy": dedup_copy
    }


def generate_tsne_figure( adata, ana, gene_symbols:list, plot_type, analysis, x_axis, y_axis, order, filters,
    flip_x, flip_y, marker_size, colorize_by, colors, colorblind_mode, center_around_median,
    expression_palette, reverse_palette, high_dpi, grid_spec, max_columns, horizontal_legend,
    skip_gene_plot=None, plot_by_group=None, two_way_palette=None
    ):
    """
    Generates a t-SNE (or other embedding) figure for single-cell gene expression data.

    This function processes the input AnnData object, applies various filters and transformations,
    and creates a plot visualizing gene expression or categorical metadata on a 2D embedding (e.g., t-SNE, UMAP).
    The plot can be customized with different color palettes, marker sizes, colorization by metadata, and more.
    The resulting figure is encoded as a base64 image string.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape n_obs × n_vars.
    ana : object
        Analysis object containing metadata such as dataset_id.
    gene_symbols : list
        List of gene symbols (or a single gene symbol) to plot.
    plot_type : str
        Type of embedding plot (e.g., 'tsne', 'umap').
    analysis : dict
        Analysis parameters or configuration.
    x_axis : str
        Name of the variable to use for the x-axis.
    y_axis : str
        Name of the variable to use for the y-axis.
    order : dict or None
        Optional ordering for categorical variables.
    filters : dict or None
        Optional filters to apply to the data.
    flip_x : bool
        Whether to flip the x-axis.
    flip_y : bool
        Whether to flip the y-axis.
    marker_size : int or None
        Size of the markers in the plot.
    colorize_by : str or None
        Name of the categorical variable to colorize by.
    colors : dict or None
        Dictionary mapping category names to colors.
    colorblind_mode : bool
        Whether to use a colorblind-friendly palette.
    center_around_median : bool
        Whether to center the color scale around the median expression.
    expression_palette : str
        Name of the color palette for gene expression.
    reverse_palette : bool
        Whether to reverse the color palette.
    high_dpi : bool
        Whether to generate a high-DPI image.
    grid_spec : str
        Grid specification for the plot layout, as a string (e.g., "0/0/10/10").
    max_columns : int or None
        Maximum number of columns in the plot grid.
    horizontal_legend : bool
        Whether to display the legend horizontally.
    skip_gene_plot : bool or None, optional
        If True, skips plotting the gene expression plot (single-gene mode).
    plot_by_group : str or None, optional
        Name of the group to split plots by.
    two_way_palette : str or None, optional
        Name of a two-way color palette for special sorting.

    Returns
    -------
    dict
        A dictionary with the following keys:
        - 'success': int, 1 if successful, -1 if an error occurred.
        - 'message': str, error or status message.
        - 'image': str, base64-encoded image of the plot (if successful).
    """
    # --- Validation and preparation ---
    try:
        selected = filter_adata_genes(adata, gene_symbols)
        selected = validate_analysis_inputs(selected, analysis, ana.dataset_id, x_axis, y_axis)
    except ValueError as ve:
        return {'success': -1, 'message': str(ve)}

    selected = prepare_adata(selected, flip_x=flip_x, flip_y=flip_y, order=order, filters=filters)
    selected = modify_genes_found_as_obs_columns(selected, list(gene_symbols))

    dedup = dedup_genes(selected, ana)
    selected = dedup['adata']
    success = dedup['success']
    message = dedup['message']
    dedup_copy = dedup['dedup_copy']

    try:
        basis = PLOT_TYPE_TO_BASIS[plot_type]
    except ValueError:
        return {'success': -1, 'message': f"{plot_type} was not a valid plot type"}

    if marker_size:
        marker_size = int(marker_size)

    if reverse_palette:
        expression_palette += "_r"

    plot_sort_order = True
    plot_vcenter = None

    if center_around_median:
        plot_vcenter = np.median(selected[:, list(gene_symbols)].X.squeeze())

    if two_way_palette:
        selected = create_two_way_sorting(selected, gene_symbols[0])
        plot_sort_order = False

    if expression_palette:
        if expression_palette.startswith("bluered"):
            create_bluered_colorscale()
        elif expression_palette.startswith("bublrd"):
            create_bublrd_colorscale()
        elif expression_palette.startswith("multicolor_diverging"):
            create_projection_colorscale()

    expression_color = create_colorscale_with_zero_gray("cividis_r" if colorblind_mode else expression_palette)

    # --- Plot setup ---
    columns = []
    titles = []
    num_plots = len(gene_symbols)

    # convert max_columns to int
    if max_columns:
        max_columns = int(max_columns)

    # Colorize logic (shared, but with some differences for single/multi)
    if colorize_by:
        num_plots += 1
        color_category = is_categorical(selected.obs[colorize_by])
        if color_category:
            color_idx_name = f"{colorize_by}_colors"
            if colors is not None and len(colors) > 2:
                selected.uns[color_idx_name] = [colors[idx] for idx in selected.obs[colorize_by].cat.categories]
            elif color_idx_name in selected.obs:
                grouped = selected.obs.groupby([colorize_by, color_idx_name], observed=False)
                if len(selected.obs[colorize_by].unique()) == len(grouped):
                    color_hex = selected.obs[color_idx_name].unique().tolist()
                    if re.search(COLOR_HEX_PTRN, color_hex[0]):
                        color_map = {name[0]: name[1] for name, group in grouped}
                        selected.uns[color_idx_name] = [color_map[k] for k in selected.obs[colorize_by].cat.categories]
            if colorblind_mode:
                cb_colors = get_colorblind_scale(len(selected.obs[colorize_by].unique()))
                color_map = {name: cb_colors[idx] for idx, name in enumerate(selected.obs[colorize_by].cat.categories)}
                selected.uns[color_idx_name] = [color_map[k] for k in selected.obs[colorize_by].cat.categories]
            num_cols = calculate_num_legend_cols(len(selected.obs[colorize_by].unique()))
            colorize_by_order = selected.obs[colorize_by].unique()
        else:
            colorize_by_order = []

        # Multi-gene: always add all gene symbols, then colorize_by
        # Single-gene: may have skip_gene_plot or plot_by_group
        if plot_by_group:
            column_order = selected.obs[plot_by_group].unique()
            num_plots += len(column_order)
            max_cols = num_plots

            if max_columns:
                max_cols = min(max_columns, num_plots)
            selected.obs["gene_expression"] = [float(x) for x in selected[:, selected.var.index.isin([gene_symbols[0]])].X]
            max_expression = max(selected.obs["gene_expression"].tolist())
            if order and plot_by_group in order:
                column_order = order[plot_by_group]
            for _, name in enumerate(column_order):
                group_name = name + "_split_by_group"
                selected.obs[group_name] = selected.obs.apply(lambda row: row["gene_expression"] if row[plot_by_group] == name else 0, axis=1)
                columns.append(group_name)
                titles.append(name)
            kwargs_ncols = max_cols
            kwargs_vmax = max_expression
        else:
            kwargs_ncols = max_columns or num_plots
            kwargs_vmax = None

        # Only add gene plot if not skipped (single-gene)
        if skip_gene_plot is None or not skip_gene_plot:
            columns.extend(gene_symbols)
            titles.extend(gene_symbols)
        columns.append(colorize_by)
        titles.append(colorize_by)
    else:
        columns.extend(gene_symbols)
        titles.extend(gene_symbols)
        kwargs_ncols = max_columns or num_plots
        kwargs_vmax = None

    # --- Plotting ---
    kwargs = {
        "basis": basis,
        "color": columns,
        "color_map": expression_color,
        "show": False,
        "use_raw": False,
        "title": titles,
        "size": marker_size,
        "sort_order": plot_sort_order,
        "vcenter": plot_vcenter,
        "return_fig": True,
        "ncols": kwargs_ncols,
    }
    if kwargs_vmax is not None:
        kwargs["vmax"] = kwargs_vmax

    io_fig = sc.pl.embedding(selected, **kwargs)
    ax = io_fig.get_axes()

    # Grid/figsize logic (shared)
    grid_spec = [int(x) for x in grid_spec.split('/')]
    row_span = grid_spec[2] - grid_spec[0]
    col_span = ceil((grid_spec[3] - grid_spec[1]) / 3)
    num_plots_wide = kwargs_ncols
    num_plots_high = ceil(len(columns) / num_plots_wide)
    io_fig.set_figwidth(calculate_figure_width(num_plots_wide, col_span))
    io_fig.set_figheight(calculate_figure_height(num_plots_high, row_span))

    # Axes/legend logic (shared)
    if isinstance(ax, list):

        # Rename axes labels for each subplot
        for f in ax:
            if f.get_label() == "<colorbar>":
                continue
            rename_axes_labels(f, x_axis, y_axis)

        last_ax = ax[-1]    # color axes
        if colorize_by and color_category:
            """
            NOTE: Quick note about legend "loc" and "bbox_to_anchor" attributes:

            bbox_to_anchor is the location of the legend relative to the plot frame.
            If x and y are 0, that is the lower-left corner of the plot.
            If bbox_to_anchor has 4 options, they are x, y, width, and height.  The last two are ratios relative to the plot. And x and y are the lower corner of the bounding box

            loc is the portion of the legend that will be at the bbox_to_anchor point.
            So, if x=0, y=0, and loc = "lower_left", the lower left corner of the legend will be anchored to the lower left corner of the plot
            """

            num_horizontal_cols = 2 * num_cols # Number of columns in horizontal legend

            (handles, labels) = sort_legend(last_ax, colorize_by_order, num_horizontal_cols, horizontal_legend)
            last_ax.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
            if horizontal_legend:
                last_ax.get_legend().remove() # Remove legend added by scanpy
                last_ax.legend(loc="upper right", bbox_to_anchor=[1, -0.05, 0, 0], frameon=False, ncol=num_horizontal_cols, handles=handles, labels=labels)
    else:
        rename_axes_labels(ax, x_axis, y_axis)

    # Clean up
    if selected.isbacked:
        selected.file.close()
    if os.path.exists(dedup_copy):
        os.remove(dedup_copy)

    with io.BytesIO() as io_pic:
        if high_dpi:
            dpi = max(150, int(selected.shape[0] / 100))
            sc.settings.set_figure_params(dpi_save=dpi)
            io_fig.set_figwidth(num_plots_wide * 10)
            io_fig.set_figheight(num_plots_high * 10)
            io_fig.savefig(io_pic, format='png', bbox_inches="tight")
        else:
            sc.settings.set_figure_params(dpi_save=150)
            io_fig.savefig(io_pic, format='webp', bbox_inches="tight")
        io_pic.seek(0)
        plt.close()
        image = base64.b64encode(io_pic.read()).decode("utf-8")

    return {
        "success": success,
        "message": message,
        "image": image
    }

# --- Resource classes ---

class MGTSNEData(Resource):
    def post(self, dataset_id):
        session_id = request.cookies.get("gear_session_id", "")
        args = multi_gene_parser.parse_args()

        gene_symbols = args.get('gene_symbols', [])
        validation = validate_args(
            dataset_id, gene_symbols, args.get('analysis'), session_id, args.get('projection_id')
        )
        if validation['success'] != 1:
            return validation
        return generate_tsne_figure(
            validation['adata'], validation['ana'], gene_symbols,
            args.get('plot_type', "tsne_static"),
            args.get('analysis', None),
            args.get('x_axis', 'tSNE_1'),
            args.get('y_axis', 'tSNE_2'),
            args.get('order', {}),
            args.get('obs_filters', {}),
            args.get('flip_x', False),
            args.get('flip_y', False),
            args.get('marker_size', None),
            args.get('colorize_legend_by', None),
            args.get('colors', {}),
            args.get('colorblind_mode', False),
            args.get('center_around_median', False),
            args.get('expression_palette', "YlOrRd"),
            args.get('reverse_palette', False),
            args.get('high_dpi', False),
            args.get('grid_spec', "1/1/2/2"),
            args.get('max_columns', None),
            args.get('horizontal_legend', False),
        )

class TSNEData(Resource):
    def post(self, dataset_id):
        session_id = request.cookies.get("gear_session_id", "")
        args = single_gene_parser.parse_args()
        gene_symbol = args.get('gene_symbol', None)
        validation = validate_args(
            dataset_id, gene_symbol, args.get('analysis'), session_id, args.get('projection_id')
        )
        if validation['success'] != 1:
            return validation
        return generate_tsne_figure(
            validation['adata'], validation['ana'], [gene_symbol],
            args.get('plot_type', "tsne_static"),
            args.get('analysis', None),
            args.get('x_axis', 'tSNE_1'),
            args.get('y_axis', 'tSNE_2'),
            args.get('order', {}),
            args.get('obs_filters', {}),
            args.get('flip_x', False),
            args.get('flip_y', False),
            args.get('marker_size', None),
            args.get('colorize_legend_by', None),
            args.get('colors', {}),
            args.get('colorblind_mode', False),
            args.get('center_around_median', False),
            args.get('expression_palette', "YlOrRd"),
            args.get('reverse_palette', False),
            args.get('high_dpi', False),
            args.get('grid_spec', "1/1/2/2"),
            args.get('max_columns', None),
            args.get('horizontal_legend', False),
            args.get('skip_gene_plot', False),
            args.get('plot_by_group', None),
            args.get('two_way_palette', False)
        )
