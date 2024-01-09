
import base64
import io
import os
import re
import sys
from math import ceil
from pathlib import Path
from time import sleep

import geardb
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from flask import request
from flask_restful import Resource
from matplotlib import cm

from werkzeug.utils import secure_filename

from gear.orthology import get_ortholog_file, get_ortholog_files_from_dataset, map_single_gene

sc.settings.set_figure_params(dpi=100)
sc.settings.verbosity = 0

PLOT_TYPE_TO_BASIS = {
    "tsne_static": "tsne",
    "tsne": "tsne",             # legacy
    "umap_static": "umap",
    "pca_static": "pca"
}
COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"

NUM_LEGENDS_PER_COL = 12    # Max number of legend items per column allowed in vertical legend
NUM_HORIZONTAL_COLS = 8 # Number of columns in horizontal legend

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
PROJECTIONS_BASE_DIR = abs_path_www.joinpath('projections')

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id):
    """
    Maps a gene symbol to its corresponding orthologous gene symbol in a given dataset.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        gene_organism_id (str): The organism ID of the gene symbol.
        dataset_organism_id (str): The organism ID of the dataset.

    Returns:
        str: The mapped orthologous gene symbol, or None if no mapping is found.
    """
    if gene_organism_id and gene_organism_id != dataset_organism_id:
        ortholog_file = get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")
        return map_single_gene(gene_symbol, ortholog_file)
    else:
        for ortholog_file in get_ortholog_files_from_dataset(dataset_organism_id, "ensembl"):
            try:
                return map_single_gene(gene_symbol, ortholog_file)
            except:
                continue
    return None

def check_gene_in_dataset(adata, gene_symbols):
    """
    Check if any of the given gene symbols are present in the dataset.

    Args:
        adata (AnnData): Annotated data object.
        gene_symbols (list): List of gene symbols to check.

    Returns:
        bool: True if any of the gene symbols are present in the dataset, False otherwise.
    """
    gene_filter = adata.var.gene_symbol.isin(gene_symbols)
    return gene_filter.any()

def get_analysis(analysis, dataset_id, session_id):
    """Return analysis object based on various factors."""
    # If an analysis is posted we want to read from its h5ad
    if analysis:
        user = geardb.get_user_from_session_id(session_id)
        ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=user.id)

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

def calculate_figure_height(num_plots):
    """Determine height of tsne plot based on number of group elements."""
    return (num_plots * 4) + (num_plots -1)

def calculate_figure_width(num_plots):
    """Determine width of tsne plot based on number of group elements."""
    return (num_plots * 6) + (num_plots -1)

def calculate_num_legend_cols(group_len):
    """Determine number of columns legend should have in tSNE plot."""
    return ceil(group_len / NUM_LEGENDS_PER_COL)

def create_colorscale_with_zero_gray(colorscale):
    """Take a predefined colorscale, and change the 0-value color to gray, and return."""
    from matplotlib import cm

    # Create custom colorscale with gray at the 0.0 level
    # Src: https://matplotlib.org/tutorials/colors/colormap-manipulation.html
    cmap = cm.get_cmap(colorscale, 256)
    newcolors = cmap(np.linspace(0, 1, 256))  # split colormap into 256 parts over 0:1 range
    gray = np.array([192/256, 192/256, 192/256, 1])
    newcolors[0, :] = gray
    return mcolors.ListedColormap(newcolors)

def create_projection_adata(dataset_adata, dataset_id, projection_id):
    # Create AnnData object out of readable CSV file
    # ? Does it make sense to put this in the geardb/Analysis class?
    projection_id = secure_filename(projection_id)
    dataset_id = secure_filename(dataset_id)

    projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
    # Sanitize input to prevent path traversal
    projection_adata_path = projection_dir.joinpath("{}.h5ad".format(projection_id))

    # SAdkins - Run into issues where the h5ad present had different observation columns than the dataset adata, leading to KeyError
    #  when using certain columns. For now, we will always createa new backed h5ad file.
    #if projection_adata_path.is_file():
    #    return sc.read_h5ad(projection_adata_path)#, backed="r")

    projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
    try:
        projection_adata = sc.read_csv(projection_csv_path)
    except Exception as e:
        print(str(e), file=sys.stderr)
        raise PlotError("Could not create projection AnnData object from CSV.")
    projection_adata.obs = dataset_adata.obs
    projection_adata.obsm = dataset_adata.obsm
    # Close dataset adata so that we do not have a stale opened object
    if dataset_adata.isbacked:
        dataset_adata.file.close()
    projection_adata.var["gene_symbol"] = projection_adata.var_names
    # Associate with a filename to ensure AnnData is read in "backed" mode
    projection_adata.filename = projection_adata_path
    return projection_adata

def create_projection_pca_colorscale():
    """Create a diverging colorscale but with black in the middle range."""
    from matplotlib.colors import LinearSegmentedColormap

    # Src: https://matplotlib.org/stable/tutorials/colors/colormap-manipulation.html#directly-creating-a-segmented-colormap-from-a-list
    nodes = [0.0, 0.25, 0.4, 0.5, 0.6, 0.75, 1.0]
    colors = ["lightblue", "blue", "darkblue", "black", "darkred", "red", "lightcoral"]
    return LinearSegmentedColormap.from_list("projection_pca", list(zip(nodes, colors)))


def get_colorblind_scale(n_colors):
    """Get a colorblind friendly colorscale (Viridis). Return n colors spaced equidistantly."""
    cividis = cm.get_cmap("viridis", n_colors)
    colors = cividis.colors
    # convert to hex since I ran into some issues using rpg colors
    return [mcolors.rgb2hex(color) for color in colors]

def get_projection_algorithm(dataset_id, projection_id):
    dataset_projection_json_file = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id, "projections.json")
    # Projection already exists, so we can just return info we want to return in a message
    import json
    with open(dataset_projection_json_file) as projection_fh:
        projections_dict = json.load(projection_fh)
    for genecart_id, configs in projections_dict.items():
        for config in configs:
            if config["uuid"] == projection_id:
                # Handle legacy algorithm
                if "is_pca" in config:
                    config["algorithm"] = "pca" if config["is_pca"] == 1 else "nmf"
                return config["algorithm"]
    raise Exception("Could not find projection ID entry within projection CSV file.")

def is_categorical(series):
    """Return True if Dataframe series is categorical."""
    return series.dtype.name == 'category'

def sort_legend(figure, sort_order, horizontal_legend=False):
    """Sort legend of plot."""
    handles, labels = figure.get_legend_handles_labels()
    new_handles = [handles[idx] for idx, name in enumerate(sort_order)]
    new_labels = [labels[idx] for idx, name in enumerate(sort_order)]

    # If horizontal legend, we need to sort in a way to have labels read from left to right
    if horizontal_legend:
        leftover_cards = len(new_handles) % NUM_HORIZONTAL_COLS
        num_chunks = int(len(new_handles) / NUM_HORIZONTAL_COLS)

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

# Rename axes labels to be whatever x and y fields were passed in
# If axis labels are from an analysis, just use default embedding labels
def rename_axes_labels(ax, x_axis, y_axis):
    if not x_axis.startswith("X_"):
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)

class TSNEData(Resource):
    """Resource for retrieving tsne data from an analysis.

    Returns
    -------
    Byte stream image data
    """
    def post(self, dataset_id):
        req = request.get_json()

        gene_symbol = req.get('gene_symbol', None)
        gene_organism_id = req.get('gene_organism_id', None)
        plot_type = req.get('plot_type', "tsne_static")
        analysis = req.get('analysis', None)
        colorize_by = req.get('colorize_legend_by')
        skip_gene_plot = req.get('skip_gene_plot', False)
        plot_by_group = req.get('plot_by_group', None) # One expression plot per group
        max_columns = req.get('max_columns')   # Max number of columns before plotting to a new row
        colors = req.get('colors')
        order = req.get('order', {})
        x_axis = req.get('x_axis', 'tSNE_1')   # Add here in case old tSNE plotly configs are missing axes data
        y_axis = req.get('y_axis', 'tSNE_2')
        flip_x = req.get('flip_x', False)
        flip_y = req.get('flip_y', False)
        horizontal_legend = req.get('horizontal_legend', False)
        marker_size = req.get("marker_size", None)
        filters = req.get('obs_filters', {})    # dict of lists
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)
        projection_id = req.get('projection_id', None)    # projection id of csv output
        colorblind_mode = req.get('colorblind_mode', False)
        sc.settings.figdir = '/tmp/'

        if not gene_symbol or not dataset_id:
            return {
                "success": -1,
                "message": "Request needs both dataset id and gene symbol."
            }

        try:
            ana = get_analysis(analysis, dataset_id, session_id)
        except PlotError as pe:
            return {
                "success": -1,
                "message": str(pe)
            }

        adata = ana.get_adata(backed=True)

        if projection_id:
            try:
                adata = create_projection_adata(adata, dataset_id, projection_id)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }


        dataset = geardb.get_dataset_by_id(dataset_id)
        dataset_organism_id = dataset.organism_id

        mapped_gene_symbol = None
        gene_symbols = (gene_symbol,)

        if 'gene_symbol' not in adata.var.columns:
            return {"success": -1, "message": "The h5ad is missing the gene_symbol column."}

        if not check_gene_in_dataset(adata, gene_symbols):
            try:
                mapped_gene_symbol = get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id)
            except:
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be mapped to the dataset organism."}

            if mapped_gene_symbol:
                gene_symbols = (mapped_gene_symbol,)
                if not check_gene_in_dataset(adata, gene_symbols):
                    return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be found in the h5ad file."}
            else:
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be mapped to the dataset organism."}

        gene_filter = adata.var.gene_symbol.isin(gene_symbols)

        # Primary dataset - find tSNE_1 and tSNE_2 in obs and build X_tsne
        if analysis is None or analysis in ["null", "undefined", dataset_id]:
            # Valid analysis column names from api/resources/h5ad.py
            analysis_tsne_columns = ['X_tsne_1', 'X_tsne_2']
            analysis_umap_columns = ['X_umap_1', 'X_umap_2']
            analysis_pca_columns = ['X_pca_1', 'X_pca_2']
            # Safety check to ensure analysis was populated in adata.obsm if user selects those data series
            if x_axis in analysis_tsne_columns and y_axis in analysis_tsne_columns:
                if not 'X_tsne' in adata.obsm:
                    return {
                        'success': -1,
                        'message': 'Analysis tSNE columns were selected but values not present in adata.obsm'
                    }
            elif x_axis in analysis_umap_columns and y_axis in analysis_umap_columns:
                if not 'X_umap' in adata.obsm:
                    return {
                        'success': -1,
                        'message': 'Analysis UMAP was selected but values not present in adata.obsm'
                    }
            elif x_axis in analysis_pca_columns and y_axis in analysis_pca_columns:
                if not 'X_pca' in adata.obsm:
                    return {
                        'success': -1,
                        'message': 'Analysis PCA was selected but values not present in adata.obsm'
                    }
            else:
                # If using user-provided columns, perform safety check and create adata.obsm entry
                for ds in [x_axis, y_axis]:
                    if ds not in adata.obs:
                        return {
                            'success': -1,
                            'message': 'Dataseries {} was selected but not present in adata.obs'.format(ds)
                        }
                adata.obsm['X_tsne'] = adata.obs[[x_axis, y_axis]].values
                adata.obsm['X_umap'] = adata.obs[[x_axis, y_axis]].values
                adata.obsm['X_pca'] = adata.obs[[x_axis, y_axis]].values

        # Flip x or y axis if requested
        for key in ['X_tsne', 'X_umap', 'X_pca']:
            if key in adata.obsm:
                if flip_x:
                    adata.obsm[key][:,0] = -1 * adata.obsm[key][:,0]
                if flip_y:
                    adata.obsm[key][:,1] = -1 * adata.obsm[key][:,1]

        # We also need to change the adata's Raw var dataframe
        # We can't explicitly reset its index so we reinitialize it with
        # the newer adata object.
        # https://github.com/theislab/anndata/blob/master/anndata/base.py#L1020-L1022
        if adata.raw is not None:
            adata.raw = adata

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
                except:
                    pass

        # Filter genes and slice the adata to get a dataframe
        # with expression and its observation metadata
        try:
            selected = adata[:, gene_filter].to_memory()
        except:
            # The "try" may fail for projections as it is already in memory
            selected = adata[:, gene_filter]

        # Filter by obs filters
        if filters:
            for col, values in filters.items():
                selected_filter = selected.obs[col].isin(values)
                selected = selected[selected_filter, :]

        # Close adata so that we do not have a stale opened object
        if adata.isbacked:
            adata.file.close()

        # If selected name in adata.var is also an observation column append _orig to the column name
        selected_gene = gene_symbols[0]
        if selected_gene in selected.obs.columns:
            selected.obs["{}_orig".format(selected_gene)] = selected.obs[selected_gene]
            # delete the original column
            selected.obs.drop(selected_gene, axis=1, inplace=True)

        df = selected.to_df()
        success = 1
        message = ""
        if len(df.columns) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)

        # Drop duplicate gene symbols so that only 1 ensemble ID is used in scanpy
        selected.var = selected.var.reset_index().set_index('gene_symbol')
        # Currently the ensembl_id column is still called 'index', which could be confusing when looking at the new .index
        # Rename to end the confusion
        selected.var = selected.var.rename(columns={selected.var.columns[0]: "ensembl_id"})
        # Modify the AnnData object to not include any duplicated gene symbols (keep only first entry)
        if len(df.columns) > 1:
            scanpy_copy = ana.dataset_path().replace('.h5ad', '.scanpy_dups_removed.h5ad')
            if os.path.exists(scanpy_copy):
                os.remove(scanpy_copy)
            selected = selected[:, selected.var.index.duplicated() == False].copy(filename=scanpy_copy)

        io_fig = None
        try:
            basis = PLOT_TYPE_TO_BASIS[plot_type]
        except:
            raise("{} was not a valid plot type".format(plot_type))

        # NOTE: This may change in the future if users want plots by group w/o the colorize_by plot added
        if plot_by_group:
            skip_gene_plot = None

        if marker_size:
            marker_size = int(marker_size)

        # Reverse cividis so "light" is at 0 and 'dark' is at incresing expression
        expression_color = create_colorscale_with_zero_gray("cividis_r" if colorblind_mode else "YlOrRd")

        # In projections, reorder so that the strongest weights (positive or negative) appear in forefront of plot
        plot_sort_order = True
        plot_vcenter = None
        if projection_id:
            try:
                algo = get_projection_algorithm(dataset_id, projection_id)
                if algo == "pca":
                    median = np.median(selected[:, gene_symbol].X.squeeze())
                    sort_order = np.argsort(np.abs(median - selected[:, gene_symbol].X.squeeze()))
                    ordered_obs = selected.obs.iloc[sort_order].index
                    selected = selected[ordered_obs, :]
                    plot_sort_order = False # scanpy auto-sorts by highest value by default so we need to override that
                    plot_vcenter = median
                    expression_color = "cividis_r" if colorblind_mode else create_projection_pca_colorscale()
            except Exception as e:
                print(str(e), file=sys.stderr)

        # If colorize_by is passed we need to generate that image first, before the index is reset
        #  for gene symbols, then merge them.
        if colorize_by:
            color_category = True if is_categorical(selected.obs[colorize_by]) else False

            if color_category:
                # were custom colors passed?  the color index is the 'colorize_by' label but with '_colors' appended
                color_idx_name = "{0}_colors".format(colorize_by)

                ## why 2?  Handles the cases of a stringified "{}" or actual keyed JSON
                if colors is not None and len(colors) > 2:
                    selected.uns[color_idx_name] = [colors[idx] for idx in selected.obs[colorize_by].cat.categories]

                elif color_idx_name in selected.obs:
                    # Alternative method.  Associate with hexcodes already stored in the dataframe
                    # Making the assumption that these values are hexcodes
                    grouped = selected.obs.groupby([colorize_by, color_idx_name])
                    # Ensure one-to-one mapping between category and hexcodes
                    if len(selected.obs[colorize_by].unique()) == len(grouped):
                        # Test if names are color hexcodes and use those if applicable (if first is good, assume all are)
                        color_hex = selected.obs[color_idx_name].unique().tolist()
                        if re.search(COLOR_HEX_PTRN, color_hex[0]):
                            color_map = {name[0]:name[1] for name, group in grouped}
                            selected.uns[color_idx_name] = [color_map[k] for k in selected.obs[colorize_by].cat.categories]

                if colorblind_mode:
                    # build a cividis color map for the colorblind mode
                    cb_colors = get_colorblind_scale(len(selected.obs[colorize_by].unique()))
                    color_map = {name:cb_colors[idx] for idx, name in enumerate(selected.obs[colorize_by].cat.categories)}
                    selected.uns[color_idx_name] = [color_map[k] for k in selected.obs[colorize_by].cat.categories]

                # Calculate the number of columns in the legend (if applicable)
                num_cols = calculate_num_legend_cols(len(selected.obs[colorize_by].unique()))

                # Get for legend order.
                colorize_by_order = selected.obs[colorize_by].unique()

            """
            NOTE: Quick note about legend "loc" and "bbox_to_anchor" attributes:

            bbox_to_anchor is the location of the legend relative to the plot frame.
            If x and y are 0, that is the lower-left corner of the plot.
            If bbox_to_anchor has 4 options, they are x, y, width, and height.  The last two are ratios relative to the plot. And x and y are the lower corner of the bounding box

            loc is the portion of the legend that will be at the bbox_to_anchor point.
            So, if x=0, y=0, and loc = "lower_left", the lower left corner of the legend will be anchored to the lower left corner of the plot
            """

            # If plotting by group the plot dimensions need to be determined
            if plot_by_group:
                column_order = selected.obs[plot_by_group].unique()
                group_len = len(column_order)
                num_plots = group_len + 2

                max_cols = num_plots
                if max_columns:
                    max_cols = min(int(max_columns), num_plots)
                max_rows = ceil((num_plots) / max_cols)

                # Set up the figure specs
                figwidth = calculate_figure_width(max_cols)
                figheight = calculate_figure_height(max_rows)
                io_fig = plt.figure(figsize=(figwidth,figheight))
                spec = io_fig.add_gridspec(ncols=max_cols, nrows=max_rows)

                selected.obs["gene_expression"] = [float(x) for x in selected[:,selected.var.index.isin([gene_symbol])].X]
                max_expression = max(selected.obs["gene_expression"].tolist())

                row_counter = 0
                col_counter = 0

                # Filter expression data by "plot_by_group" group and plot each instance
                if order and plot_by_group in order:
                    column_order = order[plot_by_group]

                for _,name in enumerate(column_order):
                    # Copy gene expression dataseries to observation
                    # Filter only expression values for a particular group.
                    selected.obs["split_by_group"] = selected.obs.apply(lambda row: row["gene_expression"] if row[plot_by_group] == name else 0, axis=1)
                    f = io_fig.add_subplot(spec[row_counter, col_counter])
                    sc.pl.embedding(selected, basis=basis, color=["split_by_group"], color_map=expression_color, ax=f, show=False, use_raw=False, title=name, vmax=max_expression, size=marker_size, sort_order=plot_sort_order, vcenter=plot_vcenter)
                    rename_axes_labels(f, x_axis, y_axis)
                    col_counter += 1
                    # Increment row_counter when the previous row is filled.
                    if col_counter % max_cols == 0:
                        row_counter += 1
                        col_counter = 0
                # Add total gene plot and color plots
                if not skip_gene_plot:
                    f_gene = io_fig.add_subplot(spec[row_counter, col_counter])    # final plot with colorize-by group
                    sc.pl.embedding(selected, basis=basis, color=[gene_symbol], color_map=expression_color, ax=f_gene, show=False, use_raw=False, size=marker_size, sort_order=plot_sort_order, vcenter=plot_vcenter) # Max expression is vmax by default
                    rename_axes_labels(f_gene, x_axis, y_axis)
                    col_counter += 1
                    # Increment row_counter when the previous row is filled.
                    if col_counter % max_cols == 0:
                        row_counter += 1
                        col_counter = 0
                f_color = io_fig.add_subplot(spec[row_counter, col_counter])    # final plot with colorize-by group
                sc.pl.embedding(selected, basis=basis, color=[colorize_by], ax=f_color, show=False, use_raw=False, size=marker_size)
                rename_axes_labels(f_color, x_axis, y_axis)
                if color_category:
                    (handles, labels) = sort_legend(f_color, colorize_by_order, horizontal_legend)
                    f_color.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
                    if horizontal_legend:
                            io_fig.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                            f_color.get_legend().remove()  # Remove legend added by scanpy

            else:
                # If 'skip_gene_plot' is set, only the colorize_by plot is printed, otherwise print gene symbol and colorize_by plots
                if skip_gene_plot:
                    # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
                    io_fig = plt.figure(figsize=(6, 4))
                    if color_category and len(selected.obs[colorize_by].cat.categories) > 10:
                        io_fig = plt.figure(figsize=(13, 4))
                    spec = io_fig.add_gridspec(ncols=1, nrows=1)
                    f1 = io_fig.add_subplot(spec[0,0])
                    sc.pl.embedding(selected, basis=basis, color=[colorize_by], ax=f1, show=False, use_raw=False, size=marker_size)
                    rename_axes_labels(f1, x_axis, y_axis)
                    if color_category:
                        (handles, labels) = sort_legend(f1, colorize_by_order, horizontal_legend)
                        f1.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
                        if horizontal_legend:
                            io_fig.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                            f1.get_legend().remove()  # Remove legend added by scanpy

                else:
                    # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
                    io_fig = plt.figure(figsize=(13, 4))
                    spec = io_fig.add_gridspec(ncols=2, nrows=1, width_ratios=[1.1, 1])
                    f1 = io_fig.add_subplot(spec[0,0])
                    f2 = io_fig.add_subplot(spec[0,1])
                    sc.pl.embedding(selected, basis=basis, color=[gene_symbol], color_map=expression_color, ax=f1, show=False, use_raw=False, size=marker_size, sort_order=plot_sort_order, vcenter=plot_vcenter)
                    # BUG: the line below throws error with stacktrace
                    # ValueError: To copy an AnnData object in backed mode, pass a filename: `.copy(filename='myfilename.h5ad')`. To load the object into memory, use `.to_memory()
                    sc.pl.embedding(selected, basis=basis, color=[colorize_by], ax=f2, show=False, use_raw=False, size=marker_size)
                    rename_axes_labels(f1, x_axis, y_axis)
                    rename_axes_labels(f2, x_axis, y_axis)
                    if color_category:
                        (handles, labels) = sort_legend(f2, colorize_by_order, horizontal_legend)
                        f2.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
                        if horizontal_legend:
                            io_fig.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                            f2.get_legend().remove()  # Remove legend added by scanpy

        else:
            io_fig = sc.pl.embedding(selected, basis=basis, color=[gene_symbol], color_map=expression_color, return_fig=True, use_raw=False, size=marker_size, sort_order=plot_sort_order, vcenter=plot_vcenter)
            rename_axes_labels(io_fig.axes[0], x_axis, y_axis)

        # Close adata so that we do not have a stale opened object
        if selected.isbacked:
            selected.file.close()

        io_pic = io.BytesIO()
        io_fig.tight_layout()   # This crops out much of the whitespace around the plot. The next line does this with the legend too
        io_fig.savefig(io_pic, format='png', bbox_inches="tight")
        io_pic.seek(0)
        plt.close() # Prevent zombie plots, which can cause issues

        return {
            "success": success,
            "message": message,
            "mapped_gene_symbol": mapped_gene_symbol,
            "image": base64.b64encode(io_pic.read()).decode("utf-8")
        }
