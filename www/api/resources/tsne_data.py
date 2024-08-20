
import base64
import io
import os
import re
from math import ceil, log2

import geardb
import matplotlib as mpl
mpl.use("Agg")  # Prevents the need for a display when plotting, also thread-safe
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from flask import request
from flask_restful import Resource

from gear.plotting import PlotError
from .common import create_projection_adata

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

class TSNEData(Resource):
    """Resource for retrieving tsne data from an analysis.

    Returns
    -------
    Byte stream image data
    """
    def post(self, dataset_id):
        req = request.get_json()

        gene_symbol = req.get('gene_symbol', None)
        plot_type = req.get('plot_type', "tsne_static")
        analysis = req.get('analysis', None)
        colorize_by = req.get('colorize_legend_by')
        skip_gene_plot = req.get('skip_gene_plot', False)
        plot_by_group = req.get('plot_by_group', None) # One expression plot per group
        max_columns = req.get('max_columns')   # Max number of columns before plotting to a new row
        expression_palette = req.get('expression_palette', "YlOrRd")
        reverse_palette = req.get('reverse_palette', False)
        two_way_palette = req.get('two_way_palette', False) # If true, data extremes are in the forefront
        colors = req.get('colors')
        order = req.get('order', {})
        x_axis = req.get('x_axis', 'tSNE_1')   # Add here in case old tSNE plotly configs are missing axes data
        y_axis = req.get('y_axis', 'tSNE_2')
        flip_x = req.get('flip_x', False)
        flip_y = req.get('flip_y', False)
        horizontal_legend = req.get('horizontal_legend', False)
        marker_size = req.get("marker_size", None)
        center_around_median = req.get("center_around_median", False)
        filters = req.get('obs_filters', {})    # dict of lists
        session_id = request.cookies.get('gear_session_id')
        projection_id = req.get('projection_id', None)    # projection id of csv output
        colorblind_mode = req.get('colorblind_mode', False)
        high_dpi = req.get('high_dpi', False)
        grid_spec = req.get('grid_spec', "1/1/2/2") # start_row/start_col/end_row/end_col (end not inclusive)
        sc.settings.figdir = '/tmp/'

        # convert max_columns to int
        if max_columns:
            max_columns = int(max_columns)

        if not dataset_id:
            return {
                "success": -1,
                "message": "Request needs a dataset id."
            }

        if not gene_symbol:
            return {
                "success": -1,
                "message": "Request needs a gene symbol."
            }

        try:
            ana = geardb.get_analysis(analysis, dataset_id, session_id)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return {
                "success": -1,
                "message": "Could not retrieve analysis."
            }

        try:
            adata = ana.get_adata(backed=True)
        except Exception as e:
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

        gene_symbols = (gene_symbol,)

        if 'gene_symbol' not in adata.var.columns:
            return {
                "success": -1,
                "message": "The h5ad is missing the gene_symbol column."
            }

        # Filter genes and slice the adata to get a dataframe
        # with expression and its observation metadata
        gene_filter = adata.var.gene_symbol.isin(gene_symbols)
        if not gene_filter.any():
            return {
                'success': -1,
                'message': 'The searched gene symbol could not be found in the dataset.'
            }

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

        success = 1
        message = ""

        # Drop duplicate gene symbols so that only 1 ensemble ID is used in scanpy
        selected.var = selected.var.reset_index().set_index('gene_symbol')
        # Currently the ensembl_id column is still called 'index', which could be confusing when looking at the new .index
        # Rename to end the confusion
        selected.var = selected.var.rename(columns={selected.var.columns[0]: "ensembl_id"})
        # Modify the AnnData object to not include any duplicated gene symbols (keep only first entry)
        if (selected.var.index.duplicated(keep="first") == True).any():
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(selected_gene)

            dedup_copy = ana.dataset_path().replace('.h5ad', '.dups_removed.h5ad')
            if os.path.exists(dedup_copy):
                os.remove(dedup_copy)
            selected = selected[:, selected.var.index.duplicated() == False].copy(filename=dedup_copy)

        io_fig = None
        try:
            basis = PLOT_TYPE_TO_BASIS[plot_type]
        except ValueError:
            raise ValueError("{} was not a valid plot type".format(plot_type))

        # NOTE: This may change in the future if users want plots by group w/o the colorize_by plot added
        if plot_by_group:
            skip_gene_plot = None

        if marker_size:
            marker_size = int(marker_size)

        # Reverse cividis so "light" is at 0 and 'dark' is at incresing expression
        if reverse_palette:
            expression_palette += "_r"

        plot_sort_order = True   # scanpy auto-sorts by highest value by default
        plot_vcenter = None

        if center_around_median:
            plot_vcenter = np.median(selected[:, selected_gene].X.squeeze())

        if two_way_palette:
            selected = create_two_way_sorting(selected, selected_gene)
            plot_sort_order = False

        if expression_palette and expression_palette.startswith("bluered"):
            create_bluered_colorscale()

        if expression_palette.startswith("bublrd"):
            create_bublrd_colorscale()

        if expression_palette.startswith("multicolor_diverging"):
            create_projection_colorscale()

        expression_color = create_colorscale_with_zero_gray("cividis_r" if colorblind_mode else expression_palette)


        # These will be passed into the sc.pl.embedding function
        columns = []
        titles = []

        # Variables that may be set based on various parameters
        # basis=basis, color=columns, color_map=expression_color, show=False, use_raw=False, title=titles, ncols=max_cols, vmax=max_expression, size=marker_size, sort_order=plot_sort_order, vcenter=plot_vcenter, return_fig=True
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
            "return_fig": True
        }

        num_plots = 1   # single plot

        # If colorize_by is passed we need to generate that image first, before the index is reset
        #  for gene symbols, then merge them.
        if colorize_by:
            num_plots = 2   # gene expression and colorize_by plot

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
                    grouped = selected.obs.groupby([colorize_by, color_idx_name], observed=False)
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

            # If plotting by group the plot dimensions need to be determined
            if plot_by_group:
                column_order = selected.obs[plot_by_group].unique()
                group_len = len(column_order)
                num_plots = group_len + 2

                max_cols = num_plots
                if max_columns:
                    max_cols = min(max_columns, num_plots)

                selected.obs["gene_expression"] = [float(x) for x in selected[:,selected.var.index.isin([selected_gene])].X]
                max_expression = max(selected.obs["gene_expression"].tolist())

                # Filter expression data by "plot_by_group" group and plot each instance
                if order and plot_by_group in order:
                    column_order = order[plot_by_group]

                for _,name in enumerate(column_order):
                    # Copy gene expression dataseries to observation
                    # Filter only expression values for a particular group.
                    group_name = name + "_split_by_group"
                    selected.obs[group_name] = selected.obs.apply(lambda row: row["gene_expression"] if row[plot_by_group] == name else 0, axis=1)
                    columns.append(group_name)
                    titles.append(name)

                kwargs["ncols"] = max_cols
                kwargs["vmax"] = max_expression

            # If 'skip_gene_plot' is set, only the colorize_by plot is printed, otherwise print gene symbol and colorize_by plots
            if not skip_gene_plot:
                columns.append(selected_gene)
                titles.append(selected_gene)

            columns.append(colorize_by)
            titles.append(colorize_by)

        else:
            columns.append(selected_gene)
            titles.append(selected_gene)

        io_fig = sc.pl.embedding(selected, **kwargs)
        ax = io_fig.get_axes()

        # break grid_spec into spans
        grid_spec = grid_spec.split('/')
        grid_spec = [int(x) for x in grid_spec]
        row_span = grid_spec[2] - grid_spec[0]
        col_span = ceil((grid_spec[3] - grid_spec[1]) / 3)    # Generally these plots span columns in multiples of 4.

        # Set the figsize (in inches)
        dpi = io_fig.dpi    # default dpi is 100, but will be saved as 150 later on
        # With 2 plots as a default (gene expression and colorize_by), we want to grow the figure size slowly based on the number of plots

        num_plots_wide = max_columns if max_columns else num_plots
        num_plots_high = ceil(num_plots / num_plots_wide)

        # set the figsize based on the number of plots
        io_fig.set_figwidth(calculate_figure_width(num_plots_wide, col_span))
        io_fig.set_figheight(calculate_figure_height(num_plots_high, row_span))

        # rename axes labels
        if type(ax) == list:
            for f in ax:
                # skip colorbar
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

                (handles, labels) = sort_legend(last_ax, colorize_by_order, horizontal_legend)
                last_ax.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
                if horizontal_legend:
                    last_ax.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                    last_ax.get_legend().remove() # Remove legend added by scanpy
        else:
            rename_axes_labels(ax, x_axis, y_axis)

        # Close adata so that we do not have a stale opened object
        if selected.isbacked:
            selected.file.close()

        with io.BytesIO() as io_pic:
            # Set the saved figure dpi based on the number of observations in the dataset after filtering
            if high_dpi:
                dpi = max(150, int(selected.shape[0] / 100))
                sc.settings.set_figure_params(dpi_save=dpi)
                # Double the height and width of the figure to maintain the same size
                io_fig.set_figwidth(num_plots_wide * 10)
                io_fig.set_figheight(num_plots_high * 10)

                io_fig.savefig(io_pic, format='png', bbox_inches="tight")
            else:
                # Moved this to the end to prevent any issues with the dpi setting
                sc.settings.set_figure_params(dpi_save=150)
                io_fig.savefig(io_pic, format='webp', bbox_inches="tight")
            io_pic.seek(0)
            plt.close() # Prevent zombie plots, which can cause issues

            image = base64.b64encode(io_pic.read()).decode("utf-8")

        return {
            "success": success,
            "message": message,
            "image": image
        }
