from flask import request
from flask_restful import Resource

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap

import json, os, re
import sys
import geardb
import base64
import io
from math import ceil

sc.settings.set_figure_params(dpi=100)

PLOT_TYPE_TO_BASIS = {
    "tsne_static": "tsne",
    "tsne": "tsne",             # legacy
    "umap_static": "umap",
    "pca_static": "pca"
}
COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"

NUM_LEGENDS_PER_COL = 12

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
    # Create custom colorscale with gray at the 0.0 level
    # Src: https://matplotlib.org/tutorials/colors/colormap-manipulation.html
    ylorrd = cm.get_cmap(colorscale, 256)
    newcolors = ylorrd(np.linspace(0, 1, 256))  # split colormap into 256 parts over 0:1 range
    gray = np.array([192/256, 192/256, 192/256, 1])
    newcolors[0, :] = gray
    return ListedColormap(newcolors)

def sort_legend(figure, sort_order):
    """Sort legend of plot."""
    handles, labels = figure.get_legend_handles_labels()
    new_handles = [handles[idx] for idx, name in enumerate(sort_order)]
    new_labels = [labels[idx] for idx, name in enumerate(sort_order)]
    return (new_handles, new_labels)

class TSNEData(Resource):
    """Resource for retrieving tsne data from an analysis.

    Returns
    -------
      Byte stream image data
    """
    # This endpoint would be a get request to
    # '/api/plot/<SOME_DATASET_ID>/tsne?gene=<SOME_GENE>&analysis=<SOME_ANALYSIS_ID>
    def get(self, dataset_id):
        gene_symbol = request.args.get('gene')
        plot_type = request.args.get('plot_type', "tsne_static")
        analysis_id = request.args.get('analysis')
        colorize_by = request.args.get('colorize_by')
        skip_gene_plot = request.args.get('skip_gene_plot')
        plot_by_group = request.args.get('plot_by_group') # One expression plot per group
        max_columns = request.args.get('max_columns')   # Max number of columns before plotting to a new row
        colors = request.args.get('colors')
        order = request.args.get('order', {})
        x_axis = request.args.get('x_axis', 'tSNE_1')   # Add here in case old tSNE plotly configs are missing axes data
        y_axis = request.args.get('y_axis', 'tSNE_2')
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)
        analysis_owner_id = request.args.get('analysis_owner_id')
        sc.settings.figdir = '/tmp/'

        if not gene_symbol or not dataset_id:
            return {
                "success": -1,
                "message": "Request needs both dataset id and gene symbol."
            }

        dataset = geardb.get_dataset_by_id(dataset_id)

        if analysis_id and analysis_id not in ["null", "undefined"]:
            # need analysis_type here, but can discover it
            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id,
                                  session_id=session_id,
                                  user_id=analysis_owner_id)
            ana.discover_type()
        else:
            ana = geardb.Analysis(type='primary', dataset_id=dataset_id)

        adata = ana.get_adata(backed=True)

        gene_symbols = (gene_symbol,)
        if 'gene_symbol' in adata.var.columns:
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            if not gene_filter.any():
                return {
                    'success': -1,
                    'message': 'Gene not found',
                }

        else:
            return {
                'success': -1,
                'message': 'Missing gene_symbol in adata.var'
            }

        # Primary dataset - find tSNE_1 and tSNE_2 in obs and build X_tsne
        if analysis_id is None or analysis_id in ["null", "undefined", dataset_id]:
            for ds in [x_axis, y_axis]:
                if ds not in adata.obs:
                    return {
                        'success': -1,
                        'message': 'Dataseries {} was selected but not present in adata.obs'.format(ds)
                    }
            adata.obsm['X_tsne'] = adata.obs[[x_axis, y_axis]].values
            adata.obsm['X_umap'] = adata.obs[[x_axis, y_axis]].values
            adata.obsm['X_pca'] = adata.obs[[x_axis, y_axis]].values

        # We also need to change the adata's Raw var dataframe
        # We can't explicitly reset its index so we reinitialize it with
        # the newer adata object.
        # https://github.com/theislab/anndata/blob/master/anndata/base.py#L1020-L1022
        if adata.raw is not None:
            adata.raw = adata

        # Reorder the categorical values in the observation dataframe
        # Currently in UI only "plot_by_group" has reordering capabilities
        if order:
            # order is passed as JSON str.  In plotly_data.py it is passed as Dict though
            order = json.loads(order)
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

        selected = adata[:, gene_filter]
        df = selected.to_df()
        success = 1
        message = ""
        if len(df.columns) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)

        # Drop duplicate gene symbols so that only 1 ensemble ID is used in scanpy
        adata.var = adata.var.reset_index().set_index('gene_symbol')
        # Currently the ensembl_id column is still called 'index', which could be confusing when looking at the new .index
        # Rename to end the confusion
        adata.var = adata.var.rename(columns={adata.var.columns[0]: "ensembl_id"})
        # Modify the AnnData object to not include any duplicated gene symbols (keep only first entry)
        if len(df.columns) > 1:
            scanpy_copy = ana.dataset_path().replace('.h5ad', '.scanpy_dups_removed.h5ad')
            if os.path.exists(scanpy_copy):
                os.remove(scanpy_copy)
            adata = adata[:, adata.var.index.duplicated() == False].copy(filename=scanpy_copy)

        io_fig = None
        try:
            basis = PLOT_TYPE_TO_BASIS[plot_type]
        except:
            raise("{} was not a valid plot type".format(plot_type))

        # NOTE: This may change in the future if users want plots by group w/o the colorize_by plot added
        if plot_by_group and plot_by_group not in ["null", "undefined"]:
            skip_gene_plot = None

        new_YlOrRd = create_colorscale_with_zero_gray("YlOrRd")

        # If colorize_by is passed we need to generate that image first, before the index is reset
        #  for gene symbols, then merge them.
        if colorize_by is not None and colorize_by != 'null':
            # were custom colors passed?  the color index is the 'colorize_by' label but with '_colors' appended
            color_idx_name = "{0}_colors".format(colorize_by)

            ## why 2?  Handles the cases of a stringified "{}" or actual keyed JSON
            if colors is not None and len(colors) > 2:
                colors = json.loads(colors)
                adata.uns[color_idx_name] = [colors[idx] for idx in adata.obs[colorize_by].cat.categories]

            elif color_idx_name in adata.obs:
                # Alternative method.  Associate with hexcodes already stored in the dataframe
                # Making the assumption that these values are hexcodes
                grouped = adata.obs.groupby([colorize_by, color_idx_name])
                # Ensure one-to-one mapping between category and hexcodes
                if len(adata.obs[colorize_by].unique()) == len(grouped):
                    # Test if names are color hexcodes and use those if applicable (if first is good, assume all are)
                    color_hex = adata.obs[color_idx_name].unique().tolist()
                    if re.search(COLOR_HEX_PTRN, color_hex[0]):
                        color_map = {name[0]:name[1] for name, group in grouped}
                        adata.uns[color_idx_name] = [color_map[k] for k in adata.obs[colorize_by].cat.categories]

            # Calculate the number of columns in the legend (if applicable)
            num_cols = calculate_num_legend_cols(len(adata.obs[colorize_by].unique()))

            # Get for legend order.
            colorize_by_order = adata.obs[colorize_by].unique()

            # If plotting by group the plot dimensions need to be determined
            if plot_by_group and plot_by_group not in ["null", "undefined"]:
                column_order = adata.obs[plot_by_group].unique()
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

                adata.obs["gene_expression"] = [float(x) for x in adata[:,adata.var.index.isin([gene_symbol])].X]
                max_expression = max(adata.obs["gene_expression"].tolist())

                row_counter = 0
                col_counter = 0

                # Filter expression data by "plot_by_group" group and plot each instance
                if order and plot_by_group in order:
                    column_order = order[plot_by_group]

                for _,name in enumerate(column_order):
                    # Copy gene expression dataseries to observation
                    # Filter only expression values for a particular group.
                    adata.obs["split_by_group"] = adata.obs.apply(lambda row: row["gene_expression"] if row[plot_by_group] == name else 0, axis=1)
                    f = io_fig.add_subplot(spec[row_counter, col_counter])
                    sc.pl.embedding(adata, basis=basis, color=["split_by_group"], color_map=new_YlOrRd, ax=f, show=False, use_raw=False, title=name, vmax=max_expression)
                    col_counter += 1
                    # Increment row_counter when the previous row is filled.
                    if col_counter % max_cols == 0:
                        row_counter += 1
                        col_counter = 0
                # Add total gene plot and color plots
                if not skip_gene_plot or skip_gene_plot == "false":
                    f_gene = io_fig.add_subplot(spec[row_counter, col_counter])    # final plot with colorize-by group
                    sc.pl.embedding(adata, basis=basis, color=[gene_symbol], color_map=new_YlOrRd, ax=f_gene, show=False, use_raw=False) # Max expression is vmax by default
                    col_counter += 1
                    # Increment row_counter when the previous row is filled.
                    if col_counter % max_cols == 0:
                        row_counter += 1
                        col_counter = 0
                f_color = io_fig.add_subplot(spec[row_counter, col_counter])    # final plot with colorize-by group
                sc.pl.embedding(adata, basis=basis, color=[colorize_by], ax=f_color, show=False, use_raw=False)
                (handles, labels) = sort_legend(f_color, colorize_by_order)
                f_color.legend(bbox_to_anchor=[1, 1], ncol=num_cols, handles=handles, labels=labels)
            else:
                # If 'skip_gene_plot' is set, only the colorize_by plot is printed, otherwise print gene symbol and colorize_by plots
                if skip_gene_plot == 'true':
                    # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
                    io_fig = plt.figure(figsize=(6, 4))
                    if len(adata.obs[colorize_by].cat.categories) > 10:
                        io_fig = plt.figure(figsize=(13, 4))
                    spec = io_fig.add_gridspec(ncols=1, nrows=1)
                    f1 = io_fig.add_subplot(spec[0,0])
                    sc.pl.embedding(adata, basis=basis, color=[colorize_by], ax=f1, show=False, use_raw=False)
                    (handles, labels) = sort_legend(f1, colorize_by_order)
                    f1.legend(bbox_to_anchor=[1,1], ncol=num_cols, handles=handles, labels=labels)
                else:
                    # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
                    io_fig = plt.figure(figsize=(13, 4))
                    spec = io_fig.add_gridspec(ncols=2, nrows=1, width_ratios=[1.1, 1])
                    f1 = io_fig.add_subplot(spec[0,0])
                    f2 = io_fig.add_subplot(spec[0,1])
                    sc.pl.embedding(adata, basis=basis, color=[gene_symbol], color_map=new_YlOrRd, ax=f1, show=False, use_raw=False)
                    sc.pl.embedding(adata, basis=basis, color=[colorize_by], ax=f2, show=False, use_raw=False)
                    (handles, labels) = sort_legend(f2, colorize_by_order)
                    f2.legend(bbox_to_anchor=[1, 1], ncol=num_cols, handles=handles, labels=labels)

        else:
            io_fig = sc.pl.embedding(adata, basis=basis, color=[gene_symbol], color_map=new_YlOrRd, return_fig=True, use_raw=False)

        io_pic = io.BytesIO()
        io_fig.savefig(io_pic, format='png', bbox_inches="tight")
        io_pic.seek(0)
        plt.close() # Prevent zombie plots, which can cause issues

        return {
            "success": success,
            "message": message,
            "image": base64.b64encode(io_pic.read()).decode("utf-8")
        }







