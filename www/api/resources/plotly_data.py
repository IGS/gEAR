import copy
import json
import numbers
import re
import sys

import geardb
import pandas as pd
import plotly.express.colors as pxc
from flask import request
from flask_restful import Resource
from plotly.utils import PlotlyJSONEncoder

from gear.plotting import PlotError, generate_plot, plotly_color_map
from .common import create_projection_adata, order_by_time_point

COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"

class PlotlyData(Resource):
    """Resource for retrieving data from h5ad to be used to draw charts on UI.
    Parameters
    ----------
    gene_symbol: str
        Gene to search in adata.
    index: list[str]
        List of column names to index on.
    colors: dict
        Dictionary mapping levels to colors.
    order: dict
        Dictionary mapping order levels.
    plot_type: str
        Plot type (bar, violin, scatter or line)
    Returns
    -------
    dict
        Plot data
    """

    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        gene_symbol = req.get('gene_symbol', None)
        plot_type = req.get('plot_type')

        # tsne/umap_dynamic is just a symlink to scatter (support legacy tsne_dynamic)
        if plot_type in [ 'tsne_dynamic', 'tsne/umap_dynamic']:
            plot_type = 'scatter'

        x_axis = req.get('x_axis')
        y_axis = req.get('y_axis', "raw_value")
        z_axis = req.get('z_axis')
        label = req.get('point_label')
        hide_x_labels = req.get('hide_x_labels', False)
        hide_y_labels = req.get('hide_y_labels', False)
        hide_legend = req.get('hide_legend', False)
        color_name = req.get('color_name')
        color_map = req.get('colors')
        palette = req.get('color_palette')
        reverse_palette = req.get('reverse_palette')
        facet_row = req.get('facet_row')
        facet_col = req.get('facet_col')
        order = req.get('order', {})
        analysis = req.get('analysis', None)
        size_by_group = req.get('size_by_group')
        markersize = req.get('marker_size', 3)  # 3 is default in lib/gear/plotting.py
        jitter = req.get('jitter', False)
        x_min = req.get('x_min')
        y_min = req.get('y_min')
        x_max = req.get('x_max')
        y_max = req.get('y_max')
        x_title = req.get('x_title')
        y_title = req.get('y_title')    # Will set later if not provided
        vlines = req.get('vlines', [])    # Array of vertical line dict properties
        filters = req.get('obs_filters', {})   # Dict of lists
        projection_id = req.get('projection_id', None)    # projection id of csv output
        colorblind_mode = req.get('colorblind_mode', False)
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        # Returning initial values in case plotting errors.
        # This is to help with unpredictable issues in the Vuex config where it expects to set values after the API call is done
        return_dict = {
            "success": None,
            "message": None,
            'gene_symbol': gene_symbol,
            "mapped_gene_symbol": None,
            'plot_json': None,
            "x_axis": x_axis,
            "y_axis": y_axis,
            "z_axis": z_axis,
            "x_min": x_min,
            "y_min": y_min,
            "x_max": x_max,
            "y_max": y_max,
            "x_title": x_title,
            "y_title": y_title,
            "vlines": vlines,
            "point_label": label,
            "size_by_group": size_by_group,
            "marker_size": markersize,
            "jitter": jitter,
            "hide_x_labels": hide_x_labels,
            "hide_y_labels": hide_y_labels,
            "hide_legend": hide_legend,
            "color_name": color_name,
            "facet_row": facet_row,
            "facet_col": facet_col,
            # only send back colormap for categorical color dimension
            "plot_colors": color_map if isinstance(color_map, dict) else None,
            "plot_palette": palette,
            "reverse_palette":reverse_palette,
            "obs_filters": filters,
            "plot_order": None,
        }

        if not gene_symbol or not dataset_id:
            return_dict["success"] = -1
            return_dict["message"] = "Request needs both dataset id and gene symbol."
            return return_dict

        try:
            ana = geardb.get_analysis(analysis, dataset_id, session_id)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return_dict["success"] = -1
            return_dict["message"] = "Could not retrieve analysis."
            return return_dict

        try:
            adata = ana.get_adata(backed=True)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return_dict["success"] = -1
            return_dict["message"] = "Could not retrieve AnnData."
            return return_dict

        # quick check to ensure x, y, color, facet columns are in the adata.obs
        if x_axis not in adata.obs.columns:
            if not x_axis == "raw_value":
                return_dict["success"] = -1
                return_dict["message"] = f"X-axis arg '{x_axis}' not found in observation metadata for dataset. Please update curation."
                return return_dict

        if y_axis not in adata.obs.columns:
            if not y_axis == "raw_value":
                return_dict["success"] = -1
                return_dict["message"] = f"Y-axis arg '{y_axis}' not found in observation metadata for dataset. Please update curation."
                return return_dict

        if color_name and color_name not in adata.obs.columns:
            if not color_name == "raw_value":
                return_dict["success"] = -1
                return_dict["message"] = f"Color arg '{color_name}' not found in observation metadata for dataset. Please update curation."
                return return_dict

        if facet_row and facet_row not in adata.obs.columns:
            return_dict["success"] = -1
            return_dict["message"] = f"Facet row arg '{facet_row}' not found in observation metadata for dataset. Please adupdatejust curation."
            return return_dict

        if facet_col and facet_col not in adata.obs.columns:
            return_dict["success"] = -1
            return_dict["message"] = f"Facet column arg '{facet_col}' not found in observation metadata for dataset. Please update curation."
            return return_dict

        if label and label not in adata.obs.columns:
            if not label == "raw_value":
                return_dict["success"] = -1
                return_dict["message"] = f"Label arg '{label}' not found in observation metadata for dataset. Please update curation."
                return return_dict


        if projection_id:
            try:
                adata = create_projection_adata(adata, dataset_id, projection_id)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

        adata.obs = order_by_time_point(adata.obs)

        # Reorder the categorical values in the observation dataframe
        try:
            if order:
                obs_keys = order.keys()
                for key in obs_keys:
                    if key not in adata.obs:
                        raise PlotError(f"Sort order series '{key}' not found in observation metadata for dataset ID {dataset_id}.")

                    col = adata.obs[key]
                    try:
                        # Some columns might be numeric, therefore
                        # we don't want to reorder these
                        reordered_col = col.cat.reorder_categories(
                            order[key], ordered=True)
                        adata.obs[key] = reordered_col
                    except:
                        pass

            # get a map of all levels for each column
            columns = adata.obs.select_dtypes(include="category").columns.tolist()

            if 'replicate' in columns:
                columns.remove('replicate')


            gene_symbols = (gene_symbol,)

            if 'gene_symbol' not in adata.var.columns:
                return_dict["success"] = -1
                return_dict["message"] = "The h5ad is missing the gene_symbol column."
                return return_dict

            # Filter genes and slice the adata to get a dataframe
            # with expression and its observation metadata
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            if not gene_filter.any():
                return_dict["success"] = -1
                return_dict["message"] = "The searched gene symbol could not be found in the dataset."
                return return_dict

            try:
                selected = adata[:, gene_filter].to_memory()
            except:
                # The "try" may fail for projections as it is already in memory
                selected = adata[:, gene_filter]

            # convert adata.X to a dense matrix if it is sparse
            # This prevents potential downstream issues
            try:
                selected.X = selected.X.todense()
            except:
                pass

            # Filter by obs filters
            if filters:
                for col, values in filters.items():
                    if col not in selected.obs:
                        raise PlotError(f"Filter series '{col}' not found in observation metadata for dataset ID {dataset_id}.")

                    # if there is an "NA" value in the filters but no "NA" in the dataframe
                    # check if it is a missing value, and if so, impute it
                    if "NA" in values and "NA" not in selected.obs[col].cat.categories:
                        values.remove("NA")
                        selected.obs[col].cat.add_categories("NA")
                        selected.obs[col].fillna("NA", inplace=True)

                    selected_filter = selected.obs[col].isin(values)
                    selected = selected[selected_filter, :]
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        order_res = dict()
        for col in columns:
            try:
                # Some columns might be numeric, therefore
                # we don't want to reorder these
                if col in [x_axis, color_name, facet_col, facet_row]:
                    order_res[col] = selected.obs[col].cat.categories.tolist()
            except:
                pass

        # Close adata so that we do not have a stale opened object
        if adata.isbacked:
            adata.file.close()

        success = 1
        message = ""

        # If there are multiple rows with the same gene symbol, we will only use the first one
        # But throw a warning message
        if len(selected.var) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)
            selected = selected[:, 0]

        # Rename the single selected.var index label to "raw_value"
        # This resolves https://github.com/IGS/gEAR/issues/878 where the gene_symbol index may be the same as a observation column (i.e. projections)
        selected.var.index = pd.Index(["raw_value"])

        df = selected.to_df()
        df = pd.concat([df,selected.obs], axis=1)

        # fill any missing adata.obs values with "NA"
        # The below line gives the error - TypeError: Cannot setitem on a Categorical with a new category (NA), set the categories first
        #df = df.fillna("NA")

        # Valid analysis column names from api/resources/h5ad.py
        analysis_tsne_columns = ['X_tsne_1', 'X_tsne_2']
        analysis_umap_columns = ['X_umap_1', 'X_umap_2']
        analysis_pca_columns = ['X_pca_1', 'X_pca_2']
        X, Y = (0, 1)

        # If analysis was performed on dataset, attempt to get coordinates from the obsm table
        if hasattr(selected, 'obsm'):
            if 'X_tsne' in selected.obsm:
                if x_axis in analysis_tsne_columns and y_axis in analysis_tsne_columns:
                    # A filtered AnnData object is an 'ArrayView' object and must
                    # be accessed as selected.obsm["X_tsne"] rather than selected.obsm.X_tsne
                    df[x_axis] = selected.obsm["X_tsne"].transpose()[X]
                    df[y_axis] = selected.obsm["X_tsne"].transpose()[Y]
            elif 'X_umap' in selected.obsm:
                if x_axis in analysis_umap_columns and y_axis in analysis_umap_columns:
                    df[x_axis] = selected.obsm["X_umap"].transpose()[X]
                    df[y_axis] = selected.obsm["X_umap"].transpose()[Y]
            elif 'X_pca' in selected.obsm:
                if x_axis in analysis_pca_columns and y_axis in analysis_pca_columns:
                    df[x_axis] = selected.obsm["X_pca"].transpose()[X]
                    df[y_axis] = selected.obsm["X_pca"].transpose()[Y]

        # Close adata so that we do not have a stale opened object
        if selected.isbacked:
            selected.file.close()

        if color_map and color_name:
            # Validate if all color map keys are in the dataframe columns
            # Ran into an issue where the color map keys were truncated compared to the dataframe column values
            col_values = set(df[color_name].unique())
            diff = col_values.difference(color_map.keys())
            if diff:
                message =  "WARNING: Color map has values not in the dataframe column '{}': {}\n".format(color_name, diff)
                message += "Will set color map key values to the unique values in the dataframe column."
                print(message, file=sys.stderr)
                # If any element in diff is nan and color_map contains a valid missing value key like "NA", change the value in the dataframe to match the color_map key
                for key in list(diff):
                    if pd.isna(key) and "NA" in color_map.keys():
                        df[color_name] = df[color_name].replace({key: "NA"})
                        col_values.remove(key)  # Remove the nan value from the set
                        col_values = col_values.union({"NA"})
                        break

                # Sort both the colormap and dataframe column alphabetically
                sorted_column_values = sorted(col_values)
                updated_color_map = {}
                # Replace all the colormap values with the dataframe column values
                # There is a good chance that the dataframe column values will be in the same order as the colormap values
                for idx, val in enumerate(sorted(color_map.keys())):
                    col_val = sorted_column_values[idx]
                    updated_color_map[col_val] = color_map[val]

                color_map = updated_color_map

        if color_name and not (color_map or palette):
            # For numerical color dimensions, we want to use
            # one of plotly's baked in scales.
            if isinstance(df[color_name].iloc[0], numbers.Number):
                purples = [
                    [0, 'rgb(218, 183, 193)'],
                    [0.35, 'rgb(194, 137, 166)'],
                    [0.5, 'rgb(169, 98, 151)'],
                    [0.6, 'rgb(145, 66, 143)'],
                    [0.7, 'rgb(105, 39, 122)'],
                    [1, 'rgb(63, 19, 98)']
                ]
                color_map = purples
            else:
                names = df[color_name].unique().tolist()
                color_map = plotly_color_map(names)

                # Check if color hexcodes exist and use if validated
                color_code = "{}_colors".format(color_name)
                if color_code in df.columns:
                    grouped = df.groupby([color_name, color_code])
                    # Ensure one-to-one mapping of color names to codes
                    if len(grouped) == len(names):
                        # Test if names are color hexcodes and use those if applicable
                        color_hex = df[color_code].unique().tolist()
                        if re.search(COLOR_HEX_PTRN, color_hex[0]):
                            color_map = {name[0]:name[1] for name, group in grouped}

        # Save original passed-in colormap or palette, so that it is not written by the colorblind version
        chosen_color_map = color_map if color_map else None
        chosen_palette = palette if palette else None

        # NOTE: If no color_name category, just leave color as "purple"
        # Using the reversed cividis scale so that higher expression values are darker
        if colorblind_mode:
            # Discrete scales = Viridis
            # Continuous scales = Reversed Cividis
            if palette:
                palette = "cividis_r"
            elif color_map:
                if isinstance(color_map, list):
                    color_map = pxc.get_colorscale("cividis_r")
                elif isinstance(color_map, dict):
                    num_entries = len(color_map)
                    viridis_colors =  pxc.get_colorscale("viridis")
                    sampled_colors = pxc.sample_colorscale(viridis_colors, num_entries)
                    color_map = {key: value for key, value in zip(color_map.keys(), sampled_colors)}

        if 'replicate' in df and plot_type == 'scatter':
            df = df.drop(['replicate'], axis=1)

        # kwargs will equal various options that should be passed to various plotly update functions
        # keys for kwargs: 'annotations', 'coloraxes', 'layout', 'traces', 'xaxes', 'yaxes'
        for k in ['annotations', 'coloraxes', 'layout', 'traces', 'xaxes', 'yaxes']:
            kwargs.setdefault(k, {})
        kwargs['traces']['marker'] = dict()

        # Assign custom marker size for scatter plot if passed
        if plot_type == "scatter":
            # Strip plots can only accept integer sizes
            if not jitter and size_by_group:
                # Dataframe series-based sizes
                kwargs['traces']['marker']['size'] = size_by_group
                kwargs['traces']['marker']['sizemin'] = markersize
            elif markersize:
                # Constant size
                kwargs['traces']['marker']['size'] = int(markersize)


        x_title = x_title if x_title else x_axis

        # Modification to y_axis title if title not provided
        if not y_title:
            y_title = y_axis
            if y_axis == "raw_value":
                y_title = "expression of {}".format(gene_symbol)
                if projection_id is not None:
                    y_title = "relative " + y_title

        if plot_type == "contour" and not z_axis:
            z_axis = "raw_value"
        elif not plot_type == "contour":
            z_axis = None   # Safeguard against unintended effects

        # Create plot
        try:
            fig = generate_plot(
                df,
                x=x_axis,
                y=y_axis,
                z=z_axis,
                color_name=color_name,
                facet_row=facet_row,
                facet_col=facet_col,
                # function has side effects and mutates colormap,
                # so we make a deepcopy so we can return the original
                # in response
                colormap=copy.deepcopy(color_map),
                palette=palette,    # NOTE: Maybe pass in colormap option and determine based on type?
                reverse_palette=True if reverse_palette else False,
                category_orders=order_res,
                plot_type=plot_type,
                text_name=label,
                jitter=jitter,
                hide_x_labels=hide_x_labels,
                hide_y_labels=hide_y_labels,
                hide_legend=hide_legend,
                x_range=[x_min, x_max],
                y_range=[y_min, y_max],
                x_title=x_title,
                y_title=y_title,
                vlines=vlines,
                is_projection=projection_id is not None,
                **kwargs
            )
        except PlotError as pe:
            return_dict["success"] = -1
            return_dict["message"] = str(pe)
            return return_dict
        except Exception as e:
            # print stacktrace to stderr
            import traceback
            traceback.print_exc()
            return_dict["success"] = -1
            return_dict["message"] = "Encountered error: {}".format(str(e))
            return return_dict

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # Modify y-title so that gene display results plot is not misleading
        # This only affects JSON return value, not the plot itself
        if "expression of {}".format(gene_symbol) in y_title:
            y_title = None

        return {
            "success": success,
            "message": message,
            'gene_symbol': gene_symbol,
            'plot_json': json.loads(plot_json),
            "x_axis": x_axis,
            "y_axis": y_axis,
            "z_axis": z_axis,
            "x_min": x_min,
            "y_min": y_min,
            "x_max": x_max,
            "y_max": y_max,
            "x_title": x_title,
            "y_title": y_title,
            "vlines": vlines,
            "point_label": label,
            "size_by_group": size_by_group,
            "marker_size": markersize,
            "jitter": jitter,
            "hide_x_labels": hide_x_labels,
            "hide_y_labels": hide_y_labels,
            "hide_legend": hide_legend,
            "color_name": color_name,
            "facet_row": facet_row,
            "facet_col": facet_col,
            # only send back colormap for categorical color dimension
            "plot_colors": chosen_color_map if isinstance(chosen_color_map, dict) else None,
            "plot_palette": chosen_palette,
            "reverse_palette":reverse_palette,
            "obs_filters": filters,
            "plot_order": order_res,
        }
