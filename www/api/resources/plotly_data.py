from flask import request
from flask_restful import Resource
import pandas as pd
import copy, json, os, re
import geardb
import numbers
from gear.plotting import generate_plot, plotly_color_map, PlotError
from plotly.utils import PlotlyJSONEncoder

COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"

def create_projection_adata(dataset_adata, projection_csv):
    # Create AnnData object out of readable CSV file
    # ? Does it make sense to put this in the geardb/Analysis class?
    try:
        import scanpy as sc
        projection_adata = sc.read_csv("/tmp/{}".format(projection_csv))
    except:
        raise PlotError("Could not create projection AnnData object from CSV.")

    projection_adata.obs = dataset_adata.obs
    projection_adata.var["gene_symbol"] = projection_adata.var_names
    return projection_adata

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
        user = geardb.get_user_from_session_id(session_id)
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
        analysis_owner_id = req.get('analysis_owner_id', None)
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
        projection_id = req.get('projection_id', None)    # projection id of csv output
        colorblind_mode = req.get('colorblind_mode', False)
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        # Returning initial values in case plotting errors.
        # This is to help with unpredictable issues in the Vuex config where it expects to set values after the API call is done
        return_dict = {
            "success": None,
            "message": None,
            'gene_symbol': gene_symbol,
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
            "plot_order": None,
        }

        if not gene_symbol or not dataset_id:
            return_dict["success"] = -1
            return_dict["message"] = "Request needs both dataset id and gene symbol."
            return return_dict

        try:
            ana = get_analysis(analysis, dataset_id, session_id, analysis_owner_id)
        except PlotError as pe:
            return_dict["success"] = -1
            return_dict["message"] = str(pe)
            return return_dict

        adata = ana.get_adata(backed=True)
        adata.obs = order_by_time_point(adata.obs)

        # Reorder the categorical values in the observation dataframe
        if order:
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

        # get a map of all levels for each column
        columns = adata.obs.select_dtypes(include="category").columns.tolist()

        if 'replicate' in columns:
            columns.remove('replicate')

        order_res = dict()
        for col in columns:
            try:
                # Some columns might be numeric, therefore
                # we don't want to reorder these
                if col in [x_axis, color_name, facet_col, facet_row]:
                    order_res[col] = adata.obs[col].cat.categories.tolist()
            except:
                pass

        if projection_id:
            projection_csv = "{}.csv".format(projection_id)
            try:
                adata = create_projection_adata(adata, projection_csv)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

        gene_symbols = (gene_symbol,)

        if 'gene_symbol' in adata.var.columns:
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            if not gene_filter.any():
                return_dict["success"] = -1
                return_dict["message"] = 'Gene not found in dataset'
                return return_dict
        else:
            return_dict["success"] = -1
            return_dict["message"] = 'Missing gene_symbol in adata.var'
            return return_dict

        # Filter genes and slice the adata to get a dataframe
        # with expression and its observation metadata
        selected = adata[:, gene_filter]

        df = selected.to_df()

        success = 1
        message = ""
        if len(df.columns) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)
            df = df.iloc[:,[0]] # Note, put the '0' in a list to return a DataFrame.  Not having in list returns DataSeries instead

        df2 = pd.concat([df,selected.obs], axis=1)
        df = df2.rename(columns={df.columns[0]: "raw_value"})

        # Valid analysis column names from api/resources/h5ad.py
        analysis_tsne_columns = ['X_tsne_1', 'X_tsne_2']
        analysis_umap_columns = ['X_umap_1', 'X_umap_2']
        analysis_pca_columns = ['X_pca_1', 'X_pca_2']
        X, Y = (0, 1)

        # If analysis was performed on dataset, attempt to get coordinates from the obsm table
        if hasattr(adata, 'obsm'):
            if 'X_tsne' in adata.obsm:
                if x_axis in analysis_tsne_columns and y_axis in analysis_tsne_columns:
                    # A filtered AnnData object is an 'ArrayView' object and must
                    # be accessed as selected.obsm["X_tsne"] rather than selected.obsm.X_tsne
                    df[x_axis] = selected.obsm["X_tsne"].transpose()[X]
                    df[y_axis] = selected.obsm["X_tsne"].transpose()[Y]
            elif 'X_umap' in adata.obsm:
                if x_axis in analysis_umap_columns and y_axis in analysis_umap_columns:
                    df[x_axis] = selected.obsm["X_umap"].transpose()[X]
                    df[y_axis] = selected.obsm["X_umap"].transpose()[Y]
            elif 'X_pca' in adata.obsm:
                if x_axis in analysis_pca_columns and y_axis in analysis_pca_columns:
                    df[x_axis] = selected.obsm["X_pca"].transpose()[X]
                    df[y_axis] = selected.obsm["X_pca"].transpose()[Y]

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
                **kwargs
            )
        except PlotError as pe:
            return_dict["success"] = -1
            return_dict["message"] = str(pe)
            return return_dict
        except Exception as e:
            return_dict["success"] = -1
            return_dict["message"] = "ERROR: {}. Please contact the gEAR team if you need help resolving this".format(str(e))
            return return_dict

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # Modify y-title so that gene display results plot is not misleading
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
            "plot_colors": color_map if isinstance(color_map, dict) else None,
            "plot_palette": palette,
            "reverse_palette":reverse_palette,
            "plot_order": order_res,
        }
