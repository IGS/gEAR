from flask import request
from flask_restful import Resource
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse
import json
import os
import sys
import copy
import geardb
import numbers
from gear.plotting import create_facet_grid, get_config, plotly_color_map
from plotly.utils import PlotlyJSONEncoder


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
        gene_symbol = req.get('gene_symbol')
        plot_type = req.get('plot_type')

        # tsne_dynamic is just a symlink to scatter
        if plot_type == 'tsne_dynamic':
            plot_type = 'scatter'

        x_axis = req.get('x_axis')
        y_axis = req.get('y_axis')
        label = req.get('point_label')
        hide_x_labels = req.get('hide_x_labels', False)
        hide_y_labels = req.get('hide_y_labels', False)
        color_name = req.get('color_name')
        color_map = req.get('colors')
        facet_row = req.get('facet_row')
        facet_col = req.get('facet_col')
        order = req.get('order', {})
        analysis = req.get('analysis', None)
        analysis_owner_id = req.get('analysis_owner_id', None)
        showlegend = req.get('showlegend', True)
        markersize = req.get('marker_size', 3)  # 3 is default in lib/gear/plotting.py
        jitter = req.get('jitter', 0)

        # If an analysis is posted we want to read from its h5ad
        if analysis:
            ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                  session_id=session_id, user_id=analysis_owner_id)

            if 'type' in analysis:
                ana.type = analysis['type']
            else:
                ana.discover_type(current_user_id=user.id)

            adata = sc.read_h5ad(ana.dataset_path(), backed='r')
        else:
            ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
            h5_path = ds.get_file_path()

            # Let's not fail if the file isn't there
            if not os.path.exists(h5_path):
                return {
                    'success': -1,
                    'message': "No h5 file found for this dataset"
                }

            adata = sc.read_h5ad(h5_path, backed='r')

        # check if time point order is intially provided in h5ad
        time_point_order = adata.obs.get('time_point_order')
        if (time_point_order is not None and 'time_point' in adata.obs.columns):
            sorted_df = adata.obs.drop_duplicates().sort_values(by='time_point_order')
            # Safety check. Make sure time point is categorical before
            # calling .cat
            adata.obs['time_point'] = pd.Categorical(adata.obs['time_point'])
            col = adata.obs['time_point'].cat
            adata.obs['time_point'] = col.reorder_categories(
                sorted_df.time_point.drop_duplicates(), ordered=True)
            adata.obs = adata.obs.drop(['time_point_order'], axis=1)

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
        columns = adata.obs.columns.tolist()
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

        gene_symbols = (gene_symbol,)

        if 'gene_symbol' in adata.var.columns:
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            if not gene_filter.any():
                return {
                    'success': -1,
                    'message': 'Gene not found',
                }

            if gene_filter.sum() > 1:
                return {
                    "success": -1,
                    "message": f"Multiple Ensembl IDs matched the gene: {gene_symbols[0]}"
                }
        else:
            return {
                'success': -1,
                'message': 'Missing gene_symbol in adata.var'
            }

        # Filter genes and slice the adata to get a dataframe
        # with expression and its observation metadata
        selected = adata[:, gene_filter]

        # which of these is faster?
        # option 1:
        df = selected.to_df()
        df2 = pd.concat([df,selected.obs], axis=1)
        df = df2.rename(columns={df.columns[0]: "raw_value"})

        # option 2:
        #df = pd.DataFrame({"raw_value": selected.X, **selected.obs})

        analysis_tsne_columns = ['X_tsne_1', 'X_tsne_2']
        if x_axis in analysis_tsne_columns and y_axis in analysis_tsne_columns:
            if hasattr(adata, "obsm") and hasattr(adata.obsm, 'X_tsne'):
                X, Y = (0, 1)
                # A filtered AnnData object is an 'ArrayView' object and must
                # be accessed as selected.obsm["X_tsne"] rather than selected.obsm.X_tsne
                df["X_tsne_1"] = selected.obsm["X_tsne"].transpose()[X]
                df["X_tsne_2"] = selected.obsm["X_tsne"].transpose()[Y]

        if not color_map and color_name:
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
                pass

        if 'replicate' in df and plot_type == 'scatter':
            df = df.drop(['replicate'], axis=1)

        kwargs = {'marker':dict()}

        # Assign custom marker size for scatter plot if passed
        if markersize and plot_type == "scatter":
            kwargs['marker']['size'] = int(markersize)

        # Create plot
        fig = create_facet_grid(
            df,
            x=x_axis,
            y=y_axis if y_axis else "raw_value",
            color_name=color_name,
            facet_row=facet_row,
            facet_col=facet_col,
            # function has side effects and mutates colormap,
            # so we make a deepcopy so we can return the original
            # in response
            colormap=copy.deepcopy(color_map),
            trace_type=plot_type,
            text_name=label,
            jitter=jitter,
            hide_x_labels=hide_x_labels,
            hide_y_labels=hide_y_labels,
            **kwargs
        )

        # If categorial, order the x-axis.  Ensure it is consistent in all faceted columns
        # TODO: A similar thing for y-axis
        if x_axis in order_res:
            for key in fig['layout']:
                if key.startswith('xaxis'):
                    # Some plot types should have empty columns, and others should not
                    if plot_type == "line":
                        fig['layout'][key]['categoryarray'] = order_res[x_axis]
                    else:
                        # Find the correct set of categories within the data traces
                        anchor = fig['layout'][key]['anchor']
                        axis_cats = []
                        for trace in fig['data']:
                            if trace['yaxis'] == anchor:
                                axis_cats.extend(trace['x'])    # Do not break loop.  Sometimes x may be an empty list and categories are in another trace
                        axis_cats_set = set(axis_cats)

                        # At this point, data may be unsorted.  Use "official" x-axis order as an index to sort this particular x-axis
                        # Source: https://yuji.wordpress.com/2014/05/06/how-to-sort-a-python-list-by-the-order-of-matching-values-in-another-list/
                        sorted_cats = sorted(list(axis_cats_set), key=lambda x: order_res[x_axis].index(x))
                        fig['layout'][key]['categoryarray'] = sorted_cats

        # Add titles to axes if expression data is the y-axis
        # If axis has no facets, then modify axis title,
        # else overwrite annot.text in plotting._axis_title_annotation
        if y_axis == "raw_value":
            if facet_col:
                for annot in fig['layout']['annotations']:
                    # The annotations for each faceted axis will have text.
                    if annot['name'] == 'x-title':
                        annot['text'] = x_axis
            else:
                fig['layout']['xaxis']['title']['text'] = x_axis

            y_text = "expression of {}".format(gene_symbol)
            if facet_row:
                for annot in fig['layout']['annotations']:
                    # The annotations for each faceted axis will have text.
                    if annot['name'] == 'y-title':
                        annot['text'] = y_text
            else:
                fig['layout']['yaxis']['title']['text'] = y_text

        fig.update_layout(
            margin=dict(
                b=100,
                r=50,
            ),
        )

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        return {
            "success": 1,
            'gene_symbol': gene_symbol,
            'plot_json': json.loads(plot_json),
            "x_axis": x_axis,
            "y_axis": y_axis,
            "point_label": label,
            "marker_size": markersize,
            "jitter": jitter,
            "hide_x_labels": hide_x_labels,
            "hide_y_labels": hide_y_labels,
            "color_name": color_name,
            "facet_row": facet_row,
            "facet_col": facet_col,
            "showlegend": showlegend,
            # only send back colormap for categorical color dimension
            "plot_colors": color_map if isinstance(color_map, dict) else None,
            "plot_config": get_config(),
            "plot_order": order_res
        }
