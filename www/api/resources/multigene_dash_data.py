import json

import gear.mg_plotting as mg
import geardb
import pandas as pd
from flask import request
from flask_restful import Resource
from plotly.utils import PlotlyJSONEncoder

from gear.mg_plotting import PlotError
from .common import create_projection_adata, order_by_time_point


# SAdkins - 2/15/21 - This is a list of datasets already log10-transformed where if selected will use log10 as the default dropdown option
# This is meant to be a short-term solution until more people specify their data is transformed via the metadata
LOG10_TRANSFORMED_DATASETS = [
"320ca057-0119-4f32-8397-7761ea084ed1"
, "df726e89-b7ac-d798-83bf-2bd69d7f3b52"
, "bad48d04-db27-26bc-2324-e88506f751fd"
, "dbd715bf-778a-4923-6fe7-c587987cdb00"
, "c8d99d13-394f-a87f-5d3a-395968fdb619"
, "bee735e5-d180-332c-7892-dd751dd76bb8"
, "c4f16a12-9e98-47be-4335-b8321282919e"
, "17a07bf4-b41a-d9c3-9aa7-b4729390f57a"
, "6a0a2bca-0f86-59d0-4e3d-4457be3a71ff"
, "39e01b71-415f-afa7-0c64-f0e996be0fb7"
, "6482c608-a6bd-d8b1-6bc1-5b53c34ed61c"
, "0c5a4c18-c2a9-930c-6e52-ef411f54eb67"
, "3c02d449-61ab-4bcd-f100-5f5937b1794e"
, "23e3797f-3016-8142-cbe8-69b03131ad95"
, "b16eeb8d-d68e-c7c9-9dc9-a3f4821e9192"
, "b96f448a-315d-549d-6e8a-83cdf1ce1b5c"
, "b0420910-a0fa-e920-152d-420b6275d3af"
, "f1ce4e63-3577-8020-8307-e88f1fb98953"
, "2f79f784-f7f7-7dc3-9b3e-4c87a4346d91"
, "c32835d3-cac4-bb0e-a90a-0b41dec6617a"
, "48bab518-439e-4a17-b868-6b225abf2c73"    # Carlo dataset
, "1b12dde9-1762-7564-8fbd-1b07b750505f"
, "a2dd9f06-5223-0779-8dfc-8dce7a3897e1"
, "f7de7db2-b4cb-ebe3-7f1f-b278f46f1a7f"
, "e34fa5c6-1083-cacb-eedf-23f59f2e005f"
, "0c5fb6b0-31ab-6bfc-075d-76756ccd56b4"
, "a183b2e6-ab38-458a-52a6-5eb014d073da"
, "c4f16a12-9e98-47be-4335-b8321282919e"
, "2a25e445-2776-8913-076f-9a147a43e8b4"
, "2786d849-f11c-2de6-b22e-12c940aafe07"
, "2e3423b3-74db-d436-8357-abb3031d47e9"
, "4cb2ac62-c283-86a9-83cb-2c1b381948f2"
, "d0659d69-1a33-8b84-252c-f7ded46aa3d6"
, "cee5325d-434f-fefe-d2e6-e0be39421951"    # Unsure if log2 or log10 but necessary for correct output
, "34f8f131-8158-db83-7df9-db9003797dff"    # same as above
, "7ddb4965-e710-faf7-ee26-4ce95d7602a8"    # same as above
, "f122cac5-c79f-8ea2-166e-42415916db11"    # ""
, "173ab634-a2b1-87bc-f1ef-d288de0bcd1a"    # ""
, "80eadbe6-49ac-8eaf-f2fb-e07706cf117b"    # HRP dataset
]

CLUSTER_LIMIT = 5000

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
    return df.obs[columns].apply(lambda x: ';'.join(map(str, x)), axis=1)

class MultigeneDashData(Resource):
    """Resource for retrieving data from h5ad to be used to draw charts on UI.
    Parameters
    ----------
    gene_symbols: str
        Genes to search in adata.
    plot_type: str
        Plot type (heatmap, mg_violin, volcano)
    Returns
    -------
    dict
        Plot data
    """

    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        analysis = req.get('analysis', None)
        plot_type = req.get('plot_type')
        gene_symbols = req.get('gene_symbols', [])
        filters = req.get('obs_filters', {})    # Dict of lists
        primary_col = req.get('primary_col', None)
        secondary_col = req.get('secondary_col', None)
        sort_order = req.get('sort_order', {})
        colorscale = req.get('colorscale', None)    # If None, this can override the "reverse_colorscale" param by using defaults.
        reverse_colorscale = req.get('reverse_colorscale', False)
        # Heatmap opts
        clusterbar_fields = req.get('clusterbar_fields', [])
        subsample_limit = req.get('subsample_limit', 0)
        matrixplot = req.get('matrixplot', False)
        center_around_zero = req.get('center_around_zero', False)
        cluster_obs = req.get('cluster_obs', False)
        cluster_genes = req.get('cluster_genes', False)
        flip_axes = req.get('flip_axes', False)
        distance_metric = req.get('distance_metric', "euclidean")
        hide_obs_labels = req.get('hide_obs_labels', False)
        hide_gene_labels = req.get('hide_gene_labels', False)
        # Quadrant plot options
        compare_group1 = req.get("compare1_condition", None)
        compare_group2 = req.get("compare2_condition", None)
        fc_cutoff = float(req.get("fold_change_cutoff", 2))
        fdr_cutoff = float(req.get("fdr_cutoff", 0.05))
        include_zero_fc = req.get("include_zero_fc", True)
        # Volcano plot options
        pval_threshold = float(req.get("pvalue_threshold", 0.05))
        lower_logfc_threshold = float(req.get("lower_logfc_threshold", -2))
        upper_logfc_threshold = float(req.get("upper_logfc_threshold", 2))
        query_condition = req.get('query_condition', None)
        ref_condition = req.get('ref_condition', None)
        de_test_algo = req.get("de_test_algo", "t-test")
        use_adj_pvals = req.get('adj_pvals', True)
        annotate_nonsignificant = req.get('annotate_nonsignificant', True)
        # Violin plot options
        stacked_violin = req.get('stacked_violin', False)
        violin_add_points = req.get('violin_add_points', False)
        # Misc options
        title = req.get('plot_title', None)
        legend_title = req.get('legend_title', None)
        projection_id = req.get('projection_id', None)    # projection id of csv output
        colorblind_mode = req.get('colorblind_mode', False)
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        try:
            ana = geardb.get_analysis(analysis, dataset_id, session_id)
            # Using adata with "backed" mode does not work with volcano plot
            adata = ana.get_adata(backed=False)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return {
                'success': -1,
                'message': str(e),
            }

        adata.obs = order_by_time_point(adata.obs)

        # get a map of all levels for each column
        columns = adata.obs.select_dtypes(include="category").columns.tolist()

        if 'replicate' in columns:
            columns.remove('replicate')

        # remove _colors columns
        columns = [col for col in columns if not col.endswith('_colors')]

        # Remove any columns that are have more than 50 unique values
        # These are likely not categorical columns (i.e. barcodes) and can cause issues with the composite index
        for col in columns:
            if len(adata.obs[col].unique()) > 50:
                columns.remove(col)

        if not columns:
            return {
                "success": -1,
                "message": "There are no categorical-datatype conditions found in this dataset."
            }

        # Ensure datasets are not doubly log-transformed
        # In the case of projection inputs, we don't want to log-transform either
        is_log10 = False
        if dataset_id in LOG10_TRANSFORMED_DATASETS or projection_id:
            is_log10 = True

        success = 1
        message = ""

        if projection_id:
            try:
                adata = create_projection_adata(adata, dataset_id, projection_id)
                # For plots where var.index is used, we need to assign that a name (currently empty)
                adata.var.index = adata.var.index.rename("index")
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

        # Success levels
        # -1 Failure
        # 1 Success
        # 2 Warning - duplicate genes found
        # 3 Warning - One or more genes could not be processed
        # NOTE: The success level in a warning can be overridden by another warning or error

        # delete all "None" values from the gene_symbols list
        # These are genes that did not map in the orthology mapping
        gene_symbols = [gene for gene in gene_symbols if gene]

        # TODO: How to deal with a gene mapping to multiple Ensemble IDs
        try:
            if not gene_symbols and plot_type in ["dotplot", "heatmap", "mg_violin"]:
                raise PlotError('Must pass in some genes before creating a plot of type {}'.format(plot_type))

            if len(gene_symbols) == 1 and plot_type == "heatmap":
                raise PlotError('Heatmaps require 2 or more genes as input')

            # Some datasets have multiple ensemble IDs mapped to the same gene.
            # Drop dups to prevent out-of-bounds index errors downstream
            gene_filter, success, message = mg.create_dataframe_gene_mask(adata.var, gene_symbols)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # ADATA - Observations are rows, genes are columns
        try:
            selected = adata.to_memory()
        except:
            # The "try" may fail for projections as it is already in memory
            selected = adata

        # These plot types filter to only the specific genes.
        # The other plot types use all genes and rather annotate the specific ones.
        if plot_type in ['dotplot', 'heatmap', 'mg_violin'] and gene_filter is not None:
            selected = selected[:, gene_filter]

            try:
                if plot_type == "heatmap" and len(selected.var) == 1:
                    raise PlotError("Only one gene from the searched gene symbols was found in dataset.  The heatmap option require at least 2 genes to plot.")
            except PlotError as pe:
                return {
                    "success": -1,
                    "message": str(pe)
                }
            # Get a list of sorted ensembld IDs based on the specified gene symbol order
            ensm_to_gene = selected.var.to_dict()["gene_symbol"]
            # Since genes were restricted to one Ensembl ID we can invert keys and vals
            gene_to_ensm = {y:x for x,y in ensm_to_gene.items()}

            # Collect all genes from the unfiltered dataset
            dataset_genes = adata.var['gene_symbol'].unique().tolist()
            # Gene symbols list may have genes not in the dataset.
            normalized_genes_list, _found_genes = mg.normalize_searched_genes(dataset_genes, gene_symbols)

            # deduplicate normalized_genes_list
            normalized_genes_list = list(dict.fromkeys(normalized_genes_list))

            # Sort ensembl IDs based on the gene symbol order
            sorted_ensm = map(lambda x: gene_to_ensm[x], normalized_genes_list)

        # Reorder the categorical values in the observation dataframe
        sort_fields = []
        if sort_order:
            # Ensure selected primary and secondary columns are in the correct order
            if primary_col:
                sort_fields.append(primary_col)
            if secondary_col and secondary_col != primary_col:
                sort_fields.append(secondary_col)

            # Now reorder the dataframe
            for key in sort_fields:
                col = selected.obs[key]
                try:
                    # Some columns might be numeric, therefore
                    # we don't want to reorder these
                    reordered_col = col.cat.reorder_categories(
                        sort_order[key], ordered=True)
                    selected.obs[key] = reordered_col
                except:
                    pass

        # Make a composite index of all categorical types
        selected.obs['composite_index'] = create_composite_index_column(selected, columns)
        selected.obs['composite_index'] = selected.obs['composite_index'].astype('category')
        columns.append("composite_index")

        import sys

        # Filter dataframe on the chosen observation filters
        if filters:
            # reorder filter key the same order as sort_order if key exists
            # Mostly for fixing the order of the heatmap clusterbars
            for field in filters.keys():
                values = filters[field]
                # if there is an "NA" value in the filters but no "NA" in the dataframe
                # check if it is a missing value, and if so, impute it
                if "NA" in values and "NA" not in selected.obs[col].cat.categories:
                    values.remove("NA")
                    selected.obs[col].cat.add_categories("NA")
                    selected.obs[col].fillna("NA", inplace=True)

                mask = selected.obs[field].isin(values)
                selected = selected[mask, :]

                if sort_order and field in sort_order:
                    filters[field] = sort_order[field]

            # if the filters are empty, return an empty dataframe
            if selected.shape[0] == 0:
                return {
                    "success": -1,
                    "message": "No data found for the selected filters."
                }

            # For the remaining data, create a special composite index for the specified filters
            selected.obs['filters_composite'] = create_composite_index_column(selected, filters.keys())
            selected.obs['filters_composite'] = selected.obs['filters_composite'].astype('category')
            columns.append("filters_composite")

            # sort by the filters
            for key in sort_fields:
                col = selected.obs[key]
                reordered_col = col.cat.reorder_categories(
                    filters[key], ordered=True)
                selected.obs[key] = reordered_col

        var_index = selected.var.index.name

        if plot_type == "volcano":
            try:
                key, query_val, ref_val = mg.validate_volcano_conditions(query_condition, ref_condition)
                df = mg.prep_volcano_dataframe(selected
                    , key
                    , query_val
                    , ref_val
                    , de_test_algo
                    , is_log10
                    )
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

            # Build a dictionary to easily move gene_syms to "text" property and ensembl ids to "customdata" property
            ensm2genesymbol = pd.Series(df["gene_symbol"].values, index=df["ensm_id"]).to_dict()

            # Volcano plot expects specific parameter names (unless we wish to change the options)
            fig = mg.create_volcano_plot(df
                , query_val
                , ref_val
                , pval_threshold
                , [lower_logfc_threshold, upper_logfc_threshold]
                , use_adj_pvals
                )

            downcolor = None
            upcolor = None
            # Use reversed "cividis" colorscheme
            if colorblind_mode:
                downcolor = "rgb(254, 232, 56)"
                upcolor = "rgb(0, 34, 78)"

            mg.modify_volcano_plot(fig, query_val, ref_val, ensm2genesymbol, downcolor, upcolor)

            if gene_symbols:
                dataset_genes = df['gene_symbol'].unique().tolist()
                normalized_genes_list, _found_genes = mg.normalize_searched_genes(dataset_genes, gene_symbols)
                mg.add_gene_annotations_to_volcano_plot(fig, normalized_genes_list, annotate_nonsignificant)

        elif plot_type == "quadrant":
            # Get list of normalized genes before dataframe filtering takes place
            if gene_symbols:
                dataset_genes = adata.var['gene_symbol'].unique().tolist()
                normalized_genes_list, _found_genes = mg.normalize_searched_genes(dataset_genes, gene_symbols)
            try:
                key, control_val, compare1_val, compare2_val = mg.validate_quadrant_conditions(ref_condition, compare_group1, compare_group2)
                df = mg.prep_quadrant_dataframe(selected
                        , key
                        , control_val
                        , compare1_val
                        , compare2_val
                        , de_test_algo
                        , fc_cutoff
                        , fdr_cutoff
                        , include_zero_fc
                        , is_log10
                        )
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

            colorscale = "viridis" if colorblind_mode else None

            fig = mg.create_quadrant_plot(df, control_val, compare1_val, compare2_val, colorscale)
            # Annotate selected genes
            if gene_symbols:
                genes_not_found, genes_none_none = mg.add_gene_annotations_to_quadrant_plot(fig, normalized_genes_list)
                if genes_not_found:
                    success = 2
                    message += "<li>One or more genes did not pass cutoff filters to be in the plot: {}</li>".format(', '.join(genes_not_found))
                if genes_none_none:
                    success = 2
                    message += "<li>One or more genes had no fold change in both comparisons and will not be annotated: {}</li>".format(', '.join(genes_none_none))

        elif plot_type == "dotplot":
            df = selected.to_df()

            # Sort Ensembl ID columns by the gene symbol order
            df = df[sorted_ensm]

            if not primary_col:
                return {
                    'success': -1,
                    'message': "The 'primary_col' option required for dot plots."
                }

            groupby_filters = [primary_col]
            if secondary_col and not primary_col == secondary_col:
                groupby_filters.append(secondary_col)

            for gb in groupby_filters:
                df[gb] = selected.obs[gb]

            # 1) Flatten to long-form
            # 2) Create a gene symbol column by mapping to the Ensembl IDs
            df = df.melt(id_vars=groupby_filters)

            # Add "gene_symbol" as a column, make it categorical to ensure the sort order is preserved when melted
            df["gene_symbol"] = df[var_index].map(ensm_to_gene).astype('category')
            df["gene_symbol"] = df["gene_symbol"].cat.reorder_categories(
                        normalized_genes_list, ordered=True)
            df = df.sort_values(by=["gene_symbol"])

            # Percent of all cells in this group where the gene has expression
            percent = lambda row: round(len([num for num in row if num > 0]) / len(row) * 100, 2)
            groupby = ["gene_symbol"]
            groupby.extend(groupby_filters)

            # drop Ensembl ID index since it may not aggregate and throw warnings
            df.drop(columns=[var_index], inplace=True)

            grouped = df.groupby(groupby, observed=True)
            df = grouped.agg(['mean', 'count', ('percent', percent)]) \
                .fillna(0) \
                .reset_index()

            # Reverse Cividis so that dark is higher expression
            if colorblind_mode:
                colorscale = "cividis_r"

            fig = mg.create_dot_plot(df, groupby_filters, is_log10, title, colorscale, reverse_colorscale)

        elif plot_type == "heatmap":
            # Filter genes and slice the adata to get a dataframe
            # with expression and its observation metadata
            df = selected.to_df()

            # Sort Ensembl ID columns by the gene symbol order
            df = df[sorted_ensm]

            # Enabling subsampling to deal with potential memory issues for clustering.
            # If clustering on observations, limit samples to 10,000 or fewer
            # If a subsampling limit was set, sample based on the min of these two values
            if subsample_limit > len(df) or subsample_limit == 0:
                subsample_limit = len(df)
            if cluster_obs and len(df) > CLUSTER_LIMIT:
                subsample_limit = min(subsample_limit, CLUSTER_LIMIT)
            df = df.sample(subsample_limit, random_state=1)

            groupby_index = "composite_index"
            groupby_fields = columns
            # Create a composite to groupby
            if filters and matrixplot:
                union_fields = mg.union(list(filters.keys()), clusterbar_fields)
                selected.obs['groupby_composite'] = create_composite_index_column(selected, union_fields)
                selected.obs['groupby_composite'] = selected.obs['groupby_composite'].astype('category')
                columns.append("groupby_composite")
                groupby_index = "groupby_composite"
                union_fields.extend([groupby_index, "filters_composite"])   # Preserve filters index for downstream labeling
                groupby_fields = union_fields

            # Only add the fields that will be used downstream
            for cat in groupby_fields:
                df[cat] = selected.obs[cat]

            # Groupby to remove the replicates
            # Ensure the composite index is used as the index for plot labeling
            if matrixplot:
                grouped = df.groupby(groupby_fields, observed=False)
                df = grouped.mean() \
                    .dropna() \
                    .reset_index() \
                    .set_index(groupby_index)

            # Since this is the new index in the matrixplot, it does not exist as a droppable series
            # These two statements ensure that the current columns are the same if that option is set or not
            groupby_fields.remove(groupby_index)
            if not matrixplot:
                df = df.drop(columns=groupby_index)

            # Sort based on the specified sort order
            if primary_col:
                sortby_fields = [primary_col]
                if secondary_col and not primary_col == secondary_col:
                    sortby_fields.append(secondary_col)
                df = df.sort_values(by=sortby_fields)

            # Drop the obs metadata now that the dataframe is sorted
            # They cannot be in there when the clustergram is made
            # But save it to add back in later
            df_cols = pd.concat([df.pop(cat) for cat in groupby_fields], axis=1)

            # Reverse Cividis so that dark is higher expression
            if colorblind_mode:
                colorscale = "cividis_r"

            # "df" must be obs label for rows and genes for cols only
            fig = mg.create_clustergram(df
                , normalized_genes_list
                , is_log10
                , cluster_obs
                , cluster_genes
                , flip_axes
                , center_around_zero
                , distance_metric
                , colorscale
                , reverse_colorscale
                , hide_obs_labels
                , hide_gene_labels
                )

            # Need the obs metadata again for mapping clusterbars to indexes
            df = pd.concat([df, df_cols], axis=1)

            clusterbar_indexes = mg.build_obs_group_indexes(df, filters, clusterbar_fields)

            # Create labels based only on the included filters
            obs_labels = None
            if "filters_composite" in df and matrixplot:
                obs_labels = mg.create_clustergram_observation_labels(df, fig, "filters_composite", flip_axes)

            mg.add_clustergram_cluster_bars(fig, clusterbar_indexes, obs_labels, is_log10, flip_axes)

        elif plot_type == "mg_violin":
            df = selected.to_df()

            # Sort Ensembl ID columns by the gene symbol order
            df = df[sorted_ensm]

            if not primary_col:
                return {
                    'success': -1,
                    'message': "The 'primary_col' option required for violin plots."
                }

            groupby_filters = [primary_col]
            if secondary_col and not primary_col == secondary_col:
                groupby_filters.append(secondary_col)

            for gb in groupby_filters:
                df[gb] = selected.obs[gb]

            # 1) Flatten to long-form
            # 2) Create a gene symbol column by mapping to the Ensembl IDs
            df = df.melt(id_vars=groupby_filters)

            # Add "gene_symbol" as a column, make it categorical to ensure the sort order is preserved when melted
            df["gene_symbol"] = df[var_index].map(ensm_to_gene).astype('category')
            df["gene_symbol"] = df["gene_symbol"].cat.reorder_categories(
                        normalized_genes_list, ordered=True)

            # This is a bit redundant for regular violin plots, which sort correctly without sorting by the "groupby_fields"
            # But the stacked violin plot does not respect the sorted order of the categories after melting, so we need to sort again
            sortby_fields = ["gene_symbol"]
            sortby_fields.extend(groupby_filters)
            df = df.sort_values(by=sortby_fields)

            violin_func = mg.create_stacked_violin_plot if stacked_violin else mg.create_violin_plot

            # I think Viridis lends itself to quantitative plots than Cividis.
            if colorblind_mode:
                colorscale = "viridis"

            fig = violin_func(df
                , groupby_filters
                , is_log10
                , colorscale
                , reverse_colorscale
                )

            # Add jitter-based args (to make beeswarm plot)
            if violin_add_points:
                fig.update_traces(
                    jitter=0.25
                    , points="all"
                    , pointpos=0
                    , marker=dict(
                        color="#000000"
                        , size = 6 if len(groupby_filters) == 1 else 3
                    )
                )
        else:
            return {
                'success': -1,
                'message': "Plot type {} is not a valid multi-gene plot option".format(plot_type)
            }

        # Close adata so that we do not have a stale opened object
        if adata.isbacked:
            adata.file.close()

        # If figure is actualy a JSON error message, send that instead
        if "success" in fig and fig["success"] == -1:
            return fig

        fig.update_layout(autosize=True)

        # change background to pure white
        # Heatmap/Clustergram already does this, but this option adds extra ticks
        if not plot_type == "heatmap":
            fig.update_layout(
                template="simple_white"
            )

        # Title is addressed in the creation of dotplot subplots
        # But we can add it here for other plots
        if title and not plot_type == "dotplot":
            fig.update_layout(
                title={
                    "text":title
                    ,"x":0.5
                    ,"xref":"paper"
                    ,"y":0.9
                }
            )

        if legend_title:
            fig.update_layout(
                legend={
                    "title":{
                        "text":legend_title
                    }
                }
            )

        # Change plot elements to indicuate projections instead of genes
        if projection_id:
            fig.update_layout(
                legend_title_text=fig.layout.legend.title.text.replace("gene", "projection").replace("Gene", "Projection") if fig.layout.legend.title.text else None
                , title_text=fig.layout.title.text.replace("gene", "projection").replace("Gene", "Projection") if fig.layout.title.text else None
            )
            fig.for_each_xaxis(
                lambda a: a.update(
                    title={
                        "text": a.title.text.replace("gene", "projection").replace("Gene", "Projection") if a.title.text else None

                    }
                )
            )
            fig.for_each_yaxis(
                lambda a: a.update(
                    title={
                        "text": a.title.text.replace("gene", "projection").replace("Gene", "Projection") if a.title.text else None

                    }
                )
            )

        # Pop any default height and widths being added
        fig["layout"].pop("height", None)
        fig["layout"].pop("width", None)

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # NOTE: With volcano plots, the Chrome "devtools" cannot load the JSON response occasionally
        return {
            "success": success
            , "message": message
            , 'plot_json': json.loads(plot_json)
        }
