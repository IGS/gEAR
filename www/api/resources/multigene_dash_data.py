from flask import request
from flask_restful import Resource
import pandas as pd
import geardb

import json, os
from gear.mg_plotting import PlotError
import gear.mg_plotting as mg

from plotly.utils import PlotlyJSONEncoder


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
#, "fbe1296e-572c-d388-b9d1-6e2a6bf10b0a"
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

def get_analysis(analysis, dataset_id, session_id, analysis_owner_id):
    """Return analysis object based on various factors."""
    # If an analysis is posted we want to read from its h5ad
    if analysis:
        ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=analysis_owner_id)

        try:
            ana.type = analysis['type']
        except:
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
        analysis_owner_id = req.get('analysis_owner_id', None)
        plot_type = req.get('plot_type')
        gene_symbols = req.get('gene_symbols', [])
        filters = req.get('obs_filters', {})    # Dict of lists
        primary_col = req.get('primary_col', None)
        secondary_col = req.get('secondary_col', None)
        sort_order = req.get('sort_order', {})
        # Heatmap opts
        clusterbar_fields = req.get('clusterbar_fields', [])
        matrixplot = req.get('matrixplot', False)
        center_around_zero = req.get('center_around_zero', False)
        cluster_obs = req.get('cluster_obs', False)
        cluster_genes = req.get('cluster_genes', False)
        flip_axes = req.get('flip_axes', False)
        distance_metric = req.get('distance_metric', "euclidean")
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
        projection_csv = req.get('projection_csv', None)    # as CSV file
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        try:
            ana = get_analysis(analysis, dataset_id, session_id, analysis_owner_id)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # Using adata with "backed" mode does not work with volcano plot
        adata = ana.get_adata(backed=False)

        adata.obs = order_by_time_point(adata.obs)

        # Reorder the categorical values in the observation dataframe
        if sort_order:
            obs_keys = sort_order.keys()
            for key in obs_keys:
                col = adata.obs[key]
                try:
                    # Some columns might be numeric, therefore
                    # we don't want to reorder these
                    reordered_col = col.cat.reorder_categories(
                        sort_order[key], ordered=True)
                    adata.obs[key] = reordered_col

                    # Ensure filter order aligns with sort order
                    # NOTE: Off-chance sort order may have more elements than filter key
                    # if filters are changed after sort order list is generated.
                    filters[key] = sort_order[key]
                except:
                    pass

        # get a map of all levels for each column
        columns = adata.obs.select_dtypes(include="category").columns.tolist()

        if 'replicate' in columns:
            columns.remove('replicate')

        if not columns:
            return {
                "success": -1,
                "message": "There are no categorical-datatype conditions found in this dataset."
            }

        # Ensure datasets are not doubly log-transformed
        # In the case of projection inputs, we don't want to log-transform either
        is_log10 = False
        if dataset_id in LOG10_TRANSFORMED_DATASETS or projection_csv:
            is_log10 = True

        success = 1
        message = ""

        if projection_csv:
            try:
                adata = create_projection_adata(adata, projection_csv)
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

        # TODO: How to deal with a gene mapping to multiple Ensemble IDs
        try:

            if not gene_symbols and plot_type in ["dotplot", "heatmap", "mg_violin"]:
                raise PlotError('Must pass in some genes before creating a plot of type {}'.format(plot_type))

            if len(gene_symbols) == 1 and plot_type == "heatmap":
                raise PlotError('Heatmaps require 2 or more genes as input')

            # Some datasets have multiple ensemble IDs mapped to the same gene.
            # Drop dups to prevent out-of-bounds index errors downstream
            #var = adata.var.drop_duplicates(subset=['gene_symbol'])
            gene_filter, success, message = mg.create_dataframe_gene_mask(adata.var, gene_symbols)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # ADATA - Observations are rows, genes are columns
        selected = adata

        # These plot types filter to only the specific genes.
        # The other plot types use all genes and rather annotate the specific ones.
        if plot_type in ['dotplot', 'heatmap', 'mg_violin'] and gene_filter is not None:
            selected = selected[:, gene_filter]

            if plot_type == "heatmap" and len(selected.var) == 1:
                raise PlotError("Only one gene from the searched gene symbols was found in dataset.  The heatmap option require at least 2 genes to plot.")

        # Make a composite index of all categorical types
        selected.obs['composite_index'] = selected.obs[columns].apply(lambda x: ';'.join(map(str,x)), axis=1)
        selected.obs['composite_index'] = selected.obs['composite_index'].astype('category')
        columns.append("composite_index")

        # Filter dataframe on the chosen observation filters
        if filters:
            # Create a special composite index for the specified filters
            selected.obs['filters_composite'] = selected.obs[filters.keys()].apply(lambda x: ';'.join(map(str,x)), axis=1)
            selected.obs['filters_composite'] = selected.obs['filters_composite'].astype('category')
            columns.append("filters_composite")
            unique_composite_indexes = selected.obs["filters_composite"].unique()

            # Only want to keep indexes that match chosen filters
            # However if no filters were chosen, just use everything
            filtered_composite_indexes = mg.create_filtered_composite_indexes(filters, unique_composite_indexes.tolist())
            if filtered_composite_indexes:
                condition_filter = selected.obs["filters_composite"].isin(filtered_composite_indexes)
                selected = selected[condition_filter, :]

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

            # Volcano plot expects specific parameter names (unless we wish to change the options)
            fig = mg.create_volcano_plot(df
                , query_val
                , ref_val
                , pval_threshold
                , [lower_logfc_threshold, upper_logfc_threshold]
                , use_adj_pvals
                )
            mg.modify_volcano_plot(fig, query_val, ref_val)

            if gene_symbols:
                dataset_genes = df['gene_symbol'].unique().tolist()
                normalized_genes_list, _found_genes = mg.normalize_searched_genes(dataset_genes, gene_symbols)
                mg.add_gene_annotations_to_volcano_plot(fig, normalized_genes_list, annotate_nonsignificant)

        elif plot_type == "dotplot":
            df = selected.to_df()

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
            ensm_to_gene = selected.var.to_dict()["gene_symbol"]
            df["gene_symbol"] = df[var_index].map(ensm_to_gene)

            # Percent of all cells in this group where the gene has expression
            percent = lambda row: round(len([num for num in row if num > 0]) / len(row) * 100, 2)
            groupby = ["gene_symbol"]
            groupby.extend(groupby_filters)
            grouped = df.groupby(groupby)
            df = grouped.agg(['mean', 'count', ('percent', percent)]) \
                .fillna(0) \
                .reset_index()

            fig = mg.create_dot_plot(df, groupby_filters, is_log10, title)

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

            fig = mg.create_quadrant_plot(df, control_val, compare1_val, compare2_val)
            # Annotate selected genes
            if gene_symbols:
                genes_not_found, genes_none_none = mg.add_gene_annotations_to_quadrant_plot(fig, normalized_genes_list)
                if genes_not_found:
                    success = 2
                    message += "<li>One or more genes did not pass cutoff filters to be in the plot: {}</li>".format(', '.join(genes_not_found))
                if genes_none_none:
                    success = 2
                    message += "<li>One or more genes had no fold change in both comparisons and will not be annotated: {}</li>".format(', '.join(genes_none_none))


        elif plot_type == "heatmap":
            # Filter genes and slice the adata to get a dataframe
            # with expression and its observation metadata
            df = selected.to_df()

            # Sort alphabetically. Other plot types don't need sorting or sort via groupby
            gene_symbols.sort()
            ensm_to_gene = selected.var.to_dict()["gene_symbol"]
            # Since genes were restricted to one Ensembl ID we can invert keys and vals
            gene_to_ensm = {y:x for x,y in ensm_to_gene.items()}

            # Reorder the dataframe columns based on sorted gene symbols
            dataset_genes = adata.var['gene_symbol'].unique().tolist()
            normalized_genes_list, _found_genes = mg.normalize_searched_genes(dataset_genes, gene_symbols)
            sorted_ensm = map(lambda x: gene_to_ensm[x], normalized_genes_list)
            df = df[sorted_ensm]

            groupby_index = "composite_index"
            groupby_fields = columns
            # Create a composite to groupby
            if filters and matrixplot:
                union_fields = mg.union(list(filters.keys()), clusterbar_fields)
                selected.obs['groupby_composite'] = selected.obs[union_fields].apply(lambda x: ';'.join(map(str,x)), axis=1)
                selected.obs['groupby_composite'] = selected.obs['groupby_composite'].astype('category')
                columns.append("groupby_composite")
                groupby_index = "groupby_composite"
                union_fields.extend([groupby_index, "filters_composite"])   # Preserve filters index for downstream labeling
                groupby_fields = union_fields

            # Remove composite index so that the label is not duplicated.
            for cat in columns:
                df[cat] = selected.obs[cat]

            # Groupby to remove the replicates
            # Ensure the composite index is used as the index for plot labeling
            if matrixplot:
                grouped = df.groupby(groupby_fields)
                df = grouped.agg('mean') \
                    .dropna() \
                    .reset_index() \
                    .set_index(groupby_index)

            # Since this is the new index in the matrixplot, it does not exist as a droppable series
            # These two statements ensure that the current columns are the same if that option is set or not
            groupby_fields.remove(groupby_index)
            if not matrixplot:
                df = df.drop(columns=groupby_index)

            sort_fields = []
            if primary_col:
                sort_fields.append(primary_col)
            if secondary_col and not primary_col == secondary_col:
                sort_fields.append(secondary_col)

            # Sort the dataframe before plotting
            sortby = groupby_fields
            if sort_fields:
                sortby = sort_fields

            sorted_df = df.sort_values(by=sortby)
            df = df.reindex(sorted_df.index.tolist())

            # Drop the obs metadata now that the dataframe is sorted
            # They cannot be in there when the clustergram is made
            # But save it to add back in later
            df_cols = pd.concat([df.pop(cat) for cat in groupby_fields], axis=1)

            # "df" must be obs label for rows and genes for cols only
            fig = mg.create_clustergram(df
                , gene_symbols
                , is_log10
                , cluster_obs
                , cluster_genes
                , flip_axes
                , center_around_zero
                , distance_metric
                )

            # The current clustergram palette used when centering values around zero should be reversed
            fig.data[-1]["reversescale"] = center_around_zero
            if center_around_zero:
                fig.data[-1]["zmid"] = 0
            else:
                # both zmin and zmax are required
                fig.data[-1]["zmin"] = 0
                fig.data[-1]["zmax"] = max(map(max, fig.data[-1]["z"])) # Highest z-value in 2D array

            # Clustergram has a bug where the space where dendrograms should appear is still whitespace
            # Need to adjust the domains of those subplots if no clustering is required
            if not (cluster_genes or cluster_obs):
                fig.layout["xaxis4"]["domain"] = [0,0]
                fig.layout["xaxis5"]["domain"] = [0,0.95]
                fig.layout["yaxis2"]["domain"] = [1,1]
                fig.layout["yaxis5"]["domain"] = [0,1]

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
            ensm_to_gene = selected.var.to_dict()["gene_symbol"]
            df["gene_symbol"] = df[var_index].map(ensm_to_gene)

            violin_func = mg.create_violin_plot
            if stacked_violin:
                violin_func = mg.create_stacked_violin_plot

            fig = violin_func(df
                , groupby_filters
                , is_log10
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
                    "text":legend_title
                }
            )

        # Pop any default height and widths being added
        fig["layout"].pop("height", None)
        fig["layout"].pop("width", None)

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        # NOTE: With volcano plots, the Chrome "devtools" cannot load the JSON response occasionally
        return {
            "success": success
            , "message": message
            , 'gene_symbols': gene_symbols
            , 'plot_json': json.loads(plot_json)
        }
