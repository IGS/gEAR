#!/opt/bin/python3

import cgi, json
import math
import sys
import os
import statistics
import pandas as pd

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from geardb import Dataset

import scanpy as sc
sc.settings.verbosity = 0

def main():
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset1_id')
    condition1 = form.getvalue('dataset1_condition')
    condition2 = form.getvalue('dataset2_condition')
    std_dev_num_cutoff = float(form.getvalue('std_dev_num_cutoff'))
    fold_change_cutoff = float(form.getvalue('fold_change_cutoff'))
    log_transformation = form.getvalue('log_transformation')
    statistical_test = form.getvalue('statistical_test')

    dataset = Dataset(id=dataset_id, has_h5ad=1)
    h5_path = dataset.get_file_path()

    # This needs to match the list in get_condition_list.cgi
    # this is for the NeMO demo do handle the datasets overloaded with columns
    skip_columns = [
        'Amp_Date',
        'Amp_Name',
        'Amp_PCR_cyles',
        'ATAC_cluster_label',
        'ATAC_cluster_color',
        'Cell_Capture',
        'DNAm_cluster_color',
        'DNAm_cluster_label',
        'Donor',
        'Gender',
        'Lib_Cells',
        'Lib_Date',
        'Lib_Name',
        'Lib_PCR_cycles',
        'Lib_PassFail',
        'Lib_type',
        'Live_Cells',
        'Live_percent',
        'Mean_Reads_perCell',
        'Median_Genes_perCell',
        'Median_UMI_perCell',
        'Region',
        'Replicate_Lib',
        'RNA_cluster_color',
        'RNA_cluster_id',
        'RNA_cluster_label',
        'Saturation',
        'Seq_batch',
        'Total_Cells',
        'age_in_weeks',
        'aggr_num',
        'barcode',
        'cellType',
        'celltype_a',
        'celltype_b',
        'celltype_colors',
        'class_id',
        'class_label',
        'class_color',
        'cluster_id',
        'cluster_label',
        'cluster_color',
        'cross_species_cluster',
        'cross_species_cluster_color',
        'cross_species_cluster_id',
        'cross_species_cluster_label',
        'donor_id',
        'donor_sex',
        'donor_sex_color',
        'donor_sex_id',
        'donor_sex_label',
        'doublet.score',
        'exp_component_name',
        'gene.counts',
        'genes_detected',
        'genes_detected_color',
        'genes_detected_id',
        'genes_detected_label',
        'library_id',
        'louvain',
        'mapped_reads',
        'method',
        'nonconf_mapped_reads',
        'sample_id',
        'species_color',
        'subclass_color',
        'subclass_id',
        #'subclass_label',  # Leave this as displayed for NeMO datasets
        'tSNE_1',
        'tSNE_2',
        'total.reads',
        'total_UMIs',
        'total_UMIs_color',
        'total_UMIs_id',
        'total_UMIs_label',
        'tsne1_combined',
        'tsne2_combined',
        'tsne_1',
        'tsne_2',
        'tsne1',
        'tsne2',
        'tSNE1',
        'tSNE2',
        'tsne_x',
        'tsne_y',
        'tube_barcode',
        'UMAP1',
        'UMAP2',
        'umap1',
        'umap2',
        'umap_1',
        'umap_2',
        'uMAP_1',
        'uMAP_2',
        'umap1_combined',
        'umap2_combined',
        'umi.counts',
        'unmapped_reads',
        # NeMOAnalytics columns to skip
        "RNANumber",
        "BRNumDigit",
        "RIN",
        "AgeYR",
        "TotalNumReads",
        "TotalNumMapped",
        "log10(AgeYR+0.8)",
        "color",
        "AgeRND",
        "SampleID",
        "BioRep",
        "TechRep",
        "color",
        "X"
        ]

    if not os.path.exists(h5_path):
        result = dict()
        result['success'] = 0
        result['error'] = "No h5 file found for this dataset"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    if condition1 == condition2:
        result = dict()
        result['success'] = 0
        result['error'] = "Selected conditions are identical. Please select different conditions."
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    # To support some options passed we have to do some stats
    perform_ranking = False
    if statistical_test:
        perform_ranking = True

    adata = sc.read(h5_path)

    # Add the composite column
    cols_x = adata.obs.columns.tolist()

    # If there are multiple columns, create a composite
    if len(cols_x) > 0:
        if "replicate" in cols_x:
            cols_x.remove("replicate")
        if 'BioRep' in cols_x:
            cols_x.remove('BioRep')
        if 'TechRep' in cols_x:
            cols_x.remove('TechRep')

        for unwanted_col in skip_columns:
            if unwanted_col in cols_x:
                cols_x.remove(unwanted_col)

        composite_index = adata.obs[cols_x].apply(lambda x: ';'.join(map(str,x)), axis=1)
        adata.obs['comparison_composite_index'] = composite_index.tolist()

        # add the composite column for ranked grouping
        if perform_ranking == True:
            sc.pp.filter_cells(adata, min_genes=10)
            sc.pp.filter_genes(adata, min_cells=1)
            sc.tl.rank_genes_groups(adata, 'comparison_composite_index', use_raw=True, groups=[condition1], reference=condition2, n_genes=0, rankby_abs=False, copy=False, method=statistical_test, corr_method='benjamini-hochberg', log_transformed=False)

        # Set the index as the composite column so we can more easily match the dataset condition being searched.
        # Our condition has the convention of <column_name>;<column_name>;...
        # In order to index our dataframe on our condition, we need to map over each
        # row and set its levels with the same convention. For example, the condition
        # cochlea;GFP+;E16 would be set as an index so it can be sliced.
        adata.obs = adata.obs.set_index('comparison_composite_index')

    # AnnData does not yet allow slices on both rows and columns
    # with boolean indices, so first filter genes and then grab
    # value with corresponding condition index from obs

    # Match the condition (ex. A1;ADULT;F) to all
    # the rows in obs, and use this to aggregate the mean
    # SAdkins - This works with and without replicates
    condition_x_repls_filter = adata.obs.index == condition1
    condition_y_repls_filter = adata.obs.index == condition2
    adata_x_subset = adata[condition_x_repls_filter, :]
    adata_y_subset = adata[condition_y_repls_filter, :]

    df_x = pd.DataFrame({
        # adata.X ends up being 2 dimensional array with gene's values going down a column.
        # We tranpose so a gene's replicate values are in a list and then we take the average
        'e1_raw': [float(replicate_values.mean()) for replicate_values in adata_x_subset.X.transpose()],
        'e2_raw': [float(replicate_values.mean()) for replicate_values in adata_y_subset.X.transpose()],
        'gene_sym': adata_x_subset.var.gene_symbol
    })

    if perform_ranking:
        df_x['pvals_adj'] = adata.uns['rank_genes_groups']['pvals_adj']

    # Adding e1_condition here instead of inside DataFrame
    # constructor so we can spread the single condition string
    # across all rows in the column.
    df_x['e1_condition'] = condition1
    df_x['e2_condition'] = condition2

    result = {
               'fold_change_std_dev': None,
               'gene_ids': list(),
               'symbols': list(),
               'success': 1,
               'values': list(),
               'pvals_adj': list(),
               'x': list(),
               'y': list(),
               'fold_changes': list()
    }

    # There's too much data to show it all, else plotting fails.  So let's filter
    #  and skip those whose differences fall within N standard deviations of the mean.
    for index, row in df_x.iterrows():
        result['fold_changes'].append(fold_change(row['e1_raw'], row['e2_raw']))
        result['values'].append([row['e1_raw'], row['e2_raw']])
        result['x'].append(row['e1_raw'])
        result['y'].append(row['e2_raw'])
        result['gene_ids'].append(index)
        result['symbols'].append(row['gene_sym'])

        if perform_ranking:
            result['pvals_adj'].append(row['pvals_adj'])

    fold_change_std_dev = statistics.stdev(result['fold_changes'])

    filtered_values = list()
    filtered_pvals_adj = list()
    filtered_x = list()
    filtered_y = list()
    filtered_gene_ids = list()
    filtered_symbols = list()
    filtered_fold_changes = list()

    if std_dev_num_cutoff > 0:
        cutoff_diff = fold_change_std_dev * std_dev_num_cutoff
        idx = 0

        for (e1_raw, e2_raw) in result['values']:
            if fold_change(e1_raw, e2_raw) > cutoff_diff:
                filtered_values.append([e1_raw, e2_raw])
                filtered_x.append(e1_raw)
                filtered_y.append(e2_raw)
                filtered_gene_ids.append(result['gene_ids'][idx])
                filtered_symbols.append(result['symbols'][idx])
                filtered_fold_changes.append(fold_change(e1_raw, e2_raw))

                if perform_ranking:
                    filtered_pvals_adj.append(result['pvals_adj'][idx])

            idx += 1

        result['values'] = filtered_values
        result['pvals_adj'] = filtered_pvals_adj
        result['gene_ids'] = filtered_gene_ids
        result['symbols'] = filtered_symbols
        result['x'] = filtered_x
        result['y'] = filtered_y
        result['fold_changes'] = filtered_fold_changes
        filtered_values = list()
        filtered_pvals_adj = list()
        filtered_x = list()
        filtered_y = list()
        filtered_gene_ids = list()
        filtered_symbols = list()
        filtered_fold_changes = list()

    if fold_change_cutoff > 0:
        idx = 0

        for (e1_raw, e2_raw) in result['values']:
            if fold_change(e1_raw, e2_raw) >= fold_change_cutoff:
                filtered_values.append([e1_raw, e2_raw])
                filtered_x.append(e1_raw)
                filtered_y.append(e2_raw)
                filtered_gene_ids.append(result['gene_ids'][idx])
                filtered_symbols.append(result['symbols'][idx])
                filtered_fold_changes.append(fold_change(e1_raw, e2_raw))

                if perform_ranking:
                    filtered_pvals_adj.append(result['pvals_adj'][idx])

            idx += 1

        result['values'] = filtered_values
        result['pvals_adj'] = filtered_pvals_adj
        result['x'] = filtered_x
        result['y'] = filtered_y
        result['gene_ids'] = filtered_gene_ids
        result['symbols'] = filtered_symbols
        result['fold_changes'] = filtered_fold_changes
        filtered_values = list()
        filtered_x = list()
        filtered_y = list()

    # Is there a transformation to apply?
    if log_transformation == "2":
        for (e1_raw, e2_raw) in result['values']:
            transformed_e1 = get_log2(e1_raw)
            transformed_e2 = get_log2(e2_raw)

            if transformed_e1 is not None and transformed_e2 is not None:
                filtered_x.append(transformed_e1)
                filtered_y.append(transformed_e2)
                filtered_values.append([transformed_e1,transformed_e2])

        result['values'] = filtered_values
        result['x'] = filtered_x
        result['y'] = filtered_y

    elif log_transformation == "10":
        for (e1_raw, e2_raw) in result['values']:
            transformed_e1 = get_log10(e1_raw)
            transformed_e2 = get_log10(e2_raw)

            if transformed_e1 is not None and transformed_e2 is not None:
                filtered_x.append(transformed_e1)
                filtered_y.append(transformed_e2)
                filtered_values.append([transformed_e1,transformed_e2])

        result['values'] = filtered_values
        result['x'] = filtered_x
        result['y'] = filtered_y

    result['fold_change_std_dev'] = "{0:.2f}".format(fold_change_std_dev)
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

def fold_change(x, y):
    if x >= y:
        if y == 0:
            return x
        else:
            return x / y
    else:
        if x == 0:
            return y
        else:
            return y / x

def get_log10(val):
    if val == 0:
        return 0
    else:
        try:
            return math.log10(val)
        except ValueError:
            return None

def get_log2(val):
    if val == 0:
        return 0
    else:
        try:
            return math.log2(val)
        except ValueError:
            return None

if __name__ == '__main__':
    main()
