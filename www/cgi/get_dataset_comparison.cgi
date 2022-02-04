#!/opt/bin/python3

import cgi, json
import math
import sys
import os
import statistics
import pandas as pd
from itertools import product

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
    dataset_id = form.getvalue('dataset_id')
    condition1 = form.getvalue('condition_x')
    condition2 = form.getvalue('condition_y')
    std_dev_num_cutoff = form.getvalue('std_dev_num_cutoff')
    std_dev_num_cutoff = float(std_dev_num_cutoff) if std_dev_num_cutoff else None
    fold_change_cutoff = form.getvalue('fold_change_cutoff')
    fold_change_cutoff = float(fold_change_cutoff) if fold_change_cutoff else None
    log_transformation = form.getvalue('log_transformation')
    statistical_test = form.getvalue('statistical_test')

    dataset = Dataset(id=dataset_id, has_h5ad=1)
    h5_path = dataset.get_file_path()

    if not os.path.exists(h5_path):
        msg = "No h5 file found for this dataset"
        return_error_response(msg)

    if condition1 == condition2:
        msg = "Selected conditions are identical. Please select different conditions."
        return_error_response(msg)

    # To support some options passed we have to do some stats
    perform_ranking = False
    if statistical_test:
        perform_ranking = True

    adata = sc.read(h5_path)

    # Get all categorical columns
    #cols_x = [col for col in adata.obs.columns if adata.obs[col].dtype.name == 'category']

    # If there are multiple columns, create a composite
    # Only keep columns that had a group selected in the UI
    # NOTE: There may be an edge case where the user selected a group for cond. 1 and nothing for cond. 2 which would leave the column being omitted from the composite
    cols_to_keep = []
    cond1 = json.loads(condition1)
    cols_to_keep.extend([col for col in cond1 if not col in cols_to_keep])
    cond2 = json.loads(condition2)
    cols_to_keep.extend([col for col in cond2 if not col in cols_to_keep])

    # Add new column to combine various groups into a single index
    adata.obs['comparison_composite_index'] = adata.obs[cols_to_keep].apply(lambda x: ';'.join(map(str,x)), axis=1)
    adata.obs['comparison_composite_index'] = adata.obs['comparison_composite_index'].astype('category')
    unique_composite_indexes = adata.obs["comparison_composite_index"].unique()

    # Only want to keep indexes that match chosen filters
    cond1_composite_idx = create_filtered_composite_indexes(cond1, unique_composite_indexes.tolist())
    cond2_composite_idx = create_filtered_composite_indexes(cond2, unique_composite_indexes.tolist())

    # Exit with error if there is no valid composite index mask created
    # Example would be a duplicated obs column where only A was selected in col1 and only B was selected in col2
    if not cond1_composite_idx:
        msg = "The X-axis condition selected combination does not exist for this dataset"
        return_error_response(msg)

    if not cond2_composite_idx:
        msg = "The Y-axis condition selected combination does not exist for this dataset"
        return_error_response(msg)

    # add the composite column for ranked grouping
    if perform_ranking == True:
        sc.pp.filter_cells(adata, min_genes=10)
        sc.pp.filter_genes(adata, min_cells=1)

        # Scanpy.rank_genes_groups can handle multiple groups, but our output is designed for just one gropu
        if len(cond1_composite_idx) > 1:
            msg = "Detected multiple possible conditions for the X-axis condition (the query condition)." + \
                "In order to perform a significance test, please ensure that your set conditions are such so that only 1 possible combination of conditions can be used as the query condition." + \
                "<br />Curated set conditions: {}".format(cond1_composite_idx)
            return_error_response(msg)

        if len(cond2_composite_idx) > 1:
            msg = "Detected multiple possible conditions for the Y-axis condition (the reference condition)." + \
                "In order to perform a significance test, please ensure that your set conditions are such so that only 1 possible combination of conditions can be used as the reference condition." + \
                "<br />Curated set conditions: {}".format(cond2_composite_idx)
            return_error_response(msg)

        try:
            sc.tl.rank_genes_groups(adata, 'comparison_composite_index', use_raw=True, groups=cond1_composite_idx, reference=cond2_composite_idx[0], n_genes=0, rankby_abs=False, copy=False, method=statistical_test, corr_method='benjamini-hochberg', log_transformed=False)
        except Exception as e:
            msg = "scanpy.rank_genes_groups failed.\n{}".format(str(e))
            return_error_response(msg)

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
    condition_x_repls_filter = adata.obs.index.isin(cond1_composite_idx)
    condition_y_repls_filter = adata.obs.index.isin(cond2_composite_idx)
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
    result['condition_x_idx'] = cond1_composite_idx
    result['condition_y_idx'] = cond2_composite_idx
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

def create_filtered_composite_indexes(filters, composite_indexes):
    """Create an index based on the 'comparison_composite_index' column."""
    all_vals = [v for k, v in filters.items()]  # List of lists
    # Remove empty nested lists (list should now have nested lists equal to number of categories with selected groups)
    all_vals = [x for x in all_vals if x]

    # itertools.product returns a combation of every value from every list
    # Essentially  ((x,y) for x in A for y in B)
    filter_combinations = product(*all_vals)
    string_filter_combinations = [";".join(v) for v in filter_combinations]

    # This contains combinations of indexes that may not exist in the dataframe.
    # Use composite indexes from dataframe to return valid filtered indexes
    return intersection(string_filter_combinations, composite_indexes)

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

def intersection(lst1, lst2):
    """Intersection of two lists."""
    return list(set(lst1) & set(lst2))

def return_error_response(msg):
    result = dict()
    result['success'] = 0
    result['error'] = msg
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))
    sys.exit()

if __name__ == '__main__':
    main()
