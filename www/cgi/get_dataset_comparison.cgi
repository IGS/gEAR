#!/opt/bin/python3

import cgi
import json
import math
import os
import statistics
import sys
import traceback

import pandas as pd
import scanpy as sc

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import get_analysis

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

sc.settings.verbosity = 0

def main():
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    filters = form.getvalue('obs_filters', "")    # Dict of lists

    compare_key = form.getvalue('compare_key')
    x_compare = form.getvalue('condition_x')    # list of conditions
    y_compare = form.getvalue('condition_y')
    std_dev_num_cutoff = form.getvalue('std_dev_num_cutoff')
    std_dev_num_cutoff = float(std_dev_num_cutoff) if std_dev_num_cutoff else None
    fold_change_cutoff = form.getvalue('fold_change_cutoff')
    fold_change_cutoff = float(fold_change_cutoff) if fold_change_cutoff else None
    log_transformation = form.getvalue('log_transformation')
    statistical_test = form.getvalue('statistical_test')

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        return {
            "success": -1,
            'message': "No dataset found with that ID"
        }
    is_spatial = ds.dtype == "spatial"

    if not x_compare or not y_compare:
        msg = "Please select a condition for both the X and Y axis"
        return_error_response(msg)

    if x_compare == y_compare:
        msg = "Selected conditions are identical. Please select different conditions."
        return_error_response(msg)

    try:
        ana = get_analysis(None, dataset_id, None, is_spatial=is_spatial)
    except Exception:
        traceback.print_exc()
        return_error_response("Analysis for this dataset is unavailable.")

    try:
            args = {}
            if is_spatial:
                args['include_images'] = False
            adata = ana.get_adata(**args)
    except Exception:
        traceback.print_exc()
        return_error_response("Could not create dataset object using analysis.")

    # To support some options passed we have to do some stats
    perform_ranking = False
    if statistical_test:
        perform_ranking = True

    filters = json.loads(filters)
    # Filter by obs filters
    if filters:
        for col, values in filters.items():
            selected_filter = adata.obs[col].isin(values)
            adata = adata[selected_filter, :]

    x_compare = json.loads(x_compare)
    y_compare = json.loads(y_compare)

    # Error if any condition in x matches any condition in y
    intersection_conditions = intersection(x_compare, y_compare)
    if intersection_conditions:
        msg = f"The follwing conditions were found in both X and Y: {intersection_conditions}. Please select unique conditions."
        return_error_response(msg)

    # add the composite column for ranked grouping
    if perform_ranking:
        try:
            # TODO: Had the tool crash here and apache restart resolved it.  Need to investigate.
            sc.pp.filter_cells(adata, min_genes=10)
            sc.pp.filter_genes(adata, min_cells=1)
        except Exception as e:
            msg = "scanpy.pp.filter_cells or scanpy.pp.filter_genes failed.\n{}".format(str(e))
            return_error_response(msg)

    condition_x_repls_filter = adata.obs[compare_key].isin(x_compare)
    condition_y_repls_filter = adata.obs[compare_key].isin(y_compare)

    if perform_ranking:
        # Assign categorical values if sample is in x_compare or y_compare
        adata.obs['compare'] = 'neither'
        adata.obs.loc[condition_x_repls_filter, 'compare'] = 'x'
        adata.obs.loc[condition_y_repls_filter, 'compare'] = 'y'

        try:
            sc.tl.rank_genes_groups(adata, "compare", groups=["x"], reference="y", n_genes=0, rankby_abs=False, copy=False, method=statistical_test, corr_method='benjamini-hochberg', log_transformed=False)
        except Exception as e:
            msg = "scanpy.rank_genes_groups failed.\n{}".format(str(e))
            return_error_response(msg)

    # AnnData does not yet allow slices on both rows and columns
    # with boolean indices, so first filter genes and then grab
    # value with corresponding condition index from obs

    # Match the condition (ex. A1;ADULT;F) to all
    # the rows in obs, and use this to aggregate the mean
    # SAdkins - This works with and without replicates

    adata_x_subset = adata[condition_x_repls_filter, :]
    adata_y_subset = adata[condition_y_repls_filter, :]

    df_x = pd.DataFrame({
        # adata.X ends up being 2 dimensional array with gene's values going down a column.
        # We tranpose so a gene's replicate values are in a list and then we take the average
        'e1_raw': [float(replicate_values.mean()) for replicate_values in adata_x_subset.X.transpose()], # type: ignore
        'e2_raw': [float(replicate_values.mean()) for replicate_values in adata_y_subset.X.transpose()], # type: ignore
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

    if std_dev_num_cutoff and std_dev_num_cutoff > 0:
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

    if fold_change_cutoff and fold_change_cutoff > 0:
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
    log_base = None
    if log_transformation == "2":
        log_base = 2
    elif log_transformation == "10":
        log_base = 10

    if log_base:
        for (e1_raw, e2_raw) in result['values']:
            transformed_e1 = get_log(e1_raw, log_base)
            transformed_e2 = get_log(e2_raw, log_base)

            if transformed_e1 is not None and transformed_e2 is not None:
                filtered_x.append(transformed_e1)
                filtered_y.append(transformed_e2)
                filtered_values.append([transformed_e1,transformed_e2])

        result['values'] = filtered_values
        result['x'] = filtered_x
        result['y'] = filtered_y

    result['fold_change_std_dev'] = "{0:.2f}".format(fold_change_std_dev)
    result["compare_key"] = compare_key
    result['condition_x'] = x_compare
    result['condition_y'] = y_compare
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

def fold_change(x, y):
    if x >= y:
        if y == 0:
            return x
        return x / y
    if x == 0:
        return y
    return y / x

def get_log(val, base):
    if val == 0:
        return 0
    try:
        return math.log(val, base)
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
