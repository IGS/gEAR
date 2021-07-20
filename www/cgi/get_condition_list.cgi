#!/opt/bin/python3

"""
For a given dataset_id, returns a unique list of the conditions within that dataset.

These are omitted:

   replicate, barcode, tSNE_1, tSNE_2

Data structure returned (example):

{
   has_replicates: 1,
   conditions: [ {'class_label': 'A1_ADULT_F_rep1', 'formatted_class_label: 'A1, ADULT, F, rep1'},  ...  ]
}

"""

import cgi
import json
import os
import sys
from datetime import datetime

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from geardb import Dataset
import scanpy as sc

# The select boxes are auto-generated based on aggregates of the obs columns.  Having
#  too many in the dataset though kills the interface.  Put a limit here and then
#  just return an error if there are too many (report to user why)
MAX_ALLOWED_CONDITIONAL_COLUMNS = 5

def main():
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    dataset = Dataset(id=dataset_id, has_h5ad=1)

    h5_path = dataset.get_file_path()
    result = {'has_replicates': 0}

    # These columns should not be included in any composite indexes or labels
    # This needs to match the list in get_dataset_comparison.cgi
    #skip_columns = ['barcode', 'tSNE_1', 'tSNE_2', 'tsne_1', 'tsne_2', 'tsne_x', 'tsne_y']

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
        'cell_type_colors',
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
        'UMAP_1',
        'UMAP_2',
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
        result['success'] = 0
        result['error'] = f"No h5 file found for {dataset_id}."
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    adata = sc.read(h5_path)

    if adata.obs.empty and len(adata.obs.index) > 50:
        result['success'] = -1
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    df = adata.obs
    columns = df.columns.tolist()
    if 'replicate' in columns:
        columns.remove('replicate')
        result['has_replicates'] = 1
    if 'BioRep' in columns:
        columns.remove('BioRep')
        result['has_replicates'] = 1
    if 'TechRep' in columns:
        columns.remove('TechRep')
        result['has_replicates'] = 1

    # omit any skipped columns
    for unwanted_col in skip_columns:
        if unwanted_col in columns:
            columns.remove(unwanted_col)

    if len(columns) > MAX_ALLOWED_CONDITIONAL_COLUMNS:
        result['success'] = 0
        result['error'] = "This dataset has too many observation columns to be used in the comparison tool.  We are working to improve this so you can simply choose column groupings in the future."
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    # some datasets may have an index but no
    # columns specificed
    if len(columns) > 0:
        df['index'] = df[columns].apply(lambda x: ';'.join(map(str,x)), axis=1)
        df = df.set_index('index')

    conditions = [
        dict(
            class_label=index,
            # Todo:
            # We are using underscores to separate attributes, e.g., cellType_Condition,
            # to aid us in splitting the condition and indexing appropriately.
            # This should be formatted with commas rather than underscores when printing
            # inside the select box, however, some values have underscores in their names and would get
            # replaced. For example, Mes_Cells_2, is a cell type level. Need to possibly replace all
            # underscores with whitespace before joining columns on the underscore.
            formatted_class_label=index.replace(';', ", ")
        )
        for index in df.index.unique().tolist()
    ]

    result['conditions'] = conditions

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
