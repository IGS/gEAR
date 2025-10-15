#!/usr/bin/env python

"""
find_duplicated_dataset_genes.py - Find genes that are duplicated (multiple Ensembl IDs) in every dataset.

By: Shaun Adkins (sadkins@som.umaryland.edu)
"""

import anndata
import os

with open("/tmp/dup_genes_by_dataset.tsv", 'w') as ofh:

    DATASETS_ROOT = "/var/www/datasets"
    for f in os.listdir(DATASETS_ROOT):
        if f.endswith(".h5ad"):
            filepath = os.path.join(DATASETS_ROOT, f)
            adata = anndata.read_h5ad(filepath)
            # keep="False"... mark all dups as true (even first instance)
            # subset="gene_symbol"... only look at this column for duplicates
            try:
                duplicated_genes = adata.var[adata.var.duplicated(keep=False, subset=["gene_symbol"])]
                dataset = f
                ensembl_ids = list(duplicated_genes.index)
                genes = list(duplicated_genes.gene_symbol)

                for i in range(0, len(ensembl_ids)):
                    ofh.write("{}\t{}\t{}\n".format(dataset, ensembl_ids[i], genes[i]))
            except Exception as e:
                print("Error with dataset {}".format(f))
                print(e)