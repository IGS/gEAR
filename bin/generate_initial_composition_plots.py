#!/opt/bin/python3

"""

Generates composition plots for all datasets and stores them statically so the
analyzer interface is more responsive.  Reads datasets from the database and
skips any for which image files are already present.  Does this only for
primary datasets which are of types:

   single-cell-rnaseq
   singlecell-h5ad

"""

import os
import shutil
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

# this is needed so that we don't get TclError failures in the underlying modules
import matplotlib
matplotlib.use('Agg')

import numpy as np

import scanpy as sc
sc.settings.verbosity = 0
sc.settings.autosave = True
sc.settings.figdir = '/tmp'

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    dataset_ids = get_dataset_id_list(cursor)

    for dataset_id in dataset_ids:
        print("Processing dataset ID: {0}".format(dataset_id))
        dataset = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5ad_path = dataset.get_file_path()

        if os.path.exists(h5ad_path):
            print("\tFile found: {0}".format(h5ad_path))
        else:
            print("\tFile not found, skipping: {0}".format(h5ad_path))
            continue

        violin_image_path = h5ad_path.replace('.h5ad', '.prelim_violin.png')
        scatter_image_path = h5ad_path.replace('.h5ad', '.prelim_n_genes.png')

        # Skip if the images are already there
        if os.path.exists(violin_image_path) and os.path.exists(scatter_image_path):
            print("\tImages already found, skipping")
            continue

        ana = geardb.Analysis(dataset_id=dataset_id, type='primary')
        adata = ana.get_adata()

        print("\tGenerating figures")

        # the ".A1" is only necessary, as X is sparse - it transform to a dense array after summing
        try:
            # add the total counts per cell as observations-annotation to adata
            adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1

        except AttributeError:
            # add the total counts per cell as observations-annotation to adata
            adata.obs['n_counts'] = np.sum(adata.X, axis=1)

        sc.pp.filter_cells(adata, min_genes=3)
        sc.pp.filter_genes(adata, min_cells=300)

        axs = sc.pl.violin(adata, ['n_genes', 'n_counts'],
                       jitter=0.4, multi_panel=True, save="_prelim_violin.png")

        ax = sc.pl.scatter(adata, x='n_counts', y='n_genes', save="_prelim_n_genes.png")

        # move files written to tmp
        shutil.move("/tmp/violin_prelim_violin.png", violin_image_path)
        shutil.move("/tmp/scatter_prelim_n_genes.png", scatter_image_path)

    cursor.close()
    cnx.close()

    sys.exit()

def get_dataset_id_list(cursor):
    dataset_ids = list()

    qry = """
          SELECT id
            FROM dataset
           WHERE dtype IN ('single-cell-rnaseq', 'singlecell-h5ad')
             AND has_h5ad = 1
    """
    cursor.execute(qry)
    for row in cursor:
        dataset_ids.append(row[0])

    return dataset_ids

if __name__ == '__main__':
    main()
