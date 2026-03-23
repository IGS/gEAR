import argparse

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.rinterface_lib.callbacks as r_cbs
import rpy2.robjects.packages as rpackages
import sys
import mygene
import pandas as pd
import scanpy
import os
import argparse


def silent_handler(s:str) -> None:
    # way to bypass the R stderr output
    pass

def argument_parser():
    parser = argparse.ArgumentParser(usage="%(prog)s -r [RDS Object] -s [Share ID]",add_help=True)
    parser.add_argument('-r', '--rds', required=True, type=str)
    parser.add_argument('-s', '--share-id', required=True, type=str)
    args = vars(parser.parse_args())
    return args

def r_package_installer() -> None:
    utils = rpackages.importr('utils')
    # Install BiocManager if not installed
    if not rpackages.isinstalled('BiocManager'):
        utils.install_packages('BiocManager')
    # Import BiocManager
    BiocManager = importr('BiocManager')
    # Install Seurat, anndataR and rhdf5
    if not rpackages.isinstalled('reticulate'):
        utils.install_packages('reticulate')
    if not rpackages.isinstalled('Seurat'):
        utils.install_packages('Seurat')
    if not rpackages.isinstalled('anndataR'):
        BiocManager.install('anndataR')
    if not rpackages.isinstalled('rhdf5'):
        BiocManager.install('rhdf5')


def r_package_importer(package_name:str):
    """
    Import installed package, if not installed return message
    Input:
        package_name: R package name to import
    Output:
        The R package that was imported or if there's an error the message will be returned
    """
    importErrorMessage = ""
    try:
        pkg = importr(package_name)
        return pkg
    except:
        importErrorMessage += f"{package_name} not installed or can not be imported"
        sys.exit(importErrorMessage)


def seurat_to_anndata(file_path:str,share_name:str):
    """
    file_path: path to rds or rdata file
    share_name: final h5ad string name to be expected (without h5ad)

    return:
        h5ad or name or False
    """
    # Use R's readRDS to load the object.
    # The result is an R object within the Python environment.
    r_seurat_obj = base.readRDS(file_path)
    ro.globalenv['seurat_obj'] = r_seurat_obj
    # Using anndataR write out a converted h5ad
    ro.r('adata <- as_AnnData(seurat_obj)')
    try:
        ro.r(f'write_h5ad(adata, "{share_name}.h5ad")')
        return f'{share_name}.h5ad'
    # In cases where the write fails we will assume the h5ad already exists
    except:
        print(f"h5ad name already exists {share_name}.h5ad")
        return False


def openh5ad(h5ad_name):
    """Just open the supplied h5ad file"""
    adata = scanpy.read_h5ad(h5ad_name)
    return adata

def genes_to_ensembl(adata):
    # We are calling an external API for genes to ensembl mapping
    # Potentially problematic down the road if this shuts down
    genes = adata.var.index.tolist()
    mg = mygene.MyGeneInfo()
    mg_genes = mg.querymany(genes, scopes="symbol", fields="ensembl.gene", species="mouse")
    ensembl_mapping_dict = {}
    for mg_gene in mg_genes:
        gene_name = mg_gene['query']
        if 'ensembl' in mg_gene.keys():
            if isinstance(mg_gene['ensembl'],list):
                # Currently taking first value, not sure of a better way to handle one gene having multiple ensembl IDs
                ensembl_mapping_dict[gene_name] = mg_gene['ensembl'][0]['gene']
            else:
                ensembl_mapping_dict[gene_name] = mg_gene['ensembl']['gene']
    count = 0
    # We still need an ensembl id for the genes that do not actually have them.
    # So here we create a FAKE# for each one so that it can be searchable in gEAR
    for gene in genes:
        if gene not in ensembl_mapping_dict.keys():
            ensembl_mapping_dict[gene] = f"Fake{count}"
            count += 1
    # Overwrite the current adata.var
    adata.var = pd.DataFrame(
        index=list(ensembl_mapping_dict.values()), data={"gene_names": list(ensembl_mapping_dict.keys())}
    )
    return adata


def reduction_to_metadata(adata):
    # Discussion with Carlo and Brian resulted in us determining we would like to
    # take the first 2 values of each reduction
    # PCA in the future, and potentially other reductions may need more
    for reduction in adata.obsm:
        if 'pca' in reduction:
            continue
        elif 'tsne' in reduction.lower() or 'umap' in reduction.lower():
            adata.obs[f'{reduction}_1'] = adata.obsm[reduction][:,0]
            adata.obs[f'{reduction}_2'] = adata.obsm[reduction][:,1]
        elif adata.obsm[reduction].shape[1] > 1:
            for i in range(2):
                adata.obs[f'{reduction}_{i}'] = adata.obsm[reduction][:,i]
    return adata


def main():
    global anndataR, rhdf5, seurat, base
    arguments = argument_parser()
    # Args
    rds_path = arguments['rds']
    share_name = arguments['share_id']
    # remove the R output from writing to stderr
    r_cbs.consolewrite_print = silent_handler
    r_cbs.consolewrite_warnerror = silent_handler
    # Begin R package importing
    r_package_installer()
    base = rpackages.importr('base')
    seurat = r_package_importer('Seurat')
    rhdf5 = r_package_importer("rhdf5")
    anndataR = r_package_importer('anndataR')
    # Take the RDS and output the most basic h5ad
    h5ad_name = seurat_to_anndata(rds_path,share_name)
    # Below are some changes and checks to the h5ad to correctly format for gEAR
    if h5ad_name:
        adata = openh5ad(f'tmp{h5ad_name}')
        adata = genes_to_ensembl(adata)
        adata = reduction_to_metadata(adata)
        adata.write({h5ad_name})
        os.remove(f'tmp{h5ad_name}')


if __name__ == "__main__":
    main()