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



def seurat_to_anndata(file_path: str, share_name: str, output_dir: str = "."):
    """
    file_path: path to rds or rdata file
    share_name: final h5ad string name to be expected (without h5ad)
    output_dir: directory to write the temporary h5ad file into

    return:
        absolute path to tmp h5ad, or False on failure
    """
    # Suppress R console output and ensure required packages are loaded,
    # since this function may be called as a module in cgi script (not via main()).
    r_cbs.consolewrite_print = silent_handler
    r_cbs.consolewrite_warnerror = silent_handler
    # Import required R packages
    base = rpackages.importr('base')
    r_package_importer('Seurat')
    r_package_importer('rhdf5')
    r_package_importer('anndataR')
    # Use R's readRDS to load the object.
    # The result is an R object within the Python environment.
    r_seurat_obj = base.readRDS(file_path)
    ro.globalenv['seurat_obj'] = r_seurat_obj
    # Using anndataR write out a converted h5ad
    ro.r('adata <- as_AnnData(seurat_obj)')
    output_path = os.path.join(output_dir, f'tmp_{share_name}.h5ad')
    try:
        ro.r(f'write_h5ad(adata, "{output_path}")')
        return output_path
    # In cases where the write fails we will assume the h5ad already exists
    except:
        print(f"h5ad name already exists {output_path}")
        return False


def openh5ad(h5ad_name):
    """Just open the supplied h5ad file"""
    adata = scanpy.read_h5ad(h5ad_name)
    return adata

def genes_to_ensembl(adata, taxid=None):
    # We are calling an external API for genes to ensembl mapping
    # Potentially problematic down the road if this shuts down
    if taxid is None:
        return None
    genes = adata.var.index.tolist()
    mg = mygene.MyGeneInfo()
    mg_genes = mg.querymany(genes, scopes="symbol", fields="ensembl.gene", species=f"{taxid}")
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
        if adata.obsm[reduction].shape[1] > 1:
            for i in range(2):
                adata.obs[f'{reduction}_{i+1}'] = adata.obsm[reduction][:,i]
    return adata


def layer_to_X(adata, layer_name):
    # Possibility for Seurat -> Anndata conversion doesn not create the X matrix.
    # Use adata.layers['data'] as X
    adata.X = adata.layers[layer_name]
    return adata

def main():
    arguments = argument_parser()
    # Args
    rds_path = arguments['rds']
    share_name = arguments['share_id']
    r_package_installer()
    # Take the RDS and output the most basic h5ad
    h5ad_name = seurat_to_anndata(rds_path,share_name)
    # Below are some changes and checks to the h5ad to correctly format for gEAR
    if h5ad_name:
        adata = openh5ad(f'tmp_{h5ad_name}')
        adata = genes_to_ensembl(adata)
        if adata is None:
            sys.exit("TaxID not supplied")
        adata = reduction_to_metadata(adata)
        adata.write({h5ad_name.replace('tmp_','')})
        os.remove(f'tmp_{h5ad_name}')


if __name__ == "__main__":
    main()
