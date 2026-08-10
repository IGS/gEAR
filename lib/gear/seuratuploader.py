import argparse
import os
import sys

import mygene
import pandas as pd
import rpy2.rinterface_lib.callbacks as r_cbs
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import scanpy
from rpy2.robjects.packages import importr
import time


def silent_handler(s:str) -> None:
    # way to bypass the R stderr output
    pass

def argument_parser():
    parser = argparse.ArgumentParser(usage="%(prog)s -r [RDS Object] -s [Share ID]",add_help=True)
    parser.add_argument('-r', '--rds', required=True, type=str)
    parser.add_argument('-s', '--share-id', required=True, type=str)
    parser.add_argument('-t','--tax-id',required=False,type=str)
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
    except Exception:
        importErrorMessage += f"{package_name} not installed or can not be imported"
        raise ImportError(importErrorMessage)


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
    try:
        r_seurat_obj = base.readRDS(file_path)
    except Exception as e:
        print(f"Error reading RDS file using readRDS: {e}", file=sys.stderr)
        print("ERROR (readRDS): Perhaps the file is not a valid RDS file (from saveRDS), or the path is incorrect.", file=sys.stderr)
        raise ValueError("Error reading RDS file")
    ro.globalenv['seurat_obj'] = r_seurat_obj
    # Using anndataR write out a converted h5ad
    ro.r('adata <- as_AnnData(seurat_obj)')
    output_path = os.path.join(output_dir, f'tmp_{share_name}.h5ad')
    try:
        ro.r(f'write_h5ad(adata, "{output_path}")')
        return output_path
    # In cases where the write fails we will assume the h5ad already exists
    except Exception:
        print(f"h5ad name already exists {output_path}", file=sys.stderr)
        raise ValueError("Error writing h5ad file to output path")


def openh5ad(h5ad_name):
    """Just open the supplied h5ad file"""
    adata = scanpy.read_h5ad(h5ad_name)
    return adata


def genes_to_ensembl(adata, taxid=None):
    if taxid is None:
        return None

    genes = adata.var.index.tolist()
    max_retries = 5
    base_delay = 2  # seconds

    for attempt in range(max_retries):
        try:
            mg = mygene.MyGeneInfo()
            mg_genes = mg.querymany(genes, scopes="symbol", fields="ensembl.gene", species=f"{taxid}")
            break  # success — exit retry loop
        except Exception as e:
            is_server_error = "500" in str(e) or "Internal Server Error" in str(e).lower()
            if is_server_error and attempt < max_retries - 1:
                delay = base_delay ** (attempt + 1)  # 2, 4, 8, 16, 32 seconds
                print(
                    f"MyGene API returned a 500 error (attempt {attempt + 1}/{max_retries}). "
                    f"Retrying in {delay}s...",
                )
                time.sleep(delay)
            else:
                # Non-500 error, or all retries exhausted
                print(f"Error occurred while querying MyGene: {e}")
                raise

    ensembl_mapping_dict = {}
    for mg_gene in mg_genes:
        gene_name = mg_gene['query']
        if 'ensembl' in mg_gene.keys():
            if isinstance(mg_gene['ensembl'], list):
                ensembl_mapping_dict[gene_name] = mg_gene['ensembl'][0]['gene']
            else:
                ensembl_mapping_dict[gene_name] = mg_gene['ensembl']['gene']

    count = 0
    for gene in genes:
        if gene not in ensembl_mapping_dict.keys():
            ensembl_mapping_dict[gene] = f"Fake{count}"
            count += 1

    adata.var = pd.DataFrame(
        index=list(ensembl_mapping_dict.values()), data={"gene_symbol": list(ensembl_mapping_dict.keys())}
    )
    return adata


def reduction_to_metadata(adata):
    # Discussion with Carlo and Brian resulted in us determining we would like to
    # take the first 2 values of each reduction
    # PCA in the future, and potentially other reductions may need more
    try:
        for reduction in adata.obsm:
            if adata.obsm[reduction].shape[1] > 1:
                for i in range(2):
                    adata.obs[f'{reduction}_{i+1}'] = adata.obsm[reduction][:,i]
    except Exception as e:
        print(f"Error processing reductions: {e}", file=sys.stderr)
        raise ValueError("Error processing reductions in AnnData object")
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
    tax_id = arguments.get('tax_id',None)
    r_package_installer()
    # Take the RDS and output the most basic h5ad
    h5ad_name = seurat_to_anndata(rds_path,share_name)
    # Below are some changes and checks to the h5ad to correctly format for gEAR
    if h5ad_name:
        print(h5ad_name)
        adata = openh5ad(f'{h5ad_name}')
        if tax_id is None:
            raise ValueError("TaxID not supplied")
        adata = genes_to_ensembl(adata,taxid=tax_id)
        adata = reduction_to_metadata(adata)
        if adata is None:
            raise ValueError("Anndata object is None after gene conversion")
        adata.write(str(h5ad_name.replace('tmp_','').replace('./','')))
        os.remove(f'{h5ad_name}')


if __name__ == "__main__":
    main()
