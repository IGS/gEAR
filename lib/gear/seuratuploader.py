import argparse
import os
import sys
import traceback
from typing import Callable, Optional

import pandas as pd
import rpy2.rinterface_lib.callbacks as r_cbs
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import scanpy
from gear.utils.gene_mapping import map_gene_symbols_via_mygene
from rpy2.robjects.packages import importr


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

# TODO: Recently switched to the pak installer for the docker image, should consider switching this too
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
    if not rpackages.isinstalled('Signac'):
        ro.r("setRepositories(ind = 1:3)") # needed to automatically install Bioconductor dependencies
        utils.install_packages('Signac')
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
    try:
        pkg = importr(package_name)
        return pkg
    except Exception:
        importErrorMessage = f"{package_name} not installed or can not be imported"
        traceback.print_exc()
        raise ImportError(importErrorMessage)


def seurat_to_anndata(
    file_path: str,
    share_name: str,
    output_dir: str = ".",
    progress_callback: Optional[Callable[[int, str], None]] = None,
):
    """
    file_path: path to rds or rdata file
    share_name: final h5ad string name to be expected (without h5ad)
    output_dir: directory to write the temporary h5ad file into
    progress_callback: optional (progress: int, message: str) callback invoked before each of
        the three R steps to log progress, which is relayed to the user

    return:
        absolute path to tmp h5ad, or False on failure
    """

    def _report(pct: int, msg: str) -> None:
        # Always print to stderr (timestamped by journald under systemd) even without a
        # callback, so this is visible in the consumer's own logs either way.
        print(f"INFO: {msg}", file=sys.stderr, flush=True)
        if progress_callback is not None:
            progress_callback(pct, msg)

    # Suppress R console output and ensure required packages are loaded,
    # since this function may be called as a module in cgi script (not via main()).
    r_cbs.consolewrite_print = silent_handler
    r_cbs.consolewrite_warnerror = silent_handler
    # Import required R packages
    base = rpackages.importr('base')
    try:
        r_package_importer('Seurat')
        r_package_importer('rhdf5')
        r_package_importer('anndataR')
        r_package_importer('Signac')    # For some extra stuff Carlo adds in
    except ImportError:
        raise
    # Use R's readRDS to load the object.
    # The result is an R object within the Python environment.
    _report(6, "Reading Seurat/RDS object into R (large files can take a while)...")
    try:
        r_seurat_obj = base.readRDS(file_path)
    except Exception as e:
        print(f"Error reading RDS file using readRDS: {e}", file=sys.stderr)
        print("ERROR (readRDS): Perhaps the file is not a valid RDS file (from saveRDS), or the path is incorrect.", file=sys.stderr)
        raise ValueError("Error reading RDS file")
    ro.globalenv['seurat_obj'] = r_seurat_obj

    # Discover the reductions present on this object (e.g. "pca", "umap", "tsne")
    reduction_names = list(ro.r('Reductions(seurat_obj)'))

    # Build obsm_mapping using each reduction's exact name as both key and value
    obsm_mapping = ro.ListVector({name: name for name in reduction_names})
    ro.globalenv['obsm_mapping'] = obsm_mapping

    # Convert directly to a file-backed HDF5AnnData object (output_class = "HDF5AnnData")
    # instead of anndataR's default InMemoryAnnData. This streams the conversion straight to
    # the target .h5ad file rather than first building a second full in-memory representation
    # of the already-loaded R Seurat object, then a third copy when using write_h5ad()
    output_path = os.path.join(output_dir, f'tmp_{share_name}.h5ad')
    _report(10, "Converting Seurat object to AnnData (this can take a while for large datasets)...")
    try:
        ro.r(
            'adata <- as_AnnData(seurat_obj, obsm_mapping = obsm_mapping, '
            f'output_class = "HDF5AnnData", file = "{output_path}", mode = "w")'
        )
        # Close to release lock on file.
        ro.r('adata$close()')
    except Exception as e:
        print(f"Error converting Seurat object to AnnData: {e}", file=sys.stderr)
        raise ValueError("Error converting Seurat object to AnnData")

    _report(14, "Finished writing H5AD file.")
    return output_path

def openh5ad(h5ad_name):
    """Just open the supplied h5ad file"""
    adata = scanpy.read_h5ad(h5ad_name)
    return adata


def genes_to_ensembl(adata, taxid=None):
    if taxid is None:
        return None

    genes = adata.var.index.tolist()
    ensembl_mapping_dict = map_gene_symbols_via_mygene(genes, taxid, verbose=True)

    count = 0
    for gene in genes:
        if gene not in ensembl_mapping_dict:
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
