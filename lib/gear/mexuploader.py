import anndata
import pandas as pd
import numpy as np
import os, sys
import scanpy as sc
sc.settings.verbosity = 0

from gear.datasetuploader import FileType


class MexUploader(FileType):
    """
    Called by datasetuploader.py (factory) when a TAR archive (mexfiles_are_inside.tar) is being uploaded.

    """

    def _read_file(self, filepath):
        """
        Reads in TAR archive that contains the 3 MEX files: matrix.mtx, barcodes.tsv, and genes.tsv.
        Output is directed to the FileUploader object.

        Input
        -----
            filepath - /path/to/your/upload.tar
                     - TAR archive must end with file extension '.tar'
                     - The must contain 3 files by the exact names as below:
                         matrix.mtx
                         barcodes.tsv
                         genes.tsv

        Output
        ------
            Output is assigned to the DatasetUploader object:
                DatasetUploader.adata = adata
                DatasetUploader.originalFile = filepath

            'adata' is an AnnData object where data of each MEX file is assigned to
            the AnnData object:
                AnnData.X = matrix.mtx
                AnnData.obs = barcodes.tsv
                AnnData.var = genes.tsv

            'filepath' is the file path of the original file

        """
        import tarfile

        #Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        # Open the archive and extract each file into the tmp directory
        with tarfile.open(filepath) as tf:
            for entry in tf:

                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

                # Read each file (.mtx as adata.X and .tsv as pandas dataframes)
                if entry.name == 'matrix.mtx' or os.path.basename(filepath)== 'matrix.mtx':
                    adata = sc.read(filepath, cache=False).transpose()
                elif entry.name == 'barcodes.tsv' or os.path.basename(filepath)== 'barcodes.tsv':
                    obs = pd.read_csv(filepath, sep='\t', index_col=0, header=None, names=['observations'])
                elif entry.name == 'genes.tsv' or os.path.basename(filepath)=='genes.tsv':
                    var = pd.read_csv(filepath, sep='\t', index_col=0, header=None, names=['genes', 'gene_symbol'])
                else:
                    raise Exception("Unexpected file name: '{0}'. Excepted 'matrix.mtx', 'genes.tsv', or 'barcodes.tsv'.".format(entry.name))

        # Assign genes and observations to AnnData object
        adata.var = var
        adata.obs = obs

        # Apply AnnData obj and filepath to uploader obj
        self.adata = adata
        self.originalFile = filepath
        return self


    def get_plot_preview_expression(self):
        #TODO: Currently no plot types support this data type.
        pass


    def _write_to_h5ad(self, filepath=None):
        #TODO this might not be used. It's only a anndata.write command
        """
        Write AnnData object to file in h5ad format.
        """

        if self.adata is None:
            raise Exception("No AnnData object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")


        # write AnnData object to file
        try:
            self.adata.write(filename=filepath)
        except Exception as err:
            raise Exception("Error occurred while writing to file:\n", err)

        return self
