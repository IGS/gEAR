import os,sys
import anndata
import pandas as pd
import scanpy as sc


class DataArchive:
    """
    This class handles:
    1. Get dataset type (mex or 3tab)
    2. Create h5ad and return
    3. Write h5ad to specfic directory
    """

    def __init__(self, base_dir=None, extracted_file_paths=None, format=None, tarball_path=None):
        # If the archive is a tarball, this is the source path for it
        self.tarball_path = tarball_path

        # If the archive is extracted into a directory, this should be
        #  that directory base, directly within which all files are found.
        self.base_dir = base_dir

        # The full path to each file in the archive
        self.extracted_file_paths = list()

        # Can be explicit or derived (mex, 3tab, etc.)
        self.format = format


    def get_archive_type(self, data_path = None):
        """
        Determines type of dataset based on file names and extensions.
        Input
        -----
            data_path - path/to/extracted/data
        Output
        -----
            return "mex" or "3tab"
        """
        archive_filenames = os.listdir(data_path)
        archive_files=''.join(archive_filenames)
        data_type = None
        if 'expression.tab' in archive_files and 'genes.tab' in archive_files and 'observations.tab' in archive_files:
            data_type = "3tab"
        if 'matrix.mtx' in archive_files and 'barcodes.tsv' in archive_files and 'genes.tsv' in archive_files:
            data_type = "mex"
        if 'DataMTX.tab' in archive_files and 'COLmeta.tab' in archive_files and 'ROWmeta.tab' in archive_files:
            data_type = "3tab"
        return data_type

    def read_mex_files(self, data_path = None):
        """
        Reads all files in the direcory and returns an Anndata object.

        Input
        -----
            data_path - path/to/extracted/archive

        Output
        ------
            Output is anndata object created using the files in directory
            'adata' is an AnnData object where data of each MEX file is assigned to the AnnData object:
                AnnData.X = matrix.mtx
                AnnData.obs = barcodes.tsv
                AnnData.var = genes.tsv

        """
        is_en=False
        archive_files_list = os.listdir(data_path)
        archive_files = [os.path.join(data_path, s) for s in archive_files_list]
        for entry in archive_files:
            if 'matrix.mtx' in entry:
                adata = sc.read(entry, cache=False).transpose()
            elif 'barcodes.tsv' in entry:
                obs = pd.read_csv(entry, sep='\t', header=None,index_col = 0, names=['observations'])
            elif 'genes.tsv' in entry:
                var = pd.read_csv(entry, sep='\t', header=None, index_col = 0, names=['genes', 'gene_symbol'])
                if var.index.str.contains('ENS').sum() >1:
                    is_en=True
            elif 'EXPmeta.json' in entry:
                continue
            else:
                raise Exception("Unexpected file name: '{0}'. Excepted 'matrix.mtx', 'genes.tsv', or 'barcodes.tsv'.".format(entry))
        adata.var = var
        adata.obs = obs
        self.adata = adata
        self.originalPath = data_path
        return self, is_en

    def read_3tab_files(self, data_path= None):
        """
        Reads files in direcoty and returns an Anndata object

        Input
        -----
            data_path - /path/to/extracted/archive

        Output
        ------
            Output is anndata object created using the files in directory
            'adata' is an AnnData object where data of each TAB file is assigned to the AnnData object:
                AnnData.X = expression file
                AnnData.obs = observations file
                AnnData.var = genes file
        """
        is_en = False
        archive_files_list = os.listdir(data_path)
        archive_files = [os.path.join(data_path, s) for s in archive_files_list]
        for entry in archive_files:
            if 'expression.tab' in entry or 'DataMTX.tab' in entry:
                # Get columns and rows of expression data in list form.
                exp = pd.read_csv(entry, sep='\t', index_col=0, header=0)
                exp_obs = list(exp.columns)
                exp_genes= list(exp.index)
                adata = sc.read(entry, cache=False).transpose()
            elif 'observations.tab' in entry or 'COLmeta.tab' in entry:
                obs = pd.read_csv(entry, sep='\t', index_col=0, header=0)
            elif 'genes.tab' in entry or 'ROWmeta.tab' in entry:
                var = pd.read_csv(entry, sep='\t', index_col=0, header=0)
                if var.index.str.contains('ENS').sum() >1:
                    is_en=True
            elif 'EXPmeta' in entry:
                continue
            ### Needs to be changed to account for other file types PCA etc
            else:
                raise Exception("Unexpected file name: '{0}'.".format(entry))

        for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
            if str_type in obs.columns:
                obs[str_type] = pd.Categorical(obs[str_type])
        for num_type in ['replicate', 'time_point_order']:
            if num_type in obs.columns:
                obs[num_type] = pd.to_numeric(obs[num_type])

        # Ensure observations and genes are sorted the same as found in the expressions file
        obs_index = list(obs.index)
        if set(obs_index) != set(exp_obs):
            raise Exception("Observation IDs from 'expressions' and 'observations' files are not the same.")
        obs = obs.reindex(exp_obs)

        genes_index = list(var.index)
        if set(genes_index) != set(exp_genes):
            raise Exception("Gene IDs from 'expressions' and 'genes' files are not the same.")
        var = var.reindex(exp_genes)

        adata.var = var
        adata.obs = obs
        self.adata = adata
        self.originalPath = data_path

        return self, is_en

    def write_h5ad(self, output_path, gear_identifier):
        """
        Writes anndata object as h5ad. Returns path to file.

        Input
        -----
            output_path - path/to/write/converted/file
            gear_identifier - unique name for output file

        Output
        ------
            path to successfully generated h5ad format file.

        """
        if self.adata is None:
            raise Exception("No AnnData object present to write to file.")
        if output_path is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            self.adata.write(filename = output_path)
        except Exception as err:
            raise Exception("Error occurred while writing to file: ", err)
        return ""



