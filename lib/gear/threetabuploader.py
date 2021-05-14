import anndata
import pandas as pd
import os, sys
import scanpy as sc
import tarfile
sc.settings.verbosity = 0

from gear.datasetuploader import FileType


class ThreeTabUploader(FileType):
    """
    Called by datasetuploader.py (factory) when excel file is going uploaded

    Standardized names for different files:
        'expression.tab' - contains the expression values matrix.
            It is assumed the 1st column is the row index and the 1st row is the column index
        'observations.tab' - contains details of the column index from 'expression'.
        'genes.tab' - contains details of the row index from 'expression' sheet.  Should at least have gene symbol

    Standardized names for columns in each file:
        'expression.tab'
            'genes' - Must be the first column on the left. This should be the Ensembl ID
            'gene_symbol' - Such as 'Rfx7' or 'TonB'
        'observations' -
            'observations' - This column must contain the identical names of the first row from the 'expression' sheet.
            'cell_type' - This column is the anatomical region or cell_type. Examples: 'utricle', 'inner_hair_cell'
            'condition' - This column is the experimental condition of the observation. Examples: 'control', 'treated'
            'replicate' - This column is the number of replicate the observation is. If this is replicate 1 of 3, put '1'.
            'time_point' - This column is the time point, or age, at which the observation was taken. Examples: 'P0', '24'
            'time_point_order' - This column is a numeric value allowing for sorting of the time_point column values.
            'time_unit' - What is the unit of time used? hours? minutes?
    """

    def _read_file(self, filepath):
        """
        Input
        -----
            filepath - /path/to/your/some.tar[.gz]

        Output
        ------
            Output is assigned to the DatasetUploader object:
                DatasetUploader.adata = adata
                DatasetUploader.originalFile = filepath

            'adata' is an AnnData object where data of each TAB file is assigned to
            the AnnData object:
                AnnData.X = expression file
                AnnData.obs = observations file
                AnnData.var = genes file

            'filepath' is the file path of the original tarball
        """

        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

                # Read each file as pandas dataframes
                # Added functionality to read files inside folder within tarball and for NEMO three tab files
                if entry.name == 'expression.tab' or os.path.basename(filepath)== 'expression.tab' or 'DataMTX.tab' in entry.name:
                    # Get columns and rows of expression data in list form.
                    exp = pd.read_table(filepath, sep='\t', index_col=0, header=0)
                    exp_obs = list(exp.columns)
                    exp_genes= list(exp.index)

                    # Read in expressions as AnnData object
                    adata = sc.read(filepath, first_column_names=True, cache=False).transpose()
                elif entry.name == 'observations.tab' or os.path.basename(filepath)== 'observations.tab' or 'COLmeta.tab' in entry.name:
                    obs = pd.read_table(filepath, sep='\t', index_col=0, header=0)
                elif entry.name == 'genes.tab' or os.path.basename(filepath)== 'genes.tab' or 'ROWmeta.tab' in entry.name:
                    var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

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

        # Assign genes and observations to AnnData object
        adata.var = var
        adata.obs = obs

        # Apply AnnData obj and filepath to uploader obj
        self.adata = adata
        self.originalFile = filepath
        return self


    def _write_to_h5ad(self, filepath=None):
        ###Added this subroutine from mexuploader as it was not writing to h5ad without it.
        if self.adata is None:
            raise Exception("No AnnData object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            self.adata.write(filename=filepath)
        except Exception as err:
            raise Exception("Error occurred while writing to file: ", err)
        return self
