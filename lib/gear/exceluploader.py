import anndata
import pandas as pd
import os, sys
import scanpy as sc
sc.settings.verbosity = 0

from gear.datasetuploader import DatasetUploader


class ExcelUploader(DatasetUploader):
    """
    Called by datasetuploader.py (factory) when excel file is going uploaded

    Standardized names for different sheets:
        'expression' - contains the expression values matrix.
            It is assumed the 1st column is the row index and the 1st row is the column index
        'observations' - contains details of the column index from 'expression'.
        'genes' - contains details of the row index from 'expression' sheet. This sheet is optional.

    Standardized names for columns in each sheet:
        'expression'
            'genes' - Must be the first column on the left. This can contain either:
                1) Ensembl ids, or;
                2) Gene symbols
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
        Currently reads excel file 3x (once for each sheet) for validation purposes. In doing so, this allows for closer
        validation and error handling.

        Input
        -----
            filepath - /path/to/your/excel.xlsx

        Output
        ------
            Output is assigned to the DatasetUploader object:
                DatasetUploader.adata = adata
                DatasetUploader.originalFile = filepath

            'adata' is an AnnData object where data of each XLSX sheet is assigned to
            the AnnData object:
                AnnData.X = expression sheet
                AnnData.obs = observations sheet
                AnnData.var = genes sheet
                    (if absent, this only contains a pandas dataframe index of gene ids)

            'filepath' is the file path of the original file
        """

        # Get the expression matrix
        try:
            exp_df = pd.read_excel(filepath, sheet_name='expression', index_col=0).transpose()
        except:
            raise Exception("No expression sheet found. Expected spreadsheet sheet named 'expression'.")

        try:
            X = exp_df.values[:, 0:].astype(float)
        except:
            raise Exception("Encountered unexpected value type. Expected float type in expression matrix.")

        # Get counts of genes and observations
        number_obs_from_exp, number_genes_from_exp = X.shape

        # Get the observations
        try:
            obs_df = pd.read_excel(filepath, sheet_name='observations', index_col='observations')
        except:
            raise Exception("No observations sheet found. Expected spreadsheet sheet named 'observations'.")

        # Verify number observations equal those found in expression sheet
        number_obs, number_cond = obs_df.shape
        if number_obs != number_obs_from_exp:
            raise Exception("Observations sheet error. Row count ({0}) in 'observations' sheet must match column count of 'expression' sheet({1}).".format(number_obs, number_obs_from_exp))

        # Verify observations index matches expression sheet col index
        if not obs_df.index.equals(exp_df.index):
            raise Exception("Observations sheet error. The names and order of the index column in 'observations' sheet must match the columns of 'expression' sheet.")

        # Get the genes (if present), else use the .var from exp_df
        try:
            genes_df = pd.read_excel(filepath, sheet_name='genes', index_col=0, converters={'gene_symbol': str})
        except:
            # With 'genes' sheet absent try to get the genes from 'expression'
            try:
                # read expression sheet. Set genes as the index and only parse 1 column
                genes_df = pd.read_excel(filepath, sheet_name='expression', index_col=0, usecols=[0,1])

                # remove the 1 parsed column. This leaves an empty dataframe
                # with an index of gene ids to match other datasets.
                genes_df = genes_df.drop(genes_df.columns[0], axis=1)
            except Exception as err:
                error = "No 'genes' sheet found. Tried using genes column from 'expression' sheet as .var, but " + str(err)
                raise Exception(error)

        # Check for numeric gene symbols
        if 'gene_symbol' in genes_df.columns:
            digit_count = genes_df['gene_symbol'].str.isnumeric().sum()
            if digit_count > 0:
                raise Exception("Genes sheet error. {0} gene symbols are listed as numbers, not gene symbols".format(str(digit_count)))
        else:
            raise Exception("Failed to find gene_symbol column in genes tab")

        for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
            if str_type in obs_df.columns:
                obs_df[str_type] = pd.Categorical(obs_df[str_type])

        for num_type in ['replicate', 'time_point_order']:
            if num_type in obs_df.columns:
                obs_df[num_type] = pd.to_numeric(obs_df[num_type])

        # Verify number genes equal those found in expression sheet
        number_genes, number_columns = genes_df.shape
        if number_genes != number_genes_from_exp:
            raise Exception("Genes sheet error. Row count ({0}) in 'genes' sheet must match row count of 'expression' sheet({1}).".format(number_genes, number_genes_from_exp))

        # Verify genes index matches expression sheet columns
        if not genes_df.index.equals(exp_df.columns):
            raise Exception("Genes sheet error. The names and order of the index column in 'genes' sheet must match the rows of 'expression' sheet.")

        # Create AnnData object and return it
        self.adata = anndata.AnnData(X=X, obs=obs_df, var=genes_df)
        self.originalFile = filepath
        return self

    def get_plot_preview_expression(self):
        """
        Searches the adata.X values for the first set of values that are all present and
        non-zero.

        Returns tuple containing (gene_symbol, expression):
            1) Gene symbol (if present) or Ensembl id of gene.
            2) Pandas dataframe containing the expression values of the first
               gene containing values for all observations and are non-zero.
               The dataframe contains all columns of observation characteristics
               from adata.obs.
        """
        expression = None
        index = None
        gene = None
        X = pd.DataFrame(self.adata.X)
        var = self.adata.var
        obs = self.adata.obs

        index_length = len(X.index)
        for i in range(0, len(X.columns)):
            # Find a gene whose values are all populated and non-zero
            if pd.notnull(X).all(1).sum() == index_length:
                index = i

                # Get the gene's expression and add the observation characteristics
                expression = pd.DataFrame(X[i])
                expression.columns = ['raw_value']

                for col in obs.columns:
                    expression[col] = obs[col].values
                break

        #Get gene symbol (if present), otherwise use index label (Ensembl id)
        if 'gene_symbol' in var.columns:
            gene = var.iloc[index]['gene_symbol']
        else:
            gene = var.index.values[index]

        # Return gene and it's expression dataframe
        return gene, expression


    def _add_calculated_values(self, adata):
        """
        TODO: calculations in DatasetStats make inaccurate assumption when determining
              the number of replicates within a condition (observation).

        **Correcting this is now hold until handling replicates observational and technical is addressed and finalized. **

        Runs the statistical calculations on AnnData.X and stores them to AnnData.uns

        Input
        -----
        AnnData object

        Output
        ------
        AnnData object with statistics. statistics are stored as numpy arrays in
        the structured AnnData dictionary:
            AnnData.uns['Xmean'] contains a numpy array martix of means

        """
        from gear.datasetstats import DatasetStats

        # Calculate replicate averages. Return them as numpy array called .Xmean
        Xmean = DatasetStats.get_replicate_averages(adata)
        adata.uns['Xmean'] = Xmean

        # Calculate replicate standard deviations
        Xstd = DatasetStats.get_replicate_std(adata)
        adata.uns['Xstd'] = Xstd

        # Calculate replicate standard errors (of the mean)
        Xsem = DatasetStats.get_replicate_sem(adata)
        adata.uns['Xsem'] = Xsem

        # Calculate replicate p-value (uses Pearson from scipy.stats)
        Xpval = DatasetStats.get_replicate_pvalue(adata)
        adata.uns['Xpval'] = Xpval

        # TODO Needs refactored now that Xpval is stored in .uns['Xpval']
        # Calculate replicate FDR
        # Xfdr = DatasetStats.get_replicate_fdr(adata)
        # adata.uns['Xfdr'] = Xfdr

        self.adata = adata
        return self


    def _add_color_values(self, adata):
        #calculate raw and absolute coloring
        """
        TODO: calculations in DatasetColoring make inaccurate assumption when determining
              the number of replicates within a condition (observation).
            - Dataset scoped coloring needs completed

        **Correcting this is now hold until handling replicates observational and technical is addressed and finalized. **

        """
        from gear.datasetcoloring import DatasetColoring

        # Calcuate gene level coloring - raw
        XColGenRaw = DatasetColoring.get_gene_coloring(adata, color_mode='raw')
        adata.uns['XColGenRaw'] = XColGenRaw

        # Calcuate gene level coloring - absolute
        XColGenAbs = DatasetColoring.get_gene_coloring(adata, color_mode='absolute')
        adata.uns['XColGenAbs'] = XColGenAbs

        # Calcuate tissue level coloring - raw
        XColTisRaw = DatasetColoring.get_tissue_coloring(adata, color_mode='raw')
        adata.uns['XColTisRaw'] = XColTisRaw

        # Calcuate tissue level coloring - absolute
        XColTisAbs = DatasetColoring.get_tissue_coloring(adata, color_mode='absolute')
        adata.uns['XColTisAbs'] = XColTisAbs

        # TODO
        # get_color_dataset()
        # get_color_abs_dataset()


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
