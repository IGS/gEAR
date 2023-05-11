import os,sys
import shutil
import uuid
import anndata
import pandas as pd
import scanpy as sc
import gzip, tarfile

from pathlib import Path

class DataArchive:
    """
    This class handles:
    1. Get dataset type (mex or 3tab)
    2. Create h5ad and return
    3. Write h5ad to specfic directory
    """

    def __init__(self, base_dir=None, extracted_file_paths=None, format=None, tarball_path=None, dataset_id=None):
        # If the archive is a tarball, this is the source path for it
        self.tarball_path = tarball_path
        if self.tarball_path:
            self.base_dir = self.extract_dataset()

        # If the archive is extracted into a directory, this should be
        #  that directory base, directly within which all files are found.
        if base_dir:
            self.base_dir = base_dir
        if not self.base_dir:
            raise Exception("No base dataset directory detected from arguments or extracted tarball")

        # The full path to each file in the archive
        if extracted_file_paths:
            self.extracted_file_paths = list(extracted_file_paths.split(","))
            # TODO: Do one Path.iterdir here instead of in all these methods

        # Can be explicit or derived (mex, 3tab, etc.)
        self.format = format if format else self.get_archive_type(data_path=self.base_dir)
        if self.format == "tabcounts":
            self.format == "3tab"

        # Dataset ID to name H5AD file
        self.dataset_id = dataset_id if dataset_id else str(uuid.uuid4())

        # Are ensembl IDs present in our file?
        self.ens_present = False

    def extract_dataset(self, save_path=None):
        """
        Input: A path to an input dataset tarball, and the base directory where output can be
            written temporarily.

        Output: This function will extract the .tar or .tar.gz file and return the path to the
            directory created.

        Assumptions:  The specification states that the tar or tarball should create a unique
            directory name within which all the files of the dataset are contained.

        Example:
            Input:  /path/to/DLPFCcon322polyAgeneLIBD.3tab.tar.gz
            Output: /path/to/DLPFCcon322polyAgeneLIBD/DLPFCcon322polyAgeneLIBD_COLmeta.tab
                                                    ./DLPFCcon322polyAgeneLIBD_DataMTX.tab
                                                    ./DLPFCcon322polyAgeneLIBD_ROWmeta.tab
                                                    ./DLPFCcon322polyAgeneLIBD_EXPmeta.json
            Returns: /path/to/DLPFCcon322polyAgeneLIBD
        """

        if not (self.tarball_path and tarfile.is_tarfile(self.tarball_path)):
            return None

        save_path = Path(save_path) if save_path else Path(self.tarball_path.parent)

        tar = tarfile.open(self.tarball_path)
        sample_file = tar.next().name   # Get path of first member so we can extract directory later
        if not sample_file:
            raise Exception("Tar file {} appears to be empty".format(self.tarball_path))
        tar.extractall(path=save_path)
        tar.close()
        return save_path

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
        data_path = data_path if data_path else self.base_dir

        archive_filepaths = Path(data_path).iterdir()
        archive_filenames = [str(f) for f in archive_filepaths]
        archive_files=''.join(archive_filenames)
        data_type = None
        if 'expression.tab' in archive_files and 'genes.tab' in archive_files and 'observations.tab' in archive_files:
            data_type = "3tab"
        if 'matrix.mtx' in archive_files and 'barcodes.tsv' in archive_files and ('genes.tsv' in archive_files or 'features.tsv' in archive_files):
            data_type = "mex"
        if 'DataMTX.tab' in archive_files and 'COLmeta.tab' in archive_files and 'ROWmeta.tab' in archive_files:
            data_type = "3tab"
        if ".h5ad" in archive_files:
            return "h5ad"
        return data_type

    def convert_to_h5ad(self, output_dir):
        """
        Input: An extracted directory containing expression data files to be converted
            to H5AD.  These can be MEX or 3tab and should be handled appropriately.

        Output: An H5AD file should be created and the path returned.  The name of the
            file created should use the ID passed, like this:

            f50c5432-e9ca-4bdd-9a44-9e1d624c32f5.h5ad

        TBD: Error with writing to file.
        """
        # If dataset directory has h5ad file, skip that step
        if self.format=="h5ad":
            # assume ENSEMBL IDs are not present if h5ad was already passed to us
            self.ens_present = False
            file_gen = Path(self.base_dir).iterdir()
            return str(next(file_gen))

        outdir_name_path = Path(output_dir).joinpath(self.dataset_id + ".h5ad")
        outdir_name = str(outdir_name_path)
        # If errors occur in gEAR's parsing steps propagate upwards
        try:
            if self.format == "3tab":
                self.read_3tab_files()
            elif self.format == "mex":
                self.read_mex_files()
            else:
                raise Exception("Undetermined Format: {0}".format(self.format))
            self.write_h5ad(output_path=outdir_name)
        except:
            raise

        return outdir_name

    def read_mex_files(self, data_path=None):
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
        data_path = data_path if data_path else self.base_dir

        archive_filepaths = Path(data_path).iterdir()
        archive_files = [str(f) for f in archive_filepaths]

        for entry in archive_files:
            if 'matrix.mtx' in entry:
                unzipped_entry = gunzip_file(entry)
                adata = sc.read(unzipped_entry, cache=False).transpose()
            elif 'barcodes.tsv' in entry:
                unzipped_entry = gunzip_file(entry)
                obs = pd.read_csv(unzipped_entry, sep='\t', header=None,index_col = 0, names=['observations'])
            elif 'genes.tsv' in entry or 'features.tsv' in entry:
                unzipped_entry = gunzip_file(entry)
                var = pd.read_csv(unzipped_entry, sep='\t', header=None, index_col = 0, names=['genes', 'gene_symbol'])
                if var.index.str.contains('ENS').sum() >1:
                    self.ens_present = True
            elif 'EXPmeta.json' in entry:
                continue
            else:
                raise Exception("Unexpected file name: '{0}'. Excepted 'matrix.mtx', 'genes.tsv', or 'barcodes.tsv'.".format(entry))
        adata.var = var
        adata.obs = obs
        self.adata = adata
        self.originalPath = data_path

    def read_3tab_files(self, data_path=None):
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
        data_path = data_path if data_path else self.base_dir

        archive_filepaths = Path(data_path).iterdir()
        archive_files = [str(f) for f in archive_filepaths]

        for entry in archive_files:
            if 'expression.tab' in entry or 'DataMTX.tab' in entry:
                unzipped_entry = gunzip_file(entry)
                # Get columns and rows of expression data in list form.
                exp = pd.read_csv(unzipped_entry, sep='\t', index_col=0, header=0)
                exp_obs = list(exp.columns)
                exp_genes= list(exp.index)
                adata = sc.read(entry, first_column_names=True, cache=False).transpose()
            elif 'observations.tab' in entry or 'COLmeta.tab' in entry:
                unzipped_entry = gunzip_file(entry)
                obs = pd.read_csv(unzipped_entry, sep='\t', index_col=0, header=0)
            elif 'genes.tab' in entry or 'ROWmeta.tab' in entry:
                unzipped_entry = gunzip_file(entry)
                var = pd.read_csv(unzipped_entry, sep='\t', index_col=0, header=0)
                if var.index.str.contains('ENS').sum() >1:
                    self.ens_present = True
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

    def write_h5ad(self, output_path):
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

def gunzip_file(gzip_file:str):
    """Run "gunzip" on a file and return the extracted filename."""
    if not gzip_file.endswith(".gz"):
        return gzip_file
    gunzip_file = gzip_file.replace(".gz", "")
    with gzip.open(gzip_file, 'rb') as f_in:
        with open(gunzip_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Now that file is extracted.  Remove file
    Path(gzip_file).unlink()

    return gunzip_file

