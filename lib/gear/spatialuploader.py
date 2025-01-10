
import tarfile, os
from abc import ABC, abstractmethod

import pandas as pd
import anndata as ad
import spatialdata as sd

import spatialdata_io as sdio
from spatialdata_io.experimental import to_legacy_anndata

#class SpatialUploader(FileType):
class SpatialUploader(metaclass=ABC):

    NORMALIZED_TABLE_NAME = "table"

    @property
    @abstractmethod
    def normalized_table_name(self):
        return self.NORMALIZED_TABLE_NAME

    @property
    @abstractmethod
    def adata(self):
        return self._adata

    @adata.setter
    @abstractmethod
    def adata(self, adata=ad.AnnData()):
        self._adata = adata

    @property
    @abstractmethod
    def sdata(self):
        return self._sdata

    @sdata.setter
    @abstractmethod
    def sdata(self, sdata=sd.SpatialData()):
        self._sdata = sdata

    @abstractmethod
    def _read_file(self, filepath):
        pass

    @abstractmethod
    def _write_to_zarr(self, filepath=None):
        if self.sdata is None:
            raise Exception("No spatial data object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            self.sdata.write(file_path=filepath)
        except Exception as err:
            raise Exception("Error occurred while writing to file: ", err)
        return self

    @abstractmethod
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

class CurioUploader(SpatialUploader):
    # TODO: Test this class
    """
    Called by datasetuploader.py (factory) when a Curio Seeker dataset is going to be uploaded

    Standardized names for different files:
    * <dataset_id>_`'anndata.h5ad'`: Counts and metadata file.
    * <dataset_id>_`'cluster_assignment.txt'`: Cluster assignment file.
    * <dataset_id>_`'Metrics.csv'`: Metrics file.
    * <dataset_id>_`'variable_features_clusters.txt'`: Variable features clusters file.
    * <dataset_id>_`'variable_features_spatial_moransi.txt'`: Variable features Moran’s I file.
    """

    def _read_file(self, filepath):
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

            sdata = sdio.curio(tmp_dir)

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata.var["gene_symbol"] = sdata.var.index

            # set the index to the ensembl id (gene_ids)
            sdata.var.set_index("gene_ids", inplace=True)

            self.sdata = sdata

            # table name should already be "table" for Visium

            adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires", table_name=self.NORMALIZED_TABLE_NAME)

            # Apply AnnData obj and filepath to uploader obj
            self.adata = adata
            self.originalFile = filepath
            return self

    def _write_to_zarr(self, filepath=None):
        return super()._write_to_zarr(filepath)

    def _write_to_h5ad(self, filepath=None):
        return super()._write_to_h5ad(filepath)

class VisiumUploader(SpatialUploader):
    # TODO: Test this class
    """
    Called by datasetuploader.py (factory) when a Visium or Visium-HD dataset is going to be uploaded

    Standardized names for different files:
    * (<dataset_id>_)`'filtered_feature_bc_matrix.h5'`: Counts and metadata file.
    * 'spatial/tissue_hires_image.png': High resolution image.
    * 'spatial/tissue_lowres_image.png': Low resolution image.
    * 'scalefactors_json.json': Scalefactors file.
    * 'tissue_positions_list.csv' (SpaceRanger 1) or 'tissue_positions.csv' (SpaceRanger 2): Spots positions file.
    * fullres_image_file: large microscopy image used as input for space ranger.
    """

    def _read_file(self, filepath):

        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

            sdata = sdio.visium(tmp_dir)#, dataset_id="spatialdata")    # Provide a name to standarize downstream usage

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata.var["gene_symbol"] = sdata.var.index

            # set the index to the ensembl id (gene_ids)
            sdata.var.set_index("gene_ids", inplace=True)

            self.sdata = sdata

            # table name should already be "table" for Visium

            adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires", table_name=self.NORMALIZED_TABLE_NAME)

            # Apply AnnData obj and filepath to uploader obj
            self.adata = adata
            self.originalFile = filepath
            return self

    def _write_to_zarr(self, filepath=None):
        return super()._write_to_zarr(filepath)

    def _write_to_h5ad(self, filepath=None):
        return super()._write_to_h5ad(filepath)

class VisiumHDUploader(SpatialUploader):
    """
    Called by datasetuploader.py (factory) when a Visium-HD dataset is going to be uploaded

    Explanation of Space Ranger v3 output here -> https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview#hd-outputs

    Required files that will break the upload if not present:
    /binned_outputs/square_008um/filtered_feature_bc_matrix.h5
    /binned_outputs/feature_slice.h5

    Special note: We have observed that bin sizes finer than 8 microns per pixel will generally have more cells, which can lead to longer analysis times.
    For now, we will attempt to use the "square_008um" binned output.

    """

    TABLE_NAME = "square_008um"


    def _read_file(self, filepath, ):
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        binned_outputs_dir = "{}/binned_outputs".format(tmp_dir)
        bin_008_dataset_path = "{}/{}/".format(binned_outputs_dir, self.TABLE_NAME)
        clustering_csv_path = "{}/analysis/clustering/gene_expression_graphclust/clusters.csv".format(bin_008_dataset_path)

        absolute_path = os.path.abspath(binned_outputs_dir)

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(binned_outputs_dir, entry.name)
                tf.extract(entry, path=binned_outputs_dir)

        if not os.path.exists("{}/feature_slice.h5".format(binned_outputs_dir)):
            raise Exception("feature_slice.h5 file not found in /binned_outputs directory in tarball.")

        # If clustering file does not exist, raise an exception
        if not os.path.exists(clustering_csv_path):
            raise Exception("clusters.csv file not found in tarball.")

        # If clustering file does not have "Barcode" and "Cluster" columns, raise an exception
        with open(clustering_csv_path, 'r') as f:
            first_line = f.readline()
            if "Barcode" not in first_line or "Cluster" not in first_line:
                raise Exception("clusters.csv file does not have 'Barcode' and 'Cluster' columns in clusters.csv file in tarball.")

        # https://github.com/scverse/spatialdata-io/issues/212

        # Create a symlink for feature_slice.h5 within "visium_dataset_path" to include a "spatialdata_" prefix
        # This is a workaround for the current implementation of the visium_hd function

        if not os.path.exists("{}/spatialdata_feature_slice.h5".format(binned_outputs_dir)):
            os.symlink("{}/feature_slice.h5".format(absolute_path), "{}/spatialdata_feature_slice.h5".format(binned_outputs_dir))

            sdata = sdio.visium_hd(tmp_dir
                                   , dataset_id="spatialdata"   # Provide a name to standarize downstream usage
                                   , bin_size=8
                                   , filtered_counts_file=True
                                   , load_all_images=False
                                   , fullres_image_file=None
                                   , bins_as_squares=True
                                   )

            # add clustering information to the vis_sdata.table.obs dataframe
            clustering = pd.read_csv(clustering_csv_path)
            # make barcode as index
            clustering.set_index('Barcode', inplace=True)
            sdata[self.TABLE_NAME].obs['clusters'] = clustering['Cluster'].astype('category')

            # To get the adata equivalent, look at sdata.tables["table"]
            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.TABLE_NAME].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata[self.TABLE_NAME].var["gene_symbol"] = sdata[self.TABLE_NAME].var.index

            # set the index to the ensembl id (gene_ids)
            sdata[self.TABLE_NAME].var.set_index("gene_ids", inplace=True)

            # Set the table name to the normalized table name
            sdata.tables[self.normalized_table_name] = sdata[self.TABLE_NAME]

            self.sdata = sdata

            adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires", table_name=self.normalized_table_name)

            # Apply AnnData obj and filepath to uploader obj
            self.adata = adata
            self.originalFile = filepath
            return self

    def _write_to_zarr(self, filepath=None):
        return super()._write_to_zarr(filepath)

    def _write_to_h5ad(self, filepath=None):
        return super()._write_to_h5ad(filepath)

class XeniumUploader(SpatialUploader):
    # TODO: Test this class
    """
    Called by datasetuploader.py (factory) when a Xenium dataset is going to be uploaded

    Standardized names for different files:
    * 'experiment.xenium': File containing specifications.
    * 'nucleus_boundaries.parquet': Polygons of nucleus boundaries.
    * 'cell_boundaries.parquet': Polygons of cell boundaries.
    * 'transcripts.parquet': File containing transcripts.
    * 'cell_feature_matrix.h5': File containing cell feature matrix.
    * 'cells.parquet': File containing cell metadata.
    * 'morphology_mip.ome.tif': File containing morphology mip.
    * 'morphology_focus.ome.tif': File containing morphology focus.

    More information on Xenium Ranger outputs can be found here: https://www.10xgenomics.com/support/software/xenium-ranger/latest/analysis/outputs/XR-output-overview
    """

    def _read_file(self, filepath):
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

            sdata = sdio.xenium(tmp_dir)

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata.var["gene_symbol"] = sdata.var.index

            # set the index to the ensembl id (gene_ids)
            sdata.var.set_index("gene_ids", inplace=True)

            self.sdata = sdata

            # table name should already be "table" for Visium

            adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires", table_name=self.NORMALIZED_TABLE_NAME)

            # Apply AnnData obj and filepath to uploader obj
            self.adata = adata
            self.originalFile = filepath
            return self

    def _write_to_zarr(self, filepath=None):
        return super()._write_to_zarr(filepath)

    def _write_to_h5ad(self, filepath=None):
        return super()._write_to_h5ad(filepath)