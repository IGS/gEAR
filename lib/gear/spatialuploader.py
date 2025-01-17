
import tarfile, os
from abc import ABC, abstractmethod

import pandas as pd
import anndata as ad
import spatialdata as sd

import spatialdata_io as sdio
from spatialdata_io.experimental import to_legacy_anndata

class SpatialUploader(ABC):

    NORMALIZED_TABLE_NAME = "table"

    @property
    def normalized_table_name(self):
        return self.NORMALIZED_TABLE_NAME

    @property
    def adata(self):
        return self._adata

    @adata.setter
    def adata(self, adata=ad.AnnData()):
        self._adata = adata

    @property
    def sdata(self):
        return self._sdata

    @sdata.setter
    def sdata(self, sdata=sd.SpatialData()):
        self._sdata = sdata

    @property
    @abstractmethod
    def has_images(self):
        pass

    @property
    @abstractmethod
    def coordinate_system(self):
        pass

    @property
    @abstractmethod
    def platform(self):
        pass

    @property
    @abstractmethod
    def img_name(self):
        pass

    @abstractmethod
    def _read_file(self, filepath):
        pass

    def _filter_sdata_by_coords(self):
        # Filter to only the hires image boundaries
        if not self.img_name:
            return self

        x = len(self.sdata.images[self.img_name].x)
        y = len(self.sdata.images[self.img_name].y)
        sdata = sd.bounding_box_query(self.sdata,
                axes=("x", "y"),
                min_coordinate=[0, 0],
                max_coordinate=[x, y],
                target_coordinate_system=self.coordinate_system,
                filter_table=True,
                )
        self.sdata = sdata
        return self

    def _convert_sdata_to_adata(self, include_images=None, table_name=None):
        if self.sdata is None:
            raise Exception("No spatial data object present to convert to AnnData object.")

        if include_images is None:
            include_images = self.has_images

        # Generally everything should already be converted to the normalized table name
        if table_name is None:
            table_name = self.NORMALIZED_TABLE_NAME

        try:
            # table name should already be "table" for Visium
            adata = to_legacy_anndata(self.sdata, include_images=include_images, coordinate_system=self.coordinate_system, table_name=table_name)
        except Exception as err:
            raise Exception("Error occurred while converting spatial data object to AnnData object: ", err)
        self.adata = adata
        return self

    def _write_to_zarr(self, filepath=None):
        if self.sdata is None:
            raise Exception("No spatial data object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            # Add platform type as metadata. Can use downstream
            self.sdata.tables[self.NORMALIZED_TABLE_NAME].uns["platform"] = self.platform

            # Will fail if file already exists
            self.sdata.write(file_path=filepath)
        except Exception as err:
            raise Exception("Error occurred while writing to file: ", err)
        return self

    def _write_to_h5ad(self, filepath=None):
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
    """
    Called by datasetuploader.py (factory) when a Curio Seeker dataset is going to be uploaded

    Standardized names for different files:
    * <dataset_id>_`'anndata.h5ad'`: Counts and metadata file.
    * <dataset_id>_`'cluster_assignment.txt'`: Cluster assignment file.
    * <dataset_id>_`'Metrics.csv'`: Metrics file.
    * <dataset_id>_`'variable_features_clusters.txt'`: Variable features clusters file.
    * <dataset_id>_`'variable_features_spatial_moransi.txt'`: Variable features Moran’s I file.

    It seems Curio Seeker only provides gene IDs (as indexes to an empty dataframe),
    so you need to modify the included h5ad file to include ensembl IDs in adata.var.
    You can use add_ensembl_id_to_h5ad_missing_release.py for that, and include the revisted h5ad file in the tarball.
    """

    @property
    def has_images(self):
        return False

    @property
    def coordinate_system(self):
        return "global"

    @property
    def platform(self):
        return "curio"

    @property
    def img_name(self):
        return None

    def _read_file(self, filepath):
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        h5ad_file = None
        spatial_moransi_file = None

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if entry.name.startswith("._") or entry.name.startswith(".DS_Store"):
                    continue
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

                if entry.name.endswith("anndata.h5ad"):
                    h5ad_file = filepath
                elif entry.name.endswith("variable_features_spatial_moransi.txt"):
                    spatial_moransi_file = filepath

            """
            The sdio.curio function below uses the Moran's I-score file to overwrite the adata.var that was previously set
            with gene IDs. We need to pull out the original h5ad file and set the adata.var back to the original Dataframe
            """

            # Read in the h5ad file
            adata = ad.read_h5ad(h5ad_file)
            # Create mapping dict.
            # Original index name (ensembl_id) was created in add_ensembl_id_to_h5ad_missing_release.py
            gene_symbol_to_ensembl_id = adata.var.reset_index().set_index("gene_symbol").to_dict()["ensembl_id"]
            var_features_moransi = pd.read_csv(spatial_moransi_file, sep="\t", header=0)

            # Replace the gene symbols with the ensembl ids
            var_features_moransi.index = var_features_moransi.index.map(gene_symbol_to_ensembl_id)
            var_features_moransi.to_csv(spatial_moransi_file, sep="\t", header=True, index=True, index_label=False)

            # Now are ready to read in to a SpatialData object
            sdata = sdio.curio(tmp_dir)

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

            # Change obs "cluster" to "clusters" to harmonize
            sdata[self.NORMALIZED_TABLE_NAME].obs.rename(columns={"cluster": "clusters"}, inplace=True)

            self.sdata = sdata

            # table name should already be "table" for Visium

            self.originalFile = filepath
            return self

    def _convert_sdata_to_adata(self, include_images=None):
        return super()._convert_sdata_to_adata(include_images)

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

    @property
    def has_images(self):
        return True

    @property
    def coordinate_system(self):
        return "downscaled_hires"

    @property
    def platform(self):
        return "visium"

    @property
    def img_name(self):
        return "spatialdata_hires_image"

    def _read_file(self, filepath):

        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if entry.name.startswith("._") or entry.name.startswith(".DS_Store"):
                    continue
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

            sdata = sdio.visium(tmp_dir)#, dataset_id="spatialdata")    # Provide a name to standarize downstream usage

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata[self.NORMALIZED_TABLE_NAME].var["gene_symbol"] = sdata[self.NORMALIZED_TABLE_NAME].var.index

            # set the index to the ensembl id (gene_ids)
            sdata[self.NORMALIZED_TABLE_NAME].var.set_index("gene_ids", inplace=True)

            self.sdata = sdata
            self.originalFile = filepath
            return self

    def _convert_sdata_to_adata(self, include_images=None):
        return super()._convert_sdata_to_adata(include_images)

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
    /binned_outputs/square_008um/spatial/scalefactors_json.json
    /binned_outputs/square_008um/spatial/tissue_hires_image.png

    Recommended tar command to create tarball:
    `tar cvf <dataset>.tar binned_outputs/feature_slice.h5 binned_outputs/square_008um spatial`

    Special note: We have observed that bin sizes finer than 8 microns per pixel will generally have more cells, which can lead to longer analysis times.
    For now, we will attempt to use the "square_008um" binned output.

    """

    table_name = "square_008um"

    @property
    def has_images(self):
        return True

    @property
    def coordinate_system(self):
        return "downscaled_hires"

    @property
    def platform(self):
        return "visium_hd"

    @property
    def img_name(self):
        return "spatialdata_hires_image"

    def _read_file(self, filepath):
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        binned_outputs_dir = "{}/binned_outputs".format(tmp_dir)
        bin_008_dataset_path = "{}/{}/".format(binned_outputs_dir, self.table_name)
        clustering_csv_path = "{}/analysis/clustering/gene_expression_graphclust/clusters.csv".format(bin_008_dataset_path)

        absolute_path = os.path.abspath(binned_outputs_dir)

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if entry.name.startswith("._") or entry.name.startswith(".DS_Store"):
                    continue

                # IF directory has "square_" but not "square_008um", skip
                if "square_" in entry.name and "square_008um" not in entry.name:
                    continue

                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

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

            sdata = sdio.visium_hd(binned_outputs_dir
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
            sdata[self.table_name].obs['clusters'] = clustering['Cluster'].astype('category')

            # To get the adata equivalent, look at sdata.tables["table"]
            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.table_name].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata[self.table_name].var["gene_symbol"] = sdata[self.table_name].var.index

            # set the index to the ensembl id (gene_ids)
            sdata[self.table_name].var.set_index("gene_ids", inplace=True)

            # Set the table name to the normalized table name
            sdata.tables[self.normalized_table_name] = sdata[self.table_name]

            self.sdata = sdata
            self.originalFile = filepath
            return self

    def _convert_sdata_to_adata(self, include_images=None):
        return super()._convert_sdata_to_adata(include_images)

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

    @property
    def has_images(self):
        return True

    @property
    def coordinate_system(self):
        return "downscaled_hires"

    @property
    def platform(self):
        return "xenium"

    @property
    def img_name(self):
        return "spatialdata_hires_image"

    def _read_file(self, filepath):
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if entry.name.startswith("._") or entry.name.startswith(".DS_Store"):
                    continue
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

            sdata = sdio.xenium(tmp_dir)

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            sdata[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            sdata[self.NORMALIZED_TABLE_NAME].var["gene_symbol"] = sdata[self.NORMALIZED_TABLE_NAME].var.index

            # set the index to the ensembl id (gene_ids)
            sdata[self.NORMALIZED_TABLE_NAME].var.set_index("gene_ids", inplace=True)

            self.sdata = sdata
            self.originalFile = filepath
            return self

    def _convert_sdata_to_adata(self, include_images=None):
        return super()._convert_sdata_to_adata(include_images)

    def _write_to_zarr(self, filepath=None):
        return super()._write_to_zarr(filepath)

    def _write_to_h5ad(self, filepath=None):
        return super()._write_to_h5ad(filepath)

### Helper constants

SPATIALTYPE2CLASS = {
    "visium": VisiumUploader,
    "visium_hd": VisiumHDUploader,
    "xenium": XeniumUploader,
    "curio": CurioUploader
}