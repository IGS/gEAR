
import os
import tarfile
import typing
from abc import ABC, abstractmethod
from pathlib import Path

import anndata as ad
import pandas as pd
import spatialdata as sd
import spatialdata_io as sdio
import xarray
from gear.utils import update_adata_with_ensembl_ids
from spatialdata_io.experimental import from_legacy_anndata, to_legacy_anndata

if typing.TYPE_CHECKING:
    from anndata import AnnData
    from spatialdata import SpatialData

class SpatialHandler(ABC):
    """
    Abstract base class for handling spatial transcriptomics data, providing a unified interface for managing and processing spatial and annotated data objects.

    This class defines properties and methods for:
    - Managing normalized table names.
    - Accessing and setting AnnData and SpatialData objects.
    - Querying metadata such as image presence, coordinate system, platform, and image names.
    - Reading spatial data files.
    - Filtering spatial data to image boundaries.
    - Converting spatial data to AnnData format.
    - Writing spatial data to Zarr and AnnData (H5AD) file formats.

    Subclasses must implement abstract properties and methods to specify platform-specific details and file reading logic.

    Attributes:
        NORMALIZED_TABLE_NAME (str): Default name for the normalized table.

    Properties:
        normalized_table_name (str): Name of the normalized table.
        adata (AnnData): The AnnData object associated with this handler.
        sdata (SpatialData): The SpatialData object associated with this handler.
        has_images (bool): Whether the handler has associated images (abstract).
        coordinate_system (str): The coordinate system used (abstract).
        platform (str): The name of the platform (abstract).
        img_name (str | None): The name of the associated image (abstract).

    Methods:
        _read_file(filepath: str) -> SpatialHandler:
            Reads and processes a spatial data file from the given filepath (abstract).

        filter_sdata_by_coords() -> SpatialHandler:
            Filters the spatial data to include only elements within the boundaries of the high-resolution image.

        _convert_sdata_to_adata(include_images: bool | None = None, table_name=None) -> SpatialHandler:
            Converts the internal spatial data object to an AnnData object.

        _write_to_zarr(filepath: str | None = None) -> SpatialHandler:

        _write_to_h5ad(filepath: str | None = None) -> SpatialHandler:
            Writes the current AnnData object to an H5AD file at the specified file path.
    """

    NORMALIZED_TABLE_NAME = "table"

    @property
    def normalized_table_name(self):
        """
        Returns the name of the normalized table associated with this instance.

        Returns:
            str: The normalized table name.
        """
        return self.NORMALIZED_TABLE_NAME

    @property
    def adata(self) -> "AnnData":
        """
        Returns the AnnData object associated with this instance.

        Returns:
            AnnData: The underlying AnnData object.
        """
        return self._adata

    @adata.setter
    def adata(self, adata: "AnnData") -> None:
        """
        Sets the AnnData object for the instance.

        Parameters:
            adata (AnnData): The annotated data matrix to be assigned to the instance.
        """
        self._adata = adata

    @property
    def sdata(self) -> "SpatialData":
        """
        Returns the associated SpatialData object.

        Returns:
            SpatialData: The spatial data instance associated with this handler.
        """
        return self._sdata

    @sdata.setter
    def sdata(self, sdata: "SpatialData") -> None:
        """
        Sets the SpatialData object for the handler.

        Parameters:
            sdata (SpatialData): The SpatialData instance to assign to the handler.
        """
        self._sdata = sdata

    @property
    @abstractmethod
    def has_images(self) -> bool:
        """
        Checks if the spatial handler has associated images.

        Returns:
            bool: True if images are present, False otherwise.
        """
        pass

    @property
    @abstractmethod
    def coordinate_system(self) -> str:
        """
        Returns the coordinate system used by the spatial handler.

        Returns:
            str: A string representing the name or type of the coordinate system (e.g., 'WGS84', 'UTM', etc.).
        """
        pass

    @property
    @abstractmethod
    def platform(self) -> str:
        """
        Returns the name of the platform as a string.

        Returns:
            str: The name of the platform.
        """
        pass

    @property
    @abstractmethod
    def img_name(self) -> str | None:
        """
        Returns the name of the image associated with the current instance.

        Returns:
            str | None: The image name if available, otherwise None.
        """
        pass

    @abstractmethod
    def _read_file(self, filepath: str) -> "SpatialHandler":
        """
        Reads and processes a spatial data file from the given filepath.

        Args:
            filepath (str): The path to the spatial data file to be read.

        Returns:
            SpatialHandler: An instance of SpatialHandler containing the loaded spatial data.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            IOError: If there is an error reading the file.
            ValueError: If the file format is invalid or unsupported.
        """
        pass

    def filter_sdata_by_coords(self) -> "SpatialHandler":
        """
        Filters the spatial data (`self.sdata`) to include only elements within the boundaries of the high-resolution image.

        The method determines the spatial extent of the image specified by `self.img_name` and restricts the spatial data to this region.
        It supports images stored as either `xarray.DataArray` or `xarray.DataTree`. If the image is not found, the method returns the object unchanged.

        Returns:
            SpatialHandler: The current instance with `self.sdata` filtered to the image boundaries.
        """
        # Filter to only the hires image boundaries
        if not self.img_name:
            return self

        # Image can be DataArray or DataTree depending on if multiple scales are present for the image
        img = self.sdata[self.img_name]
        if isinstance(img, xarray.DataTree):
            coords = sd.get_pyramid_levels(img, n=0)
        elif isinstance(img, xarray.DataArray):
            coords = img
        else:
            coords = self.sdata.images[self.img_name]  # Fallback to preserve original behavior

        # Get the coordinates of the image
        x = len(coords.x) # type: ignore
        y = len(coords.y) # type: ignore
        sdata: "SpatialData" = sd.bounding_box_query(self.sdata,
                axes=("x", "y"),
                min_coordinate=[0, 0],
                max_coordinate=[x, y],
                target_coordinate_system=self.coordinate_system,
                filter_table=True,
                ) # type: ignore

        self.sdata = sdata
        return self

    def _convert_sdata_to_adata(self, include_images: bool | None = None, table_name=None) -> "SpatialHandler":
        """
        Converts the internal spatial data object (`sdata`) to an AnnData object and assigns it to `self.adata`.

        Parameters:
            include_images (bool | None, optional): Whether to include image data in the conversion. If None, defaults to `self.has_images`.
            table_name (str, optional): The name of the table to use for conversion. If None, defaults to `self.NORMALIZED_TABLE_NAME`.

        Returns:
            SpatialHandler: The current instance with the `adata` attribute set.

        Raises:
            Exception: If `self.sdata` is None or if an error occurs during conversion.
        """
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
            print(str(err))
            raise Exception("Error occurred while converting spatial data object to AnnData object: " + str(err))
        self.adata = adata
        return self

    def _write_to_zarr(self, filepath: str | None =None) -> "SpatialHandler":
        """
        Writes the spatial data object to a Zarr file at the specified file path.

        This method performs the following steps:
        - Checks that a spatial data object (`self.sdata`) is present.
        - Ensures a destination file path is provided.
        - Adds the platform type as metadata to the normalized table.
        - Converts the data matrix to float type for compatibility.
        - Attempts to write the spatial data to the specified file path.
        - If an error occurs during writing, removes any partially created directory.

        Args:
            filepath (str | None): The destination file path where the Zarr file will be written.

        Returns:
            SpatialHandler: The current instance of the SpatialHandler.

        Raises:
            Exception: If no spatial data object is present, no file path is provided, or an error occurs during writing.
        """
        if self.sdata is None:
            raise Exception("No spatial data object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            # Add platform type as metadata. Can use downstream
            self.sdata.tables[self.NORMALIZED_TABLE_NAME].uns["platform"] = self.platform

            # See https://github.com/pydata/xarray/issues/3476#issuecomment-1115045538 for why we need to convert to unicode
            self.sdata.tables[self.NORMALIZED_TABLE_NAME].X = self.sdata.tables[self.NORMALIZED_TABLE_NAME].X.astype("float")

            # Will fail if file already exists
            self.sdata.write(file_path=filepath)
        except Exception as err:
            # remove the directory if it was created
            if os.path.exists(filepath):
                import shutil
                shutil.rmtree(filepath)
            raise Exception("Error occurred while writing to file: " + str(err))
        return self

    def _write_to_h5ad(self, filepath: str | None=None) -> "SpatialHandler":
        """
        Writes the current AnnData object (`self.adata`) to an H5AD file at the specified file path.

        Parameters:
            filepath (str | None): The destination file path where the AnnData object should be written.
                If None, an exception is raised.

        Returns:
            SpatialHandler: Returns self after successfully writing the file.

        Raises:
            Exception: If no AnnData object is present, if no file path is provided, or if an error occurs during writing.
                In case of a write error, any partially written file at the destination path is removed.
        """
        if self.adata is None:
            raise Exception("No AnnData object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            self.adata.write(filename=Path(filepath))
        except Exception as err:
            # remove the file if it was created
            if os.path.exists(filepath):
                os.remove(filepath)
            raise Exception("Error occurred while writing to file: ", err)
        return self

class CoxMxHandler(SpatialHandler):
    """
    Factory class for CoxMx dataset uploads and conversions.

    Standardized names for different files:
    * <dataset_id>_`'anndata.h5ad'`: Counts and metadata file.
    * <dataset_id>_`'cluster_assignment.txt'`: Cluster assignment file.
    * <dataset_id>_`'Metrics.csv'`: Metrics file.
    * <dataset_id>_`'variable_features_clusters.txt'`: Variable features clusters file.
    * <dataset_id>_`'variable_features_spatial_moransi.txt'`: Variable features Moran’s I file.
    """

    @property
    def has_images(self) -> bool:
        return False

    @property
    def coordinate_system(self) -> str:
        return "global"

    @property
    def platform(self) -> str:
        return "coxmx"

    @property
    def img_name(self) -> str | None:
        return None

    def _read_file(self, filepath: str, **kwargs) -> "SpatialHandler":
        return self

    def convert_sdata_to_adata(self, include_images: bool | None=None) -> "SpatialHandler":
        return super()._convert_sdata_to_adata(include_images)

    def write_to_zarr(self, filepath: str | None=None) -> "SpatialHandler":
        return super()._write_to_zarr(filepath)

    def write_to_h5ad(self, filepath: str | None=None) -> "SpatialHandler":
        return super()._write_to_h5ad(filepath)



class CurioHandler(SpatialHandler):
    """
    Factory class for Curio Seeker dataset uploads and conversions.

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
    def has_images(self) -> bool:
        return False

    @property
    def coordinate_system(self) -> str:
        return "global"

    @property
    def platform(self) -> str:
        return "curio"

    @property
    def img_name(self) -> str | None:
        return None

    def _read_file(self, filepath: str, **kwargs) -> "SpatialHandler":
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
                if ".DS_Store" in entry.name or "._" in entry.name:
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
        with gene IDs. We need to pull out the original h5ad file and modify the Moran's I-score file to use the ensembl IDs.
        """

        if h5ad_file is None:
            raise Exception("h5ad file not found in tarball.")

        from geardb import get_dataset_by_id
        organism_id = None
        if not kwargs.get("organism_id"):
            dataset = get_dataset_by_id(kwargs.get("dataset_id"))
            if dataset:
                organism_id = dataset.organism_id

        if not organism_id:
            raise Exception("Organism ID not found in dataset metadata or provided as an argument.")

        # Read in the h5ad file
        adata = ad.read_h5ad(h5ad_file)

        # Add ensemble IDs to the adata.var
        adata = update_adata_with_ensembl_ids(adata, organism_id, "UNMAPPED_")

        # Create mapping dict.
        # Original index name (ensembl_id) was created in add_ensembl_id_to_h5ad_missing_release.py
        gene_symbol_to_ensembl_id = adata.var.reset_index().set_index("gene_symbol").to_dict()["ensembl_id"]
        if spatial_moransi_file is None:
            raise Exception("Moran's I score file not found in tarball.")

        var_features_moransi = pd.read_csv(spatial_moransi_file, sep="\t", header=0)

        # Replace the gene symbols with the ensembl ids
        var_features_moransi.index = var_features_moransi.index.map(gene_symbol_to_ensembl_id)
        var_features_moransi.to_csv(spatial_moransi_file, sep="\t", header=True, index=True, index_label=False)

        # Now are ready to read in to a SpatialData object
        sdata = sdio.curio(tmp_dir)

        # To get the adata equivalent, look at sdata.tables["table"]

        # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
        sdata.tables[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

        # Change obs "cluster" to "clusters" to harmonize
        sdata.tables[self.NORMALIZED_TABLE_NAME].obs = sdata.tables[self.NORMALIZED_TABLE_NAME].obs.rename(columns={"cluster": "clusters"})

        self.sdata = sdata

        # table name should already be "table" for Visium

        self.originalFile = filepath
        return self

    def convert_sdata_to_adata(self, include_images: bool | None = None) -> "SpatialHandler":
        return super()._convert_sdata_to_adata(include_images)

    def write_to_zarr(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_zarr(filepath)

    def write_to_h5ad(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_h5ad(filepath)

class GeoMxHandler(SpatialHandler):
    """
    Code is mostly inspired by https://github.com/LiHongCSBLab/SOAPy/blob/153095a44200a07a73a6a72c9978adfa1581c853/SOAPy_st/pp/all2adata.py#L229
    I wanted to install SOAPy but ran into pip requirement compatibility issues.  For example, we use a later version of AnnData in gEAR than SOAPy does.

    Factory class for GeoMx dataset uploads and conversions.

    Required files:
    * "xlsx" file with information.
      * This Excel file must contain a sheet named "SegmentProperties" with a column named "SegmentDisplayName" which will be used as the cell ID.
      * This Excel file must contain a sheet named "TargetCountMatrix" or "BioProbeCountMatrix" with the counts matrix.

    NOT IMPLEMENTED - Polygon data from XML files

    It seems GeoMx Excel data has gene IDs but can have multiple accessions (none are ensembl IDs),
    so you need to modify add the ensembl IDs to the adata.var.
    You can use add_ensembl_id_to_h5ad_missing_release.py for that, and include the revisted h5ad file in the tarball.

    """

    @property
    def has_images(self) -> bool:
        return False

    @property
    def coordinate_system(self) -> str:
        return "global"

    @property
    def platform(self) -> str:
        return "geomx"

    @property
    def img_name(self) -> str | None:
        return None

    def _read_file(self, filepath: str, **kwargs) -> "SpatialHandler":
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        information_file = None

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if ".DS_Store" in entry.name or "._" in entry.name:
                    continue
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

                if entry.name.endswith(".xlsx"):
                    information_file = filepath

        if information_file is None:
            raise Exception("Excel file containing sample and count information not found in tarball.")

        # Validate the Excel file
        import openpyxl
        wb = openpyxl.load_workbook(information_file)
        if "SegmentProperties" not in wb.sheetnames:
            raise Exception("Excel file must contain a sheet named 'SegmentProperties' with a column named 'SegmentDisplayName'.")
        # Determine if the count matrix is in "TargetCountMatrix" or "BioProbeCountMatrix" and store sheet name as variable
        if "TargetCountMatrix" in wb.sheetnames:
            count_matrix_sheet = "TargetCountMatrix"
        elif "BioProbeCountMatrix" in wb.sheetnames:
            count_matrix_sheet = "BioProbeCountMatrix"
        else:
            raise Exception("Excel file must contain a sheet named 'TargetCountMatrix' or 'BioProbeCountMatrix' with the counts matrix.")

        from geardb import get_dataset_by_id
        organism_id = kwargs.get("organism_id", None)
        if not organism_id:
            dataset = get_dataset_by_id(kwargs.get("dataset_id"))
            if dataset:
                organism_id = dataset.organism_id

        if not organism_id:
            raise Exception("Organism ID not found in dataset metadata or provided as an argument.")

        # Get observation table data
        roi_obs_df = pd.read_excel(
            information_file,
            sheet_name='SegmentProperties',
            index_col='SegmentDisplayName',
            header=0,
        )
        obs = roi_obs_df

        # Get count matrix data
        roi_counts_df = pd.read_excel(
            information_file,
            sheet_name=count_matrix_sheet,
            index_col=0,
            header=0,
        ).T

        if count_matrix_sheet == "BioProbeCountMatrix":
            # If the count matrix is from BioProbeCountMatrix, need to set TargetName values as the column names
            # The sheet has accession IDs but some entries have multiple IDs, so we will map Ensembl IDs later
            roi_counts_df.columns = roi_counts_df.loc["TargetName"].to_numpy()
            roi_counts_df = roi_counts_df.drop("TargetName")

        counts = roi_counts_df.loc[obs.index, :]

        var = pd.DataFrame(index=counts.columns)

        adata = ad.AnnData(counts.values, obs=obs, var=var, uns={}, obsm={})
        adata.obsm['spatial'] = obs.loc[:, ['ROICoordinateX', 'ROICoordinateY']].to_numpy()

        # Sanitize adata.obs so that column names only contain alphanumeric characters, underscores, dots and hyphens.
        adata.obs.columns = adata.obs.columns.str.replace(r'+', '_Plus')    # special case for '+' character
        adata.obs.columns = adata.obs.columns.str.replace(r'[^a-zA-Z0-9_.-]', '_')
        adata.obs.columns = adata.obs.columns.str.replace(r' ', '_')

        # Add ensemble IDs to the adata.var
        adata = update_adata_with_ensembl_ids(adata, organism_id, "UNMAPPED_")

        # Convert to SpatialData object
        sdata = from_legacy_anndata(adata)

        # GeoMx is focused more on spatial bulk rather than spatial single-cell, so we will set the clusters to the index
        sdata.tables[self.NORMALIZED_TABLE_NAME].obs["clusters"] = sdata.tables[self.NORMALIZED_TABLE_NAME].obs.index

        self.sdata = sdata
        self.originalFile = filepath
        return self

    def convert_sdata_to_adata(self, include_images: bool | None = None) -> "SpatialHandler":
        return super()._convert_sdata_to_adata(include_images)

    def write_to_zarr(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_zarr(filepath)

    def write_to_h5ad(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_h5ad(filepath)

class VisiumHandler(SpatialHandler):
    # NOTE: Uploads work but it cannot be used in a spatial panel yet because clusters have not been provided.
    """
    Factory class for Visium dataset uploads and conversions.

    Standardized names for different files:
    * (<dataset_id>_)`'filtered_feature_bc_matrix.h5'`: Counts and metadata file.
    * 'analysis/clustering/gene_expression_graphclust/clusters.csv': Clustering information. Preferable if "Cluster" column has actual labels instead of numbers.
    * 'spatial/tissue_hires_image.png': High resolution image.
    * 'spatial/tissue_lowres_image.png': Low resolution image.
    * 'spatial/scalefactors_json.json': Scalefactors file.
    * 'spatial/tissue_positions_list.csv' (SpaceRanger 1) or 'spatial/tissue_positions.csv' (SpaceRanger 2): Spots positions file.
    * fullres_image_file: (NOT USED) large microscopy image used as input for space ranger.

    Recommended tar command to create tarball:
    `tar cvf <dataset>.tar <dataset>_filtered_feature_bc_matrix.h5 spatial analysis`
    """

    @property
    def has_images(self) -> bool:
        return True

    @property
    def coordinate_system(self) -> str:
        return "downscaled_hires"

    @property
    def platform(self) -> str:
        return "visium"

    @property
    def img_name(self) -> str | None:
        return "spatialdata_hires_image"

    def _read_file(self, filepath: str, **kwargs) -> "SpatialHandler":

        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if ".DS_Store" in entry.name or "._" in entry.name:
                    continue
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)


        clustering_csv_path = "{}/analysis/clustering/gene_expression_graphclust/clusters.csv".format(tmp_dir)
        # If clustering file does not exist, raise an exception
        if not os.path.exists(clustering_csv_path):
            raise Exception("clusters.csv file not found in tarball.")

        # If clustering file does not have "Barcode" and "Cluster" columns, raise an exception
        with open(clustering_csv_path, 'r') as f:
            first_line = f.readline()
            if "Barcode" not in first_line or "Cluster" not in first_line:
                raise Exception("clusters.csv file does not have 'Barcode' and 'Cluster' columns in clusters.csv file in tarball.")



        sdata = sdio.visium(path=tmp_dir, dataset_id="spatialdata")    # Provide a name to standarize downstream usage

        # add clustering information to the vis_sdata.table.obs dataframe
        clustering = pd.read_csv(clustering_csv_path)
        # make barcode as index
        clustering = clustering.set_index('Barcode')
        sdata.tables[self.NORMALIZED_TABLE_NAME].obs['clusters'] = clustering['Cluster'].astype('category')

        # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
        sdata.tables[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

        # currently gene symbols are the index, need to move them to a column
        sdata.tables[self.NORMALIZED_TABLE_NAME].var["gene_symbol"] = sdata.tables[self.NORMALIZED_TABLE_NAME].var.index

        # set the index to the ensembl id (gene_ids)
        sdata.tables[self.NORMALIZED_TABLE_NAME].var = sdata.tables[self.NORMALIZED_TABLE_NAME].var.set_index("gene_ids")

        self.sdata = sdata
        self.originalFile = filepath
        return self

    def convert_sdata_to_adata(self, include_images: bool | None = None) -> "SpatialHandler":
        return super()._convert_sdata_to_adata(include_images)

    def write_to_zarr(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_zarr(filepath)

    def write_to_h5ad(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_h5ad(filepath)

class VisiumHDHandler(SpatialHandler):
    """
    Factory class for Visium HD dataset uploads and conversions.

    Explanation of Space Ranger v3 output here -> https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview#hd-outputs

    Required files that will break the upload if not present:
    /binned_outputs/square_008um/filtered_feature_bc_matrix.h5
    /binned_outputs/square_008um/analysis/clustering/gene_expression_graphclust/clusters.csv
    /binned_outputs/feature_slice.h5
    /binned_outputs/square_008um/spatial/scalefactors_json.json
    /binned_outputs/square_008um/spatial/tissue_hires_image.png
    /binned_outputs/square_008um/spatial/tissue_lowres_image.png

    Recommended tar command to create tarball:
    `tar cvf <dataset>.tar binned_outputs/feature_slice.h5 binned_outputs/square_008um`

    Special note: We have observed that bin sizes finer than 8 microns per pixel will generally have more cells, which can lead to longer analysis times.
    For now, we will attempt to use the "square_008um" binned output.

    """

    table_name = "square_008um"

    @property
    def has_images(self) -> bool:
        return True

    @property
    def coordinate_system(self) -> str:
        return "downscaled_hires"

    @property
    def platform(self) -> str:
        return "visium_hd"

    @property
    def img_name(self) -> str | None:
        return "spatialdata_hires_image"

    def _read_file(self, filepath: str, **kwargs) -> "SpatialHandler":
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
                if ".DS_Store" in entry.name or "._" in entry.name:
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
                                , load_all_images=False  # CytAssist image is not helpful for us.
                                , fullres_image_file=None
                                , bins_as_squares=True
                                )

        # add clustering information to the vis_sdata.table.obs dataframe
        clustering = pd.read_csv(clustering_csv_path)
        # make barcode as index
        clustering = clustering.set_index('Barcode')
        sdata.tables[self.table_name].obs['clusters'] = clustering['Cluster'].astype('category')

        # To get the adata equivalent, look at sdata.tables["table"]
        # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
        sdata.tables[self.table_name].var_names_make_unique()

        # currently gene symbols are the index, need to move them to a column
        sdata.tables[self.table_name].var["gene_symbol"] = sdata.tables[self.table_name].var.index

        # set the index to the ensembl id (gene_ids)
        sdata.tables[self.table_name].var = sdata.tables[self.table_name].var.set_index("gene_ids")

        # Set the table name to the normalized table name
        sdata.tables[self.NORMALIZED_TABLE_NAME] = sdata.tables[self.table_name]

        self.sdata = sdata
        self.originalFile = filepath
        return self

    def convert_sdata_to_adata(self, include_images: bool | None = None) -> "SpatialHandler":
        return super()._convert_sdata_to_adata(include_images)

    def write_to_zarr(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_zarr(filepath)

    def write_to_h5ad(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_h5ad(filepath)

class XeniumHandler(SpatialHandler):
    """
    Factory class for Xenium dataset uploads and conversions.

    Standardized names for different files:
    * (REQ) 'experiment.xenium': File containing specifications.
    * (REQ) 'cell_feature_matrix.h5': File containing cell feature matrix.
    * (REQ) 'cells.parquet': File containing cell metadata.
    * (REQ) 'morphology_focus.ome.tif': File containing morphology focus or a "morphology_focus" directory containing multiple images.
    * 'nucleus_boundaries.parquet': Polygons of nucleus boundaries.
    * 'cell_boundaries.parquet': Polygons of cell boundaries.
    * 'transcripts.parquet': File containing transcripts.
    * 'cells.zarr.zip': Zarr file containing cell and nucleus label data (NOT USING FOR NOW)
    * 'analysis/clustering/gene_expression_graphclust/clusters.csv': Clustering information. Preferable if "Cluster" column has actual labels instead of numbers.

    Currently we are not using the cells.zarr.zip file because of memory issues when converting the resulting cell labels as part of the AnnData conversion process.

    More information on Xenium Ranger outputs can be found here: https://www.10xgenomics.com/support/software/xenium-ranger/latest/analysis/outputs/XR-output-overview
    """

    @property
    def has_images(self) -> bool:
        return True

    @property
    def coordinate_system(self) -> str:
        return "global"

    @property
    def platform(self) -> str:
        return "xenium"

    @property
    def img_name(self) -> str | None:
        return "morphology_focus"

    def _read_file(self, filepath: str, **kwargs) -> "SpatialHandler":
        # Get tar filename so tmp directory can be assigned
        tar_filename = filepath.rsplit('/', 1)[1].rsplit('.')[0]
        tmp_dir = '/tmp/' + tar_filename

        if os.path.isdir(tmp_dir):
            # Remove any existing directory
            os.system("rm -rf {}".format(tmp_dir))

        # settings to enable or disable based on if a file is present in the uploaded tarball
        include_raster_labels = False
        cell_boundaries_present = False
        nucleus_boundaries_present = False
        transcripts_present = False

        with tarfile.open(filepath) as tf:
            for entry in tf:
                # Skip any BSD tar artifacts, like files that start with ._ or .DS_Store
                if ".DS_Store" in entry.name or "._" in entry.name:
                    continue
                # Extract file into tmp dir
                filepath = "{0}/{1}".format(tmp_dir, entry.name)
                tf.extract(entry, path=tmp_dir)

                if entry.name == "cells.zarr.zip":
                    include_raster_labels = True
                if entry.name == "cell_boundaries.parquet":
                    cell_boundaries_present = True
                if entry.name == "nucleus_boundaries.parquet":
                    nucleus_boundaries_present = True
                if entry.name == "transcripts.parquet":
                    transcripts_present = True

        # If clustering file does not exist, raise an exception
        clustering_csv_path = "{}/analysis/clustering/gene_expression_graphclust/clusters.csv".format(tmp_dir)
        if not os.path.exists(clustering_csv_path):
            raise Exception("clusters.csv file not found in tarball.")

        # If clustering file does not have "Barcode" and "Cluster" columns, raise an exception
        with open(clustering_csv_path, 'r') as f:
            first_line = f.readline()
            if "Barcode" not in first_line or "Cluster" not in first_line:
                raise Exception("clusters.csv file does not have 'Barcode' and 'Cluster' columns in clusters.csv file in tarball.")

        sdata = sdio.xenium(tmp_dir
                            , cells_labels=False # Avoid adding polygons to SpatialData object (for now due to out-of-memory issues)
                            , nucleus_labels=False
                            , cell_boundaries=cell_boundaries_present
                            , nucleus_boundaries=nucleus_boundaries_present
                            , transcripts=transcripts_present
                            , cells_as_circles=True  # Table is associated with the cells instead of the nuclei (faster performance)
                            , morphology_mip=False   # Using the morphology_focus image instead
                            )

        # In code, it seems that the Xenium reader is supposed to set the index to the "barcodes" column
        # But this column is not found, so we need to manually replace with "cell_id"
        sdata.tables[self.NORMALIZED_TABLE_NAME].obs["Barcode"] = sdata.tables[self.NORMALIZED_TABLE_NAME].obs["cell_id"]
        sdata.tables[self.NORMALIZED_TABLE_NAME].obs = sdata.tables[self.NORMALIZED_TABLE_NAME].obs.set_index("Barcode")

        # Change annotation target from "cell_circles" to "cell_labels"
        #sdata["table"].obs["region"] = "cell_labels"
        #sdata.set_table_annotates_spatialelement(
        #    table_name="table", region="cell_labels", region_key="region", instance_key="cell_labels"
        #)

        # add clustering information to the vis_sdata.table.obs dataframe
        clustering = pd.read_csv(clustering_csv_path)
        # make barcode as index
        clustering = clustering.set_index('Barcode')
        sdata.tables[self.NORMALIZED_TABLE_NAME].obs['clusters'] = clustering['Cluster'].astype('category')

        # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
        sdata.tables[self.NORMALIZED_TABLE_NAME].var_names_make_unique()

        # currently gene symbols are the index, need to move them to a column
        sdata.tables[self.NORMALIZED_TABLE_NAME].var["gene_symbol"] = sdata.tables[self.NORMALIZED_TABLE_NAME].var.index

        # set the index to the ensembl id (gene_ids)
        sdata.tables[self.NORMALIZED_TABLE_NAME].var = sdata.tables[self.NORMALIZED_TABLE_NAME].var.set_index("gene_ids")

        self.sdata = sdata
        self.originalFile = filepath
        return self

    def convert_sdata_to_adata(self, include_images: bool | None = None) -> "SpatialHandler":
        return super()._convert_sdata_to_adata(include_images)

    def write_to_zarr(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_zarr(filepath)

    def write_to_h5ad(self, filepath: str | None = None) -> "SpatialHandler":
        return super()._write_to_h5ad(filepath)


### Helper constants

SPATIALTYPE2CLASS = {
    #"cosmx": CoxMxHandler,
    "curio": CurioHandler,
    "geomx": GeoMxHandler,
    "visium": VisiumHandler,
    "visium_hd": VisiumHDHandler,
    "xenium": XeniumHandler
}