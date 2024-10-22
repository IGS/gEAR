
import os, json, sys
import tarfile
import pandas as pd

from gear.datasetuploader import FileType

import spatialdata_io as sdio
from spatialdata_io.experimental import to_legacy_anndata



class VisiumUploader(FileType):
    """
    Called by datasetuploader.py (factory) when a Visium or Visium-HD dataset is going to be uploaded

    Standardized names for different files:
    * (<dataset_id>_)`'filtered_feature_bc_matrix.h5'`: Counts and metadata file.
    * 'spatial/tissue_hires_image.png': High resolution image.
    * 'spatial/tissue_lowres_image.png': Low resolution image.
    * 'scalefactors_json.json': Scalefactors file.
    * 'tissue_positions_list.csv' (SpaceRanger 1) or 'tissue_positions.csv' (SpaceRanger 2): Spots positions file.
    * fullres_image_file: large microscopy image used as input for space ranger.

    In addition, we may encounter 'tissue_positions.parquet' instead of the CSV file. If this file is found it will be converted into 'tissue_positions_list.csv' (Space Ranger 1.3 output)
    All files need to be in the same directory format as Space Ranger output. This should be tarballed in order to upload to gEAR.

    Special note: We have observed that bin sizes finer than 8 microns per pixel will generally have more cells, which can lead to longer analysis times.
    For now, we will add a warning and ask the user to downsample before uploading. But will not prevent the upload. This warning should also be added to other tools that use the same data.
    A quick way to determine bin size is to look at the 'scalefactors_json.json' file. The 'bin_size_um' key will give you the bin size in microns.

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

            # If parquet file exists, convert to Space Ranger csv output
            parquet_file = "{}/spatial/tissue_positions.parquet".format(tmp_dir)
            if os.path.exists(parquet_file):
                df = pd.read_parquet("{}/spatial/tissue_positions.parquet".format(tmp_dir))
                df.to_csv("{}/spatial/tissue_positions_list.csv".format(tmp_dir), index=False)

            # look for the "bin_size_um" key in the scalefactors_json.json file
            # if the bin size is less than 8 microns per pixel, we will assign a flag to this dataset
            scalefactors_file = "{}/scalefactors_json.json".format(tmp_dir)
            if os.path.exists(scalefactors_file):
                with open(scalefactors_file) as f:
                    data = json.load(f)
                    bin_size = data["bin_size_um"]
                    if bin_size < 8:
                        self.large_dataset_flag = True
                        print(f"Bin size for uploaded file {tar_filename} is less than 8 microns per pixel. This may lead to longer analysis times. Please consider downsampling before uploading.", file=sys.stderr)

            # TODO: have branching condition to handle Visium-HD datasets
            sdata = sdio.visium(tmp_dir, dataset_id="spatialdata")    # Provide a name to standarize downstream usage

            # To get the adata equivalent, look at sdata.tables["table"]

            # The Space Ranger h5 matrix has the gene names as the index, need to move them to a column and set the index to the ensembl id
            #sdata.var_names_make_unique()

            # currently gene symbols are the index, need to move them to a column
            #sdata.var["gene_symbol"] = sdata.var.index

            # set the index to the ensembl id (gene_ids)
            #sdata.var.set_index("id", inplace=True)

            self.sdata = sdata

            adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires")

            # Apply AnnData obj and filepath to uploader obj
            self.adata = adata
            self.originalFile = filepath
            return self

    def _write_to_zarr(self, filepath=None):
        if self.sdata is None:
            raise Exception("No spatial data object present to write to file.")
        if filepath is None:
            raise Exception("No destination file path given. Provide one to write file.")
        try:
            self.sdata.write(filepath)
        except Exception as err:
            raise Exception("Error occurred while writing to file: ", err)

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
