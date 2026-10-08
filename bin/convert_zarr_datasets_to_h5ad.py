"""
For each zarr dataset in www/datasets/spatial, do the following preprocessing:

1. Create a SpatialHandler object based on the dataset type
2. Run spatial_obj.convert_sdata_to_adata()
3. Write to h5ad file with the same name as the zarr dataset in the parent directory

"""

from pathlib import Path
import spatialdata as sd

CURRENT_DIR = Path(__file__).parent.resolve()

libdir = CURRENT_DIR.parent / "lib"
import sys
sys.path.insert(0, str(libdir))
import gear.spatialhandler as sh

ZARR_DIR = CURRENT_DIR.parent / "www" / "datasets" / "spatial"
H5AD_DIR = ZARR_DIR.parent
ZARR_EXTENSION = ".zarr"

def main():
    # Iterate through each zarr dataset in the specified directory
    for zarr_dataset in ZARR_DIR.glob(f"*{ZARR_EXTENSION}"):
        dataset_name = zarr_dataset.stem
        print(f"---Processing dataset: {dataset_name}")

        # Step 1: Create a SpatialHandler object based on the dataset type
        sdata = sd.read_zarr(zarr_dataset)
        try:
            platform = sdata.tables["table"].uns["platform"]
        except KeyError:
            raise ValueError("No platform information found in the dataset")

        # Ensure the spatial data type is supported
        if platform not in sh.SPATIALTYPE2CLASS.keys():
            raise ValueError(
                f"Invalid or unsupported spatial data type: {platform}"
            )

        spatial_obj = sh.SPATIALTYPE2CLASS[platform]()
        spatial_obj.sdata = sdata

        # Step 2: Run spatial_obj.convert_sdata_to_adata()
        try:
            spatial_obj.convert_sdata_to_adata()
        except Exception as e:
            print(f"Error converting spatial data to AnnData: {e}")
            raise

        # Step 3: Write to h5ad file with the same name as the zarr dataset in the parent directory
        h5ad_path = H5AD_DIR / f"{dataset_name}.h5ad"
        spatial_obj.adata.write_h5ad(h5ad_path)
        spatial_obj.adata.file.close()

        print(f"---Saved h5ad file to: {h5ad_path}")


if __name__ == "__main__":
    main()
