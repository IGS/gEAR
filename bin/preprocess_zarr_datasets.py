"""
For each zarr dataset in www/datasets/spatial, do the following preprocessing:

1. Create a SpatialHandler object based on the dataset type
2. Run spatial_obj.subset_sdata()
3. Run spatial_obj.scale_and_translate_images()
4. Run spatial_obj.merge_centroids_with_obs()
5. Run the single-cell workbench steps on the spatial_obj.tables["table"] (AnnData object) using default parameters:
    1. Run scanpy.pp.highly_variable_genes()
    2. Run scanpy.pp.pca()
    3. Run scanpy.pp.neighbors()
    4. Run scanpy.tl.umap()
6. Save the processed SpatialData object back to zarr with "_new" appended to the original dataset name
7. Replace the original dataset with the new processed dataset

"""

from pathlib import Path
import scanpy as sc
import spatialdata as sd

CURRENT_DIR = Path(__file__).parent.resolve()

libdir = CURRENT_DIR.parent / "lib"
import sys
sys.path.insert(0, str(libdir))
import gear.spatialhandler as sh

ZARR_PATH = CURRENT_DIR.parent / "www" / "datasets" / "spatial"
ZARR_EXTENSION = ".zarr"

def main():
    # Iterate through each zarr dataset in the specified directory
    for zarr_dataset in ZARR_PATH.glob(f"*{ZARR_EXTENSION}"):
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

        # Step 2: Run spatial_obj.subset_sdata()
        spatial_obj = sh.SPATIALTYPE2CLASS[platform]()
        spatial_obj.sdata = sdata
        spatial_obj.subset_sdata()

        # Step 3: Run spatial_obj.scale_and_translate_images()
        spatial_obj = spatial_obj.scale_and_translate_sdata()

        # Step 4: Run spatial_obj.merge_centroids_with_obs()
        # The SpatialData object table should have coordinates, but they are not translated into the image space
        # Each observation has an associated polygon "shape" in the image space, and we can get the centroid of that shape
        spatial_obj.merge_centroids_with_obs()

        # Step 5: Run the single-cell workbench steps on the spatial_obj.tables["table"] (AnnData object) using default parameters
        adata = spatial_obj.sdata.tables["table"]

        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        adata.var_names_make_unique()

        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        spatial_obj.sdata.tables["table"] = adata

        # Step 6: Save the processed SpatialData object back to zarr with "_new" appended
        new_zarr_path = zarr_dataset.parent / f"{dataset_name}_new{ZARR_EXTENSION}"
        spatial_obj.save_to_zarr(new_zarr_path)

        # Step 7: Replace the original dataset with the new processed dataset
        zarr_dataset.unlink()  # Remove original dataset
        new_zarr_path.rename(zarr_dataset)  # Rename new dataset to original name

        print(f"---Finished processing dataset: {dataset_name}")

if __name__ == "__main__":
    main()
