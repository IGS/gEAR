# Uploading Spatial Transcriptomics Datasets

This document is meant to be a guide on how to upload a spatial transcriptomics dataset before the uploader tool supports spatial datasets.

## Currently supported platforms

* Curio Seeker
  * No image support
* Visium
  * Have not tested yet
* Visium HD
* Xenium
  * Have not tested yet

## Metadata prep

Use the metadata_template.xlsx to store metadata.  Where the dataset type dropdown is located, delete the dropdown cell and add "spatial" in order to bypass the dropdown requirements and get it into the right type.  Right now, I haven't looked into adding "spatial" into the list of choices yet for validation.

To add the Excel metadata to the database, run the following command:

`cd <gear-root>/bin; <python_path>/python add_excel_metadata_to_db.py -i <metadata>.xlsx -oi <owner_id>`

Note that the python path and that the owner ID for the dataset will both be variable depending on your environment.  When you run this command, you shall get both a dataset_id value and a share_id value.  Take note of the dataset ID as you will need it for later.

## File prep

Our data format type of choice for spatial data is ZARR, which is a data storage specification for storage of large N-dimensional typed arrays (tensors). Spatial data can come from a variety of platforms, and we have utilities in `<gear-path>/lib/gear/spatialuploader.py` to harmomize some of these by converting into a .zarr file. This module also has a dictionary that maps supported formats to the right uploader class, and this can be used to see what platforms are supported.

If the dataset is already in .zarr format (converted previously and copying to a new server, for instance), you can just move the dataset by doing `mv <dataset>.zarr <gear-path>/www/datasets/spatial/<dataset-id>.zarr` where the dataset_id is the uuid4 you got from uploading the metadata.

### Getting Ensembl IDs

UPDATED: This step is now integrated into the Curio and GeoMx dataset uploaders. Either pass in the organism ID directly to `upload_spatial_dataset.py` as an argument or let it be read from existing uploaded metadata.

It seems that the 10x platforms (i.e. Visium and Xenium) have Ensembl IDs in the output that can be converted into the format we need... as the index for that DataFrame.  However for other platforms, like Curio Seeker, only gene symbols are present in the output files.  So we need to use other scripts to determine the Ensembl ID annotations.

Run the following to map gene symbols to Ensembl IDs stored in the database (this example is for Curio which outputs an H5AD file):

`cd <gear-root>/bin; <python_path>/python add_ensembl_id_to_h5ad_missing_release.py -i <unmapped-dataset>.h5ad -o <new-dataset>.h5ad -org <organism_id> -idp UNMAPPED`

For this example, you need the organism_id for the correct organism, found in the database.

### Uploading the file

Before uploading the spatial dataset, ensure that it is in a tarball format and it has the requisitie files.  You can check the docuemntation for each class in `<gear-path>/lib/gear/spatialuploader.py` to learn what files are required at minimum and where they should be in the tarball. In general though, the output contents and structure for the platform should be preserved within the tarball.

This it the command to convert and upload the spatial dataset to the correct location in gEAR.:

`cd <gear-root>/bin; <python_path>/python upload_spatial_dataset.py -i <dataset>.tar -t <platform> -d <dataset-id>`

When this is finished, the final location will be printed to screen, meaning things are complete. If the upload fails, the extract tarball content directory will be deleted upon rerun, so that a previous set of contents will not influence the current run.  You can view the extracted tarball contents in your `/tmp` area.