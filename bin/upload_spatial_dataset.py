#!/opt/bin/python3

"""
upload_spatial_dataset.py - Upload spatial data to the correct location in gEAR dataset storage.

Assumes metadata is already in the database (use add_excel_metadata_to_db.py).

Also assumes that the spatial data has ensembl IDs present for the given platform. Certain platforms, like Curio Seeker,
may only provide gene IDs, so you need to modify the included h5ad file to include ensembl IDs in adata.var.
You can use add_ensembl_id_to_h5ad_missing_release.py for that.

"""

import argparse
import os
import sys
from pathlib import Path

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
from gear import spatial_handler

DEST_DIRPATH = Path(__file__).resolve().parent.parent / "www" / "datasets" / "spatial"
OUTPUT_SUFFIX = ".zarr"

def main():

    parser = argparse.ArgumentParser( description='Loads both spatial data and metadata information' )

    # Spatial file
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input spatial tarball to be read' )
    parser.add_argument('-t', '--type', type=str, required=True, help='Type of spatial data being uploaded')
    parser.add_argument('-d', '--dataset_id', type=str, required=True, help='Numerical dataset ID to be used' )
    parser.add_argument('-org', '--organism_id', type=int, required=False, help='Organism ID of the spatial data. Only needed if data requires remapping of gene symbols to ensembl IDs and metadata is not present in database.' )
    args = parser.parse_args()

    # Ensure the spatial data type is supported
    if args.type not in spatial_handler.SPATIALTYPE2CLASS.keys():
        print("Invalid or unsupported spatial data type")
        print("Supported types: {0}".format(spatial_handler.SPATIALTYPE2CLASS.keys()))
        sys.exit(1)
    print("Processing spatial data of type: {0}".format(args.type))
    print("Input file: {0}".format(args.input_file))
    sp_class = spatial_handler.SPATIALTYPE2CLASS[args.type]()

    sp_class._read_file(args.input_file, organism_id=args.organism_id, dataset_id=args.dataset_id)

    output_filename = args.dataset_id

    output_path = DEST_DIRPATH / (output_filename + OUTPUT_SUFFIX)
    print("Writing to {0}".format(output_path))
    sp_class._write_to_zarr(filepath=output_path)

if __name__ == '__main__':
    main()

