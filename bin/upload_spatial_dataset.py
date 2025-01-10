#!/opt/bin/python3

import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear import spatialuploader
from add_excel_metadata_to_db import add_metadata_to_db

SPATIALTYPE2CLASS = {
    "visium": spatialuploader.VisiumUploader,
    "visiumhd": spatialuploader.VisiumHDUploader,
    "xenium": spatialuploader.XeniumUploader,
    "curio": spatialuploader.CurioUploader
}

def main():

    parser = argparse.ArgumentParser( description='Loads both spatial data and metadata information' )

    # Spatial file
    parser.add_argument('-f', '--input_file', type=str, required=True, help='Path to an input spatial tarball to be read' )
    parser.add_argument('-t', '--type', type=str, required=True, help='Type of spatial data being uploaded')
    # metadata file
    parser.add_argument('-i', '--metadata_file', type=str, required=False, help='Path to an metadata XLS file to be read' )
    parser.add_argument('-oi', '--owner_id', type=int, required=False, help='Numerical user ID who will own the dataset' )
    parser.add_argument('-s', '--skip_validation', action='store_true', help='Path to an output file to be created' )
    args = parser.parse_args()

    # Current spatial data types supported are : "visium", "visiumhd", "xenium", "curio"
    if args.type not in ["visium", "visiumhd", "xenium", "curio"]:
        print("Invalid or unsupported spatial data type")
        sys.exit(1)


    if args.metadata_file:
        dataset_uid, share_uid = add_metadata_to_db(args.metadata_file, args.owner_id, args.skip_validation)

        print("Added database entry for new dataset: {0}".format(dataset_uid))
        print("\tShare ID: {0}".format(share_uid))

if __name__ == '__main__':
    main()

