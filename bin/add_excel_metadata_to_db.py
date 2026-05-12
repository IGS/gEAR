#!/usr/bin/env python3

"""
Takes Excel gEAR metadata file and inserts it as a row in the dataset table.

Example input file:

https://umgear.org/user_templates/metadata_template.xlsx
"""

import argparse
import json
import os
import sys
import uuid

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)

from gear.metadata import Metadata

def add_metadata_to_db(input_file, owner_id, skip_validation=False):
    dataset_uid = str(uuid.uuid4())
    share_uid = str(uuid.uuid4()).split('-')[0]

    metadata = Metadata(file_path=input_file)

    if not skip_validation:
        metadata.validate()

    metadata.add_field_value('dataset_uid', dataset_uid)
    metadata.add_field_value('share_uid', share_uid)
    metadata.add_field_value('owner_id', owner_id)

    metadata.save_to_mysql(status='completed')

    return dataset_uid, share_uid

def main():
    parser = argparse.ArgumentParser( description='Parses Excel gEAR metadata file and stores in DB')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-oi', '--owner_id', type=int, required=True, help='Numerical user ID who will own the dataset' )
    parser.add_argument('-s', '--skip_validation', action='store_true', help='Path to an output file to be created' )
    args = parser.parse_args()

    dataset_uid, share_uid = add_metadata_to_db(args.input_file, args.owner_id, args.skip_validation)

    print("Added database entry for new dataset: {0}".format(dataset_uid))
    print("\tShare ID: {0}".format(share_uid))


if __name__ == '__main__':
    main()







