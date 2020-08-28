#!/usr/bin/env python3

"""
Converts Excel gEAR metadata file to JSON.  Only the first two columns of the 
spreadsheet are considered and the first row is assumed to be headers.

Example input file:

https://umgear.org/user_templates/metadata_template.xlsx
"""

import argparse
import json
import os
import sys

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)

from gear.metadata import Metadata

def main():
    parser = argparse.ArgumentParser( description='Converts Excel gEAR metadata file to JSON')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-s', '--skip_validation', action='store_true', help='Path to an output file to be created' )
    args = parser.parse_args()

    metadata = Metadata(file_path=args.input_file)

    if not args.skip_validation:
        metadata.validate()

    metadata.write_json(file_path=args.output_file)


if __name__ == '__main__':
    main()







