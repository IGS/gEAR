#!/usr/bin/env python3

"""
This script intends to mimic the upload process for a gEAR dataset metadata file.


Input
-----

-d --test_directory = Path to directory containing xlsx files to be used as test metadata files.
-f --test_file      = Path to xlsx file to be used as test metadata file

Output
------

Returns STDOUT:


"""

#enable if tests/base_template.xlsx fails
import traceback

import argparse
parser = argparse.ArgumentParser( description='Put a description of your script here')
parser.add_argument('-d', '--test_directory', type=str, required=False, help='Path to directory containing xlsx test metadata files.' )
parser.add_argument('-f', '--test_file', type=str, required=False, help='Path to xlsx test metadata file.' )
args = parser.parse_args()
test_directory = args.test_directory
test_file = args.test_file


import os, sys
dolley_path = '/home/dolley/gear/lib'
jorvis_path = '/home/jorvis/git/gEAR/lib'

gear_cgi_path = os.path.dirname(os.path.realpath(__file__))
gear_www_path = os.path.dirname(gear_cgi_path)
gear_lib_path = gear_www_path.replace("www", "lib")

if os.path.isdir(dolley_path):
    sys.path.append(dolley_path)
if os.path.isdir(jorvis_path):
    sys.path.append(jorvis_path)

# Should always be a valid path
if os.path.isdir(gear_lib_path):
    sys.path.append(gear_lib_path)

from gear.metadatauploader import MetadataUploader as uploader
# from gear.metadatavalidator import MetadataValidator as validator

files = list()
if test_directory:
    import glob
    files = glob.glob(test_directory + '/*.xlsx')
else:
    files.append(test_file)

number_tests_performed = 0
number_tests_failed = 0
tests_failed = list()
for f in files:
    number_tests_performed += 1
    try:
        # Initize metadatauploader
        uploader = uploader()

        # Read in metadata file
        uploader._read_file(f)

        # Validate inputs
        #  - checks length in all fields
        #  - checks for 200 status code for pubmed and geo ids
        #  - TODO need to add more validation checks
        uploader._validate_values()

        # Prints the dataframe with 'is_valid' column added
        # print(uploader.metadata)

        #TODO: Are required fields populated?
        # What are the required fields?

        #TODO: Upload info into MySQL (dataset table)


        # Write metadata to json file.
        uploader._write_to_json(filepath='/home/dolley/gear/tests/TEST_metadata_output.json')

    except Exception as e:
        tests_failed.append(f)
        number_tests_failed += 1
        print('ERROR:\t{0}'.format(e))
        traceback.print_exc() # enable to see where base_template.xlsx failed
        continue

print("\n{0}\nTotal Number of Tests Run:\t{1}\nTotal Number of Tests Failed:\t{2}" \
        .format('-'*40, number_tests_performed, number_tests_failed))
print("\tTests Failed:\n\t\t{0}".format('\n\t\t'.join(sorted(tests_failed))))
