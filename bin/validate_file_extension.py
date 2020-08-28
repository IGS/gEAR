#!/opt/bin/python3

"""
For gEAR expression data files. Checks that the file extension matches the actual
text delimiter.

Reads the 1st line, and if delimiter and file extension agree. Script closes file
and moves to next file. However, if the delimiter and file extension disagree,
the file extension is then changed to match the delimiter.

Example: '\t'  .tab ==> No actions
         ','   .tab ==> Change extension to .csv
         ','   .csv ==> No actions
         ','   .tab ==> Change extension to .tab
"""


import argparse
import csv
import os

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description="Validates file extension agrees with file's text delimiter")
    parser.add_argument('-f', '--input_file', type=str, required=False, help='Path to an input file' )
    parser.add_argument('-l', '--input_list', type=str, required=False, help='Path to an input list file' )
    args = parser.parse_args()

    files = list()

    if args.input_file is not None:
        files.append(args.input_file)
    elif args.input_list is not None:
        files = utils.read_list_file(args.input_list)
    else:
        raise Exception("ERROR: You must pass either -i or -l options")

    changed_num = 0
    for file in files:
        if file.lower().endswith( ('.svg', '.png', '.jpeg', '.jpg') ):
            #skip image files
            continue


        print("Processing file: {0}".format(file))

        file_extension = os.path.splitext(file)[1]
        print("... file extension: {0}".format(file_extension))

        fh = open(file)

        #Check if TAB file is tab delimited
        if file_extension == '.tab':
            try:
                reader = csv.reader(fh, delimiter='\t')
                for line in reader:
                    if len(line) > 1:
                        print("... no changed needed.")
                        break
                    else:
                        raise Exception("... Exception: Incorrect file extension.")
            except:
                new_file = file.replace(file_extension, '.csv')
                os.rename(file, new_file)
                print("... extension changed to CSV")

                changed_num += 1

        else:
            try:
                reader = csv.reader(fh)
                for line in reader:
                    if len(line) > 1:
                        print("... no changed needed.")
                        break
                    else:
                        raise Exception("... Exception: Incorrect file extension.")
            except:
                new_file = file.replace(file_extension, '.tab')
                os.rename(file, new_file)
                print("... extension changed to TAB")

                changed_num += 1

        fh.close()

    print("\nFinished.\n{0} Files renamed.".format(changed_num))


if __name__ == '__main__':
    main()
