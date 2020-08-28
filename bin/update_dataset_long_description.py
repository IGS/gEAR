#!/usr/bin/env python3

"""
There are cases where we want a detailed or HTML-based long description, which isn't really
supported by the Dataset Manager currently.  

This script lets you read a file and makes its contents the dataset.ldesc of a given
dataset ID.  This file is usually an HTML snippet.

It minifies the HTML before uploading.

"""

import os
import sys
import argparse
import htmlmin

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

import geardb


def main():
    parser = argparse.ArgumentParser( description='')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-id', '--dataset_id', type=str, required=True, help='ID of the dataset to update' )
    args = parser.parse_args()

    comments = htmlmin.minify(
        open(args.input_file).read(),
        remove_empty_space=True
    )
    
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    qry = "UPDATE dataset SET ldesc = %s WHERE id = %s"
    cursor.execute(qry, (comments, args.dataset_id))

    cursor.close()
    cnx.commit()
    print("Finished.\n Dataset updated.")



if __name__ == '__main__':
    main()
