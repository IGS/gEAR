#!/usr/bin/env python3

"""

Renames the column headers of a singlecell dataset TAB file
2 Input files
    1) Header (Gene barcode) and Cell type file - 2 column tab-delimited file, example of a row:
            AAACATACACGACT	Inner Hair Cells
    2) Single Expression Data file
"""

import argparse
import csv
import json
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    # output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to file containing single cell header info' )
    parser.add_argument('-e', '--expression_file', type=str, required=True, help='Path to file containing single cell expression data' )
    parser.add_argument('-o', '--output_file_name', type=str, required=True, help='Name of output file' )
    args = parser.parse_args()

    header_info_filename = args.input_file
    expression_filename = args.expression_file
    output_name = "./{0}".format(args.output_file_name)

    line_num = 0

    f_info = open(header_info_filename)
    if header_info_filename.endswith('.tab'):
        reader = csv.reader(f_info, delimiter='\t')
    else:
        reader = csv.reader(f_info)

    gene_cell_pairs = list()
    cell_names = list()

    # Get the Gene barcodes and Cell Type Names
    for cols in reader:
        if line_num == 0:
            #Skip column headers
            line_num += 1
            continue
        else:
            # Example of columns to capture:  AAACATACACGACT	Inner Hair Cells
            gene_barcode = cols[0]
            cell_name = cols[1]

            if not cell_name in cell_names:
                cell_names.append(cell_name)
            gene_cell_pairs.append({'gene_barcode': gene_barcode, 'cell_name': cell_name})

            line_num += 1


    # Format the cell_names into gEAR headers
    for name in cell_names:
        celltype_count = 0
        for item in gene_cell_pairs:
            if name == item['cell_name']:
                celltype_count += 1
                formatted_name = name.replace(' ', '_') + '--Cell_' + str(celltype_count)
                item['cell_name'] = formatted_name

    f_info.close()

    f_expression = open(expression_filename)
    if expression_filename.endswith('.tab'):
        reader = csv.reader(f_expression, delimiter='\t')
    else:
        reader = csv.reader(f_expression)

    # Expression File....
    ex_line_num = 0
    output_file = open(output_name, "w")
    output_writer = csv.writer(output_file, delimiter='\t')
    for row in reader:
        if ex_line_num != 0:
            ex_line_num += 1
            output_writer.writerow(row)
        else:
            print(row)
            new_headers_row = ['Gene_symbol']
            for barcode in row:
                for item in gene_cell_pairs:
                    if barcode == item['gene_barcode']:
                        new_headers_row.append(item['cell_name'])

            output_writer.writerow(new_headers_row)
            ex_line_num += 1

    print('\nDONE!')

    f_expression.close()


if __name__ == '__main__':
    main()
