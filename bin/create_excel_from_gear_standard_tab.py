#!/opt/bin/python3

"""
This takes an input tab file and creates a three-tabbed Excel spreadsheet for import into gEAR.

Input tab file should have columns:

1. Ensembl_id
2. Gene symbol
3-...  Expression columns

"""

import argparse
import os
import xlsxwriter

def main():
    parser = argparse.ArgumentParser( description='Create a tabbed Excel file from gEAR standard tab')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument("-lf", "--large_file", help="Pass this to flush memory regularly and avoid running out of RAM on large files", action="store_true")
    args = parser.parse_args()

    if args.large_file:
        workbook = xlsxwriter.Workbook(args.output_file, {'constant_memory': True})
    else:
        workbook = xlsxwriter.Workbook(args.output_file)
    
    expression_sheet = workbook.add_worksheet('expression')
    genes_sheet = workbook.add_worksheet('genes')
    obs_sheet = workbook.add_worksheet('observations')

    line_number = 0

    for line in open(args.input_file):
        line_number += 1
        line = line.rstrip()
        cols = line.split("\t")

        expression_sheet.write(line_number - 1, 0, cols[0])

        col_num = 1
        for c in cols[2:]:
            expression_sheet.write(line_number - 1, col_num, c)
            col_num += 1
        
        if line_number == 1:
            obs = cols[2:]
            obs_line_num = 0
            obs_sheet.write(obs_line_num, 0, 'observations')

            for item in obs:
                obs_line_num += 1
                obs_sheet.write(obs_line_num, 0, item)

            genes_sheet.write(line_number - 1, 0, 'genes')
            genes_sheet.write(line_number - 1, 1, 'gene_symbol')
        else:
            var = cols[0:2]
            genes_sheet.write(line_number - 1, 0, var[0])
            genes_sheet.write(line_number - 1, 1, var[1])

    workbook.close()

if __name__ == '__main__':
    main()







