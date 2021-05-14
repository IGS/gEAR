#!/opt/bin/python3

"""
Current data structure returned:

result = {'filename':filename,
          'rows':[],   'missing':[],
          'success':1, 'row_count':200
         }

Each 'row' entry there are just the tabs from the data file

# Example of secondary gene symbol: Dnah2

"""

import cgi
import json
import csv
import os
import sys
import re

sys.path.append('../..')
import loaderutils

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    PREVIEW_ROW_LIMIT = 10

    print('Content-Type: application/json\n\n')

    # added buffered=True to prevent error: 'unread result found'
    # http://stackoverflow.com/a/33632767/2900840
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    gene_cache = loaderutils.cache_genes_by_gene_symbol(cursor, lower=True)
    primary_gene_cache = loaderutils.cache_genes_by_primary_gene_symbol(cursor, lower=True)

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    filename = form.getvalue('filename')

    result = {'filename':filename, 'rows':[], 'missing':[], 'secondary':[], 'success':1, 'row_count':0, 'col_count':0}
    line_num = 0
    gene_symbol_col = None
    preview_column_idxs = list()

    fh = open("../uploads/files/{0}".format(filename))

    # detect delimiter before opening file
    # source: https://docs.python.org/3.5/library/csv.html#csv.Sniffer
    dialect = csv.Sniffer().sniff( fh.read(1024) )
    fh.seek(0)
    # reader = csv.reader(fh, delimiter=dialect.delimiter)
    reader = csv.reader(fh, dialect)

    primary_genes = list()
    skipped_secondary = list()

    for cols in reader:
        # a minimum of two columns is required
        if len(cols) < 2: continue

        data_cols = list()

        if line_num == 0:
            col_num = 0

            for col in cols:
                if not col.startswith('#'):
                    if ('gene' in col.lower() or 'symbol' in col.lower()) and \
                    ('type' not in col.lower() and 'source' not in col.lower()):
                        gene_symbol_col = col_num

                    preview_column_idxs.append(col_num)
                    data_cols.append(col)

                    col_num += 1

            result['col_count'] = col_num
            result['rows'].append(data_cols)

        elif gene_symbol_col is not None:
            # This means we're on a data row, and the gene symbol column has been stored.
            #  We need to check this column against the cached gene symbols
            if cols[gene_symbol_col].lower() in gene_cache:

                if cols[gene_symbol_col].lower() in primary_gene_cache:
                    primary_genes.append(cols[gene_symbol_col]) #track the genes kept

                    if len(result['rows']) < PREVIEW_ROW_LIMIT:
                        col_idx = 0
                        for col in cols:
                            if col_idx in preview_column_idxs:
                                data_cols.append(col)

                            col_idx += 1

                        result['rows'].append(data_cols)
                else:
                    #a secondary sym - could be added if its primary is not present.
                    skipped_secondary.append(cols[gene_symbol_col])

            else:
                result['missing'].append(cols[gene_symbol_col])

        line_num += 1

    # Go through secondary genes and get their gene ids
    for gene in skipped_secondary:
        qry = """
        SELECT gene_id
        FROM gene_symbol
        WHERE label = %s
        """
        cursor.execute(qry, (gene,))

        # Use gene ids to get list of primary and secondary symbols
        for (gene_id,) in cursor:
            query = """
            SELECT label, is_primary
            FROM gene_symbol
            WHERE gene_id = %s
            """
            cursor.execute(query, (gene_id,))

            # Mark duplicate secondary gene symbols to be skipped
            for (label, is_primary) in cursor:
                if is_primary == 1:
                    if label in primary_genes:
                        # when primary symbol is found in dataset's gene list,
                        # skip the secondary gene
                        result['secondary'].append(gene)

    cursor.close()
    cnx.close()

    result['row_count'] = line_num

    print(json.dumps(result))

if __name__ == '__main__':
    main()
