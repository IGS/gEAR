#!/usr/bin/env python3

"""
Collects site statistics and outputs results as JSON file into www/stats_archive/

Statistics collected:
-- Number of Registered users
-- Number of datasets
-- Number of expression points
-- Number of searches (SELECT queries) since script was last run
    -NOTES:
    -- This query returns the number of SELECT queries performed:
            SHOW GLOBAL STATUS WHERE variable_name = 'Com_select';
            +---------------+-------+
            | Variable_name | Value |
            +---------------+-------+
            | Com_select    | 115   |
            +---------------+-------+
            1 row in set (0.00 sec)

        The issue is this is ALL SELECT queries, so the number is inflated as many UI actions involve multiple MySQL queries

"""

import gzip
import json
import os
import re
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

import gear.db

stats_file_symlink = os.path.abspath('../www/stats_archive/site_stats.json')

tables_to_count = ['guser', 'dataset', 'expression'] #MySQL tables to count


def main():
    print('Starting to gather gEAR site stats...')
    #output_file
    stats_file = os.path.abspath('../www/stats_archive/site_stats.' + datestamp() + '.json')
    log_base_dir = '/var/log/apache2'

    result = {}

    ## search stats
    result['fp_gene_search_count'] = get_fp_gene_searches(log_base_dir)

    #connect to gEAR MySQL
    mysql_cnx = gear.db.MySQLDB().connect()
    cursor = mysql_cnx.cursor()

    #get counts for each MySQL table
    for table in tables_to_count:
        if table == 'expression':
            h5ad_count = get_h5ad_count()
            raw_count += h5ad_count
        else:
            raw_count = get_count(cursor=cursor, column='id', table=table)
            
        result[table + '_count'] = raw_count

    # Get count of select commands performed and format it
    result['search_count'] = get_sql_select_count(cursor)

    # Print results to file
    with open(stats_file, 'w') as outfile:
        json.dump(result, outfile)

    try:
        #try to create symlink
        os.symlink(stats_file, stats_file_symlink)
    except FileExistsError:
        # symlink exists, remove it and create an updated one
        os.remove(stats_file_symlink)
        os.symlink(stats_file, stats_file_symlink)

    print('Completed.\n')
    print(json.dumps(result, indent=4, sort_keys=True))
    cursor.close()


def get_count(cursor=None, column=None, table=None):
    # Returns the row count of the table
    print('Getting count for MySQL table: {0}'.format(table))

    qry = """SELECT COUNT({0})
        FROM {1}
    """.format(column, table)
    cursor.execute(qry)

    count = None
    for row in cursor:
        count = row[0]
        break

    print('...Done')
    return count

def get_fp_gene_searches(base_dir):
    log_file_base_name = 'ssl_umgear_access.log'

    file_is_compressed = False
    gene_search_count = 0

    for fname in os.listdir(base_dir):
        if not fname.startswith(log_file_base_name):
            continue

        fpath = "{0}/{1}".format(base_dir, fname)

        if fname.endswith('.gz'):
            fh = gzip.open(fpath, 'rb')
            file_is_compressed = True
        else:
            fh = open(fpath, 'rU')
            file_is_compressed = False

        for line in fh:
            if file_is_compressed:
                line = line.decode()

            m = re.search('search_genes.py', line)
            if m:
                gene_search_count += 1

    return gene_search_count

def get_sql_select_count(cursor=None):
    #Returns count of select commands from MySQL global status
    print('Counting MySQL searches')

    prior_search_count = 0

    #get the last recorded count
    if os.path.isfile(stats_file_symlink):
        with open(stats_file_symlink, 'r') as f:
            site_stats = json.load(f)

        if site_stats['search_count']:
            prior_search_count = site_stats['search_count']

    qry = """SHOW GLOBAL STATUS WHERE variable_name = 'Com_select'"""
    cursor.execute(qry)

    current_count = None
    for row in cursor:
        current_count = int(row[1])
        break

    #The difference is the new count
    new_count = current_count - prior_search_count

    #Remove the queries from this script
    new_count -= len(tables_to_count)

    print('...Done')
    return new_count


def get_h5ad_count():
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.sparse import csr_matrix

    h5ad_dir = '../www/datasets/'
    h5ad_expression_count = 0
    print('Counting expression points - h5ad datasets')

    for h5_file in os.listdir(h5ad_dir):
        if h5_file.endswith('h5ad'):
            sc.settings.verbosity = 0
            adata = sc.read_h5ad(os.path.join(h5ad_dir, h5_file))
            print('... counting {0}'.format(h5_file))

            #Get only the observed data
            obs_matrix = adata.X.copy()

            # I feel this is much more accurate than prior methods
            # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.csr_matrix.nnz.html
            # "Get the count of explicitly-stored values (nonzeros)"
            sp_arr = csr_matrix(obs_matrix)
            this_count = sp_arr.nnz
            h5ad_expression_count += sp_arr.nnz

            print("\tAdding {0} expression points for a total of: {1}".format(this_count, h5ad_expression_count))
            
    print('...Done')
    return h5ad_expression_count


def datestamp():
    # Returns a date-time stamp
    import datetime
    return datetime.datetime.now().strftime('%Y%m%d') #Example: 20180207




if __name__ == '__main__':
    main()
