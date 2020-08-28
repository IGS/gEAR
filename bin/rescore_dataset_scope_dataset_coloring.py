#!/usr/bin/env python3

"""

ALTER TABLE expression ADD gene_based_abs_color INT NOT NULL DEFAULT 0 AFTER gene_based_color;
ALTER TABLE expression ADD tissue_based_abs_color INT NOT NULL DEFAULT 0 AFTER tissue_based_color;
ALTER TABLE expression ADD dataset_based_abs_color INT NOT NULL DEFAULT 0 AFTER dataset_based_color;

./bin/rescore_dataset_scope_dataset_coloring.py -id c69485b2-6f8d-c60e-7337-e7ebad89b2c0 > c69485b2-6f8d-c60e-7337-e7ebad89b2c0.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id 154b90c5-f427-f7ff-b63e-13b3a82f946a > 154b90c5-f427-f7ff-b63e-13b3a82f946a.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id 27c500eb-3879-2602-c4ee-f8ed769fe7ca > 27c500eb-3879-2602-c4ee-f8ed769fe7ca.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id 5f9b633b-3520-97b0-1fa9-13b5f66df52b > 5f9b633b-3520-97b0-1fa9-13b5f66df52b.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id 6fdd350c-4f82-07e2-3a39-408f105db16d > 6fdd350c-4f82-07e2-3a39-408f105db16d.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id 70fc21b1-5eaf-3594-36ab-f7fce58b587f > 70fc21b1-5eaf-3594-36ab-f7fce58b587f.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id 7812a487-932b-32f7-2de7-33dd3155c849 > 7812a487-932b-32f7-2de7-33dd3155c849.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id cf8272cb-57fa-e841-0b50-9198e62fe2ff > cf8272cb-57fa-e841-0b50-9198e62fe2ff.dscope.update.sql
./bin/rescore_dataset_scope_dataset_coloring.py -id e0fa50cc-ddcd-47c4-0df6-55935d8260c4 > e0fa50cc-ddcd-47c4-0df6-55935d8260c4.dscope.update.sql

cat *dscope.update.sql > foo
mv foo mass.20160325.update.abs.dataset.sql


"""

import argparse
import mysql.connector
import configparser
import os
import re
import statistics
import sys
import math

def main():
    parser = argparse.ArgumentParser( description='Rescore the coloring for a given dataset')

    ## output file to be written
    parser.add_argument('-id', '--dataset_id', type=str, required=True, help='Dataset ID to rescore' )
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('gear.ini')

    # if this is set to a value, debugging information will display for that gene
    debug_gene_id = None

    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])

    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)

    cursor = cnx.cursor()

    gene_scores = dict()

    qry = """
        SELECT id, gene_id, class_label, raw_value
          FROM expression
         WHERE dataset_id = %s
    """
    cursor.execute(qry, (args.dataset_id, ))

    dataset_min = None
    dataset_max = None
    dataset_vals = list()

    for (eid, gene_id, clabel, rval) in cursor:
        if dataset_min is None or rval < dataset_min:
            dataset_min = rval

        if dataset_max is None or rval > dataset_max:
            dataset_max = rval

        dataset_vals.append(rval)

        gene_scores[eid] = {'id': gene_id, 'tissue': clabel, 'rval': rval}

    dataset_vals = sorted(dataset_vals)
    cutoff_idx = int(len(dataset_vals) * .95)
    print("Dataset max was {0} but will be adjusted to {1} at 5% cutoff".format(dataset_max, dataset_vals[cutoff_idx]), file=sys.stderr)
    max_cutoff = dataset_vals[cutoff_idx]

    for exp_id in gene_scores:
        e = gene_scores[exp_id]
        abs_step_size = max_cutoff / 255

        if abs_step_size == 0:
            abs_step_size = 0.01

        if debug_gene_id is not None and gene_id == debug_gene_id:
            print("Debugging information for gene id: {0}".format(gene_id))
            print("-----------------------------------------------")
            print("Abs step size: {0}".format(abs_step_size))

        abs_steps = e['rval'] / abs_step_size
        abs_color_idx = abs_steps

        if abs_color_idx > 254:
            abs_color_idx = 254
        elif abs_color_idx < 0:
            abs_color_idx = 0

        if debug_gene_id is None or debug_gene_id == gene_id:
            print("UPDATE expression SET dataset_based_abs_color = {0} WHERE id = {1};".format(abs_color_idx, exp_id))

    cnx.commit()
    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()



