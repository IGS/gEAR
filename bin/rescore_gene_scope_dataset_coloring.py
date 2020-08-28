#!/usr/bin/env python3

"""

ALTER TABLE expression ADD gene_based_abs_color INT NOT NULL DEFAULT 0 AFTER gene_based_color;
ALTER TABLE expression ADD tissue_based_abs_color INT NOT NULL DEFAULT 0 AFTER tissue_based_color;
ALTER TABLE expression ADD dataset_based_abs_color INT NOT NULL DEFAULT 0 AFTER dataset_based_color;

./bin/rescore_gene_scope_dataset_coloring.py -id c69485b2-6f8d-c60e-7337-e7ebad89b2c0 > c69485b2-6f8d-c60e-7337-e7ebad89b2c0.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id 154b90c5-f427-f7ff-b63e-13b3a82f946a > 154b90c5-f427-f7ff-b63e-13b3a82f946a.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id 27c500eb-3879-2602-c4ee-f8ed769fe7ca > 27c500eb-3879-2602-c4ee-f8ed769fe7ca.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id 5f9b633b-3520-97b0-1fa9-13b5f66df52b > 5f9b633b-3520-97b0-1fa9-13b5f66df52b.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id 6fdd350c-4f82-07e2-3a39-408f105db16d > 6fdd350c-4f82-07e2-3a39-408f105db16d.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id 70fc21b1-5eaf-3594-36ab-f7fce58b587f > 70fc21b1-5eaf-3594-36ab-f7fce58b587f.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id 7812a487-932b-32f7-2de7-33dd3155c849 > 7812a487-932b-32f7-2de7-33dd3155c849.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id cf8272cb-57fa-e841-0b50-9198e62fe2ff > cf8272cb-57fa-e841-0b50-9198e62fe2ff.update.sql
./bin/rescore_gene_scope_dataset_coloring.py -id e0fa50cc-ddcd-47c4-0df6-55935d8260c4 > e0fa50cc-ddcd-47c4-0df6-55935d8260c4.update.sql

cat *.update.sql | grep _abs > foo
mv foo mass.20160325.update.abs.sql


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

    ofh = open("{0}.update.undo.sql".format(args.dataset_id), 'wt')
    
    qry = """
        SELECT id, gene_id, class_label, raw_value, gene_based_color
          FROM expression
         WHERE dataset_id = %s
    """
    cursor.execute(qry, (args.dataset_id, ))

    for (eid, gene_id, clabel, rval, gcidx) in cursor:
        if gene_id not in gene_scores:
            gene_scores[gene_id] = list()

        gene_scores[gene_id].append({'id':eid, 'val':rval, 'gcidx':gcidx})

    for gene_id in gene_scores:
        gene_min = None
        gene_max = None
        score_sum = 0
        score_count = 0
        score_mean = 0

        ## this needs to be handled
        if len(gene_scores[gene_id]) == 1:
            continue

        for score in gene_scores[gene_id]:
            if gene_min is None or score['val'] < gene_min:
                gene_min = score['val']

            if gene_max is None or score['val'] > gene_max:
                gene_max = score['val']

            score_sum += score['val']
            score_count += 1

        score_mean = score_sum / score_count
        max_diff_from_mean = gene_max - score_mean
        if (score_mean - gene_min) > max_diff_from_mean: max_diff_from_mean = score_mean - gene_min
        step_size = (max_diff_from_mean * 2) / 255
        #abs_step_size = (gene_max - gene_min) / 255
        abs_step_size = gene_max / 255

        if step_size == 0:
            step_size = 0.01

        if abs_step_size == 0:
            abs_step_size = 0.01

        if gene_min == 1:
            #print("WARN: min was 0 for gene:{0}, so had to adjust it to 1".format(gene_id))
            gene_min = 1
                
        gene_span = gene_max - gene_min
        span_adj = math.log(gene_max / gene_min, 2)
        # The adjustment can't be more than 10
        if span_adj > 10: span_adj = 10

        span_adj_perc = span_adj / 10

        if debug_gene_id is not None and gene_id == debug_gene_id:
            print("Debugging information for gene id: {0}".format(gene_id))
            print("-----------------------------------------------")
            print("Gene score mean: {0}".format(score_mean))
            print("Max diff from mean: {0}".format(max_diff_from_mean))
            print("Step size: {0}".format(step_size))
            print("Abs step size: {0}".format(abs_step_size))
            print("Gene span: {0}".format(gene_span))
            print("Span adjustment: {0}".format(span_adj))
            print("Span adjustment %: {0}".format(span_adj_perc))
            
        for score in gene_scores[gene_id]:
            steps = (score['val'] - score_mean) / step_size
            steps = steps * span_adj_perc
            color_idx = 127 + steps

            abs_steps = score['val'] / abs_step_size
            abs_color_idx = abs_steps

            if color_idx > 254:
                color_idx = 254
            elif color_idx < 0:
                color_idx = 0

            if abs_color_idx > 254:
                abs_color_idx = 254
            elif abs_color_idx < 0:
                abs_color_idx = 0

            score['gcidx_new'] = color_idx

            #print("Gene:{0}, min:{1}, max:{2}, span:{3}, log2(span_adj):{4:.1f}, raw:{6}, color_idx:{5}".format(
            #    gene_id, gene_min, gene_max, gene_span, span_adj, score['gcidx_new'], score['val']))
            if debug_gene_id is None or debug_gene_id == gene_id:
                ofh.write("UPDATE expression SET gene_based_color = {0} WHERE id = {1};\n".format(score['gcidx'], score['id']))
                print("UPDATE expression SET gene_based_color = {0} WHERE id = {1};".format(score['gcidx_new'], score['id']))
                print("UPDATE expression SET gene_based_abs_color = {0} WHERE id = {1};".format(abs_color_idx, score['id']))
            

    cnx.commit()
    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()


# before adjustment:
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:2198.75, color_idx:224
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:55.25, color_idx:3
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:2489.5, color_idx:255
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:20.0, color_idx:0 

# after adjustment:
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:2198.75, color_idx:156
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:55.25, color_idx:2
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:2489.5, color_idx:177
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:20.0, color_idx:0

# after adjustment + blue shift fix
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:2198.75, color_idx:195.0
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:55.25, color_idx:41.0
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:2489.5, color_idx:216.0
#Gene:637, min:20.0, max:2489.5, span:2469.5, log2(span_adj):7.0, raw:20.0, color_idx:39.0

#0-------------127--------------255

#before fix
#0-------------------177

#after fix Add 255-177 to all
#     0-------------------177
