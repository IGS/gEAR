#!/usr/bin/env python3

"""
Rescore dataset coloring for all 3 scopes: dataset, gene, and tissue

This script is a product of merging 3 recoloring scripts: 
    rescore_dataset_scope_dataset_coloring.py
    rescore_gene_scope_dataset_coloring.py
    rescore_tissue_scope_dataset_coloring.py


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

    ### BEGIN defining rescoring functions ###
    
    # Rescore Dataset Scope Coloring
    def rescore_dataset_scope_coloring(cursor):
        # if this is set to a value, debugging information will display for that gene
        debug_gene_id = None
        gene_scores = dict()
        
        dataset_min = None
        dataset_max = None
        dataset_vals = list()

        #for (eid, gene_id, clabel, rval) in cursor:
        for field in cursor:
            print(field)
            
            #order based on schema of create_schema.sql file
            eid = field[0]
            gene_id = field[2]
            clabel = field[3]
            rval = field[4]
            
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
            print("e = ", e)
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
    
    # Rescore Gene Scope Coloring
    def rescore_gene_scope_coloring(cursor):
        # if this is set to a value, debugging information will display for that gene
        debug_gene_id = None
        gene_scores = dict()

        #for (eid, gene_id, clabel, rval, gcidx) in cursor:
        for field in cursor:
            print(field)
            
            #order based on schema of create_schema.sql file
            eid = field[0]
            gene_id = field[2]
            clabel = field[3]
            rval = field[4]
            gcidx = field[11]
            
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
    
    # Rescore Tissue Scope Coloring
    def rescore_tissue_scope_coloring(cursor):
        # if this is set to a value, debugging information will display for that gene
        debug_gene_id = None
        gene_scores = dict()
        
        tissue = dict()

        #for (eid, gene_id, clabel, rval) in cursor:
        for field in cursor:
            print(field)
            
            #order based on schema of create_schema.sql file
            eid = field[0]
            gene_id = field[2]
            clabel = field[3]
            rval = field[4]
            
            
            if clabel not in tissue:
                tissue[clabel] = {'min': None, 'max': None, 'vals': list()}
                
            if tissue[clabel]['min'] is None or rval < tissue[clabel]['min']:
                tissue[clabel]['min'] = rval

            if tissue[clabel]['max'] is None or rval > tissue[clabel]['max']:
                tissue[clabel]['max'] = rval

            tissue[clabel]['vals'].append(rval)

            gene_scores[eid] = {'id': gene_id, 'tissue': clabel, 'rval': rval}

        for clabel in tissue:
            print("Tissue:{0}, min:{1}, max:{2}".format(clabel, tissue[clabel]['min'], tissue[clabel]['max']), file=sys.stderr)
            tissue[clabel]['vals'] = sorted(tissue[clabel]['vals'])

            cutoff_idx = int(len(tissue[clabel]['vals']) * .95)
            max_cutoff = tissue[clabel]['vals'][cutoff_idx]
            print("\t5% max cutoff adjusts to: {0}, idx: {1}/{2}".format(max_cutoff, cutoff_idx, len(tissue[clabel]['vals'])), file=sys.stderr)

            if max_cutoff < tissue[clabel]['max']:
                tissue[clabel]['max'] = max_cutoff

        for exp_id in gene_scores:
            e = gene_scores[exp_id]
            abs_step_size = tissue[ e['tissue'] ]['max'] / 255

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
                print("UPDATE expression SET tissue_based_abs_color = {0} WHERE id = {1};".format(abs_color_idx, exp_id))

    ### END of rescoring functions ###



    parser = argparse.ArgumentParser( description='Rescore the coloring for a given dataset')

    ## output file to be written
    parser.add_argument('-id', '--dataset_id', type=str, required=True, help='Dataset ID to rescore' )
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('gear.ini')

#    # if this is set to a value, debugging information will display for that gene
#    debug_gene_id = None

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

#    gene_scores = dict()
    
    ofh = open("{0}.update.undo.sql".format(args.dataset_id), 'wt')
    
    qry = """
        SELECT *
          FROM expression
         WHERE dataset_id = %s
    """
    cursor.execute(qry, (args.dataset_id, ))
    cursor.fetchone() #returns an object from query) http://www.tutorialspoint.com/python/python_database_access.htm
    print(cursor)
    
    # rescore dataset's coloring...
    rescore_dataset_scope_coloring(cursor)
    cnx.commit()
    
    rescore_gene_scope_coloring(cursor)
    cnx.commit()
    
    rescore_tissue_scope_coloring(cursor)

    cnx.commit()
    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()

  
