#!/opt/bin/python3

"""

This script takes a dataset stored in our legacy relational database and converts 
it to a three tab set of files.  Previously, the column header format convention provided
some metadata for the columns.  These conventions and how they are transformed are:

Bar and line plots:
   Condition--Replicate_[pval|sd]
   Untreated--P_0
   Treated--P_1_sd
   Treated--P_1_pval

Violin plot
   Group--Label
   OHC--Cell1
   OHC--Cell2

SVG
   Label

| ENSG00000271579 | AC078880.3     |   0.035522 |   0.0268332 | 8699        | 2h        |
| ENSG00000271579 | AC078880.3     |  0.0585463 |   0.0608098 | 8699        | 6h        |
| ENSG00000271579 | AC078880.3     |  0.0337311 |   0.0309742 | HBEC_strep  | 6h        |
| ENSG00000271579 | AC078880.3     |  0.0504333 |   0.0594605 | K79         | 2h        |
| ENSG00000271579 | AC078880.3     |  0.0267187 |   0.0158348 | K79         | 6h        |
| ENSG00000276308 | AC078880.4     |          0 |           0 | 8699        | 2h        |
| ENSG00000276308 | AC078880.4     |          0 |           0 | 8699        | 6h        |
| ENSG00000276308 | AC078880.4     |          0 |           0 | HBEC_strep  | 6h        |
| ENSG00000276308 | AC078880.4     |          0 |           0 | K79         | 2h        |
| ENSG00000276308 | AC078880.4     |          0 |           0 | K79         | 6h        |

These three files are created in the output directory:

  expression.tab
  genes.tab
  observations.tab


"""

import argparse
import os
import re
import sys

import mysql.connector

sys.path.append("{0}/../lib".format(os.path.dirname(sys.argv[0])))
import geardb

def main():
    parser = argparse.ArgumentParser( description='Annotation extractor -> H5AD')
    parser.add_argument('-i', '--input_dataset_id', type=str, required=True, help='ID of dataset to read from db')
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='Path to the directory where three output files will be created')
    
    args = parser.parse_args()

    cnx = geardb.Connection()
    cursor = geardb.Connection().get_cursor()

    genes_fh = open("{0}/genes.tab".format(args.output_directory), 'wt')
    genes_fh.write("Ensembl_ID\tgene_symbol\n")
    
    expression_fh = open("{0}/expression.tab".format(args.output_directory), 'wt')
    obs_fh = open("{0}/observations.tab".format(args.output_directory), 'wt')

    qry = """
          SELECT g.ensembl_id, gs.label AS gene_symbol, e.raw_value, 
                 e.std_dev, e.group_name, e.sec_label, e.class_label
            FROM expression e
                 JOIN dataset d ON d.id=e.dataset_id
                 JOIN gene g ON e.gene_id = g.id
                 JOIN gene_symbol gs ON gs.gene_id=g.id
           WHERE d.id = %s
                 AND gs.is_primary = 1
    """
    cursor.execute(qry, [args.input_dataset_id, ])

    id_to_symbol_idx = dict()
    class_labels = list()
    expr = dict()
    obs = dict()
    obs_order = list()

    has_cluster = False
    has_condition = False
    
    for (ensembl_id, gene_symbol, raw_val, std_dev, group, condition, class_label) in cursor:
        if ensembl_id not in id_to_symbol_idx:
            id_to_symbol_idx[ensembl_id] = gene_symbol
            genes_fh.write("{0}\t{1}\n".format(ensembl_id, gene_symbol))

        if class_label not in class_labels:
            class_labels.append(class_label)

        if ensembl_id not in expr:
            expr[ensembl_id] = dict()

        if group:
            has_cluster = True

        if condition:
            has_condition = True

        if class_label in expr[ensembl_id]:
            raise Exception("Got duplicate class_label ({0}) for ensembl_id ({1})".format(class_label, ensembl_id))
        else:
            # e=expression, g=group, c=condition
            expr[ensembl_id][class_label] = raw_val
            
            if class_label in obs:
                # make sure the g and c are the same, else drama
                if obs[class_label]['g'] != group or  obs[class_label]['c'] != condition:
                    raise Exception("Found duplicate class name ({0}) but with different group/conditions parsed".format(class_label))
            else:
                obs[class_label] = {'g':group, 'c':condition}
                obs_order.append(class_label)

    obs_fh.write("observations")
    if has_cluster: obs_fh.write("\tcluster")
    if has_condition: obs_fh.write("\tcondition")
    obs_fh.write("\n")

    for class_label in obs_order:
        obs_fh.write(class_label)
        if has_cluster: obs_fh.write("\t{0}".format(obs[class_label]['g']))
        if has_condition: obs_fh.write("\t{0}".format(obs[class_label]['c']))
        obs_fh.write("\n")

    # Expression header row
    expression_fh.write("Ensembl_ID")

    for class_label in obs_order:
        expression_fh.write("\t{0}".format(class_label))

    expression_fh.write("\n")
        
    # Expression lines
    for ensembl_id in expr:
        expression_fh.write("{0}".format(ensembl_id))

        for class_label in obs_order:
            expression_fh.write("\t{0}".format(expr[ensembl_id][class_label]))

        expression_fh.write("\n")
    
    cursor.close()

    genes_fh.close()
    expression_fh.close()
    obs_fh.close()
            
if __name__ == '__main__':
    main()







