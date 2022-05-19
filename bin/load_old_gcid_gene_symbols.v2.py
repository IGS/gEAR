#!/opt/bin/python3
"""
Loads gene symbols from the old GCID "gene_symbols" table as new entries in the "gene" table on the new GCID server

sys.argv[1] = gene_mapping.txt where the fields are "ensembl_id" and "label"
"""

import csv, os, sys

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
from geardb import Connection
import mysql.connectors

gene_select_qry = """
    SELECT id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
            molecule, start, stop, gene_symbol, product, biotype
    FROM gene
    WHERE ensembl_id = %s
    ORDER BY gene_symbol, organism_id, ensembl_release DESC
"""

gene_insert_qry = """
    INSERT INTO gene (ensembl_id, ensembl_version, ensembl_release, genbank_acc,
                        organism_id, molecule, start, stop, gene_symbol, product,
                        biotype)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
"""

def main():

    conn = Connection()
    cursor = conn.get_cursor()

    with open(sys.argv[1]) as mapping_file:
        """
        tsv_reader = csv.DictReader(mapping_file, delimiter="\t")
        line_count = 0
        for row in tsv_reader:
            # Skip header
            if line_count == 0:
                line_count += 1
                continue

            # Get our elements
            ensembl_id = row["ensembl_id"]
            label = row["label"]
        """

        for line in mapping_file:
            line = line.rstrip()
            (ensembl_id, label) = line.split('\t')


            # Attempt to see if label is the "gene_symbol" in the "gene" table for this ensembl_id
            # If label does not exist as a gene_symbol, then insert a new entry
            cursor.execute(gene_select_qry, (ensembl_id,))

            label_found = False
            for (id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
                 molecule, start, stop, gene_symbol, product, biotype) in cursor:
                if label == gene_symbol:
                    label_found = True
                    break
            else:
                # If ensembl ID was not found, insert as very empty entry
                ensembl_version = None
                genbank_acc = None
                organism_id = 11
                molecule = None
                start = None
                stop = None
                product = None
                biotype = None
            if not label_found:
                print("-- inserting {}".format(label), file=sys.stderr)
                try:
                    cursor.execute(gene_insert_qry, (ensembl_id, ensembl_version, -1, genbank_acc, organism_id,
                        molecule, start, stop, label, product, biotype))
                    conn.commit()
                except mysql.connector.errors.DataError as e:
                    print("ERROR -- skipping {}".format(label), file=sys.stderr)

if __name__ == '__main__':
    main()
