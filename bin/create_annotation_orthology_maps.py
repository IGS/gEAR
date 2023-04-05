#!/usr/bin/env python3

"""

Transforms a directory of orthology and ID mapping files into a data structure (h5) we can 
more rapidly access.

Expected files:

# Get the orthology itself
#   AGR has aggregated these:
#   PhylomeDB, OrthoFinder, ZFIN, Hieranoid, OMA, Ensembl Compara, Roundup, InParanoid, PANTHER, TreeFam, and OrthoInspector.
wget https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz

# Get the ID mappings
# Taxon IDs: NCBITaxon:9606, NCBITaxon:10116, NCBITaxon:10090, NCBITaxon:7955
# Species: Homo sapiens, Rattus norvegicus, Mus musculus, Danio rerio
wget https://www.alliancegenome.org/api/geneMap/ensembl?species=9606 -O id.map.9606.tab
wget https://www.alliancegenome.org/api/geneMap/ensembl?species=10116 -O id.map.10116.tab
wget https://www.alliancegenome.org/api/geneMap/ensembl?species=10090 -O id.map.10090.tab
wget https://www.alliancegenome.org/api/geneMap/ensembl?species=7955 -O id.map.7955.tab

Output file convention:

orthomap.{qry_organism_id}.{qry_source}.{qry_release}__{tgt_organism_id}.{tgt_source}.{tgt_release}.hdf5

Example:

orthomap.1.ensembl.93__2.ensembl.95.hdf5

The organism_id is from the gEAR MySQL table organism_id

Output format:

Dataframe where:

  - Source identifier (index)
  - Source gene symbol
  - Target identifier
  - Target gene symbol

"""

import argparse
import os
import pandas as pd
import sys


def main():
    parser = argparse.ArgumentParser( description='Create H5 dataframe from AGR orthology files')
    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to directory with input files' )
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='Path to an output directory base where files will be created' )
    args = parser.parse_args()

    # Only those taxon IDs here will be included
    # taxon_id -> organism_id
    taxon_ids = {7955: 3,
                 9606: 2,
                 10090: 1,
                 10116: 6,
                 6239: 8
    }

    # load each of the annotations
    refannot = dict()
    for taxon_id in taxon_ids:
        annotation_file_path = "{0}/id.map.{1}.tab".format(args.input_directory, taxon_id)
        refannot[taxon_id] = load_annotation(annotation_file_path)
            
    orthofile_path = "{0}/ORTHOLOGY-ALLIANCE_COMBINED.tsv".format(args.input_directory)
    orthomap = load_orthomap(orthofile_path)

    #orthomap.10090.ensembl.93__7955.ensembl.95.hdf5

    for tx1 in taxon_ids:
        for tx2 in taxon_ids:
            if tx1 == tx2:
                continue

            identifiers = list()
            annot = {'gs1': [], 'id2': [], 'gs2': []}

            print("\t\tProcessing orthomap between {0} and {1}".format(tx1, tx2), file=sys.stderr)

            for g1_id in orthomap[tx1][tx2]:
                # For some reason there are identifiers in the ortholog map which aren't found in
                #  their own annotation references, such as ZFIN:ZDB-GENE-041111-196 (as of 2022-09-12
                #  We can only skip them.
                if g1_id in refannot[tx1]:
                    g2_id = orthomap[tx1][tx2][g1_id]['winning_id']
                    
                    if g2_id in refannot[tx2]:
                        identifiers.append(refannot[tx1][g1_id]['ensembl_id'])
                        annot['gs1'].append(refannot[tx1][g1_id]['gene_symbol'])
                        annot['id2'].append(refannot[tx2][g2_id]['ensembl_id'])
                        annot['gs2'].append(refannot[tx2][g2_id]['gene_symbol'])
                    else:
                        print("WARN: feature ({0}) found in ortholog mapping but has no entry in taxon ({1}) reference file".format(g2_id, tx2))
                        continue
                else:
                    print("WARN: feature ({0}) found in ortholog mapping but has no entry in taxon ({1}) reference file".format(g1_id, tx1))
                    continue

            h5_file_path = os.path.abspath("{0}/orthomap.{1}.ensembl__{2}.ensembl.hdf5".format(
                args.output_directory, taxon_ids[tx1], taxon_ids[tx2]))
            df = pd.DataFrame(annot, index=identifiers)
            df.to_hdf(os.path.basename(h5_file_path), os.path.dirname(h5_file_path))
    
def load_annotation(fpath):
    """
    There can be duplicates in the annotation. For example:

    NCBITaxon:7955	ENSDARG00000009582	ZFIN:ZDB-GENE-050517-19	abcc6b.1	
    NCBITaxon:7955	ENSDARG00000105403	ZFIN:ZDB-GENE-050517-19	abcc6b.1

    Both gene identifiers have the  same gene symbol but map to different Ensembl IDs. Each
    instance of this I checked showed the older Ensembl ID to have been removed
    from the Ensembl database and replaced by the newer one.

    Found a few instances of weirdness. For example, in their annotation files, one 
    rat identifier is annotated with both mouse and rat ensembl id:

    NCBITaxon:10116	ENSMUSG00000018326	RGD:61998	Ywhab	56011
    NCBITaxon:10116	ENSRNOG00000010945	RGD:61998	Ywhab	56011

    They also chose to use commas to separate multiple entries vs. additional rows. Ex:

    NCBITaxon:9606	ENSG00000266967	HGNC:43946,HGNC:28417	PTGES3L-AARSD1,AARSD1	100885850,80755
    NCBITaxon:9606	ENSG00000108825	HGNC:43946,HGNC:28417	PTGES3L-AARSD1,AARSD1	100885850,80755

    Given that we're providing mapping like HGNC:43946 -> ENSG00000266967 we just do this for both keys
    """
    print("Loading annotation file: {0}".format(fpath))
    annotation = dict()
    line_num = 0

    for line in open(fpath):
        line_num += 1

        if line_num == 1:
            continue
        
        line = line.rstrip()
        cols = line.split("\t")
        feat_ids = cols[2].split(',')

        for feat_id in feat_ids:
            if feat_id in annotation:
                # Take the later annotation ID, numerical string comparison works here
                if cols[1] > annotation[feat_id]['ensembl_id']:
                    print("Found duplicate feat ID: {0} -- {1}, overwriting previous: {2}".format(
                        feat_id, cols[1], annotation[feat_id]['ensembl_id']), file=sys.stderr)

                    annotation[feat_id] = {'ensembl_id': cols[1], 'gene_symbol': cols[3]}
            else:
                annotation[feat_id] = {'ensembl_id': cols[1], 'gene_symbol': cols[3]}

    return annotation

def load_orthomap(fpath):
    """
    The Genome Alliance orthology map can contain multiple matches for any query, even
    against the same target genome. For each pairing this currently filters out the one
    with the most tools matching. In the case of a tie, the first gene is chosen.
    """
    omap = dict()
    line_num = 0
    for line in open(fpath):
        if line[0] == '#':
            continue

        line_num += 1

        if line_num == 1:
            continue
        
        line = line.rstrip()
        cols = line.split("\t")

        gene1_id  = cols[0]
        gene1_tax = int(cols[2].split(":")[1])
        gene2_id  = cols[4]
        gene2_tax = int(cols[6].split(":")[1])
        algos_match = int(cols[9])

        if gene1_tax not in omap:
            omap[gene1_tax] = dict()

        if gene2_tax not in omap[gene1_tax]:
            omap[gene1_tax][gene2_tax] = dict()

        if gene1_id not in omap:
            omap[gene1_tax][gene2_tax][gene1_id] = {'match_count': 0, 'winning_id': None}

        if algos_match > omap[gene1_tax][gene2_tax][gene1_id]['match_count']:
            omap[gene1_tax][gene2_tax][gene1_id] = {'match_count': algos_match, 'winning_id': gene2_id}

    return omap
    
        
if __name__ == '__main__':
    main()







