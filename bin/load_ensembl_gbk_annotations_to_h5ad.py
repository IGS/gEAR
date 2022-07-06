#!/opt/bin/python3

"""

This is a unified, all-purpose annotation loader for the gEAR database. It converts the 
GBK-formatted ENSEMBL releases to an H5AD file.  

TODO:

Add a mode option to have it run in two modes:

   - 'full' (default): Stores all attributes of the annotation
   - 'minimal': Stores only the identifier and gene symbol of each gene

In either case, the identifier is the index.

Currently, it exports the 'full' version.

"""

import argparse
import gzip
import os
import re
import sys

import pandas as pd
from biocode import annotation

sys.path.append("{0}/../lib".format(os.path.dirname(sys.argv[0])))
import geardb

from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser( description='Annotation loader for the gEAR')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Path to an input directory to be read')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file to be created')
    args = parser.parse_args()

    files = os.listdir(args.input_dir)

    if len(files) == 0:
        raise Exception("ERROR: No input files matching the expected naming convention found")

    skipped_biotypes = ['misc_RNA', 'pseudogene', 'ribozyme']

    identifiers = list()
    annot = {
        'gene_symbol': [],
        'biotype': [],
        'product': [],
        'version': [],
        'chromosome': [],
        'fmin': [],
        'fmax': [],
        # these are lists
        'go_ids': [],
        'dbxrefs': []
    }

    for file in sorted(files):
        if file.endswith('.dat.gz'):
            genes_inserted = 0
            
            fpath = "{0}/{1}".format(args.input_dir, file)
            print("Processing file: {0} ... ".format(fpath), file=sys.stderr, end='')

            # Each of the file names across organisms so far has something like '.chromosome.X.' in
            #  the name.  Let's use this to get a molecule name
            m = re.search("\.(chromosome|nonchromosomal)\.(\w+)", file)
            if m.group(1) == "chromosome":
                chr_name = m.group(2)
            elif m.group(1) == "nonchromosomal":
                chr_name = "nonchromosomal"
            else:
                raise Exception("Error: failed to pull chromosome name from file ({0}).  Expected to contain '.chromosome.X'".format(file))

            with gzip.open(fpath, 'rt') as fh:
                for gb_record in SeqIO.parse(fh, "genbank"):
                    annotations = process_gb_record(gb_record)

                    for ann in annotations:
                        # make sure the product name isn't too long
                        if ann.product_name is not None and len(ann.product_name) > 255:
                            ann.product_name = ann.product_name[:255]

                        # make sure the gene product isn't too long
                        if ann.gene_symbol is not None and len(ann.gene_symbol) > 20:
                            ann.gene_symbol = ann.gene_symbol[:20]

                        identifiers.append(ann.other_attributes['ensembl_id'])
                        annot['version'].append(ann.other_attributes['ensembl_version'])
                        annot['chromosome'].append(chr_name)
                        annot['fmin'].append(ann.other_attributes['loc']['fmin'])
                        annot['fmax'].append(ann.other_attributes['loc']['fmax'])
                        annot['gene_symbol'].append(ann.gene_symbol)
                        annot['product'].append(ann.product_name)
                        annot['biotype'].append(ann.other_attributes['biotype'])
                            
                        genes_inserted += 1

                        # The Ensembl records have many, many verified duplications in the annotation.  Keep this list
                        #  so we only load one of each.
                        go_ids = set()
                        for go_annot in ann.go_annotations:
                            go_ids.add(go_annot.go_id)
                        annot['go_ids'].append(list(go_ids))
                        

                        dbxrefs = set()
                        for dbxref in ann.dbxrefs:
                            label = "{0}:{1}".format(dbxref.db, dbxref.identifier)
                            dbxrefs.add(label)
                        annot['dbxrefs'].append(list(dbxrefs))

            print("{0} genes added".format(genes_inserted), file=sys.stderr)

    df = pd.DataFrame(annot, index=identifiers)
    df.to_hdf(os.path.basename(args.output_file), os.path.dirname(args.output_file))

def get_biotype(feat):
    """
    For the misc_RNA feature types the actual subtype is stored as the first /note value.  Example:
     misc_RNA        complement(2623651..2623758)
                     /gene="ENSGALG00000026164.1"
                     /db_xref="RefSeq_ncRNA:NR_105442"
                     /note="miRNA"
                     /note="transcript_id=ENSGALT00000044715.1"
    """
    biotype = None

    if 'note' in feat.qualifiers:
        biotype = feat.qualifiers['note'][0]

    return biotype

    
def get_coordinates(feat):
    fmin = int(feat.location.start)
    fmax = int(feat.location.end)

    if feat.location.strand == 1:
        strand = '+'
    elif feat.location.strand == -1:
        strand = '-'
    else:
        strand = None

    return {'fmin': fmin, 'fmax': fmax, 'strand': strand}

def get_dbxrefs(feat):
    """
    Looks for sections like this to extract DB xrefs.  Here are some unique ones within a single CDS entry:

     CDS             join(91298396..91298448,91300528..91300640, ...
                     /gene="ENSMUSG00000026307.12"
                     /db_xref="CCDS:CCDS15159"
                     /db_xref="RefSeq_peptide:NP_057926"
                     /db_xref="RefSeq_mRNA:NM_016717"
                     /db_xref="RefSeq_mRNA_predicted:XM_006529743"
                     /db_xref="RefSeq_ncRNA_predicted:XR_387181"
                     /db_xref="RefSeq_peptide_predicted:XP_006529806"
                     /db_xref="Uniprot/SPTREMBL:A0A0R4J069"
                     /db_xref="UCSC:uc007caf.1"
                     /db_xref="EMBL:AC102874"
                     /db_xref="GO:0003824"
                     /db_xref="MGI_trans_name:Scly-201"
                     /db_xref="Vega_translation:181101"
                     /db_xref="KEGG_Enzyme:00730+2.8.1.7"
                     /db_xref="Reactome:R-MMU-1430728"
                     /db_xref="UniParc:UPI0000022876"

    Duplicates from the same source are expected and allowed in the schema.  So, two references to EMBL is
    OK, but not if both have the same value.  That sort of duplicate is simply skipped.

    Sources skipped because they are handled elsewhere:  GO

    """
    dbxref_sources_to_skip = ['GO']
    dbxrefs = list()

    if 'db_xref' in feat.qualifiers:
        for dbxref in feat.qualifiers['db_xref']:
            m = re.match('(.+):(.+)', dbxref)
            if m:
                if m.group(1) not in dbxref_sources_to_skip:
                    dbxrefs.append(annotation.Dbxref(db=m.group(1), identifier=m.group(2)))

    return dbxrefs

def get_genbank_acc(feat):
    """
    Looks for sections like this:

         CDS             742516..743481
                     /gene="ENSGALG00000045828.1"
                     /db_xref="RefSeq_mRNA_predicted:XM_425093"

    And extracts out the RefSeq_mRNA_predicted accession
    """
    if 'db_xref' in feat.qualifiers:
        for dbxref in feat.qualifiers['db_xref']:
            m = re.match('RefSeq_mRNA_predicted:(.+)', dbxref)
            if m:
                return m.group(1)

    return None

def get_go_annotations(feat):
    """
    Looks for sections like this to extract GO terms

     CDS             join(3366667..3366969,3463389..3463463)
                     /gene="ENSMUSG00000040653.6"
                     /protein_id="ENSMUSP00000149688.1"
                     /db_xref="EMBL:AC162384"
                     /db_xref="EMBL:CH466562"
                     /db_xref="GO:0005737"
                     /db_xref="GO:0042325"
    """
    go_terms = list()

    if 'db_xref' in feat.qualifiers:
        for dbxref in feat.qualifiers['db_xref']:
            m = re.match('GO:(.+)', dbxref)
            if m:
                go_terms.append(annotation.GOAnnotation(go_id=m.group(1)))

    return go_terms

def get_go_index(curs):
    """
    Creates an index of exiting GO terms in the database, where the key is the 
    GO id (like 'GO:013454') and the value is simply True for all entries, allowing
    for an O(1) hash lookup
    """
    qry = "SELECT go_id FROM go"
    curs.execute(qry)

    idx = dict()
    for row in curs:
        idx[row[0]] = True

    return idx
                

def get_gene_sym(feat):
    if 'locus_tag' in feat.qualifiers:
        for locus_tag in feat.qualifiers['locus_tag']:
            # take off anything after the space, if present
            return locus_tag.split(' ')[0]

    return None

def get_product(feat):
    product = None

    if 'note' in feat.qualifiers:
        product = feat.qualifiers['note'][0]

        # skip the products like this:
        #   /note="transcript_id=ENSGALT00000058308.1"
        if product.startswith('transcript_id='):
            return None

        # process the mouse products, which are like: limb and neural patterns [Source:MGI Symbol;Acc:MGI:1918115]
        # taking off everything starting with [Source:
        m = re.match('(.+) \[Source', product)
        if m:
            product = m.group(1)

    return product

def process_gb_record(record):
    mol_id = record.name
    annotations = list()
    current_annot = None

    for feat in record.features:
        if feat.type == 'source': continue
        
        # Only do the features associated with a gene
        if 'gene' not in feat.qualifiers:
            continue

        if feat.type == 'gene':
            # append the gene we've been working on and reset.
            if current_annot is not None:
                annotations.append(current_annot)
            
            current_annot = process_gb_gene(feat)

        elif feat.type == 'mRNA':
            current_annot.other_attributes['biotype'] = 'mRNA'
        
        elif feat.type == 'misc_RNA':
            # here expect lincRNA, miRNA, snoRNA, etc.
            # annotation['biotype'] = get_biotype(feat)
            biotype = get_biotype(feat)
            if biotype is not None:
                current_annot.other_attributes['biotype'] = biotype
            else:
                # Skip biotypes that gEAR is not currently interested in
                continue

        elif feat.type == 'CDS':
            current_annot.other_attributes['genbank_acc'] = get_genbank_acc(feat)
            current_annot.go_annotations = get_go_annotations(feat)
            current_annot.dbxrefs = get_dbxrefs(feat)

    # Don't forget the last one, rookie
    if current_annot is not None:
        annotations.append(current_annot)

    return annotations

def process_gb_gene(gene):
    annot = annotation.FunctionalAnnotation()
    annot.other_attributes['ensembl_version'] = gene.qualifiers['gene'][0]

    m = re.match('(.+)\.\d+', annot.other_attributes['ensembl_version'])
    if m:
        annot.other_attributes['ensembl_id'] = m.group(1)
    else:
        annot.other_attributes['ensembl_id'] = annot.other_attributes['ensembl_version']
        annot.other_attributes['ensembl_version'] = "{0}.0".format(annot.other_attributes['ensembl_id'])

    annot.gene_symbol = get_gene_sym(gene)
    annot.product_name = get_product(gene)
    annot.other_attributes['loc'] = get_coordinates(gene)

    return annot
    
            
if __name__ == '__main__':
    main()







