#!/usr/bin/env python3

"""
Genomic data got from here:

https://www.gencodegenes.org/human/release_19.html
 -> Comprehensive gene annotation -> GFF3

wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

"""
from biocode import gff, things
import json
import re
import sys

# biocode Assembly objects
mols = dict()

# biocode Region objects with other stuff too
loci = dict()

# biocode.Gene objects by assembly
genes = dict()

for line in open('hl_loci.tab'):
    line = line.rstrip()
    cols = line.split("\t")
    cols[4] = cols[4].replace(',', '')

    m = re.match("(chr.+):(\d+)\-(\d+)", cols[4])
    if m:
        (chrom, start, stop) = m.groups()

        if chrom not in mols:
            print("Creating Assembly object for {0}".format(chrom), file=sys.stderr)
            mols[chrom] = things.Assembly(id=chrom)
        
        loci[cols[0]] = {
            'omim_link': cols[1],
            'paper_author': cols[2],
            'paper_link': cols[3],
            'chrom': chrom,
            'region': things.Region()
        }
        loci[cols[0]]['region'].locate_on(target=mols[chrom], fmin=int(start) - 1, fmax=int(stop), strand=1)

for line in open('hg19.refGene.txt'):
    line = line.rstrip()
    cols = line.split("\t")

    if len(cols) != 16:
        continue

    # skip if this is on a chromosome we don't have ranges for
    if cols[2] not in mols:
        continue

    gene_symbol = cols[12]
    gene = things.Gene(id=gene_symbol)
    gene.locate_on(target=mols[cols[2]], fmin=int(cols[4]), fmax=int(cols[5]), strand=1)

    if cols[2] not in genes:
        genes[cols[2]] = list()

    genes[cols[2]].append(gene)
        
# loop each region, then print out the genes within that region
data = dict()
for locus in loci:
    chrom = loci[locus]['chrom']
    for gene in genes[chrom]:
        if gene.overlaps_with(loci[locus]['region']):
            gene_loc = gene.location()
            #print("Gene {0} ({2}:{3}-{4}) appears to fall within locus {1}".format(
            #    gene.id, locus, gene_loc.on.id, gene_loc.fmin + 1, gene_loc.fmax))

            if gene.id.lower() not in data:
                data[gene.id.lower()] = {
                    'locus': locus,
                    'links': [{'label': 'OMIM', 'url': loci[locus]['omim_link']},
                              {'label': "Publication: {0}".format(loci[locus]['paper_author']), 'url': loci[locus]['paper_link']}
                    ]
                }
            
print(json.dumps(data, indent=3))
print("Finished.".format(), file=sys.stderr)
        




















