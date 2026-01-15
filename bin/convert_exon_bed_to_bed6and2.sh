#!/bin/bash

# For a given BED6 file of genes, convert to BED6+2 format by adding exons
# from a BED6 exon.bed file as a list of exonStarts and exonEnds

# Assuming gene and exon BED input files are sorted by chromosome and start position.
# It should have at least 6 columns: chrom, chromStart, chromEnd, name (transcript ID), score, strand.

set -euo pipefail

# Usage: convert_exon_bed_to_bed12.sh <gene_bed6> <exon_bed6> <out_bed12>
gene_bed="$1"
exon_bed="$2"
bed12_out="$3"

# Ensure exon entries are grouped by transcript name and sorted by start.
# Use process substitution to feed sorted exons as the first file to awk,
# and the gene_bed as the second file. AWK collects exons keyed by col4,
# then iterates genes and emits BED12. Genes with no exons emit a single-block BED12.

awk -F'\t' '
BEGIN { OFS = "\t" }

# First file: sorted exon_bed (via process substitution)
FNR==NR {
    # skip lines missing required fields
    if ($1 == "" || $2 == "" || $3 == "" || $4 == "") next
    # ensure start/end are integers
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) next

    name = $4
    # record chrom/strand for the transcript (first occurrence)
    if (!(name in exonChrom)) {
        exonChrom[name] = $1
        exonStrand[name] = ($6 == "" ? "+" : $6)
    }
    exonCount[name]++
    i = exonCount[name]
    exonStart[name, i] = $2
    exonEnd[name, i]   = $3
    next
}

# Second file: gene_bed
{

    # skip lines missing required fields
    if ($1 == "" || $2 == "" || $3 == "" || $4 == "") next
    # ensure start/end are integers
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) next

    g_chrom = $1
    g_start = $2
    g_end   = $3
    g_name  = $4
    g_score = ($5 == "" ? 0 : $5)
    g_strand = ($6 == "" ? "+" : $6)

    if (exonCount[g_name] > 0) {
        count = exonCount[g_name]
        chrom = exonChrom[g_name]
        chromStart = exonStart[g_name,1]
        chromEnd   = exonEnd[g_name,count]

        starts = ""
		ends = ""
        for (i = 1; i <= count; i++) {
            start = exonStart[g_name,i]
			end = exonEnd[g_name,i]
            starts = starts ((i==1) ? "" : ",") start
			ends = ends ((i==1) ? "" : ",") end
        }
        print chrom, chromStart, chromEnd, g_name, g_score, g_strand, starts, ends
    } else {
        # Fallback: emit single block using gene coordinates
        chrom = g_chrom
        chromStart = g_start
        chromEnd = g_end
        starts = g_start
		ends = g_end
        print chrom, chromStart, chromEnd, g_name, g_score, g_strand, starts, ends
    }
}
' <(sort -k4,4 -k2,2n "$exon_bed") "$gene_bed" > "$bed12_out"
# ...existing code...
