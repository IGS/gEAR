#!/usr/bin/env python3

import argparse
import os
import subprocess
import tempfile
import pandas as pd

"""
This script converts STAR SJ.out.tab files into bigInteract format for visualization in the UCSC Genome Browser.

The "star_file" must be removed of chromosome entries where the chromosome is not in the chrom-sizes file
"""

def star_to_interact(star_file, output_dir):
    # STAR columns: chrom, start, end, strand, motif, annot, unique_reads, multi_reads, max_overhang
    df = pd.read_csv(star_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'strand', 'motif', 'annot', 'uniq', 'multi', 'overhang'])

    # Filter out junctions with low read support if desired
    df = df[df['uniq'] > 0]

    # open up a tempfile BED file
    bedfile = output_dir + "/temp.bed"

    with open(bedfile, 'w') as f:
        for idx, row in df.iterrows():
            chrom = row['chrom']
            # STAR intron coordinates are 1-based, convert to 0-based BED if necessary
            start = int(row['start'])
            end = int(row['end'])
            strand = int(row['strand'])
            score = min(1000, int(row['uniq'] * 10)) # Scale score for display density
            val = row['uniq']

            # Construct the 18 columns for bigInteract
            line = f"{chrom}\t{start}\t{end}\tjunc_{idx}\t{score}\t{val}\t.\t128,0,128\t" \
                   f"{chrom}\t{start}\t{start+1}\tdonor\t{strand}\t" \
                   f"{chrom}\t{end-1}\t{end}\tacceptor\t{strand}\n"
            f.write(line)
    return bedfile

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert STAR SJ.out.tab to bigInteract BED format.")
    parser.add_argument("--star_file", help="Input STAR SJ.out.tab file")
    parser.add_argument("--chromsizes_file", help="Path to genome chrom-sizes file.")
    parser.add_argument("--output_file", help="Output bigInteract file")
    parser.add_argument("--bedtobigbed_path", help="Path to UCSC bedToBigBed utility (includes executable name)")

    args = parser.parse_args()

    bed_file = star_to_interact(args.star_file, os.path.dirname(args.output_file))

    # sort the output BED file and convert to bigInteract using bedToBigBed
    subprocess.run(f"sort -k1,1 -k2,2n {bed_file} > {bed_file}.sorted", shell=True)
    subprocess.run(f"{args.bedtobigbed_path} -as=./interact.as -type=bed5+13 {bed_file}.sorted {args.chromsizes_file} {args.output_file}", shell=True)