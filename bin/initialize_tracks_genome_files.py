#!/opt/bin/python

"""

The purpose of this script is to initialize the genome annotation files needed for
the gosling epigenome browser.  The files will be stored in ../www/tracks/genomes/<assembly>

These steps will be done by the pymart package
1. Go to https://useast.ensembl.org/biomart/martview
2. Choose organism dataset
3. Click "Attributes", and hit the Structures radio button
4. Click the following filters in this order (to simulate a BED file):
   1. Chromosome/scaffold name
   2. Gene start (bp)
   3. Gene end (bp)
   4. Gene name
   5. Transcript count (this won't be used)
   6. Strand
5. Check "Unique results only" and download as TSV. This is the genes file.
6. Remove header line from file
7. Change "1" strand numbers to "+", Change "-1" to "-"
   1. `awk -F'\t' -v OFS='\t' '{gsub(/-1/, "-", $6); gsub(/1/, "+", $6); print}' <assembly>.gene.txt > <assembly>.gene.fixed.txt`
8. Click the following filters in this order:
   1. Chromosome/scaffold name
   2. Exon start (bp)
   3. Exon end (bp)
   4. Gene name
   5. Transcript count (this won't be used)
   6. Strand
9. Check "Unique results only" and download as TSV. This is the exons file
10. Remove header line from file
11. Change "1" strand numbers to "+", Change "-1" to "-"
    1. `awk -F'\t' -v OFS='\t' '{gsub(/-1/, "-", $6); gsub(/1/, "+", $6); print}' <assembly>.exon.txt > <assembly>.exon.fixed.txt`
12. Convert both files to UCSC "chr" annotation as the next step needs that format
    1. python /Users/sadkins/ucsc_utils/chromToUcsc -s -i <file.bed> -o <file.renamed.bed> -a https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/bigZips/<assembly>.chromAlias.txt
    2. Works with any tab-file.  If binary, like BAM, just pipe in then back out using that tool (like samtools)
13. Time to sort both files
    1. Get chromosome ordered file from https://github.com/pkerpedjiev/negspy/blob/master/negspy/data
    2. bedtools sort -i <file.renamed.bed> -faidx <assembly>.chromInfo.txt > <file.sorted.bed.
    3. -faidx "chrom sizes" argument must define chromosome order in first column
    4. Alternatively if a chromInfo.txt file is not available, can do: `sort -k1V -k2n -k3n example.bed`
14. BGZIP each file
15. Create index for each file
    1. tabix -s 1 -b 2 -e 3 -p bed <file.bed.gz>
    2. -p will be "gff" if doing a gff file or "vcf" if that

In the end we need to have 4 files for each assembly
../www/tracks/genomes/<assembly>/<assembly>.gene.bed.gz
../www/tracks/genomes/<assembly>/<assembly>.gene.bed.gz.tbi
../www/tracks/genomes/<assembly>/<assembly>.exon.bed.gz
../www/tracks/genomes/<assembly>/<assembly>.exon.bed.gz.tbi

"""

import sys
from pathlib import Path
import urllib.request
import subprocess
import shlex

import apybiomart
# pymart errored on import for 3.11.12, so switched to apybiomart
# Also had to change the connection tester URL as the default one was dead.
from Bio import bgzf
import pandas as pd

GENOMES_ROOT = Path(__file__).resolve().parent.parent / "www" / "tracks" / "genomes"

assemblies = [
    "mm10",
    #"mm39",
    "danRer10",
    "galGal6",
    "hg19",
    "hg38",
    "rn6",
    #"calJac3"
    ]

assembly_2_biomart_dataset = {
    "mm10": "mmusculus_gene_ensembl",   # mouse
    #"mm39": "mmusculus_gene_ensembl",   # mouse
    "danRer10": "drerio_gene_ensembl",  # zebrafish
    "galGal6": "ggallus_gene_ensembl",  # chicken
    "hg19": "hsapiens_gene_ensembl",    # human
    "hg38": "hsapiens_gene_ensembl",    # human
    "rn6": "rnorvegicus_gene_ensembl",   # rat
    #"calJac3": "cjacchus_gene_ensembl"  # marmoset
}

# Some of these assemblies do not exist in the current verison of biomart
assembly_2_biomart_host = {
    "mm10": "http://nov2020.archive.ensembl.org/biomart/martservice",
    #"mm39": "http://useast.ensembl.org/biomart/martservice",
    "danRer10": "http://may2015.archive.ensembl.org/biomart/martservice",
    "galGal6": "http://apr2022.archive.ensembl.org/biomart/martservice",
    "hg19": "http://grch37.ensembl.org/biomart/martservice",
    "hg38": "http://useast.ensembl.org/biomart/martservice",
    "rn6": "http://may2021.archive.ensembl.org/biomart/martservice",
    #"calJac3": "http://may2015.archive.ensembl.org/biomart/martservice"
}

MART_NAME = "ENSEMBL_MART_ENSEMBL"

def build_chrom_to_ucsc(assembly) -> dict:
    """
    Convert Ensembl-style chromosome names to UCSC-style names using a mapping file from a URL.

    The chromAlias.txt URL has a header and the following columns:
    # sequenceName	alias names	UCSC database: mm10
    """

    # Get the URL
    chrom_alias_url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.chromAlias.txt"

    with urllib.request.urlopen(chrom_alias_url) as response:
        lines = response.read().decode('utf-8').strip().split('\n')

    # Build the mapping dictionary
    alias_dict = {}
    for line in lines[1:]:  # Skip header
        try:
            ucsc, _, _, ensembl = line.split('\t')
            alias_dict[ensembl] = ucsc
        except:
            # Occasionally there is missing data, so skip the line.
            continue
    return alias_dict

def chrom_to_ucsc(alias_dict, chrom: str) -> str | None:

    if chrom in alias_dict:
        return alias_dict[chrom]
    else:
        return None

def create_tabix_indexed_file(bgzipped_file, file_type="bed") -> None:
    """
    Tabix-indexes a genomic file.
    bgzipped_file: Path to the bgzipped input file.
    file_type: Type of the file (e.g., "bed", "vcf", "gff").
    """

    tabix_cmd = f"tabix -s 1 -b 2 -e 3 -p {file_type} {bgzipped_file}"
    subprocess.run(shlex.split(tabix_cmd), check=True)

def get_chrom_sizes(assembly: str, genome_dir: Path) -> list:
    """
    Get chromInfo.txt file from the negspy Github repository, and return the chromosomes (1st column) as a list.
    """
    chrom_sizes_url = f"https://raw.githubusercontent.com/pkerpedjiev/negspy/refs/heads/master/negspy/data/{assembly}/chromInfo.txt"
    # Save URl to geneome_dir as {assembly}.chromInfo.txt
    chrom_info_path = genome_dir / f"{assembly}.chromInfo.txt"
    urllib.request.urlretrieve(chrom_sizes_url, chrom_info_path)

    with open(chrom_info_path, 'r') as f_in:
        lines = f_in.read().strip().split('\n')

    # Extract chromosome names from the first column
    chroms = [line.split('\t')[0] for line in lines]
    return chroms

def initialize_genome_files(assembly: str):

    # Create genome directory if it doesn't exist
    genome_dir = GENOMES_ROOT / assembly
    genome_dir.mkdir(parents=True, exist_ok=True)

    chrom_dict = build_chrom_to_ucsc(assembly)

    sorted_chroms = get_chrom_sizes(assembly, genome_dir)

    dataset_name = assembly_2_biomart_dataset[assembly]

    # GENES
    print(f"-- Fetching gene annotations for {assembly}...")

    # Fetch gene annotations
    gene_attributes = [
        'chromosome_name',
        'start_position',
        'end_position',
        'external_gene_name',
        'strand'
    ]
    server = apybiomart.classes.Query(dataset=dataset_name, attributes=gene_attributes, filters={})
    server.host = assembly_2_biomart_host[assembly]
    #server.mart = MART_NAME

    gene_data: pd.DataFrame = server.query()
    # The query column namess are different from the attribute names
    gene_data.columns = gene_attributes

    print(f"-- Processing gene query results for {assembly}...")

    # Process gene data
    gene_bed_path = genome_dir / f"{assembly}.gene.bed"

    # Process exon data (vectorized)
    gene_data['chrom'] = gene_data['chromosome_name'].map(lambda c: chrom_to_ucsc(chrom_dict, c))
    gene_data = gene_data.dropna(subset=['chrom'])
    gene_data['strand'] = gene_data['strand'].map({1: '+', -1: '-', 0: '.'})
    gene_data['score'] = 0
    gene_data = gene_data[['chrom', 'start_position', 'end_position', 'external_gene_name', 'score', 'strand']]
    gene_data.columns = ['chrom', 'start', 'end', 'gene_name', 'score', 'strand']

    # Sort and write
    gene_data["chrom"] = pd.Categorical(gene_data["chrom"], categories=sorted_chroms, ordered=True)
    gene_data = gene_data.sort_values(["chrom", "start", "end"])
    gene_data.to_csv(gene_bed_path, sep='\t', header=False, index=False, columns=['chrom', 'start', 'end', 'gene_name', 'score', 'strand'])

    ### EXONS

    print(f"-- Fetching exon annotations for {assembly}...")

    # Fetch exon annotations
    exon_attributes = [
        'chromosome_name',
        'exon_chrom_start',
        'exon_chrom_end',
        'external_gene_name',
        'strand'
    ]
    server = apybiomart.classes.Query(dataset=dataset_name, attributes=exon_attributes, filters={})
    server.host = assembly_2_biomart_host[assembly]

    exon_data: pd.DataFrame = server.query()
    # The query column namess are different from the attribute names
    exon_data.columns = exon_attributes

    print(f"-- Processing exon query results for {assembly}...")

    # Process exon data
    exon_bed_path = genome_dir / f"{assembly}.exon.bed"

    # Process exon data (vectorized)
    exon_data['chrom'] = exon_data['chromosome_name'].map(lambda c: chrom_to_ucsc(chrom_dict, c))
    exon_data = exon_data.dropna(subset=['chrom'])
    exon_data['strand'] = exon_data['strand'].map({1: '+', -1: '-', 0: '.'})
    exon_data['score'] = 0
    exon_data = exon_data[['chrom', 'exon_chrom_start', 'exon_chrom_end', 'external_gene_name', 'score', 'strand']]
    exon_data.columns = ['chrom', 'start', 'end', 'gene_name', 'score', 'strand']

    # Sort and write
    exon_data["chrom"] = pd.Categorical(exon_data["chrom"], categories=sorted_chroms, ordered=True)
    exon_data = exon_data.sort_values(["chrom", "start", "end"])
    exon_data.to_csv(exon_bed_path, sep='\t', header=False, index=False, columns=['chrom', 'start', 'end', 'gene_name', 'score', 'strand'])


    print(f"-- Compressing and indexing files for {assembly}...")

    # BGZIP and index the files
    for bed_path in [gene_bed_path, exon_bed_path]:
        gz_path = bed_path.with_suffix(bed_path.suffix + '.gz')
        with bed_path.open('rb') as f_in, bgzf.BgzfWriter(gz_path) as f_out:
            # write the entirety of the input
            f_out.write(f_in.read())
        create_tabix_indexed_file(gz_path, file_type="bed")

for assembly in assemblies:
    print(f"Processing assembly {assembly}...")
    initialize_genome_files(assembly)