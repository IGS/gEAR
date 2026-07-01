# Preparing and adding files to HiGlass

## Acquiring Gene and Exon Annotation Files

This creates a BED file where the exons are set in blockSizes and blockStarts columns

1. Go to [https://useast.ensembl.org/biomart/martview](https://useast.ensembl.org/biomart/martview)
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
   1. `awk -F'\t' -v OFS='\t' '{gsub(/-1/, "-", $6); gsub(/1/, "+", $6); print}' mm10.gene.txt > mm10.gene.fixed.txt`
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
    1. `awk -F'\t' -v OFS='\t' '{gsub(/-1/, "-", $6); gsub(/1/, "+", $6); print}' mm10.exon.txt > mm10.exon.fixed.txt`
12. Convert both files to UCSC "chr" annotation as the next step needs that format
    1. python \~/ucsc\_utils/chromToUcsc \-i \<file.bed\> \-o \<file.renamed.bed\> \-a [https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chromAlias.txt](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chromAlias.txt)
    2. If script errors because of a chromosome that wasn't found, pass in \-s to the command.
    3. Works with any tab-file.  If binary, like BAM, just pipe in then back out using that tool (like samtools)
13. Time to sort both files
    1. Get chromosome ordered file from [https://github.com/pkerpedjiev/negspy/blob/master/negspy/data](https://github.com/pkerpedjiev/negspy/blob/master/negspy/data/mm10/chromInfo.txt)
    2. bedtools sort \-i \<file.renamed.bed\> \-faidx \<assembly\>.chromInfo.txt \> \<file.sorted.bed.
    3. \-faidx "chrom sizes" argument must define chromosome order in first column
    4. This will error if there is a chromosome (i.e. random scaffold) in the bed file but not in the chromInfo file.
    5. Genomes should have the same chromosomes per version (i.e. mm10 to mm39), so you can just update the sizes. Do not leave an empty line in the file.
       1. This is assuming you don't care about the random scaffolds, which will have different IDs. Just delete the scaffolds
       2. Sizes can be obtained from UCSC [https://hgdownload.gi.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes](https://hgdownload.gi.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes) example using mm39
    6. The chromInfo.txt file needs to also be copied to the "https://umgear.org/tracks/genomes/<assembly>"" directory
14. Combine both files into one BED12 file
    1. Run `gEAR/bin/convert_exon_bed_to_bed12.sh gene.bed exon.bed gene6and2.bed`
    2. For each gene, this adds all exons as a list of exonStart and exonEnd columns.  A lot of Gosling examples use this (see "Way 1" above)
    3. I renamed these as {genome}.annotation.bed
15. tar up all the annotation.bed files and send to the higlass-manage-prod server… unzip
    1. Or just pass the single bed file instead if copying just one.
    2. If genome (i.e. mm39) is not in the negspy list of genomes, copy the chromInfo.txt file you made with the chromosome order and sizes.
16. Turn into a beddb file
    1. scp onto higlass-manage-prod, and ssh in
    2. `cd \~/higlass/bin`
    3. `source .activate` to activate the venv
    4. If genome was not in negspy you need to ingest the chromInfo.txt you copied over
       1. I don't think this is actually needed.
       2. `higlass-manage ingest \--filetype chromsizes-tsv \--datatype chromsizes \--assembly {genome} /tmp/{genome}.chromInfo.txt`
    5. `clodius aggregate bedfile –assembly={genome} {genome}.annotation.bed`
       1. If using a custom assembly, add "--chromsizes-filename \<filepath\>" to the command.
    6. I was not able to run clodius on galGal6
    7. The "assembly" argument genomes are retrieved from negspy
17. Ingest with Higlass
    1. `higlass-manage ingest {genome}.annotation.bed.beddb`

## Hi-C file ingest (manual)

NOTE: These commands are run as part of the epigenome uploader in gEAR. You do not need to run these manually.

Example files at https://server.gosling-lang.org/api/v1/tilesets/?t=cooler

1. Convert .hic to .mcool (multi-resolution Cooler file)
   1. hic2cool.hic2cool\_convert(\<infile\>, \<outfile\>, 0\) ([https://github.com/4dn-dcic/hic2cool/](https://github.com/4dn-dcic/hic2cool/))
2. Ingest .mcool file to use in HiGlass
   1. Check [https://docs.higlass.io/higlass\_server.html\#uploading-data-post](https://docs.higlass.io/higlass_server.html#uploading-data-post)

## STAR splice junction file ingest (manual)

### Convert to BEDDB (Gosling)

The STAR splice junction tab output file is fine for this.

1. scp the file onto the "higlass-manage-prod" server, then ssh in
2. `cd \~/higlass/bin`
3. `source .activate` to activate the venv
4. `clodius aggregate bedfile –assembly={genome} --chromsizes-filename={genome}.chromInfo.txt SJ.out.tab`
5. `higlass-manage ingest SJ.out.tab.beddb`

### Convert to BigInteract (UCSC)

I wrote a script in "bin" to do this.

`python ~/git/gEAR/bin/convert_STAR_sj_tab_into_biginteract.py --star_file=<SJ_out.tab> --chromsizes_file=<genome.chromSizes.txt> --output_file=<SJ.biginteract> --bedtobigbed_path=/path/to/bedToBigBed`

* star_file = Splice junctions tab file from STAR. Remove any chromosome entries that are not in the chromSizes file
* chromsizes_file = List of sorted chromosomes with chromosome sizes. Should be available per assembly.  Ask @adkinsrs if you can't find it.
* output_file = path for the bigInteract file. Extension suggestions can be .bi or .biginteract
* bedtobigbed_path = path to the UCSC Utils bedToBigBed executable.