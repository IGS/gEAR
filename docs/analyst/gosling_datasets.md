# Dataset information for Gosling

Gosling (https://gosling-lang.org/) is a grammar-based toolkit for scalable and interactive genomics data visualization. We use this in the gEAR platform to display various types of epigenetic data.

Some of the features we have implemented in our use of Gosling include:

* Genome and regional position indicators
* Support for BED, BigWig, and VCF files
* ENSEMBL-based gene and exon annotation for a few common reference genomes
* In expanded view:
  * Side-by-side comparisons (i.e. promoter regions, gene vs gene)
  * Quick search a second gene
  * Export to UCSC Genome Browser for more advanced control.

## Currently supported reference genome codes

Subject to additions

* danRer10
* galGal6 (currently having issues getting the annotation uploaded to HiGlass)
* hg19
* hg38
* mm10
* rn6

## Manual dataset upload process

Before uploading, you need to do the following:

* Populate an upload metadata spreadsheet and scp it onto the server (I generally use /tmp)
  * Ensure the dataset type is set to "gosling"
* Create your Track Hub directory structure
  * Archive that up, and scp it onto the server.
  * Some files need to be ingested into HiGlass. @adkinsrs can help with this

On the server to upload data, do the following:

1. Add the metadata to the database to obtain a dataset ID
  a. `cd <gear_root>/bin; /opt/Python-3.10.4/bin/python3 ./add_excel_metadata_to_db.py -i <metadata_file>.xlsx -oi <user_id>`
2. In MySQL, insert an entry into the dataset_display table, so this becomes the owner's primary display
  a. `insert into dataset_display (dataset_id, user_id, label, plot_type, plotly_config) values (<dataset_id>, <user_id>, "<unused_display_label>", "gosling", '{"hubUrl":"http://umgear.org/tracks/<dataset_name>/hub.txt", "assembly":"mm10", "gene_symbol":"Atoh1"}');`
  b. `insert into dataset_preference (user_id, dataset_id, display_id) values (<user_id>, <dataset_id>, <new_display_id>);`
  c. A gene symbol must be added to the "plotly_config" JSON but mostly because UI code expects a gene.
  d. Replace "mm10" with your organism of choice
  e. It should be correct, but if the `dataset.has_h5ad` property is not set to 0, update that.
3. Move and unarchive Track Hub directory structure within `<gear_root>/www/tracks/`

## Building the Track Hub

Gosling uses the UCSC Track Hub (https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html) syntax to map Track Hub properties to those definition mappings within the Gosling interface. Typically files will be stored locally on the server, but some cases (i.e. HiC) data will need to be stored on the HiGlass server and loaded as URLs in the tracksDb.txt file.

### Dataset organization

```

<dataset_name>/
├── hub.txt
├── genomes.txt                   # Info must match directory name for genome
├── mm10/                         # Or other genome
│   └── file1.bw                  # Example track file
│   └── file2.bw
│   └── trackDb.txt                  # Contains track information for all tracks

```

### Supported trackDb.txt properties

* track
* bigDataUrl
* shortLabel
* longLabel
* color
* visibility
* type
* container multiwig
* parent (if using a container to aggregate tracks. Must indent sub-track information)

The trackDb.txt file also allows for custom metadata tags, which will be subject to future additon.

* gos_higlass_url (filled in during upload if type is hic)

### Special considerations for certain file types

For any of these cases, the "type" property in the trackDb file should point to the file type that can be used in the UCSC Genome Browser.  Both the Gosling equivalent and the UCSC equivalent need to be present.

#### BED

* Gosling only accepts BED type
  * Gosling can also accept BEDDB (which is ingested into HiGlass) but I am staying away from this for now.
* UCSC only accepts BigBed type

#### HiC

* Gosling only accepts .mcool as a HiGlass ingest
* UCSC only accepts .hic

## How to upload data to the HiGlass Server

This currently applies to HiC data
