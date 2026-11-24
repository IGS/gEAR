# Welcome to the gEAR upload documentation!

The gene Expression Analysis Resource platform currently supports the upload of:

- Bulk RNASeq data
- Single-Cell Sequencing data
- Microarray data
- Spatial transcriptomics data
- Epigenetic data

To assist users with uploading their data, we maintain documentation for
uploading data from several different starting formats including programmatically generated files or through Excel. To find the documentation
best suited for your upload needs, please see the links below. If you have any questions, or would like assistance with uploading please contact the curator team [(Contact us)](https://umgear.org/contact.html).

> [!IMPORTANT]
>
> ## Data Requirements
>
> All data uploaded to the gEAR platform must be de-identified and contain no personally identifiable information (PII), protected health information (PHI), or other data subject to regulatory obligations.
>
> By uploading data to gEAR, you represent that:
>
> - The data being uploaded was collected, processed, stored, and shared in accordance with all applicable laws and regulations
> - There is no PII or PHI present in the data upload
> - You have obtained all necessary consents, permissions, and authorizations needed for use of the data

## Using the Dataset Uploader

Expression datasets (RNA-seq, microarray and spatial) are stored and searched at the gene level, with one row per gene. Transcript- or isoform-level counts, such as those from long-read RNA-seq, must be summed to genes before upload. To show isoform-level detail alongside them, upload coverage and transcript-model tracks as an epigenome dataset; see [Long-read and isoform data](#long-read-and-isoform-data).

Open the [Dataset Uploader](https://umgear.org/upload_dataset.html) from "Dataset Uploader" in the navigation panel (you must be logged in). The uploader walks you through these steps:

1. **Enter metadata.** Fill in the form or upload the completed [metadata template](https://umgear.org/user_templates/metadata_template.xlsx) with "Or upload a metadata file". Choose the "Dataset type" (Single-cell RNA-seq, Single-nucleus RNA-seq, Bulk RNA-seq, Microarray, Spatial or Epigenome). If your data are in GEO, enter the GEO ID (like GSEnnnnn or GSMnnnnnnn) and click "Lookup" to fill in the contact, organism, platform, instrument, library and PubMed fields automatically; check the results before continuing. Confirm that the data contain no personally identifiable information and click "Submit metadata".
2. **Upload dataset.** Choose a format -- "MEX / 3-tab format", "MS Excel", "RDS / Seurat", "H5AD / Python", "Spatial (in tar format)" or "Epigenetic Data" -- then choose your file and click "Upload dataset". Each format has a "Learn more" link with examples. (Uploading from a URL instead of a file is planned but not yet available.)
   - **MEX / 3-tab:** upload a `.tar`, `.tar.gz` or `.zip` archive (the MEX example download is a `.tar` file). For 3-tab it must contain `expression.tab`, `genes.tab` and `observations.tab`. For MEX (10x Genomics Cell Ranger output) it must contain `matrix.mtx`, `barcodes.tsv` and `genes.tsv`, or the newer gzipped `matrix.mtx.gz`, `barcodes.tsv.gz` and `features.tsv.gz`; these may sit inside a folder in the archive (for example `filtered_feature_bc_matrix/`). The genes file must have Ensembl IDs in its first column. MEX files contain no cell metadata, so the dataset starts without cell groupings or tSNE/UMAP coordinates. gEAR does not cluster uploads; to add clusters and embeddings, run an analysis in the Single Cell Workbench, or upload the data in a format that includes your own metadata (3-tab, H5AD or RDS). See [Preparing cell metadata](#preparing-cell-metadata).
   - **MS Excel:** the file must be `.xlsx`. Older `.xls` files are not accepted; open them in Excel and save as `.xlsx` first.
3. **Dataset processing.** gEAR converts and checks your file on the server. The page shows the status; click "Continue" when it finishes.
4. **Dataset post-processing.** If some metadata columns could be read either as numbers or as categories, the "Resolve ambiguous column data types" step asks you to mark each one as "continuous" or "categorical". This step is optional but helps the curators offer the right plot types.
5. **Finalize submission.** Choose whether the dataset is "Public" or "Private (can be changed later)" and click "Finalize submission".
6. **Curate dataset.** Your dataset is now in the Dataset Explorer, but it has no displays yet. Click "Curate dataset" to open the Single-gene curator and create the first displays (see [Curate data/ build plots](gEARWiki.md#curate-data-build-plots)). Spatial datasets are displayed automatically, and epigenome datasets skip this step ("View dataset").

Epigenome uploads have an extra "Build track hub" step; see [Epigenetic data (Gosling)](#epigenetic-data-gosling).

### Resuming a submission

Uploads are saved as you go. When you return to the uploader, any unfinished uploads are listed under "Submissions in progress" with their share ID, status, dataset type and title. Click "View / resume" to continue where you left off, or "Delete" to discard the upload. To begin a different upload, click "Start" under "Or start a new submission".

### If processing fails

If the "Dataset processing" step stops with an error, the message says what went wrong. The most common problem is an archive (`.tar`, `.tar.gz` or `.zip`) that can't be extracted:

- **"The uploaded archive appears to be incomplete"**: the file ends before all of its contents. This usually happens when the archive was still being written, or a copy or download was interrupted, before you uploaded it. The file on your computer may look the right size and still be cut short.
- **"The uploaded archive could not be read"**: the file is damaged or isn't really an archive (for example, a file renamed to `.tar.gz`).

Before uploading again, check that the archive opens completely on your own computer. On macOS or Linux, run one of these in a terminal:

```bash
tar -tf my_dataset.tar.gz     # for .tar or .tar.gz; lists every file, with no errors at the end
unzip -t my_dataset.zip       # for .zip; should finish with "No errors detected"
```

If the check reports an error (for example "Unexpected EOF in archive" on Linux, or "Truncated input file" or "truncated gzip input" on macOS), re-create the archive from the original files, for example with `tar -czf my_dataset.tar.gz expression.tab genes.tab observations.tab`, and check it again. Then delete the failed submission under "Submissions in progress" and upload the new archive. If the archive checks out but processing still fails, [contact us](https://umgear.org/contact.html) with the share ID shown in the uploader.

For an RDS / Seurat upload, the most common problem is a message saying the object **has no normalized expression layer named "data"**. See [Preparing a Seurat object](#preparing-a-seurat-object) for how to check and fix the object before uploading it again.

### Direct H5AD Uploads

gEAR supports uploads of H5AD datasets (those created using the Anndata structure).  If you plan to upload one of these, it is strongly recommended that the `Anndata.var` has unique identifiers as the DataFrame index, and gene symbol names in a "gene_symbol" column.  If gene symbols are in the index, we can infer Ensembl IDs based on what annotations are stored in our database. However, there is a chance Ensembl IDs may not be found for some gene symbols and they will be flagged with a generic identifier.

**Write the file with pandas 2.x.** pandas 3 stores text columns with a new string type, and H5AD files written with pandas 3 (for example with anndata 0.13) save those columns as "nullable strings". gEAR can't process these yet, and the upload fails with a message mentioning `anndata.settings.allow_write_nullable_strings`. To avoid this, write the H5AD in an environment with pandas older than 3.0 (for example pandas 2.3 with anndata 0.10, 0.11 or 0.12): read your existing file there and save it again with `adata.write_h5ad()`.

### Direct RDS / Seurat Uploads

If you have been working with your data in R, you can upload a Seurat object directly with the "RDS / Seurat" format. Save the object with `saveRDS()` (the file must end in `.rds`). gEAR converts it to H5AD on the server: the normalized "data" layer is used for expression values, dimensionality reductions (such as UMAP, t-SNE and PCA) are copied into the cell metadata so they can be plotted, and gene symbols are mapped to Ensembl IDs using the organism (taxon ID) given in your metadata. Genes whose symbols cannot be matched are given a generic identifier.

#### Preparing a Seurat object

Check these before saving the object with `saveRDS()`:

- **A normalized "data" layer.** gEAR uses the `data` layer of the object's default assay for expression values. Run `NormalizeData()` (or `SCTransform()`) before saving. The raw `counts` layer is never used in its place, so an object without a `data` layer can't be uploaded.
- **Joined layers (Seurat v5).** A Seurat v5 object that was split by sample (with `split()`, or for integration with `IntegrateLayers()`) keeps one layer per sample, named like `counts.1` and `data.1`, or `data.ctrl` and `data.treated`. Run `JoinLayers()` to merge them into single `counts` and `data` layers. Your integrated reductions and clusters are kept.
- **The default assay.** Only the default assay is converted. If it is an "integrated" assay (from Seurat v4's `IntegrateData()`), it usually holds only the 2,000 or so genes used for integration. To upload every gene, set `DefaultAssay(obj) <- "RNA"` and normalize that assay.
- **Gene symbols as feature names.** The object's feature names (`rownames(obj)`) should be gene symbols. They are mapped to Ensembl IDs using the organism in your metadata.
- **Reductions.** Each reduction stored on the object (for example `pca` and `umap`) is kept, and its first two components are added to the cell metadata as `<name>_1` and `<name>_2` (for example `umap_1` and `umap_2`) so they can be plotted.
- **Cell metadata.** The object's `meta.data` becomes the dataset's cell metadata. See [Preparing cell metadata](#preparing-cell-metadata) for the column names gEAR treats specially, such as `cluster` and `cell_type`.

To check an object in R, and fix it if needed:

```r
library(Seurat)
obj <- readRDS("my_object.rds")

DefaultAssay(obj)       # the assay gEAR converts
Layers(obj)             # needs a layer named "data"
Reductions(obj)         # for example "pca" "umap"
head(rownames(obj))     # should be gene symbols

# Only if needed:
obj <- JoinLayers(obj)      # layers are split by sample ("data.1", "data.2", ...)
obj <- NormalizeData(obj)   # there is no "data" layer

saveRDS(obj, "my_object.rds")
```

## Bulk RNAseq data

Bulk RNAseq data can be uploaded to gEAR through three paths (from Recount3 data, through Excel or via tab-deliminated files). In addition to the data files, every submission requires a standard metatdata template to be filled out with basic information about your dataset. A blank template form can be found on the upload page or linked [here](https://umgear.org/user_templates/metadata_template.xlsx).

> [!NOTE]
> Count files must be normalized prior to upload, any common normalization method is accepted

> [!TIP]
> To open the links in the flowchart, right click and choose open in new window otherwise link will show as blocked

```mermaid
graph TB

A["How do you wish to format your data?"]
A-->B["In Excel (Link)"]
A-->D["From count matrix (Link)"]
A-->E["From Recount 3 (Link to instructions)"]
A-->G["From H5AD"]
E-->F["Upload to gEAR (Link)"]
B-->F
D-->F
G-->F


click B "https://docs.google.com/presentation/d/1lU7wqWmeW907GBGfK0oBi06hcrvra-PY/edit#slide=id.p1" "Excel upload instructions"
click D "https://github.com/IGS/gEAR/wiki/Prepare-3-tab-format-for-RNAseq-dataset-uploading" "Preparing count matrix for upload"
click E "https://github.com/songeric1107/Host_data_on_gEAR/blob/225e61631db4f0a60acb3abdb90ba55ace814e87/script/prepare_from_recount3.md" "Preparing Recount3 data for upload"
click F "https://umgear.org/upload_dataset.html" "Uploading to gEAR"
style F fill:gray,stroke:black,stroke-width:2px,color:blue
style D color:blue
style E color:blue
style B color:blue
```

- [Preparing bulk RNAseq data using code (overview slides)](https://docs.google.com/presentation/d/1lYbgACVi-931EHTGNIw1bSZWIIUBcg3o/edit?usp=sharing&ouid=102015920709954238045&rtpof=true&sd=true)
  - [Preparing bulk RNAseq data (R code)](https://github.com/IGS/gEAR/wiki/Prepare-3-tab-format-for-RNAseq-dataset-uploading)
  - [Preparing RNAseq data from Recount 3 (R code)](https://github.com/songeric1107/Host_data_on_gEAR/blob/225e61631db4f0a60acb3abdb90ba55ace814e87/script/prepare_from_recount3.md)
- [Preparing bulk RNAseq data via Excel (slides)](https://docs.google.com/presentation/d/1lU7wqWmeW907GBGfK0oBi06hcrvra-PY/edit?usp=sharing&ouid=102015920709954238045&rtpof=true&sd=true)
- [Example upload files](https://drive.google.com/drive/folders/1OYZ7-FjgTBwNrZDqBI3dAA7bv8QmwFE-?usp=sharing)

## Single Cell RNAseq

Single cell sequencing data can be uploaded to gEAR through multiple paths. The most common path for data upload is creating three tab-deliminated text files which are compressed together for upload, but other upload methods are avialable.  In addition to the data files, every submission requires a standard metatdata template to be filled out with basic information about your dataset. A blank template form can be found on the upload page or linked [here](https://umgear.org/user_templates/metadata_template.xlsx).

> [!NOTE]
> Count files must be normalized prior to upload, any common normalization method is accepted

> [!TIP]
> To open the links in the flowchart, right click and choose open in new window otherwise link will show as blocked

```mermaid
graph TB

A[How do you wish to format your data?]
A-->B["Excel (Link to instructions)"]
A-->D["From Seurat object (Link)"]
A-->E["From 10X raw matrix (Link)"]
A-->G["From SingleCellExperiment object (Link)"]
G-->F["Count matrix (3 tab format)"]
A-->I["Other (3 tab format)"]
A-->J["From H5AD"]
E-->F
D-->F
F-->H["Upload to gEAR (Link)"]
B-->H
I-->F
J-->H


click B "https://docs.google.com/presentation/d/1ptk78OJAQJnyRKe3Gkejqa43Vh-cP21h/edit?usp=sharing&ouid=102015920709954238045&rtpof=true&sd=true" "Excel upload instructions"
click D "https://github.com/IGS/gEAR/wiki/Uploading-Single-Cell-Sequencing-data#b-prepare-the-files-based-on-seurat-processed-output-all-the-command-below-are-r-script" "Preparing Seurat object for upload"
click E "https://github.com/IGS/gEAR/wiki/Uploading-Single-Cell-Sequencing-data#a-prepare-the-files-using-the-raw-matrix-from-10x-platform" "Preparing 10X matrix for upload"
click F "https://umgear.org/upload_dataset.html" "Uploading to gEAR"
click G "https://github.com/IGS/gEAR/wiki/Uploading-Single-Cell-Sequencing-data#c-prepare-the-files-from-singlecellexperiment-object" "Preparing data from a SingleCellExperiment Object"
click I "https://docs.google.com/presentation/d/1wZcAy3QLEIndXXB82mBcgUzWNNUvwII_/edit?usp=sharing&ouid=102015920709954238045&rtpof=true&sd=true"

style H fill:gray,stroke:black,stroke-width:2px,color:blue
style B color:blue
style D color:blue
style E color:blue
style G color:blue
style I color:blue
```

- [Preparing single cell data via code (slides)](https://docs.google.com/presentation/d/1_YlLlQCXobkfjtIQddSZBuOo8iz3Tgnf/edit?usp=sharing&ouid=102015920709954238045&rtpof=true&sd=true)
  - [Preparing single cell data (R code)](https://github.com/IGS/gEAR/wiki/Uploading-Single-Cell-Sequencing-data)
- [Preparing single cell data via Excel (slides)](https://docs.google.com/presentation/d/1ptk78OJAQJnyRKe3Gkejqa43Vh-cP21h/edit?usp=sharing&ouid=102015920709954238045&rtpof=true&sd=true)
- Example files
  - [Processed scRNA files](https://drive.google.com/drive/folders/1LHhhCIV5LmYspjfHccYr-gD1kW-bswut?usp=sharing)
  - [Raw scRNA matrix](https://drive.google.com/drive/folders/1c6pjqj-oruNeSsYDoZtJbZv-nEmF0bcT?usp=sharing)


### Preparing cell metadata

A few metadata (observation) column names get special treatment. Add them to your observations table, the Excel "observations" sheet, `adata.obs`, or the Seurat metadata before uploading.

- **Clusters and cell types.** At upload, gEAR looks for a cluster column named `cluster`, `cell_type`, `cluster_label` or `subclass_label` and records it as the dataset's clustering in the [primary analysis](https://github.com/IGS/gEAR/blob/main/docs/wiki/gEARWiki.md#primary-vs-de-novo-analyses). gEAR does not compute clusters itself.
- **tSNE and UMAP coordinates.** Pairs of columns named `tSNE_1`/`tSNE_2` (or `tSNE1`/`tSNE2`) and `UMAP_1`/`UMAP_2` (or `UMAP1`/`UMAP2`, `uMAP_1`/`uMAP_2`) are recognized as embeddings, so tSNE/UMAP displays can be made straight away. For H5AD uploads, coordinates already stored in `adata.obsm["X_tsne"]` or `adata.obsm["X_umap"]` are used too.
- **Category colors (`_colors` columns).** To choose the color of each category, add a column named after the category column plus `_colors`, holding a hex color such as `#1f77b4` for every cell or sample. Each category must map to exactly one color. For example:

  | observations | cell_type | cell_type_colors |
  | --- | --- | --- |
  | cell1 | IHC | #1f77b4 |
  | cell2 | OHC | #ff7f0e |
  | cell3 | IHC | #1f77b4 |

  gEAR uses these colors in dynamic and static plots and as the starting colors in the curators, and doesn't list the `_colors` column as a separate category. H5AD files written by Scanpy may instead carry a color list in `adata.uns["cell_type_colors"]` (one color per category, in category order); that list is used for static tSNE/UMAP/PCA plots only. Colors picked in the curator, and colorblind mode, take precedence. See [Category colors](https://github.com/IGS/gEAR/blob/main/docs/wiki/gEARWiki.md#category-colors).
- **Replicate columns.** Columns named `replicate`, `BioRep` or `TechRep` are hidden from the filter and plot category lists, so use them for replicate IDs you don't need to plot by. The comparison tool doesn't need these columns: its statistical tests treat each sample or cell as one observation. See [Comparison tool](https://github.com/IGS/gEAR/blob/main/docs/wiki/gEARWiki.md#comparison-tool).

## Microarray data

- [Prepare GEO microarray data for upload](https://github.com/songeric1107/Host_data_on_gEAR/blob/225e61631db4f0a60acb3abdb90ba55ace814e87/script/GEO_microarray_data_to_gEAR.R)

## Spatial data

Spatial transcriptomic data from various platforms can be uploaded to gEAR. Set the "Dataset type" to "Spatial" in the metadata step. Then, in the upload step, choose "Spatial (in tar format)" and select your platform from the dropdown so gEAR knows what type of data to expect. Uploads are a tarball, either plain (".tar", as made by the `tar cvf` commands shown in the requirements) or gzipped (".tar.gz"); after selecting a platform, click "Show requirements" to see the files the tarball must contain and how they must be named.

Currently supported platforms:

- CosMx SMI
- Curio Seeker
- GeoMx DSP
- 10x Genomics Visium
- 10x Genomics Visium HD
- 10x Genomics Xenium

If you are interested in uploading spatial transcriptomic data from a platform that is not in the above list, it is recommended to follow the instructions for the Single Cell RNAseq section. As long as your metadata includes the spatial row and column (X and Y) coordinates, then those can be selected in various dataset curation plots.

Once uploaded, spatial datasets are displayed automatically when you search a gene; no curation is needed (see [Spatial displays](gEARWiki.md#spatial-displays)).

## Epigenetic data (Gosling)

Epigenetic data are uploaded and displayed using integration with the Gosling epigenome viewer. We use the [UCSC Track Hub](https://genome.ucsc.edu/goldenpath/help/hubQuickStart.html) system to allow for easy mapping to the Gosling viewer spec as well as allow for exporting to the UCSC Genome Browser for more advanced viewings.

You can upload your trackhub data from the [gEAR uploader page](http://umgear.org/upload_dataset.html).

1. Upload or enter your metadata. Make sure "Epigenome" is selected for "Dataset type".
2. In the Upload Dataset step, choose the "Epigenetic Data" type.
    - Optionally, prepopulate the track hub builder by entering a remote hub.txt URL ("Provide UCSC track hub URL to preload tracks") or by uploading a local file with "Upload hub.txt" (use one or the other, not both). Select the assembly, then click "Configure Trackhub".
    - If an uploaded hub.txt refers to track files by relative or local paths, gEAR warns you ("Hub file contains relative file references"); you will need to upload or link those track files yourself in the next step.
3. Create your hub data and your tracks data. For each track choose the "Track Type" (BigBed, BigWig, Hi-C or VCF) and either link a remote URL or upload the file.
    - To overlay several BigWig tracks (for example replicates), create a track and check "Check to make this track a multiWig container", then assign the other BigWig tracks to it with "Multiwig Group". Tracks in the group share the container's visibility setting.
4. Click "Build track hub", which starts the staging process. Some file types will also be converted to a format that Gosling can accept.
    - You can leave the page and resume later from "Submissions in progress".
5. When complete, you have the option to make this a public or a private dataset. After everything is finalized, uploading is complete!

Track hub files will be copied over to the gEAR server. This is to ensure that files are always available to be streamed by Gosling.

### Long-read and isoform data

Expression datasets in gEAR are gene-level only, but you can show isoform-level detail from long-read (or short-read) RNA-seq as tracks in an epigenome dataset:

- **Coverage (BigWig).** Convert each sample's or condition's read alignments to a BigWig coverage track. Differences in exon coverage between tracks show where isoform use changes. Group the tracks in a multiWig container to overlay them.
- **Transcript models (BigBed).** Convert your transcript models (for example a GTF from your isoform caller) to BigBed. Gosling draws each BigBed entry as a single bar from its start to its end, without separate exons. To see full exon-intron structures, use "View in UCSC Genome Browser" on the epigenome display, which opens the same track hub in the UCSC Genome Browser.

BAM files are not supported. Convert them to BigWig first.

### Currently supported UCSC Track Hub fields

Any settings beyond the following will be ignored.

- track
- bigDataUrl
- shortLabel
- longLabel
- color
- visibility
- type
- container (multiWig only)
- parent (for tracks inside a multiWig container)

For the hub.txt file, we also accept the newer "useOneFile=on" single-file format. If your hub does not use this setting, ensure that the trackDb.txt file for the assembly you wish to use is present and listed in genomes.txt, as gEAR will only read the tracks for the selected assembly. Tracks for other assemblies are ignored. Multi-file hubs (without "useOneFile=on") can only be prepopulated from a hub URL, not from an uploaded hub.txt file.

### Currently supported dataset formats (Track Hub "type")

- bigBed
- bigWig
- hic
- vcfTabix

BigWig tracks can also be grouped in multiWig containers (the "container multiWig" and "parent" settings), which are drawn overlaid.

### Currently supported reference genome codes

- danRer10
- galGal6 (COMING SOON)
- hg19
- hg38
- mm10
- mm39
- rn6

If you want to upload epigenome data from one of the human references, because of concerns with personally-identifiable information, Hi-C and VCF file types will not be processed.

