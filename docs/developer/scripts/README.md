# gEAR Utility Scripts Documentation

The `/bin` directory contains 110 utility scripts for various tasks related to data management, conversion, validation, and administration. This guide categorizes these scripts and provides usage scenarios.

## Quick Reference

- [Data Conversion & Format](#data-conversion--format) - Converting between data formats
- [H5AD Manipulation](#h5ad-manipulation) - Modifying H5AD files
- [Database & Loading](#database--loading) - Loading data into the system
- [Gene Annotations](#gene-annotations) - Working with gene symbols and annotations
- [Dataset Management](#dataset-management) - Managing datasets in gEAR
- [Validation & Testing](#validation--testing) - Validating data and testing
- [Visualization & SVG](#visualization--svg) - Display and visualization utilities
- [Profiling & Statistics](#profiling--statistics) - Performance and usage metrics
- [Administration & Maintenance](#administration--maintenance) - System administration

## General Usage Notes

- Most scripts require the gEAR library path to be accessible
- Scripts typically include help text: `script_name.py --help`
- Many scripts require database access (configure `gear.ini` first)
- Always test on non-production data first
- Some scripts may be obsolete or deprecated (marked where known)

---

## Data Conversion & Format

Scripts for converting between different data formats (H5AD, 3-tab, MEX, Excel, etc.)

### Core Conversion Scripts

#### `convert_3tab_to_h5ad.py`
Converts 3-tab format (expression.tab, genes.tab, observations.tab) to H5AD format.

**Use case**: Converting uploaded datasets to gEAR's H5AD format
```bash
./bin/convert_3tab_to_h5ad.py -i /path/to/3tab/directory -o output.h5ad
```

#### `convert_3tab_to_h5ad_lowmem.py`
Low-memory version of the above, processes data in chunks.

**Use case**: Converting very large datasets that don't fit in memory
```bash
./bin/convert_3tab_to_h5ad_lowmem.py -i /path/to/3tab/directory -o output.h5ad -r 500
```
- `-r`: Row chunk size (default 500)

#### `h5ad_convert_from_MEX.py`
Converts Matrix Exchange (MEX) format to H5AD.

**Use case**: Converting 10X Genomics MEX output to gEAR format
```bash
./bin/h5ad_convert_from_MEX.py -i /path/to/mex/directory -o output.h5ad
```

#### `h5ad_convert_to_3tab.py`
Converts H5AD back to 3-tab format for export or manual editing.

**Use case**: Exporting datasets for external analysis
```bash
./bin/h5ad_convert_to_3tab.py -i input.h5ad -o /output/directory
```

#### `convert_zarr_datasets_to_h5ad.py`
Converts Zarr spatial datasets to H5AD format for analysis.

**Use case**: Converting preprocessed spatial datasets to H5AD for downstream analysis
```bash
./bin/convert_zarr_datasets_to_h5ad.py
```
**Note**: Processes all Zarr datasets in `www/datasets/spatial/` directory

### Format Translation

#### `create_MEX_from_gear_standard_tab.py`
Creates MEX format from gEAR standard tab format.

**Use case**: Preparing data for tools that require MEX format

#### `create_excel_from_gear_standard_tab.py`
Creates Excel files from gEAR standard tab format.

**Use case**: Creating human-readable exports

#### `convert_legacy_dataset_to_h5ad.py` / `convert_legacy_dataset_to_3tab.py`
Converts old gEAR dataset formats to current formats.

**Use case**: Migrating old datasets to new format
**Status**: May be obsolete if all datasets already migrated

### Metadata Handling

#### `add_excel_metadata_to_db.py`
Loads dataset metadata from Excel template into database.

**Use case**: Batch uploading dataset metadata
```bash
./bin/add_excel_metadata_to_db.py -i metadata_template.xlsx
```

#### `convert_excel_metadata_to_json.py`
Converts Excel metadata to JSON format.

**Use case**: Preparing metadata for API consumption

---

## H5AD Manipulation

Scripts for modifying existing H5AD files without full conversion.

### Adding/Modifying Data

#### `add_analysis_tSNE_to_h5ad.py`
Adds tSNE coordinates to H5AD file's obsm layer.

**Use case**: Adding pre-computed tSNE results
```bash
./bin/add_analysis_tSNE_to_h5ad.py -i input.h5ad -t tsne_coords.tab
```
- Expects columns named tSNE1 and tSNE2

#### `h5ad_add_colors_column.py`
Adds color information for categorical data.

**Use case**: Preparing datasets for visualization with custom colors

#### `add_ensembl_ids_to_h5ad.py`
Adds Ensembl gene IDs to H5AD var (gene) annotations.

**Use case**: Enriching gene annotations with Ensembl IDs
```bash
./bin/add_ensembl_ids_to_h5ad.py -i input.h5ad -o output.h5ad
```

#### `add_ensembl_id_to_h5ad_missing_release.py`
Finds best matching Ensembl release and adds IDs, creating fake IDs for unmapped genes.

**Use case**: Handling datasets with only gene symbols
```bash
./bin/add_ensembl_id_to_h5ad_missing_release.py -i input.h5ad -o output.h5ad
```

#### `add_ensembl_id_to_h5ad_missing_release__file_based.py`
File-based version that avoids memory issues with large datasets.

**Use case**: Processing very large datasets
**Note**: Requires `mapped/`, `unmapped/`, `merged/` directories in cwd

### Modifying Observations/Variables

#### `h5ad_rename_obs_column.py`
Renames columns in the observations (cells) metadata.

**Use case**: Standardizing column names
```bash
./bin/h5ad_rename_obs_column.py -i input.h5ad -o output.h5ad --old "OldName" --new "NewName"
```

#### `h5ad_rename_obs_index.py`
Renames the observation index (cell IDs).

**Use case**: Fixing cell barcode formatting

#### `h5ad_replace_column_value.py`
Replaces values in observation columns.

**Use case**: Fixing categorical labels or typos

#### `h5ad_duplicate_column.py`
Duplicates a column in observations metadata.

**Use case**: Creating backup of column before modification

### Fixing Issues

#### `h5ad_fix_numeric_headers.py`
Fixes issues with numeric column headers.

**Use case**: Resolving upload errors from numeric headers

### Testing & Analysis

#### `h5ad_test_read.py`
Tests if H5AD file can be read successfully.

**Use case**: Validating H5AD file integrity
```bash
./bin/h5ad_test_read.py -i dataset.h5ad
```

#### `h5ad_summarize_observations.py`
Prints summary statistics about observations.

**Use case**: Quick dataset overview
```bash
./bin/h5ad_summarize_observations.py -i dataset.h5ad
```

---

## Database & Loading

Scripts for loading data into the gEAR database.

### Annotation Loading

#### `load_genbank_gbk_annotations.py` / `load_genbank_gbk_annotations_v2.py`
Loads gene annotations from GenBank GBK files.

**Use case**: Initial setup - loading organism annotations
```bash
./bin/load_genbank_gbk_annotations.py -i /path/to/genbank -id 1 -r 94
```
- `-id`: Organism ID
- `-r`: Release number

#### `load_ensembl_gbk_annotations.py`
Loads annotations from Ensembl GBK files.

**Use case**: Loading/updating Ensembl annotations

#### `load_gene_ontology.py`
Loads Gene Ontology terms and relationships.

**Use case**: Initial setup - populating GO database
```bash
./bin/load_gene_ontology.py -i go.obo
```

#### `load_gene_synonyms.py`
Loads gene synonym mappings.

**Use case**: Enabling gene symbol search with synonyms

### External Data Loading

#### `load_komp_repository_data.py`
Loads data from KOMP (Knockout Mouse Project) repository.

**Use case**: Integrating KOMP knockout data

#### `load_mgi_mirna_aliases.py` / `load_mirna_fam_data.py`
Loads miRNA annotation data from MGI.

**Use case**: Supporting miRNA datasets

#### `load_zfin_id_urls.py`
Loads ZFIN (zebrafish) database URLs.

**Use case**: Supporting zebrafish datasets with external links

### Dataset Operations

#### `upload_spatial_dataset.py`
Uploads spatial transcriptomics datasets.

**Use case**: Loading Visium or other spatial datasets
```bash
./bin/upload_spatial_dataset.py -i spatial_data/ -o dataset_id
```
See also: `docs/uploading_spatial_dataset.md`

#### `FACS_RNAseq_Data.loader.py`
Specialized loader for FACS RNA-seq data.

**Use case**: Loading sorted cell population RNA-seq data
**Note**: May be dataset-specific

### Utility Loaders

#### `load_old_gcid_gene_symbols.py` / `load_old_gcid_gene_symbols.v2.py`
Loads gene symbols from old GCID system.

**Use case**: Migration from legacy system
**Status**: Likely obsolete

---

## Gene Annotations

Scripts for working with gene symbols, Ensembl IDs, and other annotations.

### Finding Best Ensembl Release

#### `find_best_ensembl_release_from_h5ad.py`
Analyzes H5AD file to find which Ensembl release best matches the genes.

**Use case**: Identifying correct annotation version for dataset
```bash
./bin/find_best_ensembl_release_from_h5ad.py -i dataset.h5ad
```

#### `find_best_ensembl_release_match.py`
Similar to above but works with gene lists.

**Use case**: Checking annotation compatibility

### Gene Symbol Operations

#### `add_gene_symbol_to_file.py`
Adds gene symbols to tab-delimited file (assumes Ensembl ID in first column).

**Use case**: Enriching data files with gene symbols
```bash
./bin/add_gene_symbol_to_file.py -i input.tab -o output.tab
```

#### `add_ensembl_ids_to_tab_file.py`
Adds Ensembl IDs to tab file with gene symbols.

**Use case**: Reverse of above - adding IDs to symbol-only files

#### `check_h5ad_for_missing_gene_symbol.py`
Checks if H5AD files have gene_symbol column in .var.

**Use case**: Dataset validation
```bash
./bin/check_h5ad_for_missing_gene_symbol.py -d /path/to/h5ad/directory
```

### Annotation Processing

#### `create_annotation_orthology_maps.py`
Creates orthology mappings between species.

**Use case**: Cross-species analysis support

#### `create_hdf5_from_ensembl_gbk.py`
Creates HDF5 annotation files from Ensembl GenBank.

**Use case**: Preparing annotation data for gEAR

#### `parse_genbank_annotations.py`
Parses GenBank annotation files.

**Use case**: Processing annotation data before loading

#### `update_gene_coordinates_from_genbank.py`
Updates gene coordinates in database from GenBank.

**Use case**: Refreshing coordinate data for new genome releases

### Scanning & Reporting

#### `scan_datasets_for_ensembl_and_gene_symbols.py`
Scans all datasets to report which have Ensembl IDs and symbols.

**Use case**: Auditing dataset annotations across system
```bash
./bin/scan_datasets_for_ensembl_and_gene_symbols.py > annotation_report.txt
```

---

## Dataset Management

Scripts for managing datasets within gEAR.

### Dataset Operations

#### `add_datasets_to_group.py`
Adds multiple datasets to a dataset collection (profile).

**Use case**: Bulk organization of datasets into collections
```bash
./bin/add_datasets_to_group.py -p profile_id -d dataset1,dataset2,dataset3
```

#### `remove_marked_datasets.py`
Removes datasets marked for deletion.

**Use case**: Cleanup of flagged datasets
**Warning**: Ensure datasets are backed up before running

#### `change_ownership.py`
Changes ownership of datasets or profiles.

**Use case**: Transferring datasets between users
```bash
./bin/change_ownership.py --dataset dataset_id --new_owner user_id
# or
./bin/change_ownership.py --profile profile_id --new_owner user_id --include_displays
```
**Warning**: Use with care - affects access permissions

### Reporting

#### `report_dataset_dimensions.py`
Reports dimensions (genes x observations) of datasets.

**Use case**: Quick inventory of dataset sizes
```bash
./bin/report_dataset_dimensions.py > dataset_sizes.txt
```

#### `report_dataset_usage_stats.py`
Reports usage statistics for datasets.

**Use case**: Understanding which datasets are most accessed
```bash
./bin/report_dataset_usage_stats.py -o stats_report.csv
```

#### `list_datasets.py`
Lists all datasets in the system.

**Use case**: Getting dataset inventory
```bash
./bin/list_datasets.py --format csv > datasets.csv
```

### Updates

#### `update_dataset_long_description.py`
Updates dataset long descriptions in bulk.

**Use case**: Maintaining dataset documentation

#### `update_dataset_table_schematic_images.py`
Updates schematic images for datasets.

**Use case**: Refreshing dataset visualizations

---

## Validation & Testing

Scripts for validating data integrity and testing functionality.

### File Validation

#### `validate_tab_file.py`
Validates tab-delimited files for common issues.

**Use case**: Pre-upload validation
```bash
./bin/validate_tab_file.py -i data.tab
```

#### `validate_3tab_xls_gene_symbols.py`
Validates gene symbols in 3-tab/Excel files.

**Use case**: Ensuring gene symbols match database

#### `validate_file_extension.py`
Checks file extensions match content.

**Use case**: Detecting renamed files

#### `byte_validator.py`
Checks for non-UTF8 characters that can break upload.

**Use case**: Pre-upload validation for encoding issues
```bash
./bin/byte_validator.py -i data.tab
```

### Data Validation

#### `validate_expression_observations_correlation.py`
Validates that expression matrix matches observations.

**Use case**: Ensuring matrix dimensions are correct

#### `validate_observation_ids.py`
Validates observation (cell) IDs are unique and valid.

**Use case**: Pre-upload validation

#### `check_singlecell_column_count.py`
Returns total columns and any rows with mismatched counts.

**Use case**: Finding corrupted rows in single-cell data
```bash
./bin/check_singlecell_column_count.py -i expression.tab
```

### System Testing

#### `test_database_connection.py`
Tests database connectivity.

**Use case**: Setup validation
```bash
./bin/test_database_connection.py
```

#### `test_mysql_stringtype.py`
Tests MySQL string type handling.

**Use case**: Debugging MySQL configuration issues

---

## Visualization & SVG

Scripts for managing visualizations, displays, and SVG anatomy diagrams.

### SVG Manipulation

#### `copy_svg_ids_to_class_names.py` / `copy_svg_group_ids_to_class_names.py`
Copies SVG element IDs to class names for styling.

**Use case**: Preparing SVG anatomical diagrams
```bash
./bin/copy_svg_ids_to_class_names.py -i input.svg -o output.svg
```

#### `set_svg_class_names.py`
Sets class names on SVG elements.

**Use case**: Bulk SVG element classification

#### `replace_anatomy_ids_with_class_names.py`
Replaces anatomy IDs with standardized class names.

**Use case**: Standardizing anatomy SVG files

### Display Management

#### `generate_static_display_images.py`
Generates static preview images for displays.

**Use case**: Creating dataset thumbnails
```bash
./bin/generate_static_display_images.py -d dataset_id
```

#### `remove_duplicate_layout_displays.py`
Removes duplicate display entries.

**Use case**: Database cleanup after duplicates created

#### `fix_missing_display_config_options.py`
Fixes displays missing configuration options.

**Use case**: Migration after schema changes

### Color & Styling

#### `generate_color_gradient_css.py`
Generates CSS for color gradients.

**Use case**: Creating color schemes for heatmaps

#### `rescore_dataset_coloring_v3.py`
Rescores dataset coloring based on new algorithms.

**Use case**: Updating to improved color schemes

#### `rescore_dataset_scope_dataset_coloring.py` / `rescore_gene_scope_dataset_coloring.py` / `rescore_tissue_scope_dataset_coloring.py`
Rescores coloring at different scopes.

**Use case**: Optimizing color mapping for specific views

### Epigenome Browser (Gosling)

#### `migrate_epiviz_to_gosling.py`
Migrates saved Epiviz curations to Gosling format.

**Use case**: One-time migration from legacy Epiviz displays to Gosling
```bash
./bin/migrate_epiviz_to_gosling.py
```
**Note**: Builds hub, genomes, tracksdb, and groups files for Gosling browser
**Status**: Used during Epiviz to Gosling transition

#### `initialize_tracks_genome_files.py`
Initializes genome annotation files for Gosling epigenome browser.

**Use case**: Setting up genome annotations for Gosling visualization
```bash
./bin/initialize_tracks_genome_files.py --assembly <assembly_name>
```
**Note**: Downloads and formats genome data from Ensembl via BioMart

#### `aggregate_gosling_bigwig_groups.py`
Aggregates bigWig files for Gosling group tracks.

**Use case**: Creating group-level tracks by averaging bigWig signals
```bash
./bin/aggregate_gosling_bigwig_groups.py --hub <hub_file>
```
**Note**: I/O intensive; run as background job for large genomes

---

## Profiling & Statistics

Scripts for performance profiling and generating statistics.

### Performance Profiling

#### `profile_single_heatmap_run.py`
Profiles a single heatmap rendering for performance.

**Use case**: Debugging slow heatmap rendering
```bash
./bin/profile_single_heatmap_run.py -d dataset_id -g gene_list
```

#### `profile_single_projectr_tsne_run.py`
Profiles ProjectR tSNE computation.

**Use case**: Debugging ProjectR performance issues

### Usage Statistics

#### `get_site_stats.py`
Generates overall site usage statistics.

**Use case**: Reporting to stakeholders
```bash
./bin/get_site_stats.py -o site_stats.json
```

#### `itemize_profile_metrics.py`
Generates detailed metrics for profiles (dataset collections).

**Use case**: Understanding profile usage patterns

### ID Generation

#### `generate_profile_share_ids.py`
Generates share IDs for profiles missing them.

**Use case**: Migration or fixing missing share IDs
```bash
./bin/generate_profile_share_ids.py
```

#### `generate_genecart_share_ids.py` / `generate_share_ids.py`
Generates share IDs for gene carts and other entities.

**Use case**: Enabling sharing functionality

---

## Administration & Maintenance

Miscellaneous administrative and maintenance scripts.

### Data Fixes

#### `fix_array_x_axis.py` / `fix_x.py`
Fixes specific issues with array data.

**Use case**: Addressing known data issues
**Status**: May be dataset-specific, possibly obsolete

#### `fix_double_encoded_json.py`
Fixes double-encoded JSON in database.

**Use case**: Fixing import errors from legacy system

#### `nullify_na_columns.py`
Converts "NA" strings to NULL in database.

**Use case**: Database cleanup

#### `preprocess_zarr_datasets.py`
Preprocesses spatial Zarr datasets with standard single-cell workflow.

**Use case**: Automated preprocessing of spatial transcriptomics datasets
```bash
./bin/preprocess_zarr_datasets.py
```
**Note**: Runs scanpy workflow (highly variable genes, PCA, neighbors, UMAP) on all Zarr datasets in `www/datasets/spatial/`

### User & Permission Management

#### `change_gene_column_case.py`
Changes case of gene symbols (upper/lower/initcap).

**Use case**: Standardizing gene symbol formatting
```bash
./bin/change_gene_column_case.py -i input.tab -o output.tab --mode initcap
```

#### `reassign_curations_to_owner.py`
Reassigns curations to dataset owners.

**Use case**: Fixing ownership after bulk curation
**Note**: Originally created for specific conference scenario

#### `get_email_list.py`
Gets email list of users (for announcements).

**Use case**: User communication
```bash
./bin/get_email_list.py --active_only > emails.txt
```

### External Data Integration

#### `fetch_komp_repository_data.py`
Fetches data from KOMP repository.

**Use case**: Updating KOMP integration

#### `query_geo.py`
Queries GEO (Gene Expression Omnibus) database.

**Use case**: Finding datasets to import

#### `get_shield_pages.py` / `process_shield_pages.py`
Downloads and processes pages from SHIELD database.

**Use case**: Integrating SHIELD gene expression data
**Status**: May be deprecated

### Utilities

#### `calculate_averages_stddev.py`
Calculates averages and standard deviations from Excel files.

**Use case**: Statistical processing of tabular data
```bash
./bin/calculate_averages_stddev.py -i input.xlsx -c config.ini -o output.xlsx
```

#### `tab_file_info.py`
Reports information about tab-delimited files.

**Use case**: Quick file inspection
```bash
./bin/tab_file_info.py -i data.tab
```

#### `export_gene_cart_sql.py`
Exports gene cart (gene list) as SQL.

**Use case**: Backing up or migrating gene lists

#### `parse_apache_error_logs.py`
Parses Apache error logs for patterns.

**Use case**: Debugging production issues
```bash
./bin/parse_apache_error_logs.py -i /var/log/apache2/error.log > summary.txt
```

### Display Helpers

#### `generate_help_ids.py`
Generates help IDs for UI elements.

**Use case**: Maintaining help system

#### `generate_initial_composition_plots.py`
Generates initial composition plots for datasets.

**Use case**: Pre-computing common visualizations

#### `update_plotly_json_db.py`
Updates Plotly JSON configurations in database.

**Use case**: Migrating to new Plotly versions

### Data Manipulation

#### `rename_single_cell_dataset_headers.py`
Renames headers in single-cell datasets.

**Use case**: Standardizing column names

#### `find_duplicated_dataset_genes.py`
Finds duplicate genes within datasets.

**Use case**: Quality control
```bash
./bin/find_duplicated_dataset_genes.py -i dataset.h5ad
```

### Misc

#### `convert_layout_member_datasets_to_displays.py`
Converts old layout structure to new display structure.

**Use case**: Database migration
**Status**: Likely obsolete after migration complete

#### `get_profile_dataset_download_commands.py`
Generates gcloud commands to download profile datasets.

**Use case**: Local testing with production data
```bash
./bin/get_profile_dataset_download_commands.py -p profile_share_id
```

---

## Notes on Script Status

### Potentially Obsolete Scripts

These scripts may no longer be needed (verify before removal):
- `convert_layout_member_datasets_to_displays.py` - Migration script
- `load_old_gcid_gene_symbols.py` - Legacy system migration
- `get_shield_pages.py` / `process_shield_pages.py` - SHIELD integration may be deprecated
- `fix_x.py` / `fix_array_x_axis.py` - Specific data fixes
- Any script referencing "old" or "legacy" in name

### Scripts Needing Review

- Scripts in Administration & Maintenance may be one-off fixes
- Check with @adkinsrs about current status of migration scripts

## Getting Help

- Most scripts include `--help` flag for detailed usage
- Check script header comments for docstrings
- Ask in team channels about script purpose if unclear
- Create an issue on GitHub for script bugs

## Contributing

When creating new scripts:
1. Include clear docstring explaining purpose
2. Add `--help` argument with argparse
3. Follow naming conventions (verb_noun_object.py)
4. Update this documentation
5. Consider if script should be integrated into API instead
