DROP TABLE dataset_ontology_term, ontology_term;
ALTER TABLE dataset
  DROP COLUMN repository_accession,
  DROP COLUMN device,
  DROP COLUMN assay,
  DROP COLUMN library_cons_method,
  DROP COLUMN rna_extraction_method,
  DROP COLUMN microarray_hybridization_protocol,
  DROP COLUMN normalization_method,
  DROP COLUMN log_transformation,
  DROP COLUMN primary_analysis_completed,
  DROP COLUMN file_format,
  DROP COLUMN metadata_version,
  DROP COLUMN metadata_source,
  DROP COLUMN curation_flags,
  DROP COLUMN pi_name;
ALTER TABLE dataset MODIFY annotation_release VARCHAR(25);