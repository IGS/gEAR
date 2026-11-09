-- Emergency rollback of bin/migrations/390_dataset_metadata_v2.sql.
-- The utf8mb4 columns are left as they are: the code from before 390 works with them, and converting
--  back to latin1 would fail (or lose characters) once any text latin1 can't store has been saved.
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