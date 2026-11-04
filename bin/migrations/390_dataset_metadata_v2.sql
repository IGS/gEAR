-- Schema changes for the v2 dataset metadata template (GitHub issue #390, #421, #422).
-- Brings an existing gear_portal database in line with create_schema.sql.
--
-- Usage: mysql -u <user> -p gear_portal < bin/migrations/390_dataset_metadata_v2.sql
--
-- Only adds columns/tables (plus widening annotation_release), so existing code keeps working.
-- Columns are appended in the same order as create_schema.sql so that every server ends up with
--  the same column order (bin/migrate_datasets_across_vms.sh copies rows with SELECT *).

-- GenBank/RefSeq assembly accessions (e.g. GCF_000001635.27) are not integers
ALTER TABLE dataset MODIFY annotation_release VARCHAR(50);

ALTER TABLE dataset
    ADD COLUMN repository_accession VARCHAR(100),
    ADD COLUMN device VARCHAR(30),
    ADD COLUMN assay VARCHAR(30),
    ADD COLUMN library_cons_method TEXT,
    ADD COLUMN rna_extraction_method TEXT,
    ADD COLUMN microarray_hybridization_protocol TEXT,
    ADD COLUMN normalization_method VARCHAR(100),
    ADD COLUMN log_transformation VARCHAR(30),
    ADD COLUMN primary_analysis_completed TINYINT,
    ADD COLUMN file_format VARCHAR(30),
    ADD COLUMN metadata_version VARCHAR(10),
    ADD COLUMN metadata_source VARCHAR(30),
    ADD COLUMN curation_flags TEXT;

CREATE TABLE ontology_term (
       term_id                  VARCHAR(30) PRIMARY KEY,
       ontology                 VARCHAR(20) NOT NULL,
       label                    VARCHAR(255),
       date_fetched             DATETIME DEFAULT CURRENT_TIMESTAMP
) ENGINE=INNODB;

-- Foreign keys need matching character sets, and the existing dataset table is latin1
CREATE TABLE dataset_ontology_term (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       dataset_id               VARCHAR(50) NOT NULL,
       field                    VARCHAR(30) NOT NULL,
       term_id                  VARCHAR(30) NOT NULL,
       UNIQUE KEY uk_dataset_field_term (dataset_id, field, term_id),
       INDEX idx_term (term_id),
       FOREIGN KEY (dataset_id)
          REFERENCES dataset(id)
          ON DELETE CASCADE
) ENGINE=INNODB DEFAULT CHARSET=latin1;
