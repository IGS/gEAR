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
    ADD COLUMN curation_flags TEXT,
    ADD COLUMN pi_name VARCHAR(100),
    ADD INDEX idx_pi_name (pi_name);

-- Term labels can hold Greek letters and other characters latin1 can't store
CREATE TABLE ontology_term (
       term_id                  VARCHAR(30) PRIMARY KEY,
       ontology                 VARCHAR(20) NOT NULL,
       label                    VARCHAR(255),
       date_fetched             DATETIME DEFAULT CURRENT_TIMESTAMP
) ENGINE=INNODB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

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

-- Store any character users type or paste into the metadata, such as the β in "β cell" (GitHub #59, #633).
--  The dataset table is latin1, which can't, so the upload failed when it was saved, after all the processing.
--  Only the columns holding user text are converted. dataset.id stays latin1 because the foreign keys of
--  other (latin1) tables point at it. TEXT columns become MEDIUMTEXT, since a utf8mb4 character can take
--  up to 4 bytes and TEXT holds 65535 bytes. The columns of each FULLTEXT index (title, ldesc, geo_id,
--  pubmed_id) must share a character set, so they are converted in one statement. Safe to run on servers
--  where #633 already converted title.
ALTER TABLE dataset
    MODIFY title VARCHAR(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
    MODIFY ldesc MEDIUMTEXT CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY geo_id VARCHAR(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY pubmed_id VARCHAR(20) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY platform_id VARCHAR(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY instrument_model VARCHAR(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY library_selection VARCHAR(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY library_source VARCHAR(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY library_strategy MEDIUMTEXT CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY contact_email VARCHAR(100) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY contact_institute VARCHAR(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY contact_name VARCHAR(100) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY library_cons_method MEDIUMTEXT CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY rna_extraction_method MEDIUMTEXT CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY microarray_hybridization_protocol MEDIUMTEXT CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY normalization_method VARCHAR(100) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci,
    MODIFY pi_name VARCHAR(100) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;

-- Keywords are stored as tags
ALTER TABLE tag MODIFY label VARCHAR(55) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
