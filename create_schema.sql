-- drop database gear_portal; create database gear_portal; use gear_portal;
-- source /home/jorvis/git/gEAR/create_schema.sql

CREATE TABLE organism (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       label          VARCHAR(255) NOT NULL,
       genus          VARCHAR(255),
       species        VARCHAR(255),
       strain         VARCHAR(255),
       taxon_id       INT
) ENGINE=INNODB;

-- DO NOT change these values without making corresponding changes in the annotation loading scripts
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (1, 'Mouse', 'Mus', 'musculus', NULL, 10090);
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (2, 'Human', 'Homo', 'sapiens', 'sapiens', 9606);
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (3, 'Zebrafish', 'Danio', 'rerio', NULL, 7955);
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (5, 'Chicken', 'Gallus', 'gallus', NULL, 9031);
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (6, 'Rat', 'Rattus', 'norvegicus', NULL, 10116);
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (7, 'Marmoset', 'Callithrix', 'jacchus', 'jacchus', 9483);
INSERT INTO organism (id, label, genus, species, strain, taxon_id)
       VALUES (8, 'Roundworm', 'Caenorhabditis', 'elegans', 'WBcel235', 6239);

-- https://docs.python.org/3.4/library/http.cookies.html
CREATE TABLE guser (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       user_name      VARCHAR(255),
       email          VARCHAR(255),
       institution    VARCHAR(255),
       pass           VARCHAR(50),
       colorblind_mode TINYINT(1) DEFAULT 0,
       updates_wanted TINYINT(1),
       is_admin       TINYINT(1) DEFAULT 0,
       help_id        VARCHAR(50),
       date_created   DATETIME DEFAULT CURRENT_TIMESTAMP,
       default_org_id INT NOT NULL DEFAULT 1,
       layout_id      INT,
       is_curator     TINYINT(1) DEFAULT 0,
       FOREIGN KEY fk_guser_doi(default_org_id) REFERENCES organism(id)
) ENGINE=INNODB;

-- password is a hashlib md5 hexdigest
INSERT INTO guser (id, user_name, email, institution, pass, updates_wanted, is_admin)
       VALUES (0, 'gEAR admin', 'admin@localhost', 'UMaryland', 'fcdf1dc2c1ef7dec3dbb1a6e2c5e3c8a', 0, 1);
INSERT INTO guser (email, user_name, institution, pass, updates_wanted, is_admin)
       VALUES('jorvis@gmail.com', 'Joshua Orvis', 'IGS', 'e81e78d854d86edc38ba45c443662aee', 0, 1);

-- Group is a reserved word, so we get gEAR Group (ggroup)
CREATE TABLE ggroup (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       creator_id     INT NOT NULL,
       label          VARCHAR(255) NOT NULL,
       FOREIGN KEY (creator_id) REFERENCES guser(id) ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE user_group_membership (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       user_id        INT NOT NULL,
       group_id       INT NOT NULL,
       FOREIGN KEY (user_id) REFERENCES guser(id) ON DELETE CASCADE,
       FOREIGN KEY (group_id) REFERENCES ggroup(id) ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE user_session (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       user_id        INT,
       session_id     VARCHAR(255),
       FOREIGN KEY (user_id)
          REFERENCES guser(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE gene (
       id               INT PRIMARY KEY AUTO_INCREMENT,
       ensembl_id       VARCHAR(40),
       ensembl_version  VARCHAR(40),
       ensembl_release  INT NOT NULL,
       genbank_acc      VARCHAR(20),
       organism_id      INT NOT NULL,
       molecule         varchar(100),
       start            int(11),
       stop             int(11),
       gene_symbol      VARCHAR(20),
       product          VARCHAR(255),
       biotype          VARCHAR(100),
       INDEX            org_idx (organism_id),
       INDEX            org_sym (organism_id, gene_symbol),
       INDEX            gene_sym (gene_symbol),
       FOREIGN KEY (organism_id) REFERENCES organism(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

-- Changes here mean bin/export_gene_cart_sql.py need to be updated too
CREATE TABLE gene_cart (
       id              INT PRIMARY KEY AUTO_INCREMENT,
       user_id         INT NOT NULL,
       organism_id     INT NOT NULL,
       gctype          VARCHAR(50) NOT NULL DEFAULT 'unweighted-list',
       label           VARCHAR(255) NOT NULL,
       ldesc           TEXT,
       share_id        VARCHAR(50),
       is_public       TINYINT DEFAULT 0,
       is_domain       TINYINT DEFAULT 0,
       date_added      DATETIME DEFAULT CURRENT_TIMESTAMP,
       FULLTEXT        text_idx (label, ldesc),
       -- INDEX           share_id (share_id),
       FOREIGN KEY (user_id) REFERENCES guser(id) ON DELETE CASCADE,
       FOREIGN KEY (organism_id) REFERENCES organism(id) ON DELETE CASCADE
) ENGINE=INNODB;

-- Changes here mean bin/export_gene_cart_sql.py need to be updated too
CREATE TABLE gene_cart_member (
       id              INT PRIMARY KEY AUTO_INCREMENT,
       gene_cart_id    INT NOT NULL,
       gene_symbol     VARCHAR(20) NOT NULL,
       FOREIGN KEY (gene_cart_id)
          REFERENCES gene_cart(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE gene_cart_group_membership (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       gene_cart_id      INT NOT NULL,
       group_id       INT NOT NULL,
       FOREIGN KEY (gene_cart_id) REFERENCES gene_cart(id) ON DELETE CASCADE,
       FOREIGN KEY (group_id) REFERENCES ggroup(id) ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE gene_symbol (
       id               INT PRIMARY KEY AUTO_INCREMENT,
       gene_id          INT NOT NULL,
       label            VARCHAR(30),
       is_primary       TINYINT(1) DEFAULT 0,
       INDEX            gene_symbol_label_idx (label),
       FOREIGN KEY (gene_id) REFERENCES gene(id)
) ENGINE=INNODB;
CREATE INDEX idx_gene_symbol__gene_id ON gene_symbol (gene_id);
CREATE INDEX idx_gene_symbol__label ON gene_symbol (label);

CREATE TABLE go (
       go_id          VARCHAR(20) PRIMARY KEY,
       name           VARCHAR(255) NOT NULL,
       namespace      VARCHAR(30) NOT NULL,
       def            TEXT
) ENGINE=INNODB;

CREATE TABLE gene_go_link (
       id             INT UNIQUE KEY AUTO_INCREMENT,
       gene_id        INT NOT NULL,
       go_id          VARCHAR(20) NOT NULL,
       PRIMARY KEY (gene_id, go_id),
       FOREIGN KEY (gene_id) REFERENCES gene(id)
) ENGINE=INNODB;

CREATE TABLE gene_dbxref (
       id             INT UNIQUE KEY AUTO_INCREMENT,
       gene_id        INT NOT NULL,
       dbxref         VARCHAR(100) NOT NULL,
       PRIMARY KEY (gene_id, dbxref),
       FOREIGN KEY (gene_id) REFERENCES gene(id)
) ENGINE=INNODB;

CREATE TABLE mirna_family (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       stem_loop_id   INT NOT NULL,
       mature_id      INT NOT NULL,
       family_id      VARCHAR(20),
       family_label   VARCHAR(20),
       FOREIGN KEY (stem_loop_id) REFERENCES gene(id),
       FOREIGN KEY (mature_id) REFERENCES gene(id)
) ENGINE=INNODB;

CREATE TABLE anatomy (
       id            INT PRIMARY KEY AUTO_INCREMENT,
       organism_id   INT NOT NULL,
       label         VARCHAR(100),
       parent_id     INT,
       FOREIGN KEY (parent_id)
          REFERENCES anatomy(id)
          ON DELETE CASCADE,
       FOREIGN KEY (organism_id)
          REFERENCES organism(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

-- Any changes to lengths here need to be also changed in form validation of upload_dataset.js
CREATE TABLE dataset (
       id                       VARCHAR(50) PRIMARY KEY,
       owner_id                 INT NOT NULL,
       title                    VARCHAR(255) NOT NULL,
       organism_id              INT NOT NULL,
       pubmed_id                VARCHAR(20),
       geo_id                   VARCHAR(50),
       is_public                TINYINT DEFAULT 0,
       is_downloadable          TINYINT DEFAULT 1,
       ldesc                    TEXT,
       date_added               DATETIME,
       dtype                    VARCHAR(50) NOT NULL DEFAULT 'svg-expression',
       -- paths are relative to the root, so probably like datasets_uploaded/{dataset_id}.jpg
       schematic_image          VARCHAR(255),
       share_id                 VARCHAR(50),
       math_default             VARCHAR(50) NOT NULL DEFAULT 'raw', -- options: 'raw', 'log2', 'log10'
       marked_for_removal       TINYINT DEFAULT 0,
       load_status              VARCHAR(20), -- options: 'pending', 'loading', 'completed', 'failed',
       has_h5ad                 TINYINT DEFAULT 0,
       platform_id              VARCHAR(255),
       instrument_model         VARCHAR(255),
       library_selection        VARCHAR(255),
       library_source           VARCHAR(255),
       library_strategy         TEXT,
       contact_email            VARCHAR(100),
       contact_institute        VARCHAR(255),
       contact_name             VARCHAR(100),
       annotation_source        VARCHAR(20),
       plot_default             VARCHAR(50), -- options: 'bar', 'line', 'violin'
       annotation_release       INT,
       FULLTEXT                 text_idx (title, ldesc),
       FULLTEXT                 text_with_geo_idx (title, ldesc, geo_id),
       FULLTEXT                 text_with_geo_pubmed_idx (title, ldesc, geo_id, pubmed_id),
       FOREIGN KEY (owner_id)
          REFERENCES guser(id)
          ON DELETE CASCADE,
       FOREIGN KEY (organism_id)
          REFERENCES organism(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE dataset_display (
       id            INT PRIMARY KEY AUTO_INCREMENT,
       dataset_id    VARCHAR(50) NOT NULL,
       user_id       INT NOT NULL,
       label         VARCHAR(255),
       plot_type     VARCHAR(20),
       plotly_config TEXT,

       FOREIGN KEY (dataset_id)
          REFERENCES dataset(id)
          ON DELETE CASCADE,
       FOREIGN KEY (user_id)
          REFERENCES guser(id)
          ON DELETE CASCADE,
       INDEX user_dataset (user_id, dataset_id)
) ENGINE=INNODB;

CREATE TABLE dataset_preference (
   user_id        INT NOT NULL,
   dataset_id     VARCHAR(50) NOT NULL,
   display_id     INT NOT NULL,
   is_multigene   TINYINT(1) DEFAULT 0,
   primary key (user_id, dataset_id, is_multigene),

   FOREIGN KEY (dataset_id)
      REFERENCES dataset(id)
      ON DELETE CASCADE,
   FOREIGN KEY (user_id)
      REFERENCES guser(id)
      ON DELETE CASCADE,
   FOREIGN KEY (display_id)
      REFERENCES dataset_display(id)
      ON DELETE CASCADE
) ENGINE=INNODB;

-- Stores custom external URLs to be displayed with each dataset
CREATE TABLE dataset_link (
      id                        INT PRIMARY KEY AUTO_INCREMENT,
      dataset_id                VARCHAR(50) NOT NULL,
      resource                  VARCHAR(100) NOT NULL,
      label                     VARCHAR(100) NOT NULL,
      url                       VARCHAR(255) NOT NULL
);

CREATE TABLE dataset_shares (
      id                        INT PRIMARY KEY AUTO_INCREMENT,
      dataset_id                VARCHAR(50) NOT NULL,
      user_id                   INT NOT NULL,
      is_allowed                TINYINT(1) DEFAULT 0,
      FOREIGN KEY (dataset_id)
          REFERENCES dataset(id)
          ON DELETE CASCADE,
      FOREIGN KEY (user_id)
          REFERENCES guser(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE event (
      id                        INT PRIMARY KEY AUTO_INCREMENT,
      label                     VARCHAR(255) NOT NULL,
      max_attendees             INT NOT NULL,
      waitlist_size             INT NOT NULL DEFAULT 0,
      date_added                DATETIME DEFAULT CURRENT_TIMESTAMP
);

INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Introduction I', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Explore and analyze basics', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Explore and customize I', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Explore and customize II', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Introduction redo', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Analyze scRNA-seq data', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Transfer learning', 50, 10);
INSERT INTO event (label, max_attendees, waitlist_size) VALUES ('ARO 2023 - Data upload', 50, 10);

CREATE TABLE event_registration (
      id                        INT PRIMARY KEY AUTO_INCREMENT,
      event_id                  INT NOT NULL,
      user_id                   INT NOT NULL,
      date_added                DATETIME DEFAULT CURRENT_TIMESTAMP,
      FOREIGN KEY (user_id) REFERENCES guser(id) ON DELETE CASCADE,
      FOREIGN KEY (event_id) REFERENCES event(id) ON DELETE CASCADE
);

-- Recursive organizational table allowing item types (like gene carts
--  or profiles) to be grouped into 'folders'
CREATE TABLE folder (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       parent_id                INT,
       label                    VARCHAR(100) NOT NULL,
       FOREIGN KEY (parent_id) REFERENCES folder(id) ON DELETE CASCADE
) ENGINE=INNODB;

-- The label for this one is not actually displayed.  It is set in tree.js
INSERT INTO folder (id, parent_id, label) VALUES (101, NULL, 'Highlighted profiles');
INSERT INTO folder (id, parent_id, label) VALUES (102, NULL, 'Your profiles');
INSERT INTO folder (id, parent_id, label) VALUES (103, NULL, 'Group profiles');
INSERT INTO folder (id, parent_id, label) VALUES (104, NULL, 'Profiles shared with you');
INSERT INTO folder (id, parent_id, label) VALUES (105, NULL, 'Other public profiles');
INSERT INTO folder (id, parent_id, label) VALUES (106, NULL, 'Highlighted gene carts');
INSERT INTO folder (id, parent_id, label) VALUES (107, NULL, 'Your gene carts');
INSERT INTO folder (id, parent_id, label) VALUES (108, NULL, 'Group gene carts');
INSERT INTO folder (id, parent_id, label) VALUES (109, NULL, 'Gene carts shared with you');
INSERT INTO folder (id, parent_id, label) VALUES (110, NULL, 'Other public carts');

CREATE TABLE folder_member (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       folder_id                INT NOT NULL,
       item_id                  INT NOT NULL,
       item_type                VARCHAR(20) NOT NULL, -- options: 'layout', 'genecart'
       FOREIGN KEY (folder_id) REFERENCES folder(id) ON DELETE CASCADE,
       UNIQUE KEY uk_folder_item (folder_id, item_id, item_type)
) ENGINE=INNODB;

CREATE TABLE layout (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       user_id                  INT NOT NULL,
       label                    VARCHAR(255),
       is_current               TINYINT(1) DEFAULT 0,
       is_domain                TINYINT(1) DEFAULT 0,
       is_public                TINYINT(1) DEFAULT 0,
       share_id                 VARCHAR(24),
       CONSTRAINT idx_layout_share_id UNIQUE (share_id),
       FOREIGN KEY (user_id)
          REFERENCES guser(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

-- Adding this here so we miss a chicken/egg problem, since both reference each other
ALTER TABLE guser ADD CONSTRAINT FOREIGN KEY fk_guser_layout(layout_id) REFERENCES layout(id) ON DELETE CASCADE;

INSERT INTO layout VALUES (0, 0, NULL, "Hearing (default)", 1);
INSERT INTO layout VALUES (10000, 0, NULL, "Brain development (default)", 0);
INSERT INTO layout VALUES (10001, 0, NULL, "Huntingtons disease (default)", 0);

/* Soom to delete */
CREATE TABLE layout_members (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       layout_id                INT NOT NULL,
       dataset_id               VARCHAR(50) NOT NULL,
       grid_position            INT NOT NULL,
       mg_grid_position         INT NOT NULL,
       start_col                INT NOT NULL DEFAULT 1,
       mg_start_col             INT NOT NULL DEFAULT 1,
       grid_width               INT NOT NULL DEFAULT 4,
       mg_grid_width            INT NOT NULL DEFAULT 12,
       start_row                INT NOT NULL DEFAULT 1,
       mg_start_row             INT NOT NULL DEFAULT 1,
       grid_height              INT NOT NULL DEFAULT 1, -- height is number of rows spanned, which is not based on a grid
       mg_grid_height           INT NOT NULL DEFAULT 1,
       math_preference          VARCHAR(50), -- options: 'raw', 'log2', 'log10'
       plot_preference          VARCHAR(50), -- options: 'bar', 'line', 'violin'
       FOREIGN KEY (layout_id)
          REFERENCES layout(id)
          ON DELETE CASCADE,
       FOREIGN KEY (dataset_id)
          REFERENCES dataset(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

/* Soon to change to layout_members */
CREATE TABLE layout_displays (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       layout_id                INT NOT NULL,
       display_id               INT NOT NULL,
       grid_position            INT NOT NULL,
       start_col                INT NOT NULL DEFAULT 1,
       grid_width               INT NOT NULL DEFAULT 4,
       start_row                INT NOT NULL DEFAULT 1,
       grid_height              INT NOT NULL DEFAULT 1, -- height is number of rows spanned, which is not based on a grid
       math_preference          VARCHAR(50), -- options: 'raw', 'log2', 'log10'
       FOREIGN KEY (layout_id)
          REFERENCES layout(id)
          ON DELETE CASCADE,
       FOREIGN KEY (display_id)
          REFERENCES dataset_display(id)
          ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE layout_group_membership (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       layout_id      INT NOT NULL,
       group_id       INT NOT NULL,
       FOREIGN KEY (layout_id) REFERENCES layout(id) ON DELETE CASCADE,
       FOREIGN KEY (group_id) REFERENCES ggroup(id) ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE supplemental_images (
       id                       INT PRIMARY KEY AUTO_INCREMENT,
       -- label is a meant to identify the class of image
       label                    VARCHAR(50),
       ensembl_id               VARCHAR(20),
       gene_symbol              VARCHAR(20),
       image_url                VARCHAR(200),
       index                    gene_sym_idx(gene_symbol)
) ENGINE=INNODB;

CREATE TABLE tag (
	id		INT PRIMARY KEY AUTO_INCREMENT,
	label	VARCHAR(55)
) ENGINE=INNODB;

CREATE TABLE comment (
    id int PRIMARY KEY AUTO_INCREMENT,
    first_name varchar(255) DEFAULT NULL,
    last_name varchar(255) DEFAULT NULL,
    user_id int NOT NULL,
    email varchar(255) DEFAULT NULL,
    title varchar(255) DEFAULT NULL,
    message varchar(1020) DEFAULT NULL,
    is_read tinyint DEFAULT 0,
    date_added datetime DEFAULT NULL,
    FOREIGN KEY comment_ibfk_1 (user_id) REFERENCES guser(id)
) ENGINE=INNODB;

-- multiple tags to multiple comments
CREATE TABLE comment_tag (
	id          INT PRIMARY KEY AUTO_INCREMENT,
	tag_id      INT,
	comment_id  INT,
	FOREIGN KEY (tag_id) REFERENCES tag(id),
	FOREIGN KEY (comment_id) REFERENCES comment(id) ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE dataset_tag (
    id                   INT PRIMARY KEY AUTO_INCREMENT,
    tag_id               INT,
    dataset_id           VARCHAR(50) NOT NULL,
    FOREIGN KEY (tag_id)
      REFERENCES tag(id),
    FOREIGN KEY (dataset_id)
      REFERENCES dataset(id)
      ON DELETE CASCADE
) ENGINE=INNODB;

CREATE TABLE note (
  id                  INT PRIMARY KEY AUTO_INCREMENT,
  title               VARCHAR(100),
  ldesc               TEXT,
  user_id             INT,
  dataset_id          VARCHAR(50),
  is_public           TINYINT DEFAULT 0,
  date_added          DATETIME,
  date_last_change    DATETIME,
  FOREIGN KEY (user_id)
     REFERENCES guser(id),
  FOREIGN KEY (dataset_id)
     REFERENCES dataset(id)
) ENGINE=INNODB;

CREATE TABLE `dataset_epiviz` (
  `id` varchar(50) NOT NULL,
  `owner_id` int(11) NOT NULL,
  `annotation` text,
  `type` varchar(10) NOT NULL,
  `url` text NOT NULL,
  `title` text,
  `is_public` tinyint(4) DEFAULT '0',
  `description` text,
  `share_id` varchar(50) NOT NULL,
  `organism` varchar(100) DEFAULT NULL
) ENGINE=InnoDB;

CREATE TABLE user_history (
    id			    INT PRIMARY KEY AUTO_INCREMENT,
    user_id         INT NOT NULL,
    entry_date      DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
    entry_category  VARCHAR(100) NOT NULL,
    label           VARCHAR(255) NOT NULL,
    url             VARCHAR(255) NOT NULL,
    FOREIGN KEY (user_id)
     REFERENCES guser(id)
) ENGINE=InnoDB;

/* Restrictd datasets can only be accessed by users in a specific group */
/* NOTE: Not inserting at this time, as I am currently using "dataset_shares"
       to address a shared dataset
CREATE TABLE dataset_group_membership (
       id             INT PRIMARY KEY AUTO_INCREMENT,
       dataset_id     VARCHAR(50) NOT NULL,
       group_id       INT NOT NULL,
       FOREIGN KEY (dataset_id) REFERENCES dataset(id) ON DELETE CASCADE,
       FOREIGN KEY (group_id) REFERENCES ggroup(id) ON DELETE CASCADE
) ENGINE=INNODB;
*/

CREATE TABLE submission (
       id               VARCHAR(50) PRIMARY KEY,
       user_id                     INT NOT NULL,
       layout_id                   INT,
       is_finished                 TINYINT DEFAULT 0,
       is_restricted               TINYINT DEFAULT 0, /* if one dataset is restricted, then the whole submission must be */
       date_added                DATETIME DEFAULT CURRENT_TIMESTAMP,
       email_updates               TINYINT DEFAULT 0,
       FOREIGN KEY (user_id) REFERENCES guser(id),
       FOREIGN KEY (layout_id) REFERENCES layout(id)
) ENGINE=INNODB;

CREATE TABLE submission_dataset (
       id                          INT PRIMARY KEY AUTO_INCREMENT,
       dataset_id                  VARCHAR(50) NOT NULL,
       nemo_identifier             VARCHAR(20) NOT NULL, /* from nemoarchive (should we do UUID here and make new one?) */
       pulled_to_vm_status         VARCHAR(20) default "pending", /*options: 'pending', 'loading', 'completed', 'canceled', 'failed',*/
       convert_metadata_status   VARCHAR(20) default "pending", /*options: 'pending', 'loading', 'completed', 'canceled', 'failed',*/
       convert_to_h5ad_status      VARCHAR(20) default "pending", /*options: 'pending', 'loading', 'completed', 'canceled', 'failed',*/
       make_tsne_status      VARCHAR(20) default "pending", /*options: 'pending', 'loading', 'completed', 'canceled', 'failed',*/
       log_message                 TEXT,
       is_restricted               TINYINT DEFAULT 0,
       FOREIGN KEY (dataset_id) REFERENCES dataset(id) ON DELETE CASCADE
) ENGINE=INNODB;
/* For some reason the collation of "dataset" table is "latin1_swedish_ci", and modifying it to "latin1_general_ci" would need to cascade elsewhere
So for now I am just removing the COLLATE part of the ENGINE syntax */

CREATE TABLE submission_member (
       id                          INT PRIMARY KEY AUTO_INCREMENT,
       submission_id               VARCHAR(50) NOT NULL,
       submission_dataset_id       INT NOT NULL,
       FOREIGN KEY (submission_id) REFERENCES submission(id) ON DELETE CASCADE,
       FOREIGN KEY (submission_dataset_id) REFERENCES submission_dataset(id) ON DELETE CASCADE
) ENGINE=INNODB;
