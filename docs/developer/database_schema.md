# Database Schema

gEAR stores users, dataset metadata, curated displays, collections and gene lists in MySQL; expression matrices themselves live on disk as H5AD/Zarr files (see [upload pipeline](./upload_pipeline.md#on-disk-layout)). This page summarizes every `CREATE TABLE` in [`create_schema.sql`](../../create_schema.sql), grouped by domain.

Connection settings come from the `[database]` section of `gear.ini` (see [configuration](./configuration.md)); Python code connects through `geardb.Connection` in `lib/geardb.py` or `gear.db.MySQLDB` in `lib/gear/db.py`. For installation see [MySQL setup](./setup/mysql.md) and [Docker MySQL](./setup/docker_mysql.md).

## Terminology

The schema predates the current UI wording. Code and tables use the old names:

| UI term | Table / code term | Notes |
|---------|-------------------|-------|
| Dataset collection | `layout` (also "profile" in older code and in `folder` labels, e.g. `default_profile_share_id`, `[folders] profile_*`) | `geardb.Layout`, `geardb.LayoutCollection` |
| Gene list | `gene_cart` ("gene cart", "genecart") | `geardb.GeneCart`, `geardb.GeneCartCollection` |
| Display / curation | `dataset_display` | A saved plot configuration for one dataset |
| Group | `ggroup` | "group" is a MySQL reserved word |
| User | `guser` | `geardb.User` |

## Users and authentication

| Table | Purpose |
|-------|---------|
| `organism` | Supported organisms (Mouse=1, Human=2, Zebrafish=3, Chicken=5, Rat=6, Marmoset=7, Roundworm=8). IDs are hard-coded in annotation loaders; do not renumber. |
| `guser` | User accounts: name, email, institution, password hash (`pass`), `is_admin`, `is_curator`, `colorblind_mode`, `default_org_id` -> `organism`, `layout_id` -> `layout` (user's current collection), `help_id`. |
| `user_session` | Session IDs issued at login (`www/cgi/login.v2.cgi`); `user_id` -> `guser`. Looked up by `geardb.get_user_from_session_id()`. |
| `ggroup` | User-created groups; `creator_id` -> `guser`. |
| `user_group_membership` | User <-> group membership. |
| `event` | Workshop/event definitions with capacity and waitlist size (seeded with ARO 2023 sessions). |
| `event_registration` | User registrations for events (`www/cgi/set_user_event_registration.cgi`, `get_event_registration_list.cgi`). |

## Datasets

| Table | Purpose |
|-------|---------|
| `dataset` | One row per dataset. `id` is a UUID string that also names the files on disk (`www/datasets/<id>.h5ad`). Key columns: `owner_id` -> `guser`, `organism_id` -> `organism`, `title`, `ldesc`, `dtype` (values used by the uploader: `single-cell-rnaseq`, `bulk-rnaseq`, `microarray`, `spatial`, `gosling`; legacy default `svg-expression`), `share_id`, `is_public`, `is_downloadable`, `load_status` (`pending`/`loading`/`completed`/`failed`), `has_h5ad`, `marked_for_removal`, GEO/PubMed IDs, sequencing metadata, contact fields, `annotation_source`/`annotation_release`, `user_pii_affirmed`. FULLTEXT indexes on title/ldesc(/geo_id/pubmed_id) back dataset search. Column lengths are mirrored in `www/js/upload_dataset.js` validation. |
| `dataset_shares` | Per-user access grants to private datasets (`is_allowed`). |
| `dataset_link` | Custom external URLs shown with a dataset (`resource`, `label`, `url`); no foreign key declared. `geardb.DatasetLink`. |
| `dataset_tag` | Dataset <-> `tag` link. |
| `note` | User notes attached to a dataset (`title`, `ldesc`, `is_public`). |
| `dataset_epiviz` | Legacy Epiviz track metadata. No longer written by application code (Gosling replaced Epiviz; see `add_gosling_display_curation()` in `lib/geardb.py`). |
| `supplemental_images` | Gene-level image URLs by `gene_symbol`/`ensembl_id`. Not referenced by current code. |
| `anatomy` | Organism-specific anatomy hierarchy (`parent_id` self-reference). Not referenced by current application code. |

`dataset_group_membership` is present in the file but commented out; dataset sharing uses `dataset_shares` instead.

## Displays and collections (layouts)

| Table | Purpose |
|-------|---------|
| `dataset_display` | Saved plot configuration ("curation") for a dataset: `plot_type`, `label`, `plotly_config` (JSON text of plot options). `dataset_id` -> `dataset`, `user_id` -> `guser`. `geardb.DatasetDisplay`. |
| `dataset_preference` | A user's default display per dataset, separately for single-gene and multigene (`is_multigene`); PK (`user_id`, `dataset_id`, `is_multigene`), `display_id` -> `dataset_display`. |
| `layout` | A dataset collection: `user_id` owner, `label`, `share_id` (unique), `is_public`, `is_domain` (site-wide default collections), `is_current`. Seeded with IDs 0, 10000, 10001 owned by admin user 0. `geardb.Layout`. |
| `layout_displays` | Current collection membership: which `dataset_display` appears in which `layout`, with grid placement (`grid_position`, `start_col`, `grid_width`, `start_row`, `grid_height`) and `math_preference`. `geardb.LayoutDisplay`. |
| `layout_members` | Older membership model (dataset-level, separate single-gene `grid_*` and multigene `mg_*` placement). Marked "soon to delete" in the schema but still read/written by `geardb.LayoutMember` and `www/cgi/get_users_layout_members.cgi`. `bin/convert_layout_member_datasets_to_displays.py` migrates rows to `layout_displays`. |
| `layout_group_membership` | Shares a collection with a `ggroup`. |

## Gene lists (gene carts)

| Table | Purpose |
|-------|---------|
| `gene_cart` | A gene list: `user_id`, `organism_id`, `gctype` (`unweighted-list` default; also `weighted-list`, etc.), `label`, `ldesc`, `share_id`, `is_public`, `is_domain`. FULLTEXT on label/ldesc. Unweighted lists keep their genes in `gene_cart_member`; weighted lists (`gctype = 'weighted-list'`) store genes and weights in a file on disk (`www/carts/cart.<share_id>.tab`; see [upload pipeline](./upload_pipeline.md#on-disk-layout)). Changes must be mirrored in `bin/export_gene_cart_sql.py`. |
| `gene_cart_member` | Gene symbols belonging to an unweighted cart. |
| `gene_cart_group_membership` | Shares a cart with a `ggroup`. |

## Organization (folders)

| Table | Purpose |
|-------|---------|
| `folder` | Recursive folder tree (`parent_id` self-reference) used to organize collections and gene lists. Root folders 101-110 are seeded ("Highlighted profiles", "Your gene carts", ...); `[folders]` in `gear.ini` points at some of them. |
| `folder_member` | Places an item in a folder: `item_id` + `item_type` (`layout` or `genecart`); unique per (folder, item, type). |

## Genes, annotation and orthology

| Table | Purpose |
|-------|---------|
| `gene` | Annotated genes per organism and Ensembl release: `ensembl_id`, `ensembl_version`, `ensembl_release`, `genbank_acc`, `gene_symbol`, `product`, `biotype`, coordinates (`molecule`, `start`, `stop`). Loaded by `bin/load_ensembl_gbk_annotations.py` and related loaders (see [scripts README](../misc/scripts/README.md)). `geardb.Gene`. |
| `gene_symbol` | Primary and alternate symbols (synonyms) per gene (`is_primary`). |
| `go` | Gene Ontology terms (`go_id`, `name`, `namespace`, `def`); loaded by `bin/load_gene_ontology.py`. |
| `gene_go_link` | Gene <-> GO term annotations. |
| `gene_dbxref` | Cross-references to external databases per gene. |
| `mirna_family` | miRNA stem-loop / mature pairs (`stem_loop_id`, `mature_id` -> `gene`) and family IDs; loaded by `bin/load_mirna_fam_data.py`. |

Orthology mapping between organisms is file-based (`lib/gear/orthology.py` reads `orthomap.*.hdf5` files from `www/feature_mapping/`; source inputs are under `orthomap_inputs/`), not stored in these tables.

## Comments, history and submissions

| Table | Purpose |
|-------|---------|
| `tag` | Free-text tags shared by comments and datasets. |
| `comment` | Contact / feedback messages (`title`, `message`, `is_read`, submitter name/email); `user_id` -> `guser`. Read by `www/cgi/load_comment.cgi`. |
| `comment_tag` | Comment <-> tag link. |
| `user_history` | Recent-activity log per user (`entry_category`, `label`, `url`), written by `lib/gear/userhistory.py`. |
| `submission` | NeMO Archive import batch: `user_id`, optional target `layout_id`, `is_finished`, `is_restricted`, `email_updates`. |
| `submission_dataset` | Per-dataset import status for a submission: `nemo_identifier`, step statuses (`pulled_to_vm_status`, `convert_metadata_status`, `convert_to_h5ad_status`, `make_tsne_status`), `log_message`. |
| `submission_member` | Submission <-> submission_dataset link. |

The three `submission*` tables are not referenced by current code in `lib/`, `www/cgi/` or `www/api/`.

## Key relationships

```
organism 1--* gene 1--* gene_symbol / gene_go_link / gene_dbxref
guser 1--* user_session
guser 1--* dataset 1--* dataset_display 1--* layout_displays *--1 layout *--1 guser
                   \--* dataset_shares, dataset_tag, note, dataset_preference
guser 1--* gene_cart 1--* gene_cart_member
ggroup *--* guser (user_group_membership), layout (layout_group_membership), gene_cart (gene_cart_group_membership)
folder 1--* folder_member --> layout | gene_cart (by item_type)
guser.layout_id --> layout   (added by ALTER TABLE after layout is created)
```

Most child tables use `ON DELETE CASCADE`, so deleting a `guser` removes their datasets, displays, collections and gene lists. Exceptions without cascade include `note`, `comment`, `user_history`, `submission`, `gene_symbol`, `gene_go_link` and `gene_dbxref`.

## Schema drift

- `gene_urls` (columns `gene_id`, `label`, `url`) is inserted into by `lib/loaderutils.py`, `bin/load_zfin_id_urls.py` and `bin/load_mirna_fam_data.py` but has no `CREATE TABLE` in `create_schema.sql`.
- Several `bin/` scripts (`remove_marked_datasets.py`, `rescore_*_dataset_coloring.py`, `convert_legacy_dataset_to_3tab.py`) query a legacy `expression` table that is not in the schema.
- The seed `guser` rows use MD5 hashes while the column comment mentions SHA3-256.

## Related documentation

- [Code map](./code_map.md) - `lib/geardb.py` classes that wrap these tables
- [API reference](./api_reference.md)
- [Configuration](./configuration.md)

---

Last updated: September 2026
