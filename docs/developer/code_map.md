# gEAR Code Map

This page maps the gEAR source tree to what each part does. It was put together by reading the code, so every path below exists in the repository. Use it to find the file behind a feature, then see the linked topic pages for detail.

Related pages: [API reference](./api_reference.md), [upload pipeline](./upload_pipeline.md), [configuration](./configuration.md), [database schema](./database_schema.md), [testing](./testing.md), [services](./services/README.md), [setup](./setup/README.md), [maintenance scripts](../misc/scripts/README.md).

## Top-level layout

| Path | Contents |
|------|----------|
| `lib/` | Shared Python library, imported by the CGIs, the Flask API, the consumers and `bin/` scripts |
| `www/` | Apache web root: HTML pages, `js/`, `css/`, `include/` partials, `cgi/`, `api/` (Flask), `plugins/` |
| `listeners/` | RabbitMQ consumer daemons and their Dockerfiles |
| `services/` | Separate services: `projectr/`, `spatial/`, `nemoanalytics-import-builder/` |
| `systemd/` | Unit files for the consumers and the spatial Panel server |
| `docker/` | Docker and Compose files for local/dev deployments, Apache config, sample SQL dumps |
| `bin/` | Maintenance, loader and one-off migration scripts |
| `tests/` | UI tests: Mocha + Playwright and pytest + SeleniumBase suites, plus older legacy runners (see [testing](./testing.md)) |
| `gear.ini.template` | Template for the server config that `gear.serverconfig.ServerConfig` reads (see [configuration](./configuration.md)) |
| `create_schema.sql` | MySQL schema (see [database schema](./database_schema.md)) |

---

## Backend

### `lib/` top-level modules

**`lib/geardb.py`** is the main data-access layer: about 4,300 lines of MySQL-backed model classes and lookup helpers. It adds `lib/` to `sys.path` and uses `gear.db.MySQLDB` for connections and `gear.serverconfig.ServerConfig` for settings.
- **Connection:** `Connection` wraps a `mysql.connector` connection and provides `get_cursor(use_dict=False)`, `commit()` and `close()`.
- **Model classes:** `Organism`, `Layout` (a dataset collection), `LayoutMember`, `LayoutDisplay`, `Folder`, `Dataset`, `DatasetDisplay`, `DatasetLink`, `Gene`, `GeneCart` (a gene list) and `User`.
- **Collection classes:** `OrganismCollection`, `LayoutCollection`, `FolderCollection`, `DatasetCollection`, `GeneCollection` and `GeneCartCollection`.
- **Module-level lookups:**
  - `get_user_from_session_id()` and `get_user_id_from_session_id()`
  - `get_dataset_by_id()`, `get_dataset_by_share_id()`, `get_dataset_by_title()` and `get_dataset_id_from_share_id()`
  - `get_layout_by_share_id()`, `get_gene_cart_by_share_id()` and `get_display_by_id()`
  - `get_default_display()` and `get_gene_by_gene_symbol()`
  - `add_spatial_panel_curation()` and `add_gosling_display_curation()`
- **Site-domain readers:** `_read_site_domain_config()` and related helpers read `www/site_domain_prefs.json`.

**`lib/gearqueue.py`** is the RabbitMQ wrapper that application code uses.
- `Connection` works as a context manager for publishers and consumers. `AsyncConnection` subclasses it for asynchronous consumers.
- Publishers:
  - `www/cgi/process_uploaded_expression_dataset.cgi` publishes to the `anndata_upload_jobs` and `spatial_upload_jobs` queues.
  - `www/api/resources/projectr.py` publishes to `projectr`.
  - `www/api/resources/track_hub.py` publishes to `trackhub_copy_jobs`.
- Consumers: all four `listeners/*_consumer.py` scripts.

**`lib/loaderutils.py`** holds cursor-level helpers for the annotation loader scripts, such as `add_gene`, `add_gene_symbol`, `add_gene_url` and the `cache_genes_by_*` / `cache_gene_aliases` lookup builders. The `bin/load_*` scripts and `bin/update_gene_coordinates_from_genbank.py` import it, and so does `www/cgi/get_tag_list.cgi`.

**`lib/setup.py`** is a minimal setuptools file that installs `lib/` as the `gear` package. Some Docker images use it so that code can run without appending `lib` to `PYTHONPATH`.

### `lib/gear/` package

| Module | Purpose / main contents | Main users |
|--------|------------------------|------------|
| `analysis.py` | `Analysis` and `SpatialAnalysis` (path conventions for primary/public/user analyses), `AnalysisCollection`, `H5adAdapter` / `ZarrAdapter` (one interface over `.h5ad` and `.zarr`), and factories `get_analysis()` and `get_primary_analysis()` | many analysis CGIs (`h5ad_generate_*`, `get_stored_analysis.cgi`, `save_dataset_analysis.cgi`, ...), API resources |
| `anndata_processor.py` | `AnndataProcessor` converts uploads (H5AD, 3-tab, Excel, MEX, Seurat RDS) to H5AD. `process_anndata_synchronously()` and `write_status()` | `process_uploaded_expression_dataset.cgi`, `listeners/anndata_upload_consumer.py`, `spatial_processor.py` |
| `cosmx_reader.py` | `read_cosmx()` reads CosMx data in chunks, adapted from `spatialdata_io` so large expression matrices do not run out of memory | `spatialhandler.py` (`CosMxHandler`) |
| `db.py` | `MySQLDB.connect()` builds a `mysql.connector` connection from the `[database]` section of `gear.ini` | `geardb.Connection`, a few `bin/` scripts |
| `fromgeo.py` | `FromGeo.get_geo_data()` fetches GEO (GSE/GSM) metadata | `metadata.py`, `get_metadata_from_geo.cgi` |
| `metadata.py` | `Metadata` reads an uploaded metadata XLSX/JSON, validates it and saves it to MySQL | `upload_expression_metadata.cgi`, `finalize_uploaded_expression_dataset.cgi`, `bin/add_excel_metadata_to_db.py` |
| `metadatavalidator.py` | `MetadataValidator`: lists of required attributes and per-field validators | `metadata.py` only |
| `mg_plotting.py` | Multi-gene Plotly builders: dot plot, heatmap (with dendrograms and cluster bars), quadrant, (stacked) violin and volcano, plus the matching `prep_*` / `validate_*` helpers | `www/api/resources/mg_plotly_data.py` |
| `orthology.py` | Maps genes between organisms using the HDF5 ortholog files (`map_single_gene`, `map_multiple_genes`, `map_dataframe_genes`, `get_ortholog_file`) | `resources/orthologs.py`, `resources/projectr.py` |
| `plotting.py` | Single-gene Plotly figures: `generate_plot()` (bar, box, histogram, line, scatter, strip, contour and violin) and `plotly_color_map()` | `resources/plotly_data.py`, `svg_data.py`, `tsne_data.py`, `spatial_scanpy_data.py`, `common.py`, `get_h5ad_plot.cgi` |
| `primary_analysis.py` | `add_primary_analysis_to_dataset()` detects or adds clustering, tSNE and UMAP, and creates composition plots | `anndata_processor.py` |
| `queue.py` | Low-level pika helpers: `get_connection()`, `get_async_connection()` and `RabbitMQQueue` | `lib/gearqueue.py` only (see the overlap note below) |
| `serverconfig.py` | `ServerConfig().parse()` returns a `ConfigParser` for `<repo>/gear.ini` | used almost everywhere |
| `seuratuploader.py` | Converts a Seurat RDS file to AnnData through rpy2 (`seurat_to_anndata()`). Also runs as a CLI (`main()`) | `anndata_processor.py` (imported only when processing an RDS file) |
| `spatial_processor.py` | `process_spatial_synchronously()` picks the handler from `SPATIALTYPE2CLASS`, cleans up obs and writes Zarr | `process_uploaded_expression_dataset.cgi`, `listeners/spatial_upload_consumer.py` |
| `spatialhandler.py` | `SpatialHandler` abstract base class, with `CosMxHandler`, `CurioHandler`, `GeoMxHandler`, `VisiumHandler`, `VisiumHDHandler` and `XeniumHandler`. `SPATIALTYPE2CLASS` maps format names to handlers | `spatial_processor.py`, `store_expression_dataset.cgi`, `resources/spatialpanel.py`, `bin/*zarr*` scripts |
| `trackhub.py` | Track hub parsing and validation, and `TrackHubProcessor` (downloads files, converts bigBed to BED and `.hic` to mcool, ingests into HiGlass). Shared by the API and a consumer | `resources/track_hub.py`, `resources/gosling_spec.py`, `listeners/gosling_upload_consumer.py` |
| `userhistory.py` | `UserHistory` records and reads a user's activity history | `add_to_user_history.cgi`, `get_user_history_entries.cgi`, the search and save CGIs |

**Two RabbitMQ modules: `lib/gearqueue.py` and `lib/gear/queue.py`.** Only `lib/gearqueue.py` imports `gear.queue`, and it wraps `gear.queue.RabbitMQQueue().connect(...)` inside `gearqueue.Connection`. All application code (CGI, API, listeners) imports `gearqueue`. Treat `gear.queue` as an internal detail of `gearqueue` and do not import it directly.

### `lib/gear/utils/`

These are cross-cutting helpers. Import the specific submodule, for example `from gear.utils.job_coordination import log_line`.

| Module | Purpose | Main users |
|--------|---------|------------|
| `gene_mapping.py` | Gene symbol to Ensembl ID mapping for AnnData and var DataFrames (`update_adata_with_ensembl_ids`, `update_var_with_ensembl_ids`, `map_gene_symbols_via_mygene`) | `anndata_processor.py`, `seuratuploader.py`, `bin/add_ensembl_id_to_h5ad_missing_release.py` |
| `job_coordination.py` | Consumer job coordination: `log_line` (timestamped log lines), `check_and_record_attempt` / `clear_attempt_count` (retry limits), and `try_acquire_lock_file` / `release_lock_file` | `listeners/*_consumer.py` |
| `obs.py` | Cleanup and categorisation of obs columns (`standardize_and_sanitize_obs`, `flag_ambiguous_obs_columns`, `apply_obs_dtype_choices`) | `anndata_processor.py`, `spatial_processor.py`, `apply_obs_dtype_choices.cgi` |
| `resource_limits.py` | `set_memory_limit_from_cgroup()` sets a memory limit from the cgroup limit so that running out of memory raises `MemoryError` instead of an OOM kill. `catch_memory_error()` is a decorator | `anndata_upload_consumer.py`, `spatial_upload_consumer.py`, `process_uploaded_expression_dataset.cgi`, API resources `orthologs.py`, `projectr.py` and `mg_plotly_data.py` |

### `www/api/` (Flask REST API)

- `www/api/api.wsgi` is the mod_wsgi entry point. It adds its own directory to `sys.path` and imports `app` from `api.py` as `application`.
- `www/api/api.py` creates the Flask app and a `flask_restful.Api`, adds `lib/` to `sys.path` and registers every resource with `api.add_resource(...)`. Route groups:
  - `/plot/<dataset_id>/...` (plotly, mg_plotly, svg, tsne, mg_tsne, gosling, spatialpanel, spatial_scanpy)
  - `/h5ad/<dataset_id>/...` (genes, analyses, aggregations, orthologs, available display types)
  - `/projectr/...`
  - `/import/trackhub/<share_uid>/copy` and `/import/dataset/<share_uid>/status`
  - `/higlass/genes/<gene_symbol>`
  - `/analysis/plotTopGenesPCA`
  - `/displays/<display_id>`
- When run directly (`__main__`), it serves the same resources under an `/api` prefix on Flask's development server.
- `www/api/resources/` has one module per resource class (for example `plotly_data.py`, `mg_plotly_data.py`, `projectr.py`, `spatialpanel.py`, `track_hub.py`). `resources/common.py` holds shared helpers and is not a registered route.

Full endpoint documentation is in the [API reference](./api_reference.md).

### `www/cgi/` conventions

`www/cgi/` contains **102 `.cgi` scripts**, including a `test.cgi` stub. It also holds a few non-CGI files: `process_contact.py`, three `.R` helpers (`plotCode1gene.R`, `plot_r_analysis.R`, `run_projectR.R`), `test.php` and `test.pl`. Most are Python 3 scripts that start with `#!/opt/bin/python3`. The catalog of endpoints is in the [API reference, CGI endpoints section](./api_reference.md#cgi-endpoints).

Three representative scripts:
- `www/cgi/get_organism_list.cgi` is read-only and takes no parameters. It uses `geardb.OrganismCollection().get_all()` and serialises with `json.dumps(result, default=lambda o: o.__dict__)`.
- `www/cgi/rename_layout.cgi` reads form parameters and checks the session. It is the standard pattern shown below.
- `www/cgi/save_new_genecart_json.cgi` reads a JSON body with `json.load(sys.stdin)` instead of `cgi.FieldStorage()`. On error it prints `Status: 500 Internal Server Error` before the header.

Conventions:
- **Library path.** About 85 CGIs use `sys.path.append(os.path.abspath(os.path.join('..', '..', 'lib')))`. This depends on Apache running CGIs with `www/cgi` as the working directory. About 14 newer scripts use `Path(__file__).resolve().parents[2] / 'lib'` instead.
- **Parameters.** Almost every script (97) reads form or query parameters with `cgi.FieldStorage()` and `form.getfirst(...)`, which always returns a single string (unlike `getvalue()`, which returns a list if a parameter is repeated). Very few read a JSON body from stdin.
- **Session.** No CGI reads cookies. The browser stores the session ID in the `gear_session_id` cookie, and `apiCallsMixin` in `common.v2.js` sends it as a `session_id` form field. CGIs resolve it with `geardb.get_user_from_session_id(session_id)` (45 scripts) or `geardb.get_user_id_from_session_id` (7 scripts). A result of `None` means the user is not logged in.
- **Output.** Scripts print `Content-Type: application/json` followed by blank lines, then `json.dumps(result)`. Errors usually come back as HTTP 200 with an `error` or `success: 0` field in the JSON body. Only a handful of scripts (6) set a non-200 `Status:` header.

Minimal skeleton, based on `rename_layout.cgi`:

```python
#!/opt/bin/python3
"""Short description of what this CGI does."""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    layout_share_id = form.getfirst('layout_share_id')

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        print(json.dumps({'error': "User must be logged in."}))
        return

    layout = geardb.get_layout_by_share_id(layout_share_id)
    if not layout:
        print(json.dumps({'error': "Dataset Collection not found."}))
        return

    # ... do work, e.g. layout.save() ...
    print(json.dumps({'layout_share_id': layout.share_id}))

if __name__ == '__main__':
    main()
```

`www/p` is also a Python CGI (no extension). It redirects short permalinks (`?s=`, `?c=`, `?l=`, `?p=` ...) to the full page URL.

### `listeners/`

This directory holds the RabbitMQ consumers. Each one is started by a `systemd/*-consumer@.service` template unit (one instance per worker) or runs in its own container.

| Consumer | Queue | Work done |
|----------|-------|-----------|
| `anndata_upload_consumer.py` | `anndata_upload_jobs` | `gear.anndata_processor` |
| `spatial_upload_consumer.py` | `spatial_upload_jobs` | `gear.spatial_processor` |
| `projectr_consumer.py` | `projectr` | projectR jobs published by `resources/projectr.py` |
| `gosling_upload_consumer.py` | `trackhub_copy_jobs` | `gear.trackhub.TrackHubProcessor` |

`Dockerfile.python_base`, `Dockerfile.anndata_upload`, `Dockerfile.spatial_upload`, `Dockerfile.projectr`, `Dockerfile.gosling_upload` and `requirements.txt` build the consumer images. See [RabbitMQ consumers](./services/rabbitmq_consumers.md).

### `services/`

- **`services/projectr/`** is a small Flask app (`main.py` with `GET /status` and `POST /`) that runs projectR through rpy2 (`rfuncs.py`). It has its own `Dockerfile`, and `[projectR_service]` in `gear.ini` controls whether it is called through Cloud Run. See [projectR service](./services/projectr.md).
- **`services/spatial/`** is a HoloViz Panel app for interactive spatial viewing:
  - `panel_app.py` and `panel_app_expanded.py` are the app entry points.
  - `panel_common.py` and `common.py` hold the viewer logic.
  - `download_plugin.py` is the download plugin.
  - The `Dockerfile` builds the image, and `systemd/spatial-panel.service` serves it on port 5006. See [spatial service](./services/spatial.md).
- **`services/nemoanalytics-import-builder/`** contains only `cors-config.json.dev` and `cors-config.json.prod`: bucket CORS rules allowing `GET` from a given origin. It has no code, and it is currently untracked in git. It probably relates to the NeMO Archive import (`[nemoarchive_import]` in `gear.ini`, `www/nemoarchive_import/`).

### `systemd/`

This directory has template units and targets for each consumer (`anndata-upload-consumer@.service` / `.target`, `spatial-upload-consumer@...`, `projectr-consumer@...`, `gosling-upload-consumer@...`), a grouping `gear-consumers.target` and `gear-consumers.slice` (shared resource limits), and `spatial-panel.service`. `ExecStart` paths use `<gear_root>` placeholders. See [setup: systemd](./setup/systemd.md).

### `docker/`

This directory holds the local and dev container setup:
- **Images:** `Dockerfile`, `Dockerfile.python` and `Dockerfile.r`, built through `docker-bake.hcl`.
- **Compose:** `docker-compose.yml` and its `.template`.
- **Apache:** `000-default.conf`, `apache2.conf`, `umgear.conf` and `wsgi.*.load`.
- **Config:** `gear.ini.docker(.template)` and `php.ini`.
- **R setup:** `install_R.sh` and `install_packages.R`.
- **Database:** sample dumps (`gear-mini.sql`, `gear-devel-*.sql`) and the `mysql/` directory.

See [setup: Docker](./setup/docker.md).

### `bin/`

`bin/` holds about 120 maintenance, annotation-loading, conversion and migration scripts (for example `load_ensembl_gbk_annotations.py`, `convert_3tab_to_h5ad.py`, `upload_spatial_dataset.py` and `create_annotation_orthology_maps.py`). They are cataloged in [maintenance scripts](../misc/scripts/README.md).

---

## Frontend

The v2 UI uses Bulma CSS and ES modules. Each page is assembled with Apache server-side includes: `<!--#include virtual="/include/..." -->` pulls in `primary_nav.html`, `header_bar.html` and any selector or tile partials. An inline `<script type="module">` then imports helpers from `js/common.v2.js`. It calls `insertVersionedJS('js/<page>.js', prefs.cache_version)` and `insertVersionedCSS(...)`, so page assets are cache-busted with `www/cache_version.json` (see [cache busting guide](./misc/cache_busting_guide.md)).

### Pages (`www/*.html`)

"Imports" lists the ES-module imports of the page's main JS file.

| Page | Main JS module | Imports | Purpose |
|------|---------------|---------|---------|
| `index.html` | `js/index.js` | `common.v2.js`, dataset-collection and gene-collection selectors | Dashboard / home and gene search entry point |
| `expression.html` | `js/expression.js` | `common.v2.js`, dataset and gene selectors, `classes/tilegrid.js` | Gene expression results in a grid of dataset tiles |
| `dataset_explorer.html` | `js/dataset_explorer.js` | `common.v2.js`, dataset-collection selector | Search datasets and manage dataset collections |
| `gene_list_manager.html` | `js/gene_list_manager.js` | `common.v2.js`, `classes/genecart.v2.js` | Search, create and edit gene lists |
| `dataset_curator.html` | `js/dataset_curator.js` | `common.v2.js`, `curator_common.js`, `helpers/plot-display-config.js`, `helpers/dataset-svg-fxns.js` | Build and save single-gene displays. Uses `include/curator_common.html` and `include/plot_config/*` |
| `multigene_curator.html` | `js/multigene_curator.js` | `common.v2.js`, `curator_common.js`, `classes/gene.js`, `classes/genecart.v2.js`, gene selector | Build and save multi-gene displays |
| `compare_datasets.html` | `js/compare_datasets.js` | `common.v2.js`, `classes/facets.js`, `gene.js`, `genecart.v2.js`, `tree.js`, gene selector | Condition comparison within a dataset |
| `projection.html` | `js/projection.js` | `common.v2.js`, dataset and pattern selectors, `classes/tilegrid.js` | projectR pattern projection onto datasets |
| `sc_workbench.html` | `js/sc_workbench.js` | `common.v2.js`, `classes/analysis.js`, `analysis-ui.js`, `dataset.js`, `gene.js`, `genecart.v2.js`, `tree.js`, `helpers/stepper-fxns.js` | Single-cell analysis workbench (QC, PCA, tSNE, clustering, markers) |
| `upload_dataset.html` | `js/upload_dataset.js` | `common.v2.js`, `classes/trackhub.js` | Dataset and track hub upload wizard. Loads `include/trackhub/*.html` at runtime (see [upload pipeline](./upload_pipeline.md)) |
| `user_profile.html` | `js/user_profile.js` | `common.v2.js` | Account settings |
| `create_account.html` | `js/create_account.js` | `common.v2.js` | Account registration |
| `forgot_password.html` | `js/forgot_password.js` | `common.v2.js` | Password recovery |
| `dominoSignal.html` | inline module only | `common.v2.js`, `classes/tree.js`, pattern-collection selector | Early prototype of a cell-cell communication page (linked from `primary_nav.html`). Not committed to git; the cell-cell communication feature is being developed separately |
| `new_page_template.html` | `js/index.js` (placeholder) | `common.v2.js` | Boilerplate for starting a new v2 page |
| `contact.html` | `js/common.js`, `js/classes/user.js`, `js/comment.js` | none (jQuery / Bootstrap 4) | **Legacy v1** contact form. Submits through `cgi/create_github_issue.cgi` (and uses `get_tag_list.cgi`). `cgi/process_contact.py` is not referenced |

Other HTML lives under `www/landing/<name>/index.html` (per-project landing pages with their own `index.js`), `www/workshop/index.html` and `www/plugins/`.

### `www/js/common.v2.js`

`common.v2.js` is the shared ES module that every v2 page imports. It imports `classes/user.v2.js`. Its responsibilities:

- **Common UI setup (`initCommonUI`)**:
  - Loads site preferences with `getDomainPreferences()` (`/site_domain_prefs.json` merged with `/cache_version.json`).
  - Rebrands the page title and logos per domain, renders domain citations and injects Google Analytics.
  - Highlights the active item in `#primary-nav` using the page's `data-nav-link`, collapses or expands the sidebar (remembered in the `gear_sidebar_collapsed` cookie), and wires the close controls for Bulma notifications and modals (including Escape).
  - Calls `loadPlugins()`.
- **Session and login**:
  - `checkForLogin()` reads the `gear_session_id` cookie, validates it through `apiCallsMixin.getSessionInfo()` and builds `CURRENT_USER`.
  - `doLogin()` posts the login form and sets the cookie. Logout removes the cookie.
  - Other exports: `getCurrentUser()` and `registerPageSpecificLoginUIUpdates()`, plus internal show/hide helpers for logged-in and logged-out elements.
- **`apiCallsMixin`**: about 60 async methods wrapping axios calls to `cgi/*.cgi` and `/api/...`. Examples are `fetchDatasets`, `fetchPlotlyData`, `fetchGeneCarts`, `saveDatasetDisplay`, `fetchProjection` and `pollProjectRStatus`. They pass `apiCallsMixin.sessionId` and `colorblindMode`. This mixin is the main map from UI actions to backend endpoints.
- **Plugins**: `loadPlugin()` and `loadPlugins()` inject `plugins/<name>/<page>.html`, `.css` and `.js` for the plugins that `SITE_PREFS.enabled_plugins` lists for the current page.
- **Utilities**:
  - Notifications: `createToast()` (toasts) and `logErrorInConsole()`, plus a global `unhandledrejection` handler that shows a toast.
  - Modals: `openModal()` and `closeModal()`.
  - URL parameters: `getUrlParameter()` and `rebindUrlParam()`.
  - Other helpers: `escapeHtml()`, `convertToFormData()`, `copyToClipboard()`, `guid()`, `trigger()`, `commonDateTime()`, `disableAndHideElement()` / `enableAndShowElement()`, `insertVersionedJS()` / `insertVersionedCSS()`, `loadDomainFunding()` (fetches `include/by_domain/<domain>/funding.html`) and `getRootUrl()`.

`www/js/curator_common.js` is a second shared module for the two curator pages. It covers the dataset tree, analysis and plot type selection, the facet widget, loading plot-config partials (`includeHtml()`), plot creation and cloning or saving displays, with `register*` hooks for page-specific behaviour.

### `www/js/classes/`

| File | Contents |
|------|----------|
| `analysis.js` | `Analysis` plus one class per sc_workbench step (`AnalysisStepPrimaryFilter`, `QCByMito`, `SelectVariableGenes`, `PCA`, `tSNE`, `Clustering`, `MarkerGenes`, `CompareGenes`, `LabeledTsne`) |
| `analysis-ui.js` | `AnalysisUI` DOM bindings and step-blocking helpers for the workbench |
| `citation.js` | `Citation`, which formats dataset citations for tiles |
| `dataset.js` | `Dataset`, a lightweight client-side dataset model |
| `facets.js` | `FacetWidget`, a faceted filter widget (compare_datasets) |
| `gene.js` | `Gene` and `WeightedGene` |
| `genecart.v2.js` | `GeneCart`, `WeightedGeneCart` and `LabeledGeneCart`, which save through `apiCallsMixin` |
| `genecart.js` | **Legacy v1** `GeneCart` / `WeightedGeneCart`. Not loaded by any page (only mentioned in CGI docstrings) |
| `tilegrid.js` | `TileGrid` and `DatasetTile`, which lay out and render the display tiles on expression.html and projection.html |
| `trackhub.js` | `Hub`, `HubContainer`, `Track` and `TrackContainer` for the track hub upload UI |
| `tree.js` | `Tree` base class with `ProjectionSourceTree`, `GeneCartTree`, `ProfileTree` and `DatasetTree` (built on wunderbaum) |
| `user.v2.js` | `User` model used by `common.v2.js` |
| `user.js` | **Legacy v1** `User`. Loaded only by `contact.html` |

### `www/js/helpers/`

| File | Contents |
|------|----------|
| `dataset-svg-fxns.js` | Colours SVG anatomical displays by expression (`colorSVG`, `drawSVGLegend`). Imports d3 v7 from jsDelivr |
| `plot-display-config.js` | Per-plot-type display config, plus `setHeatmapHeightBasedOnGenes`, `adjustStackedViolinHeight` and `attachAxisLabelTooltips` |
| `stepper-fxns.js` | Step/wizard state helpers (`openNextStepWithHrefs`, `passStepWithHref`, `failStepWithHref`, `resetStepperWithHrefs`) |

### `www/include/` partials

| Partial | Loaded by | Notes |
|---------|-----------|-------|
| `primary_nav.html`, `header_bar.html` | SSI in every v2 page | Left nav and top bar (login form, user menu, help links) |
| `dataset-collection-selector/`, `gene-collection-selector/`, `pattern-collection-selector/` | SSI for the `.html`. The page JS imports the `.js` | Each directory has `.html`, `.js` and `.css`: reusable selector widgets |
| `tile-grid/tile.html` | SSI (expression.html, projection.html) | Tile template used by `tilegrid.js` |
| `curator_common.html` | SSI (dataset_curator.html, multigene_curator.html) | Shared curator layout |
| `plot_config/pre_plot/*.html`, `plot_config/post_plot/*.html` | Runtime `includeHtml()` in the curator JS | Per-plot-type option forms |
| `trackhub/hub.html`, `trackhub/track.html` | Runtime `includeHtml()` in `upload_dataset.js` | Track hub form templates |
| `by_domain/<domain>/` (`funding.html`, `footer.html`, `site_label_bar.html`, `page_title_root.html`, `index_highlighted_dataset.html`) | `loadDomainFunding()` (v2). jQuery `.load()` in `common.js` (v1) | Per-site branding for the gear, nemo, sengear, gcid, inflammation, cancergear and node-cmtrf domains |
| `navigation_bar.html` | jQuery `.load()` in legacy `common.js` | **Legacy v1** navbar |
| `create_account.html` | No current references found | Probably a legacy v1 fragment |

### `www/plugins/`

Site-specific add-ons (`aro_2026_tab/`, `deafness_gene_annotation/`, `hrp_landing_tab/`). They are loaded per page by `loadPlugins()` according to `enabled_plugins` in `site_domain_prefs.json`. See [plugins](./services/plugins.md).

### Vendor libraries

- **Local, `www/js/vendor/`:** v2 pages use only `js.cookie.js`, `jsrender.20181003.min.js` (index, expression, create_account, forgot_password) and `snap.svg-min.js` (curator, expression, projection). The rest of the directory (jQuery 1.11, Bootstrap plugins, jquery.fileupload*, mCustomScrollbar, freewall and similar) is v1-era.
- **CSS:** `www/css/bulma/` holds Bulma, and `www/css/vendor/` holds other vendor CSS. The v2 theme is in `common.v2.css` and `gear-theme-purple.(s)css`.
- **CDN (v2 pages):**
  - axios and accessibility-widgets (unpkg), on every v2 page.
  - Plotly 3.5.0 (cdn.plot.ly).
  - wunderbaum 0.10.0, @floating-ui, html5sortable, nice-select2 and autoComplete.js (jsDelivr).
  - interactjs and intro.js (unpkg).
  - d3 v5 (d3js.org, curator pages) and d3 v7 ESM (`dataset-svg-fxns.js`).
- **CDN (legacy pages):** jQuery 3.3.1, jQuery UI, Bootstrap 4.1.3, popper.js, select2, bootstrap-select and Plotly 1.54.1.

See [webpage dependencies](./misc/webpage_dependencies.md) for more.

### Legacy v1 files still present

| File | Status |
|------|--------|
| `www/js/common.js` | v1 shared script (jQuery). Loaded by `contact.html` only. Its User Guide link opens the GitHub wiki |
| `www/js/classes/user.js`, `www/js/classes/genecart.js` | v1 classes. `user.js` is used only by `contact.html`. `genecart.js` is unused |
| `www/js/comment.js` | v1 comment/contact form logic, used by `contact.html` |
| `www/contact.html` | v1 page, but **still linked from v2** (`include/header_bar.html`, `create_account.html` and the `include/by_domain/*/footer.html` files) |
| `www/js/demo_timecourse.js`, `www/js/workshop.js` | Not referenced by any page (unused) |
| `www/include/navigation_bar.html`, `www/include/create_account.html` | v1 partials (see the partials table) |

No v2 page loads `common.js`, `user.js` or `genecart.js`.

---

Last updated: September 2026
