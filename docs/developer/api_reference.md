# gEAR REST API Reference

This page documents the Flask REST API that serves plot data, h5ad metadata, projectR jobs and import status to the gEAR front end. Every route, parameter and response key below comes from the code in `www/api/`.

For the wider codebase layout see [code_map.md](./code_map.md). Configuration keys referenced here are described in [configuration.md](./configuration.md).

## Deployment

### Application object

`www/api/api.py` builds the application:

- It sets the matplotlib backend to `Agg` so that importing scanpy in resource modules does not try to open a display.
- It lowers the process soft `RLIMIT_DATA` to 95% of its current value, so that runaway allocations raise `MemoryError` instead of taking down the host. Some handlers catch that error with `@catch_memory_error()` from `gear.utils.resource_limits` (`Orthologs.post`, `MGPlotlyData.post`, and `projectr_callback` in `resources/projectr.py`). The decorator returns `{"success": -1, ...}` with HTTP 500.
- It puts `<gear root>/lib` on `sys.path`, so resources can import `geardb`, `gear.*`, `gearqueue` and similar modules.
- It creates `app = Flask(__name__)` and `api = flask_restful.Api(app)`, imports each resource class from `www/api/resources/` and registers it with `api.add_resource(...)`.
- It reads the `DEBUG` environment variable (`1`, `true`, `True`, `TRUE`). The value is only used by `app.run(debug=debug, threaded=True)` when the file runs directly as `python api.py`.

### Serving under Apache

`www/api/api.wsgi` puts its own directory on `sys.path` and exposes `from api import app as application` for mod_wsgi. The Apache config in [setup/apache.md](./setup/apache.md) mounts it at the `/api` prefix:

```
WSGIDaemonProcess api user=www-data group=www-data processes=16 threads=2
WSGIScriptAlias /api /var/www/api/api.wsgi
```

The routes in `api.py` are registered **without** the `/api` prefix. Apache adds it, so the client calls `/api/plot/<dataset_id>` and Flask sees `/plot/<dataset_id>`. The Flask development server has no prefix. `api.py` contains a commented-out `PrefixMiddleware` example for local development without Apache.

### Authentication convention

There is no token auth layer. Handlers identify the user by the gEAR session ID and resolve it with `geardb.get_user_from_session_id(session_id)` where they need a user. The session ID is read in one of three ways:

| Source | Used by |
|---|---|
| Cookie `gear_session_id` (`request.cookies.get("gear_session_id")`) | Most resources: `PlotlyData`, `MGPlotlyData`, `TSNEData`, `MGTSNEData`, `SpatialScanpyData`, `H5ad`, `GeneSymbols`, `Aggregations`, `Analyses`, `Orthologs`, `ProjectR`, `TrackHubCopy`, `DatasetProcessingStatus` |
| JSON body key `session_id` | `AvailableDisplayTypes`, `MGAvailableDisplayTypes` |
| Form field `session_id` | `TopPCAGenes` |

Most read endpoints do not reject anonymous callers. They pass the session ID (possibly `None` or `""`) to `gear.analysis.get_analysis`, which decides which analyses the caller can see. Only `TrackHubCopy` and `DatasetProcessingStatus` return 401 when the session does not resolve to a user.

### Response conventions

- Most handlers return HTTP 200 with a `success` key: `1` means success, `2` means success with a warning in `message` (seen in `PlotlyData`, `MGPlotlyData`, `SvgData`), and `-1` or `0` means failure. Callers must check `success`, not only the HTTP status.
- Some newer handlers set real HTTP status codes (`Orthologs`, `GoslingSpec`, `TrackHubCopy`, `DatasetProcessingStatus`, `DatasetDisplay`, `ProjectRStatus`).
- `flask_restful.abort(...)` is used for 400, 403 and 404 in `resources/projectr.py`.

### Shared helpers (`www/api/resources/common.py`)

| Name | Purpose |
|---|---|
| `get_adata_from_analysis(analysis, dataset_id, session_id, backed=False)` | Loads a full `AnnData` for a non-spatial analysis (or the primary dataset when `analysis` is empty). Raises `ValueError` if the analysis is spatial. |
| `get_adata_shadow(analysis, dataset_id, session_id)` / `get_adata_shadow_from_analysis(...)` | Returns a lightweight `shadows.AnnDataShadow` over the analysis h5ad, for metadata-only reads. |
| `get_spatial_adata(analysis, dataset_id, session_id)` | Loads `AnnData` from a `SpatialAnalysis`. Raises `ValueError` if the analysis is not spatial. |
| `create_projection_adata(dataset_adata, dataset_id, projection_id)` | Builds an `AnnData` from `www/projections/by_dataset/<dataset_id>/<projection_id>.csv`, copying `obs`, `obsm` and `uns` from the source. Raises `PlotError` on failure. |
| `clip_expression_values(adata, min_clip, max_clip)` | Clips `adata.X`. Backs the `expression_min_clip` parameter. |
| `order_by_time_point(obs_df)` | Orders the `time_point` category by the `time_point_order` column, then drops that column. |
| `ANNOTATION_BEDDB_UID`, `HIGLASS_URL` | Assembly-to-HiGlass-tileset map and HiGlass server base URL, used by `higlass.py`. |
| `PROJECTIONS_BASE_DIR` | `www/projections`. |

All helpers normalize the `analysis` argument with `gear.analysis.normalize_analysis_input`, so it can be a dict (`{"id": ..., "type": ...}`), a string ID or `None`.

## Route summary

`www/api/api.py` makes 24 `api.add_resource` calls (25 URL patterns, because `PlotlyData` has two). Paths are relative to the `/api` prefix. The Methods column lists the HTTP methods the Resource class defines.

| Route | Methods | Handler class | Module |
|---|---|---|---|
| `/plot/<dataset_id>` | POST | `PlotlyData` | `resources/plotly_data.py` |
| `/plot/<dataset_id>/plotly` | POST | `PlotlyData` | `resources/plotly_data.py` |
| `/plot/<dataset_id>/mg_plotly` | POST | `MGPlotlyData` | `resources/mg_plotly_data.py` |
| `/plot/<dataset_id>/svg` | GET | `SvgData` | `resources/svg_data.py` |
| `/plot/<dataset_id>/tsne` | POST | `TSNEData` | `resources/tsne_data.py` |
| `/plot/<dataset_id>/mg_tsne` | POST | `MGTSNEData` | `resources/tsne_data.py` |
| `/plot/<dataset_id>/gosling` | GET (POST, PUT, DELETE are stubs) | `GoslingSpec` | `resources/gosling_spec.py` |
| `/plot/<dataset_id>/spatialpanel` | POST | `SpatialPanel` | `resources/spatialpanel.py` |
| `/plot/<dataset_id>/spatial_scanpy` | POST | `SpatialScanpyData` | `resources/spatial_scanpy_data.py` |
| `/higlass/genes/<gene_symbol>` | GET | `HiGlassGene` | `resources/higlass.py` |
| `/projectr/<dataset_id>` | POST | `ProjectR` | `resources/projectr.py` |
| `/projectr/<dataset_id>/output_file` | POST | `ProjectROutputFile` | `resources/projectr.py` |
| `/projectr/<projection_id>/status` | GET | `ProjectRStatus` | `resources/projectr.py` |
| `/h5ad/<dataset_id>` | GET | `H5ad` | `resources/h5ad.py` |
| `/h5ad/<share_uid>/availableAnalysisTools` | GET | `AvailableAnalysisTools` | `resources/available_analysis_tools.py` |
| `/h5ad/<dataset_id>/availableDisplayTypes` | POST | `AvailableDisplayTypes` | `resources/available_display_types.py` |
| `/h5ad/<dataset_id>/mg_availableDisplayTypes` | POST | `MGAvailableDisplayTypes` | `resources/available_display_types.py` |
| `/h5ad/<dataset_id>/aggregations` | POST | `Aggregations` | `resources/aggregations.py` |
| `/h5ad/<dataset_id>/analyses` | GET | `Analyses` | `resources/analyses.py` |
| `/h5ad/<dataset_id>/orthologs` | POST | `Orthologs` | `resources/orthologs.py` |
| `/h5ad/<dataset_id>/genes` | GET | `GeneSymbols` | `resources/gene_symbols.py` |
| `/import/trackhub/<share_uid>/copy` | POST | `TrackHubCopy` | `resources/track_hub.py` |
| `/import/dataset/<share_uid>/status` | POST | `DatasetProcessingStatus` | `resources/dataset_processing.py` |
| `/analysis/plotTopGenesPCA` | POST | `TopPCAGenes` | `resources/top_pca_genes.py` |
| `/displays/<display_id>` | GET | `DatasetDisplay` | `resources/dataset_display.py` |

## Plot resources

### PlotlyData: `POST /plot/<dataset_id>` and `/plot/<dataset_id>/plotly`

Builds a single-gene Plotly figure with `gear.plotting`. Supported `plot_type` values include `bar`, `violin`, `scatter`, `line` and `contour`. `tsne_dynamic` and `tsne/umap_dynamic` are aliases for `scatter`.

JSON body keys:

| Key | Notes |
|---|---|
| `gene_symbol` | Required. |
| `plot_type` | See above. |
| `analysis` | Analysis dict or ID. Omit to use the primary dataset. |
| `x_axis`, `y_axis` (default `raw_value`), `z_axis` | Obs columns. `z_axis` is required for `contour`. |
| `point_label`, `color_name`, `colors`, `color_palette`, `reverse_palette`, `colorblind_mode` | Color settings. |
| `facet_row`, `facet_col`, `order`, `size_by_group`, `marker_size` (default 3), `jitter` | Layout settings. |
| `x_min`, `x_max`, `y_min`, `y_max`, `x_title`, `y_title`, `vlines` | Axis settings. |
| `hide_x_labels`, `hide_y_labels`, `hide_legend` | Booleans. |
| `obs_filters` | Dict of obs column to list of allowed values. |
| `projection_id` | Plot a projectR pattern instead of a gene (see [ProjectR](#projectr-resources)). |
| `expression_min_clip` | Lower clip for expression values. |
| `return_image` | If true, return a base64 PDF instead of JSON. |
| `custom_props` | Extra kwargs passed to the plotting function. |

Response: `success`, `message`, `gene_symbol`, `plot_json` (the Plotly figure), the plot settings echoed back (`x_axis`, `y_axis`, `z_axis`, axis limits and titles, `vlines`, `point_label`, `size_by_group`, `marker_size`, `jitter`, `hide_*`, `color_name`, `facet_row`, `facet_col`), plus `plot_colors`, `plot_palette`, `reverse_palette`, `obs_filters` and `plot_order`. With `return_image` the response is `{success, message, image, image_format: "pdf"}`.

Errors (HTTP 200, `success: -1`): no JSON body, missing gene symbol, dataset not found, h5ad file not found, `PlotError` from plotting, or any other exception (as `"Encountered error: ..."`).

### MGPlotlyData: `POST /plot/<dataset_id>/mg_plotly`

Builds a multi-gene Plotly figure. The `plot_type` values handled are `dotplot`, `heatmap`, `mg_violin`, `volcano` and `quadrant`. Any other value returns `"Plot type ... is not a valid multi-gene plot option"`.

| Key group | JSON body keys |
|---|---|
| Common | `analysis`, `plot_type`, `gene_symbols`, `obs_filters`, `projection_id`, `expression_min_clip`, `colorblind_mode`, `return_image`, `custom_props`, `plot_title`, `legend_title` |
| Grouping and sorting | `primary_col` (required for `dotplot`, `mg_violin` and heatmap matrixplots), `secondary_col`, `sort_order`, `clusterbar_fields`, `subsample_limit` |
| Heatmap | `matrixplot`, `center_around_zero`, `cluster_obs`, `cluster_genes`, `flip_axes`, `distance_metric` (default `euclidean`), `hide_obs_labels`, `hide_gene_labels`, `colorscale`, `reverse_colorscale` |
| Volcano | `query_condition`, `ref_condition`, `de_test_algo` (default `t-test`), `pvalue_threshold`, `lower_logfc_threshold`, `upper_logfc_threshold`, `adj_pvals`, `annotate_nonsignificant` |
| Quadrant | `compare1_condition`, `compare2_condition`, `fold_change_cutoff` (default 2), `fdr_cutoff` (default 0.05), `include_zero_fc` |
| Violin | `stacked_violin`, `violin_add_points` |

Response: `{success, message, plot_json}`, or `{success, message, image, image_format: "pdf"}` when `return_image` is set. A `MemoryError` returns HTTP 500.

### TSNEData and MGTSNEData: `POST /plot/<dataset_id>/tsne` and `/plot/<dataset_id>/mg_tsne`

These render static matplotlib/scanpy embedding plots and return them as base64 images. Parameters come from `reqparse` parsers defined at the top of `resources/tsne_data.py`.

- Shared keys: `plot_type` (default `tsne_static`; valid values are the keys of `PLOT_TYPE_TO_BASIS`: `tsne_static`, `tsne`, `umap_static`, `pca_static`, `mg_tsne_static`, `mg_umap_static`, `mg_pca_static`), `analysis`, `x_axis` (default `tSNE_1`), `y_axis` (default `tSNE_2`), `colorize_legend_by`, `max_columns`, `expression_palette` (default `YlOrRd`), `reverse_palette`, `colors`, `order`, `flip_x`, `flip_y`, `horizontal_legend`, `marker_size`, `center_around_median`, `vmin`, `vmax`, `make_zero_gray`, `enforce_equal_aspect`, `obs_filters`, `projection_id`, `expression_min_clip`, `colorblind_mode`, `high_dpi`.
- `TSNEData` only: `gene_symbol`, `plot_by_group`, `hide_group_nonmembers`, `skip_gene_plot`, `two_way_palette`.
- `MGTSNEData` only: `gene_symbols` (list).

Response: `{success, message, image, image_format}`. The format is `webp`, or `png` if the figure exceeds the 16383 px WebP limit, or `pdf` when `high_dpi` is true. Validation failures return `success: -1`, for example a missing gene, an unavailable analysis, a missing `gene_symbol` column, or inconsistent `vmin`/`vmax`/median-centering settings.

### SvgData: `GET /plot/<dataset_id>/svg`

Returns expression values for coloring SVG anatomy images. This endpoint reads the primary dataset h5ad directly and ignores analyses.

Query args: `gene` (required), `projection_id`, `expression_min_clip`, `vmin`, `vmax`.

Response: `success` (`2` if the symbol maps to several Ensembl IDs), `message`, `scores` (keys `dataset`, `gene`, `tissue`, plus `user_defined` when `vmin` or `vmax` is given), and the expression dataframe spread as top-level keys (`data` plus obs columns). When a `cell_type` column exists, the dataframe gains `<cell_type>--mean` rows.

### GoslingSpec: `GET /plot/<dataset_id>/gosling`

Fetches a UCSC-style track hub and builds a Gosling visualization spec (BAM, BED, bigWig, bigInteract, VCF and HiC tracks).

Query args: `gene`, `assembly` and `hub_url` (all required), and `zoom` (`"true"` or `"false"`).

Response: `{success, spec, position, message, hic_found}`. `position` is `chr:start-end`, or `NA` if the gene was not located. Errors: 400 for a missing parameter or a hub fetch/parse failure, 404 if the dataset is not found. `post`, `put` and `delete` exist on the class but are empty stubs.

### SpatialPanel: `POST /plot/<dataset_id>/spatialpanel`

Prepares cached inputs for the HoloViz Panel spatial viewer and returns the script that embeds it. See [services/spatial.md](./services/spatial.md).

JSON body keys: `gene_symbol` (a gene, or a pattern name when projecting), `projection_id`, `is_zoomed`, `min_genes`, `expression_min_clip`, `disable_save` (default true), and the optional initial view range `x_range_start`, `x_range_end`, `y_range_start`, `y_range_end`. The range is used only when all four are given.

Side effects: writes `<gene_symbol>.csv` (or `<projection_id>_<gene_symbol>.csv`), `spatial_img_<channel>.npy`, `spatial_props.json` and, on image failure, `image_extraction_error.log` under `www/cache/spatial_panel/<dataset_id>/`.

Response: `{script, success, message}`. Failures return `success: 0` with the error in `message`.

### SpatialScanpyData: `POST /plot/<dataset_id>/spatial_scanpy`

Runs a quick scanpy pipeline (filter, normalize, log1p, PCA, neighbors, UMAP) on a spatial analysis and renders a comparison figure.

JSON body keys: `gene_symbols` (required), `analysis`, `projection_id`, `colorblind_mode`, `high_dpi`.

Response: `{success, message, image}`. The image is base64 WebP, or PNG when `high_dpi` is true.

## HiGlass

### HiGlassGene: `GET /higlass/genes/<gene_symbol>`

Query arg: `assembly`, which must be a key of `ANNOTATION_BEDDB_UID` (`danRer10`, `hg19`, `hg38`, `mm10`, `mm39`, `rn6`). The handler queries `HIGLASS_URL/suggest/` and returns the matching entry as-is (`chr`, `txStart`, `txEnd`, `score`, `geneName`). It prefers an exact case-insensitive match and falls back to the first hit. It returns `null` when HiGlass has no results. An unsupported assembly raises `ValueError`, which surfaces as HTTP 500.

## ProjectR resources

These endpoints project a weighted (pattern) gene cart onto a dataset. For the Cloud Run service, the RabbitMQ consumer and the `projectR_service` settings, see [services/projectr.md](./services/projectr.md). Results are stored under `www/projections/`:

- `by_dataset/<dataset_id>/<projection_id>.csv` (and `<projection_id>_pval.csv` for full output)
- `by_dataset/<dataset_id>/projections.json` and `by_genecart/<genecart_id>/projections.json` (config registries)
- `job_status/job_<projection_id>.json` (polling state)

Shared `reqparse` arguments: `genecart_id` (required), `scope`, `algorithm`, `zscore` (boolean, default false), `full_output` (boolean), and `analysis` (not used).

### ProjectROutputFile: `POST /projectr/<dataset_id>/output_file`

Checks whether a projection already exists for this dataset, gene cart, `algorithm` and `zscore` combination. It creates empty `projections.json` files if they are missing. Response: `{"projection_id": <uuid or null>}`.

### ProjectR: `POST /projectr/<dataset_id>`

Also accepts `projection_id`. When omitted, the handler derives a deterministic UUID from `(dataset_id, genecart_id, algorithm, zscore)`. `full_output` falls back to `[projectR_service] full_output` in `gear.ini` and is forced to false unless `algorithm == "nmf"`.

Flow:

1. If the CSV (and the p-value CSV, when full output is requested) already exists, return `{status: "complete", result: {success, message, projection_id, num_common_genes, num_genecart_genes, num_dataset_genes}, error}`.
2. If `job_<projection_id>.json` shows `pending`, `running` or `complete`, return it unchanged. If it shows `failed`, delete it and rerun.
3. If `<csv>.lock` exists, return `status: "running"`.
4. Otherwise write a `pending` status. If `[projectR_service] queue_enabled` is true, publish `{dataset_id, session_id, projection_id, genecart_id, scope, algorithm, zscore, full_output}` to the RabbitMQ queue `projectr` and return the pending status. If not, run `projectr_callback` synchronously and return its final status.

Status shape: `{status, result, error}`, where `status` is one of `pending`, `running`, `complete` or `failed`. Errors: 400 (unsafe path components) and 403 (resolved path outside `www/projections`) via `abort`.

### ProjectRStatus: `GET /projectr/<projection_id>/status`

Returns the contents of `job_<projection_id>.json`. When the status is `complete`, the file is deleted after it is read, so the next poll gets 404. Errors: 403 for an invalid ID or path, 404 if the status file is missing.

## h5ad metadata resources

These back the dataset curator, the multigene viewer and the expression pages. Most return `success: -1` with a `message` when the dataset or its file is missing. Spatial datasets (`dtype == "spatial"`) go through `get_spatial_adata`. Other datasets use `get_adata_shadow`.

| Resource | Input | Response keys |
|---|---|---|
| `H5ad` `GET /h5ad/<dataset_id>` | query `analysis_id` | `success`, `num_obs`, `obs_columns` (includes `X_tsne_1/2`, `X_umap_1/2`, `X_pca_1/2` when present in `obsm`; excludes replicate columns and `time_point_order`), `obs_levels` (categorical columns with at most 50 levels, including `NA` if nulls exist), `obs_levels_truncated` (categorical columns with more than 50 levels), `has_replicates` |
| `GeneSymbols` `GET /h5ad/<dataset_id>/genes` | query `analysis` | `success`, `gene_symbols` |
| `Analyses` `GET /h5ad/<dataset_id>/analyses` | cookie session | `success`, `public`, `private`: analyses (from `AnalysisCollection.get_all_by_dataset_id`) with tSNE or UMAP computed. `private` is empty for anonymous users. |
| `Aggregations` `POST /h5ad/<dataset_id>/aggregations` | reqparse `analysis_id`, `filters` (dict of column to values) | `success`, `aggregations` (list of `{name, count, items: [{name, count}]}` per categorical obs column, with zero-count categories kept and `nan` renamed `Data not available`), `total_count` |
| `AvailableDisplayTypes` `POST /h5ad/<dataset_id>/availableDisplayTypes` | JSON `dataset_id`, `session_id`, `analysis_id` | Flat dict of booleans: `scatter`, `tsne_static`, `umap_static`, `pca_static`, `tsne/umap_dynamic`, `bar`, `violin`, `line`, `svg` |
| `MGAvailableDisplayTypes` `POST /h5ad/<dataset_id>/mg_availableDisplayTypes` | JSON `dataset_id`, `session_id`, `analysis_id` | Flat dict of booleans: `dotplot`, `heatmap`, `mg_violin`, `volcano` (needs 2 or more categorical columns), `quadrant` (needs 3 or more), `mg_pca_static`, `mg_tsne_static`, `mg_umap_static` |
| `AvailableAnalysisTools` `GET /h5ad/<share_uid>/availableAnalysisTools` | path share UID | `success`, `share_id`, `available_analysis_tools`: `{dataset-curator, multigene-viewer, compare-tool, sc-workbench}` mapped to booleans from the module-level `tools` table keyed on dataset type |

Notes:

- Both display-type endpoints read `dataset_id` from the JSON body and ignore the path segment. The body must repeat it.
- `AvailableDisplayTypes` returns `success: -1` if the dataset has no obs columns, unless it is `svg-expression` with an uploaded SVG at `datasets_uploaded/<dataset_id>.svg`. `MGAvailableDisplayTypes` returns `success: -1` if there are no categorical columns.

### Orthologs: `POST /h5ad/<dataset_id>/orthologs`

Maps the requested gene symbols to the symbols present in the dataset, using ortholog files for the gene and dataset organisms.

JSON body keys: `gene_symbols` (required list), `analysis` (dict, its `id` is used), `gene_organism_id`, and `exclusive_org` (`"true"` or `"false"`, default `"false"`). When false, the handler also checks other organisms.

Response: `{"success": 1, "mapping": {<input symbol>: [<dataset symbols>]}}` with HTTP 200. Errors return HTTP 400 as `{"error": ...}`: invalid JSON, no gene symbols, dataset not found, or dataset has no organism. File-load errors also return 400 but in the form `{"success": -1, "message": ...}`. A `MemoryError` returns 500.

## Import resources

These endpoints support the dataset upload workflow. See [upload_pipeline.md](./upload_pipeline.md). Both use the per-user staging area `www/uploads/files/<session_id>/<share_uid>/`, and both return 401 if the `gear_session_id` cookie does not resolve to a user.

### TrackHubCopy: `POST /import/trackhub/<share_uid>/copy`

Imports a track hub as a Gosling-format dataset. The request is **multipart form data**, not JSON:

| Form field | Notes |
|---|---|
| `hub_json` | Required. JSON string. |
| `tracks` | Required. JSON list of track stanzas, each with an `id`. |
| `assembly` | Required. |
| `dry_run` | `"true"` skips saving files. |
| `<...>[<track_id>][file]` | Optional uploaded files. Saved into the staging area, and the filename is recorded on the stanza as `uploadedFileName`. |

The handler writes `status.json` in the staging area, sets `dataset_format: "gosling"` in the staging `metadata.json`, and builds `hub_url` as `<domain>/tracks/<dataset_uid>`. The domain is `http://localhost:8080` when `ENVIRONMENT=development`. If `[dataset_uploader] queue_enabled` is true, the job is published to the RabbitMQ queue `trackhub_copy_jobs` and the handler returns HTTP 202 `{success: true, message, job_id}`. If not, it runs `gear.trackhub.process_trackhub_synchronously` (using the `[higlass]` settings) and returns 200 or 500. See [services/README.md](./services/README.md) for the consumers.

### DatasetProcessingStatus: `POST /import/dataset/<share_uid>/status`

The single status-polling endpoint for expression, spatial and track hub uploads. The optional JSON body key is `dataset_format` (for example `expression`, `spatial` or `gosling`).

The handler returns the staging `status.json` with HTTP 200 and sets `progress` to 100 when the status is complete. While the status is `processing` (and the format is not `gosling`), it checks that the recorded `job_id` exists, or else that the legacy `process_id` is still alive (`ps -p`). If both checks fail, it returns HTTP 400 with `status: "error"`. Other errors are returned as `{success: false, message, status: "error", progress: 0}` with HTTP 401 (invalid session) or 404 (no status file).

## Other resources

### TopPCAGenes: `POST /analysis/plotTopGenesPCA`

Used by the single-cell workbench. The request is **form-encoded**: `dataset_id`, `session_id`, `analysis_id`, `analysis_type`, and `pcs` (comma-separated, 2 to 5 components). The handler runs `sc.pl.pca_loadings` and saves the PNG into the `figures/` subdirectory next to the analysis h5ad. Response: `{"success": 1}`. Failures return `success: 0` (validation) or `success: -1` (load or plot errors) with a `message`.

### DatasetDisplay: `GET /displays/<display_id>`

Returns the display record from `geardb.get_display_by_id` with HTTP 200, or `{"message": "Display not found"}` with 404.

## CGI endpoints

Legacy endpoints in `www/cgi/`, served directly by Apache at `/cgi/<name>`. Most read form or JSON parameters, look up the user from `session_id`, and print a JSON body. Each script's module docstring lists its inputs and output in more detail. 105 entries, grouped by area.

### Account/auth

| Script | Purpose | Key parameters |
|---|---|---|
| `add_to_user_history.cgi` | Record an entry in the user's activity history | session_id, entry_category, label (+ any extra form fields) |
| `check_existing_email.cgi` | Check whether an email address already has an account | email |
| `create_account.cgi` | Create a new user account and return a session ID | first-last, email, institution, password, email_updates, colorblind_mode, verification_code_long, verification_code_short |
| `get_session_info.cgi` | Return basic user info (email, name) for a session | session_id |
| `get_session_info.v2.cgi` | Return full user record for a session | session_id |
| `get_user_history_entries.cgi` | Return the user's most recent history entries | session_id, num_entries |
| `login.cgi` | Log in with email/password and return a session_id (0/-1 on failure) | user_email, user_pass |
| `login.v2.cgi` | Log in (v2 form field names) and return a session_id | user-email, user-password |
| `save_user_account_changes.cgi` | Change password (via help_id) or account settings | scope, help_id, session_id, new_password, email, institution, colorblind_mode, want_updates |
| `save_user_default_organism.cgi` | Save the user's default organism | session_id, default_org_id |
| `send_email.cgi` | Send forgot-password or user-verification email | email, scope, destination_page, verification_code_long |
| `update_password.cgi` | Update a user's password | help_id, password |
| `validate_help_id.cgi` | Validate a password-reset help_id and return the user name | help_id |

### Datasets

| Script | Purpose | Key parameters |
|---|---|---|
| `apply_note_changes.cgi` | Create, edit, or remove dataset notes | session_id, scope (new/edit/remove), note_id, dataset_id, title, ldesc, access_level |
| `download_source_file.cgi` | Download a dataset's tarball, H5AD, or metadata file | dataset_id or share_id, type (tarball/h5ad/metadata), analysis_id, session_id |
| `get_condition_list.cgi` | List unique conditions (obs class labels) within a dataset | dataset_id |
| `get_dataset_info.cgi` | Return metadata for a single dataset | dataset_id, include_shape |
| `get_dataset_list.cgi` | Return the datasets for a layout, search, permalink, or default domain profile | session_id, layout_share_id, permalink_share_id, search_terms, scope, order, only_types, default_domain |
| `get_h5ad_obs_columns.cgi` | List observation columns in a dataset's H5AD file | dataset_id |
| `get_shared_info.cgi` | Return gene/dataset/layout/owner info for a dataset or layout share | dataset_share_id, layout_share_id |
| `get_shared_users_list.cgi` | List users a dataset has been shared with (owner only) | session_id, dataset_id |
| `manage_dataset_shares.cgi` | Share/unshare a dataset with a user and list shares | session_id, dataset_id, to_share |
| `mark_dataset_for_removal.cgi` | Mark an owned dataset for later removal | session_id, dataset_id |
| `remove_dataset.cgi` | Remove an owned dataset (soft delete, drop shares and layout displays) | session_id, dataset_id |
| `save_datasetinfo_changes.cgi` | Edit metadata (title, visibility, PubMed/GEO, description) of an owned dataset | session_id, dataset_id, visibility, is_downloadable, title, pubmed_id, geo_id, ldesc |

### Displays

| Script | Purpose | Key parameters |
|---|---|---|
| `delete_dataset_display.cgi` | Delete a saved dataset display owned by the user | session_id, id |
| `get_dataset_display.cgi` | Fetch a single dataset display by ID | display_id |
| `get_dataset_display_image.cgi` | Return the static preview image URL for a display | dataset_id, display_id |
| `get_dataset_displays.cgi` | List the user's and the dataset owner's saved displays for a dataset | dataset_id, session_id |
| `get_default_display.cgi` | Get the default display ID for a dataset (user, else owner) | dataset_id, session_id, is_multigene |
| `get_embedded_tsne_display.cgi` | Build a tSNE plotly_config from columns embedded in the dataset | dataset_id |
| `get_h5ad_plot.cgi` | Render available plot types for a gene from a dataset's H5AD | dataset_id, gene_symbol, session_id, group_by, colors |
| `save_dataset_display.cgi` | Create/update a saved display and regenerate its static preview PNG | id, session_id, dataset_id, label, plot_type, plotly_config, is_local |
| `save_default_display.cgi` | Set the user's default single/multi-gene display for a dataset | session_id, dataset_id, display_id, is_multigene |

### Collections (layouts)

| Script | Purpose | Key parameters |
|---|---|---|
| `add_display_to_layout.cgi` | Add a dataset display as a member of a layout | session_id, layout_share_id, display_id |
| `add_layout.cgi` | Create a new, empty layout for the user | session_id, layout_name |
| `get_user_layouts.cgi` | Return layouts visible to the user, grouped by category with folders | session_id, layout_share_id, no_domain, include_members |
| `get_users_layout_members.cgi` | Return the display members of a layout | session_id, layout_share_id |
| `remove_dataset_from_layout.cgi` | Remove a dataset from a layout | session_id, dataset_id, layout_share_id |
| `remove_display_from_layout.cgi` | Remove a display from a layout | session_id, display_id, layout_share_id |
| `remove_layout.cgi` | Delete an owned layout (not the site default) | session_id, layout_share_id |
| `rename_layout.cgi` | Rename a layout | session_id, layout_share_id, layout_name |
| `save_layout_arrangement.cgi` | Save grid positions/sizes of displays in a layout | session_id, layout_share_id, layout_arrangement |
| `save_user_chosen_layout.cgi` | Save the user's selected layout (guser.layout_share_id) | session_id, layout_share_id |
| `update_layout_visibility.cgi` | Toggle a layout between public and private | layout_share_id, visibility |

### Gene lists (carts)

| Script | Purpose | Key parameters |
|---|---|---|
| `download_weighted_gene_cart.cgi` | Download a weighted gene cart file | share_id |
| `get_gene_cart_members.cgi` | Return the gene symbols in a gene cart | share_id, session_id |
| `get_unweighted_gene_cart_preview.cgi` | Preview genes (symbol, product) in an unweighted gene cart | share_id |
| `get_user_gene_carts.cgi` | Return gene carts visible to the user, grouped by category | session_id, share_id, cart_type, group_by_type, include_members |
| `get_weighted_gene_cart_preview.cgi` | Summarize a weighted gene cart (gene count, weight labels) | share_id |
| `remove_gene_cart.cgi` | Delete an owned gene cart | session_id, share_id |
| `save_genecart_changes.cgi` | Edit metadata of an owned gene cart | session_id, gc_id, visibility, organism_id, title, ldesc |
| `save_new_genecart_form.cgi` | Create a new gene cart from form data (pasted or uploaded genes) | session_id, new_cart_label, new_cart_organism_id, new_cart_ldesc, is_public, new_cart_upload_type, new_cart_pasted_genes |
| `save_new_genecart_json.cgi` | Create a new gene cart from JSON on stdin (weighted carts saved to file) | JSON body: GeneCart fields (session_id, label, genes, gctype, weight_labels, ...) |

### Search

| Script | Purpose | Key parameters |
|---|---|---|
| `search_datasets.cgi` | Search/filter datasets with paging and extended attributes | session_id, custom_list, search_terms, organism_ids, dtypes, date_added, ownership, layout_share_id, include_public_collection_membership, page, limit, sort_by |
| `search_gene_carts.cgi` | Search/filter gene carts with paging | session_id, search_terms, organism_ids, date_added, ownership, page, limit, sort_by |
| `search_genes.cgi` | Look up genes and their annotations by symbol | session_id, search_gene_symbol, exact_match, is_multi, layout_share_id |

### Analysis/workbench

| Script | Purpose | Key parameters |
|---|---|---|
| `copy_dataset_analysis.cgi` | Copy or move a stored analysis between analysis classes (unsaved/saved/public) | session_id, dataset_id, source_analysis_id, source_analysis_type, dest_analysis_id, dest_analysis_type |
| `delete_dataset_analysis.cgi` | Delete a user's stored (non-primary) analysis | session_id, dataset_id, analysis_id, analysis_type |
| `get_PCs_from_anndata.cgi` | Return principal component loadings (varm) for an analysis as JSON | session_id, dataset_id, analysis_id, analysis_type |
| `get_analysis_image.cgi` | Stream a plot image generated by an analysis step | session_id, dataset_id, analysis_id, analysis_type, analysis_name |
| `get_dataset_comparison.cgi` | Compare expression between two conditions (fold change, p-values, filters) | dataset_id, compare_key, condition_x, condition_y, obs_filters, fold_change_cutoff, std_dev_num_cutoff, log_transformation, statistical_test |
| `get_h5ad_dataset_list.cgi` | List H5AD datasets viewable by the user (public/user/shared) for the workbench | session_id, share_id |
| `get_stored_analysis.cgi` | Return the JSON for one stored analysis | session_id, dataset_id, analysis_id, analysis_type |
| `get_stored_analysis_list.cgi` | List stored analyses for a dataset (primary/public/user_saved/user_unsaved) | session_id, dataset_id |
| `h5ad_apply_primary_filter.cgi` | Apply cell/gene count filters to an H5AD (copies primary into user space first) | analysis_id, analysis_type, dataset_id, session_id, filter_cells_lt_n_genes, filter_cells_gt_n_genes, filter_genes_lt_n_cells, filter_genes_gt_n_cells |
| `h5ad_compare_genes.cgi` | Compare marker genes between a query cluster and a reference cluster (or all others) | analysis_id, analysis_type, dataset_id, session_id, query_cluster, reference_cluster, n_genes, method, corr_method, group_labels |
| `h5ad_find_marker_genes.cgi` | Rank marker genes per cluster and return the top-gene table | analysis_id, analysis_type, dataset_id, session_id, n_genes, compute_marker_genes |
| `h5ad_generate_clusters.cgi` | Run Leiden/Louvain clustering and rename/merge/drop clusters | analysis_id, analysis_type, dataset_id, session_id, resolution, compute_clusters, cluster_info, plot_tsne, plot_umap |
| `h5ad_generate_marker_gene_visualization.cgi` | Plot dotplot and stacked violin for chosen marker genes | analysis_id, analysis_type, dataset_id, session_id, marker_genes |
| `h5ad_generate_pca.cgi` | Compute and plot PCA and PCA variance ratio | analysis_id, analysis_type, dataset_id, session_id, compute_pca, genes_to_color |
| `h5ad_generate_tsne.cgi` | Compute neighbors/tSNE/UMAP and render plots | analysis_id, analysis_type, dataset_id, session_id, n_pcs, n_neighbors, random_state, genes_to_color, use_scaled, compute_neighbors, compute_tsne, compute_umap, plot_tsne, plot_umap |
| `h5ad_identify_variable_genes.cgi` | Normalize/log and flag highly variable genes | analysis_id, analysis_type, dataset_id, session_id, norm_counts_per_cell, flavor, n_top_genes, min_mean, max_mean, min_dispersion, regress_out, scale_unit_variance, save_dataset |
| `h5ad_preview_primary_filter.cgi` | Check that pre-generated pre-filter composition plots exist | analysis_id, analysis_type, dataset_id, session_id |
| `h5ad_qc_by_mito.cgi` | QC cells by mitochondrial gene content and optionally filter | analysis_id, analysis_type, dataset_id, session_id, genes_prefix, filter_mito_perc, filter_mito_count, save_dataset |
| `save_dataset_analysis.cgi` | Save the JSON state of an analysis pipeline | session_id, dataset_id, analysis_id, analysis_type, analysis_vetting, state, label |

### Uploads

| Script | Purpose | Key parameters |
|---|---|---|
| `apply_obs_dtype_choices.cgi` | Apply user-chosen categorical/continuous dtypes to flagged obs columns of a staged upload | session_id, share_uid, choices (JSON) |
| `check_dataset_processing_status.cgi` | Report processing status/progress of an uploading dataset | session_id, share_uid |
| `delete_upload_in_progress.cgi` | Delete a user's in-progress upload directory | session_id, share_uid, dataset_id |
| `finalize_uploaded_expression_dataset.cgi` | Finalize an upload: load metadata to MySQL and migrate H5AD/source files | session_id, share_uid, dataset_uid, dataset_format, dataset_visibility, user_pii_affirmed, perform_analysis_migration |
| `get_ambiguous_obs_columns.cgi` | Return obs columns flagged as possibly mis-typed for the upload review step | session_id, share_uid |
| `get_metadata_from_geo.cgi` | Fetch series metadata from GEO for a GEO ID | geo_id |
| `get_uploads_in_progress.cgi` | List the user's incomplete uploads | session_id |
| `process_uploaded_expression_dataset.cgi` | Queue an uploaded expression dataset for processing via RabbitMQ | session_id, share_uid, dataset_uid, dataset_format, spatial_format, dataset_type |
| `store_expression_dataset.cgi` | Save the uploaded expression data file for processing | session_id, share_uid, dataset_format, spatial_format, expected_size |
| `store_expression_metadata.cgi` | Save uploader metadata form to metadata.json | session_id, share_uid, dataset_uid, title, summary, dataset_type, taxon_id, organism, contact_*, geo_id, pubmed_id, ... |
| `upload_expression_metadata.cgi` | Parse an uploaded metadata file and return it as JSON | session_id, metadata-dataset-id, metadata-file-input |

### Projection

| Script | Purpose | Key parameters |
|---|---|---|
| `download_projection.cgi` | Download projection coefficient/p-value CSVs as a zip | dataset_id or share_id, projection_id |
| `get_pattern_element_list.cgi` | List pattern labels (and top weighted genes) for a pattern gene cart | source_id, scope |
| `get_pattern_weighted_genes.cgi` | Return all genes and weights for one pattern of a weighted gene cart | source_id, pattern_id |

### Misc

| Script | Purpose | Key parameters |
|---|---|---|
| `create_github_issue.cgi` | Create a GitHub issue from a site comment/feedback form | comment_title, comment, comment_tag, submitter_firstname, submitter_lastname, submitter_email, private_check, screenshot |
| `get_citation_from_pubmed_id.cgi` | Fetch a formatted citation for a PubMed ID (cached, throttled) | pubmed_id |
| `get_event_registration_list.cgi` | Report event attendance and the user's registration/waitlist status | session_id, min_event_id, max_event_id |
| `get_organism_list.cgi` | List organisms in the database | (none) |
| `get_tag_list.cgi` | Return tags for comment-form autocomplete | (none) |
| `load_comment.cgi` | Store a user question/comment and tags in the database | submitter_firstname, submitter_lastname, submitter_email, comment_title, comment, comment_tag, super_impressive_security_check |
| `process_contact.py` | Save a contact-form submission to a file | submitterName, InputEmail, keep_updated, super_impressive_security_check, comments |
| `set_user_event_registration.cgi` | Register or unregister a user for an event | session_id, event_id, registration_status |
| `stats.cgi` | Report site usage statistics from Google Analytics 4 | query string: days, top_n, realtime |
| `test.cgi` | Test stub (prints a hello message) | none |
| `test.php` | Test stub (non-Python) | none |
| `test.pl` | Test stub (non-Python) | none |
| `update_share_id.cgi` | Change the share ID (permalink) of an owned dataset, layout, or gene cart | session_id, share_id, new_share_id, scope |
| `validate_share_id.cgi` | Validate a layout or dataset share ID (permalink) | session_id, share_id, scope |

---

Last updated: September 2026
