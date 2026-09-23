# gEAR Configuration Reference

This page covers every server-side setting in `gear.ini`, what each key does, and which code reads it. It also covers the two JSON files that the browser reads for site-wide settings, `www/site_domain_prefs.json` and `www/cache_version.json`.

See also: [Setup overview](./setup/README.md), [RabbitMQ setup](./setup/rabbitmq.md), [HiGlass setup](./setup/higlass.md), [projectR service](./services/projectr.md), [Code map](./code_map.md), [Upload pipeline](./upload_pipeline.md).

## How configuration is loaded

### `gear.ini`

`gear.ini` sits at the repository root (for example `/opt/gEAR/gear.ini`). Git ignores it (see `.gitignore`). Make it by copying `gear.ini.template`, or `docker/gear.ini.docker.template` for Docker (which `docker/Dockerfile` copies in as `gear.ini.docker`), then fill in the placeholders. Never commit a filled-in `gear.ini`.

**Canonical loader: `lib/gear/serverconfig.py`**

```python
from gear.serverconfig import ServerConfig

servercfg = ServerConfig().parse()          # configparser.ConfigParser
host = servercfg["dataset_uploader"]["queue_host"]
enabled = servercfg.getboolean("projectR_service", "queue_enabled", fallback=False)
```

- `ServerConfig.parse()` builds the path from the module's own location, `os.path.dirname(__file__) + "/../../gear.ini"`. That always points at `<repo>/gear.ini`, whatever the current working directory is.
- **No caching.** Every `parse()` call creates a new `ConfigParser` and reads the file again. Most modules call it once at import time and keep the result:
  - `lib/geardb.py` stores it as module attribute `geardb.servercfg`. Other code reuses it, for example `www/cgi/get_session_info.v2.cgi`.
  - `www/api/resources/projectr.py` and `bin/profile_*.py` store it as `this.servercfg`.
  - `listeners/*_consumer.py` store it as module-level `servercfg`.
  - `lib/gear/db.py` (`Connection.connect()`) and `lib/gear/userhistory.py` (`UserHistory.__init__`) call `ServerConfig().parse()` each time they run.
- If `gear.ini` is missing, `ConfigParser.read()` fails silently. The error only shows up later, as a `KeyError` or `NoSectionError` when a key is accessed.

**Other loaders.** These modules build their own `configparser.ConfigParser()` and read `gear.ini` directly:

| Loader | Path used | Notes |
|---|---|---|
| `www/api/resources/track_hub.py` | `Path(__file__).parents[3] / 'gear.ini'` | Module-level `_config` |
| `www/cgi/process_uploaded_expression_dataset.cgi` | `gear_root / 'gear.ini'` | Module-level `_config` |
| `www/cgi/send_email.cgi`, `www/cgi/get_users_layout_members.cgi`, `www/cgi/get_shared_info.cgi` | `'../../gear.ini'` (relative to the CWD, which is `www/cgi`) | `get_shared_info.cgi` reads the file but uses no keys |
| `tests/*.py` (Selenium/SeleniumBase tests) | `'../gear.ini'` | Must be run from `tests/` |
| Many `bin/*.py` loaders (e.g. `load_gene_ontology.py`, `rescore_*`, `add_datasets_to_group.py`) | `'gear.ini'` or `'../gear.ini'` | Depends on the CWD. Only `[database]` is read |

**Docker.** `listeners/Dockerfile.{anndata_upload,gosling_upload,spatial_upload,projectr}` copy the host's `gear.ini` into the consumer image. So when `gear.ini` changes, the listener images must be rebuilt, or the file mounted into them.

### Environment variables

A few settings come from environment variables, not `gear.ini`:

| Variable | Read by | Purpose |
|---|---|---|
| `ENVIRONMENT` | `lib/gear/trackhub.py`, `www/api/resources/{gosling_spec,spatialpanel,track_hub}.py`, `www/cgi/finalize_uploaded_expression_dataset.cgi` | `development` changes the URL/host handling for local Docker. Set in `docker/docker-compose.yml.template` |
| `DEBUG` | `www/api/api.py`, `services/projectr/main.py` | Flask debug mode |
| `PORT` | `services/projectr/main.py` | Cloud Run listen port (default 8080) |
| `GITHUB_ACCESS_TOKEN` | `www/cgi/create_github_issue.cgi` | Token used to open GitHub issues |
| `GEAR_STATS_ALLOWED_ORIGIN` | `www/cgi/stats.cgi` | Optional CORS origin. The GA4 property ID is hard-coded in that file |

The standalone Cloud Run service (`services/projectr/`) does not read `gear.ini`.

## `gear.ini` sections

The example values below come from `gear.ini.template`. Passwords and other secrets are shown as `<redacted>`. `docker/gear.ini.docker.template` has the same keys but different defaults: `database.host = db`, every `queue_host = queue`, `projectR_service.auth_enabled = 0`, `nemoarchive_import.queue_enabled = 0`, `dataset_uploader.queue_enabled = 0`.

Boolean keys are read with `ConfigParser.getboolean`, so `1/0`, `true/false`, `yes/no` and `on/off` all work.

### `[database]`

MySQL/MariaDB connection. See [MySQL setup](./setup/mysql.md).

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `user` | DB user name | `gear` | `lib/gear/db.py`, `bin/create_test_mysql_dump.py`, about 20 `bin/*.py` loaders |
| `password` | DB password | `<redacted>` | same as above |
| `host` | DB host name | `localhost` (`db` in Docker) | same as above |
| `name` | Database (schema) name | `gear_portal` | same as above |
| `port` *(not in template)* | DB port | `3306` (default) | `lib/gear/db.py`, via `config['database'].get('port', 3306)` |

`lib/gear/db.py` (`Connection`) is the connection that `geardb` and nearly every CGI/API module use.

### `[email_sender]`

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `address` | Sender address for outgoing mail (SMTP via `smtp.gmail.com:587`) | `[gear_email]` | `www/cgi/send_email.cgi` |
| `password` | SMTP password (app password) for that address | `<redacted>` | `www/cgi/send_email.cgi` |

### `[content]`

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `default_layout_share_id` | Share ID of the layout (dataset collection) shown to anonymous users and to users with no saved layout | `[default_layout_share_id]` | `www/cgi/get_session_info.v2.cgi` (via `geardb.servercfg`), `www/cgi/get_users_layout_members.cgi` |

### `[test]`

Credentials for the browser-automation tests in `tests/`. See [Testing](./testing.md).

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `host` | Base URL of the portal under test | `http://localhost/` | `tests/accounts__create_account.py`, `tests/accounts__log_in.py`, `tests/datasets__upload_bulk-rnaseq.py` |
| `user_name` | Display name for the test account | `Gear Tester` | `tests/accounts__create_account.py` |
| `user_email` | Test account e-mail (can be fake) | `testing@testgear.org` | `tests/accounts__*.py`, `tests/datasets__upload_bulk-rnaseq.py`, `tests/test_{compare_datasets,front_page,multigene_curator,sc_workbench}.py` |
| `user_institution` | Institution for the test account | `Institute of Things Not Breaking` | `tests/accounts__create_account.py` |
| `password` | Test account password | `<redacted>` | same files as `user_email` |

### `[folders]`

IDs of the top-level ("master") folder rows for profiles (layouts) and gene carts. These rows must exist in the `folder` table.

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `profile_domain_master_id` | Root folder for domain (site-curated) profiles | `101` | `lib/geardb.py` `FolderCollection.get_root_folders()` |
| `profile_user_master_id` | Root folder for a user's own profiles | `102` | same |
| `profile_group_master_id` | Root folder for group profiles | `103` | same |
| `profile_shared_master_id` | Root folder for profiles shared with the user | `104` | same |
| `profile_public_master_id` | Root folder for public profiles | `105` | same |
| `cart_domain_master_id` | Root folder for domain gene carts | `106` | same |
| `cart_user_master_id` | Root folder for a user's own gene carts | `107` | same |
| `cart_group_master_id` | Root folder for group gene carts | `108` | same |
| `cart_shared_master_id` | Root folder for gene carts shared with the user | `109` | same |
| `cart_public_master_id` | Root folder for public gene carts | `110` | same |

The key name is built at runtime as `f"{folder_type}_{scope}_master_id"`, so a plain grep for these names finds nothing. `get_root_folders()` is only called from `FolderCollection.get_tree_by_folder_ids()`, and nothing in `lib/`, `www/` or `bin/` calls that method. In practice the whole section is probably unused.

### `[history]`

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `enable_user_history` | Records the user's own actions (searches, saved layouts and carts) for the dashboard | `1` | `lib/gear/userhistory.py` (`getboolean`, no fallback) |
| `history_count` | How many recent entries of each type to show | `10` | `lib/gear/userhistory.py` (`getint`, no fallback) |

There is no fallback, so a missing `[history]` section raises an error in every CGI that imports `UserHistory` (`search_genes.cgi`, `search_datasets.cgi`, `add_layout.cgi`, and others).

### `[projectR_service]`

Settings for projectR pattern projection. See [projectR service](./services/projectr.md) and [RabbitMQ setup](./setup/rabbitmq.md).

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `hostname` | Cloud Run service URL. Used both as the POST URL and as the OIDC token audience | `[cloud_run_service_url]` | `www/api/resources/projectr.py` (`fetch_one`), `bin/profile_single_projectr_tsne_run.py` |
| `cloud_run_enabled` | `1` sends projections to Cloud Run. `0` runs projectR locally through rpy2 | `1` | `www/api/resources/projectr.py` |
| `queue_enabled` | `1` publishes projection jobs to RabbitMQ (handled by `listeners/projectr_consumer.py`). `0` runs them in the API process | `0` | `www/api/resources/projectr.py` |
| `queue_host` | RabbitMQ host for projectR jobs | `localhost` (`queue` in Docker) | `www/api/resources/projectr.py`, `listeners/projectr_consumer.py` |
| `full_output` | Default for requests that do not set `full_output`. `1` also returns the p-value matrix (NMF only) | `1` | `www/api/resources/projectr.py` |
| `auth_enabled` | `1` attaches a GCP OIDC bearer token (needed for production Cloud Run). `0` skips auth for local use | `1` (`0` in Docker) | `www/api/resources/projectr.py` (`get_auth_headers`) |

### `[nemoarchive_import]`

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `importer_id` | ID of the NeMO Archive importer | `[importer_id]` | **No code reads it** |
| `gcp_project_id` | GCP project for the NeMO import | `[gcp_project_id]` | **No code reads it** |
| `credentials_json` | Path to the GCP service-account JSON | `[path_to_credentials_json]` | **No code reads it** |
| `queue_enabled` | Enables the RabbitMQ queue for NeMO imports | `1` | **No code reads it** |
| `queue_host` | RabbitMQ host for NeMO imports | `localhost` | **No code reads it** |

No tracked code references `nemoarchive_import`. `.gitignore` still ignores `www/nemoarchive_import/*.json`. The untracked `services/nemoanalytics-import-builder/` folder holds only CORS config files. Treat the section as left over from an older NeMO Archive import feature.

### `[dataset_uploader]`

RabbitMQ settings for dataset uploads (AnnData/h5ad, spatial, and track hub/Gosling). See [Upload pipeline](./upload_pipeline.md) and [RabbitMQ consumers](./services/rabbitmq_consumers.md).

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `queue_enabled` | `1` publishes upload jobs to RabbitMQ. `0` makes the publisher raise `QueueDisabledError`, or log a message, and the caller falls back to synchronous processing | `1` | `www/cgi/process_uploaded_expression_dataset.cgi` (`queue_spatial_job`, `queue_anndata_job`), `www/api/resources/track_hub.py` (`queue_trackhub_job`) |
| `queue_host` | RabbitMQ host for upload jobs (publisher and consumers) | `localhost` (`queue` in Docker) | the same two publishers, plus `listeners/anndata_upload_consumer.py`, `listeners/spatial_upload_consumer.py`, `listeners/gosling_upload_consumer.py` |

### `[higlass]`

HiGlass server used for track hub and Gosling track ingestion. See [HiGlass setup](./setup/higlass.md) and [HiGlass file upload](./misc/higlass_file_upload.md).

| Key | Meaning | Example | Read by |
|---|---|---|---|
| `admin_user` | HiGlass (Django) admin user for tileset uploads | `[user]` | `www/api/resources/track_hub.py`, `listeners/gosling_upload_consumer.py`. Both pass it to `lib/gear/trackhub.py` as `higlass_admin_user` |
| `admin_pass` | Password for that user | `<redacted>` | same, as `higlass_admin_pass` |
| `hostname` | Base URL of the HiGlass server | `[host]` | same, as `higlass_hostname` |

All three are read with an empty-string fallback, so a missing section does not raise an error, but HiGlass uploads will fail.

### RabbitMQ connection details

There is no `[rabbitmq]` section. Each feature sets its own `queue_host`. `lib/gear/queue.py` (which `lib/gearqueue.py` wraps) connects with `pika.ConnectionParameters(host=host)` only, which means the default port 5672 and the default `guest` credentials. Changing the port or credentials needs a code change. See [RabbitMQ setup](./setup/rabbitmq.md).

## Unused / unverified keys

Checked by grepping tracked code under `lib/`, `www/`, `listeners/`, `services/`, `bin/`, `tests/`, `docker/` and `systemd/` for both `['section']['key']` and `.get*('section', 'key')` forms.

| Key | Status |
|---|---|
| `[nemoarchive_import] importer_id`, `gcp_project_id`, `credentials_json`, `queue_enabled`, `queue_host` | Not read anywhere (whole section unused) |
| `[folders]` (all 10 keys) | Read only by `FolderCollection.get_root_folders()`, which only `get_tree_by_folder_ids()` calls, and nothing calls that |

## Keys read in code but missing from the template

| Key | Read by | Notes |
|---|---|---|
| `[database] port` | `lib/gear/db.py` | Optional. Defaults to 3306. The `bin/*.py` loaders that connect with `mysql.connector` directly ignore it |

## `www/site_domain_prefs.json`

This file holds branding and feature switches for each deployment (for example umgear.org or a NeMO portal). It is tracked in git, and each deployment edits it in place. A `site_domain_prefs.json.bak` file sits next to it. There are no `site_domain_prefs.<domain>.json` variants in the repository.

| Key | Meaning | Read by |
|---|---|---|
| `domain_label` | Internal domain ID (`gear`, `nemo`, ...). Used for per-domain images under `img/by_domain/<label>/` and for domain-specific UI | `lib/geardb.py` (`geardb.domain_label`), `www/cgi/send_email.cgi`, `www/js/common.js`, `www/js/common.v2.js` |
| `domain_short_display_label` | Short brand name shown in the UI and in e-mails (for example `gEAR`) | `lib/geardb.py` (`geardb.domain_short_label`), `www/js/common*.js`, `www/index.html`, `www/include/navigation_bar.html` |
| `domain_tagline` | Tagline under the brand name | `www/js/common.v2.js` |
| `citations` | List of `{lines: [...], links: [{label, url}]}` shown in the "cite us" UI | `www/js/common.v2.js` |
| `domain_url` | Public base URL of the site | `lib/geardb.py` (`geardb.domain_url`), `www/cgi/send_email.cgi`, `www/cgi/finalize_uploaded_expression_dataset.cgi`, `www/api/resources/{gosling_spec,spatialpanel,track_hub}.py`, `bin/migrate_epiviz_to_gosling.py` |
| `enabled_plugins` | Map of plugin name to the list of pages it loads on. See [Plugins](./services/plugins.md) | `www/js/common.js`, `www/js/common.v2.js` |
| `google_analytics_4_measurement_id` | GA4 measurement ID injected with `gtag` | `www/js/common.js`, `www/js/common.v2.js` |
| `links_out` | Switches for external link-outs (`HomoloGene`, `DVD`) | `lib/geardb.py` (`geardb.links_out`) |

How it is loaded:

- **Python.** `lib/geardb.py` `_read_site_domain_config()` reads `<repo>/www/site_domain_prefs.json` when the module is imported. It returns `{}` if the file is missing, which happens in consumer containers that do not mount `www/`. The values are exposed as module attributes: `geardb.domain_url`, `domain_label`, `domain_short_label`, `links_out`.
- **Browser.** `getDomainPreferences()` in `www/js/common.v2.js` fetches `/site_domain_prefs.json` and `/cache_version.json` and merges `cache_version` into the returned object. Legacy pages use `www/js/common.js`, which loads the file with jQuery AJAX into `SITE_PREFS`.

## `www/cache_version.json`

This file has one key, `cache_version` (a timestamp such as `2026.09.03.212437`).

- **Writer:** the git pre-commit hook `.githooks/pre-commit`, registered as a local hook in `.pre-commit-config.yaml`. On each commit that stages other files, it rewrites the value with `date +"%Y.%m.%d.%H%M%S"`.
- **Reader:** `getDomainPreferences()` in `www/js/common.v2.js`. Pages pass `prefs.cache_version` to `insertVersionedCSS()` / `insertVersionedJS()`, which add `?v=<version>` to asset URLs.
- `bin/audit_cache_busting.py` also reads it, to report the current version.

See the [Cache busting guide](./misc/cache_busting_guide.md) for details.

---

Last updated: September 2026
