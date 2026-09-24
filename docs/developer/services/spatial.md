# Spatial Panel Service

The Spatial Panel service is an interactive viewer for spatial transcriptomics datasets built with [Panel](https://panel.holoviz.org) (HoloViz). The gEAR Flask API prepares per-gene data, and the Panel server renders it in the browser over a Bokeh websocket.

## Overview

- **Technology**: Panel / Bokeh / Datashader, Python 3.14
- **Code**: `services/spatial/`
- **Source data**: SpatialData Zarr stores in `www/datasets/spatial/<dataset_id>.zarr`
- **Panel input**: CSV and NumPy cache files in `www/cache/spatial_panel/<dataset_id>/`
- **Deployment**: systemd (`systemd/spatial-panel.service`) on servers, or the `panel` service in Docker Compose
- **Port**: 5006

### Files in `services/spatial/`

| File | Purpose |
| --- | --- |
| `panel_app.py` | Condensed viewer (`CondensedSpatialViewer`), used in dataset tiles |
| `panel_app_expanded.py` | Expanded/zoomed viewer (`ExpandedSpatialViewer`) |
| `panel_common.py` | Viewer classes; reads settings from `pn.state.session_args` |
| `common.py` | Cache readers (`retrieve_dataframe`, `retrieve_image_array`, `list_image_channels`, ...) |
| `download_plugin.py` | Tornado handler mounted at `/spatial_download` that exports the current view as standalone HTML |
| `Dockerfile`, `requirements.txt` | Container image for the Panel server |

If either app is opened with no URL arguments it renders "OK", which is a quick health check.

## Request Flow

```
Browser (dataset tile / expanded view)
    │  POST /api/plot/<dataset_id>/spatialpanel
    ▼
Flask API (www/api/resources/spatialpanel.py)
    │  - reads www/datasets/spatial/<dataset_id>.zarr
    │    (platform from tables/table/uns/platform, AnnData from tables/table)
    │  - writes www/cache/spatial_panel/<dataset_id>/<gene>.csv
    │    and spatial_img_<channel>.npy + spatial_props.json (if the platform has images)
    │  - returns a Bokeh server_document() <script> tag
    ▼
Browser loads <domain>/panel/ws/panel_app (or panel_app_expanded)
    │  Apache proxies /panel and /panel/ws to port 5006
    ▼
Panel server reads the cached CSV/NumPy files and renders the plot
```

The `script` returned by the API points at `panel_app_expanded` when the request sets `is_zoomed`, otherwise `panel_app`. URL arguments passed to Panel include `dataset_id`, `gene_symbol`, `projection_id`, `expression_min_clip`, `nosave`, `filename`, and optionally `x_range_start`/`x_range_end`/`y_range_start`/`y_range_end`.

A second endpoint, `POST /api/plot/<dataset_id>/spatial_scanpy` (`www/api/resources/spatial_scanpy_data.py`), produces static Scanpy images and does not use the Panel server.

## Data Storage

Spatial datasets are SpatialData Zarr stores written by `lib/gear/spatialhandler.py`. The expression table is read from `tables/table` inside the store:

```
www/datasets/spatial/<dataset_id>.zarr/
├── tables/
│   └── table/          # AnnData (obs, var, X, uns/platform, ...)
├── images/             # present for image-capable platforms
├── shapes/             # platform dependent
└── points/             # platform dependent
```

The API reads only `tables/table` for most requests (`ad.read_zarr(<store>/tables/table)`) and loads the full store with `spatialdata.read_zarr()` only the first time images need to be extracted into the cache.

Cache layout used by the Panel app:

```
www/cache/spatial_panel/<dataset_id>/
├── <gene_symbol>.csv                     # or <projection_id>_<pattern>.csv
├── spatial_img_<channel>.npy             # one per image channel
├── spatial_props.json                    # original image height/width
└── image_extraction_error.log            # only if extraction failed
```

Delete a dataset's cache directory to force regeneration.

## Setup

### Systemd (servers)

`systemd/spatial-panel.service` runs:

```
/opt/Python-3.14.4/bin/panel serve <gear_root>/services/spatial/panel_app.py \
    <gear_root>/services/spatial/panel_app_expanded.py \
    --plugins=download_plugin --address=0.0.0.0 --port=5006 --num-procs=4 \
    --allow-websocket-origin="<domain url>" --global-loading-spinner
```

with `PYTHONPATH=<gear_root>/services/spatial` and `BOKEH_RESOURCES=cdn`. Replace `<gear_root>` and `<domain url>` (and adjust the Python path) before installing:

```bash
sudo cp systemd/spatial-panel.service /etc/systemd/system/
# edit placeholders in /etc/systemd/system/spatial-panel.service
sudo systemctl daemon-reload
sudo systemctl enable --now spatial-panel.service
```

The Python used must have the packages in `services/spatial/requirements.txt`. Apache must proxy `/panel` and `/panel/ws` to port 5006 (see `docker/umgear.conf` for the directives).

### Docker Compose

The `panel` service in `docker/docker-compose.yml.template`:

```yaml
panel:
  environment:
    - BOKEH_RESOURCES=cdn
  image: adkinsrs/spatial_panel_app:latest
  networks:
    - gear
  ports:
    - "5006:5006"
  pull_policy: always
  restart: always
  volumes:
    - <spatial_panel_cache_path>:/gEAR/www/cache/spatial_panel
    - <spatial_code_path>:/gEAR/services/spatial
```

The image is built by the `panel` target in `docker/docker-bake.hcl` (context `services/spatial`). The container command serves `panel_app.py` and `panel_app_expanded.py` with `--plugins download_plugin` and allows websocket origin `localhost:8080`.

To build locally instead of pulling:

```bash
cd services/spatial
docker build -t adkinsrs/spatial_panel_app:latest .
```

### Configuration

There is no Panel-specific section in `gear.ini`. Relevant settings:

- The Panel URL is derived from the domain URL (`geardb._read_domain_url()`) plus `/panel/ws/<app>`. With `ENVIRONMENT=development`, the API uses `http://localhost:8080`.
- `--allow-websocket-origin` in the systemd unit or Dockerfile must match the public domain.

## Loading Spatial Datasets

Datasets are normally uploaded through the web uploader, which queues the conversion to the `spatial_upload_consumer` (see [rabbitmq_consumers.md](./rabbitmq_consumers.md)). To load one from the command line, add metadata first (`bin/add_excel_metadata_to_db.py`), then:

```bash
./bin/upload_spatial_dataset.py -i /path/to/data.tar.gz -t visium -d <dataset_id>
```

| Option | Description |
| --- | --- |
| `-i`, `--input_file` | Spatial tarball, `.tar` or `.tar.gz` (required) |
| `-t`, `--type` | Platform: `cosmx`, `curio`, `geomx`, `visium`, `visium_hd`/`visiumhd`, `xenium` (required) |
| `-d`, `--dataset_id` | Dataset ID (required) |
| `-org`, `--organism_id` | Organism ID; needed for `cosmx`, `curio`, `geomx` when gene symbols must be mapped to Ensembl IDs and metadata is not in the database |
| `--h5ad` | Write an `.h5ad` (table only) to `www/datasets/` instead of a Zarr store (testing) |
| `--overwrite` | Overwrite an existing Zarr store |

Output goes to `www/datasets/spatial/<dataset_id>.zarr`. See also [Uploading a Spatial Dataset](../../analyst/uploading_spatial_dataset.md).

## Monitoring and Troubleshooting

```bash
sudo systemctl status spatial-panel.service
sudo journalctl -u spatial-panel.service -f
docker compose logs panel -f            # Docker
```

Health check: open `http://localhost:5006/panel_app` with no arguments; it should display "OK".

Common problems:

- **Blank plot / websocket errors**: `--allow-websocket-origin` does not match the page's host, or Apache is not proxying `/panel/ws`.
- **Stale data after re-upload**: remove `www/cache/spatial_panel/<dataset_id>/`.
- **No background image**: check `image_extraction_error.log` in the cache directory; some platforms do not provide images.
- **Port 5006 in use**: `sudo lsof -i :5006`.
- **Memory**: `--num-procs=4` starts four worker processes; reduce it on small hosts.

## Development

Run locally without Docker:

```bash
cd services/spatial
PYTHONPATH=. panel serve panel_app.py panel_app_expanded.py \
    --plugins=download_plugin --port 5006 \
    --allow-websocket-origin=localhost:8080 --dev
```

`--dev` enables auto-reload. The apps expect cache files produced by the API, so drive them through the gEAR UI (or call `/api/plot/<dataset_id>/spatialpanel` first).

### Updating dependencies

```bash
# Docker: edit services/spatial/requirements.txt, then rebuild
cd docker
DATE=$(date +%Y-%m-%d) docker buildx bake --allow=fs.read=.. panel
docker compose pull panel && docker compose up -d panel

# Server: install into the Python used by the unit, then restart
/opt/Python-3.14.4/bin/pip install -r services/spatial/requirements.txt
sudo systemctl restart spatial-panel.service
```

## Related Documentation

- [Uploading a Spatial Dataset](../../analyst/uploading_spatial_dataset.md)
- [RabbitMQ Consumers](./rabbitmq_consumers.md) (spatial upload consumer)
- [Docker Setup](../setup/docker.md)
- [New Server Setup](../setup/new_server.md)
- Panel: <https://panel.holoviz.org>, SpatialData: <https://spatialdata.scverse.org>
