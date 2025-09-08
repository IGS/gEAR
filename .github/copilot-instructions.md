# Copilot Instructions for gEAR Portal

## Project Overview
- **gEAR** is a LAMP-stack portal for gene expression data visualization and analysis, supporting microarray, bulk/single-cell RNA-Seq, and epigenetic data. Key technologies: Python, MySQL, H5AD, D3 (SVG/Snap), Plot.ly, and (legacy) Epiviz.
- Functions:
    - Dataset and gene exploration and visualization
    - Dataset management (upload, delete, etc.)
    - Dataset collection management
    - Gene list management
    - Dataset curation (creating custom visualizations)
    - Projection analysis
    - Single-cell RNA-Seq analysis (similar to Seurat pipeline)
    - Comparison tool for gene expression data within a dataset
    - Integration with external data sources
- Epiviz is being sunsetted in favor of the Gosling epigenome browser.
- Major production instances: [UMgEAR](https://umgear.org) (hearing research), [NeMO Analytics](https://nemoanalytics.org) (brain research).

## Architecture & Key Components
- **Backend**:
	- Utility scripts in `bin/` (one-off/manual tasks)
	- Main API in `www/api/` (Flask)
	- Legacy CGI endpoints in `www/cgi/` (still used for single-purpose scripts)
	- Microservices in `services/` (e.g., ProjectR)
	- Libraries in `lib/`
	- MySQL database, HDF5/H5AD for large matrix storage
    - RabbitMQ for message brokering; consumers in `listeners/`
    - Datasets are read into `anndata.AnnData` objects.
    - Scanpy (Python package) for tSNE/UMAP plots and single-cell workbench.
    - Spatial data: backend dashboard service uses Panel (Holoviz suite); spatial data stored in Zarr file stores
        - Zarr stores are read in as `spatialdata.SpatialData` objects, and convert to `AnnData` for downstream analysis.
- **Frontend**:
	- JavaScript in `www/js/`
        - ES modules
        - any file with `.v2` listed are newer UI code replacing legacy "classic" JS modules
        - any code that uses jQuery is legacy and may not be maintained
	- D3.js (SVG rendering only) and Snap.svg for graphics
	- Plot.ly for interactive plots
	- Bulma for CSS (`www/css/bulma/`)
	- jQuery is being phased out; any code using jQuery is legacy and may not be maintained
- **Services**: Containerized microservices (e.g., `services/projectr/` for ProjectR, deployable via Docker/Cloud Run)
- **Systemd**: Service templates in `systemd/` for running background workers; use `systemctl` for management

## Developer Workflows
- **Setup**: See `docs/setup.new_server.notes.md` for server install steps. Large RAM/CPU recommended.
- **Testing**: Automated UI tests in `tests/` using Mocha, Chai, Playwright. Run with `npm test` in `tests/`. The testing suite is incomplete; API tests are planned, and front-end tests mock API responses for speed and CI compatibility.
- **Service Management**: Use systemd templates for scalable worker services. Start with `sudo systemctl start <service>`, reload with `sudo systemctl daemon-reload`, enable for reboot persistence.
- **Containerization**: Build Docker images for services (see `services/projectr/README.md`). For M1 Macs, use `--platform linux/amd64`.

## Project-Specific Patterns & Conventions
- **Branching**: Follows [nvie.com git branching model](https://nvie.com/posts/a-successful-git-branching-model/).
- **Testing Selectors**: Use data attributes in HTML for Playwright tests (see Playwright docs for `data-testid`).
- **Mocking**: Front-end tests mock API/database responses; update mocks if server-side changes.
- **Visualization**: D3.js (SVG only) and Snap.svg in `www/js/`; Plot.ly for interactive plots; Bulma for layout.
- **UI Migration**: `.v2` files in `www/js/` are ES modules and represent the newer UI code. Classic JS modules and jQuery are legacy.
- **Configuration**: Main config in `gear.ini`.
- **Documentation**: End-user docs at [GitHub Wiki](https://github.com/IGS/gEAR/wiki); developer notes in `README.md` and `docs/`.

### Code Guidelines
- Follow established coding conventions (e.g., naming, file structure).
- Write modular, reusable code with clear separation of concerns.
- Include comments and documentation for complex logic.
- Ensure accessibility and responsiveness in UI components.
- Integration and end-to-end testing are important, but we may be rushed and unfortunately skip them in favor of feature development.
- Ensure security best practices are followed (e.g., input validation, authentication).
- Javascript code should be written in a functional style where possible, using ES6+ features.
- For Javascript, prefer "const" over "function" for arrow functions.
- Python is using Ruff for linting and formatting.

## Integration Points
- **Gosling (future)**: Epigenome browser replacing Epiviz for epigenetic data visualization.
- **Epiviz (legacy)**: Embedded for epigenetic data visualization (being sunsetted).
- **ProjectR**: Microservice for matrix projection, deployed via Docker/Cloud Run; configure endpoint in `gear.ini`.
- **External Data**: H5AD files for dataset storage; MySQL for metadata.

## Examples
- To run UI tests: `cd tests && npm test`
- To build ProjectR service on M1 Mac: `docker build --platform linux/amd64 --no-cache -t projectr_service .`
- To start a systemd worker: `sudo systemctl start worker@1.service`

## Key Files & Directories
- `bin/`: One-off Python utility scripts
- `www/api/`: Main Flask API endpoints
- `www/cgi/`: Legacy CGI endpoints (still used for single-purpose scripts)
- `services/`: Microservices (ProjectR, etc.)
- `systemd/`: Service templates
- `tests/`: Automated UI tests (incomplete)
- `www/js/`: Front-end JS (classic modules, ES modules in `.v2` files)
- `www/css/bulma/`: Bulma CSS framework
- `listeners/`: RabbitMQ consumer scripts
- `gear.ini`: Main configuration
- `docs/`: Setup and developer notes

---
For questions, see [GitHub Wiki](https://github.com/IGS/gEAR/wiki) or open an issue on GitHub.
