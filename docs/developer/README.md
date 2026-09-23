# gEAR Developer Documentation

Welcome to the gEAR developer documentation. This guide is intended for developers working on the gEAR codebase.

## Quick Links

- [Setup Guides](./setup/README.md) - Server installation and configuration
- [Services](./services/README.md) - Microservices documentation (ProjectR, Spatial Panel, RabbitMQ consumers, plugins)
- [Architecture Overview](#architecture-overview)
- [Development Workflow](#development-workflow)

### Reference

- [API Reference](./api_reference.md) - Flask API endpoints
- [Configuration](./configuration.md) - `gear.ini` sections and keys
- [Testing](./testing.md) - Mocha/Playwright and pytest/SeleniumBase suites
- [Upload Pipeline](./upload_pipeline.md) - Dataset upload flow
- [Code Map](./code_map.md) - Where things live in the codebase
- [Database Schema](./database_schema.md) - MySQL tables

### Miscellaneous

- [Cache Busting Guide](./misc/cache_busting_guide.md)
- [Webpage Dependencies](./misc/webpage_dependencies.md)
- [Release Test Plan](./misc/release_test_plan.md)
- [Utility Scripts (`bin/`)](../misc/scripts/README.md)

## Architecture Overview

gEAR is a LAMP-stack portal with the following key components:

### Backend Components

- **Python/Flask API**: Main API endpoints in `www/api/`
- **Legacy CGI**: Single-purpose scripts in `www/cgi/` (still in use)
- **MySQL Database**: Metadata and user data storage
- **HDF5/H5AD**: Large expression matrix storage using AnnData format
- **RabbitMQ**: Message broker for async job processing
- **Libraries**: Core functionality in `lib/`

### Frontend Components

- **JavaScript**: ES modules in `www/js/`
  - `common.v2.js` holds shared page code (login, site preferences, plugins, cache-busting helpers) and is imported by every top-level page except `contact.html` and `manual.html`. The `.v2` suffix is historical; only `common.v2.js`, `classes/genecart.v2.js` and `classes/user.v2.js` (plus `css/common.v2.css`) carry it, and other current files do not.
  - jQuery-based code is legacy (being phased out)
- **D3.js & Snap.svg**: SVG-based data visualization
- **Plot.ly**: Interactive plots
- **Bulma**: CSS framework (`www/css/bulma/`)
- **Plugin Support**: in `www/plugins/`

### Services & Workers

- **ProjectR Service**: Containerized microservice for matrix projection
- **Spatial Panel**: Dashboard service using Panel (Holoviz suite)
- **RabbitMQ Consumers**: Background workers in `listeners/`
- **Systemd Services**: Service templates in `systemd/`

### Data Processing

- **Scanpy**: tSNE/UMAP plots and single-cell analysis
- **AnnData**: Primary data structure for datasets
- **SpatialData**: Zarr-based spatial data storage

## Development Workflow

### Branching Model

We follow the [nvie.com git branching model](https://nvie.com/posts/a-successful-git-branching-model/):

- `main`: Production-ready code
- `devel`: Development branch
- Feature branches: `feature/feature-name`
- Hotfix branches: `hotfix/issue-description`

### Setting Up a Development Environment

1. **Clone the Repository**

   ```bash
   git clone https://github.com/IGS/gEAR.git
   cd gEAR
   ```

2. **Follow Setup Guides**
   - See [setup/README.md](./setup/README.md) for detailed instructions
   - Recommended: 16+ cores, 100GB+ RAM for production-like environment

3. **Configure gear.ini**

   ```bash
   cp gear.ini.template gear.ini
   # Edit gear.ini with your configuration
   ```

4. **Install Git Hooks (Pre-commit)**

   This project uses [pre-commit](https://pre-commit.com/) to manage git hooks. `.pre-commit-config.yaml` currently defines one hook, `cache-version-bump`, which bumps the cache version in `www/cache_version.json` on commits with other staged changes. (Standard hooks such as trailing-whitespace are present but commented out; Ruff is not run by pre-commit.)

   **Setup:**
   ```bash
   # Install pre-commit framework
   pip install pre-commit

   # Install the git hooks defined in .pre-commit-config.yaml
   pre-commit install

   # (Optional) Run the hooks on all files to verify setup
   pre-commit run --all-files
   ```

   After this, the hooks will run automatically on every `git commit`. See the [cache busting guide](./misc/cache_busting_guide.md) for details on how the cache version bump works.

### Testing

- **UI Tests (Mocha + Playwright)**: `tests/test/*.test.js`, assertions via `expect` from `playwright/test`

  ```bash
  cd tests
  npm test
  ```

- **pytest + SeleniumBase**: End-to-end UI tests in `tests/test_*.py` that run against a live site (run from `tests/`: `pytest test_front_page.py`, add `--data=localhost` for Docker)
- **API Tests**: No dedicated suite yet
- Mocha front-end tests mock API responses for speed and CI compatibility

See [testing.md](./testing.md) and [tests/README.md](../../tests/README.md).

### Code Style Guidelines

#### Python

- Use snake_case for variables/functions, PascalCase for classes
- Follow PEP 8 style guide
- Use Ruff for linting/formatting when available
- Type hints encouraged but not required

#### JavaScript

- Functional style using ES6+ features
- Prefer `const` over `function` for arrow functions
- Prefer `for...of` loops over `forEach`
- Use camelCase for variables/functions, PascalCase for classes
- Use `getElementById`/`getElementsByClassName` over `querySelector` when possible

#### ES Module Code Order

1. Imports (must be at top)
2. Constants and variables
3. Functions and classes
4. Initialization logic
5. Event listeners (at bottom to avoid hoisting issues)

## Key Directories

```
gEAR/
├── bin/                    # Utility scripts (see docs/misc/scripts/README.md)
├── docs/                   # Documentation
│   ├── analyst/            # Analyst documentation
│   ├── developer/          # Developer documentation (you are here)
│   ├── misc/               # Miscellaneous documentation
│   ├── wiki/               # End-user wiki content
│   ├── posters/            # Historical presentations (DO NOT MODIFY)
│   └── ui-v2-design/       # UI v2 prototypes (DO NOT MODIFY)
├── docker/                 # Docker configuration
├── lib/                    # Python libraries
├── listeners/              # RabbitMQ consumer scripts
├── services/               # Microservices (ProjectR, Spatial)
├── systemd/                # Systemd service templates
├── tests/                  # Automated tests
└── www/                    # Web application
    ├── api/                # Flask API
    ├── cgi/                # Legacy CGI scripts
    ├── css/                # Stylesheets (Bulma)
    ├── datasets/           # Dataset storage (H5AD files)
    └── js/                 # JavaScript code
```

## Common Tasks

### Working with Datasets

- Datasets are stored as H5AD files (AnnData format)
- Reading datasets: `adata = anndata.read_h5ad(path)`
- Spatial data uses Zarr stores, converted to AnnData for analysis

### Adding a New Feature

1. Create a feature branch from `devel`
2. Implement changes following code guidelines
3. Write/update tests
4. Test locally
5. Create pull request to `devel`

### Running Services Locally

See [services/README.md](./services/README.md) for details on:

- Starting systemd services ([RabbitMQ consumers](./services/rabbitmq_consumers.md), [systemd](./setup/systemd.md))
- Running ProjectR locally vs Cloud Run ([projectr.md](./services/projectr.md))
- Running the spatial panel dashboard ([spatial.md](./services/spatial.md))

### Debugging

- **Backend**: Check logs in `/var/log/apache2/` or use Flask debugging
- **Frontend**: Browser developer tools
- **Services**: `sudo journalctl -u <service-name>`

## Integration Points

### External Data Sources

- **Gosling**: Epigenome browser for visualizing epigenetic data (replaced Epiviz)
- **ProjectR**: Microservice for matrix projection (Docker/Cloud Run)

### Configuration

Main configuration file: `gear.ini` (copy from `gear.ini.template`; see [configuration.md](./configuration.md))

Key sections:

- `[database]`: MySQL connection
- `[projectR_service]`: ProjectR endpoint, Cloud Run toggle, and its RabbitMQ `queue_enabled` / `queue_host`
- `[dataset_uploader]`: RabbitMQ `queue_enabled` / `queue_host` for the upload consumers
- `[nemoarchive_import]`: NeMO Archive importer settings (also has queue settings)

There is no separate `[rabbitmq]` section.

## Resources

- **GitHub Wiki**: [End-user documentation](https://github.com/IGS/gEAR/wiki)
- **Issue Tracker**: [GitHub Issues](https://github.com/IGS/gEAR/issues)
- **Production Instances**:
  - [UMgEAR](https://umgear.org) - Hearing research
  - [NeMO Analytics](https://nemoanalytics.org) - Brain research

## Getting Help

1. Check this documentation
2. Review relevant code examples in the repository
3. Ask team members (@adkinsrs maintains infrastructure docs)
4. Create an issue on GitHub

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes following the guidelines
4. Submit a pull request

All contributions should include:

- Clear commit messages
- Updated documentation
- Tests (when applicable)
- Code following style guidelines
