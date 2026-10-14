# Microservices Documentation

gEAR uses several microservices for specialized functionality. This directory consolidates documentation for all services.

## Available Services

### [ProjectR Service](./projectr.md)

Matrix projection service for dimensionality reduction and analysis.

- **Location**: `services/projectr/`
- **Deployment**: Docker container on Google Cloud Run or local
- **Purpose**: Performs ProjectR matrix projections for tSNE/UMAP analysis

### [Spatial Panel Service](./spatial.md)

Dashboard service for spatial transcriptomics data visualization.

- **Location**: `services/spatial/`
- **Deployment**: systemd (`spatial-panel.service`) or Docker container (Panel app)
- **Purpose**: Interactive visualization of spatial datasets

### [RabbitMQ Consumers](./rabbitmq_consumers.md)

Background workers that process async jobs from message queue.

- **Location**: `listeners/`
- **Deployment**: Systemd services (or Docker Compose services in development)
- **Purpose**: Process H5AD uploads, spatial uploads, Gosling track uploads, and ProjectR jobs

### [How to create Plugins](./plugins.md)

Information on how to extend the gEAR platform with a domain-specific plugins

- **Location**: `www/plugins/`

## Service Architecture

```
gEAR API (Flask)
    |
    ├──> ProjectR Service (Cloud Run or local)
    |    └── Matrix projection computations
    |
    ├──> Spatial Panel Service
    |    └── Spatial data visualization
    |
    └──> RabbitMQ
         └──> Consumer Workers (systemd)
              ├── AnnData upload consumer   (anndata_upload_jobs)
              ├── Spatial upload consumer   (spatial_upload_jobs)
              ├── Gosling upload consumer   (trackhub_copy_jobs)
              └── ProjectR consumer         (projectr)
```

## Configuration

Services are configured via `gear.ini` (see `gear.ini.template`). ProjectR:

```ini
[projectR_service]
hostname = [cloud_run_service_url]
;; 0 - disable, 1 - enable
cloud_run_enabled = 1
;; 0 - disable RabbitMQ, 1 - enable. Disabling could lead to potential server crashes if many jobs are run simultaneously
queue_enabled = 0
queue_host = localhost
```

The upload consumers read `queue_enabled` / `queue_host` from `[dataset_uploader]`. The Spatial Panel service has no `gear.ini` section; see [spatial.md](./spatial.md#configuration).

## Deployment Options

### Option 1: Cloud Run (Recommended for Production)

- Deploy ProjectR as Cloud Run service
- No local R installation needed
- Auto-scaling
- RabbitMQ is recommended (see Option 2 for configuration option)
- Configuration: `cloud_run_enabled = 1`

### Option 2: Local with RabbitMQ (Development)

- Run ProjectR consumers via systemd
- Requires local R installation
- Good for development/testing
- Configuration: `queue_enabled = 1`

### Option 3: In-Process (Quick Testing)

- Run ProjectR directly in Apache process
- Simple but blocks request thread
- Only for testing
- Configuration: `cloud run enabled = 0; queue_enabled = 0`

## Quick Start

### Setting Up ProjectR Cloud Service

1. Build Docker image (see [projectr.md](./projectr.md))
2. Push to Google Artifact Registry
3. Deploy to Cloud Run
4. Update `gear.ini` with service URL

### Setting Up Local RabbitMQ Consumers

1. Install RabbitMQ (see [rabbitmq_consumers.md](./rabbitmq_consumers.md))
2. Install R and packages for ProjectR (see [../setup/r_rpy2.md](../setup/r_rpy2.md))
3. Install the unit files, including `gear-consumers.slice`, then start all consumer groups:

   ```bash
   sudo systemctl start gear-consumers.target
   ```

### Setting Up Spatial Panel

1. Install `systemd/spatial-panel.service` or use the `panel` compose service (see [spatial.md](./spatial.md))
2. Start service:

   ```bash
   sudo systemctl start spatial-panel.service
   ```

## Monitoring Services

### Cloud Run Services

- View in Google Cloud Console
- Check logs via Cloud Logging
- Monitor request counts and latency

### Monitoring Systemd Services

```bash
# Check service status
sudo systemctl status projectr-consumer@1.service

# View logs
sudo journalctl -u projectr-consumer@1.service -f

# Restart service
sudo systemctl restart projectr-consumer@1.service
```

### RabbitMQ

```bash
# Management UI (if enabled)
http://localhost:15672

# Command line
sudo rabbitmqctl list_queues
```

## Scaling Services

### Cloud Run

- Auto-scales based on request load
- Configure max instances in Cloud Run console
- Adjust CPU/memory per instance

### Systemd Services

```bash
# Start additional workers
sudo systemctl start projectr-consumer@2.service
sudo systemctl start projectr-consumer@3.service

# Or use target to manage all at once
sudo systemctl start projectr-consumer.target
```

## Troubleshooting

### ProjectR Service Issues

- **Service not responding**: Check Cloud Run logs
- **Timeout errors**: Increase timeout in Cloud Run config (currently 1200s)
- **Memory errors**: Increase memory allocation (currently 16GB)

### RabbitMQ Consumer Issues

- **Jobs not processing**: Check systemd service status
- **Consumer crashes**: Check logs via journalctl
- **Connection errors**: Verify RabbitMQ is running and `queue_host` in gear.ini

### Spatial Panel Issues

- **Panel not loading**: Check systemd service logs
- **Port conflicts**: Verify port 5006 is available
- **Zarr data errors**: Verify `www/datasets/spatial/<dataset_id>.zarr` exists; clear `www/cache/spatial_panel/<dataset_id>/`

## Development

For local development, a Docker Compose stack is used and recommended.  See [the Docker setup documentation](../setup/docker.md) for more information.

## Service-Specific Documentation

- **ProjectR**: See [projectr.md](./projectr.md) for detailed setup and deployment
- **Spatial Panel**: See [spatial.md](./spatial.md) for configuration
- **RabbitMQ Consumers**: See [rabbitmq_consumers.md](./rabbitmq_consumers.md) for worker management

## Additional Resources

- Systemd service templates: `systemd/` (see [../setup/systemd.md](../setup/systemd.md))
- Docker configurations: `docker/`
- Setup guides: [../setup/](../setup/README.md)
