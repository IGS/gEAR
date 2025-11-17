# Spatial Panel Service

The Spatial Panel service provides an interactive dashboard for visualizing spatial transcriptomics datasets using the Panel framework (Holoviz suite).

## Overview

- **Technology**: Panel (Holoviz), Python
- **Data Format**: Zarr stores (SpatialData objects converted to AnnData)
- **Deployment**: Docker container via systemd
- **Port**: 5006 (default)

## Architecture

```
User Request
    ↓
Flask API
    ↓
Spatial Panel Service (Panel app on port 5006)
    ↓
Zarr Data Store (www/datasets/spatial/)
    ↓
Rendered Visualization
```

## Data Storage

Spatial datasets are stored as Zarr file stores:
- **Location**: `www/datasets/spatial/`
- **Format**: SpatialData Zarr stores
- **Processing**: Converted to AnnData for downstream analysis

## Setup

### Prerequisites

- Python 3.8+
- Panel (installed via requirements.txt)
- SpatialData package
- Sufficient disk space for Zarr stores

### Docker Build

```bash
cd services/spatial
docker build -t panel_app .
```

### Systemd Service

The spatial panel service is managed via systemd:

**Service File**: `systemd/spatial-panel.service`

```bash
# Copy service file
sudo cp systemd/spatial-panel.service /etc/systemd/system/

# Reload systemd
sudo systemctl daemon-reload

# Enable service (persist on reboot)
sudo systemctl enable spatial-panel.service

# Start service
sudo systemctl start spatial-panel.service
```

### Docker Compose

If using docker-compose (see `docker/docker-compose.yml.template`):

```yaml
panel:
  image: panel_app
  ports:
    - "5006:5006"
  volumes:
    - ./www/datasets/spatial:/data
```

## Configuration

Configure in `gear.ini`:

```ini
[spatial]
panel_port = 5006
panel_host = localhost
data_path = /var/www/datasets/spatial
```

## Usage

### Accessing the Dashboard

Once the service is running:
```
http://localhost:5006
```

Or through gEAR's spatial visualization interface which proxies requests.

### Loading Spatial Datasets

Use the `upload_spatial_dataset.py` script:

```bash
./bin/upload_spatial_dataset.py -i /path/to/spatial/data -o dataset_id
```

See also: `docs/uploading_spatial_dataset.md`

## Monitoring

### Service Status
```bash
sudo systemctl status spatial-panel.service
```

### Logs
```bash
# Follow live logs
sudo journalctl -u spatial-panel.service -f

# View recent logs
sudo journalctl -u spatial-panel.service -n 100
```

### Docker Logs (if using docker-compose)
```bash
docker compose logs panel -f
```

## Troubleshooting

### Service Won't Start

**Check service status:**
```bash
sudo systemctl status spatial-panel.service
```

**Common issues:**
- Port 5006 already in use
- Missing Python dependencies
- Incorrect data path
- Permission issues accessing Zarr stores

### Port Conflicts

If port 5006 is in use:
```bash
# Find process using port
sudo lsof -i :5006

# Kill process if needed
sudo kill -9 <PID>
```

Or change port in configuration.

### Data Access Issues

Verify permissions:
```bash
ls -la /var/www/datasets/spatial/
# Should be readable by service user
```

### Memory Issues

Spatial datasets can be large. Increase available memory:
- Docker: Adjust in docker-compose.yml
- Systemd: Add memory limits in service file

## Development

### Running Locally

For development without Docker:

```bash
cd services/spatial
panel serve app.py --port 5006 --dev
```

The `--dev` flag enables auto-reload on code changes.

### Testing

Test with a sample spatial dataset:
```bash
# Upload test dataset
./bin/upload_spatial_dataset.py -i test_data/ -o test_dataset

# Access via browser
http://localhost:5006?dataset=test_dataset
```

## Data Format

### Input Format

Spatial datasets should include:
- Expression matrix
- Spatial coordinates
- Cell/spot metadata
- Optional: Image data

### Zarr Store Structure

```
dataset_id.zarr/
├── expression/
├── coordinates/
├── metadata/
└── images/ (optional)
```

### Conversion to AnnData

The Panel service converts SpatialData to AnnData for analysis:
```python
import spatialdata as sd
import anndata as ad

sdata = sd.read_zarr("dataset.zarr")
adata = sdata.to_anndata()
```

## API Integration

The Panel service integrates with gEAR's Flask API:

**API Endpoint**: `/api/spatial/visualize`

**Request:**
```json
{
  "dataset_id": "abc123",
  "gene": "SOX2",
  "colormap": "viridis"
}
```

**Response:**
```json
{
  "panel_url": "http://localhost:5006/render?..."
}
```

## Performance

### Optimization Tips

1. **Pre-compute common views**: Generate static images for common genes
2. **Use appropriate chunk sizes**: Optimize Zarr chunk size for access patterns
3. **Limit simultaneous requests**: Panel can handle multiple users but has limits
4. **Cache rendered outputs**: Cache frequently accessed visualizations

### Scaling

For production with many users:
- Run multiple Panel instances behind load balancer
- Use Redis for shared caching
- Consider serverless deployment (Cloud Run)

## Updates

### Updating Panel Version

```bash
# Update requirements
pip install --upgrade panel

# Rebuild Docker image
docker build -t panel_app .

# Restart service
sudo systemctl restart spatial-panel.service
```

### Updating Service Code

```bash
# Pull latest code
git pull origin devel

# Rebuild image
docker build -t panel_app .

# Restart service
sudo systemctl restart spatial-panel.service
```

## Related Documentation

- `docs/uploading_spatial_dataset.md` - Uploading spatial datasets
- `systemd/README.md` - Systemd service management
- `docker/docker_notes.md` - Docker setup

## Support

For issues specific to spatial visualization:
1. Check Panel documentation: https://panel.holoviz.org
2. Review SpatialData docs: https://spatialdata.scverse.org
3. Create GitHub issue with "spatial" label
