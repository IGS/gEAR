# Server Setup and Configuration

This directory contains guides for setting up and configuring a gEAR server instance.

## Setup Guides

### Core Installation

- **[New Server Setup](./new_server.md)** - Complete guide for setting up a fresh gEAR instance
  - System requirements and cloud instance setup
  - Base system configuration
  - Quick reference for common setup tasks

### Component Setup

- **[MySQL Database](./mysql.md)** - MySQL installation and configuration
- **[Apache Web Server](./apache.md)** - Apache setup with mod_wsgi
- **[Python Environment](./python.md)** - Python packages and dependencies
- **[R and rpy2](./r_rpy2.md)** - R installation for ProjectR (if running locally)
- **[RabbitMQ](./rabbitmq.md)** - Message broker for async job processing

### Optional Components

- **[Epiviz](./epiviz.md)** - Epigenome browser setup (legacy, being sunsetted)
- **[Docker](./docker.md)** - Docker-based development environment
- **[Systemd Services](./systemd.md)** - Background worker services

## Quick Start

For a new server installation:

1. Start with [New Server Setup](./new_server.md)
2. Follow the component guides in order:
   - MySQL → Python → Apache → RabbitMQ (optional) → R (optional)
3. Configure `gear.ini` from template
4. Set up systemd services for background workers
5. Load annotation data and test datasets

## System Requirements

### Minimum Requirements
- **OS**: Ubuntu 22.04 LTS (recommended)
- **CPU**: 2 vCPUs (4+ recommended for production)
- **RAM**: 16GB (48GB+ recommended for production)
- **Storage**: 100GB+ SSD

### Recommended for Production
- **CPU**: 16+ cores
- **RAM**: 100GB+
- **Storage**: 300GB+ SSD (more depending on dataset sizes)

## Configuration Files

- `gear.ini.template` - Main configuration template
- `docker/gear.ini.docker.template` - Docker-specific configuration
- `create_schema.sql` - Database schema

## Data Directories

After setup, ensure these directories exist and have proper permissions:

```
www/
├── datasets/           # H5AD dataset files
├── datasets/spatial/   # Zarr spatial data stores
├── analyses/           # Analysis results
├── carts/              # Gene lists
├── projections/        # ProjectR results
├── uploads/files/      # Upload staging area
└── img/dataset_previews/  # Dataset preview images
```

### Permissions

The Apache user (usually `www-data`) needs write access:

```bash
cd /var/www  # or your gEAR www directory
chmod 777 datasets datasets/spatial analyses/* carts/ projections/ uploads/files/ img/dataset_previews/
```

## Common Issues

### Apache won't start
- Check Apache error logs: `/var/log/apache2/error.log`
- Verify mod_wsgi is installed: `apache2ctl -M | grep wsgi`
- Check gear.ini permissions and syntax

### Database connection errors
- Verify MySQL is running: `systemctl status mysql`
- Check database credentials in `gear.ini`
- Ensure database user has proper permissions

### Python import errors
- Activate virtual environment if using one
- Install missing packages: `pip install -r requirements.txt`
- Check Python path in wsgi configuration

### RabbitMQ connection issues
- Verify RabbitMQ is running: `systemctl status rabbitmq-server`
- Check credentials and virtual host in `gear.ini`
- Review RabbitMQ logs: `/var/log/rabbitmq/`

## Migration from Existing Instance

### Data Transfer

Transfer these directories from old to new instance:

```bash
$HOME/git/gEAR/www/datasets
$HOME/git/gEAR/www/analyses
```

### Database Migration

1. Export database from old instance:
   ```bash
   mysqldump -u root -p gear_db > gear_backup.sql
   ```

2. Import on new instance:
   ```bash
   mysql -u root -p gear_db < gear_backup.sql
   ```

## Docker Alternative

For development, consider using Docker instead of a full server setup:
- See [Docker Setup Guide](./docker.md)
- Faster setup but requires more disk space
- Good for local development and testing

## Maintenance

### Regular Tasks
- Update system packages: `sudo apt update && sudo apt upgrade`
- Monitor disk space (datasets can grow large)
- Review Apache/MySQL logs for errors
- Backup database regularly

### Updating gEAR
```bash
cd /path/to/gEAR
git pull origin devel
# Restart Apache if needed
sudo systemctl restart apache2
```

## Getting Help

- Check component-specific setup guides for detailed instructions
- Review logs for error messages
- Create an issue on [GitHub](https://github.com/IGS/gEAR/issues)
- Contact @adkinsrs for infrastructure questions
