# RabbitMQ Consumer Services

RabbitMQ consumers are background worker processes that handle asynchronous job processing for compute-intensive tasks.

## Overview

- **Technology**: RabbitMQ message broker, Python consumers
- **Location**: `listeners/` directory
- **Deployment**: Systemd services
- **Purpose**: Async processing of dataset uploads (H5AD, spatial, Gosling tracks) and ProjectR jobs
- **Queue client**: `lib/gearqueue.py` (`Connection`, `AsyncConnection`, built on pika)
- **Logs**: each consumer appends to `/var/log/gEAR_queue/<queue_name>.log`, plus stdout/stderr in the journal

## Architecture

```
Flask API
    ↓
RabbitMQ Queue
    ↓
Consumer Workers (systemd)
    ↓
Job Processing (uploads, ProjectR)
    ↓
Result Storage (Database, File System)
```

## Available Consumers

### Anndata Upload Consumer

Facilitates uploading of various file formats into Anndata (H5AD) format.

- **Listener**: `listeners/anndata_upload_consumer.py`
- **Queue**: `anndata_upload_jobs`
- **Service Template**: `systemd/anndata-upload-consumer@.service`
- **Service Group**: `systemd/anndata-upload-consumer.target`


### Gosling Upload Consumer

Facilitates uploading of track files for the epigenome uploader.

- **Listener**: `listeners/gosling_upload_consumer.py`
- **Queue**: `trackhub_copy_jobs`
- **Service Template**: `systemd/gosling-upload-consumer@.service`
- **Service Group**: `systemd/gosling-upload-consumer.target`

### Spatial Dataset Upload Consumer

Facilitates uploading of spatial transcriptomics datasets (Visium, VisiumHD, Curio, GeoMx, CosMx, Xenium), converting them to a SpatialData object and writing the result as a Zarr store. Kept separate from the Anndata Upload Consumer since spatial uploads produce a different output format (Zarr, not H5AD) and depend on the spatialdata/spatialdata_io stack.

- **Listener**: `listeners/spatial_upload_consumer.py`
- **Queue**: `spatial_upload_jobs`
- **Service Template**: `systemd/spatial-upload-consumer@.service`
- **Service Group**: `systemd/spatial-upload-consumer.target`

### ProjectR Consumer

Processes matrix projection jobs for dimensionality reduction.

- **Listener**: `listeners/projectr_consumer.py`
- **Queue**: `projectr`
- **Service Template**: `systemd/projectr-consumer@.service`
- **Service Group**: `systemd/projectr-consumer.target`

### Queue Summary

| Consumer | Queue | `gear.ini` section (`queue_host`) | Workers in `.target` |
| --- | --- | --- | --- |
| `anndata_upload_consumer.py` | `anndata_upload_jobs` | `[dataset_uploader]` | 2 |
| `gosling_upload_consumer.py` | `trackhub_copy_jobs` | `[dataset_uploader]` | 3 |
| `spatial_upload_consumer.py` | `spatial_upload_jobs` | `[dataset_uploader]` | 2 |
| `projectr_consumer.py` | `projectr` | `[projectR_service]` | 3 |

### Consumer Group Target

`systemd/gear-consumers.target` groups every consumer's `.target` together, so all of them can be started/stopped/enabled with one command instead of one per consumer.

- **Service Group**: `systemd/gear-consumers.target`
- **Wants**: `anndata-upload-consumer.target`, `gosling-upload-consumer.target`, `spatial-upload-consumer.target`, `projectr-consumer.target`

Add new consumers to this file's `Wants=` line as they're created.

## Setup

### Prerequisites

1. **RabbitMQ Server**

   ```bash
   sudo apt install rabbitmq-server
   sudo systemctl enable rabbitmq-server
   sudo systemctl start rabbitmq-server
   ```

   See also: [../setup/rabbitmq.md](../setup/rabbitmq.md)

2. **Python Dependencies**

   ```bash
   pip install -r listeners/requirements.txt  # includes pika
   ```

3. **R and Packages** (for ProjectR)

   See: [../setup/r_rpy2.md](../setup/r_rpy2.md)

### Configuration

Queue settings live in `gear.ini` (see `gear.ini.template`). Each section that uses RabbitMQ has the same two keys:

```ini
[dataset_uploader]
;; 0 - disable RabbitMQ, 1 - enable
queue_enabled = 1
queue_host = localhost

[projectR_service]
;; 0 - disable RabbitMQ, 1 - enable. Disabling could lead to potential server crashes if many jobs are run simultaneously
queue_enabled = 0
queue_host = localhost
```

`[dataset_uploader]` is read by the anndata, spatial and Gosling upload consumers; `[projectR_service]` by the ProjectR consumer. `[nemoarchive_import]` has the same keys for the NeMO Archive importer. There are no credential or virtual-host settings; connections use pika defaults on `queue_host`.

### Installing Services

```bash
# Copy service files to systemd directory
cd systemd
sudo cp projectr-consumer@.service /etc/systemd/system/
sudo cp projectr-consumer.target /etc/systemd/system/
sudo cp gosling-upload-consumer@.service /etc/systemd/system
sudo cp gosling-upload-consumer.target /etc/systemd/system
sudo cp anndata-upload-consumer@.service /etc/systemd/system
sudo cp anndata-upload-consumer.target /etc/systemd/system
sudo cp spatial-upload-consumer@.service /etc/systemd/system
sudo cp spatial-upload-consumer.target /etc/systemd/system
sudo cp gear-consumers.target /etc/systemd/system
sudo cp gear-consumers.slice /etc/systemd/system

# Replace <gear_root> in each *@.service file with the gEAR checkout path,
# and check the interpreter path (/opt/bin/python3).
# Size MemoryHigh/MemoryMax in gear-consumers.slice for the VM (see ../setup/rabbitmq.md).

# Reload systemd
sudo systemctl daemon-reload

# Enable services to start on boot (individually...)
sudo systemctl enable projectr-consumer.target gosling-upload-consumer.target anndata-upload-consumer.target spatial-upload-consumer.target

# ...or all at once via the group target
sudo systemctl enable gear-consumers.target
```

## Starting Services

### Using Service Target (Recommended)

Start all consumers in one group:

```bash
# Start selected consumer groups
sudo systemctl start projectr-consumer.target gosling-upload-consumer.target anndata-upload-consumer.target spatial-upload-consumer.target

# Check status
sudo systemctl status projectr-consumer.target anndata-upload-consumer.target
```

### Starting Every Consumer Group at Once

`gear-consumers.target` wants every individual consumer's `.target`, so one command brings up all consumer groups:

```bash
sudo systemctl start gear-consumers.target

# Check status of everything it started
sudo systemctl status gear-consumers.target
```

### Individual Workers

Start specific numbered workers:

Using projectr-consumer as an example.

```bash
# Start worker 1
sudo systemctl start projectr-consumer@1.service

# Start worker 2
sudo systemctl start projectr-consumer@2.service

# Check status
sudo systemctl status projectr-consumer@1.service
```

### Scaling Workers

Start multiple workers for parallel processing:

```bash
# Start 4 workers
sudo systemctl start projectr-consumer@1.service
sudo systemctl start projectr-consumer@2.service
sudo systemctl start projectr-consumer@3.service
sudo systemctl start projectr-consumer@4.service

# Or use target to manage predefined workers
sudo systemctl start projectr-consumer.target
```

## Monitoring

These use projectr-consumer as an example

### Service Status

```bash
# Check all consumers
systemctl status 'projectr-consumer@*'

# Check specific consumer
sudo systemctl status projectr-consumer@1.service
```

### Logs

```bash
# Follow logs for consumer 1
sudo journalctl -u projectr-consumer@1.service -f

# View recent logs
sudo journalctl -u projectr-consumer@1.service -n 100

# View logs for all consumers
sudo journalctl -u 'projectr-consumer@*' -f
```

### RabbitMQ Monitoring

```bash
# Command line
sudo rabbitmqctl list_queues
sudo rabbitmqctl list_consumers

# Web UI (if management plugin enabled)
http://localhost:15672
# Default credentials: guest/guest
```

Enable management plugin:

```bash
sudo rabbitmq-plugins enable rabbitmq_management
```

## Systemd Service Template

The `@` symbol in service filenames indicates a template. This allows spawning multiple instances.

**Example**: `projectr-consumer@.service`

```ini
[Unit]
Description="ProjectR Consumer for RabbitMQ - #%i"
Documentation=https://github.com/IGS/gEAR/blob/main/docs/developer/services/rabbitmq_consumers.md
After=rabbitmq-server.service
StartLimitIntervalSec=300
StartLimitBurst=5

[Service]
Type=simple
Environment=APACHE_STARTED_BY_SYSTEMD=true
ExecStart=/opt/bin/python3 <gear_root>/listeners/projectr_consumer.py
Slice=gear-consumers.slice
KillMode=mixed
PrivateTmp=true
Restart=always
RestartSec=2s

[Install]
WantedBy=multi-user.target
```

The other `*@.service` files follow the same pattern (the Gosling unit omits the `StartLimit*` settings). `Slice=gear-consumers.slice` puts every worker under the aggregate memory limit defined in `systemd/gear-consumers.slice`.

**Usage:**

- `%i` is replaced with instance number
- `projectr-consumer@1.service` → Worker 1
- `projectr-consumer@2.service` → Worker 2

## Target Files

Target files group related services together.

**Example**: `projectr-consumer.target`

```ini
[Unit]
Description=ProjectR Consumer Workers
Wants=projectr-consumer@1.service projectr-consumer@2.service projectr-consumer@3.service

[Install]
WantedBy=multi-user.target
```

**Benefits:**

- Start/stop all workers with one command
- Manage workers as a group
- Ensure dependencies are met

### Grouping Targets Together

A target's `Wants=` can list other targets, not just services — so one target can manage a group of consumer groups.

**Example**: `gear-consumers.target`

```ini
[Unit]
Description=All Upload/Job RabbitMQ Consumer Workers
Wants=anndata-upload-consumer.target gosling-upload-consumer.target spatial-upload-consumer.target projectr-consumer.target

[Install]
WantedBy=multi-user.target
```

`systemctl start gear-consumers.target` starts every listed target, which in turn starts each of their `@1`/`@2`/`@3` worker instances.

## Writing Custom Consumers

Use the existing consumers as templates rather than writing raw pika code. Each one:

1. Adds `lib/` to `sys.path` and imports `gearqueue` and `ServerConfig`.
2. Sets a module-level `queue_name` and a log file under `/var/log/gEAR_queue/`.
3. Defines an `_on_request(channel, method_frame, properties, body)` callback that decodes the JSON message, does the work, and acknowledges the message.
4. Wraps `gearqueue.AsyncConnection(host=..., publisher_or_consumer="consumer", queue_name=queue_name, on_message_callback=_on_request, ...)` in a `Consumer` class that reconnects with backoff.
5. Reads `queue_host` from the relevant `gear.ini` section in `main()`.

Publishers use `gearqueue.Connection(host=..., publisher_or_consumer="publisher")` and `publish(queue_name=..., message=...)` (see `www/api/resources/projectr.py` for an example).

### Creating Service File

1. Create the consumer script in `listeners/`
2. Add `systemd/<name>@.service` and `systemd/<name>.target` (copy an existing pair; keep `Slice=gear-consumers.slice`)
3. Add the new target to `Wants=` in `systemd/gear-consumers.target`
4. Optionally add `listeners/Dockerfile.<name>`, a `docker/docker-bake.hcl` target, and a compose service
5. Copy the unit files to `/etc/systemd/system/`, `daemon-reload`, and start the target

## Troubleshooting

### Consumer Not Processing Jobs

**Check RabbitMQ connection:**

```bash
# Verify RabbitMQ is running
sudo systemctl status rabbitmq-server

# Check queue has messages
sudo rabbitmqctl list_queues
```

**Check consumer logs:**

```bash
sudo journalctl -u projectr-consumer@1.service -n 50
```

**Verify gear.ini configuration:**

- Verify `queue_host` points at the RabbitMQ server
- Verify `queue_enabled = 1` in the relevant section (otherwise the web side does not publish to the queue)

### Consumer Crashes

**Check for errors in logs:**

```bash
sudo journalctl -u projectr-consumer@1.service -p err
```

**Common issues:**

- Python import errors (missing dependencies)
- R package errors (ProjectR consumer)
- Database connection issues
- File permission errors

**Restart consumer:**

```bash
sudo systemctl restart projectr-consumer@1.service
```

### Jobs Stuck in Queue

**Verify consumers are running:**

```bash
systemctl status 'projectr-consumer@*'
```

**Check message acknowledgment:**

- Ensure `basic_ack()` is called after processing
- Check for exceptions in processing code

**Purge queue (if needed):**

```bash
sudo rabbitmqctl purge_queue projectr
```

**Warning:** This deletes all messages!

## Performance Tuning

### Worker Count

Number of workers depends on:

- Available CPU cores
- Memory per job
- Job duration

**Example:**

- ProjectR jobs: 1-4 workers (memory intensive)
- Light jobs: Up to core count

### Prefetch Count

Controls how many messages a worker fetches at once:

```python
channel.basic_qos(prefetch_count=1)  # Process one at a time
channel.basic_qos(prefetch_count=5)  # Prefetch 5 messages
```

**Lower prefetch** = Better load balancing
**Higher prefetch** = Better throughput (if jobs are quick)

### Queue Durability

Make queues survive RabbitMQ restarts:

```python
channel.queue_declare(queue='my_jobs', durable=True)
```

Send persistent messages:

```python
channel.basic_publish(
    exchange='',
    routing_key='my_jobs',
    body=message,
    properties=pika.BasicProperties(delivery_mode=2)  # Persistent
)
```

## Best Practices

1. **Always acknowledge messages** after successful processing
2. **Handle exceptions** gracefully to avoid losing jobs
3. **Use durable queues** for important jobs
4. **Monitor queue length** to detect processing issues
5. **Set appropriate timeouts** for long-running jobs
6. **Log errors** for debugging
7. **Test locally** before deploying to production

## Maintenance

### Updating Consumer Code

```bash
# Pull latest code
git pull origin devel

# Restart all consumer workers (restarting a .target does not restart the
# services it Wants=, so address the instances directly)
sudo systemctl restart 'projectr-consumer@*' 'anndata-upload-consumer@*' \
    'spatial-upload-consumer@*' 'gosling-upload-consumer@*'
```

### Viewing Queue Statistics

```bash
# Via CLI
sudo rabbitmqctl list_queues name messages consumers

# Via management UI
http://localhost:15672/#/queues
```

### Clearing Old Messages

If messages are stuck or invalid:

```bash
# Purge specific queue
sudo rabbitmqctl purge_queue projectr
```

## Related Documentation

- [RabbitMQ setup](../setup/rabbitmq.md) - installation and memory sizing
- [systemd](../setup/systemd.md) - service management
- [R / rpy2 setup](../setup/r_rpy2.md) - R setup for ProjectR
- [ProjectR service](./projectr.md)
- RabbitMQ docs: <https://www.rabbitmq.com/documentation.html>

## Getting Help

- Check RabbitMQ logs: `/var/log/rabbitmq/`
- Review consumer logs via journalctl
- Test connection with simple producer/consumer scripts
- Create GitHub issue with "rabbitmq" label
