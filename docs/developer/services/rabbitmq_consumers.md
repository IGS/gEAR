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

## How the Units Fit Together

The unit files live in `systemd/`. They form three levels, with every worker in one shared memory slice:

```text
gear-consumers.target                      (everything)
├── anndata-upload-consumer.target         (one family)
│   ├── anndata-upload-consumer@1.service  (one worker process)
│   └── anndata-upload-consumer@2.service
├── gosling-upload-consumer.target
│   └── gosling-upload-consumer@1..3.service
├── spatial-upload-consumer.target
│   └── spatial-upload-consumer@1..2.service
└── projectr-consumer.target
    └── projectr-consumer@1..3.service

gear-consumers.slice  -- every @N.service runs in it (Slice=), giving an aggregate memory cap
```

- **`<name>@.service`** is a template. The `@` means systemd can make any number of instances from it: `projectr-consumer@1.service`, `projectr-consumer@2.service`, and so on. `%i` in the file is replaced with the instance number. Each instance is one worker process.
- **`<name>.target`** is a consumer family. Its `Wants=` line lists which numbered instances to run, so the worker count is set there.
- **`gear-consumers.target`** groups every family target.
- **`gear-consumers.slice`** sets an aggregate memory limit across all workers. See [RabbitMQ setup](../setup/rabbitmq.md#preventing-runaway-memory).

### Why `Wants=` and `PartOf=` both matter

The two settings link the levels in opposite directions, and each controls different commands:

- **`Wants=`** goes on the parent (`gear-consumers.target` wants each family target; a family target wants its workers). It is used **only for start**. Starting a parent starts any wanted unit that is not already running.
- **`PartOf=`** goes on the child (each `@.service` is `PartOf=` its family target; each family target is `PartOf=gear-consumers.target`). It is used **only for stop and restart**. Stopping or restarting a parent also stops or restarts the children that are `PartOf` it.

Both links are needed at every level. Without `PartOf=`, `systemctl restart gear-consumers.target` restarts only the target itself, and the workers keep running the old code.

| Command | Effect on workers |
| --- | --- |
| `sudo systemctl start gear-consumers.target` | Starts any worker that is stopped. **Running workers are left alone**, so they do not pick up new code. |
| `sudo systemctl restart gear-consumers.target` | Restarts every worker. Use this to load new code. |
| `sudo systemctl stop gear-consumers.target` | Stops every worker. |
| `sudo systemctl restart projectr-consumer.target` | Restarts only that family's workers. |
| `sudo systemctl restart projectr-consumer@2.service` | Restarts one worker. |

> Older installed copies of the unit files may not have these `PartOf=` lines. Check with `grep PartOf /etc/systemd/system/*consumer*`. Without them, stop and restart of `gear-consumers.target` (and of `projectr-consumer.target`) do not reach the workers. Copy in the current unit files once, as described in [Changing unit files](#changing-unit-files).

When adding a new consumer, add its target to `Wants=` in `gear-consumers.target` **and** give the new target `PartOf=gear-consumers.target`. See [Creating Service File](#creating-service-file).

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

This is a one-time step on a new server. For later changes to the unit files, see [Changing unit files](#changing-unit-files).

```bash
cd /path/to/gEAR/systemd
sudo cp *-consumer@.service *-consumer.target gear-consumers.target gear-consumers.slice /etc/systemd/system/

# Replace <gear_root> with the gEAR checkout path on this server, and check the
# interpreter path (/opt/bin/python3) in each *@.service file.
sudo sed -i 's|<gear_root>|/path/to/gEAR|g' /etc/systemd/system/*-consumer@.service

# Size MemoryHigh/MemoryMax in gear-consumers.slice for the VM (see ../setup/rabbitmq.md).

sudo systemctl daemon-reload

# Start on boot, then start now. Enabling gear-consumers.target is enough: it
# pulls in every family target, which pulls in every worker.
sudo systemctl enable gear-consumers.target
sudo systemctl start gear-consumers.target
```

The log directory `/var/log/gEAR_queue/` must exist and be writable by `www-data` before the workers start. See [RabbitMQ setup](../setup/rabbitmq.md).

## Starting and Stopping

Use `gear-consumers.target` for everything, a family target for one consumer type, or a numbered instance for one worker. See the table in [Why `Wants=` and `PartOf=` both matter](#why-wants-and-partof-both-matter) for what each command does.

```bash
# All consumers
sudo systemctl start gear-consumers.target
sudo systemctl restart gear-consumers.target
sudo systemctl stop gear-consumers.target

# One family
sudo systemctl restart anndata-upload-consumer.target

# One worker
sudo systemctl restart projectr-consumer@2.service
```

`start` never restarts a worker that is already running. To load new code, use `restart`.

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

## Rolling Out Updates

Use this section when deploying new code to a server. There are two kinds of change, and they need different steps:

- **Code changes** in `listeners/`, `lib/`, or API code the consumers import (for example `www/api/resources/projectr.py`). The workers only need a restart.
- **Unit file changes** in `systemd/`. These need the files copied into `/etc/systemd/system/` as well. See [Changing unit files](#changing-unit-files).

### 1. Find out what changed

```bash
cd /path/to/gEAR
git fetch
git diff --stat HEAD origin/devel -- listeners/ lib/ www/api/ systemd/
```

If nothing under `systemd/` changed, follow the code-only steps below. If it did, go to [Changing unit files](#changing-unit-files) instead.

### 2. Check for jobs in progress

A restart kills any job a worker is running at that moment. Workers acknowledge a job only after finishing it, and they do not shut down gracefully: systemd sends SIGTERM and the worker exits right away. RabbitMQ then puts the unacknowledged job back in the queue, and a worker runs it again **from the start** once it comes back up.

For the anndata and spatial uploads, each run counts as an attempt (a `.attempt_count` file in the upload's staging directory, two attempts allowed). An upload that is killed by a restart has used one of its two attempts.

Before restarting, check whether any worker is busy:

```bash
sudo rabbitmqctl list_queues name messages_ready messages_unacknowledged consumers
```

- `messages_unacknowledged` is the number of jobs being worked on right now. If it is 0 for a queue, that family can be restarted safely.
- `messages_ready` is the number of jobs waiting. These are safe; the queues are durable, so waiting jobs stay in the queue during a restart and are picked up afterwards.

If a queue shows jobs in progress, either wait for them to finish, or restart only the idle families now (for example `sudo systemctl restart gosling-upload-consumer.target`) and come back for the busy ones. Large anndata uploads and projectR jobs can take a long time, so check the logs in `/var/log/gEAR_queue/` to see what is running.

### 3. Pull and restart (code changes only)

```bash
cd /path/to/gEAR
git pull origin devel
sudo systemctl restart gear-consumers.target
```

No `daemon-reload` is needed when only Python code changed. `daemon-reload` is only for changes to the unit files.

### 4. Verify

```bash
# Every target and worker should be "active". "Active: ... since" should show the restart time.
systemctl list-units --all 'gear-consumers.target' '*-consumer.target' '*-consumer@*'
systemctl status '*-consumer@*' | grep -E 'consumer@.*service -|Active:'

# Each queue should have as many consumers as its target has workers (see Queue Summary above)
sudo rabbitmqctl list_queues name consumers

# Look for start-up errors
sudo tail -n 20 /var/log/gEAR_queue/*.log
sudo journalctl -u '*-consumer@*' --since '5 min ago' -p warning
```

### Changing Unit Files

`sudo systemctl daemon-reload` reads the unit files again but does **not** change any worker that is already running. A running worker keeps the settings it started with (`User=`, `ExecStart=`, `Environment=`, `Slice=`, and so on) until it is restarted, and `start` does nothing to a worker that is already running. So after a unit file change, every affected worker must be stopped and started again.

Use this sequence:

```bash
# 1. Check for jobs in progress (step 2 above). Then stop every worker while the
#    old unit files are still installed.
sudo systemctl stop gear-consumers.target

# 2. Copy in the new unit files and fill in <gear_root>
cd /path/to/gEAR/systemd
sudo cp *-consumer@.service *-consumer.target gear-consumers.target gear-consumers.slice /etc/systemd/system/
sudo sed -i 's|<gear_root>|/path/to/gEAR|g' /etc/systemd/system/*-consumer@.service

# 3. Load the new unit files
sudo systemctl daemon-reload

# 4. Make any fixes the change needs (see the table below)

# 5. Start every worker with the new settings, then verify (step 4 above)
sudo systemctl start gear-consumers.target
```

Stopping **before** copying the files matters for two reasons:

- systemd still has the old definitions loaded, so the stop reaches every running worker, including any instance or unit that the new files rename or remove. If you remove a unit file first, systemd lists the unit as `not-found` but the worker keeps running.
- Nothing runs during the gap, so you can fix things like file ownership before any worker starts with the new settings.

> **First time only:** if the installed unit files predate the `PartOf=` lines (`grep PartOf /etc/systemd/system/*consumer*` shows no target lines), `stop gear-consumers.target` does not reach the workers. In step 1, stop them directly instead: `sudo systemctl stop gear-consumers.target '*-consumer.target' '*-consumer@*'`.

What each kind of change needs:

| Change | What to do |
| --- | --- |
| Edit an existing setting (`ExecStart=`, `Environment=`, `Restart*=`, `StartLimit*=`) | The sequence above. A plain `daemon-reload` + `restart gear-consumers.target` also works. |
| Change `User=` or `Group=` (for example root to www-data) | The sequence above, and in step 4 give the new user ownership of the files the old user created (see below). Otherwise the workers can't write their logs and fail to start. |
| Rename or delete a unit file, or remove an instance from a target's `Wants=` | The sequence above. Stopping first is required here. Afterwards, run `sudo systemctl reset-failed` to clear leftover entries. |
| Add an instance to a target's `Wants=` | Copy the files, `daemon-reload`, then `start` the target. Only the new instance starts. See [Scaling Workers](#scaling-workers). |
| Change memory limits in `gear-consumers.slice` | Copy the file and `daemon-reload`. Check the limits with `systemctl show gear-consumers.slice -p MemoryHigh -p MemoryMax`. |
| Change the `[Install]` section (`WantedBy=`) | Copy the files, `daemon-reload`, then `sudo systemctl reenable <unit>`. |

Per-server memory overrides in drop-in directories (`/etc/systemd/system/<family>@.service.d/memory.conf`, see [RabbitMQ setup](../setup/rabbitmq.md#preventing-runaway-memory)) are separate files, so copying in new base unit files does not overwrite them.

#### Changing the user a consumer runs as

Files created while a worker ran as root stay owned by root, and a worker running as `www-data` can't write to them. The main ones are the consumer logs and the per-upload state files in each upload staging directory (`www/uploads/files/<session>/<share_uid>/`: `.job.lock`, `.attempt_count`, `status.json`). In step 4 of the sequence above:

```bash
sudo chown www-data:www-data /var/log/gEAR_queue/*.log
sudo find /path/to/gEAR/www/uploads/files -user root -exec chown www-data:www-data {} +
```

To confirm which user each worker really runs as, and when it started:

```bash
ps -o user,pid,lstart,args -C python3 | grep listeners/
```

### Scaling Workers

The number of workers in a family is set by the `Wants=` line in its target, for example `systemd/projectr-consumer.target`. Check memory before adding workers; see [Worker Count](#worker-count) and [RabbitMQ setup](../setup/rabbitmq.md#preventing-runaway-memory).

To **add** a worker, add the instance (for example `projectr-consumer@4.service`) to `Wants=`, copy the target file into `/etc/systemd/system/`, run `sudo systemctl daemon-reload`, then `sudo systemctl start projectr-consumer.target`. Only the new instance starts.

To **remove** a worker, stop it first, then update the target:

```bash
sudo systemctl stop projectr-consumer@3.service
# remove projectr-consumer@3.service from Wants= in the target, copy it in, then:
sudo systemctl daemon-reload
```

A removed instance is still `PartOf` the target, so if it is left running it keeps running, and it is restarted along with the target.

To run an extra worker temporarily without editing any file, start the instance directly: `sudo systemctl start projectr-consumer@4.service`. It stops at the next reboot, or with `stop`.

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
2. Add `systemd/<name>@.service` and `systemd/<name>.target` (copy an existing pair). Keep `Slice=gear-consumers.slice` and `PartOf=<name>.target` in the service, and `PartOf=gear-consumers.target` in the target
3. Add the new target to `Wants=` in `systemd/gear-consumers.target`
4. Optionally add `listeners/Dockerfile.<name>`, a `docker/docker-bake.hcl` target, and a compose service
5. Copy the unit files to `/etc/systemd/system/`, fill in `<gear_root>`, `daemon-reload`, and start the target (see [Changing Unit Files](#changing-unit-files))

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

### Worker Shows `failed` with `start-limit-hit`

`systemctl status` shows `Failed with result 'start-limit-hit'`. The worker crashed 5 times within 300 seconds (`StartLimitBurst=` / `StartLimitIntervalSec=`), so systemd stopped restarting it. `start` is refused until the failure is cleared. Find and fix the cause in the logs first, then:

```bash
sudo journalctl -u projectr-consumer@1.service -n 100
sudo systemctl reset-failed 'projectr-consumer@*'
sudo systemctl start projectr-consumer.target
```

### Worker Still Running Old Code or Settings

`start` does not restart a worker that is already running, and `daemon-reload` does not change a running worker. Use `restart` after code changes, and follow [Changing Unit Files](#changing-unit-files) after unit file changes. To see when each worker started and which user it runs as:

```bash
ps -o user,pid,lstart,args -C python3 | grep listeners/
```

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

For deploying new code or unit files, see [Rolling Out Updates](#rolling-out-updates).

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
