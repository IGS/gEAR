# Setting up RabbitMQ

RabbitMQ is needed to act as a message broker, so that some of the load is taken off of the Flask instance when performing API requests.  In particular, this is relevant to the projectR API calls, where some of the dataset operations can be memory-intensive. By putting this responsibility outside of the Apache worker, we can hopefully prevent Apache from crashing, and better control the load of memory-intensive requests

To install RabbitMQ on your server, run the script here (and be sure to click the tab on this page for your specific OS version):

https://www.rabbitmq.com/docs/install-debian#apt-quick-start-cloudsmith

Test installation worked by checking `which rabbitmq-server`

## Using the messaging broker

There is a module at `<root>/lib/gearqueue.py` that contains a class to connect to the RabbitMQ messaging broker. This script uses a python module called "pika" under the hood.

## Creating a log file to view RabbitMQ logs

As root, I would create a file in `/var/log/gEAR_queue` named `<service>.log` where `<service>` is the name of the RabbitMQ consumer service (i.e. projectr). The `/var/log/gEAR_queue` directory is owned by root:adm with 750 permissions and the service log file within should be 644 permissions.

## Running a particular consumer

NOTE: This is automatically handled in the projectr_consumer system.d service file

First, make sure a directory is present under /var/log/gEAR_queue (you may have to create this as root). If you are not going to run the consumer listener as root, ensure the user has the same group-write privileges as the directory.

The consumer scripts are stored at `<root>/listeners/<files>`.  Let it run in the background (preferably with `nohup`)

Example script, run by root:
`sudo nohup /opt/bin/python3 ./listeners/projectr_consumer.py >>/var/log/gEAR_queue/projectr.log 1>/dev/null 2>>/var/log/gEAR_queue/projectr.log`

Executing a script multiple times will spawn off more workers.

## Purging a queue

Occasionally you may need to purge a queue, so that zombie jobs will not run and clog up the queue before the newer, actual jobs need to run.  To purge, run `sudo rabbitmqctl purge_queue <queue_name>`

## Making changes to the code

In most cases, the executing code is located in the callback function.  If this code is changed, the consumer daemon must be re-deployed.

## Preventing runaway memory

Both the projectR consumer and the anndata (H5AD/Seurat) upload consumer can be memory-intensive
on some inputs, so each guards against a duplicate/OOM-prone run of the same job:

- The projectR consumer refuses to run two instances of the *exact same* projection at once
  (`projectr_callback` locks on a `.lock` file next to the output CSV). If a client resubmits a
  job while the original run is still in progress (for example, clearing a stale-looking job
  status file and retrying), the resubmit detects the lock and returns a "running" status instead
  of starting a duplicate worker.
- The anndata upload consumer (`listeners/anndata_upload_consumer.py`) does the same for a
  duplicate/redelivered *queue message*: it acquires a non-blocking lock
  (`gear.utils.try_acquire_lock_file`) on a `.job.lock` file in the job's staging directory before
  processing, and if another worker already holds it, just acks and drops the duplicate delivery
  instead of starting a second conversion.
- The anndata upload consumer also self-imposes an `RLIMIT_AS` ceiling at process start
  (`set_memory_limit_from_cgroup()`, `lib/gear/utils.py`), sized to a fraction of the
  container/VM's cgroup memory limit, so an approaching OOM raises a catchable `MemoryError`
  instead of an uncatchable kernel `SIGKILL`. `process_uploaded_expression_dataset.cgi` sets the
  same guard for its synchronous fallback path (used when the queue is disabled or unreachable),
  since that runs the same processing inside an Apache CGI worker instead.

**Why the duplicate-processing guards above matter as much as the memory caps below:** if
`/etc/rabbitmq/rabbitmq.conf` doesn't set `consumer_timeout` (RabbitMQ's compiled-in default is 30
minutes), any job that legitimately runs longer than that gets its message redelivered to another
worker while the first is still processing it — neither knows about the other, nothing is freed,
and each successive redelivery adds its own full memory footprint on top of the others still in
flight. This has been observed in practice with large Seurat/RDS conversions (which can easily run
past 30 minutes): memory climbing in a staircase pattern as a second, then third worker picks up
the same job. **Set this explicitly** — create the file if it doesn't already exist:

```bash
sudo tee -a /etc/rabbitmq/rabbitmq.conf <<'EOF'
consumer_timeout = 43200000
EOF
sudo systemctl restart rabbitmq-server
```

(12 hours, comfortably above any realistic job duration; adjust down once real run durations are
known.) Confirm it took effect with
`sudo rabbitmqctl eval 'application:get_env(rabbit, consumer_timeout).'`. The per-job locks above
are still worth keeping as defense in depth (a worker crash or redeploy can cause a redelivery
too), but fixing the broker timeout removes the everyday trigger for it entirely.

As extra insurance against any worker consuming excessive memory (from these or any other cause),
add a per-VM memory cap via a systemd drop-in for each consumer family. **Size each family
independently** — they have very different memory profiles (a projectR chunk vs. a full Seurat/RDS
conversion, which has been observed to peak north of 20GB for a "only" ~3GB input file, since
`readRDS()` + `as_AnnData()` each hold a full copy in R's memory at once) — and account for every
consumer family running on the *same* VM when picking worker counts, not just one family in
isolation:

```bash
sudo mkdir -p /etc/systemd/system/projectr-consumer@.service.d
sudo tee /etc/systemd/system/projectr-consumer@.service.d/memory.conf <<'EOF'
[Service]
MemoryAccounting=true
MemoryHigh=14G
MemoryMax=18G
EOF

sudo mkdir -p /etc/systemd/system/anndata-upload-consumer@.service.d
sudo tee /etc/systemd/system/anndata-upload-consumer@.service.d/memory.conf <<'EOF'
[Service]
MemoryAccounting=true
MemoryHigh=24G
MemoryMax=28G
EOF

sudo systemctl daemon-reload
sudo systemctl restart projectr-consumer.target anndata-upload-consumer.target
```

The anndata figures above are close to the observed single-job peak (~22GB for a 3GB input), not a
guarantee — a larger upload can plausibly exceed `MemoryMax` and get killed on a *legitimate* single
run, not just a duplicate. Watch real peaks (`systemd-cgtop`, or cgroup `memory.peak`) on genuinely
large uploads and raise the cap if that happens, sized to the largest file the platform needs to
support. This lets a runaway worker get OOM-killed within its own cgroup (and restart, per each
`systemd/<unit>@.service`'s `StartLimitIntervalSec`/`StartLimitBurst`/`Restart=` settings) instead
of taking down the whole VM.

**These per-service caps stop one runaway job from taking down the whole VM — they don't reserve
headroom for every consumer family to run at full worker count simultaneously.** For example, 3
`projectr-consumer` workers at `MemoryMax=18G` (54G) plus 2 `anndata-upload-consumer` workers at
`MemoryMax=28G` (56G) already sum to well over a 61GB VM's capacity, on paper — that's fine as long
as they aren't all genuinely pegged at once (and two Seurat conversions landing simultaneously,
alone, can already approach the VM's full capacity), but if that becomes a regular pattern (heavy
uploads and heavy projections happening at the same time), or as spatial datasets grow large enough
that `spatial-upload-consumer` joins the same math, four independent per-service caps still don't
know about each other and can't stop the *combined* total from exceeding the VM. This is also why
`systemd/anndata-upload-consumer.target` defaults to fewer workers than projectR's 3 — Seurat/RDS
conversion's much larger per-job footprint doesn't leave room for the same worker count on a
shared VM.

### Aggregate cap across all consumer families: `gear-consumers.slice`

[systemd/gear-consumers.slice](../../../systemd/gear-consumers.slice) adds an aggregate ceiling on
top of the per-service caps above. Every templated consumer service (`projectr-consumer@`,
`anndata-upload-consumer@`, `spatial-upload-consumer@`, `gosling-upload-consumer@`) joins it via
`Slice=gear-consumers.slice` in its `[Service]` section. cgroup v2 then enforces the slice's
`MemoryHigh`/`MemoryMax` against the *combined* memory of every worker process across every family
in the slice, regardless of which mix happens to be busy at a given moment — something no number of
independent per-service caps can do on their own.

The two levels are complementary, not redundant: a per-service cap still catches one pathological
job before it alone can consume the *entire* shared budget and starve every other worker; the slice
cap is what catches the combined-usage case (e.g. several large, individually-legitimate jobs
landing across different families at once) that per-service caps structurally can't see. Deploy the
slice unit and reload before (or alongside) the per-service drop-ins:
```bash
sudo cp systemd/gear-consumers.slice /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl restart gear-consumers.slice
```
The `MemoryHigh=44G`/`MemoryMax=48G` in the checked-in unit are illustrative starting points for a
61GB VM (leaving headroom for Apache/MySQL/OS) — verify against the real non-consumer baseline
usage on the actual VM and adjust; not a final answer. Confirm it's active and grouping workers with
`systemctl status gear-consumers.slice` and `systemd-cgtop`.

## Troubleshooting

### (406, "PRECONDITION_FAILED - inequivalent arg 'durable' for queue 'projectr' in vhost '/': received 'false' but current is 'true'")

This error is probably popping up in the consumer. With this error, after attempting to start the consumer, you will probably see another error along the lines of `pika.exceptions.ChannelWrongStateError: Channel is closed` in the logs. This probably means that the queue was created in one context (durable=True) and is now attempted to be run in another context (durable=False).  Just run `sudo rabbitmqctl delete_queue <queue_name>` and then restart the consumer.
