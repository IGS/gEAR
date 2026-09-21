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

- The projectR consumer already refuses to run two instances of the *exact same* projection at
  once (`projectr_callback` locks on a `.lock` file next to the output CSV). If a client resubmits
  a job while the original run is still in progress (for example, clearing a stale-looking job
  status file and retrying), the resubmit now detects the lock and returns a "running" status
  instead of starting a duplicate worker.
- The anndata upload consumer (`listeners/anndata_upload_consumer.py`) self-imposes an `RLIMIT_AS`
  ceiling at process start (`set_memory_limit_from_cgroup()`, `lib/gear/utils.py`), sized to a
  fraction of the container/VM's cgroup memory limit, so an approaching OOM raises a catchable
  `MemoryError` instead of an uncatchable kernel `SIGKILL`. `process_uploaded_expression_dataset.cgi`
  sets the same guard for its synchronous fallback path (used when the queue is disabled or
  unreachable), since that runs the same processing inside an Apache CGI worker instead.

As extra insurance against any worker consuming excessive memory (from these or any other cause),
add a per-VM memory cap via a systemd drop-in for each consumer family, sized to the VM's RAM and
worker count (e.g. for a 61 GB VM running 3 workers of each):

```bash
for unit in projectr-consumer anndata-upload-consumer; do
  sudo mkdir -p "/etc/systemd/system/${unit}@.service.d"
  sudo tee "/etc/systemd/system/${unit}@.service.d/memory.conf" <<'EOF'
[Service]
MemoryAccounting=true
MemoryHigh=14G
MemoryMax=18G
EOF
done
sudo systemctl daemon-reload
sudo systemctl restart projectr-consumer.target anndata-upload-consumer.target
```

This lets a runaway worker get OOM-killed within its own cgroup (and restart, per each
`systemd/<unit>@.service`'s `StartLimitIntervalSec`/`StartLimitBurst`/`Restart=` settings) instead
of taking down the whole VM.

It's also worth setting an explicit `consumer_timeout` in `/etc/rabbitmq/rabbitmq.conf` (RabbitMQ
defaults to 30 minutes) if any projectR jobs are expected to legitimately run longer than that —
otherwise the broker will close the channel and requeue the message to another worker while the
first is still running, which is exactly the kind of duplicate-processing scenario described
above:

```
consumer_timeout = 43200000
```

## Troubleshooting

### (406, "PRECONDITION_FAILED - inequivalent arg 'durable' for queue 'projectr' in vhost '/': received 'false' but current is 'true'")

This error is probably popping up in the consumer. With this error, after attempting to start the consumer, you will probably see another error along the lines of `pika.exceptions.ChannelWrongStateError: Channel is closed` in the logs. This probably means that the queue was created in one context (durable=True) and is now attempted to be run in another context (durable=False).  Just run `sudo rabbitmqctl delete_queue <queue_name>` and then restart the consumer.
