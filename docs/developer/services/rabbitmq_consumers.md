# RabbitMQ Consumer Services

RabbitMQ consumers are background worker processes that handle asynchronous job processing for compute-intensive tasks.

## Overview

- **Technology**: RabbitMQ message broker, Python consumers
- **Location**: `listeners/` directory
- **Deployment**: Systemd services
- **Purpose**: Async processing of ProjectR, analysis, and other long-running tasks

## Architecture

```
Flask API
    ↓
RabbitMQ Queue
    ↓
Consumer Workers (systemd)
    ↓
Job Processing (ProjectR, Analysis, etc.)
    ↓
Result Storage (Database, File System)
```

## Available Consumers

### ProjectR Consumer

Processes matrix projection jobs for dimensionality reduction.

- **Listener**: `listeners/projectr_consumer.py`
- **Queue**: `projectr_jobs`
- **Service Template**: `systemd/projectr-consumer@.service`
- **Service Group**: `systemd/projectr-consumer.target`

### Additional Consumers

More consumers can be added following the same pattern:
- Analysis pipeline consumer
- Dataset processing consumer
- Export generation consumer

## Setup

### Prerequisites

1. **RabbitMQ Server**
   ```bash
   sudo apt install rabbitmq-server
   sudo systemctl enable rabbitmq-server
   sudo systemctl start rabbitmq-server
   ```
   
   See also: `docs/developer/setup/rabbitmq.md`

2. **Python Dependencies**
   ```bash
   pip install pika  # RabbitMQ Python client
   ```

3. **R and Packages** (for ProjectR)
   
   See: `docs/developer/setup/r_rpy2.md`

### Configuration

Configure in `gear.ini`:

```ini
[projectr_service]
;; 0 - disable RabbitMQ, 1 - enable. Disabling could lead to potential server crashes if many jobs are run simultaneously
queue_enabled = 0
queue_host = localhost
```

### Installing Services

```bash
# Copy service files to systemd directory
cd systemd
sudo cp projectr-consumer@.service /etc/systemd/system/
sudo cp projectr-consumer.target /etc/systemd/system/

# Reload systemd
sudo systemctl daemon-reload

# Enable services to start on boot
sudo systemctl enable projectr-consumer.target
```

## Starting Services

### Using Service Target (Recommended)

Start all consumers in the group:

```bash
# Start all ProjectR consumers
sudo systemctl start projectr-consumer.target

# Check status
sudo systemctl status projectr-consumer.target
```

### Individual Workers

Start specific numbered workers:

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
Description=ProjectR Consumer Worker %i
After=network.target rabbitmq-server.service

[Service]
Type=simple
User=www-data
WorkingDirectory=/var/www/gEAR
ExecStart=/usr/bin/python3 listeners/projectr_consumer.py
Restart=always
RestartSec=10

[Install]
WantedBy=projectr-consumer.target
```

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
Wants=projectr-consumer@1.service projectr-consumer@2.service

[Install]
WantedBy=multi-user.target
```

**Benefits:**
- Start/stop all workers with one command
- Manage workers as a group
- Ensure dependencies are met

## Writing Custom Consumers

### Basic Consumer Template

```python
#!/usr/bin/env python3
import pika
import json
import sys

def callback(ch, method, properties, body):
    """Process incoming job"""
    job_data = json.loads(body)
    
    # Process job
    result = process_job(job_data)
    
    # Acknowledge message
    ch.basic_ack(delivery_tag=method.delivery_tag)
    
    return result

def process_job(data):
    """Your job processing logic"""
    # Implement your processing here
    pass

def main():
    # Connect to RabbitMQ
    credentials = pika.PlainCredentials('guest', 'guest')
    parameters = pika.ConnectionParameters(
        host='localhost',
        credentials=credentials
    )
    
    connection = pika.BlockingConnection(parameters)
    channel = connection.channel()
    
    # Declare queue
    channel.queue_declare(queue='my_jobs', durable=True)
    
    # Set QoS (prefetch count)
    channel.basic_qos(prefetch_count=1)
    
    # Start consuming
    channel.basic_consume(
        queue='my_jobs',
        on_message_callback=callback
    )
    
    print('Waiting for messages...')
    channel.start_consuming()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted')
        sys.exit(0)
```

### Creating Service File

1. Create consumer script in `listeners/`
2. Create systemd service file in `systemd/`
3. Copy to `/etc/systemd/system/`
4. Start service

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
- Check RabbitMQ credentials
- Verify `use_rabbitmq = true`

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
sudo rabbitmqctl purge_queue projectr_jobs
```
**Warning:** This deletes all messages!

### Connection Errors

**RabbitMQ not accessible:**
```bash
# Check if RabbitMQ is listening
sudo netstat -tlnp | grep 5672

# Verify firewall rules
sudo ufw status
```

**Authentication errors:**
- Verify credentials in `gear.ini`
- Check RabbitMQ user permissions:
  ```bash
  sudo rabbitmqctl list_users
  sudo rabbitmqctl set_permissions -p / guest ".*" ".*" ".*"
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

### Updating Consumer Code

```bash
# Pull latest code
git pull origin devel

# Restart consumers
sudo systemctl restart projectr-consumer.target
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
sudo rabbitmqctl purge_queue projectr_jobs
```

## Related Documentation

- `docs/developer/setup/rabbitmq.md` - RabbitMQ setup
- `docs/developer/setup/r_rpy2.md` - R setup for ProjectR
- `systemd/README.md` - Systemd service management
- RabbitMQ docs: https://www.rabbitmq.com/documentation.html

## Getting Help

- Check RabbitMQ logs: `/var/log/rabbitmq/`
- Review consumer logs via journalctl
- Test connection with simple producer/consumer scripts
- Create GitHub issue with "rabbitmq" label
