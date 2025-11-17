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

## Troubleshooting

### (406, "PRECONDITION_FAILED - inequivalent arg 'durable' for queue 'projectr' in vhost '/': received 'false' but current is 'true'")

This error is probably popping up in the consumer. With this error, after attempting to start the consumer, you will probably see another error along the lines of `pika.exceptions.ChannelWrongStateError: Channel is closed` in the logs. This probably means that the queue was created in one context (durable=True) and is now attempted to be run in another context (durable=False).  Just run `sudo rabbitmqctl delete_queue <queue_name>` and then restart the consumer.
