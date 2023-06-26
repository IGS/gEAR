# Setting up RabbitMQ

RabbitMQ is needed to act as a message broker, so that some of the load is taken off of the Flask instance when performing API requests.  In particular, this is relevant to the projectR API calls, where some of the dataset operations can be memory-intensive. By putting this responsibility outside of the Apache worker, we can hopefully prevent Apache from crashing, and better control the load of memory-intensive requests

To install RabbitMQ on your server, run this script (taken from https://www.rabbitmq.com/install-debian.html#apt-quick-start-cloudsmith).  Do note that the script was outdated so I updated "bionic" to "jammy" to reflect the version of Ubuntu being used.

```bash
#!/usr/bin/sh

sudo apt-get install curl gnupg apt-transport-https -y

## Team RabbitMQ's main signing key
curl -1sLf "https://keys.openpgp.org/vks/v1/by-fingerprint/0A9AF2115F4687BD29803A206B73A36E6026DFCA" | sudo gpg --dearmor | sudo tee /usr/share/keyrings/com.rabbitmq.team.gpg > /dev/null
## Cloudsmith: modern Erlang repository
curl -1sLf https://dl.cloudsmith.io/public/rabbitmq/rabbitmq-erlang/gpg.E495BB49CC4BBE5B.key | sudo gpg --dearmor | sudo tee /usr/share/keyrings/io.cloudsmith.rabbitmq.E495BB49CC4BBE5B.gpg > /dev/null
## Cloudsmith: RabbitMQ repository
curl -1sLf https://dl.cloudsmith.io/public/rabbitmq/rabbitmq-server/gpg.9F4587F226208342.key | sudo gpg --dearmor | sudo tee /usr/share/keyrings/io.cloudsmith.rabbitmq.9F4587F226208342.gpg > /dev/null

## Add apt repositories maintained by Team RabbitMQ
sudo tee /etc/apt/sources.list.d/rabbitmq.list <<EOF
## Provides modern Erlang/OTP releases
##
deb [signed-by=/usr/share/keyrings/io.cloudsmith.rabbitmq.E495BB49CC4BBE5B.gpg] https://dl.cloudsmith.io/public/rabbitmq/rabbitmq-erlang/deb/ubuntu jammy main
deb-src [signed-by=/usr/share/keyrings/io.cloudsmith.rabbitmq.E495BB49CC4BBE5B.gpg] https://dl.cloudsmith.io/public/rabbitmq/rabbitmq-erlang/deb/ubuntu jammy main

## Provides RabbitMQ
##
deb [signed-by=/usr/share/keyrings/io.cloudsmith.rabbitmq.9F4587F226208342.gpg] https://dl.cloudsmith.io/public/rabbitmq/rabbitmq-server/deb/ubuntu jammy main
deb-src [signed-by=/usr/share/keyrings/io.cloudsmith.rabbitmq.9F4587F226208342.gpg] https://dl.cloudsmith.io/public/rabbitmq/rabbitmq-server/deb/ubuntu jammy main
EOF

## Update package indices
sudo apt-get update -y

## Install Erlang packages
sudo apt-get install -y erlang-base \
                        erlang-asn1 erlang-crypto erlang-eldap erlang-ftp erlang-inets \
                        erlang-mnesia erlang-os-mon erlang-parsetools erlang-public-key \
                        erlang-runtime-tools erlang-snmp erlang-ssl \
                        erlang-syntax-tools erlang-tftp erlang-tools erlang-xmerl

## Install rabbitmq-server and its dependencies
sudo apt-get install rabbitmq-server -y --fix-missing
```

Test installation worked by checking `which rabbitmq-server`

## Using the messaging broker

There is a module at `<root>/lib/gearqueue.py` that contains a class to connect to the RabbitMQ messaging broker. This script uses a python module called "pika" under the hood.

## Creating a log file to view RabbitMQ logs

As root, I would create a file in `/var/log/gEAR_queue` named `<service>.log` where `<service>` is the name of the RabbitMQ consumer service (i.e. projectr). The `/var/log/gEAR_queue` directory is owned by root:adm with 750 permissions and the service log file within should be 644 permissions.

## Running a particular consumer

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

### Finding the errors

Typically we are using a callback to run the process after the rabbitmq consumer receives the message and payload. It is important to remember that the callback will be under the rabbitMQ process if that is enabled, and so printing to sys.stderr will not send to Apache logs like typical API/CGI calls do. If you wish, you can pass the filehandle from the consumer as an argument to the callback and print to that so those messages go to your consumer log (and print to sys.stderr if you are not using RabbitMQ).  See the ProjectR consumer and callback for an example.

### (406, "PRECONDITION_FAILED - inequivalent arg 'durable' for queue 'projectr' in vhost '/': received 'false' but current is 'true'")

This error is probably popping up in the consumer. With this error, after attempting to start the consumer, you will probably see another error along the lines of `pika.exceptions.ChannelWrongStateError: Channel is closed` in the logs. This probably means that the queue was created in one context (durable=True) and is now attempted to be run in another context (durable=False).  Just run `sudo rabbitmqctl delete_queue <queue_name>` and then restart the consumer.
