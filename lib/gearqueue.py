from datetime import datetime
import json
import os, sys
import pika

import gear.queue

# SAdkins - Originally this module created module global variables for the connection and a single channel
# However, this conection was being shared among each of the Flask API requests, so if one threw an error, they all did

class Connection:
    '''
    Connect to RabbitMQ queue as a publisher or consumer. Publish job to and consume
    jobs from the queue.

    Examples from load_dataset queue:
        www/cgi/load_dataset_queue_publisher.cgi
        www/cgi/load_dataset_queue_consumer.cgi
    '''

    def __init__(self, host="localhost", publisher_or_consumer=None):
        '''
        Establish a connection to RabbitMQ queue as a queue publisher or consumer

        Argument:
        ---------
        publisher_or_consumer = publisher - to publish (submit) jobs to queue
                                consumer - to consume (process) job from queue
        '''
        try:
            self.connection = gear.queue.RabbitMQQueue().connect(host, publisher_or_consumer=publisher_or_consumer)
        except:
            raise
        self.channel = self.connection.channel()


    # properties to make this a context manager
    def __enter__(self):
        return "entered MQ connection!"
    def __exit__(self, exc_type, exc_value, traceback):
        print("exited MQ connection!", file=sys.stderr)

    def publish(self, queue_name=None, message=None, **kwargs):
        '''
        Submits a job to queue.

        1) Opens a channel from queue connection.
        2) Declares queue. Only queue at the moment is 'load_dataset'
        3) Submits (publishes) job to queue

        Arguments:
        ----------
        queue_name = string value of the queue being connected to.
                Example: 'load_dataset'
        message = python dict of key values.
                Example:
                    {'uploaded_expression_name': 'base_template.xlsx',
                    'uploaded_metadata_name': 'test_metadata_v3.xlsx',
                    'session_id': 'ada703dc-c38a-40fd-9f62-fd34fff2d8be',
                    'dataset_uid': '35f5b227-59e9-75bb-dd41-f26fff691506',
                    'share_uid': '6c7e7775'
                    }
        kwargs - keys to pass to pike.BasicProperties
        '''
        if queue_name is None:
            raise Exception("Error: Cannot publish. No queue_name given.")
        if message is None:
            raise Exception("Error: Cannot publish. No message given.")

        # Declare queue to use
        self.channel.queue_declare(queue=queue_name, durable=True)

        # Enabled delivery confirmations. This is REQUIRED.
        self.channel.confirm_delivery()

        # Send message (dataset ids) to job queue
        try:
            self.channel.basic_publish(exchange='',
                                routing_key=queue_name,
                                body=json.dumps(message),
                                properties=pika.BasicProperties(
                                    delivery_mode = 2 # make message persistent (disk instead of memory)
                                    , content_type="application/json"
                                    , **kwargs
                                ))
        except pika.exceptions.UnroutableError:
            print('Message was returned')
            raise

        self.connection.process_data_events(time_limit=None)

    def consume(self, queue_name=None, on_message_callback=None, num_messages=1, auto_ack=False, skip_queue_declare=False):
        '''
        With a consumer connection and channel, retrieves 1 job from a named queue and
        processes it.

        Arguments:
        ----------
        queue_name = Name of the queue to get job from.
                Example: 'load_dataset'
        on_message_callback =
            The function to call when consuming with the signature on_message_callback(channel, method, properties, body), where
                channel: pika.channel.Channel (this channel)
                method: pika.spec.Basic.Deliver
                properties: pika.spec.BasicProperties
                body: bytes
        num_messages = Number of simultaneous unacknowledged messages to received at once
        auto_ack = auto-acknowledge the reply
        skip_queue_declare = if True, skip because the queue was already declared

        '''
        if on_message_callback is None:
            raise Exception("Error: Cannot consume. No callback function given.")
        if queue_name is None:
            raise Exception("Error: Cannot consume. No queue_name given.")

        #Declare queue to use
        if not skip_queue_declare:
            self.channel.queue_declare(queue=queue_name, durable=True)
        self.channel.basic_qos(prefetch_count=num_messages) # balances worker load

        # Initiate callback
        self.channel.basic_consume(queue=queue_name
                            , on_message_callback=on_message_callback
                            , auto_ack=auto_ack
                            )

        #self.connection.process_data_events(time_limit=0)

        return self

    def replyto_consume(self, on_message_callback):
        # Sets up a direct RPC (request/reply) consumer
        # See https://www.rabbitmq.com/direct-reply-to.html
        # See https://pika.readthedocs.io/en/stable/examples/direct_reply_to.html?highlight=reply_to#direct-reply-to-example

        self.consume(
            queue_name="amq.rabbitmq.reply-to"
            , on_message_callback=on_message_callback
            , auto_ack=True
            , skip_queue_declare=True
            )

        return self

    def continue_consuming(self, queue_name=None):
        '''
        Continues to consume jobs as they appear on queue.

        If the queue broker drops the connection, restart the 'load_dataset' consumer.
        For large datasets with a long processing time, the broker may close the
        connection despite setting a long heartbeat and socket timeout when making
        the initial consumer to queue connection.

        Argument:
        ---------
        queue_name = name of queue to write logs to.
        '''

        stream = sys.stderr
        if queue_name:
            stream = '/var/log/gEAR_queue/{}.log'.format(queue_name)
            print("Queue {} ready".format(queue_name), file=open(stream, 'a'))

        try:
            print("{0}\tWaiting for messages. To exit press CTRL+C".format( str(datetime.now()) ), file=open(stream, 'a'))
            self.channel.start_consuming()
        except pika.exceptions.ConnectionClosed:
            print("{0}\tConnection to queue lost. Restarting consumer...".format( str(datetime.now()) ), file=open(stream, 'a'))
        return self


    def close(self):
        '''
        Closes queue connection
        '''
        self.connection.close()
