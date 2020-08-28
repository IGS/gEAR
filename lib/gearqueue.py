from datetime import datetime
import json
import os, sys
import pika

import gear.queue

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
this.queue_connection = None
this.queue_channel = None


class Connection:
    '''
    Connect to RabbitMQ queue as a publisher or consumer. Publish job to and consume
    jobs from the queue.

    Examples from load_dataset queue:
        www/cgi/load_dataset_queue_publisher.cgi
        www/cgi/load_dataset_queue_consumer.cgi
    '''

    def __init__(self, publisher_or_consumer=None):
        '''
        Establish a connection to RabbitMQ queue as a queue publisher or consumer

        Argument:
        ---------
        publisher_or_consumer = publisher - to publish (submit) jobs to queue
                                consumer - to consume (process) job from queue
        '''
        if this.queue_connection is None:
            this.queue_connection = gear.queue.RabbitMQQueue().connect(publisher_or_consumer=publisher_or_consumer)

    def publish(self, queue_name=None, message=None):
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
        '''
        if queue_name is None:
            raise Exception("Error: Cannot publish. No queue_name given.")
        if message is None:
            raise Exception("Error: Cannot publish. No message given.")

        channel = this.queue_connection.channel()
        #Declare queue to use
        channel.queue_declare(queue=queue_name, durable=True)

        # Send message (dataset ids) to job queue
        channel.basic_publish(exchange='',
                            routing_key=queue_name,
                            body=json.dumps(message),
                            properties=pika.BasicProperties(
                                delivery_mode = 2 # make message persistent
                            ))

        channel.close()


    def get_consumer_channel(self):
        '''
        Returns a RabbitMQ connection channel.

        This way methods consume and continue_consuming can run on the same channel.
        '''
        channel = None
        if this.queue_channel is None:
            channel = this.queue_connection.channel()
        else:
            channel = this.queue_channel

        return channel


    def consume(self, callback=None, queue_name=None, channel=None):
        '''
        With a consumer connection and channel, retrieves 1 job from a named queue and
        processes it.

        Arguments:
        ----------
        callback = Function used to process job.
                Example: See function 'callback' in www/cgi/load_dataset_queue_consumer.cgi
        queue_name = Name of the queue to get job from.
                Example: 'load_dataset'
        channel = Consumer channel opened from Connection.get_consumer_channel()
        '''
        if callback is None:
            raise Exception("Error: Cannot consume. No callback function given.")
        if channel is None:
            raise Exception("Error: Cannot consume. No channel given.")
        if queue_name is None:
            raise Exception("Error: Cannot consume. No queue_name given.")

        this.queue_connection.process_data_events(time_limit=0)

        #Declare queue to use
        channel.queue_declare(queue=queue_name, durable=True)
        channel.basic_qos(prefetch_count=1)

        # Initiate callback
        channel.basic_consume(callback,
                              queue=queue_name)

        return self

    def continue_consuming(self, channel=None):
        '''
        Continues to consume jobs as they appear on queue.

        If the queue broker drops the connection, restart the 'load_dataset' consumer.
        For large datasets with a long processing time, the broker may close the
        connection despite setting a long heartbeat and socket timeout when making
        the initial consumer to queue connection.

        Argument:
        ---------
        channel = Consumer channel opened from Connection.get_consumer_channel()
        '''
        if channel is None:
            raise Exception("Error: Cannot continue consuming. No channel given.")

        try:
            print("{0}\tWaiting for messages. To exit press CTRL+C".format( str(datetime.now()) ), file=open('/var/log/gEAR_queue/load_dataset.log', 'a'))
            channel.start_consuming()
        except pika.exceptions.ConnectionClosed:
            print("{0}\tConnection to queue lost. Restarting consumer...".format( str(datetime.now()) ), file=open('/var/log/gEAR_queue/load_dataset.log', 'a'))
            os.system('./cgi/load_dataset_queue_consumer.cgi')

        return self


    def close(self):
        '''
        Closes queue connection
        '''
        this.queue_connection.close()
