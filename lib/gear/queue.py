import os, sys
import pika


def get_connection(heartbeat_time=None, socket_time=None):
    try:
        # Connect to job queue
        # heartbeat_interval: http://stackoverflow.com/a/16155184/2900840
        connection = pika.BlockingConnection(
                                pika.ConnectionParameters(host='localhost',
                                heartbeat_interval=heartbeat_time,
                                socket_timeout=socket_time))
        return connection

    except Exception as err:
        print("Failed to connect to RabbitMQ:\n{0}".format(str(err)), file=sys.stderr)


class RabbitMQQueue:
    def __init__(self):
        connection = None

    def connect(self, publisher_or_consumer=None):
        '''
        Connect to RabbitMQ. Can connect as a publisher or consumer.
        Each connects with a different heartbeat interval and socket time.

        parameters:
            publisher_or_consumer =
                'publisher' = connect to RabbitMQ as publisher.
                                heartbeat = 0
                                socket time = 5
                'consumer'  = connect to RabbitMQ as consumer at MAX length of time
                                heartbeat = 21600
                                socket time = 21600
        '''
        if publisher_or_consumer is None:
            raise Exception("Must provide value 'publisher' or 'consumer' for 'publisher_or_consumer' argument")

        if publisher_or_consumer == 'publisher':
            heartbeat_time = 0
            socket_time = 5

        if publisher_or_consumer == 'consumer':
            heartbeat_time = 21600
            socket_time = 21600

        connection = get_connection(heartbeat_time=heartbeat_time, socket_time=socket_time)
        return connection
