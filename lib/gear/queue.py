import os, sys
import pika


def get_connection(host="localhost", heartbeat_time=None, socket_time=None):
    try:
        # Connect to job queue
        # heartbeat_interval: http://stackoverflow.com/a/16155184/2900840
        # Assumes default port, user, and password
        connection = pika.BlockingConnection(
                                pika.ConnectionParameters(host=host,
                                heartbeat=heartbeat_time,
                                socket_timeout=socket_time))
        return connection

    except Exception as err:
        print("Failed to connect to RabbitMQ:\n{0}".format(str(err)), file=sys.stderr)
        raise

def get_async_connection(host="localhost", heartbeat_time=None, socket_time=None):
    """Create and return an asynchronous connection."""
    try:
        connection = pika.SelectConnection(
                            pika.ConnectionParameters(host=host,
                                heartbeat=heartbeat_time,
                                socket_timeout=socket_time)
                            )
        return connection

    except Exception as err:
        print("Failed to connect to RabbitMQ:\n{0}".format(str(err)), file=sys.stderr)
        raise

class RabbitMQQueue:
    def __init__(self):
        connection = None

    def connect(self, host="localhost", publisher_or_consumer=None, async_connection=False):
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

        # Having "publisher_or_consumer=None" is valid, especially with reply-to connections.
        # So default to "consumer"
        if not publisher_or_consumer:
            publisher_or_consumer = "consumer"

        if publisher_or_consumer == 'publisher':
            heartbeat_time = 0
            socket_time = 5

        if publisher_or_consumer == 'consumer':
            heartbeat_time = 21600
            socket_time = 21600

        if async_connection:
            return get_async_connection(host, heartbeat_time=heartbeat_time, socket_time=socket_time)
        return get_connection(host, heartbeat_time=heartbeat_time, socket_time=socket_time)


