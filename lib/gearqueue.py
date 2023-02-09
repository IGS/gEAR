from datetime import datetime
import json
import os, sys
import functools
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

    def __init__(self, host="localhost", publisher_or_consumer=None, async_connection=False, pid=None):
        '''
        Establish a connection to RabbitMQ queue as a queue publisher or consumer

        Argument:
        ---------
        publisher_or_consumer = publisher - to publish (submit) jobs to queue
                                consumer - to consume (process) job from queue
        '''
        if not pid:
            pid = "[x]"
        self.pid=pid

        try:
            self.connection = gear.queue.RabbitMQQueue().connect(host, publisher_or_consumer=publisher_or_consumer, async_connection=async_connection)
        except:
            raise

    ### Connection context manager functions

    def __enter__(self):
        #return "entered MQ connection!"
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        #print("exited MQ connection!", file=sys.stderr)
        pass

    ### Channel publish/consume functions

    def open_channel(self):
        self.channel = self.connection.channel()

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
        self.channel.queue_declare(queue=queue_name)

        # Enabled delivery confirmations. This is REQUIRED.
        self.channel.confirm_delivery()

        # Send message (dataset ids) to job queue
        try:
            self.channel.basic_publish(exchange='message',
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
        self._consumer_tag = self.channel.basic_consume(queue=queue_name
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

        stream_fh = sys.stderr
        if queue_name:
            stream = '/var/log/gEAR_queue/{}.log'.format(queue_name)
            with open(stream, "a") as stream_fh:
                print("Queue {} ready".format(queue_name), file=stream_fh)

            try:
                print("{0}\tWaiting for messages. To exit press CTRL+C".format( str(datetime.now()) ), file=stream_fh)
                self.channel.start_consuming()
            except pika.exceptions.ConnectionClosed:
                print("{0}\tConnection to queue lost. Restarting consumer...".format( str(datetime.now()) ), file=stream_fh)
        return self


    def close(self):
        '''
        Closes queue connection
        '''
        self.connection.close()

class AsyncConnection(Connection):
    """Class to implement a callback-style SelectConnection.
    This is more performant and can handle situations such as unexpected failures."""

    def __init__(self, host="localhost", publisher_or_consumer=None, queue_name=None, on_message_callback=None,pid=None, logfile=None):
        super().__init__(host=host, publisher_or_consumer=publisher_or_consumer, async_connection=True, pid=pid)
        # At this point we should have self.channel defined

        self.log_fh = sys.stderr
        if logfile:
            self.log_fh = open(logfile, "a")

        self.queue = queue_name
        self.on_message = on_message_callback

        # Callback code snippets are from https://github.com/pika/pika/blob/main/examples/asynchronous_consumer_example.py
        # The example code makes a callback out of each individual step, which I think is a bit overkill (and does not mesh well
        # with our abstracted codebase), so I tried to eliminate some callbacks
        self.should_reconnect = False
        self.was_consuming = False
        self.exchange = 'message'

        self._closing = False
        self._consumer_tag = None
        self._consuming = False

        # No need to add on_open_error_callback since the connection is open ;-)
        self.connection.add_on_open_callback(self._on_connection_open)
        self.connection.add_on_close_callback(self._on_connection_closed)


    def __del__(self):
        try:
            # If logfile is not sys.stderr, close at this time.
            if not self.log_fh.closed:
                self.log_fh.close()
        except:
            pass

    ### Connection callback functions

    def _on_connection_open(self, _unused_connection):
        """Callback for when RabbitMQ connection opens."""
        print("{} - Connection opened".format(self.pid), flush=True, file=self.log_fh)
        self.open_channel()

    def _on_connection_closed(self, _unused_connection, reason):
        """Callback for when connection has closed unexpectedly."""
        self.connection.channel = None
        if self._closing:
            self.connection.connection.ioloop.stop()
        else:
            print("{} - Connection closed - {}".format(self.pid, reason), flush=True, file=self.log_fh)
            self.reconnect()

    ### Channel callback functions

    def _on_channel_open(self, channel):
        """Callback for when the channel opens."""
        self.channel = channel
        print("{} - Channel opened".format(self.pid), flush=True, file=self.log_fh)
        self.channel.add_on_close_callback(self._on_channel_closed)
        self.setup_exchange(self.exchange)

    def _on_channel_closed(self, channel, reason):
        """Callback for when the channel unexpectedly closes."""
        print("{} - Channel {} was closed - {}".format(self.pid, channel, reason), flush=True, file=self.log_fh)
        self.close_connection()

    def _on_consumer_cancelled(self, method_frame):
        """Invoked by pika when RabbitMQ sends a Basic.Cancel for a consumer
        receiving messages.
        :param pika.frame.Method method_frame: The Basic.Cancel frame
        """
        print("{} - Consumer was cancelled remotely, shutting down: {}".format(self.pid, method_frame), flush=True, file=self.log_fh)
        if self.channel:
            self.channel.close()

    ### SelectConnection helper functions
    def open_channel(self):
        """Open a new channel with RabbitMQ by issuing the Channel.Open RPC
        command. When RabbitMQ responds that the channel is open, the
        on_channel_open callback will be invoked by pika.
        """
        print("{} - Creating new channel".format(self.pid), flush=True, file=self.log_fh)
        self.connection.channel(on_open_callback=self._on_channel_open)

    def close_connection(self):
        if self.connection.is_closing or self.connection.is_closed:
            print("{} - Connection is closing or already closed".format(self.pid), flush=True, file=self.log_fh)
        else:
            print("{} - Connection closed".format(self.pid), flush=True, file=self.log_fh)
            self.close()

    def reconnect(self):
        """Will be invoked if the connection can't be opened or is
        closed. Indicates that a reconnect is necessary then stops the
        ioloop.
        """
        self.should_reconnect = True
        if not self._closing:
            self._closing = True
            print("{} - Stopping".format(self.pid), flush=True, file=self.log_fh)
            if self._consuming:
                self.stop_consuming()
                self.connection.ioloop.start()
            else:
                self.connection.ioloop.stop()
            print("{} - Stopped".format(self.pid), flush=True, file=self.log_fh)

    def setup_exchange(self, exchange_name):
        """Setup the exchange on RabbitMQ by invoking the Exchange.Declare RPC
        command. When it is complete, the on_exchange_declareok method will
        be invoked by pika.
        :param str|unicode exchange_name: The name of the exchange to declare
        """
        print("{} - Declaring exchange".format(self.pid, exchange_name), flush=True, file=self.log_fh)
        # Note: using functools.partial is not required, it is demonstrating
        # how arbitrary data can be passed to the callback when it is called
        cb = functools.partial(
            self.on_exchange_declareok, userdata=exchange_name)
        from pika.exchange_type import ExchangeType
        self.channel.exchange_declare(
            exchange=exchange_name,
            exchange_type=ExchangeType.topic,
            callback=cb)

    def on_exchange_declareok(self, _unused_frame, userdata):
        """Invoked by pika when RabbitMQ has finished the Exchange.Declare RPC
        command.
        :param pika.Frame.Method unused_frame: Exchange.DeclareOk response frame
        :param str|unicode userdata: Extra user data (exchange name)
        """
        print("{} - Exchange declared".format(self.pid, userdata), flush=True, file=self.log_fh)
        self.setup_queue(self.queue)

    def setup_queue(self, queue_name):
        """Setup the queue on RabbitMQ by invoking the Queue.Declare RPC
        command. When it is complete, the on_queue_declareok method will
        be invoked by pika.
        :param str|unicode queue_name: The name of the queue to declare.
        """
        print("{} - Declaring queue".format(self.pid, queue_name), flush=True, file=self.log_fh)
        cb = functools.partial(self.on_queue_declareok, userdata=queue_name)
        self.channel.queue_declare(queue=queue_name, callback=cb)

    def on_queue_declareok(self, _unused_frame, userdata):
        """Method invoked by pika when the Queue.Declare RPC call made in
        setup_queue has completed. In this method we will bind the queue
        and exchange together with the routing key by issuing the Queue.Bind
        RPC command. When this command is complete, the on_bindok method will
        be invoked by pika.
        :param pika.frame.Method _unused_frame: The Queue.DeclareOk frame
        :param str|unicode userdata: Extra user data (queue name)
        """
        queue_name = userdata
        print("{} - Binding {} to {}".format(self.pid, queue_name, self.exchange), flush=True, file=self.log_fh)
        cb = functools.partial(self.on_bindok, userdata=queue_name)
        self.channel.queue_bind(
            queue_name,
            self.exchange,
            callback=cb)

    def on_bindok(self, _unused_frame, userdata):
        """Invoked by pika when the Queue.Bind method has completed. At this
        point we will set the prefetch count for the channel.
        :param pika.frame.Method _unused_frame: The Queue.BindOk response frame
        :param str|unicode userdata: Extra user data (queue name)
        """
        print("{} - Queue bound {}".format(self.pid, userdata), flush=True, file=self.log_fh)
        self.set_qos()

    def set_qos(self):
        """This method sets up the consumer prefetch to only be delivered
        one message at a time. The consumer must acknowledge this message
        before RabbitMQ will deliver another one. You should experiment
        with different prefetch values to achieve desired performance.
        """
        self.channel.basic_qos(
            prefetch_count=1, callback=self.on_basic_qos_ok)

    def on_basic_qos_ok(self, _unused_frame):
        """Invoked by pika when the Basic.QoS method has completed. At this
        point we will start consuming messages by calling start_consuming
        which will invoke the needed RPC commands to start the process.
        :param pika.frame.Method _unused_frame: The Basic.QosOk response frame
        """
        print("{} - QOS set to {}".format(self.pid, "1"), flush=True, file=self.log_fh)
        self.start_consuming()

    def start_consuming(self):
        """Start consuming messages, setitng properties and adding some callbacks along the way."""
        # Information on why this is necessary at https://pika.readthedocs.io/en/stable/examples/connecting_async.html
        # Basically, this line allows the consumer to block consuming data to trigger callback actions
        #self.connection.ioloop.start()

        print("{} - Issuing consumer-related RPC commands".format(self.pid), flush=True, file=self.log_fh)
        self.channel.add_on_cancel_callback(self._on_consumer_cancelled)
        self._consumer_tag = self.channel.basic_consume(
            self.queue, self.on_message)
        self.was_consuming = True
        self._consuming = True

    def stop_consuming(self):
        """Tell RabbitMQ that you would like to stop consuming by sending the
        Basic.Cancel RPC command.
        """
        if self.channel:
            print("{} - Sending a Basic.Cancel RPC command to RabbitMQ".format(self.pid), flush=True, file=self.log_fh)
            cb = functools.partial(
                self.on_cancelok, userdata=self._consumer_tag)
            self.channel.basic_cancel(self._consumer_tag, cb)

    def on_cancelok(self, _unused_frame, userdata):
        """This method is invoked by pika when RabbitMQ acknowledges the
        cancellation of a consumer. At this point we will close the channel.
        This will invoke the on_channel_closed method once the channel has been
        closed, which will in-turn close the connection.
        :param pika.frame.Method _unused_frame: The Basic.CancelOk frame
        :param str|unicode userdata: Extra user data (consumer tag)
        """
        self._consuming = False
        print("{} - RabbitMQ acknowledged the cancellation of the consumer: {}".format(self.pid, userdata), flush=True, file=self.log_fh)
        print("{} - Closing channel".format(self.pid), flush=True, file=self.log_fh)
        self._channel.close()

    def run(self):
        """Run the example consumer by connecting to RabbitMQ and then
        starting the IOLoop to block and allow the SelectConnection to operate.
        """
        #self._connection = self.connect()
        # Information on why this is necessary at https://pika.readthedocs.io/en/stable/examples/connecting_async.html
        # Basically, this line allows the consumer to block consuming data to trigger callback actions
        self.connection.ioloop.start()

    def stop(self):
        """Cleanly shutdown the connection to RabbitMQ by stopping the consumer
        with RabbitMQ. When RabbitMQ confirms the cancellation, on_cancelok
        will be invoked by pika, which will then closing the channel and
        connection. The IOLoop is started again because this method is invoked
        when CTRL-C is pressed raising a KeyboardInterrupt exception. This
        exception stops the IOLoop which needs to be running for pika to
        communicate with RabbitMQ. All of the commands issued prior to starting
        the IOLoop will be buffered but not processed.
        """
        if not self._closing:
            self._closing = True
            print("{} - Stopping".format(self.pid), flush=True, file=self.log_fh)
            if self._consuming:
                self.stop_consuming()
                self._connection.ioloop.start()
            else:
                self._connection.ioloop.stop()
            print("{} - Stopped".format(self.pid), flush=True, file=self.log_fh)
