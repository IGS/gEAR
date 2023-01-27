#!/opt/bin/python3

"""
projectr_consumer.py - RabbitMQ messaging consumer
"""

import os, sys, json
import gc

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import gearqueue

www_path = os.path.abspath(os.path.join('..', 'www'))
sys.path.append(www_path)

from api.resources.projectr import projectr_callback

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

queue_name = "projectr"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
stream = '/var/log/gEAR_queue/{}.log'.format(queue_name)
pid = os.getpid()

# Structure influenced by https://pika.readthedocs.io/en/stable/examples/direct_reply_to.html?direct-reply-to-example
# and https://www.rabbitmq.com/tutorials/tutorial-six-python.html

def _on_request(channel, method_frame, properties, body):
    """Callback to handle new message. Also replies to original publisher queue."""
    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)
    dataset_id = deserialized_body["dataset_id"]
    genecart_id = deserialized_body["genecart_id"]
    projection_id = deserialized_body["projection_id"]
    session_id = deserialized_body["session_id"]
    scope = deserialized_body["scope"]
    is_pca = deserialized_body["is_pca"]

    with open(stream, "a") as fh:
        print("{} - [x] - Received request for dataset {} and genecart {}".format(pid, dataset_id, genecart_id), flush=True, file=fh)
        output_payload = projectr_callback(dataset_id, genecart_id, projection_id, session_id, scope, is_pca, fh)

        # Send the output back to the Flask API call
        try:
            import pika
            channel.basic_publish(
                    exchange=""
                    , routing_key=properties.reply_to
                    , body=json.dumps(output_payload)
                    , properties=pika.BasicProperties(delivery_mode=2, content_type="application/json")
                    )
            print("{} - [x] - Publishing response for dataset {} and genecart {}".format(pid, dataset_id, genecart_id), flush=True, file=fh)
            channel.basic_ack(delivery_tag=delivery_tag)
        except Exception as e:
            print("{} - Could not deliver response back to client".format(pid), flush=True, file=fh)
            print("{} - {}".format(pid, str(e)), flush=True, file=fh)
        finally:
            gc.collect()

class Consumer:
    """This is a RabbitMQ consumer that will reconnect if the nested
    AsyncConnection indicates that a reconnect is necessary.
    """

    # Code taken from
    # https://github.com/pika/pika/blob/74b2c3be7ee4e366aa59a0ea05571d820dd145fc/examples/asynchronous_consumer_example.py

    def __init__(self, host):
        self._reconnect_delay = 0
        self.host = host

        # This differences from the example implementation
        # I connect the channel here, instead of the example's "connect" method.
        # I also decouple "start_consuming" out of the callback stack
        self._consumer = gearqueue.AsyncConnection(host=self.host, publisher_or_consumer="consumer", queue_name=queue_name, on_message_callback=_on_request, pid=pid, logfile=stream)

    def run(self):
        while True:
            try:
                # create the consumer.
                self._consumer.run()
                self._consumer.start_consuming()
            except KeyboardInterrupt:
                self._consumer.stop()
                break
            self._maybe_reconnect()

    def _maybe_reconnect(self):
        import time
        if self._consumer.should_reconnect:
            self._consumer.stop()
            reconnect_delay = self._get_reconnect_delay()
            time.sleep(reconnect_delay)
            self._consumer = gearqueue.AsyncConnection(host=self.host, publisher_or_consumer="consumer", queue_name=queue_name, on_message_callback=_on_request, pid=pid, logfile=stream)

    def _get_reconnect_delay(self):
        if self._consumer.was_consuming:
            self._reconnect_delay = 0
        else:
            self._reconnect_delay += 1
        if self._reconnect_delay > 30:
            self._reconnect_delay = 30
        return self._reconnect_delay

def main():
    #TODO: potentially add multiprocessing to start workers via script rather than running multiple times on command line
    host = this.servercfg['projectR_service']['queue_host']
    consumer = Consumer(host=host)
    consumer.run()

if __name__ == '__main__':
    main()

