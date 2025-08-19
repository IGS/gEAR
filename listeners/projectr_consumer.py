#!/opt/bin/python3

"""
projectr_consumer.py - RabbitMQ messaging consumer
"""

import gc
import json
import os
import sys
from pathlib import Path

lib_path = str(Path(__file__).resolve().parents[1].joinpath("lib"))
sys.path.append(lib_path)

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig  # noqa: E402

this.servercfg = ServerConfig().parse() # type: ignore

queue_name = "projectr"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
stream = "/var/log/gEAR_queue/{}.log".format(queue_name)
pid = os.getpid()

# Structure influenced by https://pika.readthedocs.io/en/stable/examples/direct_reply_to.html?direct-reply-to-example
# and https://www.rabbitmq.com/tutorials/tutorial-six-python.html


def _on_request(channel, method_frame, properties, body):
    """Callback to handle new message. Also replies to original publisher queue."""

    www_path = str(Path(__file__).resolve().parents[1].joinpath("www"))
    sys.path.append(www_path)

    from api.resources.projectr import projectr_callback  # type: ignore

    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)
    dataset_id = deserialized_body["dataset_id"]
    genecart_id = deserialized_body["genecart_id"]
    projection_id = deserialized_body["projection_id"]
    session_id = deserialized_body["session_id"]
    scope = deserialized_body["scope"]
    algorithm = deserialized_body["algorithm"]
    zscore = deserialized_body["zscore"]
    full_output = deserialized_body["full_output"]

    with open(stream, "a") as fh:
        print(
            "{} - [x] - Received request for dataset {} and genecart {}".format(
                pid, dataset_id, genecart_id
            ),
            flush=True,
            file=fh,
        )

        try:
            # Run the callback function to generate the reply payload
            # We only need the output for a non-rabbitMQ implementation
            output_payload = projectr_callback(
                dataset_id,
                genecart_id,
                projection_id,
                session_id,
                scope,
                algorithm,
                zscore,
                full_output,
                fh,
            )
            channel.basic_ack(delivery_tag=delivery_tag)
        except Exception as e:
            print("{} - Caught error '{}'".format(pid, str(e)), flush=True, file=fh)
            channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
        finally:
            gc.collect()


class Consumer:
    """This is a RabbitMQ consumer that will reconnect if the nested
    AsyncConnection indicates that a reconnect is necessary.
    """

    # Code taken from
    # https://github.com/pika/pika/blob/74b2c3be7ee4e366aa59a0ea05571d820dd145fc/examples/asynchronous_consumer_example.py

    def __init__(self, host: str) -> None:
        self._reconnect_delay = 0
        self.host = host

        import gearqueue

        # This differences from the example implementation
        # I connect the channel here, instead of the example's "connect" method.
        # I also decouple "start_consuming" out of the callback stack
        self._consumer = gearqueue.AsyncConnection(
            host=self.host,
            publisher_or_consumer="consumer",
            queue_name=queue_name,
            on_message_callback=_on_request,
            pid=pid,
            logfile=stream,
            purge_queue=True,  # Purge the queue before starting
        )

    def run(self) -> None:
        while True:
            try:
                # create the consumer.
                self._consumer.run()
                self._consumer.start_consuming()
            except KeyboardInterrupt:
                self._consumer.stop()
                break
            self._maybe_reconnect()

    def _maybe_reconnect(self) -> None:
        import time

        import gearqueue

        if self._consumer.should_reconnect:
            self._consumer.stop()
            reconnect_delay = self._get_reconnect_delay()
            print("{} - Reconnecting after {} seconds".format(pid, reconnect_delay))
            time.sleep(reconnect_delay)
            self._consumer = gearqueue.AsyncConnection(
                host=self.host,
                publisher_or_consumer="consumer",
                queue_name=queue_name,
                on_message_callback=_on_request,
                pid=pid,
                logfile=stream,
            )

    def _get_reconnect_delay(self) -> int:
        if self._consumer.was_consuming:
            self._reconnect_delay = 0
        else:
            self._reconnect_delay += 1
        if self._reconnect_delay > 30:
            self._reconnect_delay = 30
        return self._reconnect_delay


def main() -> None:
    # TODO: potentially add multiprocessing to start workers via script rather than running multiple times on command line
    host = this.servercfg["projectR_service"]["queue_host"]
    consumer = Consumer(host=host)
    consumer.run()


if __name__ == "__main__":
    main()
