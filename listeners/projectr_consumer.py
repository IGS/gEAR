#!/opt/bin/python3

"""
projectr_consumer.py - RabbitMQ messaging consumer
"""

import os, sys, json

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
    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)
    dataset_id = deserialized_body["dataset_id"]
    genecart_id = deserialized_body["genecart_id"]
    projection_id = deserialized_body["projection_id"]
    session_id = deserialized_body["session_id"]
    scope = deserialized_body["scope"]
    is_pca = deserialized_body["is_pca"]

    with open(stream, "a") as fh:
        print("{} - [x] - Received request for dataset {} and genecart {}".format(pid, dataset_id, genecart_id), file=fh)
        output_payload = projectr_callback(dataset_id, genecart_id, projection_id, session_id, scope, is_pca)

        # Send the output back to the Flask API call
        try:
            import pika
            channel.basic_publish(
                    exchange=""
                    , routing_key=properties.reply_to
                    , body=json.dumps(output_payload)
                    , properties=pika.BasicProperties(delivery_mode=2, content_type="application/json")
                    )
            print("{} - [x] - Publishing response for dataset {} and genecart {}".format(pid, dataset_id, genecart_id), file=fh)
            channel.basic_ack(delivery_tag=delivery_tag)
        except Exception as e:
            print("{} - Could not deliver response back to client".format(pid), file=fh)
            print("{} - {}".format(pid, str(e)), file=fh)

host = this.servercfg['projectR_service']['queue_host']

try:
    connection = gearqueue.Connection(host=host, publisher_or_consumer="consumer")
except:
    raise

with connection:
    # create the consumer.
    connection.consume(
        queue_name=queue_name
        , on_message_callback=_on_request
        , num_messages=1
    )
    # Consume indefinitely
    connection.continue_consuming(queue_name=queue_name)
