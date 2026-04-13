#!/opt/bin/python3

"""
anndata_upload_consumer.py - RabbitMQ consumer for expression dataset upload jobs.

Processes H5AD, 3-tab, Excel, and MEX format datasets.
"""

import gc
import json
import os
import traceback
from pathlib import Path

from gear.serverconfig import ServerConfig  # noqa: I001

servercfg = ServerConfig().parse()

queue_name = "anndata_upload_jobs"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
logfile = f"/var/log/gEAR_queue/{queue_name}.log"
pid = os.getpid()

gear_root = Path(__file__).resolve().parents[1]
user_upload_base = gear_root / 'www' / 'uploads' / 'files'


def _on_request(channel, method_frame, properties, body) -> None:
    """Callback to handle new anndata upload job message."""
    from gear.anndata_processor import AnndataProcessor  # noqa: E402

    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)

    job_id = deserialized_body["job_id"]
    share_uid = deserialized_body["share_uid"]
    dataset_uid = deserialized_body["dataset_uid"]
    dataset_format = deserialized_body["dataset_format"]
    perform_primary_analysis = deserialized_body.get("perform_primary_analysis", False)

    with open(logfile, "a") as fh:
        print(
            f"{pid} - [x] Received request for anndata job {job_id}",
            flush=True,
            file=fh,
        )

        if not user_upload_base.is_dir():
            print(
                f"{pid} - ERROR: User upload base directory {user_upload_base} does not exist",
                flush=True,
                file=fh,
            )
            channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
            return

        try:
            # Infer staging_area from share_uid directory structure
            staging_area = None
            for session_dir in user_upload_base.iterdir():
                candidate = session_dir / share_uid
                if candidate.is_dir():
                    staging_area = candidate
                    break

            if not staging_area:
                raise FileNotFoundError(f"Could not find staging area for {share_uid}")

            status_file = staging_area / "status.json"

            # Process the job
            processor = AnndataProcessor(
                job_id=job_id,
                share_uid=share_uid,
                staging_area=staging_area,
                status_file=status_file,
                dataset_uid=dataset_uid,
            )

            result = processor.process(
                dataset_format=dataset_format,
                perform_primary_analysis=perform_primary_analysis,
            )

            print(f"{pid} - Job {job_id}: {result['message']}", flush=True, file=fh)
            channel.basic_ack(delivery_tag=delivery_tag)

        except Exception as e:
            traceback.print_exc()
            print(f"{pid} - Caught error '{str(e)}'", flush=True, file=fh)
            channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
        finally:
            gc.collect()


class Consumer:
    """RabbitMQ consumer with automatic reconnection for anndata uploads."""

    def __init__(self, host: str) -> None:
        self._reconnect_delay = 0
        self.host = host

        import gearqueue  # noqa: F401

        self._consumer = gearqueue.AsyncConnection(
            host=self.host,
            publisher_or_consumer="consumer",
            queue_name=queue_name,
            on_message_callback=_on_request,
            pid=pid,
            logfile=logfile,
            purge_queue=False,
        )

    def run(self) -> None:
        """Run the consumer with automatic reconnection."""
        while True:
            try:
                self._consumer.run()
            except KeyboardInterrupt:
                self._consumer.stop()
                break
            self._maybe_reconnect()

    def _maybe_reconnect(self) -> None:
        """Attempt reconnection with exponential backoff."""
        import time

        import gearqueue  # noqa: F401

        if self._consumer.should_reconnect:
            self._consumer.stop()
            reconnect_delay = self._get_reconnect_delay()
            print(
                f"{pid} - Reconnecting after {reconnect_delay} seconds",
                flush=True,
            )
            time.sleep(reconnect_delay)
            self._consumer = gearqueue.AsyncConnection(
                host=self.host,
                publisher_or_consumer="consumer",
                queue_name=queue_name,
                on_message_callback=_on_request,
                pid=pid,
                logfile=logfile,
            )

    def _get_reconnect_delay(self) -> int:
        """Calculate reconnect delay with exponential backoff."""
        if self._consumer.was_consuming:
            self._reconnect_delay = 0
        else:
            self._reconnect_delay += 1
        if self._reconnect_delay > 30:
            self._reconnect_delay = 30
        return self._reconnect_delay


def main() -> None:
    """Start the anndata processing consumer."""
    host = servercfg["dataset_uploader"]["queue_host"]
    consumer = Consumer(host=host)
    consumer.run()


if __name__ == "__main__":
    main()