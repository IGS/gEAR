#!/opt/bin/python3

"""
gosling_upload_consumer.py - RabbitMQ consumer for track hub processing jobs.
"""

import gc
import json
import os
import sys
from pathlib import Path

lib_path = str(Path(__file__).resolve().parents[1] / 'lib')
sys.path.insert(0, lib_path)

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig  # noqa: E402

this.servercfg = ServerConfig().parse()  # type: ignore

queue_name = "trackhub_copy_jobs"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
logfile = f"/var/log/gEAR_queue/{queue_name}.log"
pid = os.getpid()

# Gear root for file operations
gear_root = Path(__file__).resolve().parents[1]
user_upload_base = gear_root / 'www' / 'uploads' / 'files'


def _on_request(channel, method_frame, properties, body):
    """Callback to handle new trackhub job message."""
    from gear.trackhub import TrackHubProcessor  # noqa: E402

    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)

    job_id = deserialized_body["job_id"]
    share_uid = deserialized_body["share_uid"]
    hub_json = deserialized_body["hub_json"]
    assembly = deserialized_body["assembly"]
    track_stanzas = deserialized_body["track_stanzas"]
    dry_run = deserialized_body.get("dry_run", False)

    with open(logfile, "a") as fh:
        print(
            f"{pid} - [x] - Received request for trackhub job {job_id}",
            flush=True,
            file=fh,
        )

        try:
            # Infer session_id from share_uid directory structure
            staging_area = None
            for session_dir in user_upload_base.iterdir():
                candidate = session_dir / share_uid
                if candidate.is_dir():
                    staging_area = candidate
                    break

            if not staging_area:
                raise FileNotFoundError(f"Could not find staging area for {share_uid}")

            status_file = staging_area / "trackhub_status.json"

            # Get HiGlass config
            higlass_config = None
            if this.servercfg.get("higlass"):
                higlass_config = {
                    "higlass_hostname": this.servercfg["higlass"].get("hostname", ""),
                    "higlass_admin_user": this.servercfg["higlass"].get("admin_user", ""),
                    "higlass_admin_pass": this.servercfg["higlass"].get("admin_pass", ""),
                }

            # Process the job
            processor = TrackHubProcessor(
                job_id=job_id,
                share_uid=share_uid,
                staging_area=staging_area,
                status_file=status_file,
                higlass_config=higlass_config,
            )

            result = processor.process(hub_json, assembly, track_stanzas, dry_run)
            print(f"{pid} - Job {job_id}: {result['message']}", flush=True, file=fh)
            channel.basic_ack(delivery_tag=delivery_tag)

        except Exception as e:
            print(f"{pid} - Caught error '{str(e)}'", flush=True, file=fh)
            channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
        finally:
            gc.collect()


class Consumer:
    """RabbitMQ consumer with automatic reconnection for trackhub processing."""

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
        while True:
            try:
                self._consumer.run()
            except KeyboardInterrupt:
                self._consumer.stop()
                break
            self._maybe_reconnect()

    def _maybe_reconnect(self) -> None:
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
        if self._consumer.was_consuming:
            self._reconnect_delay = 0
        else:
            self._reconnect_delay += 1
        if self._reconnect_delay > 30:
            self._reconnect_delay = 30
        return self._reconnect_delay


def main() -> None:
    """Start the trackhub processing consumer."""
    host = this.servercfg["rabbitmq"]["host"]
    consumer = Consumer(host=host)
    consumer.run()


if __name__ == "__main__":
    main()
