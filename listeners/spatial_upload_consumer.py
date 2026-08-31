#!/opt/bin/python3

"""
spatial_upload_consumer.py - RabbitMQ consumer for spatial dataset upload jobs.

Processes Visium, VisiumHD, Curio, GeoMx, CosMx, and Xenium format datasets
into a SpatialData object and writes the result as a Zarr store.

Handled separately from anndata_upload_consumer.py since spatial uploads
produce a Zarr store rather than an H5AD file, and depend on the
spatialdata/spatialdata_io stack rather than plain anndata/pandas.
"""

import gc
import json
import os
import sys
import time
import traceback
from pathlib import Path

# Gear root for file operations
gear_root = Path(__file__).resolve().parents[1]
gear_lib = gear_root / "lib"
sys.path.insert(0, str(gear_lib))

import gearqueue
from gear.serverconfig import ServerConfig  # noqa: I001

servercfg = ServerConfig().parse()

queue_name = "spatial_upload_jobs"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
logfile = f"/var/log/gEAR_queue/{queue_name}.log"
pid = os.getpid()

user_upload_base = gear_root / 'www' / 'uploads' / 'files'


def _on_request(channel, method_frame, properties, body) -> None:
    """Callback to handle new spatial upload job message."""
    from gear.anndata_processor import write_status  # noqa: E402
    from gear.spatial_processor import process_spatial_synchronously  # noqa: E402

    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)

    job_id = deserialized_body["job_id"]
    share_uid = deserialized_body["share_uid"]
    spatial_format = deserialized_body["spatial_format"]
    perform_primary_analysis = deserialized_body.get("perform_primary_analysis", False)

    with open(logfile, "a") as fh:
        print(
            f"{pid} - [x] Received request for spatial job {job_id}",
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

            # A redelivered message means this exact job was already delivered once and
            # never acked/nacked - almost always because the previous worker died mid-job
            # (e.g. a hard OOM-kill). Reprocessing would just repeat the same crash forever,
            # so fail it cleanly here instead of retrying.
            if method_frame.redelivered:
                message = (
                    "A previous attempt to process this dataset was interrupted unexpectedly "
                    "(the worker did not finish cleanly, most likely because it ran out of "
                    "available memory). Please contact the gEAR team to resolve this issue "
                    f"(share ID: {share_uid})."
                )
                print(
                    f"{pid} - Job {job_id} was redelivered; failing without retry: {message}",
                    flush=True,
                    file=fh,
                )
                write_status(
                    status_file,
                    {"job_id": job_id, "status": "error", "message": message, "progress": 0},
                )
                channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
                return

            result = process_spatial_synchronously(
                job_id=job_id,
                share_uid=share_uid,
                staging_area=staging_area,
                status_file=status_file,
                spatial_format=spatial_format,
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
    """RabbitMQ consumer with automatic reconnection for spatial uploads."""

    def __init__(self, host: str) -> None:
        self._reconnect_delay = 0
        self.host = host

        self._consumer = self._new_connection()

    def run(self) -> None:
        """Run the consumer with automatic reconnection."""
        while True:
            try:
                self._consumer.run()
            except KeyboardInterrupt:
                self._consumer.stop()
                break
            except Exception as exc:
                print(f"{pid} - Consumer loop error: {exc}", flush=True)
                traceback.print_exc()
                self._consumer.should_reconnect = True

            if not getattr(self._consumer, "should_reconnect", False):
                break

            self._maybe_reconnect()

    def _new_connection(self) -> "gearqueue.AsyncConnection":
        """Create a new AsyncConnection instance."""
        return gearqueue.AsyncConnection(
            host=self.host,
            publisher_or_consumer="consumer",
            queue_name=queue_name,
            on_message_callback=_on_request,
            pid=pid,
            logfile=logfile,
            purge_queue=False,
        )

    def _maybe_reconnect(self) -> None:
        """Attempt reconnection with exponential backoff."""
        if self._consumer.should_reconnect:
            self._consumer.stop()
            reconnect_delay = self._get_reconnect_delay()
            print(
                f"{pid} - Reconnecting after {reconnect_delay} seconds",
                flush=True,
            )
            time.sleep(reconnect_delay)
            self._consumer = self._new_connection()

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
    """Start the spatial processing consumer."""
    from gear.utils import set_memory_limit_from_cgroup

    # Spatial processing (spatialdata_io readers) can spike memory well above what a
    # clean Python exception would normally warn about. Self-impose a ceiling below the
    # container's actual cgroup limit so approaching it raises a catchable MemoryError
    # instead of the kernel OOM-killer sending an uncatchable SIGKILL.
    set_memory_limit_from_cgroup()

    host = servercfg["dataset_uploader"]["queue_host"]
    consumer = Consumer(host=host)
    consumer.run()


if __name__ == "__main__":
    main()
