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

from gear.anndata_processor import write_status  # noqa: E402
from gear.spatial_processor import process_spatial_synchronously  # noqa: E402
from gear.utils.job_coordination import (  # noqa: E402
    check_and_record_attempt,
    clear_attempt_count,
    log_line,
    release_lock_file,
    try_acquire_lock_file,
)

servercfg = ServerConfig().parse()

queue_name = "spatial_upload_jobs"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
logfile = f"/var/log/gEAR_queue/{queue_name}.log"
pid = os.getpid()

user_upload_base = gear_root / 'www' / 'uploads' / 'files'

# Original attempt + one retry. A job that still fails after this many tries is far more
# likely to genuinely need more memory than is available than to be a one-off transient
# hiccup, so it's better to fail cleanly than to let it redeliver and retry forever.
MAX_JOB_ATTEMPTS = 2


def _on_request(channel, method_frame, properties, body) -> None:
    """Callback to handle new spatial upload job message."""

    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)

    job_id = deserialized_body["job_id"]
    share_uid = deserialized_body["share_uid"]
    spatial_format = deserialized_body["spatial_format"]
    perform_primary_analysis = deserialized_body.get("perform_primary_analysis", False)

    with open(logfile, "a") as fh:
        log_line(fh, f"{pid} - [x] Received request for spatial job {job_id}")

        if not user_upload_base.is_dir():
            log_line(fh, f"{pid} - ERROR: User upload base directory {user_upload_base} does not exist")
            channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
            return

        lock_fh = None
        lockfile = None
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

            # Guard against a live duplicate: if another worker is still actively processing
            # this exact job (e.g. this delivery arrived because RabbitMQ's ack deadline was
            # exceeded while the original worker was still legitimately running -- spatial
            # processing can be slow enough for that), don't start a second concurrent run on
            # top of it -- ack-and-drop instead. Non-blocking, so it never stalls this ioloop
            # thread waiting on the lock.
            lockfile = staging_area / ".job.lock"
            lock_fh = try_acquire_lock_file(lockfile)
            if lock_fh is None:
                log_line(
                    fh,
                    f"{pid} - Job {job_id} for share {share_uid} is already being processed by "
                    "another worker (lock held); acking duplicate delivery without reprocessing.",
                )
                channel.basic_ack(delivery_tag=delivery_tag)
                return

            status_file = staging_area / "status.json"

            # Bound how many times this exact job gets retried. Recorded now, before any real
            # work starts, so the count is preserved even if this attempt itself gets killed
            # (e.g. an OOM-kill) -- otherwise a job whose memory need exceeds what's available
            # would get redelivered and retried identically, forever.
            allowed, attempt_number = check_and_record_attempt(
                staging_area, max_attempts=MAX_JOB_ATTEMPTS
            )
            if not allowed:
                message = (
                    f"This dataset failed to process after {MAX_JOB_ATTEMPTS} attempts and will "
                    "not be retried further. This usually means it needs more memory than is "
                    "available on this server, though it can also happen if several large jobs "
                    "happened to run at the same time. Please contact the gEAR team for help "
                    f"(share ID: {share_uid})."
                )
                log_line(
                    fh,
                    f"{pid} - Job {job_id} for share {share_uid} has already failed "
                    f"{attempt_number} time(s); giving up without retrying further.",
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
            if result.get("success"):
                clear_attempt_count(staging_area)

            log_line(fh, f"{pid} - Job {job_id}: {result['message']}")
            channel.basic_ack(delivery_tag=delivery_tag)
        except Exception as e:
            traceback.print_exc()
            log_line(fh, f"{pid} - Caught error '{str(e)}'")
            channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
        finally:
            if lock_fh is not None:
                release_lock_file(lock_fh, lockfile)
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
    from gear.utils.resource_limits import set_memory_limit_from_cgroup

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
