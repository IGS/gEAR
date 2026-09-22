#!/opt/bin/python3

"""
anndata_upload_consumer.py - RabbitMQ consumer for expression dataset upload jobs.

Processes H5AD, Seurat, 3-tab, Excel, and MEX format datasets.

Spatial uploads are handled differently as they create a SpatialData object
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

from gear.anndata_processor import AnndataProcessor  # noqa: E402
from gear.utils import log_line, release_lock_file, try_acquire_lock_file  # noqa: E402

servercfg = ServerConfig().parse()

queue_name = "anndata_upload_jobs"
os.makedirs("/var/log/gEAR_queue", exist_ok=True)
logfile = f"/var/log/gEAR_queue/{queue_name}.log"
pid = os.getpid()

user_upload_base = gear_root / 'www' / 'uploads' / 'files'


def _on_request(channel, method_frame, properties, body) -> None:
    """Callback to handle new anndata upload job message."""

    delivery_tag = method_frame.delivery_tag
    deserialized_body = json.loads(body)

    job_id = deserialized_body["job_id"]
    share_uid = deserialized_body["share_uid"]
    dataset_uid = deserialized_body["dataset_uid"]
    dataset_format = deserialized_body["dataset_format"]
    perform_primary_analysis = deserialized_body.get("perform_primary_analysis", False)

    with open(logfile, "a") as fh:
        log_line(fh, f"{pid} - [x] Received request for anndata job {job_id}")

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

            # Seurat/anndata conversion can legitimately run long enough to exceed RabbitMQ's
            # ack deadline, which causes the broker to redeliver this same job to another worker
            # while this one is still processing it. Guard against that duplicate run: if another
            # live process already holds the lock, just ack-and-drop this delivery rather than
            # starting a second full conversion on top of it. Non-blocking, so it never stalls
            # this ioloop thread waiting on the lock.
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

            log_line(fh, f"{pid} - Job {job_id}: {result['message']}")
            if channel.is_open:
                channel.basic_ack(delivery_tag=delivery_tag)
            else:
                # Broker likely closed the channel (e.g. ack deadline exceeded) and already
                # redelivered this message elsewhere. Acking here would raise and escape this
                # except block unhandled, so just log and move on.
                log_line(fh, f"{pid} - Channel already closed, could not ack delivery {delivery_tag}")
        except Exception as e:
            traceback.print_exc()
            log_line(fh, f"{pid} - Caught error '{str(e)}'")
            try:
                if channel.is_open:
                    channel.basic_nack(delivery_tag=delivery_tag, requeue=False)
            except Exception as nack_err:
                log_line(fh, f"{pid} - Could not nack delivery {delivery_tag}: '{str(nack_err)}'")
        finally:
            if lock_fh is not None:
                release_lock_file(lock_fh, lockfile)
            gc.collect()


class Consumer:
    """RabbitMQ consumer with automatic reconnection for anndata uploads."""

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
    """Start the anndata processing consumer."""

    from gear.utils import set_memory_limit_from_cgroup

    # Sometimes processing can spike memory well above what a
    # clean Python exception would normally warn about. Self-impose a ceiling below the
    # container's actual cgroup limit so approaching it raises a catchable MemoryError
    # instead of the kernel OOM-killer sending an uncatchable SIGKILL.
    set_memory_limit_from_cgroup()


    host = servercfg["dataset_uploader"]["queue_host"]
    consumer = Consumer(host=host)
    consumer.run()


if __name__ == "__main__":
    main()