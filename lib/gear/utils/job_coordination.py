# job_coordination.py - RabbitMQ consumer job-lock/retry/logging helpers.

import fcntl
import os
import typing
from datetime import datetime
from pathlib import Path


def log_line(fh: typing.TextIO, message: str) -> None:
    """
    Write a timestamped line to a RabbitMQ consumer log filehandle (e.g. the per-queue files
    under /var/log/gEAR_queue/). Plain print() to these files has no timestamp of its own -
    unlike stdout/stderr, which journald timestamps automatically for a systemd service - so
    without this, correlating events across workers/log lines (e.g. when a duplicate delivery
    landed relative to the original) requires cross-referencing against something else.
    """
    timestamp = datetime.now().isoformat(sep=" ", timespec="milliseconds")
    print(f"{timestamp} - {message}", flush=True, file=fh)


def check_and_record_attempt(
    state_dir: "str | Path", max_attempts: int = 2
) -> tuple[bool, int]:
    """
    Track how many times a job has been attempted, via a small counter file in the given
    directory (alongside e.g. a `.job.lock` file). Call this once per delivery -- after
    acquiring the job's lock, but before starting any real work -- so the count is recorded
    even if this attempt later gets killed (e.g. an OOM-kill), which is exactly the case this
    exists to bound: without it, a job whose memory need exceeds what's configured/available
    gets redelivered and retried forever, identically, every time.

    Because the caller is expected to hold the job's lock while calling this, the counter file
    itself needs no separate locking -- only one worker can be in that section for a given job
    at a time.

    Returns (allowed, attempt_number):
      - allowed=True: this attempt may proceed; attempt_number is 1-based (1 = first attempt,
        2 = first retry, etc).
      - allowed=False: max_attempts was already reached by a prior attempt; the caller should
        fail permanently without starting work. attempt_number is the attempt count already on
        record (i.e. how many attempts have already happened, not counting this one).
    """
    counter_file = Path(state_dir) / ".attempt_count"
    try:
        recorded_attempts = int(counter_file.read_text().strip()) if counter_file.exists() else 0
    except (ValueError, OSError):
        recorded_attempts = 0

    if recorded_attempts >= max_attempts:
        return False, recorded_attempts

    counter_file.write_text(str(recorded_attempts + 1))
    return True, recorded_attempts + 1


def clear_attempt_count(state_dir: "str | Path") -> None:
    """Remove the attempt counter -- call after a job's terminal outcome (success or a
    permanent give-up) is recorded, so a later unrelated re-upload of the same share_uid
    doesn't inherit a stale count."""
    counter_file = Path(state_dir) / ".attempt_count"
    try:
        counter_file.unlink()
    except FileNotFoundError:
        pass


def try_acquire_lock_file(filepath: "str | Path") -> typing.Optional[typing.TextIO]:
    """
    Attempt to acquire an exclusive, non-blocking lock at the given path.

    Returns the open file handle on success (caller is responsible for eventually passing it to
    release_lock_file()), or None if another live process already holds it. This never blocks
    waiting for the lock, so it's safe to call from a single-threaded event-loop callback (e.g. a
    pika/RabbitMQ consumer's on_message callback) without stalling the loop.

    The lock is released automatically by the kernel if the holding process dies for any reason
    (including an OOM-kill), so this can't be left permanently stuck the way a plain marker file
    (checked with a bare `Path.exists()`) could.
    """
    fd = open(filepath, "w+")
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        fd.close()
        return None
    fd.write(f"{os.getpid()}\n")
    fd.flush()
    return fd


def release_lock_file(fd: typing.TextIO, filepath: "str | Path") -> None:
    """Release a lock acquired via try_acquire_lock_file() and remove the lock file."""
    fd.close()
    try:
        Path(filepath).unlink()
    except FileNotFoundError:
        # This is fine, as the lock file may have been removed by another process
        pass
