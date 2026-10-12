"""Upload-worker coordination helpers: job locks, attempt counting, cancellation, archive errors and MemoryError handling."""

import io
import os
import re
import tarfile
import zipfile
import zlib

import pytest
from gear.utils import job_coordination as jc
from gear.utils.archives import ARCHIVE_READ_ERRORS, archive_error_message
from gear.utils.resource_limits import catch_memory_error

# ---------------------------------------------------------------- cancellation

def test_existing_staging_dir_is_not_cancelled(tmp_path):
    jc.raise_if_upload_deleted(tmp_path)
    jc.raise_if_upload_deleted(str(tmp_path))


def test_deleted_staging_dir_cancels(tmp_path):
    staging = tmp_path / "upload"
    with pytest.raises(jc.UploadCancelledError, match=re.escape(str(staging))):
        jc.raise_if_upload_deleted(staging)


def test_staging_path_that_is_a_file_cancels(tmp_path):
    staging = tmp_path / "upload"
    staging.write_text("")
    with pytest.raises(jc.UploadCancelledError):
        jc.raise_if_upload_deleted(staging)


# ---------------------------------------------------------------- attempts

def test_attempts_are_counted_up_to_the_limit(tmp_path):
    assert jc.check_and_record_attempt(tmp_path) == (True, 1)
    assert jc.check_and_record_attempt(tmp_path) == (True, 2)
    assert jc.check_and_record_attempt(tmp_path) == (False, 2)
    assert jc.check_and_record_attempt(tmp_path) == (False, 2)
    assert (tmp_path / ".attempt_count").read_text() == "2"


def test_max_attempts(tmp_path):
    assert jc.check_and_record_attempt(tmp_path, max_attempts=1) == (True, 1)
    assert jc.check_and_record_attempt(tmp_path, max_attempts=1) == (False, 1)
    assert jc.check_and_record_attempt(tmp_path, max_attempts=3) == (True, 2)


def test_unreadable_counter_starts_over(tmp_path):
    (tmp_path / ".attempt_count").write_text("garbage")
    assert jc.check_and_record_attempt(str(tmp_path)) == (True, 1)


def test_clear_attempt_count(tmp_path):
    jc.check_and_record_attempt(tmp_path)
    jc.clear_attempt_count(tmp_path)
    assert not (tmp_path / ".attempt_count").exists()
    jc.clear_attempt_count(tmp_path)
    assert jc.check_and_record_attempt(tmp_path) == (True, 1)


# ---------------------------------------------------------------- locks

def test_lock_is_exclusive_until_released(tmp_path):
    lockfile = tmp_path / ".job.lock"
    fd = jc.try_acquire_lock_file(lockfile)
    assert fd is not None
    assert lockfile.read_text() == f"{os.getpid()}\n"
    try:
        assert jc.try_acquire_lock_file(lockfile) is None
    finally:
        jc.release_lock_file(fd, lockfile)
    assert not lockfile.exists()

    fd = jc.try_acquire_lock_file(str(lockfile))
    assert fd is not None
    jc.release_lock_file(fd, str(lockfile))


def test_failed_lock_attempt_keeps_holder_pid(tmp_path):
    lockfile = tmp_path / ".job.lock"
    fd = jc.try_acquire_lock_file(lockfile)
    try:
        assert jc.try_acquire_lock_file(lockfile) is None
        assert lockfile.read_text() == f"{os.getpid()}\n"
    finally:
        jc.release_lock_file(fd, lockfile)


def test_stale_lock_file_is_overwritten(tmp_path):
    lockfile = tmp_path / ".job.lock"
    lockfile.write_text("99999\nleftover from a killed worker\n")
    fd = jc.try_acquire_lock_file(lockfile)
    try:
        assert lockfile.read_text() == f"{os.getpid()}\n"
    finally:
        jc.release_lock_file(fd, lockfile)


def test_lock_on_a_replaced_file_is_refused(tmp_path, monkeypatch):
    """If the holder releases (removing the file) and another worker recreates it between our open()
    and flock(), our lock is on the old file, so taking it would let two workers run at once."""
    lockfile = tmp_path / ".job.lock"
    real_flock = jc.fcntl.flock

    def racing_flock(fd, op):
        lockfile.unlink()
        lockfile.write_text("")
        real_flock(fd, op)

    monkeypatch.setattr(jc.fcntl, "flock", racing_flock)
    assert jc.try_acquire_lock_file(lockfile) is None


def test_release_tolerates_missing_lock_file(tmp_path):
    lockfile = tmp_path / ".job.lock"
    fd = jc.try_acquire_lock_file(lockfile)
    lockfile.unlink()
    jc.release_lock_file(fd, lockfile)
    assert fd.closed


def test_log_line_is_timestamped():
    fh = io.StringIO()
    jc.log_line(fh, "started job")
    assert re.fullmatch(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d\.\d{3} - started job\n", fh.getvalue())


# ---------------------------------------------------------------- archive errors

@pytest.mark.parametrize("error, filename", [
    (EOFError("Compressed file ended before the end-of-stream marker was reached"), "data.tar.gz"),
    (zlib.error("Error -3 while decompressing data"), "data.tar.gz"),
    (tarfile.ReadError("unexpected end of data"), "data.tar"),
])
def test_truncated_archive_message(error, filename):
    message = archive_error_message(error, filename)
    assert "appears to be incomplete" in message and "tar -tf <file>" in message


@pytest.mark.parametrize("error, filename", [
    (zipfile.BadZipFile("File is not a zip file"), ""),
    (tarfile.ReadError("bad"), "DATA.ZIP"),
])
def test_bad_zip_message(error, filename):
    message = archive_error_message(error, filename)
    assert "not a valid zip file" in message and "unzip -t <file>" in message


def test_bad_tar_message():
    message = archive_error_message(tarfile.ReadError("file could not be opened successfully"), "data.tar.gz")
    assert "not a valid .tar or .tar.gz file" in message and "tar -tf <file>" in message


def test_real_truncated_tar_gz_gets_incomplete_message(tmp_path):
    payload = tmp_path / "matrix.tab"
    payload.write_bytes(os.urandom(64 * 1024))
    archive = tmp_path / "upload.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(payload, arcname="matrix.tab")
    archive.write_bytes(archive.read_bytes()[: archive.stat().st_size // 2])

    with pytest.raises(ARCHIVE_READ_ERRORS) as excinfo:
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(tmp_path / "out", filter="data")
    assert "appears to be incomplete" in archive_error_message(excinfo.value, archive.name)


# ---------------------------------------------------------------- MemoryError

def test_catch_memory_error_returns_500():
    @catch_memory_error()
    def load_dataset():
        raise MemoryError("too big")

    assert load_dataset.__name__ == "load_dataset"
    assert load_dataset() == ({"message": "Exceeded memory limit", "success": -1, "error": "too big"}, 500)


def test_catch_memory_error_passes_through_results_and_other_errors():
    @catch_memory_error()
    def double(value):
        if value < 0:
            raise ValueError("negative")
        return value * 2

    assert double(4) == 8
    with pytest.raises(ValueError):
        double(-1)
