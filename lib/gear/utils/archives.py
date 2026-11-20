"""
archives.py - Friendly errors for uploaded archives (.tar, .tar.gz, .zip) that can't be extracted.

A truncated .tar.gz (e.g. one that wasn't fully created or copied) raises EOFError or zlib.error
from the gzip layer rather than tarfile.ReadError, so the upload processors catch all of these.
"""

import gzip
import tarfile
import zipfile
import zlib

# Errors raised while reading or extracting an archive that is incomplete or not a valid archive
ARCHIVE_READ_ERRORS = (tarfile.ReadError, zipfile.BadZipFile, EOFError, zlib.error, gzip.BadGzipFile)


class ArchiveReadError(Exception):
    """An uploaded archive could not be read; str() is a message for the user."""


def archive_error_message(error: BaseException, filename: str = "") -> str:
    """
    Return a message for the user explaining why the uploaded archive could not be extracted.

    Truncation (the archive ends before all of its files) gets its own message, since the fix is
    to re-create or re-copy the archive rather than to change its contents.
    """
    is_zip = str(filename).lower().endswith(".zip") or isinstance(error, zipfile.BadZipFile)
    check_command = "unzip -t <file>" if is_zip else "tar -tf <file>"
    truncated = isinstance(error, (EOFError, zlib.error)) or "unexpected end of data" in str(error)

    if truncated:
        return (
            "The uploaded archive appears to be incomplete: it ended before all of its files could "
            "be extracted. This usually means the archive was not fully created or copied before "
            f"uploading. Please re-create it, check that `{check_command}` lists every file without "
            "errors, and upload it again."
        )
    kind = "zip" if is_zip else ".tar or .tar.gz"
    return (
        f"The uploaded archive could not be read. It may be corrupted or not a valid {kind} file. "
        f"Please check that `{check_command}` lists its files without errors, and upload it again."
    )
