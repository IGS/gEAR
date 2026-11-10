"""
Stand-in for gear.anndata_processor, gear.spatial_processor and gearqueue in upload CGI tests.

process_uploaded_expression_dataset.cgi reads the repository's gear.ini. When the queue is disabled
(or there is no gear.ini) it falls back to the synchronous processors; when it is enabled it
publishes the job with gearqueue. Both stand-ins here do nothing and succeed, so tests can check
what the CGI wrote to metadata.json before handing off the job, whatever the local gear.ini says.
"""

import json


def write_status(status_file, status):
    with open(status_file, 'w') as f:
        f.write(json.dumps(status, indent=2))


def process_anndata_synchronously(**kwargs):
    return {"success": 1, "message": "processed (fake)"}


def process_spatial_synchronously(**kwargs):
    return {"success": 1, "message": "processed (fake)"}


class Connection:
    """gearqueue.Connection stand-in: publishing a job does nothing."""

    def __init__(self, host=None, publisher_or_consumer=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def open_channel(self):
        pass

    def publish(self, queue_name, message):
        pass
