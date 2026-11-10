"""
Stand-in for gear.anndata_processor and gear.spatial_processor in upload CGI tests.

process_uploaded_expression_dataset.cgi falls back to these synchronous processors when the queue
is disabled (no gear.ini in tests). They do no processing and report success, so tests can check
what the CGI wrote to metadata.json before handing off the job.
"""

import json


def write_status(status_file, status):
    with open(status_file, 'w') as f:
        f.write(json.dumps(status, indent=2))


def process_anndata_synchronously(**kwargs):
    return {"success": 1, "message": "processed (fake)"}


def process_spatial_synchronously(**kwargs):
    return {"success": 1, "message": "processed (fake)"}
