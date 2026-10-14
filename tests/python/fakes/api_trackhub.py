"""
api_trackhub.py - Stand-in for lib/gear/trackhub.py, used only when the real module can't be
imported (it needs pyBigWig and Biopython, which the test requirements don't install).

write_status matches the real function; process_trackhub_synchronously must be monkeypatched.
"""

import json


def write_status(status_file, job_id, status, progress, completed_tracks, total_tracks,
                 message="", track_statuses=None):
    status_data = {
        "job_id": job_id,
        "status": status,
        "progress": progress,
        "completed_tracks": completed_tracks,
        "total_tracks": total_tracks,
        "message": message,
        "track_statuses": track_statuses or {},
    }
    with open(status_file, "w") as f:
        json.dump(status_data, f, indent=4)


def process_trackhub_synchronously(*args, **kwargs):
    raise AssertionError("gear.trackhub stub called; monkeypatch process_trackhub_synchronously in the test")
