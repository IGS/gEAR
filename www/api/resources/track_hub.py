import configparser
import json
import sys
from pathlib import Path
from uuid import uuid4

from flask import request
from flask_restful import Resource

gear_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(gear_root / 'lib'))

_config = configparser.ConfigParser()
_config.read(gear_root / 'gear.ini')

import geardb
from gear.trackhub import process_trackhub_synchronously

user_upload_file_base = gear_root / 'www' / 'uploads' / 'files'

class QueueDisabledError(Exception):
    """Custom exception to indicate that the queue is disabled in configuration."""
    pass

def _create_initial_status_file(
    status_file: Path,
    job_id: str,
    total_tracks: int,
    message: str = "Job queued for processing",
) -> None:
    """Create initial status.json file."""
    initial_status = {
        "job_id": job_id,
        "status": "queued",
        "progress": 0,
        "completed_tracks": 0,
        "total_tracks": total_tracks,
        "message": message,
        "track_statuses": {},
    }
    with open(status_file, 'w') as f:
        json.dump(initial_status, f, indent=4)

def queue_trackhub_job(job_id: str, share_uid: str, hub_json: dict, assembly: str, track_stanzas: list, dry_run: bool = False) -> None:
    """Queue trackhub processing job to RabbitMQ."""

    # If queue is not enabled, return False
    if not _config.getboolean('dataset_uploader', 'queue_enabled', fallback=False):
        print("Queue is disabled in configuration. Cannot queue trackhub job. Falling back to synchronous processing.", file=sys.stderr)
        raise QueueDisabledError()

    import gearqueue
    host = _config["dataset_uploader"]["queue_host"]

    try:
        # Connect as a blocking RabbitMQ publisher
        connection = gearqueue.Connection(
            host=host, publisher_or_consumer="publisher"
        )
    except Exception as e:
        print(f"Error connecting to RabbitMQ: {e}", file=sys.stderr)
        raise Exception(f"Error connecting to RabbitMQ: {e}")

    with connection:
        connection.open_channel()

        payload = {
            'job_id': job_id,
            'share_uid': share_uid,
            'hub_json': hub_json,
            'assembly': assembly,
            'track_stanzas': track_stanzas,
            'dry_run': dry_run
        }

        try:
            connection.publish(
                queue_name="trackhub_copy_jobs",
                message=payload,  # method dumps JSON
            )
        except Exception as e:
            print(f"Error publishing message to RabbitMQ: {e}", file=sys.stderr)
            raise
    return



class TrackHubCopy(Resource):
    def post(self, share_uid):
        req = request.get_json()
        if req is None:
            return {"success": False, "message": "Invalid JSON body"}, 400

        session_id = request.cookies.get('gear_session_id', "")
        hub_json = req.get("hub_json")
        assembly = req.get("assembly")
        track_stanzas = req.get("tracks")
        dry_run = req.get("dry_run", False)

        result = {"success": False, "message": "", "job_id": None}

        user = geardb.get_user_from_session_id(session_id)
        if not user:
            result["message"] = "Invalid session. Please log in."
            return result, 401

        if not hub_json or not assembly or not track_stanzas:
            result["message"] = "Missing required parameters"
            return result, 400

        # Generate job ID and create status file
        job_id = str(uuid4())
        staging_area = user_upload_file_base / session_id / share_uid
        status_file = staging_area / "trackhub_status.json"

        # Create initial status
        staging_area.mkdir(parents=True, exist_ok=True)
        _create_initial_status_file(status_file, job_id, len(track_stanzas))

        result["job_id"] = job_id
        # Queue the job
        try:
            queue_trackhub_job(job_id, share_uid, hub_json, assembly, track_stanzas, dry_run)
            result["success"] = True
            result["message"] = "Track hub processing job queued"
            return result, 202  # Accepted
        except QueueDisabledError:
            higlass_config = None
            if _config.has_section('higlass'):
                higlass_config = {
                    'higlass_hostname': _config.get('higlass', 'hostname', fallback=''),
                    'higlass_admin_user': _config.get('higlass', 'admin_user', fallback=''),
                    'higlass_admin_pass': _config.get('higlass', 'admin_pass', fallback=''),
                }

            result_sync = process_trackhub_synchronously(
                job_id=job_id,
                share_uid=share_uid,
                staging_area=staging_area,
                status_file=status_file,
                hub_json=hub_json,
                assembly=assembly,
                track_stanzas=track_stanzas,
                dry_run=dry_run,
                higlass_config=higlass_config,
            )
            #

            result["success"] = result_sync["success"]
            result["message"] = result_sync["message"]
            return result, 200 if result["success"] else 500

        except Exception as e:
            result["message"] = f"Error processing track hub: {str(e)}"
            print(f"TrackHubCopy error: {str(e)}", file=sys.stderr)
            return result, 500


class TrackHubStatus(Resource):
    def post(self, share_uid):
        session_id = request.cookies.get('gear_session_id', "")
        req = request.get_json()
        if req is None:
            return {"success": False, "message": "Invalid JSON body"}, 400

        job_id = req.get("job_id")
        if not job_id:
            return {"success": False, "message": "Missing job_id"}, 400

        staging_area = user_upload_file_base / session_id / share_uid
        status_file = staging_area / "trackhub_status.json"

        if not status_file.is_file():
            return {"success": False, "message": "Job status not found"}, 404

        try:
            with open(status_file, 'r') as f:
                status_data = json.load(f)
            return status_data, 200
        except Exception as e:
            return {"success": False, "message": str(e)}, 500
