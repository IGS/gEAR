import json
import sys
from pathlib import Path
from uuid import uuid4

import pika
from flask import request
from flask_restful import Resource

gear_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(gear_root / 'lib'))

import geardb
from gear.trackhub import process_trackhub_synchronously, validate_hub_contents

user_upload_file_base = gear_root / 'www' / 'uploads' / 'files'


def queue_trackhub_job(job_id: str, share_uid: str, hub_json: dict, assembly: str, track_stanzas: list, dry_run: bool = False) -> bool:
    """Queue trackhub processing job to RabbitMQ."""
    try:
        import configparser
        config = configparser.ConfigParser()
        config.read(gear_root / 'gear.ini')

        # If queue is not enabled, return False
        if not config.getboolean('dataset_uploader', 'queue_enabled', fallback=False):
            print("Queue is disabled in configuration. Cannot queue trackhub job.", file=sys.stderr)
            return False

        rabbitmq_host = config.get('rabbitmq', 'host', fallback='localhost')
        rabbitmq_user = config.get('rabbitmq', 'user', fallback='guest')
        rabbitmq_password = config.get('rabbitmq', 'password', fallback='guest')

        credentials = pika.PlainCredentials(rabbitmq_user, rabbitmq_password)
        connection = pika.BlockingConnection(pika.ConnectionParameters(
            host=rabbitmq_host,
            credentials=credentials
        ))
        channel = connection.channel()

        queue_name = 'trackhub_copy_jobs'
        channel.queue_declare(queue=queue_name, durable=True)

        message = {
            'job_id': job_id,
            'share_uid': share_uid,
            'hub_json': hub_json,
            'assembly': assembly,
            'track_stanzas': track_stanzas,
            'dry_run': dry_run
        }

        channel.basic_publish(
            exchange='',
            routing_key=queue_name,
            body=json.dumps(message),
            properties=pika.BasicProperties(delivery_mode=2)
        )

        connection.close()
        return True
    except Exception as e:
        print(f"Error queuing trackhub job: {e}", file=sys.stderr)
        return False


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

        # Quick validation
        if not validate_hub_contents(hub_json, track_stanzas):
            result["message"] = "Hub contents failed validation"
            return result, 400

        # Generate job ID and create status file
        job_id = str(uuid4())
        staging_area = user_upload_file_base / session_id / share_uid
        status_file = staging_area / "trackhub_status.json"

        # Create initial status
        staging_area.mkdir(parents=True, exist_ok=True)
        initial_status = {
            "job_id": job_id,
            "status": "queued",
            "progress": 0,
            "completed_tracks": 0,
            "total_tracks": len(track_stanzas),
            "message": "Job queued for processing",
            "track_statuses": {}
        }
        with open(status_file, 'w') as f:
            json.dump(initial_status, f)

        # Queue the job
        if queue_trackhub_job(job_id, share_uid, hub_json, assembly, track_stanzas, dry_run):
            result["success"] = True
            result["job_id"] = job_id
            result["message"] = "Track hub processing job queued"
            return result, 202  # Accepted
        else:
            # Fall back to synchronous processing if queue is disabled
            import configparser
            config = configparser.ConfigParser()
            config.read(gear_root / 'gear.ini')

            higlass_config = None
            if config.has_section('higlass'):
                higlass_config = {
                    'higlass_hostname': config.get('higlass', 'hostname', fallback=''),
                    'higlass_admin_user': config.get('higlass', 'admin_user', fallback=''),
                    'higlass_admin_pass': config.get('higlass', 'admin_pass', fallback=''),
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

            result["success"] = result_sync["success"]
            result["job_id"] = job_id
            result["message"] = result_sync["message"]


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
