"""
track_hub.py - Import uploaded UCSC track hub tracks as a Gosling dataset.

Serves /import/trackhub/<share_uid>/copy in www/api/api.py. Jobs are queued to
RabbitMQ when enabled, otherwise processed synchronously.
"""

import configparser
import json
import os
import sys
from pathlib import Path
from uuid import uuid4

from flask import request
from flask_restful import Resource
from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(gear_root / 'lib'))

_config = configparser.ConfigParser()
_config.read(gear_root / 'gear.ini')

import geardb
from gear.trackhub import process_trackhub_synchronously, write_status

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
    write_status(
        status_file,
        job_id=job_id,
        status="queued",
        message=message,
        progress=0,
        completed_tracks=0,
        total_tracks=total_tracks,
        track_statuses={},
    )

def queue_trackhub_job(job_id: str, share_uid: str, hub_json: dict, assembly: str, track_stanzas: list, hub_url: str = "", dry_run: bool = False) -> None:
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
            'hub_url': hub_url,
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
    """
    Flask-RESTful resource that stages and processes an uploaded track hub.
    """
    def post(self, share_uid):
        """
        Save uploaded track files and queue (or run) track hub processing.

        Form keys: hub_json, assembly, tracks (JSON list of track stanzas), dry_run,
        plus "<track_id>[file]" uploads.

        Returns:
            tuple: (dict with "success", "message", and "job_id", HTTP status).
        """
        req_form = request.form
        if req_form is None:
            return {"success": False, "message": "Invalid JSON body"}, 400

        session_id = request.cookies.get('gear_session_id', "")

        # share_uid and session_id become directory names below, so reject anything that isn't a plain name
        if not share_uid or secure_filename(share_uid) != share_uid:
            return {"success": False, "message": "Invalid share ID"}, 400

        hub_json = req_form.get("hub_json")
        if not hub_json:
            return {"success": False, "message": "Missing 'hub_json' parameter"}, 400
        assembly = req_form.get("assembly")
        tracks = req_form.get("tracks")
        if not tracks:
            return {"success": False, "message": "Missing 'tracks' parameter"}, 400
        try:
            hub_json = json.loads(hub_json)
            track_stanzas: list = json.loads(tracks)
        except json.JSONDecodeError as e:
            return {"success": False, "message": f"Invalid JSON in 'hub_json' or 'tracks': {e}"}, 400
        dry_run = req_form.get("dry_run", False)
        # convert dry_run to boolean if it's a string
        if isinstance(dry_run, str):
            dry_run = dry_run.lower() == "true"

        result = {"success": False, "message": "", "job_id": None}

        user = geardb.get_user_from_session_id(session_id)
        if not user or secure_filename(session_id) != session_id:
            result["message"] = "Invalid session. Please log in."
            return result, 401

        if not hub_json or not assembly or not track_stanzas:
            result["message"] = "Missing required parameters"
            return result, 400

        # Generate job ID and create status file
        job_id = str(uuid4())
        staging_area = user_upload_file_base / session_id / share_uid
        status_file = staging_area / "status.json"
        result["job_id"] = job_id

        def fail(message, http_status):
            """Record the error in the status file (if the upload exists) and return it."""
            if staging_area.is_dir():
                write_status(
                    status_file,
                    job_id=job_id,
                    status="error",
                    message=message,
                    progress=0,
                    completed_tracks=0,
                    total_tracks=len(track_stanzas),
                    track_statuses={},
                )
            result["message"] = message
            return result, http_status

        # Check the upload's metadata, the uploaded file names and the configuration before
        #  saving anything, so a request that can't be processed leaves no files behind
        metadata_file = staging_area / 'metadata.json'
        if not metadata_file.is_file():
            return fail("Metadata file not found. Impossible to save as dataset.", 400)
        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
        except json.JSONDecodeError:
            return fail("Metadata file could not be read. Impossible to save as dataset.", 400)

        dataset_id = metadata.get("dataset_uid", "") if isinstance(metadata, dict) else ""
        if not dataset_id:
            return fail("Dataset ID not found in metadata. Impossible to save as dataset.", 400)

        uploads = []  # (track_id, file, safe_filename)
        for track_key in request.files:
            if '[file]' in track_key:
                file = request.files.get(track_key)
                if file and file.filename:
                    track_id = track_key.split('[')[1].split(']')[0]
                    # The browser-supplied name could contain "../"; keep only a safe base name
                    safe_filename = secure_filename(file.filename)
                    if not safe_filename:
                        return fail(f"Invalid file name for track '{track_id}': {file.filename!r}", 400)
                    uploads.append((track_id, file, safe_filename))

        if os.getenv("ENVIRONMENT", "production").lower() == "development":
            domain_url = "http://localhost:8080"
        else:
            domain_url = geardb._read_domain_url()
        if not domain_url:
            return fail("Domain URL not configured. Cannot process track hub.", 500)

        _create_initial_status_file(status_file, job_id, len(track_stanzas))

        # Cannot serialize File object into JSON for RabbitMQ so save immediately
        write_status(
            status_file,
            job_id=job_id,
            status="processing",
            message="First saving uploaded files and preparing track hub for processing",
            progress=0,
            completed_tracks=0,
            total_tracks=len(track_stanzas),
            track_statuses={},
        )
        uploaded_files_map = {}  # Map track_id → saved file name
        for track_id, file, safe_filename in uploads:
            # Save file to staging area
            if not dry_run:
                file.save(staging_area / safe_filename)
            uploaded_files_map[track_id] = safe_filename

        # Initialize all tracks with None, then populate from map
        for track_stanza in track_stanzas:
            track_id = track_stanza.get("id")
            if not track_id:
                print(f"Warning: Track stanza missing 'id' field. Stanza: {track_stanza}", file=sys.stderr)
                continue
            track_stanza["uploadedFileName"] = uploaded_files_map.get(track_id, None)

        # Update metadata for downstream uses
        metadata["dataset_format"] = "gosling"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=4)

        hub_url = f"{domain_url}/tracks/{dataset_id}"

        # Queue the job
        try:
            queue_trackhub_job(job_id, share_uid, hub_json, assembly, track_stanzas, hub_url, dry_run)
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
                hub_url=hub_url,
                higlass_config=higlass_config,
                dry_run=dry_run,

            )

            result["success"] = result_sync["success"]
            result["message"] = result_sync["message"]
            return result, 200 if result["success"] else 500

        except Exception as e:
            result["message"] = f"Error processing track hub: {str(e)}"
            print(f"TrackHubCopy error: {str(e)}", file=sys.stderr)
            return result, 500
