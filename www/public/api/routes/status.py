from flask import Blueprint, jsonify

# Create a blueprint for this specific noun
status_bp = Blueprint('status', __name__)

INTERNAL_API_BASE = "http://127.0.0.1:8080/api"

# Route 1: The Plural Collection
@status_bp.route('/status', methods=['GET'])
def get_status():
    # Simple heath check
    return jsonify({"status": "OK", "message": "Public API is up and running!"}), 200