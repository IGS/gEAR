import requests
from flask import Blueprint, jsonify

# Create a blueprint for this specific noun
datasets_bp = Blueprint('datasets', __name__)

INTERNAL_API_BASE = "http://127.0.0.1:8080/api"
INTERNAL_CGI_BASE = "http://127.0.0.1:8080/cgi"

# Route 1: The Plural Collection
@datasets_bp.route('/datasets', methods=['GET'])
def get_all_datasets():
    # Logic to fetch list of datasets from internal API
    response = requests.get(f"{INTERNAL_API_BASE}/datasets", timeout=3)
    return jsonify(response.json()), response.status_code

# Route 2: The Singular Item
@datasets_bp.route('/datasets/<dataset_id>', methods=['GET'])
def get_single_dataset(dataset_id):
    # Logic to fetch a single dataset and format it for the chatbot
    response = requests.get(f"{INTERNAL_API_BASE}/datasets/{dataset_id}", timeout=3)

    # Add your public "fluff" or transformation layer here
    raw_data = response.json()
    clean_data = {"id": raw_data.get("id"), "title": raw_data.get("name")}

    return jsonify(clean_data), response.status_code