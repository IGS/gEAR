#!/opt/bin/python3

"""
Generate a JupyterHub launch URL for a given dataset and language.

Accepts:
  - share_id: The dataset share ID
  - language: The notebook language ('python' or 'r')

Returns JSON with:
  - url: The full JupyterHub launch URL
  - error: Any error message (if applicable)
"""

import cgi
import json
import os
import sys
from pathlib import Path
from dotenv import dotenv_values
from itsdangerous import URLSafeTimedSerializer

# Add lib path
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


def print_json(obj, status=None):
    """Output JSON response."""
    if status:
        print(f"Status: {status}")
    print("Content-Type: application/json; charset=utf-8\n")
    print(json.dumps(obj))


def load_env_file(env_path):
    """Load environment variables from .env file."""
    if not os.path.exists(env_path):
        return None
    return dotenv_values(env_path)


def get_hub_base_url():
    """Get the JupyterHub base URL from .env or default."""
    jupyter_root = os.path.abspath(os.path.join('..', '..', 'jupyterhub'))
    env_file = os.path.join(jupyter_root, '.env')
    env_vars = load_env_file(env_file)

    if env_vars and 'HUB_PUBLIC_BASE_URL' in env_vars:
        base_url = env_vars['HUB_PUBLIC_BASE_URL'].strip('"')
        return base_url

    # Default fallback
    return "http://localhost:8000"


def get_launch_secret():
    """Get the GEAR_LAUNCH_SECRET from .env file."""
    jupyter_root = os.path.abspath(os.path.join('..', '..', 'jupyterhub'))
    env_file = os.path.join(jupyter_root, '.env')
    env_vars = load_env_file(env_file)

    if env_vars and 'GEAR_LAUNCH_SECRET' in env_vars:
        return env_vars['GEAR_LAUNCH_SECRET'].strip('"')

    return None


def main():
    # Parse form data
    form = cgi.FieldStorage()
    share_id = form.getfirst("share_id", "").strip()
    language = form.getfirst("language", "python").strip().lower()

    # Validate inputs
    if not share_id:
        print_json({"error": "No share_id provided"}, status="400 Bad Request")
        return

    if language not in ["python", "r"]:
        print_json({"error": "Language must be 'python' or 'r'"}, status="400 Bad Request")
        return

    # Get the secret
    secret = get_launch_secret()
    if not secret:
        print_json({"error": "Could not load JupyterHub secret from .env"}, status="500 Internal Server Error")
        return

    # Look up the dataset
    try:
        dataset = geardb.get_dataset_by_share_id(share_id)
        if not dataset:
            print_json({"error": f"Dataset with share_id '{share_id}' not found"}, status="404 Not Found")
            return

        host_datasets_root = os.environ.get("HOST_DATASETS_ROOT", "/var/www/datasets")
        file_path = os.path.join(host_datasets_root, dataset.id + ".h5ad")

    except Exception as e:
        print_json({"error": f"Error looking up dataset: {str(e)}"}, status="500 Internal Server Error")
        return

    # Try to get the current user; fallback to 'anonymous'
    # This would need to be enhanced to get the actual logged-in user
    username = "user"

    try:
        # Create the serializer
        s = URLSafeTimedSerializer(secret, salt="gear-launch")

        # Create the payload
        payload = {
            "username": username,
            "datasets": [file_path],
            "selected_dataset": file_path,
            "notebook_env": language,
        }

        # Generate the token
        token = s.dumps(payload)

        # Build the URL
        hub_base = get_hub_base_url().rstrip('/')
        next_url = f"/hub/user-redirect/lab/tree/gear_starters/{language}_notebook_template.ipynb"
        launch_url = f"{hub_base}/hub/gear-login?launch_token={token}&next={next_url}"

        print_json({"url": launch_url, "success": 1})

    except Exception as e:
        print_json({"error": f"Error generating launch token: {str(e)}"}, status="500 Internal Server Error")
        return


if __name__ == '__main__':
    main()
