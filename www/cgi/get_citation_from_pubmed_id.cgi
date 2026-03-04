#!/opt/bin/python3

"""
Get a pubmedId and return the citation as JSON.

Fixes: use correct datetime and sleep functions, portable shebang,
and safer requests call with timeout.
"""

import cgi
import cgitb
from datetime import datetime, timezone
import time
import json
import os
import sys
import requests

lib_path = os.path.abspath(os.path.join("..", "..", "lib"))
sys.path.append(lib_path)
import geardb

# set up cache object for pubmeds (1000 entries)
pubmed_cache = {}

def _now():
    return datetime.now(timezone.utc)


def add_to_cache(pubmed_id, data):
    pubmed_cache[pubmed_id] = {
        "data": data,
        "timestamp": _now(),
    }

    # If cache exceeds 1000 entries, remove the oldest one
    if len(pubmed_cache) > 1000:
        oldest_key = min(
            pubmed_cache.keys(),
            key=lambda k: pubmed_cache[k]["timestamp"],
        )
        del pubmed_cache[oldest_key]


def get_from_cache(pubmed_id):
    cached = pubmed_cache.get(pubmed_id)
    if cached is not None:
        return cached["data"]
    return None


# Throttle API calls to 1 every 0.4 seconds
last_api_call_time = None


def throttle_api_calls():
    global last_api_call_time
    now = _now()
    if last_api_call_time is not None:
        elapsed_time = (now - last_api_call_time).total_seconds()
        if elapsed_time < 0.4:
            # Use time.sleep from the time module; this is correct for throttling
            time.sleep(0.4 - elapsed_time)
    last_api_call_time = _now()


def print_json(obj, status=None):
    if status:
        print(f"Status: {status}")
    print("Content-Type: application/json; charset=utf-8")
    print()
    print(json.dumps(obj))


def main():
    # Parse form data
    form = cgi.FieldStorage()
    pubmed_id = form.getfirst("pubmed_id", "").strip()

    if not pubmed_id:
        print_json({"error": "No pubmed_id provided"}, status="400 Bad Request")
        return

    # Check cache first
    cached_data = get_from_cache(pubmed_id)
    if cached_data is not None:
        print_json(cached_data)
        return

    # If not in cache, fetch from API
    api_url = (
        "https://pmc.ncbi.nlm.nih.gov/api/ctxp/v1/pubmed/"
        f"?format=citation&id={pubmed_id}"
    )

    try:
        throttle_api_calls()
        response = requests.get(
            api_url,
            timeout=10,
            headers={"User-Agent": "pubmed-cgi/1.0"},
        )
        response.raise_for_status()

        try:
            data = response.json()
        except ValueError:
            data = {"raw": response.text}

        add_to_cache(pubmed_id, data)
        print_json(data)
    except requests.RequestException as e:
        print_json({"error": str(e)}, status="502 Bad Gateway")


if __name__ == "__main__":
    main()