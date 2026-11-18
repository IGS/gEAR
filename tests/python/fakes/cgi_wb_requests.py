"""
cgi_wb_requests.py - Fake `requests` for get_citation_from_pubmed_id.cgi tests (no network).

get() appends the requested URL to $FAKE_REQUESTS_LOG and answers according to $FAKE_REQUESTS_MODE:
    json (default)  200 with the body in $FAKE_REQUESTS_BODY
    text            200 with a body that isn't JSON
    http_error      503 (raise_for_status() raises HTTPError)
    unreachable     raises ConnectionError
"""

import json
import os


class RequestException(Exception):
    pass


class HTTPError(RequestException):
    pass


class ConnectionError(RequestException):   # same name as in requests
    pass


class _Response:
    def __init__(self, status_code, text):
        self.status_code = status_code
        self.text = text

    def raise_for_status(self):
        if self.status_code >= 400:
            raise HTTPError(f"{self.status_code} Server Error")

    def json(self):
        return json.loads(self.text)


def get(url, timeout=None, headers=None, **kwargs):
    with open(os.environ["FAKE_REQUESTS_LOG"], "a") as fh:
        fh.write(json.dumps({"url": url, "timeout": timeout}) + "\n")
    mode = os.environ.get("FAKE_REQUESTS_MODE", "json")
    if mode == "unreachable":
        raise ConnectionError("connection refused")
    if mode == "http_error":
        return _Response(503, "unavailable")
    if mode == "text":
        return _Response(200, "not json")
    return _Response(200, os.environ.get("FAKE_REQUESTS_BODY", "{}"))
