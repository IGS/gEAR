"""
display_requests.py - Fake "requests" for the display-saving CGIs.

post() records each plotting request in the fake database log (event "requests.post") and fails
like an unreachable API, so no image is written and no network is used.
"""

import geardb  # the fake geardb, already loaded by the test runner


class ConnectionError(Exception):
    pass


def post(url, json=None, **kwargs):
    geardb._log("requests.post", url=url, plot_type=(json or {}).get("plot_type"))
    raise ConnectionError(f"fake requests: {url} is not reachable in tests")
