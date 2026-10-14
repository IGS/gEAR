"""Fake requests (also registered as requests.exceptions): post() logs instead of calling Github."""

import geardb  # the fake geardb installed by run_cgi.py


class HTTPError(Exception):
    pass


class _Response:
    status_code = 201

    def raise_for_status(self):
        pass

    def json(self):
        return {}


def post(url, **kwargs):
    geardb._log("requests.post", url=url, json=kwargs.get("json"))
    return _Response()
