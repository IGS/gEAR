"""Fake gear.userhistory: UserHistory.add_record logs to the fake database log instead of writing."""

import geardb  # the fake geardb installed by run_cgi.py


class UserHistory:
    def __init__(self, *args, **kwargs):
        pass

    def add_record(self, **kwargs):
        geardb._log("userhistory.add_record", **kwargs)
