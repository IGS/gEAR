"""
Fake gear.userhistory for account-page CGIs.

UserHistory.get_latest_entries is the real method from lib/gear/userhistory.py, so the query it builds
(including the LIMIT clause) is logged by the fake geardb. Only __init__ is replaced, because the real one
reads gear.ini. add_record logs to the fake database log instead of writing.
"""

import importlib.util
from pathlib import Path

import geardb  # the fake geardb installed by run_cgi.py

_REAL_PATH = Path(__file__).resolve().parents[3] / "lib" / "gear" / "userhistory.py"
_spec = importlib.util.spec_from_file_location("_real_gear_userhistory", _REAL_PATH)
_real = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_real)


class UserHistory(_real.UserHistory):
    def __init__(self, connection=False, user_id=None, entry_category=None, entry_date=None, label=None, url=None):
        self.connection = connection or geardb.Connection()
        self.user_id = user_id
        self.entry_category = entry_category
        self.entry_date = entry_date
        self.label = label
        self.url = url

    def _serialize_json(self):
        return {k: v for k, v in vars(self).items() if k != "connection"}

    def add_record(self, **kwargs):
        geardb._log("userhistory.add_record", **kwargs)


# get_latest_entries() builds its entries with the module-level name, so point it at this class
_real.UserHistory = UserHistory
