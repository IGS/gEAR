"""Stand-in for gear.analysis in copy_dataset_analysis.cgi tests: records any use to stderr."""

import sys


def get_analysis(*args, **kwargs):
    print(f"FAKE get_analysis called: {args} {kwargs}", file=sys.stderr)
    raise RuntimeError("fake get_analysis should not be reached")


class Analysis:
    def __init__(self, *args, **kwargs):
        print(f"FAKE Analysis created: {args} {kwargs}", file=sys.stderr)
        self.base_path = "/nonexistent"
        self.settings_path = "/nonexistent/settings.json"
