"""Fake gear.analysis: an Analysis whose dataset_path is always missing (no analysis files in tests)."""


class Analysis:
    def __init__(self, id=None, dataset_id=None, session_id=None, **kwargs):
        self.id = id
        self.dataset_id = dataset_id
        self.session_id = session_id
        self.dataset_path = f"/nonexistent/analysis/{id}.h5ad"

    def discover_type(self):
        pass
