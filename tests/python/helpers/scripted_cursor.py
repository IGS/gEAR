"""A stand-in database cursor for the unified search tests (#1383)."""


class ScriptedCursor:
    """
    Records each query and returns the rows for the first key of `responses` that appears in
    the query text (usually a "FROM <table> <alias>" fragment). Queries matching no key return
    no rows. Every query must have exactly one parameter per %s placeholder.
    """

    def __init__(self, responses=None):
        self.responses = responses or {}
        self.queries = []

    def execute(self, qry, params=()):
        assert qry.count("%s") == len(params), qry
        self.queries.append((qry, list(params)))

    def fetchall(self):
        qry = self.queries[-1][0]
        for text, rows in self.responses.items():
            if text in qry:
                return list(rows)
        return []

    def queries_with(self, text):
        """The (query, params) pairs whose query contains text."""
        return [(qry, params) for qry, params in self.queries if text in qry]
