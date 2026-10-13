"""Regression tests for closed GitHub issues in gear.analysis and gear.anndata_processor."""

import json
import types

import gear.anndata_processor as ap
import pandas as pd
import pytest
from gear import analysis

# ---------------------------------------------------------------- #1339 stored analysis JSON


@pytest.fixture
def pipeline_json(tmp_path, monkeypatch):
    """Write a <dataset_id>.pipeline.json and point Analysis.settings_path at it."""
    path = tmp_path / "DS1.pipeline.json"
    monkeypatch.setattr(analysis.Analysis, "settings_path", property(lambda self: path))

    def _write(data: dict) -> dict:
        path.write_text(json.dumps(data))
        return json.loads(repr(analysis.Analysis(dataset_id="DS1")))

    return _write


def test_analysis_json_group_labels_dict_becomes_list(pipeline_json):
    """Regression for #1339: group_labels stored as a dict comes back as a list of its values."""
    data = pipeline_json({"group_labels": {"0": "hair cells", "1": "support cells"}})
    assert data["group_labels"] == ["hair cells", "support cells"]
    assert data["success"] == 1


def test_analysis_json_group_labels_list_is_unchanged(pipeline_json):
    """Regression for #1339: group_labels already stored as a list is returned as is (it used to raise, breaking "retrieve stored analysis")."""
    data = pipeline_json({"group_labels": ["hair cells", "support cells"], "user_session_id": "abc"})
    assert data["group_labels"] == ["hair cells", "support cells"]
    # The saving user's session ID is never sent back (it would let anyone viewing a public analysis log in as them)
    assert "analysis_session_id" not in data and "user_session_id" not in data
    assert "abc" not in json.dumps(data)


def test_analysis_json_without_group_labels(pipeline_json):
    """Regression for #1339: a stored analysis without group_labels loads."""
    assert pipeline_json({"tsne": {"tsne_calculated": True}}) == {"tsne": {"tsne_calculated": True}, "success": 1}


# ---------------------------------------------------------------- #1090 invisible characters in expression values


def test_clean_chunk_strips_whitespace_and_zero_fills_unparseable_values():
    """Regression for #1090: clean_chunk strips surrounding whitespace (including non-breaking spaces) before
    converting to numbers; values that still don't parse, such as one with an escape sequence inside, become 0."""
    chunk = pd.DataFrame(
        {"c1": [" 4.662665156", "1 "], "c2": ["5.903697\x1b[O318", "\t2.5 "], "c3": [6.455984369, "abc"]},
        index=pd.Index(["ENSG00000068781", "ENSG2"], name="gene"),
    )
    cleaned = ap.clean_chunk(chunk)
    assert cleaned.to_dict("index") == {
        "ENSG00000068781": {"c1": 4.662665156, "c2": 0.0, "c3": 6.455984369},
        "ENSG2": {"c1": 1.0, "c2": 2.5, "c3": 0.0},
    }
    assert all(dtype.kind == "f" for dtype in cleaned.dtypes)


def test_expression_matrix_with_invisible_characters_is_read_with_cleanup(tmp_path, monkeypatch):
    """Regression for #1090: an expression matrix whose values contain invisible characters still loads
    (the reader retries with clean_chunk) instead of failing the upload."""
    monkeypatch.setattr(ap, "UPLOADS_BASE_DIR", tmp_path)
    staging = tmp_path / "sess" / "SHARE"
    staging.mkdir(parents=True)
    matrix = staging / "expression.tab"
    matrix.write_text(
        "gene\tcellA\tcellB\tcellC\n"
        "ENSG00000068781\t4.662665156\t5.903697\x1b[O318\t6.455984369\n"
        "ENSG2\t 1.5 \t2\t0\n"
    )
    processor = ap.AnndataProcessor("job-1", "SHARE", staging, staging / "status.json", "DS")
    x = processor._read_expression_matrix_chunks(matrix, chunk_size=1, total_rows=2).toarray()
    assert x.tolist() == [[4.662665156, 0.0, 6.455984369], [1.5, 2.0, 0.0]]


# ---------------------------------------------------------------- UserHistory.get_latest_entries

class _Cursor:
    def __init__(self):
        self.calls = []

    def execute(self, qry, args):
        self.calls.append((" ".join(qry.split()), list(args)))

    def __iter__(self):
        return iter([("2026-10-01", "gene_search", "Sox2", "/p?g=Sox2")])


#@pytest.mark.parametrize("category, expected_where, expected_args", [
#    (None, "WHERE user_id = %s ORDER BY", [7, 5]),
#    ("gene_search", "WHERE user_id = %s AND entry_category = %s ORDER BY", [7, "gene_search", 5]),
#])
@pytest.mark.skip(reason="need to mock the 'history' section of the ini file")
def test_latest_entries_query(category, expected_where, expected_args):
    """The category filter goes before ORDER BY (it was appended after it, which is invalid SQL),
    filters on the requested category, and LIMIT is a query parameter."""
    from gear.userhistory import UserHistory

    cursor = _Cursor()
    history = UserHistory.__new__(UserHistory)
    history.user_id = 7
    history.entry_category = None
    history.connection = types.SimpleNamespace(get_cursor=lambda: cursor)

    entries = history.get_latest_entries(entry_category=category)

    ((qry, args),) = cursor.calls
    assert expected_where in qry and qry.endswith("ORDER BY entry_date DESC LIMIT %s")
    assert args == expected_args
    assert [e.label for e in entries] == ["Sox2"]
