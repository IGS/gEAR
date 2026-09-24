"""get_dataset_comparison.cgi: error responses, obs filters and per-gene list alignment."""

import json
from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
DB = {"datasets": [{"id": "DS1", "share_id": "S1", "owner_id": 1, "dtype": "single-cell-rnaseq"}]}
FAKE_MODULES = {"gear.analysis": FAKES / "comparison_analysis.py"}

N_GENES = 12


def compare(**overrides):
    """Run the CGI comparing cond A (x) against cond B (y) on DS1, with query overrides."""
    query = {
        "dataset_id": "DS1",
        "compare_key": "cond",
        "condition_x": json.dumps(["A"]),
        "condition_y": json.dumps(["B"]),
        "obs_filters": "",
    }
    query.update(overrides)
    query = {k: v for k, v in query.items() if v is not None}
    result = run_cgi("get_dataset_comparison.cgi", query=query, db=DB, fake_modules=FAKE_MODULES)
    assert result.returncode == 0, result.stderr[-2000:]
    return result.json()


def test_unknown_dataset():
    body = compare(dataset_id="NOPE")
    assert body["success"] == 0
    assert body["error"] == body["message"] == "No dataset found with that ID"


def test_no_filters_and_no_transformation():
    body = compare()
    assert body["success"] == 1, body
    assert len(body["gene_ids"]) == N_GENES
    assert body["gene_ids"][0] == "Gene1" and body["symbols"][0] == "G1"
    # G1 is higher in condition A (x)
    assert body["x"][0] > body["y"][0]
    assert body["pvals_adj"] == []
    assert body["compare_key"] == "cond"
    assert body["condition_x"] == ["A"] and body["condition_y"] == ["B"]


@pytest.mark.parametrize("overrides, fragment", [
    ({"obs_filters": "{not json"}, "Invalid JSON"),
    ({"condition_x": "[A"}, "Invalid JSON"),
    ({"compare_key": "nope"}, "Comparison column 'nope' was not found"),
    ({"obs_filters": json.dumps({"tissue": ["x"]})}, "Filter column 'tissue'"),
    # The filter keeps only cond A cells, so condition Y (B) has no cells left
    ({"obs_filters": json.dumps({"cond": ["A"]})}, "No observations match"),
    ({"condition_x": json.dumps(["Z"])}, "No observations match"),
    ({"condition_y": None}, "Please select a condition"),
    ({"condition_y": json.dumps(["A"])}, "identical"),
])
def test_errors_are_single_json_documents(overrides, fragment):
    body = compare(**overrides)
    assert body["success"] == 0
    assert fragment in body["error"]
    assert body["message"] == body["error"]


def test_filter_keeping_both_conditions():
    body = compare(obs_filters=json.dumps({"cond": ["A", "B"]}))
    assert body["success"] == 1, body
    assert len(body["gene_ids"]) == N_GENES


PER_GENE_KEYS = ("gene_ids", "symbols", "fold_changes", "x", "y", "values")


def test_log2_drops_negative_gene_from_every_list():
    body = compare(log_transformation="2")
    assert body["success"] == 1, body
    lengths = {key: len(body[key]) for key in PER_GENE_KEYS}
    assert set(lengths.values()) == {N_GENES - 1}, lengths
    assert "Gene12" not in body["gene_ids"]
    assert "G12" not in body["symbols"]
    # Lists stay aligned: G1 is still first and is still higher in X after the log
    assert body["gene_ids"][0] == "Gene1" and body["symbols"][0] == "G1"
    assert body["x"][0] > body["y"][0]


@pytest.mark.parametrize("statistical_test, log_transformation", [
    ("t-test", None),
    ("wilcoxon", "2"),
])
def test_statistical_tests(statistical_test, log_transformation):
    body = compare(statistical_test=statistical_test, log_transformation=log_transformation)
    assert body["success"] == 1, body
    # G12 has no expressing cells, so filter_genes removes it before ranking
    assert "Gene12" not in body["gene_ids"]
    assert len(body["pvals_adj"]) == len(body["gene_ids"]) == N_GENES - 1
    for key in PER_GENE_KEYS:
        assert len(body[key]) == len(body["gene_ids"]), key
    assert all(0.0 <= p <= 1.0 for p in body["pvals_adj"])


def test_pvals_follow_their_genes():
    """rank_genes_groups orders genes by score; each p-value must stay with its own gene."""
    import importlib.util
    import scanpy as sc

    spec = importlib.util.spec_from_file_location("comparison_analysis", FAKES / "comparison_analysis.py")
    fake = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fake)
    adata = fake.build_adata()
    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.filter_genes(adata, min_cells=1)
    adata.obs["compare"] = adata.obs["cond"].map({"A": "x", "B": "y"}).astype(str)
    sc.tl.rank_genes_groups(adata, "compare", groups=["x"], reference="y", n_genes=0, rankby_abs=False,
                            method="t-test", corr_method="benjamini-hochberg", log_transformed=False)
    ranked = sc.get.rank_genes_groups_df(adata, group="x")
    expected = dict(zip(ranked["names"], ranked["pvals_adj"]))

    body = compare(statistical_test="t-test")
    assert body["success"] == 1, body
    actual = dict(zip(body["gene_ids"], body["pvals_adj"]))
    assert actual == pytest.approx(expected)
    # Rank order differs from gene order in this data, so a positional copy would not match
    assert list(ranked["names"]) != body["gene_ids"]
