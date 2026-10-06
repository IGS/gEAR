"""POST /h5ad/<dataset_id>/aggregations: per-category counts with and without obs filters."""

import importlib
import importlib.util
import json
import sys
import types
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
from flask import Flask
from flask_restful import Api

FAKES = Path(__file__).resolve().parent / "fakes"


def _load_stub(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def aggregations_module():
    """
    Import resources.aggregations. If resources.common can't be imported here (it needs plotly
    and anndata shadows), a stub stands in for it while this module's tests run.
    """
    with pytest.MonkeyPatch.context() as mp:
        try:
            importlib.import_module("resources.common")
        except ImportError:
            mp.setitem(sys.modules, "resources.common", _load_stub("resources.common", FAKES / "api_common.py"))
        mp.delitem(sys.modules, "resources.aggregations", raising=False)
        yield importlib.import_module("resources.aggregations")


def make_adata():
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["IHC", "IHC", "OHC", "OHC", "OHC", "Pillar"]),
            "age": pd.Categorical(["P1", "P7", "P1", "P7", "P7", "P1"]),
            "n_genes": [10, 20, 30, 40, 50, 60],     # not categorical: never aggregated or filtered on
        },
        index=[f"cell{i}" for i in range(6)],
    )
    return anndata.AnnData(X=np.ones((6, 3), dtype=np.float32), obs=obs)


@pytest.fixture
def client(aggregations_module, monkeypatch):
    """A test client plus a record of which adata loader the resource used."""
    calls = []
    datasets = {"DS1": types.SimpleNamespace(id="DS1", dtype="single-cell-rnaseq"),
                "SP1": types.SimpleNamespace(id="SP1", dtype="spatial")}

    def loader(name):
        def _load(analysis_id, dataset_id, session_id, *args, **kwargs):
            calls.append(name)
            return make_adata()
        return _load

    mod = aggregations_module
    monkeypatch.setattr(mod.geardb, "get_dataset_by_id", lambda dataset_id, *a, **k: datasets.get(dataset_id))
    for name in ("get_adata_shadow", "get_adata_from_analysis", "get_spatial_adata"):
        monkeypatch.setattr(mod, name, loader(name))

    app = Flask(__name__)
    api = Api(app)
    api.add_resource(mod.Aggregations, "/h5ad/<dataset_id>/aggregations")
    test_client = app.test_client()
    test_client.loader_calls = calls
    return test_client


def post(client, filters=None, dataset_id="DS1"):
    body = {} if filters is None else {"filters": filters}
    resp = client.post(f"/h5ad/{dataset_id}/aggregations", json=body)
    assert resp.status_code == 200, resp.get_data(as_text=True)
    return resp.get_json()


def counts(body, column):
    (agg,) = [a for a in body["aggregations"] if a["name"] == column]
    return {item["name"]: item["count"] for item in agg["items"]}


def assert_all_ints(body):
    assert isinstance(body["total_count"], int)
    for agg in body["aggregations"]:
        assert isinstance(agg["count"], int)
        assert all(isinstance(item["count"], int) for item in agg["items"])
    json.dumps(body)


def test_no_filters(client):
    body = post(client)
    assert body["success"] == 1
    assert body["total_count"] == 6
    assert {a["name"] for a in body["aggregations"]} == {"cell_type", "age"}
    assert counts(body, "cell_type") == {"IHC": 2, "OHC": 3, "Pillar": 1}
    assert counts(body, "age") == {"P1": 3, "P7": 3}
    assert_all_ints(body)
    assert client.loader_calls == ["get_adata_shadow"]


def test_filter_by_age_counts_filtered_cells_and_keeps_zero_categories(client):
    body = post(client, {"age": ["P7"]})
    assert body["total_count"] == 3
    assert counts(body, "cell_type") == {"OHC": 2, "IHC": 1, "Pillar": 0}
    assert counts(body, "age") == {"P7": 3, "P1": 0}
    assert_all_ints(body)
    assert client.loader_calls == ["get_adata_from_analysis"]


def test_filter_by_cell_type(client):
    body = post(client, {"cell_type": ["OHC"]})
    assert body["total_count"] == 3
    assert counts(body, "cell_type") == {"OHC": 3, "IHC": 0, "Pillar": 0}
    assert counts(body, "age") == {"P1": 1, "P7": 2}
    assert_all_ints(body)


def test_filters_on_two_columns_are_combined(client):
    body = post(client, {"cell_type": ["OHC", "IHC"], "age": ["P1"]})
    assert body["total_count"] == 2
    assert counts(body, "cell_type") == {"IHC": 1, "OHC": 1, "Pillar": 0}


@pytest.mark.parametrize("filters", [{"n_genes": [10]}, {"no_such_column": ["x"]}])
def test_filter_on_non_categorical_column_is_ignored(client, filters):
    body = post(client, filters)
    assert body["success"] == 1
    assert body["total_count"] == 6
    assert counts(body, "cell_type") == {"IHC": 2, "OHC": 3, "Pillar": 1}
    assert_all_ints(body)


def test_spatial_dataset_uses_spatial_loader(client):
    body = post(client, dataset_id="SP1")
    assert body["total_count"] == 6
    assert client.loader_calls == ["get_spatial_adata"]


def test_unknown_dataset(client):
    body = post(client, dataset_id="NOPE")
    assert body == {"success": -1, "message": "No dataset found with that ID"}
    assert client.loader_calls == []
