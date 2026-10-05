"""POST /analysis/plotTopGenesPCA: the loadings figure goes to the user's unsaved analysis, never a shared directory."""

import importlib
import importlib.util
import sys
import types
from pathlib import Path

import pytest
from flask import Flask
from flask_restful import Api

FAKES = Path(__file__).resolve().parent / "fakes"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def client(tmp_path, monkeypatch):
    scanpy = pytest.importorskip("scanpy")
    try:
        importlib.import_module("resources.common")
    except ImportError:
        monkeypatch.setitem(sys.modules, "resources.common", _load("resources.common", FAKES / "api_common.py"))
    monkeypatch.delitem(sys.modules, "resources.top_pca_genes", raising=False)
    top_pca_genes = importlib.import_module("resources.top_pca_genes")

    workbench = _load("fake_workbench_analysis", FAKES / "workbench_analysis.py")
    monkeypatch.setenv("GEAR_TEST_ANALYSIS_DIR", str(tmp_path))
    monkeypatch.setattr(top_pca_genes, "get_analysis", workbench.get_analysis)
    monkeypatch.setattr(top_pca_genes.geardb, "get_dataset_by_id",
                        lambda dataset_id: types.SimpleNamespace(id=dataset_id, dtype="single-cell-rnaseq"))

    saved_to = []
    monkeypatch.setattr(scanpy.pl, "pca_loadings", lambda adata, **kwargs: saved_to.append(scanpy.settings.figdir))

    app = Flask(__name__)
    Api(app).add_resource(top_pca_genes.TopPCAGenes, "/analysis/plotTopGenesPCA")
    return app.test_client(), saved_to


@pytest.mark.parametrize("analysis_type", ["primary", "public", "user_unsaved"])
def test_loadings_saved_to_unsaved_analysis(client, tmp_path, analysis_type):
    test_client, saved_to = client
    response = test_client.post("/analysis/plotTopGenesPCA", data={
        "dataset_id": "DS1", "analysis_id": "A1", "analysis_type": analysis_type, "session_id": "owner", "pcs": "1,2"})
    assert response.get_json()["success"] == 1
    assert [Path(p) for p in saved_to] == [tmp_path / "user_unsaved" / "figures"]
