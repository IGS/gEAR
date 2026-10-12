"""Helpers in resources.projectr: path guards, chunking, projection IDs and cached results."""

import hashlib
import importlib
import json
import sys
import types
import uuid

import pandas as pd
import pytest
from werkzeug.exceptions import BadRequest


def _sliced(seq, n):
    return (seq[i:i + n] for i in range(0, len(seq), n))


# Production-only dependencies of resources.projectr, stubbed only if they aren't installed here
OPTIONAL_DEPS = {
    "aiohttp": {},
    "aiohttp_retry": {"ExponentialRetry": object, "RetryClient": object},
    "google": {},
    "google.auth": {},
    "google.auth.transport": {},
    "google.auth.transport.requests": {"Request": object},
    "google.oauth2": {},
    "google.oauth2.id_token": {"fetch_id_token": lambda *a, **k: "token"},
    "more_itertools": {"sliced": _sliced},
}


@pytest.fixture(scope="module")
def projectr_module():
    """Import resources.projectr, with stubs for missing optional dependencies while this module's tests run."""
    with pytest.MonkeyPatch.context() as mp:
        for name, attrs in OPTIONAL_DEPS.items():
            try:
                importlib.import_module(name)
            except ImportError:
                stub = types.ModuleType(name)
                stub.__dict__.update(attrs)
                mp.setitem(sys.modules, name, stub)
                parent, _, child = name.rpartition(".")
                if parent:
                    mp.setattr(sys.modules[parent], child, stub, raising=False)
        mp.delitem(sys.modules, "resources.projectr", raising=False)
        module = importlib.import_module("resources.projectr")
        mp.setitem(sys.modules, "resources.projectr", module)
        yield module


@pytest.fixture
def mod(projectr_module, monkeypatch, tmp_path):
    """resources.projectr with its projection and gene list directories under tmp_path."""
    monkeypatch.setattr(projectr_module, "PROJECTIONS_BASE_DIR", tmp_path / "projections")
    monkeypatch.setattr(projectr_module, "CARTS_BASE_DIR", tmp_path / "carts")
    (tmp_path / "carts").mkdir()
    return projectr_module


BAD_IDS = ["..", "../x", "a/b", "/etc/passwd", "a\\b", ""]


# build_projection_csv_path

@pytest.mark.parametrize(
    "scope, expected",
    [
        ("dataset", "by_dataset/DS1/P1.csv"),
        ("genecart", "by_genecart/DS1/P1.csv"),
        ("pval", "by_dataset/DS1/P1_pval.csv"),
    ],
)
def test_csv_path_for_valid_ids(mod, tmp_path, scope, expected):
    assert mod.build_projection_csv_path("DS1", "P1", scope) == tmp_path / "projections" / expected


@pytest.mark.parametrize("bad", BAD_IDS)
@pytest.mark.parametrize("position", [0, 1, 2])
def test_csv_path_rejects_unsafe_components(mod, bad, position):
    args = ["DS1", "P1", "dataset"]
    args[position] = bad
    with pytest.raises(BadRequest):
        mod.build_projection_csv_path(*args)


# build_projection_json_path

@pytest.mark.parametrize("scope", ["dataset", "genecart"])
def test_json_path_for_valid_id(mod, tmp_path, scope):
    expected = (tmp_path / "projections" / f"by_{scope}" / "DS1" / "projections.json").resolve()
    assert mod.build_projection_json_path("DS1", scope) == expected


@pytest.mark.parametrize("bad", BAD_IDS)
def test_json_path_rejects_unsafe_ids(mod, bad):
    with pytest.raises(ValueError):
        mod.build_projection_json_path(bad, "dataset")


# calculate_chunk_size

@pytest.mark.parametrize(
    "num_genes, num_samples, expected",
    [(20000, 5000, 50), (1000, 10, 1000), (3, 7, 333333), (1_000_000, 4, 1)],
)
def test_chunk_size_keeps_payload_under_limit(mod, num_genes, num_samples, expected):
    size = mod.calculate_chunk_size(num_genes, num_samples)
    assert size == expected
    assert size * num_genes <= 1_000_000


@pytest.mark.parametrize("num_genes", [1_000_001, 2_500_000, 50_000_000])
def test_chunk_size_is_at_least_one_for_huge_gene_counts(mod, num_genes):
    assert mod.calculate_chunk_size(num_genes, 10) == 1


@pytest.mark.parametrize(
    "num_genes, num_samples, match",
    [(0, 10, "no genes"), (-1, 10, "no genes"), (0, 0, "no genes"), (10, 0, "no samples")],
)
def test_chunk_size_rejects_empty_dimensions(mod, num_genes, num_samples, match):
    with pytest.raises(ValueError, match=match):
        mod.calculate_chunk_size(num_genes, num_samples)


# create_new_uuid

def test_uuid_is_deterministic(mod):
    first = mod.create_new_uuid("DS1", "cart1", "nmf", False)
    assert isinstance(first, uuid.UUID)
    assert first == mod.create_new_uuid("DS1", "cart1", "nmf", False)


def test_uuid_matches_ids_already_on_disk(mod):
    # Projection IDs written by earlier runs must keep resolving to the same files
    expected = uuid.UUID(hashlib.md5(b"DS1-cart1-nmf-False").hexdigest())
    assert mod.create_new_uuid("DS1", "cart1", "nmf", False) == expected


@pytest.mark.parametrize(
    "other",
    [("DS2", "cart1", "nmf", False), ("DS1", "cart2", "nmf", False),
     ("DS1", "cart1", "pca", False), ("DS1", "cart1", "nmf", True)],
)
def test_uuid_differs_for_different_inputs(mod, other):
    assert mod.create_new_uuid("DS1", "cart1", "nmf", False) != mod.create_new_uuid(*other)


# chunk_dataframe

def test_chunk_dataframe_splits_columns(mod):
    df = pd.DataFrame([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]], columns=[f"c{i}" for i in range(5)])
    chunks = list(mod.chunk_dataframe(df, 2))
    assert [idx for idx, _ in chunks] == [0, 1, 2]
    assert [list(chunk.columns) for _, chunk in chunks] == [["c0", "c1"], ["c2", "c3"], ["c4"]]
    pd.testing.assert_frame_equal(pd.concat([chunk for _, chunk in chunks], axis=1), df)


def test_chunk_dataframe_larger_than_frame(mod):
    df = pd.DataFrame({"a": [1], "b": [2]})
    chunks = list(mod.chunk_dataframe(df, 10))
    assert len(chunks) == 1
    pd.testing.assert_frame_equal(chunks[0][1], df)


# get_existing_projection_result

def write_dataset_projections(tmp_path, content):
    path = tmp_path / "projections" / "by_dataset" / "DS1" / "projections.json"
    path.parent.mkdir(parents=True)
    path.write_text(content if isinstance(content, str) else json.dumps(content))


BASE_RESULT = {"success": 1, "message": "", "projection_id": "P1"}
CONFIG = {"algorithm": "nmf", "zscore": False, "num_common_genes": 5,
          "num_genecart_genes": 8, "num_dataset_genes": 100}


def test_existing_result_without_projections_file(mod):
    assert mod.get_existing_projection_result("DS1", "cart1", "nmf", False, "P1") == BASE_RESULT


def test_existing_result_with_corrupt_projections_file(mod, tmp_path):
    write_dataset_projections(tmp_path, "{not json")
    assert mod.get_existing_projection_result("DS1", "cart1", "nmf", False, "P1") == BASE_RESULT


def test_existing_result_reports_gene_counts(mod, tmp_path):
    write_dataset_projections(tmp_path, {"cart1": [{**CONFIG, "algorithm": "pca"}, CONFIG]})
    result = mod.get_existing_projection_result("DS1", "cart1", "nmf", False, "P1")
    assert result["num_common_genes"] == 5
    assert result["num_genecart_genes"] == 8
    assert result["num_dataset_genes"] == 100
    assert result["message"] == ("Found 5 common genes between the target dataset (100 genes) "
                                 "and the pattern (8 genes).")


def test_existing_result_treats_missing_zscore_as_false(mod, tmp_path):
    config = {k: v for k, v in CONFIG.items() if k != "zscore"}
    write_dataset_projections(tmp_path, {"cart1": [config]})
    assert mod.get_existing_projection_result("DS1", "cart1", "nmf", False, "P1")["num_common_genes"] == 5
    assert mod.get_existing_projection_result("DS1", "cart1", "nmf", True, "P1") == BASE_RESULT


@pytest.mark.parametrize("genecart_id, algorithm", [("cart2", "nmf"), ("cart1", "fixednmf")])
def test_existing_result_without_matching_config(mod, tmp_path, genecart_id, algorithm):
    write_dataset_projections(tmp_path, {"cart1": [CONFIG]})
    assert mod.get_existing_projection_result("DS1", genecart_id, algorithm, False, "P1") == BASE_RESULT


def test_existing_result_with_no_common_genes_has_no_message(mod, tmp_path):
    write_dataset_projections(tmp_path, {"cart1": [{**CONFIG, "num_common_genes": 0}]})
    result = mod.get_existing_projection_result("DS1", "cart1", "nmf", False, "P1")
    assert result["num_common_genes"] == 0
    assert result["message"] == ""


# init_job_status

def test_init_job_status(mod):
    assert mod.init_job_status("P1") == {"status": "pending", "result": {"projection_id": "P1"}, "error": None}


# create_weighted_loading_df

def test_weighted_loading_df_reads_cart_file(mod, tmp_path):
    (tmp_path / "carts" / "cart.abc123.tab").write_text("gene_id\tgene_sym\tP1\nENSG1\tA\t0.5\n")
    df = mod.create_weighted_loading_df("abc123")
    assert list(df.columns) == ["gene_id", "gene_sym", "P1"]
    assert df.iloc[0].tolist() == ["ENSG1", "A", 0.5]


def test_weighted_loading_df_missing_file(mod):
    with pytest.raises(FileNotFoundError, match="cart.nope.tab"):
        mod.create_weighted_loading_df("nope")


@pytest.mark.parametrize("bad", BAD_IDS + ["x/../../outside"])
def test_weighted_loading_df_rejects_unsafe_ids(mod, tmp_path, bad):
    # Make "x/../../outside" resolve to a real file outside the carts directory
    (tmp_path / "carts" / "cart.x").mkdir()
    (tmp_path / "outside.tab").write_text("gene_id\tgene_sym\tP1\nENSG1\tA\t1\n")
    with pytest.raises(ValueError, match="Invalid gene list identifier"):
        mod.create_weighted_loading_df(bad)


# projectr_callback

def test_callback_fails_cleanly_for_dataset_without_samples(mod, monkeypatch, tmp_path):
    """A dataset with no samples fails the job and releases the lock instead of raising in chunking."""
    anndata = pytest.importorskip("anndata")
    import numpy as np

    job_status_dir = tmp_path / "projections" / "job_status"
    job_status_dir.mkdir(parents=True)
    (tmp_path / "projections" / "by_dataset" / "DS1").mkdir(parents=True)
    monkeypatch.setattr(mod, "JOB_STATUS_DIR", job_status_dir)
    (tmp_path / "carts" / "cart.GC1.tab").write_text("gene_id\tgene_sym\tP1\nENSG1\tA\t0.5\nENSG2\tB\t0.1\n")

    genecart = types.SimpleNamespace(num_genes=2, organism_id=1, get_gene_counts=lambda: None)
    dataset = types.SimpleNamespace(organism_id=1, dtype="single-cell-rnaseq")
    empty_adata = anndata.AnnData(X=np.zeros((0, 2)), var=pd.DataFrame(index=["ENSG1", "ENSG2"]))
    analysis = types.SimpleNamespace(dataset_path=tmp_path / "DS1.h5ad", get_adata=lambda **kwargs: empty_adata)
    monkeypatch.setattr(mod.geardb, "get_gene_cart_by_share_id", lambda share_id: genecart)
    monkeypatch.setattr(mod.geardb, "get_dataset_by_id", lambda dataset_id: dataset)
    monkeypatch.setattr(mod, "get_analysis", lambda *args, **kwargs: analysis)
    monkeypatch.setattr(mod, "calculate_chunk_size", lambda *args: pytest.fail("chunking should not be reached"))

    status = mod.projectr_callback("DS1", "GC1", "P1", "session", "weighted-list", "nmf", False, False)

    assert status["status"] == "failed"
    assert "no samples" in status["error"]
    assert status["result"]["num_common_genes"] == 2
    assert not (tmp_path / "projections" / "by_dataset" / "DS1" / "P1.csv.lock").exists()
    assert json.loads((job_status_dir / "job_P1.json").read_text())["status"] == "failed"
