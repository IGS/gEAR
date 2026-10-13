"""
cgi_wb_analysis.py - Fake gear.analysis for the workbench step CGI tests (test_cgi_workbench_steps.py).

Analyses live under $GEAR_TEST_ANALYSIS_DIR, laid out like the real ones:

    primary       <dir>/primary/<dataset_id>.h5ad
    public        <dir>/public/<analysis_id>/<dataset_id>.h5ad
    user_saved    <dir>/user_saved/<user_id>/<analysis_id>/<dataset_id>.h5ad   (user from the session)
    user_unsaved  <dir>/user_unsaved/<session_id>/<analysis_id>/<dataset_id>.h5ad

Like the real get_analysis(), a non-primary analysis whose h5ad doesn't exist raises FileNotFoundError.
get_adata() reads the analysis's h5ad if there is one, else (primary only) builds a small AnnData, so
steps can be chained: a step that writes its result is read back by the next one.
"""

import os

import anndata
import geardb
import numpy as np
import pandas as pd

N_CELLS, N_GENES = 60, 30


def build_adata():
    """60 cells in two louvain groups (g1/g2) x 30 genes; G0-G4 are higher in g1, mt-0..mt-2 are mitochondrial."""
    rng = np.random.default_rng(0)
    groups = pd.Categorical(["g1"] * (N_CELLS // 2) + ["g2"] * (N_CELLS // 2))
    X = rng.poisson(3, size=(N_CELLS, N_GENES)).astype(np.float32)
    X[: N_CELLS // 2, :5] += rng.poisson(20, size=(N_CELLS // 2, 5))
    symbols = [f"G{i}" for i in range(N_GENES - 3)] + ["mt-0", "mt-1", "mt-2"]
    obs = pd.DataFrame({"louvain": groups, "n_genes": (X > 0).sum(axis=1)},
                       index=[f"cell{i}" for i in range(N_CELLS)])
    var = pd.DataFrame({"gene_symbol": symbols}, index=[f"ENSG{i:04d}" for i in range(N_GENES)])
    return anndata.AnnData(X=X, obs=obs, var=var)


class Analysis:
    def __init__(self, id=None, dataset_id=None, session_id=None, user_id=None, type=None):
        self.id = id or dataset_id
        self.dataset_id = dataset_id
        self.session_id = session_id
        self.user_id = user_id
        self.type = type or "primary"

    @property
    def base_path(self):
        root = os.environ["GEAR_TEST_ANALYSIS_DIR"]
        if self.type == "primary":
            return os.path.join(root, "primary")
        if self.type == "public":
            return os.path.join(root, "public", self.id)
        if self.type == "user_saved":
            if self.user_id is None:
                raise Exception("no user for a user_saved analysis")
            return os.path.join(root, "user_saved", str(self.user_id), self.id)
        if self.type == "user_unsaved":
            if not self.session_id:
                raise Exception("no session for a user_unsaved analysis")
            return os.path.join(root, "user_unsaved", self.session_id, self.id)
        raise Exception(f"Invalid type {self.type!r}")

    @property
    def dataset_path(self):
        return os.path.join(self.base_path, f"{self.dataset_id}.h5ad")

    def get_adata(self, **kwargs):
        if os.path.exists(self.dataset_path):
            return anndata.read_h5ad(self.dataset_path)
        return build_adata()


class SpatialAnalysis(Analysis):
    pass


def get_analysis(analysis_obj, dataset_id, session_id, is_spatial=False):
    if not analysis_obj:
        return Analysis(dataset_id=dataset_id, type="primary")
    ana = Analysis(id=analysis_obj["id"], dataset_id=dataset_id, session_id=session_id,
                   user_id=geardb.get_user_id_from_session_id(session_id), type=analysis_obj.get("type"))
    if ana.type != "primary" and not os.path.exists(ana.dataset_path):
        raise FileNotFoundError(f"No h5ad file found for analysis {ana.id}")
    return ana


def test_analysis_for_zarr(analysis):
    return False
