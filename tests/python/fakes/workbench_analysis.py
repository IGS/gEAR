"""
Fake gear.analysis for workbench CGI tests.

Each analysis type keeps its files under $GEAR_TEST_ANALYSIS_DIR/<type>/, so a test can see which
directory a script read from or wrote to. get_adata() returns a small AnnData with two cell groups
stored under every cluster column the workbench scripts use, plus tSNE and UMAP embeddings.
"""

import os
import sys

import anndata
import numpy as np
import pandas as pd


class _Analysis:
    def __init__(self, analysis_obj, dataset_id):
        self.id = (analysis_obj or {}).get("id") or "A1"
        self.type = (analysis_obj or {}).get("type") or "primary"
        self.dataset_id = dataset_id

    @property
    def dataset_path(self):
        return os.path.join(os.environ["GEAR_TEST_ANALYSIS_DIR"], self.type, f"{self.dataset_id}.h5ad")

    def get_adata(self):
        print(f"FAKE get_adata from type {self.type}", file=sys.stderr)
        rng = np.random.default_rng(0)
        n_cells, n_genes = 40, 12
        groups = pd.Categorical(["g1"] * (n_cells // 2) + ["g2"] * (n_cells // 2))
        obs = pd.DataFrame(
            {"louvain": groups, "joint_cluster_round4_annot": groups, "subclass_label": groups},
            index=[f"cell{i}" for i in range(n_cells)],
        )
        var = pd.DataFrame({"gene_symbol": [f"G{i}" for i in range(n_genes)]},
                           index=[f"ENSG{i}" for i in range(n_genes)])
        X = rng.poisson(2, size=(n_cells, n_genes)).astype(np.float32)
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        # Precomputed embeddings, so plotting steps can run without computing them
        adata.obsm["X_tsne"] = rng.normal(size=(n_cells, 2))
        adata.obsm["X_umap"] = rng.normal(size=(n_cells, 2))
        return adata


def get_analysis(analysis_obj, dataset_id, session_id, is_spatial=False):
    return _Analysis(analysis_obj, dataset_id)
