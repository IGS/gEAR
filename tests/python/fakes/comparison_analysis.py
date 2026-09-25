"""
comparison_analysis.py - Fake gear.analysis for get_dataset_comparison.cgi tests.

get_analysis() returns an object whose get_adata() builds a small, deterministic AnnData:
12 genes (Gene1..Gene12, symbols G1..G12) x 8 cells, with obs "cond" = A (4 cells) / B (4 cells).
G1 is much higher in A; G12 is negative everywhere (so its log is undefined and it has no
expressing cells); the other genes are positive noise.
"""

import anndata
import numpy as np
import pandas as pd

N_CELLS = 8
N_GENES = 12


def build_adata():
    rng = np.random.default_rng(1234)
    X = rng.uniform(1.0, 5.0, size=(N_CELLS, N_GENES))
    X[:4, 0] += 20.0                                   # G1 higher in condition A
    X[:, N_GENES - 1] = -rng.uniform(1.0, 2.0, size=N_CELLS)   # G12 negative everywhere

    obs = pd.DataFrame(
        {"cond": pd.Categorical(["A"] * 4 + ["B"] * 4)},
        index=[f"cell{i + 1}" for i in range(N_CELLS)],
    )
    var = pd.DataFrame(
        {"gene_symbol": [f"G{i + 1}" for i in range(N_GENES)]},
        index=[f"Gene{i + 1}" for i in range(N_GENES)],
    )
    return anndata.AnnData(X=X.astype(np.float32), obs=obs, var=var)


class _FakeAnalysis:
    def __init__(self, dataset_id):
        self.dataset_id = dataset_id

    def get_adata(self, **kwargs):
        return build_adata()


def get_analysis(analysis_data, dataset_id, session_id, is_spatial=False):
    return _FakeAnalysis(dataset_id)
