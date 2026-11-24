"""gear.primary_analysis: recognizing clusters and embeddings already in an upload's obs."""

import anndata
import numpy as np
import pandas as pd
import pytest
from gear import primary_analysis as pa


def make_adata(columns):
    obs = pd.DataFrame({c: np.arange(3, dtype=float) for c in columns}, index=["c1", "c2", "c3"])
    return anndata.AnnData(obs=obs)


@pytest.mark.parametrize("columns, expected", [
    (["cell_clusters"], "cell_clusters"),
    (["seurat_clusters", "RNA_snn_res.0.8"], "seurat_clusters"),
    (["RNA_snn_res.0.8", "SCT_snn_res.1.2"], "RNA_snn_res.0.8"),
    (["leiden"], "leiden"),
    # Seurat's orig.ident is the sample or project name, not a cell grouping
    (["orig.ident"], None),
    # Existing names keep priority over the new ones
    (["seurat_clusters", "cell_type"], "cell_type"),
    (["condition", "RNA_snn_res"], None),
])
def test_find_cluster_column(columns, expected):
    assert pa.find_cluster_column(columns) == expected


def test_seurat_clusters_become_louvain():
    adata = make_adata(["RNA_snn_res.0.8"])
    assert pa.detect_clustering(adata)
    pa.add_clustering_analysis(adata)
    assert isinstance(adata.obs["louvain"].dtype, pd.CategoricalDtype)


@pytest.mark.parametrize("columns", [["PCA_1", "PCA_2"], ["PC_1", "PC_2"], ["PC1", "PC2"], ["pca_1", "pca_2"]])
def test_pca_columns_copied_to_obsm(columns):
    adata = make_adata(columns)
    assert pa.detect_pca(adata) and not pa.has_pca(adata)
    pa.add_pca_analysis(adata)
    assert adata.obsm["X_pca"].shape == (3, 2)


def test_seurat_uploader_reduction_columns_are_recognized():
    """seuratuploader.reduction_to_metadata names columns after the reduction, e.g. umap_1."""
    assert pa.detect_umap(make_adata(["umap_1", "umap_2"]))
    assert pa.detect_tsne(make_adata(["tsne_1", "tsne_2"]))
    assert pa.detect_pca(make_adata(["pca_1", "pca_2"]))


@pytest.mark.parametrize("columns, expected", [
    (["UMAP_1", "UMAP_2"], True),
    (["cell_clusters"], True),
    (["UMAP_1"], False),
    (["condition"], False),
])
def test_detect_primary_analysis(columns, expected):
    assert pa.detect_primary_analysis(make_adata(columns)) is expected
