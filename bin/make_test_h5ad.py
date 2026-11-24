"""
Makes small synthetic mouse cochlea H5AD files for testing the revised uploader (#390, #422).

    sample_metadata_full.h5ad     every sample metadata column, with missing values and one stand-in
    sample_metadata_none.h5ad     no sample metadata, to see the missing-column and primary analysis warnings

Genes (Ensembl IDs and symbols) are copied from a local mouse dataset so they match gEAR's gene table.
"""

import sys

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

SOURCE = sys.argv[1]      # a local mouse H5AD with a gene_symbol var column
OUT_DIR = sys.argv[2]

rng = np.random.default_rng(390)

source = anndata.read_h5ad(SOURCE, backed="r")
var = source.var[["gene_symbol"]].copy()
source.file.close()
var = var[~var.index.duplicated() & ~var["gene_symbol"].duplicated()]

# Cell groups: (cell_label, cell_ontology, tissue_label, tissue_ontology, marker genes, cells)
GROUPS = [
    ("IHC", "CL:0000589", "OC", "UBERON:0002227", ["Myo7a", "Pou4f3", "Otof", "Cabp2"], 120),
    ("OHC", "CL:0000601", "OC", "UBERON:0002227", ["Myo7a", "Pou4f3", "Slc26a5", "Strc"], 160),
    ("SC", "CL:0002315", "OC", "UBERON:0002227", ["Sox2", "Prox1", "Fgfr3", "S100a1"], 140),
    ("SGN", "CL:0011113", "SG", "UBERON:0000395", ["Tubb3", "Snap25", "Nefl", "Calb2"], 120),
    # An unassigned group whose cell_ontology is a stand-in, which should get a warning
    ("unassigned", "N/A", "OC", "UBERON:0002227", [], 60),
]
markers = sorted({g for group in GROUPS for g in group[4]})
missing = [g for g in markers if g not in set(var["gene_symbol"])]
if missing:
    sys.exit(f"Marker genes not in {SOURCE}: {missing}")

other = var[~var["gene_symbol"].isin(markers)].sample(n=1500, random_state=390)
var = pd.concat([var[var["gene_symbol"].isin(markers)], other])
symbol_to_col = {s: i for i, s in enumerate(var["gene_symbol"])}

rows, obs_rows, umap = [], [], []
centers = [(-6, 4), (-2, 6), (2, 2), (6, -4), (0, -6)]
for (cell_label, cell_ont, tissue_label, tissue_ont, genes, n), center in zip(GROUPS, centers):
    counts = rng.poisson(0.3, size=(n, len(var))).astype(np.float32)
    for g in genes:
        counts[:, symbol_to_col[g]] += rng.poisson(12, size=n)
    rows.append(counts)
    umap.append(rng.normal(center, 0.8, size=(n, 2)))
    for _ in range(n):
        obs_rows.append({
            "cell_label": cell_label,
            "cell_ontology": cell_ont,
            "tissue_label": tissue_label,
            "tissue_ontology": tissue_ont,
            "seurat_clusters": len(rows) - 1,      # numbers, so the column type review step asks about it
        })

X = sp.csr_matrix(np.vstack(rows))
obs = pd.DataFrame(obs_rows)
n_cells = len(obs)
obs.index = [f"cell{i:04d}" for i in range(n_cells)]
umap = np.vstack(umap)

obs.insert(0, "sample_id", obs.index)
obs["treatment"] = rng.choice(["control", "noise exposed"], size=n_cells)
obs["treatment_label"] = obs["treatment"].map({"control": "ctl", "noise exposed": "noise"})
obs["UMAP_1"], obs["UMAP_2"] = umap[:, 0], umap[:, 1]
# Accepted missing values filling whole recommended columns: no warnings, hidden in the curators
obs["age_value"] = "not collected"
obs["age_unit"] = "Not Collected"
# Real values mixed with an accepted missing value: shown in the curators
obs["sex_assigned_at_birth"] = rng.choice(["female", "male", "unavailable"], size=n_cells, p=[0.45, 0.45, 0.1])
# Not a recommended column, so the stand-in isn't warned about, but it's hidden in the curators
obs["notes"] = "N/A"
obs["replicate"] = rng.integers(1, 4, size=n_cells)
obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()

full = anndata.AnnData(X=X, obs=obs, var=var)
full.write_h5ad(f"{OUT_DIR}/sample_metadata_full.h5ad", compression="gzip")

bare = anndata.AnnData(X=X, obs=obs[["replicate", "n_counts"]].copy(), var=var)
bare.write_h5ad(f"{OUT_DIR}/sample_metadata_none.h5ad", compression="gzip")

print(full)
print(bare)
