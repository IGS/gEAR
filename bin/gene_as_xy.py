#!/opt/bin/python

"""
gene_as_xy.py - Prototype: scatter plot of one gene's expression against another's.

Plots each cell with the first gene on the x axis and the second on the y axis, colored
by cell type, with a linear regression line per cell type, and saves a PNG. Written for
one dataset (the Kelley P1 cochlea dataset) on a single gEAR instance: the dataset path,
gene symbols, cell types, colors and output path are all hard-coded. Run it from bin/
after editing those values.
"""

import anndata
import numpy as np

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

#plt.style.use('ggplot')

# Load an AnnData object using Muon

path = "../www/datasets/c4f16a12-9e98-47be-4335-b8321282919e.h5ad"  # Kelley P1 cochlea dataset

adata = anndata.read_h5ad(path)

# Filter the adata object to our two genes of interest
gene1 = "Cib2" #"Fgf8"
gene2 = "Oc90" #"Mpc1"
# find the index of the gene in the adata object
gene1_idx = adata.var[adata.var["gene_symbol"] == gene1].index[0]
gene2_idx = adata.var[adata.var["gene_symbol"] == gene2].index[0]

# filter the adata object to only include the two genes
adata = adata[:, [gene1_idx, gene2_idx]]

# convert to a DataFrame
df = adata.to_df()

# Add cell type information to the DataFrame
df["cell_type"] = adata.obs["cell_type"]

# filter to only IHC, OHC, eIHC, and eOHC cell types
#df = df[df["cell_type"].isin(["IHC", "OHC", "eIHC", "eOHC"])]
df = df[df["cell_type"].isin(["IHC", "OHC", "Oc90"])]


# replace column names (adata.var.index) with gene symbols
df.rename(columns={gene1_idx: gene1, gene2_idx: gene2}, inplace=True)

# Create a plot where each axis is a different gene
# and the cell type is the color


# Create separate groups for each "cell type"
grouped = df.groupby("cell_type", observed=True)

colors = {
    "IHC": "darkgreen",
    "OHC": "purple",
    "eIHC": "lightgreen",
    "eOHC": "thistle",
    "IPC": "blue",
    "OPC": "red",
    "Oc90": "orange"
}

# create color dataframe row for each cell type
df["color"] = df["cell_type"].apply(lambda x: colors.get(x, None))

for name, group in grouped:
    print(name)

    legend_label = name + " (n=" + str(len(group)) + ")"

    plt.plot(group[gene1].tolist(), group[gene2].tolist(), label=legend_label, linestyle='', marker='o', markersize=3, color=colors[name])

    # add regression lines
    model = LinearRegression()
    model.fit(group[[gene1]], group[[gene2]])
    r2 = model.score(group[[gene1]], group[[gene2]])
    coefficients = model.coef_
    intercept = model.intercept_

    # Plot the line
    plt.plot(group[[gene1]], intercept+coefficients[0]*group[[gene1]], color=colors[name], lw=1.5)

plt.xlabel(gene1)
plt.ylabel(gene2)

plt.legend(ncol=2, loc='upper right', fontsize=8)

# Save the plot to a file
plt.savefig("/Users/sadkins/gene_as_xy.png")


