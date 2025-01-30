import base64
import io

import geardb
import matplotlib as mpl
mpl.use("Agg")  # Prevents the need for a display when plotting, also thread-safe
import matplotlib.pyplot as plt
import scanpy as sc

import plotly.express as px


from flask import request
from flask_restful import Resource

from gear.plotting import PlotError
from .common import create_projection_adata, get_spatial_adata

sc.settings.verbosity = 0
plt.style.use('default')


"""
This API request will take a selection of genes and create a matplotlib figure with the following plots:

1. Gene expression plots per gene using scanpy.pl.spatial.
2. Cluster annotation plots using scanpy.pl.spatial.
3. (if applicable) Bare image plot.
4. Gene expression plots per gene using scanpy.pl.umap
5. Cluster annotation plots using scanpy.pl.umap.
6. Stacked violin plot of genes x clusters

"""

def map_colors(adata):
    # Assuming df is your DataFrame and it has a column "clusters"
    unique_clusters = adata.obs["clusters"].unique()
    # sort unique clusters by number
    unique_clusters = sorted(unique_clusters, key=lambda x: int(x))
    color_map = None
    if "colors" in adata.obs.columns:
        color_map = {cluster: adata.obs[adata.obs["clusters"] == cluster]["colors"].values[0] for cluster in unique_clusters}
    else:
        color_map = {cluster: px.colors.qualitative.Alphabet[i % len(px.colors.qualitative.Alphabet)] for i, cluster in enumerate(unique_clusters)}

    # Scanpy checks for the existence of adata.uns["cluster_colors"] to determine whether to use the default color map
    adata.uns["clusters_colors"] = [color_map[cluster] for cluster in unique_clusters]

    return adata

def normalize_gene_symbols(dataset_genes, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [g for cg in chosen_genes for g in dataset_genes if cg.lower() == g.lower()]
    return case_insensitive_genes

def create_plot(adata, gene_symbols, colorblind_mode):
    # Set color maps
    spatial_color_map = "YlGn"
    umap_color_map = "YlOrRd"
    if colorblind_mode:
        spatial_color_map = "viridis"
        umap_color_map = "inferno"

    # Split into 2 top/bottom subfigures.  Top side will be spatial and umap plots, and bottom side will be "stacked violin" plot
    io_fig = plt.figure(figsize=(20, 10))
    subfigs = io_fig.subfigures(2,1, height_ratios=[1, 0.5], hspace=0.05)

    extra_plots = 1
    library_id = None
    img_key = None
    spot_size = 10   # Must be provided if no image found.  Otherwise auto-calculated.
    # The spatial function was meant for Visium data loaded through the spatialdata package. It requires a library_id and img_key to show images.
    if (adata.uns["has_images"]):
        extra_plots = 2
        library_id = adata.uns["img_name"]
        img_key = "hires"
        spot_size = None

    # Spatial and UMAP plots
    ax0 = subfigs[0].subplots(2, len(gene_symbols)+extra_plots)
    ax_col = 0

    for gene in gene_symbols:
        if adata.uns["has_images"]:
            sc.pl.spatial(adata, img_key=img_key, color=gene, size=2, library_id=library_id, ax=ax0[0][ax_col], color_map=spatial_color_map, show=False)
        else:
            sc.pl.embedding(adata, basis="spatial", color=gene, size=5, ax=ax0[0][ax_col], color_map=spatial_color_map, show=False)
            # set background color to light blue so that low expression (yellow) points are more visible
            ax0[0][ax_col].set_facecolor("#99FFFF")

        sc.pl.umap(adata, color=gene, ax=ax0[1][ax_col], color_map=umap_color_map, na_color="gray", show=False)

        # remove umap title (using title=None doesn't work)
        ax0[1][ax_col].set_title("")
        ax_col +=1

    # clusters
    if adata.uns["has_images"]:
        sc.pl.spatial(adata, img_key=img_key, color="clusters", size=2, spot_size=spot_size, library_id=library_id, legend_loc=None, ax=ax0[0][ax_col], show=False)
    else:
        sc.pl.embedding(adata, basis="spatial", color="clusters", size=5, legend_loc=None, ax=ax0[0][ax_col], show=False)
        # set background color to light blue so that low expression (yellow) points are more visible
        ax0[0][ax_col].set_facecolor("#99FFFF")

    sc.pl.umap(adata, color="clusters", ax=ax0[1][ax_col], show=False)

    # remove umap title (using title=None doesn't work)
    ax0[1][ax_col].set_title("")

    if adata.uns["has_images"]:
        ax_col +=1

        # blank image
        sc.pl.spatial(adata, img_key=img_key, color=None, size=2, spot_size=spot_size, library_id=library_id, ax=ax0[0][ax_col], show=False)

        # remove axes for ax0[ax_row][1] (no umap blank image)
        ax0[1][ax_col].axis("off")

    # Stacked Violin plot
    ax1 = subfigs[1].subplots(nrows=1, ncols=1)

    if len(gene_symbols) > 1:

        # Currently there is a bug with retrieving categorical x-labels (clusters) after swapping axes.


        violin_fig = sc.pl.stacked_violin(adata, gene_symbols, swap_axes=True, title="Marker gene expression per cluster", groupby=clusters, ax=ax1, row_palette="Blue", show=False, return_fig=True)

        # Remove the existing legend and add a new vertically-oriented one
        # Need to make the figure to have it generate axes.
        violin_fig.add_totals()
        #violin_fig.swap_axes()
        violin_fig.make_figure()


        """
        # At this point, the figure has a spacer above the plot that we need to remove
        # For some reason, deleting all axes and remaking the figure removes the spacer above the plot (which was in ax[2] I think)
        violin_axes = violin_fig.fig.get_axes()
        import sys
        print(violin_axes, file=sys.stderr)

        for ax in violin_axes:
            violin_fig.fig.delaxes(ax)

        violin_fig.swap_axes()
        violin_fig.make_figure()
        """

    else:
        sc.pl.violin(adata, gene_symbols[0], groupby="clusters", ax=ax1, show=False)
        # hide legend
        ax1.get_legend().remove()

    return io_fig

class SpatialScanpyData(Resource):
    """Resource for retrieving spatial comparison data from an analysis.

    Returns
    -------
    Byte stream image data (to be converted back to an image on the client side)
    """
    def post(self, dataset_id):
        req = request.get_json()
        gene_symbols = req.get('gene_symbols', [])
        analysis = req.get('analysis', None)
        projection_id = req.get('projection_id', None)    # projection id of csv output
        colorblind_mode = req.get('colorblind_mode', False)
        high_dpi = req.get('high_dpi', False)

        session_id = request.cookies.get('gear_session_id')

        sc.settings.figdir = '/tmp/'

        success = 1
        message = ""

        if not dataset_id:
            return {
                "success": -1,
                "message": "Request needs a dataset id."
            }

        if not gene_symbols:
            return {
                "success": -1,
                "message": "Request needs at least one gene symbol."
            }

        """ CURRENTLY UNUSED
        try:
            ana = geardb.get_analysis(analysis, dataset_id, session_id)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return {
                "success": -1,
                "message": "Could not retrieve analysis."
            }
        """

        try:
            adata = get_spatial_adata(analysis, dataset_id, session_id)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return {
                "success": -1,
                "message": "Could not retrieve AnnData object."
            }

        if projection_id:
            try:
                adata = create_projection_adata(adata, dataset_id, projection_id)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

        # Filter out cells that overlap with the blank space of the image.
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)

        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        # Reset index to gene symbol to ensure gene names are the plot labels
        adata.var.reset_index(inplace=True)
        adata.var.set_index('gene_symbol', inplace=True)

        # Deduplicate gene_symbols
        adata = adata[:, adata.var.index.duplicated() == False]

        adata.var_names_make_unique()
        dataset_gene_symbols = adata.var.index.tolist()
        gene_symbols = normalize_gene_symbols(dataset_gene_symbols, gene_symbols)

        adata = map_colors(adata)

        io_fig = create_plot(adata, gene_symbols, colorblind_mode)

        with io.BytesIO() as io_pic:
            # Set the saved figure dpi based on the number of observations in the dataset after filtering
            if high_dpi:
                dpi = max(150, int(adata.shape[0] / 100))
                sc.settings.set_figure_params(dpi_save=dpi)
                io_fig.savefig(io_pic, format='png', bbox_inches="tight")
            else:
                # Moved this to the end to prevent any issues with the dpi setting
                sc.settings.set_figure_params(dpi_save=150)
                io_fig.savefig(io_pic, format='webp', bbox_inches="tight")
            io_pic.seek(0)
            plt.close() # Prevent zombie plots, which can cause issues

            image = base64.b64encode(io_pic.read()).decode("utf-8")

        return {
            "success": success,
            "message": message,
            "image": image
        }