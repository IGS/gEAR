import sys

import panel as pn

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PIL import Image
import base64
from io import BytesIO

import pandas as pd
import scanpy as sc

import spatialdata as sd
from spatialdata import bounding_box_query
from spatialdata_io.experimental import to_legacy_anndata

from pathlib import Path

spatial_path = Path("/datasets/spatial")

pn.extension('plotly'
            , loading_indicator=True
            , nthreads=4
            )

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]

def normalize_searched_gene(gene_set, chosen_gene):
    """Convert to case-insensitive version of gene.  Returns None if gene not found in dataset."""
    chosen_gene_lower = chosen_gene.lower()
    for gene in gene_set:
        try:
            if chosen_gene_lower == gene.lower():
                return gene
        except Exception:
            print(gene, file=sys.stderr)
            raise
    return None

def make_expression_scatter(df, color, size):
    return go.Scattergl(x=df["spatial1"], y=df["spatial2"], mode="markers", marker=dict(
                color=df["raw_value"],
                colorscale=color,  # You can choose any colorscale you like
                size=size,  # Adjust the marker size as needed
                colorbar=dict(
                    len=1,  # Adjust the length of the colorbar
                    thickness=15,  # Adjust the thickness of the colorbar (default is 30)
                    #title=f"{gene_symbol} expression",  # Title for the colorbar
                    x=0.2575
                ),
                symbol="square",
                ),
                showlegend=False,
                hovertemplate="Expression: %{marker.color:.2f}<extra></extra>"
            )

def make_zoom_fig_callback(event):
    if event and "range" in event:
        # determine if first or second plot
        x = "x" if "x" in event["range"] else "x2" if "x2" in event["range"] else "x3"
        y = "y" if "y" in event["range"] else "y2" if "y2" in event["range"] else "y3"

        range_x1 = event["range"][x][0]
        range_x2 = event["range"][x][1]
        range_y1 = event["range"][y][0]
        range_y2 = event["range"][y][1]

        # Filter the data based on the selected range
        selected_data = df[(df["spatial1"] >= range_x1) & (df["spatial1"] <= range_x2) & (df["spatial2"] >= range_y1) & (df["spatial2"] <= range_y2)]

        return make_zoom_fig(selected_data)
    # else return original figure
    orig_styles = True
    return make_zoom_fig(df, orig_styles)

def make_zoom_fig(df, orig_styles=False):
    range_x1 = df["spatial1"].min()
    range_x2 = df["spatial1"].max()
    range_y1 = df["spatial2"].min()
    range_y2 = df["spatial2"].max()

    # Make a 3-column subplot
    fig = make_subplots(rows=1, cols=3, column_titles=(f"Local Expression", "Clusters", "Image Only"), horizontal_spacing=0.1)

    marker_size = 2 if orig_styles else 7

    fig.update_xaxes(range=[range_x1, range_x2], title_text="spatial1", showticklabels=False)
    fig.update_yaxes(range=[range_y2, range_y1], title_text="spatial2", showticklabels=False, title_standoff=0)

    ### Image only
    # Add earlier so that the hover will show the expression or cluster traces on top
    image_trace = go.Image(source=base64_string)
    fig.add_trace(image_trace, row=1, col=1)
    fig.add_trace(image_trace, row=1, col=2)
    fig.add_trace(image_trace, row=1, col=3)

    ### Expression plot
    fig.add_trace(make_expression_scatter(df, "YlOrRd", marker_size), row=1, col=1)

    ### Cluster plot
    for cluster in unique_clusters:
        cluster_data = df[df["Clusters"] == cluster]
        fig.add_trace(
            go.Scattergl(
                x=cluster_data["spatial1"],
                y=cluster_data["spatial2"],
                mode="markers",
                marker=dict(color=color_map[cluster], size=marker_size, symbol="square"),
                name=str(cluster),
                text=cluster_data["Clusters"]
            ), row=1, col=2
        )

    # make legend markers bigger
    fig.update_layout(legend=dict(indentation=-15, itemsizing='constant', x=0.585))

    # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
    # NOTE: It seems setting a fixed width and height will override the panel sizing mode after updating, so don't set it here
    fig.update_layout(
        xaxis=dict(domain=[0, 0.26]),
        xaxis2=dict(domain=[0.32, 0.58]),
        xaxis3=dict(domain=[0.74, 1.00]),
        margin=dict(l=20, r=0, t=50, b=0, autoexpand=False),
        # For the zoomed figure, we don't want to allow interaction
        dragmode=False
    )
    return fig

### MAIN APP ###

# If not params passed, just show OK as a way to test the app
if not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()
else:
    if not "dataset_id" in pn.state.location.query_params:
        raise ValueError("Please provide a dataset_id")

    if not "gene_symbol" in pn.state.location.query_params:
        raise ValueError("Please provide a gene_symbol")

    #if not "image_str" in pn.state.location.query_params:
    #    raise ValueError("Please provide an encoded image")

    # Read in zarr file into spatialdata object
    dataset_id = pn.state.location.query_params["dataset_id"]
    zarr_path = spatial_path / f"{dataset_id}.zarr"
    if not zarr_path.exists():
        raise ValueError(f"Dataset {dataset_id} not found")
    sdata = sd.read_zarr(zarr_path)

    # Filter to only the hires image boundaries
    img_to_use = "spatialdata_hires_image"
    x = len(sdata.images[img_to_use].x)
    y = len(sdata.images[img_to_use].y)

    sdata = bounding_box_query(sdata,
                            axes=("x", "y"),
                            min_coordinate=[0, 0],
                            max_coordinate=[x, y],
                            target_coordinate_system="downscaled_hires",
                            filter_table=True,
                            )

    # Create AnnData object
    # Need to include image since the bounding box query does not filter the image data by coordinates
    # Each Image is downscaled (or upscaled) during rendering to fit a 2000x2000 pixels image (downscaled_hires)
    adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires", table_name="table")

    # Filter out cells that overlap with the blank space of the image.
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    adata.var_names_make_unique()

    # Convert image array to base64 string for performance in loading
    # https://plotly.com/python/imshow/#passing-image-data-as-a-binary-string-to-goimage
    spatial_img = adata.uns["spatial"]["spatialdata_hires_image"]["images"]["hires"]
    pil_img = Image.fromarray(spatial_img) # PIL image object
    prefix = "data:image/png;base64,"
    with BytesIO() as stream:
        pil_img.save(stream, format="png")
        base64_string = prefix + base64.b64encode(stream.getvalue()).decode("utf-8")

    gene_symbol = pn.state.location.query_params["gene_symbol"]

    dataset_genes = set(adata.var['gene_symbol'].unique())
    norm_gene_symbol = normalize_searched_gene(dataset_genes, pn.state.location.query_params["gene_symbol"])

    gene_filter = adata.var.gene_symbol == norm_gene_symbol
    selected = adata[:, gene_filter]
    selected.var.index = pd.Index(["raw_value"])
    df = selected.to_df()

    # Add spatial coords from adata.obsm
    X, Y = (0, 1)
    df["spatial1"] = selected.obsm["spatial"].transpose()[X].tolist()
    df["spatial2"] = selected.obsm["spatial"].transpose()[Y].tolist()

    df["Clusters"] = selected.obs["clusters"].astype("category")
    # drop NaN clusters
    df = df.dropna(subset=["Clusters"])

    # Assuming df is your DataFrame and it has a column "Clusters"
    unique_clusters = df["Clusters"].unique()
    # sort unique clusters by number
    unique_clusters = sorted(unique_clusters, key=lambda x: int(x))
    color_map = {cluster: px.colors.qualitative.Alphabet[i % len(px.colors.qualitative.Alphabet)] for i, cluster in enumerate(unique_clusters)}
    # Map the colors to the clusters
    df["color"] = df["Clusters"].map(color_map)

    # Make a 3-column subplot
    fig = make_subplots(rows=1, cols=3, column_titles=(f"{norm_gene_symbol} Expression", "Clusters", "Image Only"), horizontal_spacing=0.1)

    fig.update_xaxes(range=[0, spatial_img.shape[1]], title_text="spatial1", showticklabels=False)
    fig.update_yaxes(range=[spatial_img.shape[0], 0], title_text="spatial2", showticklabels=False, title_standoff=0)

    ### Image only
    # Add earlier so that the hover will show the expression or cluster traces on top
    image_trace = go.Image(source=base64_string)
    fig.add_trace(image_trace, row=1, col=1)
    fig.add_trace(image_trace, row=1, col=2)
    fig.add_trace(image_trace, row=1, col=3)

    ### Expression plot
    fig.add_trace(make_expression_scatter(df, "YlGn", 2), row=1, col=1)

    ### Cluster plot
    for cluster in unique_clusters:
        cluster_data = df[df["Clusters"] == cluster]
        fig.add_trace(
            go.Scattergl(
                x=cluster_data["spatial1"],
                y=cluster_data["spatial2"],
                mode="markers",
                marker=dict(color=color_map[cluster], size=2, symbol="square"),
                name=str(cluster),
                text=cluster_data["Clusters"]
            ), row=1, col=2
        )

    # make legend markers bigger
    fig.update_layout(legend=dict(indentation=-15, itemsizing='constant', x=0.585))

    # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
    fig.update_layout(
        xaxis=dict(domain=[0, 0.26]),
        xaxis2=dict(domain=[0.32, 0.58]),
        xaxis3=dict(domain=[0.74, 1.00]),
        margin=dict(l=20, r=0, t=50, b=0, autoexpand=False),
        width=1920, height=1080,
        # Set dragmode properties on the main figure
        dragmode="select", selectdirection="d"
    )
    fig_pane = pn.pane.Plotly(fig
                              , config={"doubleClick":"reset","displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove}
                              , height=350
                              , sizing_mode="stretch_both"
                              )

    # Create a row for the zoomed in view
    bound_fig = pn.bind(make_zoom_fig_callback, fig_pane.param.selected_data, watch=True)
    zoom_pane = pn.pane.Plotly(bound_fig
                               , config={'displayModeBar': False, 'doubleClick': None}
                               , height=350
                               , sizing_mode="stretch_both"
                               )

    # Create the app
    # TODO: Defer loading of the images
    # TODO: explore Datashader
    # TODO: Add some loading indicators for plot drawing
    layout = pn.Column(
        pn.pane.Markdown('## Select a region to modify zoomed in view in the bottom panel', height=30),
        fig_pane,
        pn.layout.Divider(height=5),    # default margins
        pn.pane.Markdown('## Zoomed in view',height=30),
        zoom_pane
        , width=1168, height=800
    ).servable(title="Spatial Data Viewer", location=True)
