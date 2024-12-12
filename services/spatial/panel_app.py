import sys

import panel as pn

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PIL import Image
import base64
from io import BytesIO


import scanpy as sc

import spatialdata as sd
from spatialdata import bounding_box_query
from spatialdata_io.experimental import to_legacy_anndata

from pathlib import Path

spatial_path = Path("/datasets/spatial")

pn.extension('plotly'
            , loading_indicator=True
            , sizing_mode="stretch_both"
            , nthreads=4
            )

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]

orig_x1 = None
orig_x2 = None
orig_y1 = None
orig_y2 = None

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

def make_zoom_fig_callback(event):
    # Fallback to original values if no event
    range_x1 = orig_x1
    range_x2 = orig_x2
    range_y1 = orig_y1
    range_y2 = orig_y2

    if event and "range" in event:
        range_x1 = event["range"]["x"][0]
        range_x2 = event["range"]["x"][1]
        range_y1 = event["range"]["y"][0]
        range_y2 = event["range"]["y"][1]

    #error_row[0].object = f"Selected range: x1={range_x1}, x2={range_x2}, y1={range_y1}, y2={range_y2}"

    # Filter the data based on the selected range
    selected_data = df[(df["spatial1"] >= range_x1) & (df["spatial1"] <= range_x2) & (df["spatial2"] >= range_y1) & (df["spatial2"] <= range_y2)]

    make_zoom_fig(selected_data)

def make_zoom_fig(df):

    range_x1 = df["spatial1"].min()
    range_x2 = df["spatial1"].max()
    range_y1 = df["spatial2"].min()
    range_y2 = df["spatial2"].max()

    image_trace = go.Image(source=base64_string)

    # Make a 3-column subplot
    fig = make_subplots(rows=1, cols=3, subplot_titles=("Gene Expression", "Clusters", "Image Only"), horizontal_spacing=0.15)

    ### Expression plot
    fig.add_trace(image_trace, row=1, col=1)
    fig.add_trace(make_expression_scatter(df, norm_gene_symbol, "YlGn", 2), row=1, col=1)

    ### Cluster plot
    fig.add_trace(image_trace, row=1, col=2)
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
    fig.update_layout(legend=dict(itemsizing='constant', x=0.635), dragmode="select")

    ### Image only
    fig.add_trace(image_trace, row=1, col=3)

    fig.update_xaxes(range=[range_x1, range_x2], title_text="spatial1")
    fig.update_yaxes(range=[range_y2, range_y1], title_text="spatial2")

    # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
    fig.update_layout(
        xaxis=dict(domain=[0, 0.25]),
        xaxis2=dict(domain=[0.4, 0.65]),
        xaxis3=dict(domain=[0.75, 1]),
    )

def make_expression_scatter(df, gene_symbol, color, size):
    return go.Scattergl(x=df["spatial1"], y=df["spatial2"], mode="markers", marker=dict(
                color=df[gene_symbol],
                colorscale=color,  # You can choose any colorscale you like
                size=size,  # Adjust the marker size as needed
                colorbar=dict(
                    len=1,  # Adjust the length of the colorbar
                    thickness=20,  # Adjust the thickness of the colorbar (default is 30)
                    #title=f"{gene_symbol} expression",  # Title for the colorbar
                    x=0.25
                ),
                symbol="square",
                ),
                showlegend=False,
                hovertemplate="Expression: %{marker.color:.2f}<extra></extra>"
            )

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
    adata = to_legacy_anndata(sdata, include_images=True, coordinate_system="downscaled_hires")

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

    dataset_genes = set(adata.var.index.unique())
    norm_gene_symbol = normalize_searched_gene(dataset_genes, pn.state.location.query_params["gene_symbol"])

    selected = adata[:, norm_gene_symbol]
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

    image_trace = go.Image(source=base64_string)

    # Make a 3-column subplot
    fig = make_subplots(rows=1, cols=3, subplot_titles=("Gene Expression", "Clusters", "Image Only"), horizontal_spacing=0.15)

    ### Expression plot
    fig.add_trace(image_trace, row=1, col=1)
    fig.add_trace(make_expression_scatter(df, norm_gene_symbol, "YlGn", 2), row=1, col=1)

    ### Cluster plot
    fig.add_trace(image_trace, row=1, col=2)
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
    fig.update_layout(legend=dict(itemsizing='constant', x=0.635), dragmode="select")

    ### Image only
    fig.add_trace(image_trace, row=1, col=3)

    fig.update_xaxes(range=[0, spatial_img.shape[1]], title_text="spatial1")
    fig.update_yaxes(range=[spatial_img.shape[0], 0], title_text="spatial2")

    # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
    fig.update_layout(
        xaxis=dict(domain=[0, 0.25]),
        xaxis2=dict(domain=[0.4, 0.65]),
        xaxis3=dict(domain=[0.75, 1]),
    )


    fig_pane = pn.pane.Plotly(fig, config={"doubleClick":"reset","displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove})

    # Create a row for the zoomed in view
    zoom_row = pn.Row(
            pn.pane.Plotly(go.Figure(), config={'displayModeBar': False}),
            pn.pane.Plotly(go.Figure(), config={'displayModeBar': False}),
            pn.pane.Plotly(go.Figure(), config={'staticPlot': True}),
            height=400
        )

    # get x and y range for whole figure
    orig_x1, orig_x2 = fig.layout.xaxis.range
    orig_y1, orig_y2 = fig.layout.yaxis.range

    # First time, just show the whole image
    #make_zoom_figs(None)

    error_row = pn.Row(
        pn.pane.Str("Any errors will go here", styles={"color": "red"})
    )

    #pn.bind(make_zoom_figs, expression_plot_pane.param.selected_data, watch=True)
    #pn.bind(make_zoom_figs, cluster_plot_pane.param.selected_data, watch=True)

    # Create the app
    # TODO: Defer loading of the images
    # TODO: explore Datashader
    # TODO: Add some loading indicators for plot drawing
    layout = pn.Column(
        '## Select a region in expression or cluster plot to show zoomed in view',
        pn.Row(
            fig_pane,
            height=400
        ),
        pn.layout.Divider(),    # default margins
        '## Zoomed in view',
        zoom_row
        #, error_row
    ).servable(title="Spatial Data Viewer", location=True)
