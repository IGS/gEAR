import panel as pn

import plotly.express as px
import plotly.graph_objects as go

from PIL import Image
import base64
from io import BytesIO


import scanpy as sc

import spatialdata as sd
from spatialdata import bounding_box_query
from spatialdata_io.experimental import to_legacy_anndata

from pathlib import Path

spatial_path = Path("/datasets/spatial")

pn.extension('plotly')

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]


def make_zoom_figs(event):
        if not event:
            return
        if not "range" in event:
            return

        try:
            range_x1 = event["range"]["x"][0]
            range_x2 = event["range"]["x"][1]
            range_y1 = event["range"]["y"][0]
            range_y2 = event["range"]["y"][1]

            error_row[0].object = f"Selected range: x1={range_x1}, x2={range_x2}, y1={range_y1}, y2={range_y2}"

            # Filter the data based on the selected range
            selected_data = df[(df["spatial1"] >= range_x1) & (df["spatial1"] <= range_x2) & (df["spatial2"] >= range_y1) & (df["spatial2"] <= range_y2)]

            ### Update expression plot
            zoomed_expression_figure = go.Figure()
            zoomed_expression_figure.add_trace(image_trace)
            zoomed_expression_figure.add_trace(make_expression_scatter(selected_data, gene_symbol, 10))

            zoomed_expression_figure.update_xaxes(range=[range_x1, range_x2], title_text="spatial1")
            zoomed_expression_figure.update_yaxes(range=[range_y2, range_y1], title_text="spatial2")

            ### Cluster plot
            zoomed_cluster_figure = go.Figure()
            zoomed_cluster_figure.add_trace(image_trace)
            for cluster in unique_clusters:
                cluster_data = selected_data[selected_data["Clusters"] == cluster]
                zoomed_cluster_figure.add_trace(
                    go.Scattergl(
                        x=cluster_data["spatial1"],
                        y=cluster_data["spatial2"],
                        mode="markers",
                        marker=dict(color=color_map[cluster], size=10, symbol="square"),
                        name=str(cluster),
                        text=cluster_data["Clusters"]
                    )
                )

            # make legend markers bigger
            zoomed_cluster_figure.update_layout(legend=dict(itemsizing='constant'))

            zoomed_cluster_figure.update_xaxes(range=[range_x1, range_x2], title_text="spatial1")
            zoomed_cluster_figure.update_yaxes(range=[range_y2, range_y1], title_text="spatial2")

            ### Image only
            zoomed_image_figure = go.Figure()
            zoomed_image_figure.add_trace(image_trace)
            zoomed_image_figure.update_xaxes(range=[range_x1, range_x2], title_text="spatial1")
            zoomed_image_figure.update_yaxes(range=[range_y2, range_y1], title_text="spatial2")

            # Update the zoom_row figures
            # The object attribute is the plotly figure passed to the pane
            zoom_row[0].object = zoomed_expression_figure
            zoom_row[1].object = zoomed_cluster_figure
            zoom_row[2].object = zoomed_image_figure
        except Exception as e:
            error_row[0].object = f"Error: {e}"


def make_expression_scatter(df, gene_symbol, size):
    return go.Scattergl(x=df["spatial1"], y=df["spatial2"], mode="markers", marker=dict(
                color=df[gene_symbol],
                colorscale='YlGn',  # You can choose any colorscale you like
                size=size,  # Adjust the marker size as needed
                colorbar=dict(
                    title=f"{gene_symbol} expression",  # Title for the colorbar
                ),
                symbol="square",
                ),
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
    sc.pp.filter_cells(adata, min_genes=300)
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

    selected = adata[:, gene_symbol]
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

    ### Expression plot
    expression_plot_figure = go.Figure()
    expression_plot_figure.add_trace(image_trace)
    expression_plot_figure.add_trace(make_expression_scatter(df, gene_symbol, 2))

    expression_plot_figure.update_xaxes(range=[0, spatial_img.shape[1]], title_text="spatial1")
    expression_plot_figure.update_yaxes(range=[spatial_img.shape[0], 0], title_text="spatial2")
    expression_plot_figure.update_layout(title="Gene Expression", title_x=0.5, dragmode="select")

    expression_plot_pane = pn.pane.Plotly(expression_plot_figure, height=400, sizing_mode="stretch_width", config={"doubleClick":"reset","displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove})

    ### Cluster plot
    cluster_plot_figure = go.Figure()
    cluster_plot_figure.add_trace(image_trace)
    for cluster in unique_clusters:
        cluster_data = df[df["Clusters"] == cluster]
        cluster_plot_figure.add_trace(
            go.Scattergl(
                x=cluster_data["spatial1"],
                y=cluster_data["spatial2"],
                mode="markers",
                marker=dict(color=color_map[cluster], size=2, symbol="square"),
                name=str(cluster),
                text=cluster_data["Clusters"]
            )
        )

    # make legend markers bigger
    cluster_plot_figure.update_layout(legend=dict(itemsizing='constant'))

    cluster_plot_figure.update_xaxes(range=[0, spatial_img.shape[1]], title_text="spatial1")
    cluster_plot_figure.update_yaxes(range=[spatial_img.shape[0], 0], title_text="spatial2")
    cluster_plot_figure.update_layout(title="Clusters", title_x=0.5, dragmode="select")

    cluster_plot_pane = pn.pane.Plotly(cluster_plot_figure, height=400, sizing_mode="stretch_width", config={"doubleClick":"reset", "displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove})

    ### Image only
    image_only_figure = go.Figure()
    image_only_figure.add_trace(image_trace)

    image_only_figure.update_xaxes(range=[0, spatial_img.shape[1]], title_text="spatial1")
    image_only_figure.update_yaxes(range=[spatial_img.shape[0], 0], title_text="spatial2")
    image_only_figure.update_layout(title="Image Only", title_x=0.5)

    image_plot_pane = pn.pane.Plotly(image_only_figure, height=400, sizing_mode="stretch_width", config={'staticPlot': True})

    # Create a row for the zoomed in view
    zoom_row = pn.Row(
            pn.pane.Plotly(go.Figure(), height=400, sizing_mode="stretch_width", config={'displayModeBar': False}),
            pn.pane.Plotly(go.Figure(), height=400, sizing_mode="stretch_width", config={'displayModeBar': False}),
            pn.pane.Plotly(go.Figure(), height=400, sizing_mode="stretch_width", config={'staticPlot': True}),
        )

    error_row = pn.Row(
        pn.pane.Str("Any errors will go here", styles={"color": "red"})
    )

    pn.bind(make_zoom_figs, expression_plot_pane.param.selected_data, watch=True)
    pn.bind(make_zoom_figs, cluster_plot_pane.param.selected_data, watch=True)

    # Create the app
    # TODO: Defer loading of the images
    # TODO: Add some loading indicators for plot drawing
    layout = pn.Column(
        '## Select a region in expression or cluster plot to show zoomed in view',
        pn.Row(
            expression_plot_pane,
            cluster_plot_pane,
            image_plot_pane
        ),
        pn.layout.Divider(margin=(-20, 0, 0, 0)),
        '## Zoomed in view',
        zoom_row,
        error_row
    ).servable()
