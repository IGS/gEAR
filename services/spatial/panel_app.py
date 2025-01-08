import sys, logging

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

class SpatialPlot():

    def __init__(self, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, dragmode):
        self.df = df
        self.spatial_img = spatial_img
        self.color_map = color_map
        self.gene_symbol = gene_symbol
        self.fig = go.Figure()
        self.expression_name = expression_name
        self.expression_color = expression_color
        self.dragmode = dragmode
        self.marker_size = 2
        self.base64_string = None
        self.range_x1 = 0
        self.range_x2 = spatial_img.shape[1]
        self.range_y1 = 0
        self.range_y2 = spatial_img.shape[0]

    def make_expression_scatter(self):
        df = self.df
        return go.Scattergl(x=df["spatial1"], y=df["spatial2"], mode="markers", marker=dict(
                    color=df["raw_value"],
                    colorscale=self.expression_color,  # You can choose any colorscale you like
                    size=self.marker_size,  # Adjust the marker size as needed
                    colorbar=dict(
                        len=1,  # Adjust the length of the colorbar
                        thickness=15,  # Adjust the thickness of the colorbar (default is 30)
                        #title=f"{gene_symbol} expression",  # Title for the colorbar
                        x=0.2575
                    ),
                    symbol="square",
                    ),
                    showlegend=False,
                    unselected=dict(marker=dict(opacity=1)),    # Helps with viewing unselected regions
                    hovertemplate="Expression: %{marker.color:.2f}<extra></extra>"
                )

    def add_cluster_traces(self):
        fig = self.fig
        df = self.df
        color_map = self.color_map
        for cluster in color_map:
            cluster_data = df[df["Clusters"] == cluster]
            fig.add_trace(
                go.Scattergl(
                    x=cluster_data["spatial1"],
                    y=cluster_data["spatial2"],
                    mode="markers",
                    marker=dict(color=color_map[cluster], size=self.marker_size, symbol="square"),
                    name=str(cluster),
                    text=cluster_data["Clusters"],
                    unselected=dict(marker=dict(opacity=1)),
                ), row=1, col=2
            )

    def make_fig(self, static_size=False):
        self.create_subplots()
        self.update_axes()
        self.add_image_trace()
        self.fig.add_trace(self.make_expression_scatter(), row=1, col=1)
        self.add_cluster_traces()

        # make legend markers bigger
        self.fig.update_layout(legend=dict(indentation=-15, itemsizing='constant', x=0.585))

        # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
        self.fig.update_layout(
            xaxis=dict(domain=[0, 0.26]),
            xaxis2=dict(domain=[0.32, 0.58]),
            xaxis3=dict(domain=[0.74, 1.00]),
            margin=dict(l=20, r=0, t=50, b=0),
            width=1920 if static_size else None,
            height=1080 if static_size else None,
            # Set dragmode properties on the main figure
            dragmode=self.dragmode,
            selectdirection="d"
        )
        return self.fig

    def convert_image(self):
        # Convert image array to base64 string for performance in loading
        # https://plotly.com/python/imshow/#passing-image-data-as-a-binary-string-to-goimage
        pil_img = Image.fromarray(self.spatial_img)  # PIL image object
        # Crop the image to the selected range
        # https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.crop
        pil_img = pil_img.crop((self.range_x1, self.range_y1, self.range_x2, self.range_y2))
        prefix = "data:image/png;base64,"
        with BytesIO() as stream:
            pil_img.save(stream, format="png")
            self.base64_string = prefix + base64.b64encode(stream.getvalue()).decode("utf-8")

    def create_subplots(self):
        # Make a 3-column subplot
        fig = make_subplots(rows=1, cols=3, column_titles=(f"{self.expression_name} Expression", "Clusters", "Image Only"), horizontal_spacing=0.1)
        self.fig = fig

    def update_axes(self):
        self.fig.update_xaxes(range=[self.range_x1, self.range_x2], title_text="spatial1", showticklabels=False)
        self.fig.update_yaxes(range=[self.range_y2, self.range_y1], title_text="spatial2", showticklabels=False, title_standoff=0)

    def add_image_trace(self):
        self.convert_image()
        image_trace = go.Image(source=self.base64_string, x0=self.range_x1, y0=self.range_y1)
        self.fig.add_trace(image_trace, row=1, col=1)
        self.fig.add_trace(image_trace, row=1, col=2)
        self.fig.add_trace(image_trace, row=1, col=3)


class SpatialNormalSubplot(SpatialPlot):

    def create_pane(self):
        return pn.pane.Plotly(self.make_fig(static_size=True)
                    , config={"doubleClick":"reset","displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

class SpatialZoomSubplot(SpatialPlot):

    def __init__(self, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params):
        super().__init__(df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params)

        # Preserve the original dataframe for filtering
        self.orig_df = df

    def create_pane(self):
        return pn.pane.Plotly(self.make_fig(static_size=False)
                    , config={'displayModeBar': False, 'doubleClick': None}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

    def make_zoom_fig_callback(self, event):
        if event and "range" in event:
            # determine if first or second plot
            x = "x" if "x" in event["range"] else "x2" if "x2" in event["range"] else "x3"
            y = "y" if "y" in event["range"] else "y2" if "y2" in event["range"] else "y3"

            range_x1 = event["range"][x][0]
            range_x2 = event["range"][x][1]
            range_y1 = event["range"][y][0]
            range_y2 = event["range"][y][1]

            df = self.orig_df

            # Filter the data based on the selected range
            self.df = df[(df["spatial1"] >= range_x1) & (df["spatial1"] <= range_x2) & (df["spatial2"] >= range_y1) & (df["spatial2"] <= range_y2)]

            self.range_x1 = range_x1
            self.range_x2 = range_x2
            self.range_y1 = range_y1
            self.range_y2 = range_y2
            self.expression_color = "YlOrRd"
            self.marker_size = 7

            # TODO: Either a) clear selected points on "not this" plot or b) mirror selection on all plots

        return self.create_pane()

class SpatialPanel(pn.viewable.Viewer):

    def __init__(self, dataset_id, gene_symbol, **params):
        super().__init__(**params)
        self.dataset_id = dataset_id
        self.gene_symbol = gene_symbol

        self.prep_adata()
        self.spatial_img = self.adata.uns["spatial"]["spatialdata_hires_image"]["images"]["hires"]
        self.create_gene_df()
        self.map_colors()

        self.fig_subplot = SpatialNormalSubplot(self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, self.norm_gene_symbol, "YlGn", dragmode="select")
        self.zoom_subplot = SpatialZoomSubplot(self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, "Local", "YlOrRd", dragmode=False)

    def __panel__(self):
        self.fig_pane = self.fig_subplot.create_pane()
        self.zoom_pane = pn.bind(self.zoom_subplot.make_zoom_fig_callback, self.fig_pane.param.selected_data, watch=False)

        return pn.Column(
            pn.pane.Markdown('## Select a region to modify zoomed in view in the bottom panel', height=30),
            self.fig_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown('## Zoomed in view',height=30),
            self.zoom_pane,
            width=1168, height=800
        )

    def prep_adata(self):
        zarr_path = spatial_path / f"{self.dataset_id}.zarr"
        if not zarr_path.exists():
            raise ValueError(f"Dataset {self.dataset_id} not found")
        sdata = sd.read_zarr(zarr_path)

        # Filter to only the hires image boundaries
        img_to_use = "spatialdata_hires_image"
        x = len(sdata.images[img_to_use].x)
        y = len(sdata.images[img_to_use].y)

        self.sdata = bounding_box_query(sdata,
                                axes=("x", "y"),
                                min_coordinate=[0, 0],
                                max_coordinate=[x, y],
                                target_coordinate_system="downscaled_hires",
                                filter_table=True,
                                )

        # Create AnnData object
        # Need to include image since the bounding box query does not filter the image data by coordinates
        # Each Image is downscaled (or upscaled) during rendering to fit a 2000x2000 pixels image (downscaled_hires)
        adata = to_legacy_anndata(self.sdata, include_images=True, coordinate_system="downscaled_hires", table_name="table")

        # Filter out cells that overlap with the blank space of the image.
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)

        adata.var_names_make_unique()
        self.adata = adata

    def create_gene_df(self):
        adata = self.adata
        dataset_genes = set(adata.var['gene_symbol'].unique())
        self.norm_gene_symbol = normalize_searched_gene(dataset_genes, self.gene_symbol)

        gene_filter = adata.var.gene_symbol == self.norm_gene_symbol
        selected = adata[:, gene_filter]
        selected.var.index = pd.Index(["raw_value"])
        df = selected.to_df()

        # Add spatial coords from adata.obsm
        X, Y = (0, 1)
        df["spatial1"] = selected.obsm["spatial"].transpose()[X].tolist()
        df["spatial2"] = selected.obsm["spatial"].transpose()[Y].tolist()

        df["Clusters"] = selected.obs["clusters"].astype("category")
        # drop NaN clusters
        self.df = df.dropna(subset=["Clusters"])

    def map_colors(self):
        df = self.df
        # Assuming df is your DataFrame and it has a column "Clusters"
        unique_clusters = df["Clusters"].unique()
        # sort unique clusters by number
        self.unique_clusters = sorted(unique_clusters, key=lambda x: int(x))
        self.color_map = {cluster: px.colors.qualitative.Alphabet[i % len(px.colors.qualitative.Alphabet)] for i, cluster in enumerate(self.unique_clusters)}
        # Map the colors to the clusters
        df["color"] = df["Clusters"].map(self.color_map)
        self.df = df


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

### MAIN APP ###

# Reset logging level to "error" to suppress a bokeh "dropping patch" info message
# https://github.com/bokeh/bokeh/issues/13229
logging.getLogger().setLevel(logging.ERROR)

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
    gene_symbol = pn.state.location.query_params["gene_symbol"]

    # Create the app
    # TODO: Defer loading of the images
    # TODO: explore Datashader
    # TODO: Add some loading indicators for plot drawing
    sp_panel = SpatialPanel(dataset_id, gene_symbol).servable(title="Spatial Data Viewer", location=True)