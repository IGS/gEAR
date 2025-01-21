import sys, logging

import panel as pn
import param

import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PIL import Image
import base64
from io import BytesIO

import pandas as pd
import scanpy as sc

import spatialdata as sd

from pathlib import Path

lib_path = Path(__file__).resolve().parent.parent.parent / 'lib'
sys.path.append(str(lib_path))
from gear import spatialuploader

spatial_path = Path("/datasets/spatial")

pio.templates.default = "simple_white"  # no gridlines, white background

pn.extension('plotly'
            , loading_indicator=True
            , defer_load=True
            , nthreads=4
            )

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]

class Settings(param.Parameterized):
    gene_symbol = param.String(doc="Gene symbol to display")
    dataset_id = param.String(doc="Dataset ID to display")
    selection_x1 = param.Integer(doc="left selection range")
    selection_x2 = param.Integer(doc="right selection range")
    selection_y1 = param.Integer(doc="upper selection range")
    selection_y2 = param.Integer(doc="lower selection range")

class SpatialPlot():
    """
    Generalized class for creating a spatial plot with a gene expression heatmap, cluster markers, and an image.
    """

    def __init__(self, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, dragmode):
        self.df = df
        self.spatial_img = spatial_img  # None or numpy array
        self.color_map = color_map
        self.gene_symbol = gene_symbol
        self.fig = go.Figure()
        self.expression_name = expression_name
        self.expression_color = expression_color
        self.dragmode = dragmode
        self.marker_size = 2
        self.base64_string = None
        self.range_x1 = 0 if self.spatial_img is not None else min(df["spatial1"])
        self.range_x2 = self.spatial_img.shape[1] if self.spatial_img is not None else max(df["spatial1"])
        self.range_y1 = 0 if self.spatial_img is not None else min(df["spatial2"])
        self.range_y2 = self.spatial_img.shape[0] if self.spatial_img is not None else max(df["spatial2"])

    def make_expression_scatter(self):
        df = self.df
        return go.Scattergl(x=df["spatial1"], y=df["spatial2"], mode="markers", marker=dict(
                    color=df["raw_value"],
                    colorscale=self.expression_color,  # You can choose any colorscale you like
                    size=self.marker_size,  # Adjust the marker size as needed
                    colorbar=dict(
                        len=1,  # Adjust the length of the colorbar
                        thickness=15,  # Adjust the thickness of the colorbar (default is 30)
                        x=self.max_x1 - 0.0025,  # Adjust the x position of the colorbar
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
        # TODO: Add click event for removing equivalent points from expression plot if legend item is clicked
        for cluster in color_map:
            cluster_data = df[df["clusters"] == cluster]
            fig.add_trace(
                go.Scattergl(
                    x=cluster_data["spatial1"],
                    y=cluster_data["spatial2"],
                    mode="markers",
                    marker=dict(color=color_map[cluster], size=self.marker_size, symbol="square"),
                    name=str(cluster),
                    text=cluster_data["clusters"],
                    unselected=dict(marker=dict(opacity=1)),
                ), row=1, col=2
            )

    def make_fig(self, static_size=False):
        self.create_subplots()
        self.update_axes()
        # domain is adjusted whether there are images or not
        if self.spatial_img is not None:
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.26]),
                xaxis2=dict(domain=[0.33, 0.59]),
                xaxis3=dict(domain=[0.74, 1.00]),
            )
        else:
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.4]),
                xaxis2=dict(domain=[0.5, 0.90]),
            )

        # Get max domain for axis 1 and axis 2
        self.max_x1 = self.fig.layout.xaxis.domain[1]
        self.max_x2 = self.fig.layout.xaxis2.domain[1]

        if self.spatial_img is not None:
            self.add_image_trace()
        self.fig.add_trace(self.make_expression_scatter(), row=1, col=1)
        self.add_cluster_traces()

        # make legend markers bigger
        self.fig.update_layout(legend=dict(indentation=-15, itemsizing='constant', x=self.max_x2 + 0.01))

        # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
        self.fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=0),
            width=1920 if static_size else None,
            height=1080 if static_size else None,
            plot_bgcolor="#CCFFFF", # Mirror the background color of visium images
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
        if self.spatial_img is not None:
            self.create_subplots_three()
        else:
            self.create_subplots_two()

    def create_subplots_two(self):
        self.fig = make_subplots(rows=1, cols=2, column_titles=(f"{self.expression_name} Expression", "clusters"), horizontal_spacing=0.1)

    def create_subplots_three(self):
        self.fig = make_subplots(rows=1, cols=3, column_titles=(f"{self.expression_name} Expression", "clusters", "Image Only"), horizontal_spacing=0.1)

    def update_axes(self):
        self.fig.update_xaxes(range=[self.range_x1, self.range_x2], title_text="spatial1", showgrid=False, showticklabels=False, ticks="")
        self.fig.update_yaxes(range=[self.range_y2, self.range_y1], title_text="spatial2", showgrid=False, showticklabels=False, ticks="", title_standoff=0)

    def add_image_trace(self):
        self.convert_image()
        image_trace = go.Image(source=self.base64_string, x0=self.range_x1, y0=self.range_y1)
        self.fig.add_trace(image_trace, row=1, col=1)
        self.fig.add_trace(image_trace, row=1, col=2)
        self.fig.add_trace(image_trace, row=1, col=3)

class SpatialNormalSubplot(SpatialPlot):
    """
    Class for creating a spatial plot with a gene expression heatmap, cluster markers, and an image.
    """

    def create_pane(self):
        return pn.pane.Plotly(self.make_fig(static_size=True)
                    , config={"doubleClick":"reset","displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

class SpatialZoomSubplot(SpatialPlot):
    """
    Class for creating a zoomed-in spatial plot with a gene expression heatmap, cluster markers, and an image.
    """

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
            # TODO: Sometimes the selection does not trigger the callback, need to investigate

        return self.create_pane()

class SpatialPanel(pn.viewable.Viewer):
    """
    Class for prepping the spatial data and setting up the Panel app.
    """

    def __init__(self, dataset_id, gene_symbol, min_genes=200, **params):
        super().__init__(**params)
        self.dataset_id = dataset_id
        self.gene_symbol = gene_symbol

        self.prep_sdata()
        self.prep_adata()

        self.min_genes_slider = pn.widgets.IntSlider(name='Filter - Mininum genes per observation', start=0, end=500, step=5, value=200)
        self.refresh_dataframe(min_genes)

        pn.bind(self.refresh_dataframe, self.min_genes_slider.param.value_throttled, watch=True)


    def __panel__(self):
        return pn.Column(
            pn.pane.Markdown('## Select a region to modify zoomed in view in the bottom panel', height=30),
            self.min_genes_slider,
            self.fig_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown('## Zoomed in view',height=30),
            self.zoom_pane,
            width=1100, height=725
        )

    def prep_sdata(self):
        zarr_path = spatial_path / f"{self.dataset_id}.zarr"
        if not zarr_path.exists():
            raise ValueError(f"Dataset {self.dataset_id} not found")
        sdata = sd.read_zarr(zarr_path)

        try:
            platform = sdata.tables["table"].uns["platform"]
        except KeyError:
            raise ValueError("No platform information found in the dataset")

        # Ensure the spatial data type is supported
        if platform not in spatialuploader.SPATIALTYPE2CLASS.keys():
            print("Invalid or unsupported spatial data type")
            print("Supported types: {0}".format(spatialuploader.SPATIALTYPE2CLASS.keys()))
            sys.exit(1)

        # Use uploader class to determine correct helper functions
        self.spatial_obj = spatialuploader.SPATIALTYPE2CLASS[platform]()
        self.spatial_obj.sdata = sdata
        # Dictates if this will be a 2- or 3-plot figure
        self.has_images = self.spatial_obj.has_images

        # Filter by bounding box (mostly for images)
        self.spatial_obj._filter_sdata_by_coords()

    def prep_adata(self):
        # Create AnnData object
        # Need to include image since the bounding box query does not filter the image data by coordinates
        # Each Image is downscaled (or upscaled) during rendering to fit a 2000x2000 pixels image (downscaled_hires)
        self.spatial_obj._convert_sdata_to_adata()
        adata = self.spatial_obj.adata

        self.spatial_img = None
        if self.has_images:
            img_name = self.spatial_obj.img_name
            self.spatial_img = adata.uns["spatial"][img_name]["images"]["hires"]
        self.adata = adata

    def filter_adata(self):
        # Filter out cells that overlap with the blank space of the image.
        sc.pp.filter_cells(self.adata, min_genes=self.min_genes)
        sc.pp.normalize_total(self.adata, inplace=True)
        sc.pp.log1p(self.adata)
        self.adata.var_names_make_unique()

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

        # Add cluster info
        if "clusters" not in selected.obs:
            raise ValueError("No cluster information found in adata.obs")

        df["clusters"] = selected.obs["clusters"].astype("category")
        # drop NaN clusters
        self.df = df.dropna(subset=["clusters"])

    def refresh_dataframe(self, min_genes):
        self.min_genes = min_genes
        self.filter_adata()
        self.create_gene_df()
        self.map_colors()
        self.fig_subplot = SpatialNormalSubplot(self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, self.norm_gene_symbol, "YlGn", dragmode="select")
        self.zoom_subplot = SpatialZoomSubplot(self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, "Local", "YlOrRd", dragmode=False)

        self.fig_pane = self.fig_subplot.create_pane()
        self.zoom_pane = pn.bind(self.zoom_subplot.make_zoom_fig_callback, self.fig_pane.param.selected_data, watch=False)


    def map_colors(self):
        df = self.df
        # Assuming df is your DataFrame and it has a column "clusters"
        unique_clusters = df["clusters"].unique()
        # sort unique clusters by number
        self.unique_clusters = sorted(unique_clusters, key=lambda x: int(x))
        self.color_map = None
        if "colors" in df:
            self.color_map = {cluster: df[df["clusters"] == cluster]["colors"].values[0] for cluster in self.unique_clusters}
        else:
            self.color_map = {cluster: px.colors.qualitative.Alphabet[i % len(px.colors.qualitative.Alphabet)] for i, cluster in enumerate(self.unique_clusters)}
            # Map the colors to the clusters
            df["color"] = df["clusters"].map(self.color_map)
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

""" TODO
settings = Settings()
pn.state.location.sync(settings, {
    'dataset_id': 'dataset_id'
    , 'gene_symbol': 'gene_symbol'
    , 'selection_x1': 'selection_x1'
    , 'selection_x2': 'selection_x2'
    , 'selection_y1': 'selection_y1'
    , 'selection_y2': 'selection_y2'})
"""

# If not params passed, just show OK as a way to test the app
if not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()
else:
    if not "dataset_id" in pn.state.location.query_params:
        raise ValueError("Please provide a dataset_id")

    if not "gene_symbol" in pn.state.location.query_params:
        raise ValueError("Please provide a gene_symbol")

    dataset_id = pn.state.location.query_params["dataset_id"]
    gene_symbol = pn.state.location.query_params["gene_symbol"]

    # Create the app
    # TODO: Defer loading of the images (https://panel.holoviz.org/how_to/callbacks/defer_load.html)
    # TODO: explore Datashader
    # TODO: Add some loading indicators for plot drawing (https://panel.holoviz.org/how_to/state/busy.html)
    sp_panel = SpatialPanel(dataset_id, gene_symbol).servable(title="Spatial Data Viewer", location=True)