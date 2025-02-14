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

gear_root = Path(__file__).resolve().parents[2]
lib_path = gear_root.joinpath('lib')

sys.path.append(str(lib_path))
from gear import spatialuploader

spatial_path = gear_root.joinpath('www/datasets/spatial')

pio.templates.default = "simple_white"  # no gridlines, white background

SECS_IN_DAY = 86400
CACHE_EXPIRATION = SECS_IN_DAY * 7  # 7 days

# TODO: explore Datashader for large datasets

# Ignore warnings about plotly GUI events, which propagate to the browser console
import warnings
warnings.filterwarnings('ignore', 'plotly.*unrecognized gui edit.*')

pn.extension('plotly'
            , loading_indicator=True
            , defer_load=True
            , nthreads=4
            )

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]

class Settings(param.Parameterized):
    """
    Settings class for configuring parameters related to gene display and selection ranges.

    Attributes:
        gene_symbol (param.String): Gene symbol to display.
        dataset_id (param.String): Dataset ID to display.
        min_genes (param.Integer): Minimum number of genes per observation, with a default of 200 and bounds between 0 and 500.
        selection_x1 (param.Integer): Left selection range.
        selection_x2 (param.Integer): Right selection range.
        selection_y1 (param.Integer): Upper selection range.
        selection_y2 (param.Integer): Lower selection range.

    Info on Param module can be found at https://param.holoviz.org/
    """

    gene_symbol = param.String(doc="Gene symbol to display")
    dataset_id = param.String(doc="Dataset ID to display")
    min_genes = param.Integer(doc="Minimum number of genes per observation", default=200, bounds=(0, 500))
    selection_x1 = param.Number(doc="left selection range", allow_None=True)
    selection_x2 = param.Number(doc="right selection range", allow_None=True)
    selection_y1 = param.Number(doc="upper selection range", allow_None=True)
    selection_y2 = param.Number(doc="lower selection range", allow_None=True)

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
        # if image is provided, the 0,0 (origin) point is the top left corner
        self.flip = False if self.spatial_img is not None else True

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
        self.create_subplots()  # sets self.fig
        self.update_axes()
        # domain is adjusted whether there are images or not
        if self.spatial_img is not None:
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.26]),
                xaxis2=dict(domain=[0.33, 0.59]),  # Leave room for cluster annotations
                xaxis3=dict(domain=[0.74, 1.00]),
            )
        else:
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.4]),
                xaxis2=dict(domain=[0.5, 0.90]),  # Leave room for cluster annotations
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

        # Mirror the background color of visium images
        # Do not set if image is present as there is a slight padding between data and axes ticks
        plot_bgcolor = None
        if self.spatial_img is None:
            plot_bgcolor = "#99FFFF"

        # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
        self.fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=10),
            width=1920 if static_size else None,
            height=1080 if static_size else None,
            plot_bgcolor=plot_bgcolor,
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
        # flip the y-axis if an image is not present
        if self.flip:
            self.fig.update_yaxes(range=[self.range_y1, self.range_y2], title_text="spatial2", showgrid=False, showticklabels=False, ticks="", title_standoff=0)
        else:
            self.fig.update_yaxes(range=[self.range_y2, self.range_y1], title_text="spatial2", showgrid=False, showticklabels=False, ticks="", title_standoff=0)

    def add_image_trace(self):
        self.convert_image()
        # Note that passing the spatial img array to the "z" property instead of the base64-encoded string causes a big performance hit
        image_trace = go.Image(source=self.base64_string
            , x0=self.range_x1
            , y0=self.range_y1
            )
        self.fig.add_trace(image_trace, row=1, col=1)
        self.fig.add_trace(image_trace, row=1, col=2)
        self.fig.add_trace(image_trace, row=1, col=3)

class SpatialNormalSubplot(SpatialPlot):
    """
    Class for creating a spatial plot with a gene expression heatmap, cluster markers, and an image.
    """

    def __init__(self, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params):
        super().__init__(df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params)

        self.selections_dict = {}

        # if all four selection range values are the same, they were not set in the query params, so use the default values
        if settings.selection_x1 != settings.selection_x2 or settings.selection_y1 != settings.selection_y2:
            # Occasionally, if a selection goes to the edge of the plot,
            # the selection may not be modified and added to query_params,
            # so use the default values if that happens
            self.selections_dict = dict(x0=settings.selection_x1, x1=settings.selection_x2, y0=settings.selection_y1, y1=settings.selection_y2)

    def refresh_spatial_fig(self):
        self.make_fig(static_size=False)
        if self.selections_dict:
            self.fig.add_selection(self.selections_dict, row="all", col="all")
        return self.fig.to_dict()

    def make_umap_plots(self):
        df = self.df
        # Make "clusters" a category so that a legend will be created
        df["clusters"] = df["clusters"].astype("category")
        fig = make_subplots(rows=1, cols=2, column_titles=(f"{self.expression_name} UMAP", "clusters UMAP"), horizontal_spacing=0.1)

        fig.update_layout(
            xaxis=dict(domain=[0, 0.4]),
            xaxis2=dict(domain=[0.5, 0.90]),  # Leave room for cluster annotations
        )

        max_x1 = fig.layout.xaxis.domain[1]
        max_x2 = fig.layout.xaxis2.domain[1]

        fig.add_scatter(col=1, row=1
                        , x=df["UMAP1"]
                        , y=df["UMAP2"]
                        , mode="markers"

                        , marker=dict(color=df["raw_value"]
                            , colorscale="cividis_r"
                            , size=2
                            , colorbar=dict(
                                len=1  # Adjust the length of the colorbar
                                , thickness=15  # Adjust the thickness of the colorbar (default is 30)
                                , x=max_x1  # Adjust the x position of the colorbar
                                )
                            )
                        , showlegend=False
                        , hovertemplate="Expression: %{marker.color:.2f}<extra></extra>"
                        )

        # Process clusters as individual traces so that all show on the legend
        df["clusters"] = df["clusters"].astype("category")
        sorted_clusters = sorted(df["clusters"].unique(), key=lambda x: int(x))
        for cluster in sorted_clusters:
            cluster_data = df[df["clusters"] == cluster]
            fig.add_trace(go.Scattergl(
                x=cluster_data["UMAP1"]
                , y=cluster_data["UMAP2"]
                , mode="markers"
                , marker=dict(
                    color=cluster_data["colors"]
                    , size=2
                    )
                , name=str(cluster)
                , text=cluster_data["clusters"]
                ), col=2, row=1)

        fig.update_xaxes(showgrid=False, showticklabels=False, ticks="", title_text="UMAP1")
        fig.update_yaxes(showgrid=False, showticklabels=False, ticks="", title_text="UMAP2", title_standoff=0)

        # make legend markers bigger
        fig.update_layout(legend=dict(indentation=-15, itemsizing='constant', x=max_x2))

        fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=False,
            selectdirection="d"
        )

        return fig.to_dict()

    def make_violin_plot(self):
        df = self.df
        # Make "clusters" a category so that the violin plot will sort the clusters in order
        df["clusters"] = df["clusters"].astype("category")

        fig = go.Figure()
        fig.update_layout(xaxis=dict(domain=[0, 0.90])) # Leave room for cluster annotations

        xmax = fig.layout.xaxis.domain[1]

        # Process clusters as individual traces so that the violin widths are scaled correctly
        sorted_clusters = sorted(df["clusters"].unique(), key=lambda x: int(x))
        for cluster in sorted_clusters:
            cluster_data = df[df["clusters"] == cluster]
            fig.add_trace(go.Violin(
                x=cluster_data["clusters"]
                , y=cluster_data["raw_value"]
                , fillcolor=self.color_map[cluster]
                , line_color=self.color_map[cluster]
                , marker=dict(
                    color="#000000"
                    , size=1
                    )
                , name=str(cluster)
                ))

        # xaxis needs to be categorical, even for numerical values
        fig.update_xaxes(type="category", title_text="Clusters")
        fig.update_yaxes(title_text="Expression", title_standoff=0)

        fig.update_layout(legend=dict(indentation=-15, itemsizing='constant', x=xmax))


        fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=False,
            selectdirection="d"
        )

        return fig.to_dict()

    def mirror_selection_callback(self, event):
        """
        For a selection event, mirror the selection across all plots.
        """

        self.selections_dict = {}

        if event and "range" in event:
            # determine if first or second plot
            x = "x" if "x" in event["range"] else "x2" if "x2" in event["range"] else "x3"
            y = "y" if "y" in event["range"] else "y2" if "y2" in event["range"] else "y3"

            range_x1 = event["range"][x][0]
            range_x2 = event["range"][x][1]
            range_y1 = event["range"][y][0]
            range_y2 = event["range"][y][1]

            self.selections_dict = dict(x0=range_x1, x1=range_x2, y0=range_y1, y1=range_y2)

        return self.refresh_spatial_fig()


class SpatialZoomSubplot(SpatialPlot):
    """
    Class for creating a zoomed-in spatial plot with a gene expression heatmap, cluster markers, and an image.
    """

    def __init__(self, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params):
        super().__init__(df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params)

        # Preserve the original dataframe for filtering
        self.orig_df = df

        # if all four selection range values are the same, they were not set in the query params, so use the default values
        if settings.selection_x1 != settings.selection_x2 or settings.selection_y1 != settings.selection_y2:
            # Occasionally, if a selection goes to the edge of the plot,
            # the selection may not be modified and added to query_params,
            # so use the default values if that happens
            if settings.selection_x1:
                self.range_x1 = settings.selection_x1
            if settings.selection_x2:
                self.range_x2 = settings.selection_x2
            if settings.selection_y1:
                self.range_y1 = settings.selection_y1
            if settings.selection_y2:
                self.range_y2 = settings.selection_y2

            # Viewing a selection, so increase the marker size
            self.calculate_marker_size()

    def calculate_marker_size(self):
        """Dynamically calculate the marker size based on the range of the selection."""
        # Calculate the range of the selection
        x_range = self.range_x2 - self.range_x1
        y_range = self.range_y2 - self.range_y1

        # Calculate the marker size based on the range of the selection
        # The marker size will scale larger as the range of the selection gets more precise
        self.marker_size = int(2 + 4000 / (x_range + y_range))

    def refresh_spatial_fig(self):
        self.fig =  self.make_fig(static_size=False)
        return self.fig.to_dict()

    def make_zoom_fig_callback(self, event):
        if event and "range" in event:
            # determine if first or second plot
            x = "x" if "x" in event["range"] else "x2" if "x2" in event["range"] else "x3"
            y = "y" if "y" in event["range"] else "y2" if "y2" in event["range"] else "y3"

            self.range_x1 = event["range"][x][0]
            self.range_x2 = event["range"][x][1]
            self.range_y1 = event["range"][y][0]
            self.range_y2 = event["range"][y][1]

            # Viewing a selection, so increase the marker size
            self.calculate_marker_size()

            # If no event, use the default values
            df = self.orig_df

            # Filter the data based on the selected range
            self.df = df[(df["spatial1"] >= self.range_x1) & (df["spatial1"] <= self.range_x2) & (df["spatial2"] >= self.range_y1) & (df["spatial2"] <= self.range_y2)]

            settings.selection_x1 = self.range_x1
            settings.selection_x2 = self.range_x2
            settings.selection_y1 = self.range_y1
            settings.selection_y2 = self.range_y2

            # TODO: Either a) clear selected points on "not this" plot or b) mirror selection on all plots
            # TODO: Sometimes the selection does not trigger the callback, need to investigate

        return self.refresh_spatial_fig()

class SpatialPanel(pn.viewable.Viewer):
    """
    Class for prepping the spatial data and setting up the Panel app.
    """

    def __init__(self, dataset_id, gene_symbol, min_genes, **params):
        super().__init__(**params)

        self.dataset_id = dataset_id
        self.gene_symbol = gene_symbol
        self.min_genes = min_genes

        self.min_genes_slider = pn.widgets.IntSlider(name='Filter - Mininum genes per observation', start=0, end=500, step=25, value=min_genes)

        self.normal_fig = dict(data=[], layout={})
        self.zoom_fig = dict(data=[], layout={})
        self.umap_fig = dict(data=[], layout={})
        self.violin_fig = dict(data=[], layout={})

        # the "figs" are dicts which allow for patching the figure objects
        # See -> https://panel.holoviz.org/reference/panes/Plotly.html#patching
        # You can also replace the Figure directly, but I had occasional issues with returning a figure ih the "bind" function

        self.normal_pane = pn.pane.Plotly(self.normal_fig
                    , config={"doubleClick":"reset","displayModeBar":True, "modeBarButtonsToRemove": buttonsToRemove}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.zoom_pane = pn.pane.Plotly(self.zoom_fig
                    , config={'displayModeBar': False, 'doubleClick': None}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.umap_pane = pn.pane.Plotly(self.umap_fig
                    , config={"displayModeBar": False}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.violin_pane = pn.pane.Plotly(self.violin_fig
                    , config={"displayModeBar": False}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.layout = pn.Column(
            pn.bind(self.init_data)
        )

        self.plot_layout = pn.Column(
            pn.pane.Markdown('## Select a region to modify zoomed in view in the bottom panel', height=30),
            self.min_genes_slider,
            self.normal_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown('## Zoomed in view',height=30),
            self.zoom_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown("## UMAP plots", height=30),
            self.umap_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown("## Violin plot", height=30),
            self.violin_pane,
            width=1100, height=1520
        )

        # SAdkins - Have not quite figured out when to use "watch" but I think it mostly applies when a callback does not return a value
        def refresh_dataframe_callback(value):
            with self.normal_pane.param.update(loading=True) \
                , self.zoom_pane.param.update(loading=True) \
                , self.umap_pane.param.update(loading=True) \
                , self.violin_pane.param.update(loading=True) \
                , self.min_genes_slider.param.update(disabled=True):
                self.refresh_dataframe(value)

        pn.bind(refresh_dataframe_callback, self.min_genes_slider.param.value_throttled, watch=True)

        def selection_callback(event):
            self.zoom_pane.object = self.zoom_fig_obj.make_zoom_fig_callback(event)
            self.normal_pane.object = self.normal_fig_obj.mirror_selection_callback(event)

        pn.bind(selection_callback, self.normal_pane.param.selected_data, watch=True)


    def __panel__(self):
        # This is run when the app is loaded. Return the final layout of the app
        return self.layout

    def loading_indicator(self, label):
        return pn.indicators.LoadingSpinner(
            value=True, name=label, align="center", color="info"
        )

    def init_data(self):
        yield self.loading_indicator("Processing data file...")

        def create_adata_pkg():
            self.prep_sdata()
            self.prep_adata()

            adata_pkg = {"adata": self.adata, "img_name": self.spatial_obj.img_name}
            return adata_pkg

        adata_cache_label = f"{self.dataset_id}_adata"

        # Load the Anndata object (+ image name) from cache or create it if it does not exist, with a 24-hour time-to-live
        adata_pkg = pn.state.as_cached(adata_cache_label, create_adata_pkg, ttl=CACHE_EXPIRATION)

        self.adata = adata_pkg["adata"]
        self.spatial_img = None
        if adata_pkg["img_name"]:
            self.spatial_img = self.adata.uns["spatial"][adata_pkg["img_name"]]["images"]["hires"]

        self.orig_adata = self.adata.copy()  # Preserve the original adata for filtering

        yield self.loading_indicator("Processing data to create plots. This may take a minute...")

        yield self.refresh_dataframe(self.min_genes)

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
        self.adata = self.spatial_obj.adata

    def filter_adata(self):
        # Need to make a copy of the original adata to avoid modifying the original (via inplace=True)
        adata = self.orig_adata.copy()

        # Filter out cells that overlap with the blank space of the image.
        sc.pp.filter_cells(adata, min_genes=self.min_genes)
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)

        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        adata.var_names_make_unique()
        self.adata = adata
        return self.adata

    def create_gene_df(self):
        adata = self.adata
        dataset_genes = set(self.adata.var['gene_symbol'].unique())
        self.norm_gene_symbol = normalize_searched_gene(dataset_genes, self.gene_symbol)
        gene_filter = adata.var.gene_symbol == self.norm_gene_symbol
        selected = adata[:, gene_filter]
        selected.var.index = pd.Index(["raw_value"])
        df = selected.to_df()

        # Add spatial coords from adata.obsm
        X, Y = (0, 1)
        df["spatial1"] = selected.obsm["spatial"].transpose()[X].tolist()
        df["spatial2"] = selected.obsm["spatial"].transpose()[Y].tolist()

        # Add UMAP coords
        df["UMAP1"] = selected.obsm["X_umap"].transpose()[X].tolist()
        df["UMAP2"] = selected.obsm["X_umap"].transpose()[Y].tolist()

        # Add cluster info
        if "clusters" not in selected.obs:
            raise ValueError("No cluster information found in adata.obs")

        df["clusters"] = selected.obs["clusters"].astype("category")

        # drop NaN clusters
        self.df = df.dropna(subset=["clusters"])

    def refresh_dataframe(self, min_genes):
        """
        Refresh the dataframe based on the selected gene and min_genes value.
        Updates Plotly dicts in the Panel app in-place
        """
        self.min_genes = min_genes

        # If the min_genes value has changed... sync to URL
        if self.min_genes != settings.min_genes:
            settings.min_genes = self.min_genes

        adata_subset_cache_label = f"{self.dataset_id}_{self.min_genes}_adata"

        # Load the subset Anndata object from cache or create it if it does not exist, with a 24-hour time-to-live
        self.adata = pn.state.as_cached(adata_subset_cache_label, self.filter_adata, ttl=CACHE_EXPIRATION)

        self.create_gene_df()   # creating the Dataframe is generally fast
        self.map_colors()

        # destroy the old figure objects (to free up memory)
        self.normal_fig_obj = None
        self.zoom_fig_obj = None

        self.normal_fig_obj = SpatialNormalSubplot(self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, self.norm_gene_symbol, "YlGn", dragmode="select")
        self.zoom_fig_obj = SpatialZoomSubplot(self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, "Local", "YlOrRd", dragmode=False)

        self.normal_fig = self.normal_fig_obj.refresh_spatial_fig()

        # The pn.bind function for the zoom callback will not trigger when the normal_fig is refreshed.
        self.zoom_fig = self.zoom_fig_obj.refresh_spatial_fig()

        self.umap_fig = self.normal_fig_obj.make_umap_plots()
        self.violin_fig = self.normal_fig_obj.make_violin_plot()

        self.normal_pane.object = self.normal_fig
        self.zoom_pane.object = self.zoom_fig
        self.umap_pane.object = self.umap_fig
        self.violin_pane.object = self.violin_fig

        return self.plot_layout

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
            df["colors"] = df["clusters"].map(self.color_map)
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
settings = Settings()

# If not params passed, just show OK as a way to test the app
if not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()

else:
    pn.state.location.sync(settings, {
        'dataset_id': 'dataset_id'
        , 'gene_symbol': 'gene_symbol'
        , 'min_genes': 'min_genes'
        , 'selection_x1': 'selection_x1'
        , 'selection_x2': 'selection_x2'
        , 'selection_y1': 'selection_y1'
        , 'selection_y2': 'selection_y2'})

    # Create the app
    sp_panel = SpatialPanel(settings.dataset_id, settings.gene_symbol, settings.min_genes)
    sp_panel.servable(title="Spatial Data Viewer", location=True)
