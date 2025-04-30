import sys, logging, tempfile, shutil

# Holoviz imports
import panel as pn
import param
import colorcet as cc
import datashader as ds
#import datashader.transfer_functions as tf

# Plotly imports
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PIL import Image
import base64
from io import BytesIO

import numpy as np
import pandas as pd
import scanpy as sc

import spatialdata as sd
import anndata

from pathlib import Path

from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[2]
lib_path = gear_root.joinpath('lib')

www_path = gear_root.joinpath('www')
PROJECTIONS_BASE_DIR = www_path.joinpath('projections')

sys.path.append(str(lib_path))
from gear.spatialuploader import SPATIALTYPE2CLASS

spatial_path = gear_root.joinpath('www/datasets/spatial')

pio.templates.default = "simple_white"  # no gridlines, white background

SECS_IN_DAY = 86400
CACHE_EXPIRATION = SECS_IN_DAY * 7  # 7 days

# TODO: explore Datashader for large datasets
# TODO: Encountering a weird bug where if I have multiple gEAR pages open, and I interact with a Pane from one tab, plots in the other tabs disappear. Doesn't affect simulatenous users
# ? Maybe make a button to create UMAP instead of creating at the start

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
zoomButtonsToRemove = buttonsToRemove + ["select2d"]

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
    min_genes = param.Integer(doc="Minimum number of genes per observation", default=0, bounds=(0, 500))
    projection_id = param.String(doc="Projection ID to display", allow_None=True)
    selection_x1 = param.Number(doc="left selection range", allow_None=True)
    selection_x2 = param.Number(doc="right selection range", allow_None=True)
    selection_y1 = param.Number(doc="upper selection range", allow_None=True)
    selection_y2 = param.Number(doc="lower selection range", allow_None=True)
    save = param.Boolean(doc="If true, save this configuration as a new display.", default=False)
    display_name = param.String(doc="Display name for the saved configuration", allow_None=True)
    make_default = param.Boolean(doc="If true, make this the default display.", default=False)
    expanded = param.Boolean(doc="If true, show the full range of plots.", default=False)

class SpatialPlot():
    """
    Generalized class for creating a spatial plot with a gene expression heatmap, cluster markers, and an image.
    """

    # TODO: Change "expression" to "relative expression" for projections

    def __init__(self, settings, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, dragmode):
        self.settings = settings
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

        # https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.experimental.to_legacy_anndata.html
        # (see section Matching of spatial coordinates and pixel coordinates)
        # The downscaled image (x,y) coordinate matches (spatial2, spatial1) in the dataframe
        # We traditionally use the 1-coord as x and the 2-coord as y, so we need to address these in various plots.
        self.range_x1 = 0 if self.spatial_img is not None else min(df["spatial1"])
        self.range_x2 = self.spatial_img.shape[1] if self.spatial_img is not None else max(df["spatial1"])
        self.range_y1 = 0 if self.spatial_img is not None else min(df["spatial2"])
        self.range_y2 = self.spatial_img.shape[0] if self.spatial_img is not None else max(df["spatial2"])


    def make_expression_scatter(self):
        df = self.df

        cvs = ds.Canvas(plot_width=300, plot_height=300)

        # If multiple data points are in the same location, use the max gene expression value (raw_value)
        # Oddly, we need to swap the x and y coordinates to match the image, but this is not a problem for the scatter plot
        agg = cvs.points(df, x='spatial2', y='spatial1', agg=ds.max('raw_value'))

        # NOTE: I attempted to use tf.Shade for an image (as I've seen in various examples) but the image was not appearing in the plot

        # There is an agg.to_dataframe(name="raw_value") method, but I elected to use what Copilot suggested
        # Extract x and y coordinates
        x_coords = agg.coords['spatial1'].values
        y_coords = agg.coords['spatial2'].values

        # Flatten the data and create a mask for non-NaN values
        agg_values = agg.values.flatten()
        x_flat = np.repeat(x_coords, len(y_coords))
        y_flat = np.tile(y_coords, len(x_coords))

        # Filter out NaN values (if any)
        mask = ~np.isnan(agg_values)
        x_filtered = x_flat[mask]
        y_filtered = y_flat[mask]
        values_filtered = agg_values[mask]

        return go.Scatter(
                x=x_filtered,
                y=y_filtered,
                mode="markers",
                marker=dict(
                    cauto=False,  # Do not adjust the color scale
                    cmin=0,  # No expression
                    cmax=df["raw_value"].max(),  # Max expression value
                    # ? This would not apply if data is log-transformed
                    color=values_filtered,
                    colorscale=self.expression_color,  # You can choose any colorscale you like
                    size=self.marker_size,  # Adjust the marker size as needed
                    colorbar=dict(
                        len=1,  # Adjust the length of the colorbar
                        thickness=15,  # Adjust the thickness of the colorbar (default is 30)
                        x=self.max_x1 - 0.005,  # Adjust the x position of the colorbar
                    )
                ),
                showlegend=False,
                unselected=dict(marker=dict(opacity=1)),    # Helps with viewing unselected regions
                hovertemplate="Expression: %{marker.color:.2f}<extra></extra>"
            )

    def add_cluster_traces(self, df):
        fig = self.fig
        df = df.sort_values(by="raw_value")

        cvs = ds.Canvas(plot_width=300, plot_height=300)
        agg = cvs.points(df, x='spatial2', y='spatial1', agg=ds.by('clusters', ds.any()))

        # ? Is there a more efficient way to do this?
        agg_df = agg.to_dataframe(name="clusters")
        # Drop rows where "clusters" is False
        agg_df = agg_df[agg_df["clusters"] != False]

        # Move spatial coordinates as columns
        agg_df = agg_df.reset_index(level=["spatial1", "spatial2"])
        # Drop the True/False column and set clusters as a column
        df = agg_df.drop("clusters", axis=1).reset_index()

        color_map = self.color_map

        # TODO: Add click event for removing equivalent points from expression plot if legend item is clicked
        for cluster in color_map:
            cluster_data = df[df["clusters"] == cluster]
            fig.add_trace(
                go.Scattergl(
                    x=cluster_data["spatial1"],
                    y=cluster_data["spatial2"],
                    mode="markers",
                    marker=dict(color=color_map[cluster], size=self.marker_size),
                    name=str(cluster),
                    text=cluster_data["clusters"],
                    unselected=dict(marker=dict(opacity=1)),
                    hovertemplate="Cluster: %{text}<extra></extra>"
                ), row=1, col=self.cluster_col
            )

    def make_fig(self, static_size=False):
        self.create_subplots()  # sets self.fig
        self.update_axes()

        # domain is adjusted whether there are images or not
        if self.spatial_img is not None:
            self.expression_col = 2
            self.cluster_col = 3
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.29]),
                xaxis2=dict(domain=[0.33, 0.62]),  # Leave room for cluster annotations
                xaxis3=dict(domain=[0.71, 1]),
            )
        else:
            self.expression_col = 1
            self.cluster_col = 2
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.45]),
                xaxis2=dict(domain=[0.55, 1]),  # Leave room for cluster annotations
            )

        # Get max domain for axis 1 and axis 2
        if self.expression_col == 1:
            self.max_x1 = self.fig.layout.xaxis.domain[1]
        else:
            self.max_x1 = self.fig.layout.xaxis2.domain[1]

        if self.spatial_img is not None:
            self.add_image_trace()
        self.fig.add_trace(self.make_expression_scatter(), row=1, col=self.expression_col)

        df = self.df
        self.add_cluster_traces(df)

        # Get longest cluster name for legend
        longest_cluster = max(df["clusters"].astype(str).apply(len))
        font_size = 12
        if longest_cluster > 10:
            font_size = 10

        # make legend markers bigger
        self.fig.update_layout(legend=dict(
            font=dict(size=font_size)
            , itemsizing='constant'
            ))

        # Mirror the background color of visium images
        # Do not set if image is present as there is a slight padding between data and axes ticks
        plot_bgcolor = None
        if self.spatial_img is None:
            plot_bgcolor = "#000000"

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

    def encode_pil_image(self, pil_img):
        prefix = "data:image/png;base64,"
        with BytesIO() as stream:
            pil_img.save(stream, format="png")
            return prefix + base64.b64encode(stream.getvalue()).decode("utf-8")

    def convert_background_image(self):
        # Convert image array to base64 string for performance in loading
        # https://plotly.com/python/imshow/#passing-image-data-as-a-binary-string-to-goimage
        pil_img = Image.fromarray(self.spatial_img)  # PIL image object
        # Crop the image to the selected range
        # https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.crop
        pil_img = pil_img.crop((self.range_x1, self.range_y1, self.range_x2, self.range_y2))

        self.base64_string = self.encode_pil_image(pil_img)

    def create_subplots(self):
        if self.spatial_img is not None:
            self.create_subplots_three()
        else:
            self.create_subplots_two()

    def create_subplots_two(self):
        self.fig = make_subplots(rows=1, cols=2, column_titles=(f"{self.expression_name} Expression", "Clusters"), horizontal_spacing=0.1)

    def create_subplots_three(self):
        self.fig = make_subplots(rows=1, cols=3, column_titles=("Image Only", f"{self.expression_name} Expression", "Clusters"), horizontal_spacing=0.1)

    def update_axes(self):
        self.fig.update_xaxes(range=[self.range_x1, self.range_x2], title_text="spatial1", showgrid=False, showticklabels=False, ticks="")
        # y-axis needs to be flipped. 0,0 is the top left corner
        self.fig.update_yaxes(range=[self.range_y2, self.range_y1], title_text="spatial2", showgrid=False, showticklabels=False, ticks="", title_standoff=0)

    def add_image_trace(self):
        self.convert_background_image()
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

    def __init__(self, settings, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params):
        super().__init__(settings, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params)

        self.selections_dict = {}

        # if all four selection range values are the same, they were not set in the query params, so use the default values
        if self.settings.selection_x1 != self.settings.selection_x2 or self.settings.selection_y1 != self.settings.selection_y2:
            # Occasionally, if a selection goes to the edge of the plot,
            # the selection may not be modified and added to query_params,
            # so use the default values if that happens
            self.selections_dict = dict(x0=self.settings.selection_x1, x1=self.settings.selection_x2, y0=self.settings.selection_y1, y1=self.settings.selection_y2)

    def refresh_spatial_fig(self):
        self.make_fig(static_size=False)
        if self.selections_dict:
            self.fig.add_selection(self.selections_dict, row="all", col="all")
        return self.fig.to_dict()

    def make_umap_plots(self):
        # NOTE: Axes are mirrored from previous versions but the coordinates do not matter... clustering is the focus
        df = self.df.sort_values(by="raw_value")

        # Make "clusters" a category so that a legend will be created
        df["clusters"] = df["clusters"].astype("category")

        fig = make_subplots(rows=1, cols=2, column_titles=(f"{self.expression_name} UMAP", "clusters UMAP"), horizontal_spacing=0.1)

        fig.update_layout(
            xaxis=dict(domain=[0, 0.45]),
            xaxis2=dict(domain=[0.55, 1]),  # Leave room for cluster annotations
        )

        max_x1 = fig.layout.xaxis.domain[1]

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

        # Assuming df is your DataFrame and it has a column "clusters"
        unique_clusters = df["clusters"].unique()
        # sort unique clusters by number if numerical, otherwise by name
        try:
            sorted_clusters = sorted(unique_clusters, key=lambda x: int(x))
        except:
            sorted_clusters = sorted(unique_clusters, key=lambda x: str(x))

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

        # Get longest cluster name for legend
        longest_cluster = max(df["clusters"].astype(str).apply(len))
        font_size = 12
        if longest_cluster > 10:
            font_size = 10

        # make legend markers bigger
        fig.update_layout(legend=dict(
            font=dict(size=font_size)
            , itemclick="toggleothers"
            , itemsizing='constant'
            ))

        fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=False,
            selectdirection="d"
        )

        return fig.to_dict()

    def make_violin_plot(self):
        df = self.df.sort_values(by="raw_value")

        # Make "clusters" a category so that the violin plot will sort the clusters in order
        df["clusters"] = df["clusters"].astype("category")

        fig = go.Figure()

        # Assuming df is your DataFrame and it has a column "clusters"
        unique_clusters = df["clusters"].unique()
        # sort unique clusters by number if numerical, otherwise by name
        try:
            sorted_clusters = sorted(unique_clusters, key=lambda x: int(x))
        except:
            sorted_clusters = sorted(unique_clusters, key=lambda x: str(x))

        # Process clusters as individual traces so that the violin widths are scaled correctly
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
        fig.update_yaxes(title_text="Expression", rangemode="tozero")

        # Get longest cluster name for legend
        longest_cluster = max(df["clusters"].astype(str).apply(len))
        font_size = 12
        if longest_cluster > 10:
            font_size = 10

        fig.update_layout(legend=dict(
            font=dict(size=font_size)
            , itemsizing='constant'
            ))

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

    def __init__(self, settings, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params):
        super().__init__(settings, df, spatial_img, color_map, gene_symbol, expression_name, expression_color, **params)

        # Preserve the original dataframe for filtering
        self.orig_df = df

        # if all four selection range values are the same, they were not set in the query params, so use the default values
        if self.settings.selection_x1 != self.settings.selection_x2 or self.settings.selection_y1 != self.settings.selection_y2:
            # Occasionally, if a selection goes to the edge of the plot,
            # the selection may not be modified and added to query_params,
            # so use the default values if that happens
            if self.settings.selection_x1:
                self.range_x1 = self.settings.selection_x1
            if self.settings.selection_x2:
                self.range_x2 = self.settings.selection_x2
            if self.settings.selection_y1:
                self.range_y1 = self.settings.selection_y1
            if self.settings.selection_y2:
                self.range_y2 = self.settings.selection_y2

            # Viewing a selection, so increase the marker size
            self.calculate_marker_size()

    def calculate_marker_size(self):
        """Dynamically calculate the marker size based on the range of the selection."""
        # Calculate the range of the selection
        x_range = self.range_x2 - self.range_x1
        y_range = self.range_y2 - self.range_y1

        # Calculate the marker size based on the range of the selection
        # The marker size will scale larger as the range of the selection gets more precise
        self.marker_size = int(1 + 2500 / (x_range + y_range))

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

            self.settings.selection_x1 = self.range_x1
            self.settings.selection_x2 = self.range_x2
            self.settings.selection_y1 = self.range_y1
            self.settings.selection_y2 = self.range_y2

        return self.refresh_spatial_fig()

class SpatialPanel(pn.viewable.Viewer):
    """
    Class for prepping the spatial data and setting up the Panel app.
    """

    def __init__(self, settings, **params):
        super().__init__(**params)

        self.settings = settings

        pn.state.location.sync(self.settings, {
            'dataset_id': 'dataset_id'
            , 'gene_symbol': 'gene_symbol'
            , 'min_genes': 'min_genes'
            , 'projection_id': 'projection_id'
            , 'selection_x1': 'selection_x1'
            , 'selection_x2': 'selection_x2'
            , 'selection_y1': 'selection_y1'
            , 'selection_y2': 'selection_y2'
            , 'save': 'save'
            , 'display_name': 'display_name'
            , 'make_default': 'make_default'
            , 'expanded': 'expanded'})

        self.dataset_id = self.settings.dataset_id
        self.gene_symbol = self.settings.gene_symbol
        self.min_genes = self.settings.min_genes
        self.projection_id = self.settings.projection_id

        # ? This can be useful for filtering datasets even for projections, but how best to word it?
        self.min_genes_slider = pn.widgets.IntSlider(name='Filter - Mininum genes per observation', start=0, end=500, step=25, value=self.min_genes)

        self.display_name = pn.widgets.TextInput(name='Display name', placeholder='Name this display to save...', width=250)
        self.save_button = pn.widgets.Button(name="Save settings", button_type="primary", width=100, align="end")
        self.make_default = pn.widgets.Checkbox(name='Make this the default display', value=False)

        # NOTE: We will perform the entire computation whether expanded or not, to avoid having to recompute the data
        self.expanded = self.settings.expanded

        self.normal_fig = dict(data=[], layout={})
        self.zoom_fig = dict(data=[], layout={})
        self.umap_fig = dict(data=[], layout={})
        self.violin_fig = dict(data=[], layout={})

        # the "figs" are dicts which allow for patching the figure objects
        # See -> https://panel.holoviz.org/reference/panes/Plotly.html#patching
        # You can also replace the Figure directly, but I had occasional issues with returning a figure ih the "bind" function

        self.normal_pane = pn.pane.Plotly(self.normal_fig
                    , config={"doubleClick":"reset","displayModeBar": True, "modeBarButtonsToRemove": buttonsToRemove}
                    , height=350 if self.expanded else None
                    , sizing_mode="stretch_width" if self.expanded else "stretch_both"
                    )

        self.zoom_pane = pn.pane.Plotly(self.zoom_fig
                    , config={"displayModeBar": True, "modeBarButtonsToRemove": zoomButtonsToRemove, 'doubleClick': None}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.umap_pane = pn.pane.Plotly(self.umap_fig
                    , config={"displayModeBar": True, "modeBarButtonsToRemove": zoomButtonsToRemove}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.violin_pane = pn.pane.Plotly(self.violin_fig
                    , config={"displayModeBar": True, "modeBarButtonsToRemove": zoomButtonsToRemove}
                    , height=350
                    , sizing_mode="stretch_width"
                    )

        self.layout = pn.Column(
            pn.bind(self.init_data)
        )

        layout_height = 312 # 360px - tile header height
        if self.expanded:
            layout_height = 1520

        self.nonexpanded_pre = pn.pane.Markdown(
            "## Click the Expand icon in the top right corner to see all plots",
            height=30, visible=True if not self.expanded else False
        )

        self.expanded_pre = pn.Column(
            pn.Row(
                pn.pane.Markdown('## Select a region to modify zoomed in view in the bottom panel', height=30, width=700),
                self.display_name,
                self.save_button,
            ),
            pn.Row(
                self.min_genes_slider,
                pn.Spacer(width=400),
                self.make_default
            ),
            visible=True if self.expanded else False,
        )

        self.expanded_post = pn.Column(
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown('## Zoomed in view',height=30),
            self.zoom_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown("## UMAP plots", height=30),
            self.umap_pane,
            pn.layout.Divider(height=5),    # default margins
            pn.pane.Markdown("## Violin plot", height=30),
            self.violin_pane,
            visible=True if self.expanded else False,
        )

        # A condensed layout should only have the normal pane.
        self.plot_layout = pn.Column(
            self.expanded_pre,
            self.nonexpanded_pre,
            self.normal_pane,
            self.expanded_post,

            width=1100, height=layout_height
        )

        # SAdkins - Have not quite figured out when to use "watch" but I think it mostly applies when a callback does not return a value
        def refresh_dataframe_callback(value):
            self.dataset_adata = self.dataset_adata_orig.copy()
            self.adata = self.adata_orig.copy()

            with self.normal_pane.param.update(loading=True) \
                , self.zoom_pane.param.update(loading=True) \
                , self.umap_pane.param.update(loading=True) \
                , self.violin_pane.param.update(loading=True) \
                , self.min_genes_slider.param.update(disabled=True):
                self.refresh_dataframe(value)

            with self.umap_pane.param.update(loading=True) \
                , self.min_genes_slider.param.update(disabled=True):
                self.add_umap()

        pn.bind(refresh_dataframe_callback, self.min_genes_slider.param.value_throttled, watch=True)

        def selection_callback(event):
            self.zoom_pane.object = self.zoom_fig_obj.make_zoom_fig_callback(event)
            self.normal_pane.object = self.normal_fig_obj.mirror_selection_callback(event)

        pn.bind(selection_callback, self.normal_pane.param.selected_data, watch=True)

        def save_settings_callback(event):
            self.settings.save = True
            self.settings.display_name = self.display_name.value
            self.settings.make_default = self.make_default.value
            # Make button "is-loading" while saving
            self.save_button.button_type = "success"
            self.save_button.name = "Saving..."
            self.save_button.disabled = True
            print(self.settings, file=sys.stderr)

        def reset_save_button_callback(event):
            if not event.name == "save":
                return
            self.settings.save = False
            self.save_button.button_type = "primary"
            self.save_button.name = "Save settings"
            self.save_button.disabled = False

        pn.bind(save_settings_callback, self.save_button, watch=True)
        self.settings.param.watch(reset_save_button_callback, "save", onlychanged=True)


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
            try:
                self.prep_sdata()
                adata = self.prep_adata()
            except ValueError as e:
                raise

            adata_pkg = {"adata": adata, "img_name": self.spatial_obj.img_name}
            return adata_pkg

        adata_cache_label = f"{self.dataset_id}_adata"

        # Load the Anndata object (+ image name) from cache or create it if it does not exist, with a 1-week time-to-live
        try:
            adata_pkg = pn.state.as_cached(adata_cache_label, create_adata_pkg, ttl=CACHE_EXPIRATION)
        except ValueError as e:
            yield pn.pane.Alert(f"Error: {e}", alert_type="danger")
            raise


        self.dataset_adata_orig = adata_pkg["adata"] # Original dataset adata
        self.dataset_adata = self.dataset_adata_orig.copy() # Copy we manipulate (filtering, etc.)
        self.adata = self.dataset_adata.copy()  # adata object to use for plotting

        # Modify the adata object to use the projection ID if it exists
        if self.projection_id:
            self.adata = self.create_projection_adata()

        self.adata_orig = self.adata.copy() # This is to restore when the min_genes slider is changed

        self.spatial_img = None
        if adata_pkg["img_name"]:
            # In certain conditions, the image multi-array may need to be squeezed so that PIL can read it
            self.spatial_img = self.adata.uns["spatial"][adata_pkg["img_name"]]["images"]["hires"].squeeze()

        yield self.loading_indicator("Processing data to create plots. This may take a minute...")

        yield self.refresh_dataframe(self.min_genes)

        with self.umap_pane.param.update(loading=True) \
            , self.min_genes_slider.param.update(disabled=True):
            self.add_umap()

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
        if platform not in SPATIALTYPE2CLASS.keys():
            print("Invalid or unsupported spatial data type")
            print("Supported types: {0}".format(SPATIALTYPE2CLASS.keys()))
            sys.exit(1)

        # Use uploader class to determine correct helper functions
        self.spatial_obj = SPATIALTYPE2CLASS[platform]()
        self.spatial_obj.sdata = sdata
        # Dictates if this will be a 2- or 3-plot figure
        self.has_images = self.spatial_obj.has_images

        # Filter by bounding box (mostly for images)
        self.spatial_obj._filter_sdata_by_coords()

    def prep_adata(self):
        # Create AnnData object
        # Need to include image since the bounding box query does not filter the image data by coordinates
        # Each Image is downscaled (or upscaled) during rendering to fit a 2000x2000 pixels image (downscaled_hires)
        try:
            self.spatial_obj._convert_sdata_to_adata()
        except Exception as e:
            print("Error converting sdata to adata:", str(e), file=sys.stderr)
            raise ValueError("Could not convert sdata to adata. Check the spatial data type and bounding box coordinates.")
        return self.spatial_obj.adata

    def create_projection_adata(self):
        dataset_adata = self.adata
        dataset_id = self.dataset_id
        projection_id = self.projection_id
        # Create AnnData object out of readable CSV file
        projection_id = secure_filename(projection_id)
        dataset_id = secure_filename(dataset_id)

        projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
        # Sanitize input to prevent path traversal
        projection_adata_path = projection_dir.joinpath("{}.h5ad".format(projection_id))
        projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
        try:
            # READ CSV to make X and var
            df = pd.read_csv(projection_csv_path, sep=',', index_col=0, header=0)
            X = df.to_numpy()
            var = pd.DataFrame(index=df.columns)
            obs = dataset_adata.obs
            obsm = dataset_adata.obsm
            uns = dataset_adata.uns
            # Create the anndata object and write to h5ad
            # Associate with a filename to ensure AnnData is read in "backed" mode
            projection_adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns, filemode='r')
        except Exception as e:
            print(str(e), file=sys.stderr)
            raise ValueError("Could not create projection AnnData object from CSV.")
        # Close dataset adata so that we do not have a stale opened object
        if dataset_adata.isbacked:
            dataset_adata.file.close()

        # For some reason the gene_symbol is not taken in by the constructor
        projection_adata.var["gene_symbol"] = projection_adata.var_names

        # write to projection_adata_path. This ensures that the file is created and up to date with latest projection results
        projection_adata.write(projection_adata_path)
        return projection_adata

    def filter_adata(self):
        # Filter out cells that overlap with the blank space of the image.

        sc.pp.filter_cells(self.dataset_adata, min_genes=self.min_genes)

        # If adata is empty, raise an error
        if self.dataset_adata.n_obs == 0:
            raise ValueError(f"No cells found with at least {self.min_genes} genes. Choose a different gene or lower the filter.")

        sc.pp.normalize_total(self.dataset_adata, inplace=True)
        sc.pp.log1p(self.dataset_adata)

        self.dataset_adata.var_names_make_unique()

    def add_umap(self):
        # Add UMAP information to the adata object. This is a slow process, so we want to show other plots while this is processing

        def create_umap():
            # We need to use the original dataset for UMAP clustering instead of the projection one
            # However we only want to use the cells that have clusters
            adata = self.dataset_adata
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            return adata

        adata_subset_cache_label = f"{self.dataset_id}_{self.min_genes}_adata"

        try:
            # Load the subset Anndata object from cache or create it if it does not exist, with a 1-week time-to-live
            adata = pn.state.as_cached(adata_subset_cache_label, create_umap, ttl=CACHE_EXPIRATION)

            X, Y = (0, 1)
            self.df["UMAP1"] = adata.obsm["X_umap"].transpose()[X].tolist()
            self.df["UMAP2"] = adata.obsm["X_umap"].transpose()[Y].tolist()

            self.umap_fig = self.normal_fig_obj.make_umap_plots()
        except ValueError as e:
            print("Error creating UMAP:", str(e), file=sys.stderr)
            layout = {
                "annotations": [
                    {
                        "text": "Something went wrong with UMAP clustering.",
                        "font": {"size": 20, "color": "red"},
                        "showarrow": False,
                        "x": 0.5,
                        "y": 1.3,
                        "xref": "paper",
                        "yref": "paper"
                    }
                ]
            }
            self.umap_fig = {"data": [], "layout": layout}  # reset the umap figure

        self.umap_pane.object = self.umap_fig

        # Add X_umap to self.adata
        self.adata.obs = adata.obs
        self.adata.obsm["X_umap"] = adata.obsm["X_umap"]

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
        df["spatial1"] = adata.obsm["spatial"].transpose()[X].tolist()
        df["spatial2"] = adata.obsm["spatial"].transpose()[Y].tolist()

        # Add cluster info
        if "clusters" not in selected.obs:
            raise ValueError("No cluster information found in adata.obs")

        df["clusters"] = selected.obs["clusters"].astype("category")

        # Drop any NA clusters
        df = df.dropna(subset=["clusters"])
        self.df = df

    def refresh_dataframe(self, min_genes):
        """
        Refresh the dataframe based on the selected gene and min_genes value.
        Updates Plotly dicts in the Panel app in-place
        """
        self.min_genes = min_genes

        # If the min_genes value has changed... sync to URL
        if self.min_genes != self.settings.min_genes:
            self.settings.min_genes = self.min_genes

        self.filter_adata()

        # self.adata should have the same subset as self.dataset_adata
        self.adata = self.adata[self.dataset_adata.obs.index]

        self.create_gene_df()   # creating the Dataframe is generally fast
        self.map_colors()
        # drop indexes from self.adata not in self.df (since clustering may have removed some cells)
        self.dataset_adata = self.dataset_adata[self.df.index]
        self.adata = self.adata[self.df.index]

        # destroy the old figure objects (to free up memory)
        self.normal_fig_obj = None
        self.zoom_fig_obj = None

        # CET_L19 is "WhBuPrRd"
        # CET_L18 is "YlOrRd"
        #normal_color = cc.cm.CET_L19
        #zoom_color = cc.cm.CET_L18

        self.normal_fig_obj = SpatialNormalSubplot(self.settings, self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, self.norm_gene_symbol, "YlGn", dragmode="select")
        self.zoom_fig_obj = SpatialZoomSubplot(self.settings, self.df, self.spatial_img, self.color_map, self.norm_gene_symbol, "Local", "YlOrRd", dragmode=False)

        self.normal_fig = self.normal_fig_obj.refresh_spatial_fig()

        # The pn.bind function for the zoom callback will not trigger when the normal_fig is refreshed.
        self.zoom_fig = self.zoom_fig_obj.refresh_spatial_fig()

        # Add annotation to UMAP plot to indicate that it is processing
        layout = {
            "annotations": [
                {
                    "text": "Performing UMAP clustering. This may take a minute...",
                    "font": {"size": 20},
                    "showarrow": False,
                    "x": 0.5,
                    "y": 1.3,
                    "xref": "paper",
                    "yref": "paper"
                }
            ]
        }

        self.umap_fig = {"data": [], "layout": layout}  # reset the umap figure
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
        # sort unique clusters by number if numerical, otherwise by name
        try:
            unique_clusters = sorted(unique_clusters, key=lambda x: int(x))
        except:
            unique_clusters = sorted(unique_clusters, key=lambda x: str(x))

        self.unique_clusters = unique_clusters

        if "colors" in df:
            self.color_map = {cluster: df[df["clusters"] == cluster]["colors"].values[0] for cluster in self.unique_clusters}
        else:
            # Some glasbey_bw_colors may not show well on a dark background so use "light" colors if images are not present
            swatch_color = cc.glasbey_bw if self.has_images else cc.glasbey_light

            self.color_map = {cluster: swatch_color[i % len(swatch_color)] for i, cluster in enumerate(self.unique_clusters)}
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
    # Create the app
    sp_panel = SpatialPanel(settings)
    sp_panel.servable(title="Spatial Data Viewer", location=True)
