import gc
import logging
import sys
import traceback
import warnings
from pathlib import Path

import anndata
import colorcet as cc
import numpy as np
import pandas as pd
import panel as pn
import param
import plotly.graph_objects as go
import plotly.io as pio
import scanpy as sc
import spatialdata as sd
from common import (
    Settings,
    SpatialFigure,
    clip_expression_values,
    normalize_searched_gene,
    sort_clusters,
)
from plotly.subplots import make_subplots
from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[2]
www_path = gear_root.joinpath("www")
PROJECTIONS_BASE_DIR = www_path.joinpath("projections")

spatial_path = gear_root.joinpath("www/datasets/spatial")

pio.templates.default = "simple_white"  # no gridlines, white background

SECS_IN_DAY = 86400
CACHE_EXPIRATION = SECS_IN_DAY * 7  # 7 days

# Ignore warnings about plotly GUI events, which propagate to the browser console
warnings.filterwarnings("ignore", "plotly.*unrecognized gui edit.*")

pn.extension("plotly", loading_indicator=True, defer_load=True, nthreads=4) # type: ignore

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]
zoomButtonsToRemove = buttonsToRemove + ["select2d"]

class ExpandedSettings(Settings):
    save = param.Boolean(
        doc="If true, save this configuration as a new display.", default=False
    )
    display_name = param.String(
        doc="Display name for the saved configuration", allow_None=True
    )
    make_default = param.Boolean(
        doc="If true, make this the default display.", default=False
    )

    nosave = param.Boolean(
        doc="If true, do not show the contents related to saving.", default=False
    )


class SpatialNormalSubplot(SpatialFigure):
    """
    Class for creating a spatial plot with a gene expression heatmap, cluster markers, and an image.

    Args:
        settings (ExpandedSettings): Configuration and state for the figure.
        df (pd.DataFrame): DataFrame containing spatial coordinates and expression/cluster data.
        spatial_img (np.ndarray | None): Optional background spatial image.
        color_map (dict): Mapping of cluster/expression values to colors.
        cluster_map (dict): Mapping of cluster IDs to cluster names or metadata.
        gene_symbol (str): Gene symbol for expression visualization.
        expression_name (str): Name of the expression metric to display.
        expression_color (str): Color to use for expression visualization.
        **params: Additional keyword arguments for customization.

    Methods:
        refresh_spatial_fig():
            Regenerates the figure and applies current selections, returning the figure as a dictionary.
        setup_cluster_plot(...):
            Sets up a cluster plot on the given Plotly figure using the provided DataFrame and plot type.
        make_umap_plots():
            Creates UMAP plots for expression and clusters.
        make_violin_plot():
            Creates a violin plot for expression by cluster.
        mirror_selection_callback(event):
            Mirrors a selection event across all plots and refreshes the figure.
    """

    def __init__(
        self,
        settings: ExpandedSettings,
        df: pd.DataFrame,
        spatial_img: np.ndarray | None,
        color_map: dict,
        cluster_map: dict,
        gene_symbol: str,
        expression_name: str,
        expression_color: str,
        **params,
    ):
        super().__init__(
            settings,
            df,
            spatial_img,
            color_map,
            cluster_map,
            gene_symbol,
            expression_name,
            expression_color,
            **params,
        )

        self.selections_dict: dict[str, float] = {}

    def refresh_spatial_fig(self) -> dict:
        """
        Refreshes the spatial figure by regenerating it and applying current selections.

        Returns:
            dict: The updated figure represented as a dictionary.
        """

        self.fig = self.make_fig()
        return self.fig.to_dict()

    def make_fig(self) -> go.Figure:
        self._setup_fig()
        if self.update_selection_ranges():
            self.mirror_selection()
        return self._add_traces_to_fig()

    def setup_cluster_plot(
        self,
        fig: go.Figure,
        dataframe: pd.DataFrame,
        x_col: str,
        y_col: str,
        plot_type: str = "scatter",
        col: int = 1,
        row: int = 1,
        showlegend: bool = True
    ) -> go.Figure:
        """
        Sets up a cluster plot on the given Plotly figure using the provided DataFrame and plot type.

        Parameters:
            fig (go.Figure): The Plotly figure to which the cluster plot will be added.
            dataframe (pd.DataFrame): The data containing cluster assignments and plotting values.
            x_col (str): The name of the column to use for the x-axis.
            y_col (str): The name of the column to use for the y-axis.
            plot_type (str, optional): The type of plot to generate ('scatter' or 'violin'). Defaults to "scatter".
            col (int, optional): The subplot column position. Defaults to 1.
            row (int, optional): The subplot row position. Defaults to 1.
            showlegend (bool, optional): Whether to display the legend. Defaults to True.

        Returns:
            go.Figure: The updated Plotly figure with the cluster plot added.

        Notes:
            - The function expects a 'clusters' column in the DataFrame.
            - Additional plot types can be added as needed.
            - The function updates layout settings such as legend font size and margins.
        """

        # Add traces
        if plot_type == "scatter":
            self.add_cluster_scatter(fig, dataframe, self.color_map, x_col, y_col, col, row)
        elif plot_type == "violin":
            # Ensure clusters are categorical and sorted
            dataframe["clusters"] = dataframe["clusters"].astype("category")
            unique_clusters = dataframe["clusters"].unique()
            sorted_clusters = sort_clusters(unique_clusters)

            for cluster in sorted_clusters:
                cluster_data = dataframe[dataframe["clusters"] == cluster]
                fig.add_trace(
                    go.Violin(
                        x=cluster_data[x_col],
                        y=cluster_data[y_col],
                        name=str(cluster),
                        fillcolor=self.color_map[cluster],
                        line_color=self.color_map[cluster],
                        marker=dict(color="#000000", size=1),
                        showlegend=showlegend,
                    ),
                )
        # Add more plot types as needed

        # Update axes and legend
        font_size = self.get_legend_font_size()
        fig.update_layout(
            legend=dict(font=dict(size=font_size), itemsizing="constant"),
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=False,
            selectdirection="d",
        )
        return fig

    def make_umap_plots(self) -> dict:
        fig = make_subplots(
            rows=1,
            cols=2,
            column_titles=(f"{self.expression_name} UMAP", "clusters UMAP"),
            horizontal_spacing=0.1,
        )

        fig.update_layout(
            xaxis=dict(domain=[0, 0.45]),
            xaxis2=dict(domain=[0.55, 1]),  # Leave room for cluster annotations
        )

        agg = self.create_datashader_agg(x="UMAP2", y="UMAP1")

        ### Expression UMAP

        expression_agg_df = agg["expression"].to_dataframe(name="raw_value")
        # Drop missing values
        expression_agg_df = expression_agg_df.dropna()
        expression_df = expression_agg_df.reset_index()

        # NOTE: I attempted to use tf.Shade for an image (as I've seen in various examples) but the image was not appearing in the plot

        max_x1 = fig.layout.xaxis.domain[1] # type: ignore

        fig.add_scatter(
            col=1,
            row=1,
            x=expression_df["UMAP1"],
            y=expression_df["UMAP2"],
            mode="markers",
            marker=dict(
                color=expression_df["raw_value"],
                colorscale="cividis_r",
                size=self.marker_size,
                colorbar=dict(
                    len=1,  # Adjust the length of the colorbar
                    thickness=15,  # Adjust the thickness of the colorbar (default is 30)
                    x=max_x1,  # Adjust the x position of the colorbar
                ),
            ),
            showlegend=False,
            hovertemplate="Expression: %{marker.color:.2f}<extra></extra>",
        )

        ### Clusters UMAP

        clusters_agg_df = agg["clusters"].to_dataframe(name="clusters_cat_codes")
        # Drop rows where "clusters" is False
        clusters_agg_df = clusters_agg_df[
            clusters_agg_df["clusters_cat_codes"]
        ]
        # The columns we want are in the multi-index, so we need to make them into a dataframe
        clusters_df = clusters_agg_df.index.to_frame(index=False)

        # map cluster cat_codes to cluster names
        clusters_df["clusters"] = (
            clusters_df["clusters_cat_codes"].map(self.cluster_map)
        )
        try:
            # If clusters are float, convert to int (for later sorting)
            clusters_df["clusters"] = clusters_df["clusters"].astype(int).astype("category")
        except ValueError:
            # If clusters are not numeric, keep as category
            clusters_df["clusters"] = clusters_df["clusters"].astype("category")

        self.setup_cluster_plot(fig, clusters_df, "UMAP1", "UMAP2", plot_type="scatter", col=2, row=1)

        fig.update_xaxes(
            showgrid=False, showticklabels=False, ticks="", title_text="UMAP1"
        )
        fig.update_yaxes(
            showgrid=False,
            showticklabels=False,
            ticks="",
            title_text="UMAP2",
            title_standoff=0,
        )

        font_size = self.get_legend_font_size()

        # make legend markers bigger
        fig.update_layout(
            legend=dict(
                font=dict(size=font_size),
                itemclick="toggleothers",
                itemsizing="constant",
            )
        )

        fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=False,
            selectdirection="d",
        )

        return fig.to_dict()

    def make_violin_plot(self) -> dict:
        """
        Creates a violin plot for expression by cluster.

        Returns:
            dict: The violin plot figure represented as a dictionary.
        """
        dataframe = self.df.sort_values(by="raw_value")

        # Make "clusters" a category so that the violin plot will sort the clusters in order
        dataframe["clusters"] = dataframe["clusters"].astype("category")

        fig = go.Figure()

        self.setup_cluster_plot(fig, dataframe, "clusters", "raw_value", plot_type="violin")

        # xaxis needs to be categorical, even for numerical values
        fig.update_xaxes(type="category", title_text="Clusters")
        fig.update_yaxes(title_text="Expression", rangemode="tozero")

        font_size = self.get_legend_font_size()

        fig.update_layout(legend=dict(font=dict(size=font_size), itemsizing="constant"))

        fig.update_layout(
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=False,
            selectdirection="d",
        )
        return fig.to_dict()

    def mirror_selection_callback(self, event: dict) -> dict:
        """
        For a selection event, mirror the selection across all plots.

        Args:
            event (dict): The selection event containing range information for the plot axes.

        Returns:
            dict: The updated spatial figure after applying the selection.
        """

        if event and "range" in event:
            # determine if first or second plot
            x = (
                "x"
                if "x" in event["range"]
                else "x2"
                if "x2" in event["range"]
                else "x3"
            )
            y = (
                "y"
                if "y" in event["range"]
                else "y2"
                if "y2" in event["range"]
                else "y3"
            )

            range_x1 = event["range"][x][0]
            range_x2 = event["range"][x][1]
            range_y1 = event["range"][y][0]
            range_y2 = event["range"][y][1]

            self.selections_dict = dict(
                x0=range_x1, x1=range_x2, y0=range_y1, y1=range_y2
            )

        return self.refresh_spatial_fig()

    def mirror_selection(self):
        """
        Mirrors the current selection across all rows and columns in the figure.

        This method adds the selections stored in `self.selections_dict` to the figure (`self.fig`)
        so that the selection is applied to every subplot (all rows and columns).

        Returns:
            None
        """

        linecolor = "black"
        if self.spatial_img is None or self.platform == "xenium":
            # If there is no spatial image, use white line color for visibility
            linecolor = "white"

        self.fig.add_selection(self.selections_dict,
                            line=dict(
                                color=linecolor,
                                width=3,
                                dash="solid"
                            ),
                            row="all", col="all")

class SpatialZoomSubplot(SpatialFigure):
    """
    Class for creating a zoomed-in spatial plot with a gene expression heatmap, cluster markers, and an image.

    Args:
        settings (ExpandedSettings): Configuration and state for the figure.
        df (pd.DataFrame): DataFrame containing spatial coordinates and expression/cluster data.
        spatial_img (np.ndarray | None): Optional background spatial image.
        color_map (dict): Mapping of cluster/expression values to colors.
        cluster_map (dict): Mapping of cluster IDs to cluster names or metadata.
        gene_symbol (str): Gene symbol for expression visualization.
        expression_name (str): Name of the expression metric to display.
        expression_color (str): Color to use for expression visualization.
        **params: Additional keyword arguments for customization.

    Methods:
        calculate_marker_size():
            Calculates and sets the marker size for the zoomed-in plot.
        refresh_spatial_fig():
            Regenerates the zoomed-in figure and applies current selections, returning the figure as a dictionary.
        make_zoom_fig_callback(event):
            Handles zoom/selection events, filters data accordingly, and refreshes the figure.
    """

    def __init__(
        self,
        settings: ExpandedSettings,
        df: pd.DataFrame,
        spatial_img: np.ndarray | None,
        color_map: dict,
        cluster_map: dict,
        gene_symbol: str,
        expression_name: str,
        expression_color: str,
        **params,
    ):
        super().__init__(
            settings,
            df,
            spatial_img,
            color_map,
            cluster_map,
            gene_symbol,
            expression_name,
            expression_color,
            **params,
        )

        # Preserve the original dataframe for filtering
        self.orig_df = df

        # Set ranges to be the selection.
        # This will fix image cropping issues on the intiial render
        if (
            self.settings.selection_x1 != self.settings.selection_x2
            or self.settings.selection_y1 != self.settings.selection_y2
        ):
            self.range_x1 = self.settings.selection_x1
            self.range_x2 = self.settings.selection_x2
            self.range_y1 = self.settings.selection_y1
            self.range_y2 = self.settings.selection_y2

    def calculate_marker_size(self) -> None:
        """
        Calculates and sets the marker size for the zoomed-in plot.

        Returns:
            None
        """
        # Calculate the range of the selection
        x_range = self.settings.selection_x2 - self.settings.selection_x1 # type: ignore
        y_range = self.settings.selection_y2 - self.settings.selection_y1 # type: ignore

        # Calculate the marker size based on the range of the selection
        # The marker size will scale larger as the range of the selection gets more precise
        self.marker_size = int(1 + 2500 / (x_range + y_range))

        if self.platform == "visium":
            self.marker_size += 3

    def refresh_spatial_fig(self) -> dict:
        """
        Refreshes the zoomed-in spatial figure by regenerating it and applying current selections.

        Returns:
            dict: The updated figure represented as a dictionary.
        """

        self.fig = self.make_fig()
        return self.fig.to_dict()

    def make_fig(self) -> go.Figure:
        self._setup_fig()
        if self.update_selection_ranges():
            # Ensure the axes are set to the selection ranges
            self.setup_zoom_fig_params()
        return self._add_traces_to_fig()

    def make_zoom_fig_callback(self, event: dict) -> dict:
        """
        Handles zoom/selection events, filters data accordingly, and refreshes the figure.

        Args:
            event (dict): The selection event containing range information for the plot axes.

        Returns:
            dict: The updated spatial figure after applying the zoom selection.
        """
        if event and "range" in event:
            # determine if first or second plot
            x = (
                "x"
                if "x" in event["range"]
                else "x2"
                if "x2" in event["range"]
                else "x3"
            )
            y = (
                "y"
                if "y" in event["range"]
                else "y2"
                if "y2" in event["range"]
                else "y3"
            )

            self.range_x1 = event["range"][x][0]
            self.range_x2 = event["range"][x][1]
            self.range_y1 = event["range"][y][0]
            self.range_y2 = event["range"][y][1]


            self.settings.selection_x1 = self.range_x1
            self.settings.selection_x2 = self.range_x2
            self.settings.selection_y1 = self.range_y1
            self.settings.selection_y2 = self.range_y2

        return self.refresh_spatial_fig()

    def setup_zoom_fig_params(self):
        """
        Configures the figure parameters for zooming into a selected region.

        This method performs the following actions:
        - Adjusts the marker size for better visibility when a selection is being viewed.
        - Updates the x-axis range of the figure to match the selected region, using `selection_x1` and `selection_x2` from the settings.
        - Updates the y-axis range of the figure to match the selected region, flipping the axis so that (0,0) is at the top-left corner, using `selection_x2` and `selection_y1` from the settings.
        - Applies these axis updates to all subplot columns

        """

        if not self.fig:
            raise ValueError("Figure is not initialized. Cannot set zoom parameters.")

        if not self.settings.selection_x1 or not self.settings.selection_x2: # type: ignore
            raise ValueError("Selection coordinates are not set. Cannot set zoom parameters.")

        # Preserve the original dataframe for filtering
        dataframe = self.orig_df

        # Filter the data based on the selected range
        self.df = dataframe[
            (dataframe["spatial1"] >= self.settings.selection_x1)
            & (dataframe["spatial1"] <= self.settings.selection_x2)
            & (dataframe["spatial2"] >= self.settings.selection_y1)
            & (dataframe["spatial2"] <= self.settings.selection_y2)
        ]

        self.calculate_marker_size()

        self.fig.update_xaxes(
            range=[self.settings.selection_x1, self.settings.selection_x2], row="all", col="all"
        )
        # y-axis needs to be flipped. 0,0 is the top left corner
        self.fig.update_yaxes(
            range=[self.settings.selection_y2, self.settings.selection_y1], row="all", col="all"
        )


class SpatialPanel(pn.viewable.Viewer):
    """
    Panel view class for displaying spatial transcriptomics visualizations in a Panel dashboard.
    """

    dataset_id: str
    gene_symbol: str
    min_genes: int
    projection_id: str | None = None
    selection_x1: float | None = None
    selection_x2: float | None = None
    selection_y1: float | None = None
    selection_y2: float | None = None

    def __init__(self, settings: ExpandedSettings, **params):
        super().__init__(**params)

        self.settings = settings

        if pn.state.location is not None:
            pn.state.location.sync(
                self.settings,
                {
                    "dataset_id": "dataset_id",
                    "gene_symbol": "gene_symbol",
                    "min_genes": "min_genes",
                    "projection_id": "projection_id",
                    "selection_x1": "selection_x1",
                    "selection_x2": "selection_x2",
                    "selection_y1": "selection_y1",
                    "selection_y2": "selection_y2",
                    "display_height": "height",
                    "display_width": "width",
                    "expression_min_clip": "expression_min_clip",
                    "save": "save",
                    "display_name": "display_name",
                    "make_default": "make_default",
                    "nosave": "nosave",
                },
            )

        self.dataset_id = self.settings.dataset_id # type: ignore
        self.gene_symbol = self.settings.gene_symbol # type: ignore
        self.min_genes = self.settings.min_genes # type: ignore
        self.projection_id = self.settings.projection_id # type: ignore

        self.nosave: bool = self.settings.nosave # type: ignore

        self.platform = None # Will be set in prep_sdata()

        layout_height = 1520
        if self.settings.display_height and self.settings.display_height > 0: # type: ignore
            layout_height: int = self.settings.display_height # type: ignore

        layout_width = 1100  # Default width
        if self.settings.display_width and self.settings.display_width > 0: # type: ignore
            layout_width: int = self.settings.display_width # type: ignore

        # When width is too short, things break down.
        if layout_width < 1100:
            layout_width = 1100
            layout_height = 1520

        self.normal_fig = dict(data=[], layout={})
        self.zoom_fig = dict(data=[], layout={})
        self.umap_fig = dict(data=[], layout={})
        self.violin_fig = dict(data=[], layout={})

        # the "figs" are dicts which allow for patching the figure objects
        # See -> https://panel.holoviz.org/reference/panes/Plotly.html#patching
        # You can also replace the Figure directly, but I had occasional issues with returning a figure ih the "bind" function

        self.normal_pane = pn.pane.Plotly(
            self.normal_fig,
            config={
                "doubleClick": "reset",
                "displayModeBar": True,
                "modeBarButtonsToRemove": buttonsToRemove,
            },
            height=350,
            sizing_mode="stretch_width",
        )

        self.zoom_pane = pn.pane.Plotly(
            self.zoom_fig,
            config={
                "displayModeBar": True,
                "modeBarButtonsToRemove": zoomButtonsToRemove,
                "doubleClick": None,
            },
            height=350,
            sizing_mode="stretch_width",
        )

        self.umap_pane = pn.pane.Plotly(
            self.umap_fig,
            config={
                "displayModeBar": True,
                "modeBarButtonsToRemove": zoomButtonsToRemove,
            },
            height=350,
            sizing_mode="stretch_width",
        )

        self.violin_pane = pn.pane.Plotly(
            self.violin_fig,
            config={
                "displayModeBar": True,
                "modeBarButtonsToRemove": zoomButtonsToRemove,
            },
            height=350,
            sizing_mode="stretch_width",
        )

        self.layout = pn.Column(pn.bind(self.init_data), width=layout_width)

       # ? This can be useful for filtering datasets even for projections, but how best to word it?
        min_slider_width = 300
        self.min_genes_slider = pn.widgets.IntSlider(
            name="Filter - Mininum genes per observation",
            start=0,
            end=500,
            step=25,
            width=min_slider_width,
            value=self.min_genes,
        )

        self.display_name = pn.widgets.TextInput(
            name="Display name",
            placeholder="Name this display to save...",
            width=250,
            visible=not self.nosave,
        )
        self.save_button = pn.widgets.Button(
            name="Save settings", button_type="primary", width=100, align="end"
            , visible=not self.nosave
        )
        self.make_default = pn.widgets.Checkbox(
            name="Make this the default display", value=False
            , visible=not self.nosave
        )

        markdown_width = 675
        spacer_width = markdown_width - min_slider_width    # Make default button should left-align with the above text input

        if spacer_width < 0:
            spacer_width = 100

        self.pre_layout = pn.Column(
            pn.Row(
                pn.pane.Markdown(
                    "## Select a region to modify zoomed in view in the bottom panel",
                    height=30,
                    width=markdown_width,
                ),
                self.display_name,
                self.save_button,
            ),
            pn.Row(self.min_genes_slider, pn.Spacer(width=spacer_width), self.make_default),
        )

        self.fig_layout = pn.Column(
            pn.layout.Divider(height=5),  # default margins
            pn.pane.Markdown("## Zoomed in view", height=30),
            self.zoom_pane,
            pn.layout.Divider(height=5),  # default margins
            pn.pane.Markdown("## UMAP plots", height=30),
            self.umap_pane,
            pn.layout.Divider(height=5),  # default margins
            pn.pane.Markdown("## Violin plot", height=30),
            self.violin_pane,
        )

        # A condensed layout should only have the normal pane.
        self.plot_layout = pn.Column(
            self.pre_layout,
            self.normal_pane,
            self.fig_layout,
            width=layout_width,
            height=layout_height,
        )

        # SAdkins - Have not quite figured out when to use "watch" but I think it mostly applies when a callback does not return a value
        def refresh_dataframe_callback(value):
            self.dataset_adata = self.dataset_adata_orig.copy()
            self.adata = self.adata_orig.copy()

            with (
                self.normal_pane.param.update(loading=True),
                self.zoom_pane.param.update(loading=True),
                self.umap_pane.param.update(loading=True),
                self.violin_pane.param.update(loading=True),
                self.min_genes_slider.param.update(disabled=True),
            ):
                self.refresh_dataframe(value)

            with (
                self.umap_pane.param.update(loading=True),
                self.min_genes_slider.param.update(disabled=True),
            ):
                self.add_umap()

            gc.collect()  # Collect garbage to free up memory


        pn.bind(
            refresh_dataframe_callback,
            self.min_genes_slider.param.value_throttled,
            watch=True,
        )

        def selection_callback(event):
            if self.zoom_fig_obj is not None:
                self.zoom_pane.object = self.zoom_fig_obj.make_zoom_fig_callback(event)
            if self.normal_fig_obj is not None:
                # If the normal pane is selected, mirror the selection across all plots
                self.normal_pane.object = self.normal_fig_obj.mirror_selection_callback(
                    event
                )
            gc.collect()  # Collect garbage to free up memory

        pn.bind(
            selection_callback, self.normal_pane.param.selected_data, watch=True
        )

        def save_settings_callback(event):
            self.settings.save = True
            self.settings.display_name = self.display_name.value
            self.settings.make_default = self.make_default.value
            # Make button "is-loading" while saving
            self.save_button.button_type = "success"
            self.save_button.name = "Saving..."
            self.save_button.disabled = True
            logging.error(self.settings)

        def reset_save_button_callback(event):
            if not event.name == "save":
                return
            self.settings.save = False
            self.save_button.button_type = "primary"
            self.save_button.name = "Save settings"
            self.save_button.disabled = False

        pn.bind(save_settings_callback, self.save_button, watch=True)
        self.settings.param.watch(
            reset_save_button_callback, "save", onlychanged=True
        )

    def __panel__(self):
        # This is run when the app is loaded. Return the final layout of the app
        return self.layout

    def loading_indicator(self, label):
        return pn.indicators.LoadingSpinner(
            value=True, name=label, align="center", color="info"
        )

    def init_data(self):
        try:
            yield self.loading_indicator("Processing data file...")

            def create_adata():
                try:
                    self.prep_sdata()
                    adata = self.prep_adata()
                except ValueError:
                    raise
                return adata

            # Load the Anndata object (+ image name) from cache or create it if it does not exist, with a 1-week time-to-live
            try:
                adata = create_adata()
            except ValueError as e:
                yield pn.pane.Alert(f"Error: {e}", alert_type="danger")
                raise

            self.dataset_adata_orig = adata  # Original dataset adata
            self.dataset_adata = (
                self.dataset_adata_orig.copy()
            )  # Copy we manipulate (filtering, etc.)
            self.adata = self.dataset_adata.copy()  # adata object to use for plotting

            # Modify the adata object to use the projection ID if it exists
            if self.projection_id:
                self.adata = self.create_projection_adata()

            if self.settings.expression_min_clip is not None:
                self.adata = clip_expression_values(self.adata, min_clip=self.settings.expression_min_clip)

            self.adata_orig = (
                self.adata.copy()
            )  # This is to restore when the min_genes slider is changed

            yield self.loading_indicator(
                "Processing data to create plots. This may take a minute..."
            )

            yield self.refresh_dataframe(self.min_genes)

            with (
                self.umap_pane.param.update(loading=True),
                self.min_genes_slider.param.update(disabled=True),
            ):
                self.add_umap()
        except Exception as e:
            logging.error(traceback.format_exc())
            yield pn.pane.Alert(
                "An unexpected error occurred. Please try again or contact support.",
                alert_type="danger",
            )
            raise

    #@pn.io.profile(name="prep_sdata_expanded", engine="pyinstrument")
    def prep_sdata(self):
        zarr_path = spatial_path / f"{self.dataset_id}.zarr"
        if not zarr_path.exists():
            raise ValueError(f"Dataset {self.dataset_id} not found")

        sdata = sd.read_zarr(zarr_path)

        try:
            self.platform = sdata.tables["table"].uns["platform"]
        except KeyError:
            raise ValueError("No platform information found in the dataset")

        lib_path = gear_root.joinpath("lib")
        sys.path.append(str(lib_path))

        from gear.spatialhandler import SPATIALTYPE2CLASS

        # Ensure the spatial data type is supported
        if self.platform not in SPATIALTYPE2CLASS.keys():
            logging.error("Invalid or unsupported spatial data type")
            logging.error("Supported types: {0}".format(SPATIALTYPE2CLASS.keys()))
            raise ValueError(
                f"Invalid or unsupported spatial data type: {self.platform}"
            )

        # Use uploader class to determine correct helper functions
        self.spatial_obj = SPATIALTYPE2CLASS[self.platform]()
        self.spatial_obj.sdata = sdata
        self.spatial_obj.subset_sdata()

        # Dictates if this will be a 2- or 3-plot figure
        self.has_images = self.spatial_obj.has_images

        # Standardize the spatial data object
        obs = self.spatial_obj.sdata.tables["table"].obs
        if not ("spatial1" in obs.columns and "spatial2" in obs.columns):
            self.spatial_obj.subset_sdata()
            self.spatial_obj.scale_and_translate_sdata()

            # The SpatialData object table should have coordinates, but they are not translated into the image space
            # Each observation has an associated polygon "shape" in the image space, and we can get the centroid of that shape
            self.spatial_obj.merge_centroids_with_obs()

        # Convert image to dataframe if it exists
        self.spatial_img = None
        try:
            self.spatial_img = self.spatial_obj.extract_img()
        except Exception as e:
            # This is not a fatal error. Some datasets do not have images
            logging.info(f"No image found or error converting image to dataframe: {e}")

    def prep_adata(self):
        # Create AnnData object
        # Need to include image since the bounding box query does not filter the image data by coordinates
        # Each Image is downscaled (or upscaled) during rendering to fit a 2000x2000 pixels image (downscaled_hires)
        try:
            self.spatial_obj.adata = self.spatial_obj.sdata.tables["table"]
        except Exception as e:
            logging.error(f"Error converting sdata to adata: {str(e)}")
            raise ValueError(
                "Could not convert sdata to adata. Check the spatial data type and bounding box coordinates."
            )
        return self.spatial_obj.adata

    def create_projection_adata(self):
        dataset_adata = self.adata
        dataset_id = self.dataset_id
        projection_id = self.projection_id
        # Create AnnData object out of readable CSV file
        if projection_id:
            projection_id = secure_filename(projection_id)
        dataset_id = secure_filename(dataset_id)

        projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
        # Sanitize input to prevent path traversal
        projection_adata_path = projection_dir.joinpath("{}.h5ad".format(projection_id))
        projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
        try:
            # READ CSV to make X and var
            dataframe = pd.read_csv(projection_csv_path, sep=",", index_col=0, header=0)
            X = dataframe.to_numpy()
            var = pd.DataFrame(index=dataframe.columns)
            obs = dataset_adata.obs
            obsm: np.ndarray = dataset_adata.obsm
            uns = dataset_adata.uns
            # Create the anndata object and write to h5ad
            # Associate with a filename to ensure AnnData is read in "backed" mode
            projection_adata = anndata.AnnData(
                X=X, obs=obs, var=var, obsm=obsm, uns=uns, filemode="r"
            )
        except Exception as e:
            logging.error(traceback.format_exc())
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
            raise ValueError(
                f"No cells found with at least {self.min_genes} genes. Choose a different gene or lower the filter."
            )

        sc.pp.normalize_total(self.dataset_adata, inplace=True)
        sc.pp.log1p(self.dataset_adata)

        self.dataset_adata.var_names_make_unique()

    #@pn.io.profile(name="add_umap_expanded", engine="pyinstrument")
    def add_umap(self):
        # Add UMAP information to the adata object. This is a slow process, so we want to show other plots while this is processing
        def create_umap():
            # We need to use the original dataset for UMAP clustering instead of the projection one
            # However we only want to use the cells that have clusters
            adata = self.dataset_adata

            VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'],
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]

            # If a pair in VALID_UMAP_PAIRS are already computed, convert to X_umap and return
            cols = adata.obs.columns.tolist()
            for umap_pair in VALID_UMAP_PAIRS:
                if all(key in cols for key in umap_pair):
                    adata.obsm["X_umap"] = adata.obs[[umap_pair[0], umap_pair[1]]].to_numpy()
                    return adata

            # If UMAP is already computed, return early
            if "X_umap" in adata.obsm.keys():
                return adata

            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            return adata

        try:
            # Load the subset AnnData object from cache or create it if it does not exist, with a 1-week time-to-live
            adata = create_umap()

            X, Y = (0, 1)
            self.df["UMAP1"] = adata.obsm["X_umap"].transpose()[X].tolist()
            self.df["UMAP2"] = adata.obsm["X_umap"].transpose()[Y].tolist()

            if self.normal_fig_obj is None:
                raise ValueError(
                    "Normal figure object is not initialized. Cannot create UMAP plots."
                )
            self.umap_fig = self.normal_fig_obj.make_umap_plots()
        except ValueError as e:
            logging.error("Error creating UMAP")
            logging.error(traceback.format_exc())
            layout = {
                "annotations": [
                    {
                        "text": "Something went wrong with UMAP clustering.",
                        "font": {"size": 20, "color": "red"},
                        "showarrow": False,
                        "x": 0.5,
                        "y": 1.3,
                        "xref": "paper",
                        "yref": "paper",
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
        dataset_genes = set(self.adata.var["gene_symbol"].unique())
        norm_gene_symbol = normalize_searched_gene(dataset_genes, self.gene_symbol)
        if norm_gene_symbol is None:
            raise ValueError(
                f"Gene symbol '{self.gene_symbol}' not found in dataset {self.dataset_id}."
            )
        self.norm_gene_symbol = norm_gene_symbol
        gene_filter = adata.var.gene_symbol == self.norm_gene_symbol
        selected = adata[:, gene_filter]
        selected.var.index = pd.Index(["raw_value"])
        dataframe = selected.to_df()

        # Add spatial coords
        dataframe["spatial1"] = selected.obs["spatial1"]
        dataframe["spatial2"] = selected.obs["spatial2"]

        # Add cluster info
        if "clusters" not in selected.obs:
            raise ValueError("No cluster information found in adata.obs")

        dataframe["clusters"] = selected.obs["clusters"].astype("category")
        dataframe["clusters_cat_codes"] = dataframe["clusters"].cat.codes.astype("category")

        self.cluster_map = {
            code: dataframe[dataframe["clusters_cat_codes"] == code]["clusters"].to_numpy()[
                0
            ]
            for code in dataframe["clusters_cat_codes"].unique()
        }

        # Drop any NA clusters
        dataframe = dataframe.dropna(subset=["clusters"])
        self.df = dataframe

    def refresh_dataframe(self, min_genes):
        """
        Refresh the dataframe based on the selected gene and min_genes value.
        Updates Plotly dicts in the Panel app in-place

        Args:
            min_genes (int): The minimum number of genes to filter observations.

        Returns:
            None
        """
        self.min_genes = min_genes

        # If the min_genes value has changed... sync to URL
        if self.min_genes != self.settings.min_genes:
            self.settings.min_genes = self.min_genes

        self.filter_adata()

        # self.adata should have the same subset as self.dataset_adata
        self.adata = self.adata[self.dataset_adata.obs.index]

        self.create_gene_df()  # creating the Dataframe is generally fast
        self.map_colors()
        # drop indexes from self.adata not in self.df (since clustering may have removed some cells)
        self.dataset_adata = self.dataset_adata[self.df.index]
        self.adata = self.adata[self.df.index]

        # destroy the old figure objects (to free up memory)
        self.normal_fig_obj = None
        self.zoom_fig_obj = None

        # CET_L19 is "WhBuPrRd"
        # CET_L18 is "YlOrRd"
        # normal_color = cc.cm.CET_L19
        # zoom_color = cc.cm.CET_L18

        self.normal_fig_obj = SpatialNormalSubplot(
            self.settings,
            self.df,
            self.spatial_img,
            self.color_map,
            self.cluster_map,
            self.norm_gene_symbol,
            self.norm_gene_symbol,
            "YlGn",
            dragmode="select",
            platform=self.platform
        )
        self.zoom_fig_obj = SpatialZoomSubplot(
            self.settings,
            self.df,
            self.spatial_img,
            self.color_map,
            self.cluster_map,
            self.norm_gene_symbol,
            "Local",
            "YlOrRd",
            dragmode=False,
            platform=self.platform
        )

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
                    "yref": "paper",
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
        dataframe = self.df
        # Assuming df is your DataFrame and it has a column "clusters"
        unique_clusters = dataframe["clusters"].unique()
        sorted_clusters = sort_clusters(unique_clusters)

        self.unique_clusters = sorted_clusters

        if "colors" in dataframe:
            self.color_map = {
                cluster: dataframe[dataframe["clusters"] == cluster][
                    "colors"
                ].to_numpy()[0]
                for cluster in self.unique_clusters
            }
        else:
            # Some glasbey_bw_colors may not show well on a dark background so use "light" colors if images are not present
            # Prepending b_ to the name will return a list of RGB colors (though glasbey_light seems to already do this)
            swatch_color = (
                cc.b_glasbey_bw if self.spatial_img is not None else cc.glasbey_light
            )

            self.color_map = {
                cluster: swatch_color[i % len(swatch_color)]
                for i, cluster in enumerate(self.unique_clusters)
            }
            # Map the colors to the clusters
            dataframe["colors"] = dataframe["clusters"].map(self.color_map)
        self.df = dataframe


### MAIN APP ###

# Reset logging level to "error" to suppress a bokeh "dropping patch" info message
# https://github.com/bokeh/bokeh/issues/13229
logging.getLogger().setLevel(logging.ERROR)
settings = ExpandedSettings()

# If not params passed, just show OK as a way to test the app
if pn.state.location is None or not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()

else:
    # Create the app
    sp_panel = SpatialPanel(settings)
    sp_panel.servable(title="Spatial Data Viewer", location=True)
