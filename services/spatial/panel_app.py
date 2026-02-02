import logging
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import panel as pn
import plotly.graph_objects as go
import plotly.io as pio
from common import (
    Settings,
    SpatialFigure,
    clip_expression_values,
)
from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[2]
www_path = gear_root.joinpath("www")
PANEL_CSV_CACHE_DIR = www_path / "cache" / "spatial_panel"
SPATIAL_IMAGE_NAME = "spatial_img.npy"

pio.templates.default = "simple_white"  # no gridlines, white background

SECS_IN_DAY = 86400
CACHE_EXPIRATION = SECS_IN_DAY * 7  # 7 days

# Ignore warnings about plotly GUI events, which propagate to the browser console
warnings.filterwarnings("ignore", "plotly.*unrecognized gui edit.*")

pn.extension("plotly", loading_indicator=True, defer_load=True, nthreads=4)  # type: ignore

# Keep only box select
buttonsToRemove = ["zoom", "pan", "zoomIn", "zoomOut", "autoScale", "lasso2d"]
zoomButtonsToRemove = buttonsToRemove + ["select2d"]


class SpatialCondensedSubplot(SpatialFigure):
    """
    SpatialCondensedSubplot is a specialized subclass of SpatialFigure designed for creating and managing spatial data visualizations with condensed subplots using Plotly.

    This class supports flexible visualization layouts depending on the presence of a spatial image and the use of clusters, enabling side-by-side comparison of raw images, cluster assignments, and gene expression data. It provides interactive selection and zooming capabilities, allowing users to focus on specific spatial regions and dynamically update the visualization based on user input.

    Key Features:
    - Dynamically generates subplots for spatial images, clusters, and gene expression.
    - Handles interactive selection events to update and zoom into regions of interest.
    - Supports both cluster-based and expression-based visualizations.
    - Manages color mapping, legends, and layout adjustments for clear presentation.
    - Integrates with a Settings object to persist and update selection ranges.

        settings (Settings): Configuration object containing selection and visualization parameters.
        df (pd.DataFrame): DataFrame containing spatial and expression data.
        spatial_img (np.ndarray | None): Optional spatial image to display as a subplot.
        color_map (dict): Mapping of cluster or expression values to colors.
        cluster_map (dict): Mapping of cluster identifiers to display names or colors.
        expression_name (str): The name of the expression metric to display.
        expression_color (str, optional): Color map for expression visualization. Defaults to "YlGn".
        zoom_expression_color (str, optional): Color map for zoomed-in expression visualization. Defaults to "YlOrRd".
        use_clusters (bool, optional): Whether to display cluster-based plots. Defaults to False.
        **params: Additional keyword arguments for customization.

    Attributes:
        orig_df (pd.DataFrame): Original DataFrame for filtering and reference.
        final_col (int): Index of the final subplot column.
        selections_dict (dict[str, float]): Stores current selection ranges for interactive updates.
        use_clusters (bool): Indicates if cluster visualization is enabled.
        zoom_expression_color (str): Color map for zoomed-in expression plots.
        annotation_col (int): Index of the subplot column used for annotation.
        marker_size (int): Size of markers in the scatter plots.

    Methods:
        make_fig(): Creates and configures the Plotly figure with appropriate subplots and traces.
        refresh_spatial_fig(): Regenerates the figure and applies current selections.
        calculate_zoom_marker_size(): Adjusts marker size based on zoom selection.
        selection_callback(event): Handles selection events and updates the figure accordingly.
        mirror_selection(): Applies the current selection to all relevant subplots.
        setup_zoom_fig_params(): Configures subplot parameters for zoomed-in views.
    """

    def __init__(
        self,
        settings: Settings,
        df: pd.DataFrame,
        spatial_img: np.ndarray | None,
        color_map: dict,
        cluster_map: dict,
        expression_name: str,
        expression_color: str = "YlGn",
        zoom_expression_color: str = "YlOrRd",
        use_clusters: bool = False,
        **params,
    ):

        super().__init__(
            settings,
            df,
            spatial_img,
            color_map,
            cluster_map,
            expression_name,
            expression_color,
            **params,
        )
        # Preserve the original dataframe for filtering
        self.orig_df = df

        self.final_col = 3 if self.spatial_img is not None else 2

        self.selections_dict: dict[str, float] = {}
        self.use_clusters = use_clusters
        self.zoom_expression_color = zoom_expression_color

        # if all four selection range values are the same, they were not set in the query params, so use the default values
        self.update_selection_ranges()

    def make_fig(self) -> go.Figure:
        """
        Creates and configures a Plotly figure with subplots for spatial data visualization.

        The figure layout and content depend on the presence of a spatial image and whether clusters are used.
        - If a spatial image is present, three subplots are created: "Image Only", "Clusters"/"Expression", and "Zoomed View".
        - If no spatial image is present, two subplots are created: "Clusters"/"Expression" and "Zoomed View".

        The method:
            - Sets up subplot domains and titles.
            - Adds image and/or cluster/expression traces as appropriate.
            - Updates axes and selection ranges.
            - Configures legend, background color, and layout properties.
            - Handles zoomed-in views and annotation columns.

        Returns:
            plotly.graph_objs._figure.Figure: The configured Plotly figure object.
        """

        # domain is adjusted whether there are images or not
        if self.spatial_img is not None:
            titles = (
                "Image Only",
                "Clusters"
                if self.use_clusters
                else f"{self.expression_name} Expression",
                "Zoomed View",
            )

            self.create_subplots(num_cols=3, titles=titles)  # sets self.fig
            self.annotation_col = 2
            if not self.use_clusters:
                self.fig.update_layout(
                    xaxis=dict(domain=[0, 0.29]),
                    xaxis2=dict(domain=[0.33, 0.62]),
                    xaxis3=dict(domain=[0.71, 1]),
                )
        else:
            titles = (
                "Clusters"
                if self.use_clusters
                else f"{self.expression_name} Expression",
                "Zoomed View",
            )

            self.create_subplots(num_cols=2, titles=titles)  # sets self.fig
            self.annotation_col = 1

            if not self.use_clusters:
                self.fig.update_layout(
                    xaxis=dict(domain=[0, 0.45]),
                    xaxis2=dict(domain=[0.55, 1]),
                )
        self.update_axes()

        # Get max domain for axis 1 or axis 2
        if self.annotation_col == 1:
            max_x = self.fig.layout.xaxis.domain[1]  # type: ignore
        else:
            max_x = self.fig.layout.xaxis2.domain[1]  # type: ignore
        cbar_loc = max_x + 0.005  # Leave space for colorbar

        if self.spatial_img is not None:
            self.add_image_trace()

        agg = self.create_datashader_agg(x="spatial2", y="spatial1")

        if self.use_clusters:
            # If using clusters, we need to add the cluster traces first
            # Legend will be placed after the final plot
            self.add_cluster_traces(
                cluster_agg=agg["clusters"],
                cluster_col=self.annotation_col,
                show_legend=True,
            )
        else:
            # If not using clusters, we can add the expression trace first
            self.fig.add_trace(
                self.make_expression_scatter(
                    agg["expression"],
                    cbar_loc=cbar_loc,
                    expression_color=self.expression_color,
                ),
                row=1,
                col=self.annotation_col,
            )

        # if all four selection range values are the same, they were not set in the query params, so use the default values
        if self.update_selection_ranges():
            self.mirror_selection()
            self.setup_zoom_fig_params()

        # Update the x and y axes ranges for the "zoom in" plot
        self.annotation_col += 1

        # Get max domain for axis 2 or axis 3
        if self.annotation_col == 2:
            max_x = self.fig.layout.xaxis2.domain[1]  # type: ignore
        else:
            max_x = self.fig.layout.xaxis3.domain[1]  # type: ignore
        cbar_loc = max_x + 0.005  # Leave space for colorbar

        # Add zoom traces
        if self.use_clusters:
            self.add_cluster_traces(
                cluster_agg=agg["clusters"],
                cluster_col=self.annotation_col,
                show_legend=False,
            )
        else:
            self.fig.add_trace(
                self.make_expression_scatter(
                    agg["expression"],
                    cbar_loc=cbar_loc,
                    expression_color=self.zoom_expression_color,
                ),
                row=1,
                col=self.annotation_col,
            )

        font_size = self.get_legend_font_size()

        # make legend markers bigger
        self.fig.update_layout(
            legend=dict(font=dict(size=font_size), itemsizing="constant"),
            margin=dict(l=20, r=0, t=50, b=10),
            width=None,
            height=None,
            dragmode=self.dragmode,
            selectdirection="d",
        )

        # Mirror the background color of visium images
        # Do not set if image is present as there is a slight padding between data and axes ticks
        plot_bgcolor = None
        if self.spatial_img is None:
            plot_bgcolor = "#000000"

        # adjust domains of all 3 plots, leaving enough space for the colorbar and legend
        self.fig.update_layout(
            plot_bgcolor=plot_bgcolor,
        )

        return self.fig

    def refresh_spatial_fig(self) -> dict:
        """
        Refreshes the spatial figure by regenerating it and applying current selections.

        Returns:
            dict: The updated figure represented as a dictionary.
        """
        # Reset some things
        self.marker_size = 2

        self.fig = self.make_fig()
        return self.fig.to_dict()

    def calculate_zoom_marker_size(self) -> None:
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
        self.marker_size = int(1 + 2500 / (x_range + y_range)) + 3

    def selection_callback(self, event: dict) -> dict:
        """
        Handles selection events from a plot and updates selection ranges.

        This method processes the selection event dictionary, determines which subplot
        the selection applies to (based on the presence of a spatial image), extracts
        the selected x and y ranges, and updates both an internal selections dictionary
        and the settings object with these values. Finally, it refreshes the spatial
        figure to reflect the new selection.

        Args:
            event (dict): A dictionary containing selection event data, expected to have
                a "range" key with x and y range information.

        Returns:
            dict: The updated spatial figure data after applying the selection.
        """

        if event and "range" in event:
            # determine if first or second plot
            if self.spatial_img is not None:
                # If there is a spatial image, we have three subplots
                x = "x" if "x" in event["range"] else "x2"
                y = "y" if "y" in event["range"] else "y2"
            else:
                x = "x"
                y = "y"

            range_x1 = event["range"][x][0]
            range_x2 = event["range"][x][1]
            range_y1 = event["range"][y][0]
            range_y2 = event["range"][y][1]

            self.selections_dict = dict(
                x0=range_x1, x1=range_x2, y0=range_y1, y1=range_y2
            )

            # update the Settings selection values
            self.settings.selection_x1 = range_x1
            self.settings.selection_x2 = range_x2
            self.settings.selection_y1 = range_y1
            self.settings.selection_y2 = range_y2

        # Selection will be mirrored across plots in make_fig()
        return self.refresh_spatial_fig()

    def mirror_selection(self):
        """
        Sets up the selection for the figure by specifying the columns to include and adding the selection to the figure.

        This method creates a list of column indices from 1 up to (but not including) `self.final_col`, and applies the selections defined in `self.selections_dict` to all rows and the specified columns of `self.fig`.

        Returns:
            None
        """

        linecolor = "black"
        #if self.spatial_img is None or self.platform == "xenium":
            # If there is no spatial image, use white line color for visibility
        #    linecolor = "white"

        cols = list(range(1, self.final_col))
        self.fig.add_selection(
            self.selections_dict,
            line=dict(
                color=linecolor,
                width=3,
                dash="solid"
            ),
            row="all",
            col=cols,
        )

    def setup_zoom_fig_params(self):
        """
        Configures the figure parameters for zooming into a selected region.

        This method performs the following actions:
        - Adjusts the marker size for better visibility when a selection is being viewed.
        - Updates the x-axis range of the figure to match the selected region, using `selection_x1` and `selection_x2` from the settings.
        - Updates the y-axis range of the figure to match the selected region, flipping the axis so that (0,0) is at the top-left corner, using `selection_x2` and `selection_y1` from the settings.
        - Applies these axis updates to the subplot column specified by `final_col`.

        Assumes that `self.fig`, `self.settings`, and `self.final_col` are properly initialized.
        """

        if not self.fig:
            raise ValueError("Figure is not initialized. Cannot set zoom parameters.")

        if not self.settings.selection_x1 or not self.settings.selection_x2: # type: ignore
            raise ValueError("Selection coordinates are not set. Cannot set zoom parameters.")

        if not self.final_col:
            raise ValueError("Final column is not set. Cannot set zoom parameters.")

        # Viewing a selection, so increase the marker size
        self.calculate_zoom_marker_size()

        self.fig.update_xaxes(
            range=[self.settings.selection_x1, self.settings.selection_x2],
            col=self.final_col,
        )
        # y-axis needs to be flipped. 0,0 is the top left corner
        self.fig.update_yaxes(
            range=[self.settings.selection_y2, self.settings.selection_y1],
            col=self.final_col,
        )


class SpatialPanel(pn.viewable.Viewer):
    """
    Panel view class for displaying spatial transcriptomics visualizations in a Panel dashboard.
    """

    dataset_id: str
    filename: str
    min_genes: int
    selection_x1: float | None = None
    selection_x2: float | None = None
    selection_y1: float | None = None
    selection_y2: float | None = None

    def __init__(self, settings: Settings, **params):
        super().__init__(**params)

        self.settings = settings

        if pn.state.location is not None:
            pn.state.location.sync(
                self.settings,
                {
                    "dataset_id": "dataset_id",
                    "filename": "filename",
                    "min_genes": "min_genes",
                    "selection_x1": "selection_x1",
                    "selection_x2": "selection_x2",
                    "selection_y1": "selection_y1",
                    "selection_y2": "selection_y2",
                    "display_height": "height",
                    "display_width": "width",
                    "expression_min_clip": "expression_min_clip",
                },
            )

        self.dataset_id = secure_filename(self.settings.dataset_id)  # type: ignore
        self.filename = secure_filename(self.settings.filename)  # type: ignore
        self.min_genes = self.settings.min_genes  # type: ignore

        self.unique_clusters = None
        self.cluster_map = None
        self.color_map = None

        layout_height = 312  # 360px - tile header height
        if self.settings.display_height and self.settings.display_height > 0: # type: ignore
            layout_height: int = self.settings.display_height # type: ignore

        layout_width = 1100  # Default width
        if self.settings.display_width and self.settings.display_width > 0: # type: ignore
            layout_width: int = self.settings.display_width # type: ignore

        # When width is too short, things break down.
        if layout_width < 1100:
            layout_width = 1100
            layout_height = 312

        self.layout = pn.Column(pn.bind(self.init_data), width=layout_width)

        self.condensed_fig = dict(data=[], layout={})

        self.use_clusters = False
        self.use_clusters_switch = pn.widgets.Switch(
            name="Feature", value=self.use_clusters, sizing_mode="fixed", margin=(10, 10)
        )

        markdown_width = 400    # Manually measured.
        markdown_padding = 20  # Total left/right padding for the markdown pane
        self.condensed_intro = pn.pane.Markdown(
            "### Click the Expand icon in the top right corner to see all plots", width=markdown_width
        )

        # Using HTML to center the labels (https://github.com/holoviz/panel/issues/1313#issuecomment-1582731241)
        switch_content_width = 250
        self.switch_layout = pn.Row(
            pn.pane.HTML("""<label><strong>Gene Expression</strong></label>""", sizing_mode="fixed", margin=(10, 10)),
            self.use_clusters_switch,
            pn.pane.HTML("""<label><strong>Clusters</strong></label>""", sizing_mode="fixed", margin=(10, 10)),
            width=switch_content_width
        )

        spacer_width = layout_width - markdown_width - markdown_padding - switch_content_width
        if spacer_width < 0:
            spacer_width = 100

        self.condensed_pre = pn.Row(
            self.condensed_intro,
            pn.Spacer(width=spacer_width),
            self.switch_layout,
            height=30
        )

        # Build a row with the following elements:
        # 1) Blank image
        # 2) Normal plot
        # 3) Zoomed in plot
        # Above these is a toggle for switching 2) and 3) between gene expression and cluster plots

        self.condensed_pane = pn.pane.Plotly(
            self.condensed_fig,
            config={
                "doubleClick": "reset",
                "displayModeBar": True,
                "modeBarButtonsToRemove": buttonsToRemove,
            },
            height=275,
            sizing_mode="stretch_width"
        )

        # A condensed layout should only have the normal pane.
        self.plot_layout = pn.Column(
            self.condensed_pre,
            self.condensed_pane,
            width=layout_width,
            height=layout_height,
        )

        def refresh_figures_callback(value) -> None:
            self.use_clusters = value
            self.refresh_figures()

        pn.bind(refresh_figures_callback, self.use_clusters_switch, watch=True)

        def selection_callback(event):
            if self.condensed_fig_obj is None:
                return

            self.condensed_pane.object = self.condensed_fig_obj.selection_callback(
                event
            )

        pn.bind(selection_callback, self.condensed_pane.param.selected_data, watch=True)

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

            spatial_img_path = PANEL_CSV_CACHE_DIR / self.dataset_id / SPATIAL_IMAGE_NAME
            self.spatial_img = None
            if spatial_img_path.is_file():
                self.spatial_img = np.load(spatial_img_path)

            df_path = PANEL_CSV_CACHE_DIR / self.dataset_id / self.filename
            if not df_path.is_file():
                raise FileNotFoundError(f"Data file not found: {df_path}")

            self.df = pd.read_csv(df_path)

            if self.settings.expression_min_clip is not None:
                self.df = clip_expression_values(self.df, min_clip=self.settings.expression_min_clip)   # type: ignore

            yield self.loading_indicator(
                "Processing data to create plots. This may take a minute..."
            )

            self.refresh_dataframe(self.min_genes)

            yield self.refresh_figures()
        except Exception as e:
            logging.error(traceback.format_exc())
            yield pn.pane.Alert(f"Error: {e}", alert_type="danger")

    def apply_gene_filter(self, min_genes):
        """
        Apply a minimum-gene filter to the panel's dataframe.

        Sets self.min_genes to the provided value and, if that value differs from
        self.settings.min_genes, synchronizes self.settings.min_genes. Then filters
        self.df in-place to retain only rows where the "n_genes_by_counts" column
        is greater than or equal to the specified minimum.

        Parameters
        ----------
        min_genes : int
            Minimum number of genes (inclusive) required for a row to be kept.

        Raises
        ------
        AttributeError
            If self, self.df, or self.settings is not correctly initialized.
        KeyError
            If the "n_genes_by_counts" column does not exist in self.df.

        Notes
        -----
        - This method mutates self.min_genes, self.settings.min_genes (conditionally),
          and self.df.
        - The method does not return a value.
        """
        self.min_genes = min_genes

        # If the min_genes value has changed... sync to URL
        if self.min_genes != self.settings.min_genes:
            self.settings.min_genes = self.min_genes

        self.df = self.df[self.df["n_genes_by_counts"] >= self.min_genes]

    def refresh_dataframe(self, value: int):
        """
        Refresh the instance dataframe and derived cluster mappings after applying a gene filter.

        This method performs the following steps:
        - Calls self.apply_gene_filter(value) to (re)populate / filter self.df based on the provided integer value.
        - Computes and stores self.unique_clusters as the unique values from the "clusters" column of the (possibly filtered) dataframe.
        - Builds self.cluster_map: a dict mapping each category code from "clusters_cat_codes" to the corresponding cluster name (the first "clusters" value found for that code).
        - Builds self.color_map: a dict mapping each cluster name to the corresponding color (the first "colors" value found for that cluster).
        - Enforces categorical dtypes for the "clusters" and "clusters_cat_codes" columns (CSV input may not preserve these dtypes).

        Parameters
        ----------
        value : int
            Integer parameter forwarded to self.apply_gene_filter. The semantics of this value depend on the implementation
            of apply_gene_filter (e.g., it may select a gene index or threshold).

        Side effects
        ------------
        - Mutates self.df (it is expected that apply_gene_filter filters or updates self.df).
        - Sets/updates the attributes: self.unique_clusters, self.cluster_map, self.color_map.
        - Casts self.df["clusters"] and self.df["clusters_cat_codes"] to pandas.Categorical dtype in-place.

        Required dataframe columns
        --------------------------
        The method expects self.df to contain at least the following columns:
        - "clusters" : cluster names/labels
        - "clusters_cat_codes" : cluster categorical codes
        - "colors" : color values associated with clusters

        Raises
        ------
        AttributeError
            If self.apply_gene_filter is not defined.
        KeyError
            If any required column ("clusters", "clusters_cat_codes", "colors") is missing from self.df.
        IndexError
            If a code or cluster has no matching rows when building cluster_map or color_map.

        Performance
        -----------
        Linear in the number of rows of self.df for the filtering and mapping steps.

        Notes
        -----
        - The implementation assumes there is at least one representative row per cluster / cluster code; missing rows will raise an IndexError when accessing the first element.
        - Because CSV files often lose dtype information, this method explicitly casts cluster-related columns to categorical. Consider switching to a binary format (e.g., Parquet) to preserve dtypes across I/O.
        """
        self.apply_gene_filter(value)
        self.unique_clusters = self.df["clusters"].unique()
        self.cluster_map = {
            code: self.df[self.df["clusters_cat_codes"] == code]["clusters"].to_numpy()[0]
            for code in self.df["clusters_cat_codes"].unique()
        }
        self.color_map = {
            cluster: self.df[self.df["clusters"] == cluster][
                "colors"
            ].to_numpy()[0]
            for cluster in self.unique_clusters
        }

        # This should have been done but the CSV file does not preserve dtypes
        # ? Should we switch to parquet files?
        self.df["clusters"] = self.df["clusters"].astype("category")
        self.df["clusters_cat_codes"] = self.df["clusters_cat_codes"].astype("category")


    def refresh_figures(self):
        self.condensed_fig_obj = None

        if self.color_map is None or self.cluster_map is None:
            raise ValueError("Color map or cluster map is not initialized.")

        self.condensed_fig_obj = SpatialCondensedSubplot(
            self.settings,
            self.df,
            self.spatial_img,
            self.color_map,
            self.cluster_map,
            "Local",
            "YlGn",
            "YlOrRd",
            self.use_clusters,
            dragmode="select",
            #platform=self.platform,
        )
        self.condensed_fig = self.condensed_fig_obj.refresh_spatial_fig()
        self.condensed_pane.object = self.condensed_fig

        return self.plot_layout


### MAIN APP ###

# Reset logging level to "error" to suppress a bokeh "dropping patch" info message
# https://github.com/bokeh/bokeh/issues/13229
logging.getLogger().setLevel(logging.ERROR)
settings = Settings()

# If not params passed, just show OK as a way to test the app
if pn.state.location is None or not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()

else:
    # Create the app
    sp_panel = SpatialPanel(settings)
    sp_panel.servable(title="Spatial Data Viewer", location=True)
