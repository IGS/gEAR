import gc
import logging
import sys
import warnings
from pathlib import Path

import anndata
import colorcet as cc
import numpy as np
import pandas as pd
import panel as pn
import plotly.graph_objects as go
import plotly.io as pio
import scanpy as sc
import spatialdata as sd
from common import (
    Settings,
    SpatialFigure,
    normalize_searched_gene,
    sort_clusters,
)
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
        gene_symbol (str): The gene symbol to visualize.
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
        gene_symbol: str,
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
            gene_symbol,
            expression_name,
            expression_color,
            **params,
        )
        # Preserve the original dataframe for filtering
        self.orig_df = df

        self.final_col = 3 if spatial_img is not None else 2

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
        self.marker_size = int(1 + 2500 / (x_range + y_range))

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
        if self.spatial_img is None or self.platform == "xenium":
            # If there is no spatial image, use white line color for visibility
            linecolor = "white"

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
    gene_symbol: str
    min_genes: int
    projection_id: str | None = None
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
                    "gene_symbol": "gene_symbol",
                    "min_genes": "min_genes",
                    "projection_id": "projection_id",
                    "selection_x1": "selection_x1",
                    "selection_x2": "selection_x2",
                    "selection_y1": "selection_y1",
                    "selection_y2": "selection_y2",
                    "display_height": "height",
                    "display_width": "width",
                },
            )

        self.dataset_id = self.settings.dataset_id  # type: ignore
        self.gene_symbol = self.settings.gene_symbol  # type: ignore
        self.min_genes = self.settings.min_genes  # type: ignore
        self.projection_id = self.settings.projection_id  # type: ignore

        self.platform = None  # Will be set in prep_sdata()

        layout_height = 312  # 360px - tile header height
        if self.settings.display_height and self.settings.display_height > 0: # type: ignore
            layout_height: int = self.settings.display_height # type: ignore

        layout_width = 1100  # Default width
        if self.settings.display_width and self.settings.display_width > 0: # type: ignore
            layout_width: int = self.settings.display_width # type: ignore

        self.layout = pn.Column(pn.bind(self.init_data), width=layout_width)

        self.condensed_fig = dict(data=[], layout={})

        self.use_clusters = False
        self.use_clusters_switch = pn.widgets.Switch(
            name="Feature", value=self.use_clusters, sizing_mode="fixed", margin=(10, 10)
        )

        markdown_width = 700
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

        spacer_width = layout_width - markdown_width - switch_content_width
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
        yield self.loading_indicator("Processing data file...")

        def create_adata_pkg():
            try:
                self.prep_sdata()
                adata = self.prep_adata()
            except ValueError:
                raise

            adata_pkg = {"adata": adata, "img_name": self.spatial_obj.img_name}
            return adata_pkg

        adata_cache_label = f"{self.dataset_id}_adata"

        # Load the Anndata object (+ image name) from cache or create it if it does not exist, with a 1-week time-to-live
        try:
            adata_pkg = pn.state.as_cached(
                adata_cache_label, create_adata_pkg, ttl=CACHE_EXPIRATION
            )
        except ValueError as e:
            yield pn.pane.Alert(f"Error: {e}", alert_type="danger")
            raise

        self.dataset_adata_orig = adata_pkg["adata"]  # Original dataset adata
        self.dataset_adata = (
            self.dataset_adata_orig.copy()
        )  # Copy we manipulate (filtering, etc.)
        self.adata = self.dataset_adata.copy()  # adata object to use for plotting

        # Modify the adata object to use the projection ID if it exists
        if self.projection_id:
            self.adata = self.create_projection_adata()

        self.adata_orig = (
            self.adata.copy()
        )  # This is to restore when the min_genes slider is changed

        self.spatial_img = None
        if adata_pkg["img_name"]:
            # In certain conditions, the image multi-array may need to be squeezed so that PIL can read it
            self.spatial_img = self.adata.uns["spatial"][adata_pkg["img_name"]][
                "images"
            ]["hires"].squeeze()

        yield self.loading_indicator(
            "Processing data to create plots. This may take a minute..."
        )

        self.refresh_dataframe(self.min_genes)
        yield self.refresh_figures()

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
            print("Invalid or unsupported spatial data type")
            print("Supported types: {0}".format(SPATIALTYPE2CLASS.keys()))
            sys.exit(1)

        # Use uploader class to determine correct helper functions
        self.spatial_obj = SPATIALTYPE2CLASS[self.platform]()
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
            raise ValueError(
                f"No cells found with at least {self.min_genes} genes. Choose a different gene or lower the filter."
            )

        sc.pp.normalize_total(self.dataset_adata, inplace=True)
        sc.pp.log1p(self.dataset_adata)

        self.dataset_adata.var_names_make_unique()

    def create_gene_df(self):
        adata = self.adata
        dataset_genes = set(self.adata.var["gene_symbol"].unique())
        norm_gene_symbol = normalize_searched_gene(dataset_genes, self.gene_symbol)
        if norm_gene_symbol is None:
            raise ValueError(
                f"Gene '{self.gene_symbol}' not found in the dataset. Please choose a different gene."
            )
        self.norm_gene_symbol = norm_gene_symbol
        gene_filter = adata.var.gene_symbol == self.norm_gene_symbol
        selected = adata[:, gene_filter]
        selected.var.index = pd.Index(["raw_value"])
        dataframe = selected.to_df()

        # Add spatial coords from adata.obsm
        X, Y = (0, 1)
        dataframe["spatial1"] = adata.obsm["spatial"].transpose()[X].tolist() # type: ignore
        dataframe["spatial2"] = adata.obsm["spatial"].transpose()[Y].tolist()

        # Add cluster info
        if "clusters" not in selected.obs:
            raise ValueError("No cluster information found in adata.obs")

        dataframe["clusters"] = selected.obs["clusters"].astype("category")
        dataframe["clusters_cat_codes"] = dataframe["clusters"].cat.codes.astype(
            "category"
        )

        self.cluster_map = {
            code: dataframe[dataframe["clusters_cat_codes"] == code][
                "clusters"
            ].to_numpy()[0]
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

        return

    def refresh_figures(self):
        self.condensed_fig_obj = None

        self.condensed_fig_obj = SpatialCondensedSubplot(
            self.settings,
            self.df,
            self.spatial_img,
            self.color_map,
            self.cluster_map,
            self.norm_gene_symbol,
            "Local",
            "YlGn",
            "YlOrRd",
            self.use_clusters,
            dragmode="select",
            platform=self.platform,
        )
        self.condensed_fig = self.condensed_fig_obj.refresh_spatial_fig()
        self.condensed_pane.object = self.condensed_fig

        gc.collect()  # Collect garbage to free up memory

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
settings = Settings()

# If not params passed, just show OK as a way to test the app
if pn.state.location is None or not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()

else:
    # Create the app
    sp_panel = SpatialPanel(settings)
    sp_panel.servable(title="Spatial Data Viewer", location=True)
