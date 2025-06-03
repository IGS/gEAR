import base64
import sys
from io import BytesIO
from typing import Literal

import datashader as ds
import numpy as np
import pandas as pd
import param
from PIL import Image, ImageEnhance
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from xarray import DataArray

### Functions


def normalize_searched_gene(gene_set, chosen_gene) -> str | None:
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


def sort_clusters(clusters) -> list:
    """
    Sort clusters by number if numerical, otherwise by name.
    """
    try:
        sorted_clusters = sorted(clusters, key=lambda x: int(x))
    except Exception:
        sorted_clusters = sorted(clusters, key=lambda x: str(x))
    return sorted_clusters


### Classes
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
        projection_id (param.String): Projection ID to display.
        save (param.Boolean): If true, save this configuration as a new display.
        display_name (param.String): Display name for the saved configuration.
        make_default (param.Boolean): If true, make this the default display.
    """

    gene_symbol = param.String(doc="Gene symbol to display")
    dataset_id = param.String(doc="Dataset ID to display")
    min_genes = param.Integer(
        doc="Minimum number of genes per observation", default=0, bounds=(0, 500)
    )
    projection_id = param.String(doc="Projection ID to display", allow_None=True)
    selection_x1 = param.Number(doc="left selection range", allow_None=True)
    selection_x2 = param.Number(doc="right selection range", allow_None=True)
    selection_y1 = param.Number(doc="upper selection range", allow_None=True)
    selection_y2 = param.Number(doc="lower selection range", allow_None=True)


class SpatialFigure:
    """
    SpatialFigure is a class for creating interactive spatial gene expression visualizations using Plotly and Datashader.

    This class provides methods to generate figures with subplots for spatial images, gene expression scatter plots, and cluster annotations. It supports overlaying spatial transcriptomics data on tissue images, visualizing gene expression levels, and displaying cluster assignments. The class handles image encoding, axis configuration, color mapping, and efficient aggregation of large datasets.

    Attributes:
        settings (Settings): Configuration object containing selection and plot settings.
        df (pd.DataFrame): DataFrame containing spatial coordinates, expression values, and cluster assignments.
        spatial_img (np.ndarray | None): Optional background image as a NumPy array.
        color_map (dict): Mapping from cluster identifiers to colors.
        cluster_map (dict): Mapping from cluster category codes to cluster names.
        gene_symbol (str): The gene symbol being visualized.
        expression_name (str): Name of the expression metric (e.g., "expression", "relative expression").
        expression_color (str): Color scale for expression visualization.
        dragmode (str): Plotly drag mode (e.g., "select").
        marker_size (int): Size of scatter plot markers.
        base64_string (str | None): Base64-encoded PNG string of the background image.
        PLOT_WIDTH (int): Width of the plot canvas for aggregation.
        PLOT_HEIGHT (int): Height of the plot canvas for aggregation.
        range_x1, range_x2, range_y1, range_y2 (float): Ranges for spatial axes.
        fig (go.Figure): The Plotly figure object.

    Methods:
        __init__(...): Initializes the SpatialFigure with data, settings, and visualization parameters.
        make_expression_scatter(expression_agg): Creates a Plotly scatter plot for gene expression.
        add_cluster_traces(cluster_agg, cluster_col): Adds cluster scatter traces to the figure.
        add_cluster_scatter(fig, dataframe, clusters, color_map, x_col, y_col, col, row): Adds scatter traces for each cluster.
        create_datashader_agg(x, y): Aggregates data using Datashader for efficient visualization.
        make_fig(): Constructs and configures the complete Plotly figure with subplots and traces.
        encode_pil_image(pil_img): Encodes a PIL image as a base64 PNG data URI.
        convert_background_image(): Converts and crops the background image, encoding it as base64.
        create_subplots(num_cols, titles): Creates subplots with the specified number of columns.
        create_subplots_two(titles): Creates a two-column subplot layout.
        create_subplots_three(titles): Creates a three-column subplot layout.
        update_axes(): Configures axes ranges, titles, and appearance.
        add_image_trace(): Adds the background image to all subplot columns.
        get_legend_font_size(): Determines legend font size based on cluster name lengths.
        update_selection_ranges(todict): Updates selection ranges or dictionary based on current settings.
    """

    # TODO: Change "expression" to "relative expression" for projections

    def __init__(
        self,
        settings: Settings,
        df: pd.DataFrame,
        spatial_img: np.ndarray | None,
        color_map: dict,
        cluster_map: dict,
        gene_symbol: str,
        expression_name: str,
        expression_color: str,
        dragmode: str = "select",
        **kwargs: dict | None,
    ) -> None:
        self.settings = settings
        self.df = df
        self.spatial_img = spatial_img  # None or numpy array
        self.color_map = color_map
        self.cluster_map = cluster_map
        self.gene_symbol = gene_symbol
        self.fig = go.Figure()
        self.expression_name = expression_name
        self.expression_color = expression_color
        self.dragmode = dragmode
        self.marker_size = 2
        self.base64_string = None
        self.PLOT_WIDTH = 200
        self.PLOT_HEIGHT = 200

        # kwargs
        self.platform = kwargs.get("platform", None)  # e.g., "xenium"

        # https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.experimental.to_legacy_anndata.html
        # (see section Matching of spatial coordinates and pixel coordinates)
        # The downscaled image (x,y) coordinate matches (spatial2, spatial1) in the dataframe
        # We traditionally use the 1-coord as x and the 2-coord as y, so we need to address these in various plots.
        self.range_x1 = 0 if self.spatial_img is not None else min(df["spatial1"])
        self.range_x2 = (
            self.spatial_img.shape[1]
            if self.spatial_img is not None
            else max(df["spatial1"])
        )
        self.range_y1 = 0 if self.spatial_img is not None else min(df["spatial2"])
        self.range_y2 = (
            self.spatial_img.shape[0]
            if self.spatial_img is not None
            else max(df["spatial2"])
        )

    def make_expression_scatter(
        self, expression_agg: DataArray, cbar_loc: float, expression_color: str | None
    ) -> go.Scatter:
        """
        Creates a Plotly Scatter plot visualizing gene expression values in spatial coordinates.

        Parameters
        ----------
        expression_agg : pandas.Series or pandas.DataFrame
            Aggregated expression values indexed by spatial coordinates. Must be convertible to a DataFrame
            with columns 'spatial1', 'spatial2', and expression values.

        Returns
        -------
        plotly.graph_objs.Scatter
            A Plotly Scatter object with marker color representing expression values and spatial coordinates
            as x and y axes.
        """
        dataframe = self.df

        agg_df = expression_agg.to_dataframe(name="raw_value")
        # Drop missing values
        agg_df = agg_df.dropna()
        dataframe = agg_df.reset_index()

        if expression_color is None:
            expression_color = self.expression_color

        return go.Scatter(
            x=dataframe["spatial1"],
            y=dataframe["spatial2"],
            mode="markers",
            marker=dict(
                cauto=False,  # Do not adjust the color scale
                cmin=0,  # No expression
                cmax=dataframe["raw_value"].max(),  # Max expression value
                # ? This would not apply if data is log-transformed
                color=dataframe["raw_value"],
                colorscale=expression_color,  # You can choose any colorscale you like
                size=self.marker_size,  # Adjust the marker size as needed
                colorbar=dict(
                    len=1,  # Adjust the length of the colorbar
                    thickness=15,  # Adjust the thickness of the colorbar (default is 30)
                    x=cbar_loc,  # Adjust the x position of the colorbar
                ),
            ),
            showlegend=False,
            unselected=dict(
                marker=dict(opacity=1)
            ),  # Helps with viewing unselected regions
            hovertemplate="Expression: %{marker.color:.2f}<extra></extra>",
        )

    def add_cluster_traces(
        self, cluster_agg: DataArray, cluster_col: int, show_legend=True
    ) -> None:
        """
        Adds cluster scatter traces to the figure based on aggregated cluster data.

        Parameters:
            cluster_agg (DataArray): An xarray DataArray containing aggregated cluster information.
            cluster_col (int): The column index to use for plotting clusters.

        Returns:
            None
        """
        fig = self.fig

        agg_df = cluster_agg.to_dataframe(name="clusters_cat_codes")
        # Drop rows where "clusters" is False
        agg_df = agg_df[agg_df["clusters_cat_codes"]]
        # The columns we want are in the multi-index, so we need to make them into a dataframe
        dataframe = agg_df.index.to_frame(index=False)

        color_map = self.color_map

        # map cluster cat_codes to cluster names
        dataframe["clusters"] = dataframe["clusters_cat_codes"].map(self.cluster_map)
        try:
            # If clusters are float, convert to int (for later sorting)
            dataframe["clusters"] = dataframe["clusters"].astype(int).astype("category")
        except ValueError:
            # If clusters are not numeric, keep as category
            dataframe["clusters"] = dataframe["clusters"].astype("category")

        self.add_cluster_scatter(
            fig,
            dataframe,
            color_map,
            "spatial1",
            "spatial2",
            col=cluster_col,
            row=1,
            show_legend=show_legend,
        )
        # TODO: Add click event for removing equivalent points from expression plot if legend item is clicked

        # Update the scatter trace to ensure unselected points are visible
        fig.update_traces(
            unselected=dict(marker=dict(opacity=1)),
            selector=dict(type="scatter"),
            col=cluster_col,
            row=1,
        )

    def add_cluster_scatter(
        self,
        fig: go.Figure,
        dataframe: pd.DataFrame,
        color_map: dict,
        x_col: str,
        y_col: str,
        col: int,
        row: int,
        show_legend: bool = True,
    ) -> None:
        """
        Adds scatter plot traces for each cluster to a Plotly figure.

        Parameters:
            fig (go.Figure): The Plotly figure to which the scatter traces will be added.
            dataframe (pd.DataFrame): The DataFrame containing the data to plot.
            clusters (list): A list of cluster identifiers to plot.
            color_map (dict): A mapping from cluster identifiers to colors.
            x_col (str): The name of the column to use for x-axis values.
            y_col (str): The name of the column to use for y-axis values.
            col (int): The column position in the subplot grid.
            row (int): The row position in the subplot grid.

        Returns:
            None
        """

        # Ensure clusters are categorical and sorted
        dataframe["clusters"] = dataframe["clusters"].astype("category")
        unique_clusters = dataframe["clusters"].unique()
        sorted_clusters = sort_clusters(unique_clusters)

        for cluster in sorted_clusters:
            cluster_data = dataframe[dataframe["clusters"] == cluster]
            fig.add_trace(
                go.Scatter(
                    x=cluster_data[x_col],
                    y=cluster_data[y_col],
                    marker=dict(color=color_map[cluster], size=self.marker_size),
                    mode="markers",
                    name=str(cluster),
                    text=cluster_data["clusters"],
                    hovertemplate="Cluster: %{text}<extra></extra>",
                    showlegend=show_legend,  # Show legend only if specified
                ),
                col=col,
                row=row,
            )

    def create_datashader_agg(self, x: str, y: str) -> DataArray:
        """
        Aggregates data points from the DataFrame using Datashader, producing a summary for visualization.

        Parameters:
            x (str): The name of the column to use for the x-axis.
            y (str): The name of the column to use for the y-axis.

        Returns:
            xarray.DataArray: An aggregated DataArray where each pixel contains the maximum 'raw_value' expression
            and a boolean indicating the presence of each cluster (by 'clusters_cat_codes').

        Notes:
            - The x and y values are swapped for aggregation compared to plotting.
            - The aggregation computes the maximum 'raw_value' per pixel and tracks cluster membership.
        """

        # Create a canvas with the specified width and height
        cvs = ds.Canvas(plot_width=self.PLOT_WIDTH, plot_height=self.PLOT_HEIGHT)

        # NOTE: This seems to affect the expression aggregation but does not matter for the clusters
        # We need to swap the x and y values for aggregation compared to what we eventually plot.

        agg = cvs.points(
            self.df,
            x=x,
            y=y,
            agg=ds.summary(
                expression=ds.max("raw_value"),
                clusters=ds.by("clusters_cat_codes", ds.any()),  # type: ignore
            ),
        )
        return agg

    def _setup_fig(self) -> None:
        """
        Creates and configures a Plotly figure with subplots for spatial and expression data visualization.

        Returns:
            plotly.graph_objs._figure.Figure: The configured Plotly figure object.
        """

        # domain is adjusted whether there are images or not
        if self.spatial_img is not None:
            self.create_subplots(num_cols=3)  # sets self.fig
            self.expression_col = 2
            self.cluster_col = 3
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.29]),
                xaxis2=dict(domain=[0.33, 0.62]),  # Leave room for cluster annotations
                xaxis3=dict(domain=[0.71, 1]),
            )
        else:
            self.create_subplots(num_cols=2)  # sets self.fig
            self.expression_col = 1
            self.cluster_col = 2
            self.fig.update_layout(
                xaxis=dict(domain=[0, 0.45]),
                xaxis2=dict(domain=[0.55, 1]),  # Leave room for cluster annotations
            )
        self.update_axes()

        if self.spatial_img is not None:
            self.add_image_trace()

    def _add_traces_to_fig(self) -> go.Figure:
        # Get max domain for axis 1 or axis 2
        if self.expression_col == 1:
            max_x = self.fig.layout.xaxis.domain[1]  # type: ignore
        else:
            max_x = self.fig.layout.xaxis2.domain[1]  # type: ignore

        agg = self.create_datashader_agg(x="spatial2", y="spatial1")

        cbar_loc = max_x + 0.005  # Leave space for colorbar

        self.fig.add_trace(
            self.make_expression_scatter(
                agg["expression"],
                cbar_loc=cbar_loc,
                expression_color=self.expression_color,
            ),
            row=1,
            col=self.expression_col,
        )
        self.add_cluster_traces(agg["clusters"], self.cluster_col)

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

        self.fig.update_layout(
            plot_bgcolor=plot_bgcolor,
        )

        return self.fig

    def encode_pil_image(self, pil_img: Image.Image) -> str:
        """
        Encodes a PIL Image object as a base64-encoded PNG data URI string.

        Args:
            pil_img (Image.Image): The PIL Image to encode.

        Returns:
            str: A string containing the base64-encoded PNG image, prefixed with the appropriate data URI scheme.
        """
        prefix = "data:image/png;base64,"
        with BytesIO() as stream:
            pil_img.save(stream, format="png")
            return prefix + base64.b64encode(stream.getvalue()).decode("utf-8")

    def convert_background_image(self) -> None:
        """
        Converts the spatial image array to a cropped PIL image based on the selected range,
        then encodes the cropped image as a base64 string for efficient loading and storage.

        Returns:
            None
        """

        if self.spatial_img is None:
            raise ValueError("Spatial image is not provided. Should not have been called.")

        img = self.spatial_img

        # (Copilot Recommendation) - Ensure image is uint8 for consistent contrast
        if img.dtype != np.uint8:
            img = img.astype(np.float32)
            img = 255 * (img - img.min()) / (img.max() - img.min() + 1e-8)
            img = img.astype(np.uint8)

        pil_img = Image.fromarray(img)
        pil_img = pil_img.crop(
            (self.range_x1, self.range_y1, self.range_x2, self.range_y2)
        )

        # Convert to RGB.
        # For some reason, converting the mode directly in Image.fromarray() throws "Not enough image data" error
        pil_img = pil_img.convert("RGB")

        # If "xenium" platform, contrast after conversion is too high, so we reduce it
        if self.platform == "xenium":
            enhancer = ImageEnhance.Contrast(pil_img)
            pil_img = enhancer.enhance(0.33)  # Reduce contrast to 33%

        self.base64_string = self.encode_pil_image(pil_img)

    def create_subplots(self, num_cols: int = 2, titles: tuple = ()) -> None:
        """
        Creates subplots with either two or three columns, delegating to the appropriate method.

        Args:
            num_cols (int, optional): Number of columns for the subplots. Must be either 2 or 3. Defaults to 2.
            titles (tuple, optional): Titles for each subplot. Defaults to an empty tuple.

        Raises:
            ValueError: If num_cols is not 2 or 3.

        Returns:
            None
        """
        if num_cols not in (2, 3):
            raise ValueError("num_cols must be 2 or 3")
        if num_cols == 3:
            self.create_subplots_three(titles)
        else:
            self.create_subplots_two(titles)

    def create_subplots_two(self, titles: tuple) -> None:
        """
        Creates a figure with two subplots arranged in a single row.

        Args:
            titles (tuple): A tuple containing titles for the two subplots. If empty or None, default titles are used.

        Returns:
            None
        """
        if not titles:
            titles = (f"{self.expression_name} Expression", "Clusters")

        self.fig = make_subplots(
            rows=1,
            cols=2,
            column_titles=titles,
            horizontal_spacing=0.1,
        )

    def create_subplots_three(self, titles: tuple) -> None:
        """
        Creates a figure with three subplots arranged in a single row, each with a specified title.

        Args:
            titles (tuple): A tuple containing titles for the three subplots. If not provided or empty,
                            default titles are used: "Image Only", "<expression_name> Expression", and "Clusters".

        Returns:
            None
        """
        if not titles:
            titles = (
                "Image Only",
                f"{self.expression_name} Expression",
                "Clusters",
            )
        self.fig = make_subplots(
            rows=1,
            cols=3,
            column_titles=titles,
            horizontal_spacing=0.1,
        )

    def update_axes(self) -> None:
        """
        Updates the x and y axes of the figure (`self.fig`) with custom settings.

        Returns:
            None
        """
        self.fig.update_xaxes(
            range=[self.range_x1, self.range_x2],
            title_text="spatial1",
            showgrid=False,
            showticklabels=False,
            ticks="",
        )
        # y-axis needs to be flipped. 0,0 is the top left corner
        self.fig.update_yaxes(
            range=[self.range_y2, self.range_y1],
            title_text="spatial2",
            showgrid=False,
            showticklabels=False,
            ticks="",
            title_standoff=0,
        )

    def add_image_trace(self) -> None:
        """
        Adds a background image trace to each subplot column in the figure.

        Returns:
            None
        """
        self.convert_background_image()
        # Note that passing the spatial img array to the "z" property instead of the base64-encoded string causes a big performance hit
        image_trace = go.Image(
            source=self.base64_string, x0=self.range_x1, y0=self.range_y1
        )
        for col in range(1, 4):
            self.fig.add_trace(image_trace, row=1, col=col)

    def get_legend_font_size(self) -> Literal[10] | Literal[12]:
        """
        Determines the appropriate font size for the legend based on the length of the longest cluster name.

        Returns:
            Literal[10] | Literal[12]: The font size to use for the legend. Returns 12 if the longest cluster name
            is 10 characters or fewer, otherwise returns 10.
        """
        # Get longest cluster name for legend
        longest_cluster = max(self.df["clusters"].astype(str).apply(len))
        font_size = 12
        if longest_cluster > 10:
            font_size = 10
        return font_size

    def update_selection_ranges(self) -> bool:
        """
        Updates the selection ranges based on the current selection coordinates.

        Checks if the selection coordinates for x and y axes have changed. If they have,
        updates the `selections_dict` attribute with the new selection boundaries and returns True.
        If the selection coordinates have not changed, returns False.

        Returns:
            bool: True if the selection ranges were updated, False otherwise.
        """

        if (
            self.settings.selection_x1 != self.settings.selection_x2
            or self.settings.selection_y1 != self.settings.selection_y2
        ):
            self.selections_dict = dict(
                x0=self.settings.selection_x1,
                x1=self.settings.selection_x2,
                y0=self.settings.selection_y1,
                y1=self.settings.selection_y2,
            )

            return True
        return False
