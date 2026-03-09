from pathlib import Path
import sys

from bokeh.models import HoverTool
import datashader as ds
from holoviews.operation.datashader import spread
import numpy as np
import pandas as pd
import param
from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[2]
www_path = gear_root.joinpath("www")
PANEL_CSV_CACHE_DIR = www_path / "cache" / "spatial_panel"
SPATIAL_IMAGE_NAME = "spatial_img.npy"

### Functions

def create_spatial_plot(df, agg, x_col='spatial1', y_col='spatial2', color_col='raw_value', cmap='YlOrRd', is_categorical=False, width=300, height=200):
    """Generates a Datashaded spatial plot colored by expression of the specified gene."""

    plot = df.hvplot.points(
        x=x_col, y=y_col, c=color_col,
        rasterize=True, aggregator=agg,
        cmap=cmap, frame_width=width, frame_height=height,
        xaxis=None, yaxis=None,

        # 1. Kill the colorbar if it's categorical
        colorbar=not is_categorical,

        # 2. Force the legend on
       # legend='right' if is_categorical else False
    )

    tools = ['box_select']
    if color_col == "raw_value":
        label_name = "Expression"
        # @image is a special variable that datashader uses to store the aggregated value for the hovered pixel
        custom_hover = HoverTool(tooltips=[
                (label_name, "@image")
            ])
    else:
        label_name = color_col.title()
        # @image is a special variable that datashader uses to store the aggregated value for the hovered pixel
        custom_hover = HoverTool(tooltips=[
                (label_name, "@field")
            ], limit=3)

    tools.append(custom_hover)  # type: ignore

    plot = plot.opts(
            tools=tools,
            active_tools=['box_select'],
            default_tools=[]
        )

    return spread(plot, px=5)

def create_umap_plot(df, color_col, cmap, is_categorical=False, width=400, height=300):
    """Generates a Datashaded UMAP."""
    agg = ds.count_cat(color_col) if is_categorical else ds.mean(color_col)

    return df.hvplot.points(
        x='UMAP_1', y='UMAP_2', c=color_col,
        rasterize=True, aggregator=agg, dynspread=True,
        cmap=cmap, frame_width=width, frame_height=height,
        xaxis=None, yaxis=None, title=f"UMAP: {color_col}",
    )

def create_violin_plot(df, y_col, group_col='cluster', cmap='Category10', width=900, height=300):
    """Generates standard bokeh violin plots (no Datashader needed here)."""
    return df.hvplot.violin(
        y=y_col, by=group_col, c=group_col, cmap=cmap,
        ylabel='Expression', xlabel='Annotation Cluster',
        title=f"Expression Distribution: {y_col}",
        frame_width=width, frame_height=height, legend=False
    ).opts(violin_fill_alpha=0.7)

def clip_expression_values(dataframe: pd.DataFrame, min_clip: float | None=None, max_clip: float | None=None) -> pd.DataFrame:
    """
    Clip values in the DataFrame's "raw_value" column.

    Parameters
    ----------
    dataframe : pd.DataFrame
        DataFrame containing a "raw_value" column with numeric values to be clipped.
    min_clip : float | None, optional
        Minimum value to clip to. If None, no lower clipping is applied.
    dataframe["raw_value"] = dataframe["raw_value"].clip(lower=min_clip, upper=max_clip)
    max_clip : float | None, optional
        Maximum value to clip to. If None, no upper clipping is applied.

    Returns
    -------
    pd.DataFrame
        The same DataFrame instance passed in, with its "raw_value" column replaced by the clipped values.

    Notes
    -----
    - This function mutates the input DataFrame in place by assigning to dataframe["raw_value"].
    - Uses pandas.Series.clip which preserves NaNs and is vectorized for performance.
    - A KeyError will be raised if the "raw_value" column is not present in the DataFrame.

    Examples
    --------
    >>> df = pd.DataFrame({"raw_value": [-5, 0, 2.5, 10]})
    >>> clip_expression_values(df, min_clip=0.0, max_clip=5.0)
    >>> df["raw_value"].tolist()
    [0.0, 0.0, 2.5, 5.0]
    """
    return dataframe

def has_selection(settings) -> bool:
    """
    Return ``True`` when the bounds stored in ``self.settings`` describe a
    non‑degenerate rectangle.

    The legacy convention used by the panel app was that *all four*
    selection values would be equal when no box was supplied (either
    via the UI or the URL).  The boolean return value allows callers to
    apply zoom/mirroring only when there really is something to zoom to.
    """

    return not (
        settings.selection_x1 == settings.selection_x2 == settings.selection_y1 == settings.selection_y2
    )

def normalize_expression_name(filename) -> str:
        """
        Extract and normalize an expression name from a filename.

        If the filename stem follows the pattern <uuid4>_<str>, extracts the <str> part.
        Additionally, if the extracted name is "unweighted", it is replaced with "Pattern".

        Args:
            filename (str or Path): The filename to process.

        Returns:
            str: The normalized expression name.

        Examples:
            >>> normalize_expression_name("550e8400-e29b-41d4-a716-446655440000_myexpression.txt")
            'myexpression'
            >>> normalize_expression_name("550e8400-e29b-41d4-a716-446655440000_unweighted.txt")
            'Pattern'
            >>> normalize_expression_name("expression.txt")
            'expression'
        """
        expression_name = str(Path(filename).stem)
        # If expression_name has a pattern of <uuid4>_<str>, extract the <str> part
        if "_" in expression_name:
            parts = expression_name.split("_", 1)
            # Test for the UUID pattern too
            if len(parts[0]) == 36 and parts[0].count("-") == 4:
                expression_name = parts[1]
                if expression_name == "unweighted":
                    expression_name = "Pattern"

        return expression_name

def retrieve_dataframe(dataset_id, filename) -> pd.DataFrame:
    """
    Retrieve a dataframe from a CSV file in the spatial data cache.

    Args:
        dataset_id (str): The identifier of the dataset. Will be sanitized using secure_filename.
        filename (str): The name of the CSV file to retrieve. Will be sanitized using secure_filename.

    Returns:
        pd.DataFrame: The dataframe loaded from the specified CSV file.

    Raises:
        FileNotFoundError: If the CSV file does not exist at the expected cache path.

    Note:
        This function attempts to load a spatial image file if it exists at the dataset cache path,
        though the loaded image is not currently used or returned.
    """
    dataset_id = secure_filename(dataset_id)  # type: ignore
    filename = secure_filename(filename)  # type: ignore

    df_path = PANEL_CSV_CACHE_DIR / dataset_id / filename
    if not df_path.is_file():
        raise FileNotFoundError(f"Data file not found: {df_path}")

    return pd.read_csv(df_path)

def retrieve_image_array(dataset_id) -> np.ndarray | None:
    """
    Retrieve a spatial image array from the cache for a given dataset ID.

    Args:
        dataset_id (str): The identifier of the dataset. Will be sanitized using secure_filename.

    Returns:
        np.ndarray | None: The image array if it exists, otherwise None.
    """
    dataset_id = secure_filename(dataset_id)  # type: ignore
    spatial_img_path = PANEL_CSV_CACHE_DIR / dataset_id / SPATIAL_IMAGE_NAME
    if spatial_img_path.is_file():
        return np.load(spatial_img_path)
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
    """

    filename = param.String(doc="Filename for the dataframe to retrieve")
    dataset_id = param.String(doc="Dataset ID to display")
    min_genes = param.Integer(
        doc="Minimum number of genes per observation", default=0, bounds=(0, 500)
    )
    projection_id = param.String(doc="Projection ID to display", allow_None=True)
    selection_x1 = param.Number(doc="left selection range", allow_None=True)
    selection_x2 = param.Number(doc="right selection range", allow_None=True)
    selection_y1 = param.Number(doc="upper selection range", allow_None=True)
    selection_y2 = param.Number(doc="lower selection range", allow_None=True)
    display_height = param.Integer(doc="Height of the display in pixels", allow_None=True)
    display_width = param.Integer(doc="Width of the display in pixels", allow_None=True)
    expression_min_clip = param.Number(doc="Minimum expression value to clip", allow_None=True)

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