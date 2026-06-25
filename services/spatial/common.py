import sys
from pathlib import Path

import datashader as ds
import numpy as np
import pandas as pd
import param
from bokeh.models import CustomJSHover, HoverTool
from holoviews.operation.datashader import dynspread, spread
from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[2]
www_path = gear_root.joinpath("www")
PANEL_CSV_CACHE_DIR = www_path / "cache" / "spatial_panel"
SPATIAL_IMAGE_NAME = "spatial_img.npy"

### Functions

def autohide_toolbar(plot, element):
    plot.state.toolbar.autohide = True

def fix_colorbar_hook(plot, element):
    """
    Forces the Bokeh colorbar to use the colormap of the top-most (selected) layer,
    bypassing the HoloViews unselected_color hijack.
    """
    try:
        from bokeh.models import ColorBar
        # Locate the colorbar in the plot layout
        colorbars = [obj for layout in ['right', 'left', 'above', 'below']
                     for obj in getattr(plot.state, layout, [])
                     if isinstance(obj, ColorBar)]

        if colorbars:
            cb = colorbars[0]
            # Find all image renderers that possess a color_mapper
            renderers = [r for r in plot.state.renderers
                         if hasattr(r, 'glyph') and hasattr(r.glyph, 'color_mapper')]

            if renderers:
                # The selected data is drawn last, meaning it is the last renderer in the list.
                # Force the colorbar to sync with it.
                cb.color_mapper = renderers[-1].glyph.color_mapper
    except Exception:
        pass

def create_spatial_plot(df, agg, x_col='spatial1', y_col='spatial2', color_col='raw_value', cmap='YlOrRd', is_categorical=False, title=None, mode="standard"):
    """Generates a Datashaded spatial plot colored by expression of the specified gene."""

    # Fixes a Holoviews bug where the linker callback cannot find the categorical metadata dimension that Datashader named
    # Gemini suggested
    if is_categorical and isinstance(cmap, dict):
        if hasattr(df[color_col], 'cat'):
            categories = df[color_col].cat.categories
        else:
            categories = sorted(df[color_col].unique())
        cmap = [cmap.get(cat, '#CCCCCC') for cat in categories]

    plot = df.hvplot.points(
        x=x_col, y=y_col, c=color_col,
        rasterize=True, aggregator=agg,
        cmap=cmap,
        xaxis=None, yaxis=None,
        title=title,

        #Kill the colorbar if it's categorical
        colorbar=not is_categorical,
        clabel="",  # No label for the colorbar
        legend=False,    # using a ghost legend so the legend does not squish the plot

        # Responsiveness if the browser is resized
        responsive=True, # Automatically stretches to fill its container
        min_height=300,   # Set a floor so plots don't collapse to 0px
        min_width=275    # Tells Bokeh the canvas cannot drop below 275px (prevent squishing if the layout for the display is not full width)
    )

    default_tools = ["box_zoom", "wheel_zoom", "pan", "reset"]
    if mode == "expanded":
        default_tools = ['box_select', 'reset']
    if color_col == "raw_value":
        label_name = "Expression"
        # Intercept the NaN values and round the long floats
        js_code_expr = """
            if (value === "NaN") {
                return "No data";
            }
            return value.toFixed(2);
        """
        formatter = CustomJSHover(code=js_code_expr)
        custom_hover = HoverTool(
            tooltips=[(label_name, "@image{custom}")],
            formatters={"@image": formatter}
        )
    else:
        label_name = color_col.title()

        # Unfortunately the datashader array only has the aggregated counts,
        # so we need to do some JS magic (thanks Gemini) to get the hover to show the category name instead of the count
        categories = list(df[color_col].cat.categories)

        # Build the JavaScript to find the dominant category in the pixel
        js_code = f"""
            const counts = value;
            const cats = {categories};
            let max_val = -1;
            let max_idx = -1;

            // Find the index of the highest count
            for (let i = 0; i < counts.length; i++) {{
                if (counts[i] > max_val) {{
                    max_val = counts[i];
                    max_idx = i;
                }}
            }}

            // Return the category name (and optionally the cell count!)
            if (max_val > 0) {{
                return cats[max_idx]; // + " (" + max_val + " cells)";
            }}
            return "No Data";
        """

        # Create the formatter and apply it strictly to the @image field
        formatter = CustomJSHover(code=js_code)
        custom_hover = HoverTool(
            tooltips=[(label_name, "@image{custom}")],
            formatters={"@image": formatter}
        )

    plot = plot.opts(
            tools=[custom_hover],
            default_tools=default_tools,
            toolbar="below",
            colorbar_opts={"width": 12},    # thin the colorbar out.
            hooks=[autohide_toolbar, fix_colorbar_hook]
        )

    # NOTE: Cluster downsampling will have a striped appearance called the Moiré Interference Pattern.
    # I tried to use dynspread to remedy it, but the data points end up being too light.
    # Mostly an issue with very dense data, like Visium HD
    return spread(plot, px=2, shape="square")


def create_umap_plot(df, agg, color_col, cmap, is_categorical=False, title=None):
    """Generates a Datashaded UMAP."""

    # Fixes a Holoviews bug where the linker callback cannot find the categorical metadata dimension that Datashader named
    # Gemini suggested
    if is_categorical and isinstance(cmap, dict):
        if hasattr(df[color_col], 'cat'):
            categories = df[color_col].cat.categories
        else:
            categories = sorted(df[color_col].unique())
        cmap = [cmap.get(cat, '#CCCCCC') for cat in categories]

    color_col_title = title if title else color_col.title()

    plot =  df.hvplot.points(
        x='UMAP1', y='UMAP2', c=color_col,
        aggregator=agg, cmap=cmap,
        xaxis="bottom", yaxis="left",
        xlabel="UMAP_1", ylabel="UMAP_2",
        title=f"UMAP: {color_col_title}",
        #Kill the colorbar if it's categorical
        colorbar=not is_categorical,
        clabel="",  # No label for the colorbar
        rasterize=True,
        responsive=True,
        legend=False,    # using a ghost legend so the legend does not squish the plot
        height=300,
    )

    if is_categorical:
        label_name = color_col.title()

        # Unfortunately the datashader array only has the aggregated counts,
        # so we need to do some JS magic (thanks Gemini) to get the hover to show the category name instead of the count
        categories = list(df[color_col].cat.categories)

        # Build the JavaScript to find the dominant category in the pixel
        js_code = f"""
            const counts = value;
            const cats = {categories};
            let max_val = -1;
            let max_idx = -1;

            // Find the index of the highest count
            for (let i = 0; i < counts.length; i++) {{
                if (counts[i] > max_val) {{
                    max_val = counts[i];
                    max_idx = i;
                }}
            }}

            // Return the category name (and optionally the cell count!)
            if (max_val > 0) {{
                return cats[max_idx]; // + " (" + max_val + " cells)";
            }}
            return "No Data";
        """

        # Create the formatter and apply it strictly to the @image field
        formatter = CustomJSHover(code=js_code)
        custom_hover = HoverTool(
            tooltips=[(label_name, "@image{custom}")],
            formatters={"@image": formatter}
        )
    else:
        label_name = "Expression"
        # Intercept the NaN values and round the long floats
        js_code_expr = """
            if (value === "NaN") {
                return "No data";
            }
            return value.toFixed(2);
        """
        formatter = CustomJSHover(code=js_code_expr)
        custom_hover = HoverTool(
            tooltips=[(label_name, "@image{custom}")],
            formatters={"@image": formatter}
        )

    default_tools = ['box_select', 'lasso_select', 'reset']
    plot = plot.opts(
            tools=[custom_hover],
            active_tools=['box_select'],
            default_tools=default_tools,
            #data_aspect=1,  # square aspect ratio for UMAP
            xticks=0, yticks=0,    # No ticks.
            colorbar_opts={"width": 12},    # thin the colorbar out.
            hooks=[autohide_toolbar, fix_colorbar_hook]
        )

    return spread(plot)

def create_violin_plot(df, y_col, group_col='cluster', cmap='Category10', title=None):
    """Generates standard bokeh violin plots (no Datashader needed here)."""
    plot_title = title if title else f"Expression Distribution: {y_col}"

    # Bokeh will try to wrap in a list, so we need to convert to a flattened list of hex colors.
    if isinstance(cmap, dict):
        # Determine the exact order of categories Bokeh will use
        if hasattr(df[group_col], 'cat'):
            categories = df[group_col].cat.categories
        else:
            categories = sorted(df[group_col].unique())

        cmap = [cmap.get(cat, '#CCCCCC') for cat in categories]

    # Sort df by the group_col
    df = df.sort_values(by=group_col)

    plot = df.hvplot.violin(
        y=y_col, by=group_col, c=group_col, cmap=cmap,
        ylabel='Expression', xlabel='Annotation Cluster',
        title=plot_title,
        min_height=400, responsive=True, legend=False
    )

    return plot.opts(
        violin_fill_alpha=0.7,
        xrotation=45,
        default_tools = ['box_select', 'lasso_select', 'reset'],
        hooks=[autohide_toolbar] # Keep the toolbar hidden until hover
    )

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
    if "raw_value" not in dataframe.columns:
        raise KeyError("DataFrame must contain a 'raw_value' column to clip.")

    dataframe["raw_value"] = dataframe["raw_value"].clip(lower=min_clip, upper=max_clip)
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
    projection_id = param.String(doc="Projection ID to display", allow_None=True)
    selection_x1 = param.Number(doc="left selection range", allow_None=True)
    selection_x2 = param.Number(doc="right selection range", allow_None=True)
    selection_y1 = param.Number(doc="upper selection range", allow_None=True)
    selection_y2 = param.Number(doc="lower selection range", allow_None=True)
    expression_min_clip = param.Number(doc="Minimum expression value to clip", allow_None=True)

    min_genes = param.Integer(default=0, doc="Minimum number of genes expressed to include a cell observation", bounds=(0, 500))
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