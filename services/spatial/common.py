from pathlib import Path
import sys

import datashader as ds
import numpy as np
import pandas as pd
import param
from bokeh.models import CustomJS, CustomJSHover, HoverTool, Span
from holoviews.operation.datashader import spread
import holoviews as hv
from werkzeug.utils import secure_filename

gear_root = Path(__file__).resolve().parents[2]
www_path = gear_root.joinpath("www")
PANEL_CSV_CACHE_DIR = www_path / "cache" / "spatial_panel"
SPATIAL_IMAGE_NAME = "spatial_img.npy"

DEFAULT_PLOT_WIDTH = 275
DEFAULT_PLOT_HEIGHT = 300

# NOTE: When you see the word "glyph" thrown around, that is Bokeh's terminology for the actual visual representation of a data point.
# Examples are circle, square, line, etc.
# A "renderer" is the object that takes the data and draws it using a glyph.

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

def link_crosshairs(figures, color='#333333', line_dash='dashed', line_width=1):
    """
    Synchronizes a crosshair (a vertical + horizontal guide line) across multiple
    already-rendered Bokeh figures that share the same coordinate space.

    Pure client-side JS (a shared Span pair per figure, updated via a single
    CustomJS on mousemove) rather than a HoloViews PointerXY stream, which
    would round-trip to the Bokeh server on every mouse move.

    Parameters
    ----------
    figures : list[bokeh.plotting.figure]
        The concrete Bokeh figures to link. They should share the same x/y
        coordinate system for the crosshair position to line up across them.
    """
    vlines, hlines = [], []
    for fig in figures:
        vline = Span(location=0, dimension='height', line_color=color,
                    line_dash=line_dash, line_width=line_width, visible=False)
        hline = Span(location=0, dimension='width', line_color=color,
                    line_dash=line_dash, line_width=line_width, visible=False)
        fig.add_layout(vline)
        fig.add_layout(hline)
        vlines.append(vline)
        hlines.append(hline)

    move_callback = CustomJS(args=dict(vlines=vlines, hlines=hlines), code="""
        const x = cb_obj.x;
        const y = cb_obj.y;
        for (let i = 0; i < vlines.length; i++) {
            vlines[i].location = x;
            vlines[i].visible = true;
            hlines[i].location = y;
            hlines[i].visible = true;
        }
    """)
    leave_callback = CustomJS(args=dict(vlines=vlines, hlines=hlines), code="""
        for (let i = 0; i < vlines.length; i++) {
            vlines[i].visible = false;
            hlines[i].visible = false;
        }
    """)

    for fig in figures:
        fig.js_on_event('mousemove', move_callback)
        fig.js_on_event('mouseleave', leave_callback)

def link_ranges(figures):
    """
    Links pan/zoom together across multiple Bokeh figures by sharing the same
    x_range/y_range model instances, so panning or box/wheel-zooming any one
    of them moves all the others in lockstep. This is the standard Bokeh
    "linked panning" technique.

    Because it's a literal shared model (not values kept in sync via a
    callback), Bokeh keeps every figure referencing it in sync automatically,
    for both the live client and this server-side session.

    Parameters
    ----------
    figures : list[bokeh.plotting.figure]
        The concrete Bokeh figures to link. They should share the same x/y
        coordinate system -- sharing ranges across figures with different
        coordinate systems will make them all zoom to whatever the first
        figure's extent is, which is rarely what you want.
    """
    if len(figures) < 2:
        return
    reference = figures[0]
    for fig in figures[1:]:
        fig.x_range = reference.x_range
        fig.y_range = reference.y_range

def _prepare_category_renderers(cmap, registry=None):
    """
    Returns a HoloViews plot hook for categorical (NdOverlay-based) spatial
    plots that does two things once the plot is actually rendered:

    1. Re-applies the intended category->color mapping directly onto each
       category's Bokeh glyph. Works around a hvplot/HoloViews issue
       where `by=` grouping combined with a dict `cmap` does not reliably
       preserve the category->color correspondence, leading to colors
       being assigned to the wrong category.
    2. If `registry` (a dict) is provided, populates it with
       {category_name: bokeh_glyph_renderer} so a caller can wire up
       per-category show/hide from an external legend afterward
    """
    def hook(plot, element):
        for key, subplot in getattr(plot, 'subplots', {}).items():
            category = key[0] if isinstance(key, tuple) else key
            renderer = subplot.handles.get('glyph_renderer')
            if renderer is None:
                continue

            if isinstance(cmap, dict):
                hexcolor = cmap.get(category, '#CCCCCC')
                if hasattr(renderer.glyph, 'fill_color'):
                    renderer.glyph.fill_color = hexcolor
                if hasattr(renderer.glyph, 'line_color') and renderer.glyph.line_color is not None:
                    renderer.glyph.line_color = hexcolor

            if registry is not None:
                registry[category] = renderer
    return hook


def create_spatial_plot(df:pd.DataFrame, x_col:str='spatial1', y_col:str='spatial2', color_col:str='raw_value', cmap:str='fire_r',
                        is_categorical:bool=False, title:str|None=None, cbar_max:float|None=None, category_renderers:dict|None=None,
                        min_height:int=DEFAULT_PLOT_HEIGHT, min_width:int=DEFAULT_PLOT_WIDTH, radius:int=3):
    """Generates a Datashaded spatial plot colored by expression of the specified gene."""

    plot_kwargs = dict(
        x=x_col, y=y_col,
        rasterize=False, #aggregator=agg,
        cmap=cmap,
        xaxis=None, yaxis=None,
        title=title,

        #Kill the colorbar if it's categorical
        colorbar=not is_categorical,
        clabel="",  # No label for the colorbar
        clim = (0, cbar_max) if cbar_max is not None else None,
        legend=False,    # using a ghost legend so the legend does not squish the plot

        # Responsiveness if the browser is resized
        responsive=True, # Automatically stretches to fill its container
        min_height=min_height,   # Set a floor so plots don't collapse to 0px
        min_width=min_width    # Tells Bokeh the canvas cannot drop below 275px (prevent squishing if the layout for the display is not full width)
    )

    if is_categorical:
        # `by=` splits this into one renderer per category
        # instead of one combined glyph, which is what makes per-category
        # show/hide possible downstream.
        plot = df.hvplot.points(by=color_col, **plot_kwargs)
    else:
        plot = df.hvplot.points(c=color_col, **plot_kwargs)


    label_name = "Expression" if color_col == "raw_value" else color_col.title()
    custom_hover = HoverTool(
        tooltips=[(label_name, f"@{color_col}")],
    )
    default_tools=["box_zoom", "wheel_zoom", "pan", "reset"]
    active_tools=["wheel_zoom", "pan"]

    hooks = [autohide_toolbar, fix_colorbar_hook]
    if is_categorical:
        # Always run this for categorical plots -- it's also what fixes the
        # by=/cmap color assignment, not just the registry population.
        hooks.append(_prepare_category_renderers(cmap, category_renderers))

    opts_kwargs = dict(hooks=hooks)
    if not is_categorical:
        # This only applies to the continuous/expression plot.
        opts_kwargs["colorbar_opts"] = {"width": 12}    # thin the colorbar out.

    plot = plot.opts(**opts_kwargs)

    # Add some opts that are scoped to the Points element itself, not the overall plot.
    # These particular options will override any defaults that hvplot may try to add in.
    plot = plot.opts(hv.opts.Points(
        radius=radius, line_color=None,
        tools=[custom_hover], default_tools=default_tools, active_tools=active_tools,
        toolbar="below",
        ))
    return plot

def create_umap_plot(df:pd.DataFrame, color_col:str="raw_value", cmap:str="cividis_r", is_categorical:bool=False, title:str|None=None, cbar_max:float|None=None, category_renderers:dict|None=None, radius:float=0.15):
    """Generates a Datashaded UMAP."""

    color_col_title = title if title else color_col.title()

    plot_kwargs = dict(
        x='UMAP1', y='UMAP2',
        rasterize=False, #aggregator=agg,
        cmap=cmap,
        xaxis="bottom", yaxis="left",
        xlabel="UMAP_1", ylabel="UMAP_2",
        title=f"UMAP: {color_col_title}",

        #Kill the colorbar if it's categorical
        colorbar=not is_categorical,
        clabel="",  # No label for the colorbar
        clim = (0, cbar_max) if cbar_max is not None else None,
        legend=False,    # using a ghost legend so the legend does not squish the plot

        # Responsiveness if the browser is resized
        responsive=True, # Automatically stretches to fill its container
        height=DEFAULT_PLOT_HEIGHT    # Tells Bokeh the canvas cannot drop below 275px (prevent squishing if the layout for the display is not full width)
    )

    if is_categorical:
        # `by=` splits this into one renderer per category
        # instead of one combined glyph, which is what makes per-category
        # show/hide possible downstream.
        plot = df.hvplot.points(by=color_col, **plot_kwargs)
    else:
        plot = df.hvplot.points(c=color_col, **plot_kwargs)

    label_name = "Expression" if color_col == "raw_value" else color_col.title()
    custom_hover = HoverTool(
        tooltips=[(label_name, f"@{color_col}")],
    )
    default_tools=["box_zoom", "wheel_zoom", "pan", "reset"]
    active_tools=["wheel_zoom", "pan"]

    hooks = [autohide_toolbar, fix_colorbar_hook]
    if is_categorical:
        # Always run this for categorical plots -- it's also what fixes the
        # by=/cmap color assignment, not just the registry population.
        hooks.append(_prepare_category_renderers(cmap, category_renderers))

    opts_kwargs = dict(hooks=hooks)
    if not is_categorical:
        # This only applies to the continuous/expression plot.
        opts_kwargs["colorbar_opts"] = {"width": 12}    # thin the colorbar out.

    plot = plot.opts(**opts_kwargs)

    # Add some opts that are scoped to the Points element itself, not the overall plot.
    # These particular options will override any defaults that hvplot may try to add in.
    plot = plot.opts(hv.opts.Points(
        radius=radius, line_color=None,
        tools=[custom_hover], default_tools=default_tools, active_tools=active_tools,
        toolbar="right",
        ))
    return plot

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

def compute_aggregation_params(df, x_col, y_col, target_markers=80_000, min_dim=275):
    """
    Computes (width, height, radius) for create_datashader_agg + create_spatial_plot.

    width/height are proportional to the data's x/y extent, so bin spacing
    comes out equal in both dimensions -- an aggregation resolution that
    isn't proportional to the data's aspect ratio produces different spacing
    in x vs y, and a single radius can't match both, causing visible tearing.
    radius is derived from that spacing so markers tile without gaps or overlap.
    """
    x_extent = df[x_col].max() - df[x_col].min()
    y_extent = df[y_col].max() - df[y_col].min()
    aspect = x_extent / y_extent if y_extent else 1.0

    height = max(min_dim, int(round((target_markers / aspect) ** 0.5)))
    width = max(min_dim, int(round(height * aspect)))

    spacing = y_extent / height
    radius = spacing / 2 * 1.05  # slight overlap avoids hairline gaps

    return width, height, radius

def create_datashader_agg(df, x: str, y: str, width:int=DEFAULT_PLOT_WIDTH, height:int=DEFAULT_PLOT_HEIGHT) -> "DataArray":
    """
    Aggregates data points from the DataFrame using Datashader, producing a summary for visualization.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the data to aggregate.
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
    cvs = ds.Canvas(plot_width=width, plot_height=height)

    # NOTE: This seems to affect the expression aggregation but does not matter for the clusters
    # We need to swap the x and y values for aggregation compared to what we eventually plot.

    df["clusters_cat_codes"] = df["clusters_cat_codes"].astype("category")

    agg = cvs.points(
        df,
        x=x,
        y=y,
        agg=ds.summary(
            expression=ds.max("raw_value"),
            clusters=ds.by("clusters_cat_codes", ds.any()),  # type: ignore
        ),
    )
    return agg

def create_clusters_df(agg):
    agg_df = agg.to_dataframe(name="clusters_cat_codes")
    agg_df = agg_df[agg_df["clusters_cat_codes"]]
    # The columns we want are in the multi-index, so we need to make them into a dataframe
    final_df =  agg_df.index.to_frame(index=False)
    return final_df

def create_expression_df(agg):
    agg_df = agg.to_dataframe(name="raw_value")
    # Drop missing values
    agg_df = agg_df.dropna()
    final_df = agg_df.reset_index()
    return final_df

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