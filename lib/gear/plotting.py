from plotly import exceptions
from plotly.colors import unlabel_rgb
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go

from itertools import cycle

# Purely kept so we can debug with sys.stderr
import sys

px.defaults.template = "none"
blank_template = go.layout.Template()

# mapping to ensure the right function is used.
PLOT_TYPE_TO_FUNCTION = {
    "bar":px.bar,
    "box":px.box,
    "histogram":px.histogram,
    "line":px.line,
    "scatter":px.scatter,
    "strip":px.strip,
    #"violin":px.violin,
    "contour":px.density_contour    # NOTE: Potentially add as secondary option to lay over scatter plot
}

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)


def _add_kwargs_info_to_annotations(fig, annotation_info):
    """Add various annotation info.  Updates 'fig' inplace."""
    fig.update_annotations(
        patch=annotation_info
    )

def _add_kwargs_info_to_coloraxes(fig, coloraxis_info):
    """Add various coloraxis info.  Updates 'fig' inplace."""
    fig.update_coloraxes(
        patch=coloraxis_info
    )

def _add_kwargs_info_to_layout(fig, layout_info):
    """Add various layout info.  Updates 'fig' inplace."""
    fig.update_layout(
        dict1=layout_info
    )

def _add_kwargs_info_to_traces(fig, trace_info):
    """Add various trace info.  Updates 'fig' inplace."""
    fig.update_traces(
        patch=trace_info
    )

def _add_kwargs_info_to_xaxes(fig, xaxis_info):
    """Add various xaxis info.  Updates 'fig' inplace."""
    fig.update_xaxes(
        patch=xaxis_info
    )

def _add_kwargs_info_to_yaxes(fig, yaxis_info):
    """Add various yaxis info.  Updates 'fig' inplace."""
    fig.update_yaxes(
        patch=yaxis_info
    )

def _add_marker_info_to_traces(fig, marker_info, plot_type, color_name=None):
    """Add various marker info as trace properties.  Updates 'fig' inplace."""

    if plot_type in ["bar"]:
        marker_info.pop('size', None)
    elif plot_type in ["violin"]:
        # Marker will make jitter points black.
        # Normally this affects fillcolor of plot, but this is manually set
        marker_info.setdefault("color","#000000")
    else:
        # In scatter plots, the 'size' was set during initial figure generation.
        # The UI Marker size slider will be used to set the minimum marker size
        if plot_type == "scatter" and not isinstance(marker_info["size"], int):
            marker_info.setdefault("sizemode", "diameter")
            marker_info.setdefault("sizeref", 0.5)
            marker_info.pop('size', None)

    if not color_name:
        # Set default purple color on trace markers
        marker_info.setdefault("color", "#401362")

    # After adding some extra marker properties, apply all trace updates.
    fig.update_traces(
        marker=marker_info,
    )

def _add_vertical_lines(fig, vline):
    """Add a vertical line to the figure.  Updates fig in-place."""
    fig.add_shape(dict(
        type= 'line',
        yref= 'paper', y0= 0, y1= 1,    # Entire length of axis
        xref= 'x', x0= vline["vl_pos"], x1= vline["vl_pos"],
        line = dict(dash=vline["vl_style"])
    ))

def _adjust_colorscale(plotting_args, colormap=None, palette=None):
    """Adjust the colorscale used for plotting."""
    # Add colorrange
    if colormap and isinstance(colormap, dict):
        plotting_args["color_discrete_map"] = colormap
    elif colormap and isinstance(colormap, list):
        plotting_args["color_continuous_scale"] = colormap
    else:
        # purple shades
        plotting_args["color_continuous_scale"] = [
            [0, 'rgb(218, 183, 193)'],
            [0.35, 'rgb(194, 137, 166)'],
            [0.5, 'rgb(169, 98, 151)'],
            [0.6, 'rgb(145, 66, 143)'],
            [0.7, 'rgb(105, 39, 122)'],
            [1, 'rgb(63, 19, 98)']
        ]
        # Palette selection supercedes "purples"
        if palette:
            plotting_args["color_continuous_scale"] = palette

    # There is no palette option for discrete mapping
    return plotting_args

def _aggregate_dataframe(df, x, y, facet_row=None, facet_col=None, color_name=None ):
    """Aggregate dataframe information under certain conditions."""

    priority_groups = _build_priority_groups(facet_row, facet_col, color_name, x)

    # Safeguard against grouping by an empty list
    if not priority_groups:
        return df

    grouped = df.groupby(priority_groups)

    # Discrete colorscale or no colorscale
    if not color_name or _is_categorical(df[color_name]):
        df = grouped.agg({
                y: ['mean', 'std']
            }) \
            .dropna(subset=[(y, "mean")]) \
            .loc[:, y] \
            .fillna(value={"std":0}) \
            .rename(columns=dict(mean=y)) \
            .reset_index()
    else:
        # Continuous colorscale
        df = grouped.agg({
                color_name: ['mean'],
                y: ['mean', 'std']
            }) \
            .dropna() \
            .reset_index()
    return df

def _build_priority_groups(facet_row=None, facet_col=None, color_name=None, x=None, y=None):
    """Determine group priority for "groupby" functions."""
    priority_groups = []
    # Add facet row or facet_col or both
    if facet_row and facet_col and not facet_row == facet_col:
        priority_groups.extend([facet_row, facet_col])
    elif (bool(facet_row) ^ bool(facet_col)):
        # XOR condition
        priority_groups.append(facet_row if facet_row else facet_col)
    elif facet_row and facet_row == facet_col:
        # Safeguard if both facets are the same... add only one
        priority_groups.append(facet_row)

    if color_name and color_name not in priority_groups:
        priority_groups.append(color_name)
    if x and x not in priority_groups:
        priority_groups.append(x)
    # Only added for contour plots
    if y and y not in priority_groups:
        priority_groups.append(y)

    return priority_groups

def _determine_annotation_shift(ds):
    """Determine number of pixels to shift annotation based on longest entry in dataseries."""
    if _is_categorical(ds):
        ds_list = ds.unique().tolist()
        longest_entry = len(max(ds_list, key = len))
        return -(longest_entry * 3 + 30)
    else:
        return -50

def _invalid_plot_type(plot_type):
    """Check for invalid plot types."""
    raise PlotError("Plot type {} is invalid!".format(plot_type))

def _is_categorical(series):
    """Return True if Dataframe series is categorical."""
    return series.dtype.name == 'category'

def _translate_and_scale(series, x):
    """Convert and return number from a linear scale using another linear scale."""
    unscaled_min = series.min()
    unscaled_max = series.max()

    # These were arbitrarily chosen. At least 1 will be added to the final total based on marker size increase
    NEW_MIN = 0
    NEW_MAX = 9

    return ( (NEW_MAX - NEW_MIN)*(x - unscaled_min) / (unscaled_max - unscaled_min) ) + NEW_MIN

def _truncate_ticktext(group_list):
    """Truncate a group of axis ticks to a specified length."""
    TRUNCATION_LEN=7    # How much of the original text to use (followed by ellipses)
    MAX_LEN_ALLOWED=10  # Any text over this limit will be truncated

    # If only 0 or 1 datapoints in group, categoryarray was not present
    if not group_list:
        return None

    new_ticktext = []
    for val in group_list:
        if len(val) > MAX_LEN_ALLOWED:
            new_ticktext.append("{}...".format(val[0:TRUNCATION_LEN]))
        else:
            new_ticktext.append(val)
    return new_ticktext

def _update_axis_titles(fig, df, x, y, facet_col=None, facet_row=None, x_title=None, y_title=None):
    """Update axis titles.  Edits "fig" inplace."""

    # Annotation defaults based on https://community.plotly.com/t/subplots-how-to-add-master-axis-titles/13927/6
    # NOTE: I've come to the conclusion that a one-size-fits-all solution is not possible for all graphs
    if facet_col:
        fig.update_xaxes(title=None)
        fig.add_annotation(
            x=0.5,
            y=-0,
            yshift=_determine_annotation_shift(df[x]),
            #yshift=-30,
            showarrow=False,
            text=x_title,
            name="x-title",
            xref="paper",
            yref="paper",
            yanchor="top",
        )

    if facet_row:
        fig.update_yaxes(title=None)
        fig.add_annotation(
            x=0,
            xshift=-40,
            y=0.5,
            showarrow=False,
            text=y_title,
            textangle=-90,
            name="y-title",
            xref="paper",
            yref="paper",
            xanchor="right",
        )

def _update_by_plot_type(fig, plot_type, force_overlay=False, use_jitter=False):
    """Updates specific to certain plot types.  Updates 'fig' inplace."""
    if plot_type == "violin":
        fig.update_traces(
            spanmode="hard",    # Do not extend violin tails beyond the min/max values

            # Jitter-based args (to make beeswarm plot)
            jitter=0.25 if use_jitter else None,
            points="all" if use_jitter else False,
            pointpos=0 if use_jitter else None,
        )
        fig.update_layout(
            violinmode='group'
        )
        if force_overlay:
            fig.update_layout(
                # Overlay is chosen if color dataseries is same as 'facet_col' dataseries
                violinmode='overlay'
            )
    elif plot_type == "bar":
        fig.update_layout(
            barmode='group'
        )
        if force_overlay:
            fig.update_layout(
                barmode='overlay'
            )
    elif plot_type in ["box", "strip"]:
        fig.update_layout(
            boxmode='group'
        )
        if force_overlay:
            fig.update_layout(
                boxmode='overlay'
            )

def generate_plot(df, x=None, y=None, z=None, facet_row=None, facet_col=None,
                      color_name=None, colormap=None, palette=None,
                      reverse_palette=False, category_orders=None,
                      plot_type='scatter', hide_x_labels=False, hide_y_labels=False,
                      hide_legend=None, text_name=None, jitter=False,
                      x_range=None, y_range=None, vlines=[], x_title=None, y_title=None,
                      **kwargs):
    """Generates and returns figure for facet grid."""

    # If replicates are present, use mean and stdev of expression data as datapoints
    if 'replicate' in df.columns and plot_type not in ['violin', 'contour']:
        df = _aggregate_dataframe(df, x, y, facet_row, facet_col, color_name)

    # Little bit of safeguarding with kwargs
    # keys for kwargs: 'annotations', 'coloraxes', 'layout', 'traces', 'xaxes', 'yaxes'

    kwargs.setdefault("traces", {"marker": {}})         # If traces does not exist
    kwargs["traces"].setdefault("marker", {})           # If markers does not exist within
    kwargs["traces"]["marker"].setdefault("size", 3)    # If size does not exist within

    # Collect all args in a dictionary for easy passing to plotting function
    plotting_args={"x":x
        , "y":y
        , "facet_row":facet_row
        , "facet_col":facet_col
        , "color":color_name
        , "category_orders": category_orders if category_orders else {}
        , "labels": {x:x_title, y:y_title, color_name:""}
        , "hover_name": text_name if text_name else y
        }

    # Ensure label is one of the labels that is not lost from "gropuby"
    if plotting_args["hover_name"] not in [x, y, facet_row, facet_col, color_name]:
        raise PlotError("Selected label {} is not the same as one of the 'x', 'y', 'facet', or 'color' conditions".format(plotting_args["hover_name"]))

    plotting_args = _adjust_colorscale(plotting_args, colormap, palette)
    plotting_args["hover_data"] = { col: False for col in df.columns.tolist() }

    # If jitter is needed for scatter plot, convert to a strip plot
    if plot_type == "scatter" and jitter:
        plot_type = "strip"

    # For scatter plots with a lot of datapoints, use WebGL rendering
    if plot_type == "scattergl":
        plot_type == "scatter"
        plotting_args["render_mode"] = "webgl"

    # For line plots, use svg render since webgl mode does not have spline as a valid line shape (as of plotly 4.14.3)
    if plot_type == "line":
        plotting_args["render_mode"] = "svg"
        plotting_args["line_shape"] = "spline"

    # Scatter plots are the only types that let you set marker size by group
    # TODO: SAdkins - this is ugly... come up with better way to handle 'integer size' vs 'size by group'
    if plot_type == "scatter":
        if not isinstance(kwargs['traces']['marker']['size'], int):
            group_size = kwargs['traces']['marker']['size']
            size_increase = kwargs['traces']['marker']['sizemin'] if "sizemin" in kwargs['traces']['marker'] else 1
            if _is_categorical(df[group_size]):
                # If categorical, pass a number equivalent to the group
                plotting_args["size"] = [i+size_increase for i in df.groupby(group_size).ngroup()]
            else:
                # If continuous data, convert data so that range is all positive integers
                scaled_range = [_translate_and_scale(df[group_size], i) for i in df[group_size]]
                plotting_args["size"] = [i+size_increase for i in scaled_range]

    # Some plottypes cannot use a certain colorscale
    if plot_type in ["strip", "box", "violin", "contour", "line"]:
        plotting_args.pop("color_continuous_scale", None)

    if plot_type in ["bar"]:
        plotting_args["error_y"] = "std" if "std" in df.columns else None

    if plot_type == 'contour':
        plotting_args["z"] = z
        plotting_args["histfunc"] = "avg"   # For expression data
        # NOTE: Almost feel these below should be under "_update_by_plot_type()"
        kwargs["traces"]["contours_coloring"] = "fill"
        # For expression data, start scale at 0
        if z == "raw_value":
            kwargs["traces"]["zmin"] = 0.0
            kwargs["traces"]["zmax"] = max(df[z].tolist())

    ### CREATE THE PLOT

    fig = None
    # Call the right Plotly express function based on the plot type
    func = PLOT_TYPE_TO_FUNCTION.get(plot_type)
    if func:
        try:
            fig = func(df, **plotting_args)
        except exceptions.PlotlyError as e:
            raise PlotError("Encountered error in plotting {}: {}".format(plot_type, str(e)))
    else:
        # SAdkins - An aside... this following code is why I wanted to use Plotly Express to generate the plots
        # TODO: put in function

        # Map indexes for subplot ordering.  Indexes start at 1 since plotting rows/cols start at 1
        facet_row_groups = category_orders[facet_row] if facet_row and facet_row in category_orders else []
        facet_row_indexes = {group: idx for idx, group in enumerate(facet_row_groups, start=1)}
        num_rows = len(facet_row_groups) if facet_row else 1
        facet_col_groups = category_orders[facet_col] if facet_col and facet_col in category_orders else []
        facet_col_indexes = {group: idx for idx, group in enumerate(facet_col_groups, start=1)}
        num_cols = len(facet_col_groups) if facet_col else 1

        # Make faceted plot
        fig = make_subplots(rows=num_rows
                , cols=num_cols
                , row_titles=facet_row_groups if facet_row else None
                , column_titles=facet_col_groups if facet_col else None
                )

        # Because of the reliance on the premade figure, this cannot go at the top of the page
        # like to the Plotly Express plot constructors
        PLOT_TYPE_TO_SPECIAL_FUNCTION = {
            "violin":fig.add_violin,
        }

        try:
            special_func = PLOT_TYPE_TO_SPECIAL_FUNCTION.get(plot_type)
        except exceptions.PlotlyError as e:
            raise PlotError("Encountered error in plotting {}: {}".format(plot_type, str(e)))
        if not special_func:
            _invalid_plot_type(plot_type)

        # TODO clean this up in a function
        new_plotting_args = { 'x': df[x]
                , "y": df[y]
                , "text":df[text_name] if text_name else y
                }
        new_plotting_args['line'] = dict(color='#401362')
        new_plotting_args['showlegend'] = False # Only add legend with color group present

        priority_groups = _build_priority_groups(facet_row, facet_col, color_name, x)
        if priority_groups:
            # Groupby will not include combinations with missing data.  This can result in missing traces for a group
            grouped = df.groupby(priority_groups)
            names_in_legend = {}
            # Name is a tuple of groupings, or a string if grouped by only 1 dataseries
            # Group is the 'groupby' dataframe
            for name, group in grouped:
                for k, v in {'x':x, 'y':y, 'text':text_name}.items():
                    new_plotting_args[k] = group[v] if v else group[y]

                # Quick plot-specific check
                if plot_type in ["violin"]:
                    if color_name and (palette or not _is_categorical(df[color_name])):
                        raise PlotError("ERROR: Tried to call continuous colorscale on violin plot.")

                # Each individual trace is a separate scalegroup to ensure plots are scaled correctly for violin plots
                new_plotting_args['scalegroup'] = name
                if isinstance(name, tuple):
                    new_plotting_args['scalegroup'] = "_".join(name)

                # If color dataseries is present, add some special configurations
                if color_name:
                    curr_color = name
                    if isinstance(name, tuple):
                        curr_color = name[priority_groups.index(color_name)]
                    new_plotting_args['name'] = curr_color

                    # If facets are present, a legend group trace can appear multiple times.
                    # Ensure it only shows once.
                    new_plotting_args['showlegend'] = True
                    if curr_color in names_in_legend:
                        new_plotting_args['showlegend'] = False
                    names_in_legend[curr_color] = True

                    new_plotting_args['line'] = dict(color='#000000')
                    new_plotting_args["legendgroup"] = curr_color
                    new_plotting_args["offsetgroup"] = curr_color   # Cleans up some weird grouping stuff, making plots thicker

                    if colormap:
                        # Use black outlines with colormap fillcolor. Pertains mostly to violin plots
                        new_plotting_args['fillcolor'] = colormap[curr_color]

                # Now determine which plot this trace should go to.  Facet column is first if row does not exist.
                # Note the "facet_row/col_indexes" enum command started indexing at 1, so no need to increment for 1-indexed subplots
                row_idx = 1
                col_idx = 1
                if isinstance(name, tuple):
                    row_idx = facet_row_indexes[name[0]] if facet_row else 1
                    # Use index in case row and col were identical
                    col_idx = facet_col_indexes[ name[priority_groups.index(facet_col)] ] if facet_col else 1
                else:
                    row_idx = facet_row_indexes[name] if facet_row else 1
                    col_idx = facet_col_indexes[name] if facet_col else 1

                special_func(**new_plotting_args, row=row_idx, col=col_idx)

        else:
            # Safeguard against grouping by an empty list
            # use dataframe instead
            special_func(**new_plotting_args, row=1, col=1)


        #TODO: Since graph_object plots don't need 'category_orders' decide if we can drop passing that to the px functions
        fig.update_layout(template=blank_template)
        # Only tick labels from the first axis should be shown
        if facet_row:
            fig.update_xaxes(showticklabels=False)
            max_x = "xaxis{}".format(num_rows if num_rows > 1 else "")
            fig['layout'][max_x]['showticklabels'] = True
        else:
            fig['layout']['xaxis']['title']['text'] = x_title
        if facet_col:
            fig.update_yaxes(showticklabels=False)
            fig['layout']['yaxis']['showticklabels'] = True
        else:
            fig['layout']['yaxis']['title']['text'] = y_title

        # Order the columns
        if x in category_orders:
            for key in fig['layout']:
                if key.startswith('xaxis'):
                    fig['layout'][key]['categoryarray'] = category_orders[x]
        elif y in category_orders:
            for key in fig['layout']:
                if key.startswith('yaxis'):
                    fig['layout'][key]['categoryarray'] = category_orders[y]

    ### POST-PLOT ADJUSTMENTS

    # Change facet annotations to read <group> instead of <cat> = <group>
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    # Update axes layouts
    fig.update_xaxes(dict(
        automargin=True,
        zeroline=False,
        showgrid=False,
        visible=False if hide_x_labels else True,
        range=x_range if x_range and not _is_categorical(df[x]) else None,
        #tickangle=-90 if facet_col else None
    ))

    # Truncate faceted column axis labels so annotation can fit
    if facet_col and _is_categorical(df[x]):
        fig.for_each_xaxis(
            lambda a: a.update(
                ticktext= _truncate_ticktext(a.categoryarray),
                tickvals= a.categoryarray,
            )
        )

    fig.update_yaxes(dict(
        automargin=True,
        zeroline=False,
        showgrid=False,
        visible=False if hide_y_labels else True,
        range=y_range if y_range and not _is_categorical(df[y]) else None
    ))

    # More general trace updates
    fig.update_traces(dict(
        opacity=None if plot_type in ["contour"] or "color_continuous_scale" in plotting_args else 0.6,
    ))

    # If 'facet_col' and 'color_name' dataseries are equal treat as if color series is not present (for plot grouping)
    # Originally included 'x' and 'color_name' but a plotly update fixed this
    force_overlay = True
    if color_name:
        force_overlay = False
        color_eq_col = False
        if facet_col and _is_categorical(df[color_name]):
            color_eq_col = df[color_name].equals(df[facet_col])
        color_eq_x = df[color_name].equals(df[x])
        force_overlay = color_eq_x or color_eq_col

    # Plottype-specific tweaks
    _update_by_plot_type(fig, plot_type, force_overlay, jitter)

    _update_axis_titles(fig, df, x, y, facet_col, facet_row, x_title, y_title)

    # Add vertical lines
    [_add_vertical_lines(fig, vl) for vl in vlines]

    # More general layout updates
    fig.update_layout(
        autosize=True,
        hovermode='closest',
        legend = dict(itemsizing="constant"),
        showlegend=False if hide_legend else None,   # 'None' just means do whatever plotly defaults to
    )

    # Use kwargs marker information to update figure information
    # Pop off marker info since it is already loaded.
    if plot_type not in ["contour"]:
        _add_marker_info_to_traces(fig, kwargs["traces"]["marker"], plot_type, color_name)
    kwargs["traces"].pop("marker", None)

    if reverse_palette:
        kwargs["coloraxes"]["reversescale"] = reverse_palette

    # Update this particular entity with kwargs information.  This is generally custom things the user wants
    # that is not avaiable in the general plotly configuration we want nor in dataset_curator options
    if "annotations" in kwargs:
        _add_kwargs_info_to_annotations(fig, kwargs["annotations"])
    if "coloraxes" in kwargs:
        _add_kwargs_info_to_coloraxes(fig, kwargs["coloraxes"])
    if "layout" in kwargs:
        _add_kwargs_info_to_layout(fig, kwargs["layout"])
    if "traces" in kwargs:
        _add_kwargs_info_to_traces(fig, kwargs["traces"])
    if "xaxes" in kwargs:
        _add_kwargs_info_to_xaxes(fig, kwargs["xaxes"])
    if "yaxes" in kwargs:
        _add_kwargs_info_to_yaxes(fig, kwargs["yaxes"])

    return fig


def get_config():
    """Get config for Plotly chart."""
    return dict(
        showLink=False,
        displaylogo=False,
        responsive=False,
        modeBarButtonsToRemove=[
            "zoom2d",
            "pan2d",
            #"select2d",
            #"lasso2d",
            "zoomIn2d",
            "zoomOut2d",
            "autoScale2d",
            # "resetScale2d",
            "hoverClosestCartesian",
            "hoverCompareCartesian",
            "zoom3d",
            "pan3d",
            "resetCameraDefault3d",
            "resetCameraLastSave3d",
            "hoverClosest3d",
            "orbitRotation",
            "tableRotation",
            "zoomInGeo",
            "zoomOutGeo",
            "resetGeo",
            "hoverClosestGeo",
            # "toImage",
            "sendDataToCloud",
            "hoverClosestGl2d",
            "hoverClosestPie",
            "toggleHover",
            "resetViews",
            "toggleSpikelines",
            "resetViewMapbox"
        ]
 )

def plotly_color_map(names):
    """
    Private plot helper method for generating colors
    for similar groups across subplots.
    Parameters
    ----------
    names:
        List of levels in a group.
    Returns
    -------
    Dictionary assigning group names to a color.
    """
    if len(names) > 1:
        d3_category_20 = [
            "#1f77b4",
            "#aec7e8",
            "#ff7f0e",
            "#ffbb78",
            "#2ca02c",
            "#98df8a",
            "#d62728",
            "#ff9896",
            "#9467bd",
            "#c5b0d5",
            "#8c564b",
            "#c49c94",
            "#e377c2",
            "#f7b6d2",
            "#7f7f7f",
            "#c7c7c7",
            "#bcbd22",
            "#dbdb8d",
            "#17becf",
            "#9edae5"
        ]
        color_generator = cycle(d3_category_20)
        return dict(zip(names, color_generator))
    else:
        return dict(zip(names, [rgb_to_hex(*unlabel_rgb('rgb(64,19,98)'))]))

def rgb_to_hex(r,g,b):
    hex = "#{:02x}{:02x}{:02x}".format(int(r),int(g),int(b))
    return hex
