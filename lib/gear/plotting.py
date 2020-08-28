from __future__ import absolute_import

from plotly import exceptions, optional_imports
import plotly.colors as clrs
from plotly.figure_factory import utils
from plotly.subplots import make_subplots
from plotly.colors import unlabel_rgb
from itertools import cycle, chain
import math, copy
from numbers import Number

import sys

pd = optional_imports.get_module('pandas')

TICK_COLOR = '#969696'
AXIS_TITLE_COLOR = '#0f0f0f'
AXIS_TITLE_SIZE = 12
GRID_COLOR = '#ffffff'
LEGEND_COLOR = '#ffffff'
PLOT_BGCOLOR = '#ffffff'
ANNOT_RECT_COLOR = '#ffffff'
PAPER_BGCOLOR = 'rgba(255, 255, 255, 0)'
LEGEND_BORDER_WIDTH = 0
LEGEND_ANNOT_X = 1.05
LEGEND_ANNOT_Y = 0.5
MAX_TICKS_PER_AXIS = 5
THRES_FOR_FLIPPED_FACET_TITLES = 10
GRID_WIDTH = 1
PLOT_LOGGING = False

VALID_TRACE_TYPES = ['scatter', 'scattergl', 'histogram', 'bar', 'box', 'line', 'violin']

CUSTOM_LABEL_ERROR = (
    "If you are using a dictionary for custom labels for the facet row/col, "
    "make sure each key in that column of the dataframe is in your facet "
    "labels. The keys you need are {}"
)

def _adjust_plot_for_jitter(trace, jitter, curr_color):
    """Make plot have jitter effect, which requires a bit of manipulation."""
    # See 'strip' def at:
    # https://github.com/plotly/plotly.py/blob/master/packages/python/plotly/plotly/express/_chart_types.py

    trace['type'] = "box"
    trace["jitter"] = jitter
    trace["boxpoints"] = "all"
    trace["pointpos"] = 0   # centers points in box
    #trace["showlegend"] = False

    # This wipes out the actual boxplot, matching it with the plot bgcolor
    trace["line"] = {"color":PAPER_BGCOLOR, "width": 0}
    trace["fillcolor"] = PAPER_BGCOLOR
    trace["whiskerwidth"] = 0
    trace["boxmean"] = False

    # Default will hover on both points and invisible boxplot
    trace["hoveron"] = "points"

    # Boxplots have no mode
    trace.pop("mode")
    return trace

def _append_to_numerical_trace(trace, group, color_name, x, y, trace_type):
    """Append extra stuff to a trace with a numerical color label"""

    trace["marker"]["color"] =  list(group[color_name])
    trace["marker"]["coloraxis"] = "coloraxis"  # Makes all subplots share the same colorbar

    if 'replicate' in group.columns and trace_type != 'violin':
        #print(group[color_name]["mean"])
        #print(list(group[color_name]["mean"]))
        trace["marker"]["color"] =  list(group[color_name]["mean"])
        if x:
            trace['x'] = group[x]['mean']
        if y:
            trace['y'] = group[y]['mean']

    return trace

def _build_priority_groups(facet_row, facet_col, color_name, x):
    """Determine group priority for "groupby" functions."""
    priority_groups = []
    # Add facet row or facet_col or both
    if facet_row and facet_col:
        priority_groups.extend([facet_row, facet_col])
    elif bool(facet_row) ^ bool(facet_col):
        # XOR condition
        priority_groups.append(facet_row if facet_row else facet_col)
    if color_name:
        priority_groups.append(color_name)
    # Always add x
    priority_groups.append(x)
    return priority_groups

def _calculate_dtick(min_range, max_range, range_are_numbers, user_dtick):
    """Determine step interval between tick labels (dtick)."""
    dtick = 1
    if range_are_numbers:
        if user_dtick:
            return user_dtick, min_range, max_range

        min_range = math.floor(min_range)
        max_range = math.ceil(max_range)

        # extend widen frame by 5% on each side
        min_range -= 0.05 * (max_range - min_range)
        max_range += 0.05 * (max_range - min_range)

        if PLOT_LOGGING:
            print("DEBUG: Adjusted ranges are min_range:{} max_range:{}".format(min_range, max_range), file=sys.stderr)

        return math.floor(
            (max_range - min_range) / MAX_TICKS_PER_AXIS
        ), min_range, max_range
    return dtick, min_range, max_range

def _create_trace(group, x, y, trace_type, kwargs_trace, kwargs_marker):
    """Create a trace for a plot some general characteristics."""
    # Start creating the trace dictionary
    trace = dict(
        type="scatter" if trace_type == 'line' else trace_type,
        **kwargs_trace
    )

    if x:
        trace['x'] = group[x]
    if y:
        trace['y'] = group[y]

    # Add some error bars if stdev is a datapoint for a barplot
    if 'std' in group.columns and trace_type == 'bar':
        trace['error_y'] = dict(
            type='data',
            array=group['std']
        )

    # Specific trace-type adjustments
    if trace_type in ['scatter', 'scattergl']:
        trace['mode'] = 'markers'
    elif trace_type in ["line"]:
        trace['mode'] = 'lines'
        trace['line'] = dict(shape='spline')
    elif trace_type in ["violin"]:
        trace['scalemode'] = "count"

    trace['marker'] = dict(**kwargs_marker)
    return trace

def _is_categorical(series):
    """Return True if Dataframe series is categorical."""
    return series.dtype.name == 'category'

def _is_flipped(num):
    """Do the axes labels need to be flipped?"""
    if num >= THRES_FOR_FLIPPED_FACET_TITLES:
        flipped = True
    else:
        flipped = False
    return flipped


def _return_label(original_label, facet_labels, facet_var):
    """Return either the default label or a custom facet label."""
    if isinstance(facet_labels, dict):
        label = facet_labels[original_label]
    elif isinstance(facet_labels, str):
        label = '{}: {}'.format(facet_var, original_label)
    else:
        label = original_label
    return label


def _annotation_dict(text, lane, num_of_lanes, SUBPLOT_SPACING, row_col='col',
                     flipped=True):
    """Customize properties for a new annotation layer."""
    l = (1 - (num_of_lanes - 1) * SUBPLOT_SPACING) / (num_of_lanes)
    xanchor = 'center'
    yanchor = 'middle'
    name = 'y-facet'
    textangle = 0
    x = (lane - 1) * (l + SUBPLOT_SPACING) + 0.5 * l
    y = 1.03

    # if not flipped, just change a couple things
    if row_col == 'row':
        y = (lane - 1) * (l + SUBPLOT_SPACING) + 0.5 * l
        x = 1.03
        textangle = 90
        name = 'x-facet'    # Used to easily identify annotation for later modification

    # If flipped, we change everything
    if flipped:
        # row_col is 'col'
        yanchor = 'bottom'
        y = 1.0
        textangle = 270
        if row_col == 'row':
            xanchor = 'left'
            yanchor = 'middle'
            y = (lane - 1) * (l + SUBPLOT_SPACING) + 0.5 * l
            x = 1.0
            textangle = 0

    annotation_dict = dict(
        textangle=textangle,
        name=name,
        xanchor=xanchor,
        yanchor=yanchor,
        x=x,
        y=y,
        showarrow=False,
        xref='paper',
        yref='paper',
        text=str(text),
        font=dict(
            size=13,
            color=AXIS_TITLE_COLOR
        )
    )
    return annotation_dict

def _axis_title_annotation(text, x_or_y_axis='x'):
    """Customize annotation layer pertaining to an overall axis title."""
    x_pos = 0.5
    y_pos = -0.3
    textangle = 0
    name = "x-title"    # Used to easily identify annotation for later modification
    if x_or_y_axis == 'y':
        x_pos = -0.15
        y_pos = 0.5
        textangle = 270
        name = "y-title"

    if not text:
        text = ''

    annot = {'font': {'color': '#000000', 'size': AXIS_TITLE_SIZE},
             'showarrow': False,
             'text': text,
             'textangle': textangle,
             'name':name,
             'x': x_pos,
             'xref': 'paper',
             'y': y_pos,
             'yref': 'paper'}
    return annot


# TODO: Various issues moved to https://github.com/jorvis/gEAR/issues/784
def _gear_facet_grid(df, x, y, facet_row, facet_col,
                color_name, colormapping, color_type, num_of_rows,
                num_of_cols, facet_row_labels, facet_col_labels,
                trace_type, flipped_rows, flipped_cols,
                SUBPLOT_SPACING, marker_color, text_name, jitter, kwargs_trace, kwargs_marker):

    print("DEBUG check: within _facet_grid, trace type is {}".format(trace_type), file=sys.stderr)

    fig = make_subplots(rows=num_of_rows, cols=num_of_cols,
                        shared_xaxes=True, shared_yaxes=True,
                        horizontal_spacing=SUBPLOT_SPACING,
                        vertical_spacing=SUBPLOT_SPACING, print_grid=False)

    if PLOT_LOGGING:
        print("DEBUG: facet_row:{0} facet_col:{1}".format(facet_row, facet_col), file=sys.stderr)
    else:
        print("DEBUG: plot debugging appears to be off", file=sys.stderr)

    # 'color_name' will be omitted if 'None'
    # Do not add color_name for 'numerical' data since that will blow up the number of traces
    priority_groups = _build_priority_groups(facet_row, facet_col, None, x) if color_type == "numerical" \
        else _build_priority_groups(facet_row, facet_col, color_name, x)

    # If replicates are present, use mean and stdev of expression data as datapoints
    if 'replicate' in df.columns and trace_type != 'violin':
        grouped = df.groupby(priority_groups)

        if color_type == "numerical":
            df = grouped.agg({
                    color_name: ['mean'],
                    "raw_value": ['mean', 'std']
                }) \
                .dropna() \
                .reset_index()
        else:
            df = grouped.agg({
                    "raw_value": ['mean', 'std']
                }) \
                .dropna() \
                .loc[:, 'raw_value'] \
                .rename(columns=dict(mean='raw_value')) \
                .reset_index()

    # When doing groupby for iteration purposes, we do not want the "x" column grouped
    # Violins will keep x, so they will be processed in the 'groupby' block to create multiple traces, instead of one
    if trace_type != "violin" or len(priority_groups) > 1:
        priority_groups.pop()   # Assumes 'x' group is always the last item.

    # Map indexes for subplot ordering.  Indexes start at 1 since plotting rows/cols start at 1
    facet_row_groups = list(df.groupby(facet_row)) if facet_row else []
    facet_row_indexes = {group[0]: idx for idx, group in enumerate(facet_row_groups, start=1)}
    facet_col_groups = list(df.groupby(facet_col)) if facet_col else []
    facet_col_indexes = {group[0]: idx for idx, group in enumerate(facet_col_groups, start=1)}

    traces = []
    row_idxs = []
    col_idxs = []
    annotations = []
    names_in_legend = {}

    # Update kwargs_marker with default values, but higher-level assignments should supercede
    if color_type == None:
        kwargs_marker.setdefault('color', marker_color)

    # If there is nothing to group by there will only be one trace
    # Essentially no facet rows or columns, and no replicates
    if not len(priority_groups):
        trace = _create_trace(df, x, y, trace_type, kwargs_trace, kwargs_marker)

        # 'categorical' colortypes will always have the 'color_name' priority group in the list
        if color_type == "numerical":
            trace = _append_to_numerical_trace(trace, df, color_name, x, y, trace_type)

        # If applicable, add some jitter to the datapoints.
        # Attempting to replicate the "strip" plot function in plotly express
        if trace_type == "scatter" and _is_categorical(df[x]) and jitter:
            trace = _adjust_plot_for_jitter(trace, jitter, None)

        if text_name and text_name in df:
            trace['text'] = df[text_name]
        traces.append(trace)
        row_idxs.append(1)  # Do not need, but is a safeguard in case annotations are added
        col_idxs.append(1)
    else:
        # https://pandas.pydata.org/docs/user_guide/groupby.html#iterating-through-groups
        # Worth noting.  If number of groups = 1 then 'name' is a string, else is a tuple
        for name, group in df.groupby(priority_groups):
            trace = _create_trace(group, x, y, trace_type, kwargs_trace, kwargs_marker)

            curr_color = None
            if color_type == "categorical":
                # If name is a tuple, color_name is last element
                curr_color = name
                if isinstance(name, tuple):
                    curr_color = name[-1]
                # Making the assumption that kwargs_marker should not overwrite these marker attributes
                trace['marker']['color'] = colormapping[curr_color]
                trace['name'] = str(curr_color)

            elif color_type == "numerical":
                trace = _append_to_numerical_trace(trace, group, color_name, x, y, trace_type)

            # If applicable, add some jitter to the datapoints.
            # Attempting to replicate the "strip" plot function in plotly express
            if trace_type == "scatter" and _is_categorical(group[x]) and jitter:
                trace = _adjust_plot_for_jitter(trace, jitter, curr_color)

            if not facet_row and not facet_col:
                if text_name and text_name in group:
                    trace['text'] = group[text_name]

            # Plotly workaround:
            # Once a trace name has been added to the legend,
            # we don't want to include future traces in the
            # legend with that name, otherwise we get legend
            # label explosion. (mostly with 'categorical' color_type)
            legend_key = name
            if color_type == "categorical":
                legend_key = curr_color

            if legend_key in names_in_legend:
                trace['showlegend'] = False
            else:
                names_in_legend[legend_key] = True

            traces.append(trace)

            # Now determine which plot this trace should go to.  Facet column is first if row does not exist.
            if isinstance(name, tuple):
                row_idxs.append(facet_row_indexes[name[0]] if facet_row else 1)
                if facet_row:
                    col_idxs.append(facet_col_indexes[name[1]] if facet_col else 1)
                else:
                    col_idxs.append(facet_col_indexes[name[0]] if facet_col else 1)
            else:
                row_idxs.append(facet_row_indexes[name] if facet_row else 1)
                col_idxs.append(facet_col_indexes[name] if facet_col else 1)

    # Create annotations for each facet row or column
    for rowname in facet_row_indexes:
        label = _return_label(rowname, facet_row_labels, facet_row)
        annotations.append(
                _annotation_dict(
                    label,
                    (num_of_rows - facet_row_indexes[rowname]) + 1, # Order from top to bottom
                    num_of_rows,    # Adds to right of plot
                    SUBPLOT_SPACING,
                    'row',
                    flipped_rows)
            )
    for colname in facet_col_indexes:
        label = _return_label(colname, facet_col_labels, facet_col)
        annotations.append(
                _annotation_dict(
                    label,
                    facet_col_indexes[colname],
                    num_of_cols,    # Adds to bottom of plot
                    SUBPLOT_SPACING,
                    'col',
                    flipped_cols)
            )

    # Traces, row, and col lists should be 1-to-1-to-1
    fig.add_traces(traces, rows=row_idxs, cols=col_idxs)
    return fig, annotations

def create_facet_grid(df, x=None, y=None, facet_row=None, facet_col=None,
                      color_name=None, colormap=None, color_is_cat=False,
                      facet_row_labels=None, facet_col_labels=None,
                      height=None, width=None, trace_type='scatter', hide_x_labels=False,hide_y_labels=False,
                      scales='fixed', dtick_x=None, dtick_y=None, text_name=None,
                      show_boxes=True, ggplot2=False, binsize=1, jitter=0, **kwargs):
    """
    Returns figure for facet grid.
    :param (pd.DataFrame) df: the dataframe of columns for the facet grid.
    :param (str) x: the name of the dataframe column for the x axis data.
    :param (str) y: the name of the dataframe column for the y axis data.
    :param (str) facet_row: the name of the dataframe column that is used to
        facet the grid into row panels.
    :param (str) facet_col: the name of the dataframe column that is used to
        facet the grid into column panels.
    :param (str) color_name: the name of your dataframe column that will
        function as the colormap variable.
    :param (str|list|dict) colormap: the param that determines how the
        color_name column colors the data. If the dataframe contains numeric
        data, then a dictionary of colors will group the data categorically
        while a Plotly Colorscale name or a custom colorscale will treat it
        numerically. To learn more about colors and types of colormap, run
        `help(plotly.colors)`.
    :param (bool) color_is_cat: determines whether a numerical column for the
        colormap will be treated as categorical (True) or sequential (False).
            Default = False.
    :param (str|dict) facet_row_labels: set to either 'name' or a dictionary
        of all the unique values in the faceting row mapped to some text to
        show up in the label annotations. If None, labeling works like usual.
    :param (str|dict) facet_col_labels: set to either 'name' or a dictionary
        of all the values in the faceting row mapped to some text to show up
        in the label annotations. If None, labeling works like usual.
    :param (int) height: the height of the facet grid figure.
    :param (int) width: the width of the facet grid figure.
    :param (str) trace_type: decides the type of plot to appear in the
        facet grid. The options are 'scatter', 'scattergl', 'histogram',
        'bar', and 'box'.
        Default = 'scatter'.
    :param (str) scales: determines if axes have fixed ranges or not. Valid
        settings are 'fixed' (all axes fixed), 'free_x' (x axis free only),
        'free_y' (y axis free only) or 'free' (both axes free).
    :param (float) dtick_x: determines the distance between each tick on the
        x-axis. Default is None which means dtick_x is set automatically.
    :param (float) dtick_y: determines the distance between each tick on the
        y-axis. Default is None which means dtick_y is set automatically.
    :param (bool) show_boxes: draws grey boxes behind the facet titles.
    :param (bool) ggplot2: draws the facet grid in the style of `ggplot2`. See
        http://ggplot2.tidyverse.org/reference/facet_grid.html for reference.
        Default = False
    :param (int) binsize: groups all data into bins of a given length.
    :param (int) jitter: Amount to offset an individual categorical x-axis
        datapoint.  The higher the number, the more extreme the jitter
        Default: 0 (no jitter)
    :param (dict) kwargs: a dictionary of scatterplot arguments.
    Examples 1: One Way Faceting
    ```
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import pandas as pd
    mpg = pd.read_table('https://raw.githubusercontent.com/plotly/datasets/master/mpg_2017.txt')
    fig = ff.create_facet_grid(
        mpg,
        x='displ',
        y='cty',
        facet_col='cyl',
    )
    py.iplot(fig, filename='facet_grid_mpg_one_way_facet')
    ```
    Example 2: Two Way Faceting
    ```
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import pandas as pd
    mpg = pd.read_table('https://raw.githubusercontent.com/plotly/datasets/master/mpg_2017.txt')
    fig = ff.create_facet_grid(
        mpg,
        x='displ',
        y='cty',
        facet_row='drv',
        facet_col='cyl',
    )
    py.iplot(fig, filename='facet_grid_mpg_two_way_facet')
    ```
    Example 3: Categorical Coloring
    ```
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import pandas as pd
    mpg = pd.read_table('https://raw.githubusercontent.com/plotly/datasets/master/mpg_2017.txt')
    fig = ff.create_facet_grid(
        mtcars,
        x='mpg',
        y='wt',
        facet_col='cyl',
        color_name='cyl',
        color_is_cat=True,
    )
    py.iplot(fig, filename='facet_grid_mpg_default_colors')
    ```
    Example 4: Sequential Coloring
    ```
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import pandas as pd
    tips = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/tips.csv')
    fig = ff.create_facet_grid(
        tips,
        x='total_bill',
        y='tip',
        facet_row='sex',
        facet_col='smoker',
        color_name='size',
        colormap='Viridis',
    )
    py.iplot(fig, filename='facet_grid_tips_sequential_colors')
    ```
    Example 5: Custom labels
    ```
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import pandas as pd
    mtcars = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/mtcars.csv')
    fig = ff.create_facet_grid(
        mtcars,
        x='wt',
        y='mpg',
        facet_col='cyl',
        facet_col_labels={4: "$\\alpha$", 6: '$\\beta$', 8: '$\sqrt[y]{x}$'},
    )
    py.iplot(fig, filename='facet_grid_mtcars_custom_labels')
    ```
    Example 6: Other Trace Type
    ```
    import plotly.plotly as py
    import plotly.figure_factory as ff
    import pandas as pd
    mtcars = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/mtcars.csv')
    fig = ff.create_facet_grid(
        mtcars,
        x='wt',
        facet_col='cyl',
        trace_type='histogram',
    )
    py.iplot(fig, filename='facet_grid_mtcars_other_trace_type')
    ```
    """

    if not pd:
        raise exceptions.ImportError(
            "'pandas' must be installed for this figure_factory."
        )

    if not isinstance(df, pd.DataFrame):
        raise exceptions.PlotlyError(
            "You must input a pandas DataFrame."
        )

    # make sure all columns are of homogenous datatype
    utils.validate_dataframe(df)

    # the tsne_dynamic trace type is an alias for scatter
    if trace_type == 'tsne_dynamic':
        trace_type = 'scatter'

    if PLOT_LOGGING:
        print("DEBUG: trace_type is: {0}".format(trace_type), file=sys.stderr)

    if trace_type in ['scatter', 'scattergl']:
        if not x or not y:
            raise exceptions.PlotlyError(
                "You need to input 'x' and 'y' if you are you are using a "
                "trace_type of 'scatter' or 'scattergl'."
            )

    for key in [x, y, facet_row, facet_col, color_name]:
        if key is not None:
            try:
                df[key]
            except KeyError:
                raise exceptions.PlotlyError(
                    "x, y, facet_row, facet_col and color_name must be keys "
                    "in your dataframe."
                )

    if trace_type not in VALID_TRACE_TYPES:
        raise exceptions.PlotlyError(
            "'trace_type' must be in {}".format(VALID_TRACE_TYPES)
        )

    if trace_type == 'histogram' or trace_type == 'line':
        SUBPLOT_SPACING = 0.06
    else:
        SUBPLOT_SPACING = 0.015

    # seperate kwargs for marker and else
    if 'marker' in kwargs:
        kwargs_marker = kwargs['marker']
    else:
        kwargs_marker = {}
    marker_color = kwargs_marker.pop('color', None)
    kwargs.pop('marker', None)
    kwargs_trace = kwargs

    if 'size' not in kwargs_marker:
        kwargs_marker['size'] = 3

    # Bar plots do not accept size markers
    if trace_type == 'bar':
        kwargs_marker.pop('size', None)

    if 'opacity' not in kwargs_marker:
        kwargs_trace['opacity'] = 0.6

    # if 'line' not in kwargs_marker:
        # kwargs_marker['line'] = {'color': 'darkgrey', 'width': 1}

    # default marker size
    if not ggplot2:
        if not marker_color:
            marker_color = '#401362'
    else:
        marker_color = 'rgb(0, 0, 0)'

    num_of_rows = 1
    num_of_cols = 1
    flipped_rows = False
    flipped_cols = False
    if facet_row:
        num_of_rows = len(df[facet_row].unique())
        flipped_rows = _is_flipped(num_of_rows)
        if isinstance(facet_row_labels, dict):
            for key in df[facet_row].unique():
                if key not in facet_row_labels.keys():
                    unique_keys = df[facet_row].unique().tolist()
                    raise exceptions.PlotlyError(
                        CUSTOM_LABEL_ERROR.format(unique_keys)
                    )
    if facet_col:
        num_of_cols = len(df[facet_col].unique())
        flipped_cols = _is_flipped(num_of_cols)
        if isinstance(facet_col_labels, dict):
            for key in df[facet_col].unique():
                if key not in facet_col_labels.keys():
                    unique_keys = df[facet_col].unique().tolist()
                    raise exceptions.PlotlyError(
                        CUSTOM_LABEL_ERROR.format(unique_keys)
                    )

    # Set up some args to pass to _gear_facet_grid function
    show_legend = False
    colormapping = None
    color_type = None   # None, 'categorical', or 'numerical'

    # If there is a color label, use either the categorial or numerical facet grid
    if color_name:
        if isinstance(colormap, dict):
            show_legend = True
            color_type = "categorical"

            clrs.validate_colors_dict(colormap, 'rgb')

            for val in df[color_name].unique():
                if val not in colormap.keys():
                    raise exceptions.PlotlyError(
                        "If using 'colormap' as a dictionary, make sure "
                        "all the values of the colormap column are in "
                        "the keys of your dictionary."
                    )

            colormapping = colormap
            if PLOT_LOGGING:
                print("DEBUG: Color type is 'categorical' with colormap dict", file=sys.stderr)

        elif isinstance(colormap, list):
            color_type = "numerical"
            colormapping = colormap
            clrs.validate_colorscale(colormapping)
            if PLOT_LOGGING:
                print("DEBUG: Color type is 'numerical' from colormap list", file=sys.stderr)
        elif isinstance(colormap, str):
            color_type = "numerical"
            if colormap in clrs.PLOTLY_SCALES.keys():
                colormapping = clrs.PLOTLY_SCALES[colormap]
            else:
                raise exceptions.PlotlyError(
                    "If 'colormap' is a string, it must be the name "
                    "of a Plotly Colorscale. The available colorscale "
                    "names are {}".format(clrs.PLOTLY_SCALES.keys())
                )
            if PLOT_LOGGING:
                print("DEBUG: Color type is 'numerical' from colormap string", file=sys.stderr)
        else:
            if isinstance(df[color_name].iloc[0], str) or color_is_cat:
                color_type = "categorical"
                # use default plotly colors for dictionary
                default_colors = clrs.DEFAULT_PLOTLY_COLORS
                colormap = {}
                j = 0
                for val in df[color_name].unique():
                    if j >= len(default_colors):
                        j = 0
                    colormap[val] = default_colors[j]
                    j += 1
                colormapping = colormap
            else:
                color_type = "numerical"
                colormapping = [
                                [0, 'rgb(218, 183, 193)'],
                                [0.35, 'rgb(194, 137, 166)'],
                                [0.5, 'rgb(169, 98, 151)'],
                                [0.6, 'rgb(145, 66, 143)'],
                                [0.7, 'rgb(105, 39, 122)'],
                                [1, 'rgb(63, 19, 98)']
                            ]
                if PLOT_LOGGING:
                    print("DEBUG: Color type is 'numerical' with no colormap", file=sys.stderr)
    else:
        if PLOT_LOGGING:
            print("DEBUG: Color type is 'None'", file=sys.stderr)

    fig, annotations = _gear_facet_grid(
        df, x, y, facet_row, facet_col, color_name, colormapping, color_type,
        num_of_rows, num_of_cols, facet_row_labels, facet_col_labels,
        trace_type, flipped_rows, flipped_cols,
        SUBPLOT_SPACING, marker_color, text_name, jitter, kwargs_trace, kwargs_marker
    )

    ### General layout adjustments
    fig['layout'].update(title='', paper_bgcolor=PAPER_BGCOLOR)
    fig['layout']['hovermode'] = "closest"
    # Default "plotly" theme produces gray plot backgrounds
    fig['layout']['template'] = "none"

    # axis titles
    x_title_annot = _axis_title_annotation('', 'x')
    y_title_annot = _axis_title_annotation('', 'y')

    # annotations
    annotations.append(x_title_annot)
    annotations.append(y_title_annot)


    # all xaxis and yaxis labels
    axis_labels = {'x': [], 'y': []}
    for key in fig['layout']:
        if 'xaxis' in key:
            axis_labels['x'].append(key)
        elif 'yaxis' in key:
            axis_labels['y'].append(key)

    string_number_in_data = False
    for var in [v for v in [x, y] if v]:
        if isinstance(df[var].tolist()[0], str):
            for item in df[var]:
                try:
                    int(item)
                    string_number_in_data = True
                except ValueError:
                    pass

    # Iterated through 'x' or 'y' axis
    for x_y in axis_labels.keys():
        # Iterate through all faceted axes
        for axis_name in axis_labels[x_y]:
            # Common to both x and y
            if string_number_in_data:
                fig['layout'][axis_name]['type'] = 'category'
            fig['layout'][axis_name]['showgrid'] = False
            fig['layout'][axis_name]['automargin'] = True
            fig['layout'][axis_name]['zeroline'] = False
            # Specific axis only
            if x_y == 'x':
                if hide_x_labels:
                    #TODO: test with 'visible' attribute instead of 'showticklabels'
                    fig['layout'][axis_name]['showticklabels'] = False

                # Uniformity of tick angles if facet groupings are present
                if facet_col:
                    fig['layout'][axis_name]['tickangle'] = 270
            elif x_y == 'y':
                fig['layout'][axis_name]['hoverformat'] = '.2f'
                if hide_y_labels:
                    fig['layout'][axis_name]['showticklabels'] = False

    fig['layout']['autosize'] = True

    # legend
    fig['layout']['showlegend'] = show_legend
    fig['layout']['legend']['bgcolor'] = LEGEND_COLOR
    fig['layout']['legend']['borderwidth'] = LEGEND_BORDER_WIDTH
    fig['layout']['legend']['x'] = 1.05
    fig['layout']['legend']['y'] = 1
    fig['layout']['legend']['yanchor'] = 'top'

    # Colorbar adjustments
    if color_type == "numerical":
        fig['layout']['coloraxis'] = {
            "colorscale":colormapping,   # Defines the range of colors for a numerical color group
            "colorbar": {'x':1.15},
            "showscale": True,

        }


    # Violin plot settings
    if trace_type == 'violin':
        if color_name:
            fig['layout']['violinmode'] = 'group'
        else:
            fig['layout']['violinmode'] = 'overlay'

    # assign annotations to figure
    fig['layout']['annotations'] = annotations


    # autoscale histogram bars
    if trace_type not in ['scatter', 'line', 'scattergl']:
        scales = 'free'

    # validate scales
    if scales not in ['fixed', 'free_x', 'free_y', 'free']:
        raise exceptions.PlotlyError(
            "'scales' must be set to 'fixed', 'free_x', 'free_y' and 'free'."
        )
    fixed_axes = None
    if scales == 'fixed':
        fixed_axes = ['x', 'y']
    elif scales == 'free_x':
        fixed_axes = ['y']
    elif scales == 'free_y':
        fixed_axes = ['x']
    elif scales == 'free':
        fixed_axes = []
    else:
        raise("Invalid scale type provided.  Must be 'fixed', 'free_x', 'free_y', or 'free'")

    # SAdkins - Removed checks for None and length and sparse matrix check
    # since recent edits should have all traces populated with data
    if len(fig['data']):
        # fixed ranges
        for x_y in fixed_axes:
            min_range = min(chain( *(trace[x_y] for trace in fig['data']) ))
            max_range = max(chain( *(trace[x_y] for trace in fig['data']) ))
            range_are_numbers = (isinstance(min_range, Number)
                                and isinstance(max_range, Number))
            if PLOT_LOGGING:
                print("DEBUG: On axis:{0} min_range:{1} max_range:{2}".format(x_y, min_range, max_range), file=sys.stderr)

            user_dtick = None
            if x_y == 'x':
                user_dtick = dtick_x
            elif x_y == 'y':
                user_dtick = dtick_y

            dtick, min_range, max_range = _calculate_dtick(min_range, max_range, range_are_numbers, user_dtick)

            # For the given axis dimension set tick attributes
            for axis_title in axis_labels[x_y]:
                fig['layout'][axis_title]['dtick'] = dtick
                fig['layout'][axis_title]['ticklen'] = 0
                if range_are_numbers:
                    fig['layout'][axis_title]['range'] = [min_range, max_range]

    else:
        if PLOT_LOGGING:
            print("DEBUG: No trace data for current plot")

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
            "select2d",
            "lasso2d",
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

def rgb_to_hex(r,g,b):
    hex = "#{:02x}{:02x}{:02x}".format(int(r),int(g),int(b))
    return hex

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
