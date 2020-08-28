from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from collections.abc import Iterable
from typing import List, Tuple, Mapping, Any
from pandas import DataFrame, Series
from plotly import tools
# from plotly.colors import DEFAULT_PLOTLY_COLORS
from plotly.colors import n_colors, unlabel_rgb, hex_to_rgb
from plotly.utils import PlotlyJSONEncoder
from itertools import cycle
import pandas as pd
import numpy as np
import plotly.offline as py
import plotly.graph_objs as go
import json
import sys

def rgb_to_hex(r,g,b):
    hex = "#{:02x}{:02x}{:02x}".format(int(r),int(g),int(b))
    return hex

@dataclass
class Plot(ABC):
    """Abstract Base Plot Class
    """
    data: 'PlotData'
    colors: dict
    @abstractmethod
    def render_go(self, **kwargs) -> Any:
        """Render bar graph object

        Parameters
        ----------
        kwargs
            Keyword arguments specific to child render_go

        Returns
        -------
        Plotly Bar Graph Obect

        Raises
        ------
        NotImplementedError
            This abstract method is meant for child class' to override.
        """
        raise NotImplementedError

    @abstractmethod
    def get_trace_data(self) -> List[Mapping]:
        """Generate trace data for rendering Plotly graph object.

        Returns
        -------
        Return list of of dictionaries that represent arguments
        for rendering Plotly graph object.

        Raises
        ------
        NotImplementedError
            This abstract method is meant for child class' to override.
        """
        raise NotImplementedError

    def make_figure(self, traces):
        """Create plot figure.

        This base method will be called for data that has 2 or less
        indices. If the class' data has a hierarchical index of three indices,
        we create a subplot figure in order to properly display groups.
        Different plot types have a unique way of rendering its subplots and
        this method is overridden.

        Parameters
        ----------
        traces
            List of Plotly graph objects (traces)

        Returns
        -------
        Graph Object Figure
        """
        return go.Figure(data=traces)

    def render(self) -> str:
        """Render plot to JSON.

        Returns
        -------
        JSON representation of plot figure.
        """
        trace_data = self.get_trace_data()
        traces = [self.render_go(**go_args) for go_args in trace_data]
        fig = self.make_figure(traces)
        fig.layout.update(self.get_layout())

        return json.dumps(fig, cls=PlotlyJSONEncoder)

    def get_layout(self) -> Mapping:
        """Get layout option for plot figure.

        Returns
        -------
        Dictionary mapping appropriate layout options.
        """
        return go.Layout(
            # Top margin is set to 20 to allow space
            # between visualization and control bar
            margin=dict(t=50, r=50, b=50, l=50),
            yaxis = dict(hoverformat = '.2f')
        )

    @staticmethod
    def _plotly_color_map(names: List[str]) -> Mapping[str, str]:
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
            linear_purple_scale = n_colors(
                lowcolor='rgb(218,183,193)',
                highcolor='rgb(64,19,98)',
                n_colors=len(names),
                colortype='rgb'
            )

            hex_scale = [rgb_to_hex(*unlabel_rgb(rgb)) for rgb in linear_purple_scale]
            # From https://github.com/plotly/plotly.py/issues/1026
            color_generator = cycle(hex_scale)
            return dict(zip(names, color_generator))
        else:
            return dict(zip(names, ['rgb(64,19,98)']))


@dataclass
class ScatterPlot(Plot):
    """Scatter Plot"""
    def render_go(self,
        x: List[float],
        y: List[float],
        name: str,
        marker: Mapping
    ) -> go.Scatter:
        """Render scatter plot."""
        return go.Scatter(
            x=x,
            y=y,
            name=name,
            marker=marker
        )
@dataclass
class BarPlot(Plot):
    """Bar Plot"""
    def render_go(self,
        x: List[str],
        y: List[float],
        name: str,
        marker: Mapping,
        error_y: Mapping) -> go.Bar:
        """Render bar graph object

        Parameters
        ----------
        x: List[str]
            Data group levels, e.g., ['Cochlea'].
        y: List[float]
            Raw values from data.
        name: str
            Name of particular subgroup, e.g., 'Cochlea'
        error_y: Mapping
            Error values for each bar for error bars.

        Returns
        -------
        Plotly Bar Graph Obect
        """
        return go.Bar(
            x=x if isinstance(x, Iterable) else [x],
            y=y if isinstance(y, Iterable) else [y],
            name=name,
            marker=marker,
            error_y=error_y)

    def get_trace_data(self) -> List[Mapping]:
        """Generate trace data for rendering Plotly graph object.

        Returns
        -------
        Return list of of dictionaries that represent arguments
        for rendering Plotly graph object.
        """
        # indices example
        # [('Cochlea', 'P0'), ('Cochlea', 'P5'),
        # ('Cochlea', 'P1'), ('Cochlea', 'P6')]
        indices = self.data.get_unique_idx()

        # trace_data example
        # [{
        #   x: ['DMSO_GFP+'+, 'DMSO_GFP-', 'DAPT_GFP+', DAPT_GFP-'],
        #   y: [224.649719, 158.241577, 222.244827, 166.381897]
        #   name: P0, <- This is the name that appears on the legend
        #   error_y: {
        #       ...
        #       array=[4.981544, 10.478561, 1.40493, 3.360636]
        #   }
        # }, ...]
        if isinstance(self.data, DataWithOneIndex):
            # If data is indexed by a single column, there is
            # no hierarchical index that we need to worry about and
            # can just pass index and raw value as x and y
            trace_data = ([
                dict(
                    x=self.data.values.index,
                    y=self.data.values.raw_value,
                    name=None,
                    # TODO COLOR MAP HERE FOR MARKER
                    marker=dict(
                        color='rgb(64,19,98)'
                    ),
                    error_y=dict(
                        type='data',
                        # have errors bars only if there are replicates
                        # in expression data frame
                        array=(
                            self.data.stats.raw_value
                        ) if self.data.has_replicates and self.data.summarize and not self.data.stats.empty else list()
                    )
                )
            ])
        else:
            if isinstance(self.data, DataWithThreeIndices):
                names = [idx[-1] for idx in indices]
            else:
                names = indices

            if not self.colors:
                self.colors = Plot._plotly_color_map(names)

            trace_data = ([
                dict(
                    x=self.data.values.loc[idx].index,
                    y=self.data.values.loc[idx].raw_value,
                    name=name,
                    marker=dict(
                        color=self.colors.get(name)
                    ),
                    error_y=dict(
                        type='data',
                        # have errors bars only if there are replicates
                        # in expression data frame
                        array=(
                            self.data.stats.loc[idx].raw_value
                        ) if self.data.has_replicates and self.data.summarize and not self.data.stats.empty else list()
                    )
                ) for name, idx in zip(names,indices)
            ])
        return trace_data

    def make_figure(self, traces: List[Any]) -> Any:
        """Create plot figure.

        If the class' data has a hierarchical index of three indices,
        we create a subplot figure in order to properly display groups.

        Parameters
        ----------
        traces
            List of Plotly graph objects (traces)

        Returns
        -------
        Graph Object Figure
        """
        if not isinstance(self.data, DataWithThreeIndices):
            return super().make_figure(traces)
        else:
            # The group with the least amount of levels
            # is the group we will create subplots for.
            GROUP_WITH_LEAST_LEVELS = 0
            subgroup_col = (
                self.data
                    .values
                    .index
                    .names[GROUP_WITH_LEAST_LEVELS]
            )

            subgroup_col_levels = (
                self.data
                    .values
                    .index
                    .get_level_values(subgroup_col)
                    .unique()
                    .tolist()
            )

            num_cols = len(subgroup_col_levels)
            ROW_NUM = 1 # Bar plot subgroups we display across
                        # the x-axis, so we need only one row
            fig = (
                tools.make_subplots(
                    rows=ROW_NUM,
                    cols=num_cols,
                    print_grid=False,
                    subplot_titles=subgroup_col_levels,
                    shared_yaxes=True
                )
            )

            # Column map allows us to map a trace
            # to the correct subplot position
            col_map = {
                level:i
                for (i, level) in enumerate(subgroup_col_levels, 1)
            }

            names = [trace.name for trace in traces]

            if not self.colors:
                self.colors.update(**Plot._plotly_color_map(names))

            names_in_legend = { name:False for name in names }

            for (subgroup, name), trace in zip(self.data.get_unique_idx(), traces):
                trace.marker.color = self.colors.get(name)

                # Plotly workaround:
                # Once a trace name has been added to the legend,
                # we don't want to include future traces in the
                # legend with that name
                if name in names_in_legend and names_in_legend[name]:
                    trace.showlegend = False
                else:
                    names_in_legend[name] = True

                fig.append_trace(trace, ROW_NUM, col_map[subgroup])

            return fig

@dataclass
class ViolinPlot(BarPlot):
    """Violin Plot"""

    def render_go(
        self,
        x: List[str],
        y: List[float],
        name: str,
        #error_y,
        marker: Mapping,
        ) -> go.Violin:
        """Render Violin graph object

        Parameters
        ----------
        x: List[str]
            x values, e.g., ['Cochlea', 'Cochlea', 'Cochlea'].
        y: List[float]
            Raw values from data.
        name: str
            Name of particular subgroup, e.g., 'Cochlea'.
        error_y: Mapping
            Error values for each bar for error bars.

        Returns
        -------
        Plotly Bar Graph Object.
        """
        return go.Violin(
            x=x if isinstance(x, Iterable) else [x],
            y=y if isinstance(y, Iterable) else [y],
            name=name,
            marker=marker,
            scalemode='count',
            showlegend=False if isinstance(self.data, DataWithOneIndex) else True,
            # box=dict(visible=True),
            # meanline=dict(visible=True)
        )
    def get_layout(self) -> Mapping:
        """Get layout option for plot figure.

        Returns
        -------
        Dictionary mapping appropriate layout options.
        """

        return (
            go.Layout(
                violinmode="overlay" if isinstance(self.data, DataWithOneIndex) else "group",
                yaxis=dict(hoverformat = '.2f')
            )
        )

    def get_trace_data(self) -> List[Mapping]:
        """Generate trace data for rendering Plotly graph object.

        Returns
        -------
        Return list of of dictionaries that represent arguments
        for rendering Plotly graph object.
        """
        # indices example
        # [('Cochlea', 'P0'), ('Cochlea', 'P5'),
        # ('Cochlea', 'P1'), ('Cochlea', 'P6')]
        indices = self.data.get_unique_idx()

        # trace_data example
        # [{
        #   x: ['DMSO_GFP+'+, 'DMSO_GFP-', 'DAPT_GFP+', DAPT_GFP-'],
        #   y: [224.649719, 158.241577, 222.244827, 166.381897]
        #   name: P0, <- This is the name that appears on the legend
        #   error_y: {
        #       ...
        #       array=[4.981544, 10.478561, 1.40493, 3.360636]
        #   }
        # }, ...]
        if isinstance(self.data, DataWithOneIndex):
            # If data is indexed by a single column, there is
            # no hierarchical index that we need to worry about and
            # can just pass index and raw value as x and y
            # df = self.data.values
            # if not self.colors:
            #     self.colors = Plot._plotly_color_map(indices)

            trace_data = []
            for idx in indices:
                i = self.data.values.loc[idx].index
                trace = {
                        "x": self.data.values.loc[idx].index,
                        "y": self.data.values.loc[idx]['raw_value'],
                        "name": None, #iidx,
                        "marker": dict(
                            # color=self.colors.get(idx),
                            color='rgb(64,19,98)'
                        ),
                        #"points": "all",
                        #"pointpos": 0,
                        # "jitter": 1,
                        # "box": {
                        #     "visible": True
                        # },
                        # "meanline": {
                        #     "visible": True
                        # }
                    }
                trace_data.append(trace)
        #     trace_data = ([
        #       dict(
        #           x=self.data.values.index,
        #           y=self.data.values.raw_value,
        #           name=None,
        #           marker=dict(
        #             color='rgb(64,19,98)'
        #           ),
        #       )
        #   ])
        else:
            if isinstance(self.data, DataWithThreeIndices):
                names = [idx[-1] for idx in indices]
            else:
                names = indices

            if not self.colors:
                self.colors = Plot._plotly_color_map(names)

            trace_data = ([
                dict(
                    x=self.data.values.loc[idx].index,
                    y=self.data.values.loc[idx].raw_value,
                    name=name,
                    marker=dict(
                        color=self.colors.get(name)
                    )

                ) for name, idx in zip(names,indices)
            ])
        return trace_data

    # def get_trace_data(self) -> List[Mapping]:
    #     """Generate trace data for rendering Plotly graph object.

    #     Note:
    #         Assumes 'cell_type' column exists and is the first
    #         level index on the data frame.

    #     Returns
    #     -------
    #     Return list of of dictionaries that represent arguments
    #     for rendering Plotly graph object.
    #     """
    #     first_index = self.data.values.index.names[0]
    #     levels = (
    #         self.data
    #             .values
    #             .index
    #             .get_level_values(first_index)
    #             .unique()
    #     )

    #     if not self.colors:
    #         self.colors = Plot._plotly_color_map(levels)

    #     trace_data = list()
    #     for level in levels:
    #         y = (
    #             self.data
    #                 .values
    #                 .loc[level]
    #                 .raw_value
    #         )
    #         # Repeat cell type for every y value
    #         # for violin plot to properly draw
    #         # distribution
    #         x = np.repeat(level, y.size)

    #         trace_data.append(
    #             dict(
    #                 x=x,
    #                 y=y,
    #                 name=level,
    #                 marker=dict(color=self.colors.get(level))
    #             )
    #         )

    #     return trace_data


@dataclass
class LinePlot(BarPlot):
    """Line Plot

    TODO:
        Currently inherits from BarPlot, which intuitively
        doesn't make sense. Inheriting from BarPlot allows
        us to follow DRY principles and not repeat similar
        methods, but should revisit this.
    """

    def __post_init__(self):
        """TODO:
            Refactor. We are rewrapping the data and aggregating
            on hard-coded columns specific to the line plot. Possibly have an index
            argument when creating plots vs. making assumptions on what to index
            on. This will also be fixed once curation view is created and settings
            will be passed in.
        """
        # self.data.index.remove('condition')
        # self.data.index.remove('time_point')

        # self.data = PlotDataFactory.wrap_data(
        #     self.data.data,
        #     ['condition', *self.data.index, 'time_point']
        # )

    def make_figure(self, traces):
        """Create plot figure.

        If the class' data has a hierarchical index of three indices,
        we create a subplot figure in order to properly display groups.

        Parameters
        ----------
        traces
            List of Plotly graph objects (traces)

        Returns
        -------
        Graph Object Figure
        """
        if not isinstance(self.data, DataWithThreeIndices):
            return super().make_figure(traces)
        else:
            # The group with the least amount of levels
            # is the group we will create subplots for.
            GROUP_WITH_LEAST_LEVELS = 0
            subgroup_col = (
                self.data
                    .values
                    .index
                    .names[GROUP_WITH_LEAST_LEVELS]
            )

            subgroup_col_levels = (
                self.data
                    .values
                    .index
                    .get_level_values(subgroup_col)
                    .unique()
                    .tolist()
            )


            num_rows = len(subgroup_col_levels)
            COL_NUM = 1 # Line plot subgroups we display across
                        # the y-axis, so we need only one column

            fig = (
                tools.make_subplots(
                    rows=num_rows,
                    cols=COL_NUM,
                    print_grid=False,
                    subplot_titles=subgroup_col_levels
                )
            )
            # Row map allows us to map a trace
            # to the correct subplot position
            row_map = {
                level:i
                for (i, level) in enumerate(subgroup_col_levels, 1)
            }

            names = [trace.name for trace in traces]

            if not self.colors:
                self.colors.update(**Plot._plotly_color_map(names))

            names_in_legend = { name:False for name in names }

            for (subgroup, name), trace in zip(self.data.get_unique_idx(), traces):
                trace.marker.color = self.colors.get(name)

                # Plotly workaround:
                # Once a trace name has been added to the legend,
                # we don't want to include future traces in the
                # legend with that name
                if names_in_legend[name]:
                    trace.showlegend = False
                else:
                    names_in_legend[name] = True

                fig.append_trace(trace, row_map[subgroup], COL_NUM)

            return fig

    def render_go(self,
        x: List[Any],
        y: List[float],
        name: str,
        marker: Mapping,
        error_y: Mapping) -> go.Scatter:
        """Render bar graph object

        Parameters
        ----------
        x: List[Any]
            X values (time_point). Values can be linear or
            categorical, e.g., ['P0', 'P1']).
        y: List[float]
            Y values
        name: str
            Name of particular subgroup, e.g., 'Male'
        error_y: Mapping
            Error values for each bar for error bars.

        Returns
        -------
        Plotly Bar Graph Obect
        """
        return go.Scatter(
            x=x if isinstance(x, Iterable) else [x],
            y=y if isinstance(y, Iterable) else [y],
            name=name,
            marker=marker,
            error_y=error_y,
            line=dict(
                shape='spline'
            )
        )

class PlotFactory:
    """Factory class for creating plots."""

    @staticmethod
    def create_plot(type: str, group_by=None, **kwargs) -> Plot:
        """Factory method for creating plots."""
        try:
            data = kwargs['data']
        except KeyError as err:
            raise KeyError(f"{err} argument for plot is missing.")

        summarize = True
        if type == 'violin':
            summarize = False

        if group_by is not None:
            kwargs['data'] = PlotDataFactory.wrap_data(data, group_by, summarize)
        else:
            kwargs['data'] = PlotDataFactory.wrap_data(data, summarize=summarize)

        if type == 'violin':
            return ViolinPlot(**kwargs)
        elif type == 'bar':
            return BarPlot(**kwargs)
        elif type == 'line':
            return LinePlot(**kwargs)
        else:
            raise ValueError(f"'{type}'' is an unsupported plot type.")

    @staticmethod
    def get_config() -> Mapping:
        """Get config for Plotly chart."""
        return dict(
            showLink=False,
            displaylogo=False,
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

@dataclass
class PlotData(ABC):
    data: DataFrame = field(repr=False)
    index: List[str]
    values: DataFrame = field(init=False, repr=False)
    stats: DataFrame = field(init=False, repr=False)
    has_replicates: bool = field(init=False)
    summarize: bool

    def __post_init__(self) -> None:
        if self.data.columns.contains('replicate'):
            self.has_replicates = True
        else:
            self.has_replicates = False

        if self.summarize and self.has_replicates:
            self.values = self._get_statistic(calculate='mean').sort_index()
            self.stats = self._get_statistic(calculate='std').sort_index()
        else:
            self.values = (
                self.data
                    .set_index(self.index)
                    .sort_index()
            )
            self.stats = None


    def _get_statistic(self, calculate: str) -> DataFrame:
        """Performs statistical aggregation on the given expression Pandas DataFrame.

        Args:
            expression:
                Pandas DataFrame of expression values.
            groupby_cols:
                A list of column names to group by and aggregate on.
            calculate:
                String representing the statistical calculation to compute. Must be
                "mean", "std", or "sem".

        Returns:
            Multi-level indexed Pandas DataFrame. Indices are built from grouping
            the top 3 with the most unique values.

        Raises:
            ValueError:
                calcuate is not "mean", "std", or "sem", or the "replicate"
                column is missing from the expression dataframe.
        """
        if calculate.lower() not in ['mean', 'std', 'sem']:
            raise ValueError("calculate must be 'mean', 'std', or 'sem'.")
        if not self.has_replicates:
            raise ValueError("Currently only performing statistics when 'replicate' column is present.")

        df_group = (
            self.data
                .drop('replicate', axis=1)
                .groupby(self.index)
        )

        compute_stat_func = getattr(df_group, calculate.lower())

        # TODO: Aggregate all statistics on data to avoid doing multiple
        # calls to _get_statistic
        # return (
        #     self.data
        #         .drop('replicate', axis=1)
        #         .groupby(self.index)
        #         .agg({
        #             'raw_value': ['mean', 'std', 'sem']
        #         })
        #         .dropna()
        # )
        return compute_stat_func().dropna()

    @abstractmethod
    def get_unique_idx(self):
        """TODO: Documentation"""
        raise NotImplementedError

@dataclass
class DataWithOneIndex(PlotData):
    """Data container for expression data with one index."""
    def get_unique_idx(self) -> List:
        """Get list of unique indices"""
        return (
            self.values
                .index
                .unique()
                .tolist()
        )

@dataclass
class DataWithTwoIndices(PlotData):
    """Data container for expression data with two index."""

    def get_unique_idx(self) -> List:
        """Get list of unique indices"""
        # Example
        # time_point  condition  raw_value
        #         P0  DMSO_GFP+  224.649719
        #             DMSO_GFP-  158.241577
        #         P1  DAPT_GFP+  222.244827
        #             DAPT_GFP-  166.381897
        # __________________________________
        # We drop the last index because it has the most number of values.
        # This way when we are plotting we can properly group the traces.
        # The example above shows a data frame with two indices, time_point
        # and condition. By dropping condition we can then just slice the
        # data frame with df.loc['P0'] and plot the appropriate
        # x=['DMSO_GFP+', 'DMSO_GFP-] and y=[224.649719, 158.241577] with name='P0'
        high_num_vals_level = self.index[-1]
        unique_idx = (
            self.values
                .sort_index(level=high_num_vals_level)
                .index
                .droplevel(high_num_vals_level)
                .unique()
                .tolist()
        )
        return unique_idx

@dataclass
class DataWithThreeIndices(DataWithTwoIndices):
    """Data container for expression data with one index."""
    pass

class PlotDataFactory:
    """Factory class to wrap expression data."""

    @staticmethod
    def wrap_data(data: DataFrame, group_by=None, summarize=True) -> PlotData:
        """Wrap expression data in the appropriate data container
        depending on how many indices it contains."""

        if group_by is None:
            group_by = PlotDataFactory.get_groupby_conditions(data)

        data_args = (data, group_by[:]) # use group_by[:] to shallow copy index

        indices_length = len(group_by)
        if indices_length == 1:
            return DataWithOneIndex(*data_args, summarize)
        elif indices_length == 2:
            return DataWithTwoIndices(*data_args, summarize)
        elif indices_length == 3:
            return DataWithThreeIndices(*data_args, summarize)
        else:
            raise NotImplementedError("Currently only support data with 1, 2, or 3 indices.")

    @staticmethod
    def get_groupby_conditions(data: DataFrame) -> List[str]:
        """Returns the top 3 conditions containing the most unique values.

        Excludes columns: 'raw_value', 'replicate', 'time_unit', 'time_point_order'

        Args:
            expression:
                Pandas DataFrame of expression values.

        Returns:
            A list of conditions (strings) in ascending order with most
            unique values.
        """
        excluded_columns = ['raw_value', 'replicate', 'time_unit', 'time_point_order']
        filter_cols = ~data.columns.isin(excluded_columns)
        conditions = (
            data.loc[:, filter_cols]
                .nunique()
                .sort_values()
                .index
                .tolist()
        )
        top_3_conditions = conditions[-3:]
        return top_3_conditions

def get_available_plot_types(data: DataFrame, index: List) -> List[Mapping]:
    """Get available plot types for expression data based on index.

    Parameters
    ----------
    data
        Expression data frame.

    Returns
    -------
    List of available plot times, where each element is a dict, e.g.,
    [{plot_type:'Bar'}, {plot_type:'Line'}, {plot_type:'Violin'}].
    That way the serialized JSON is more explicit when handling
    plot types on the client.
    """
    available_plots = list()
    if is_violin_available(data):
        available_plots.append(dict(plot_type='violin'))
    if is_line_available(data) and index[-1] == 'time_point':
        available_plots.append(dict(plot_type='line'))
    if is_bar_available(data):
        available_plots.append(dict(plot_type='bar'))

    return available_plots

def is_violin_available(data: DataFrame) -> bool:
    """
    Check if violin plot is an available visualization
    given the expression data frame.

    A violin plot is available if expression data:
        1. As long as there is a categorical column to filter on

    Parameters
    ----------
    data: DataFrame
        Expression data frame

    Returns
    -------
    Return true if expression data follows violin plot requirements,
    false otherwise.
    """
    return False

def is_line_available(data: DataFrame) -> bool:
    """
    Check if line plot is an available visualization
    given the expression data frame.

    A line plot is available if expression data:
        1. Inludes time_point column
        2. TODO: What else?

    Parameters
    ----------
    data: DataFrame
        Expression data frame

    Returns
    -------
    Return true if expression data follows line plot requirements,
    false otherwise.
    """
    return data.columns.contains('time_point')

def is_bar_available(data: DataFrame) -> bool:
    """
    Check if bar plot is an available visualization
    given the expression data frame.

    A bar plot is available if expression data:
        1. TODO

    Parameters
    ----------
    data: DataFrame
        Expression data frame

    Returns
    -------
    Return true if expression data follows bar plot requirements,
    false otherwise.
    """
    return True

def is_svg_available(data: DataFrame) -> bool:
    """TODO"""
    pass
