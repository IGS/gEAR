import json
import os, sys
import plotly.offline as py
import plotly.graph_objs as go

import pandas as pd
import scanpy as sc
sc.settings.verbosity = 0

import geardb

def get_groupby_conditions(expression=None):
    """
    Given an expression dataframe, returns the top 3 conditions containing the most
    unique values.

    Excludes columns: 'raw_value', 'replicate', 'time_unit'
    """
    from operator import itemgetter

    excluded_columns = ['raw_value', 'replicate', 'time_unit']

    #Get conditions and count unique values
    col_counts = []
    for col in expression.columns.values:
        if col in excluded_columns:
            continue

        count = expression[col].nunique()
        col_counts.append({'column': col, 'count': count})

    # Sort conditions from largest to smallest
    sorted_col_counts = sorted(col_counts, key=itemgetter('count'), reverse=True)
    sorted_cols = [col['column'] for col in sorted_col_counts]

    #keep top 3 groups
    top_groups = sorted_cols[0:3]

    #reverse groups for pandas groupby
    return list(reversed(top_groups))


def get_expression(dataset_id=None, gene_id=None):
    """
    With a dataset_uid and gEAR gene_id, returns a pandas dataframe containing the
    expression values and observation characteristics for that gene.

    Expression values are in column 'raw_value'

    Makes 2 attempts at getting the expression values:
      1. Gene symbol, but if that fails,
      2. Ensembl ID
    """
    if dataset_id is None:
        raise Exception("Unable to get expression. No 'dataset_id' given.")
    if gene_id is None:
        raise Exception("Unable to get expression. No 'gene_id' given.")
    gene = geardb.get_gene_by_id(gene_id)
    ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
    h5_path = ds.get_file_path()

    # Read file
    adata = sc.read_h5ad(h5_path)

    #Find the gene's index position (integer)
    var_index = pd.Index(adata.var.index)

    # First attempt with gene_symbol, if that fails use ensembl_id
    try:
        gene_index_label = adata.var.index[var['gene_symbol'].lower() == gene.gene_symbol.lower()].tolist()[0] #Get 1st found
        gene_index = var_index.get_loc(gene_index_label) #Get the index label's position in the index
    except:
        gene_index = var_index.get_loc(gene.ensembl_id) #Get the index label's position in the index

    expression = pd.DataFrame(adata.X[:,gene_index], columns=['raw_value'])

    # Add observation columns to expression
    for col in adata.obs.columns:
        expression[col] = adata.obs[col].values

    return expression

def get_plot_types(expression=None):
    """
    Based on expression data, determine possible plot types
    """
    if expression is None:
        raise Exception("Unable to get plot type. Argument 'expression' is None. Provide one to continue")

    plot_types = ['bar', 'line', 'violin']
    exclude_cols = ['raw_value']
    columns = expression.columns.difference(exclude_cols).tolist()

    # TODO: Try to remove each plot_type from the plot_types list
    # Violin
    if 'cell_type' in columns:
        # Count number of celltypes
        celltype_count = expression['cell_type'].nunique()

        # Only 1 violin - not best representation for this dataset
        if celltype_count == 1:
            plot_types.remove('violin')

    # Assumes Violins require 'cell_type' column
    if 'cell_type' not in columns:
        plot_types.remove('violin')

    if 'replicate' in columns and 'violin' in plot_types:
        plot_types.remove('violin')

    # Line
    if 'time_point' not in columns:
        plot_types.remove('line')

    # Bar - TODO

    # SVG - TODO Later/Last

    return plot_types

def get_statistic(expression=None, groupby_cols=None, calculate=None):
    """
    Return the stat calculation on the given expression pandas dataframe

    Note: returned dataframe is multi-level indexed dataframe. indices are build
        from groupby() on the top 3 conditions
        stat.loc[condition3, condition2, condition1]['raw_value']
    """
    if 'replicate' not in expression.columns:
        raise Exception("Error: Currently only performing statistics when 'replicate' column is present")
    if groupby_cols is None:
        raise Exception("Error: Unable to get statistic. Please provide 'groupby_cols' list")
    if calculate.lower() not in ['mean', 'std', 'sem']:
        raise Exception('Unable to get statistic. calculation given has not "mean", "std", or "sem"')

    #drop replicate column
    df = expression.drop('replicate', 1)
    df_grouped = df.groupby(groupby_cols, sort=False) #not sorting improves speed

    #calculate stat
    stat = None
    if calculate.lower() == 'mean':
        stat = df_grouped.mean()
    if calculate.lower() == 'std':
        stat = df_grouped.std()
    if calculate.lower() == 'sem':
        stat = df_grouped.sem()

    # returns pandas dataframe
    return stat

def get_plot(plot_type=None, expression=None):
    '''
    Returns a plotly graph in HTML.

    plot_type = 'violin', 'bar', or 'line'
    expression = pandas dataframe of a gene's expression raw_values

    NOTE:
        - If 'replicate' column is present, the expression datframe is grouped
        using the top 3 observation characteristics with the most unique values.
        - Presently calculates 'mean' and 'std' if 'replicate' is present.

    violin
    ------
        - NOTE: Assumes column 'cell_type' is present
        - TODO: Are there other columns that can produce violins?
    bar
    ------
        - Creates subplots if dataframe has 3 groups. The N of smallest
            group is used to create N subplots
    line
    ------
        - NOTE: Assumes column 'time_point' is present
        - Creates subplots if dataframe has 3 groups. The N of smallest
            group is used to create N subplots
    '''
    groupby_cols = get_groupby_conditions(expression=expression)

    # Get statistics if replicates are included
    if 'replicate' in expression.columns:
        values = get_statistic(expression=expression, groupby_cols=groupby_cols, calculate='mean')
        stats = get_statistic(expression=expression, groupby_cols=groupby_cols, calculate='std')
    else:
        values = expression
        stats = None


    if plot_type == 'violin':
        #Example (Multiple Traces): https://plot.ly/python/violin/
        celltypes = expression['cell_type'].unique()
        celltype_count = len(celltypes)
        data = []
        for i in range(0, celltype_count):
            trace = {
                "type": "violin",
                "x": expression["cell_type"][expression["cell_type"] == celltypes[i]].tolist(),
                "y": expression["raw_value"][expression["cell_type"] == celltypes[i]].tolist(),
                "name": celltypes[i],
                "box": {
                    "visible": True
                },
                "meanline": {
                    "visible": True
                },
                "points": "all",
                "pointpos": expression["raw_value"][expression["cell_type"] == celltypes[i]].tolist(),
                "jitter": 0,
                "scalemode": "count",
                "meanline": {
                    "visible": True
                },
                "span": [
                    0
                ]
            }
            data.append(trace)

        fig = {
            "data": data,
            "layout" : {
                "title": "",
                "yaxis": {
                    "zeroline": False,
                },
                # "autosize": True,
                # "height": 350,
                "legend": {
                    "orientation": "h"
                },
                "margin": {
                    "l": 50,
                    "r": 50,
                    "b": 75,
                    "t": 25,
                    "pad": 4
                }
            }
        }

    if plot_type == 'bar':
        #Prepare data traces
        if len(groupby_cols) == 1:
            try:
                cond1_vals = tuple(values.index.get_level_values(groupby_cols[0]).unique())
            except:
                #'replicate' column therefore no multi-level index
                cond1_vals = tuple(values[groupby_cols[0]].unique())

            data = []
            for i in range(0, len(cond1_vals)):
                if 'cell_type' in values.columns:
                    x = cond1_vals[i]
                    y = values['raw_value'][values['cell_type'] == cond1_vals[i]].tolist()
                else:
                    x = [ cond1_vals[i] ]
                    y = [ values.loc[cond1_vals[i]]['raw_value'] ]

                if stats is not None:
                    error_y = {
                        "type": "data",
                        "array": [ stats.loc[cond1_vals[i]]['raw_value'] ]
                    }
                else:
                    error_y = {}

                trace = {
                    "type": "bar",
                    "x": x,
                    "y": y,
                    "error_y": error_y,
                    "name": cond1_vals[i]
                    }
                data.append(trace)

            fig = {
                "data": data,
                "layout": {
                    "xaxis": {
                        "type": "category"
                    },
                    # "autosize": True,
                    # "height": 350,
                    "legend": {
                        "orientation": "h"
                    },
                    "margin": {
                        "l": 50,
                        "r": 50,
                        "b": 125,
                        "t": 25,
                        "pad": 4
                    }
                }
            }

        if len(groupby_cols) == 2:
            cond1_vals = tuple(values.index.get_level_values(groupby_cols[0]).unique())
            cond2_vals = tuple(values.index.get_level_values(groupby_cols[1]).unique())

            data = []
            for cond2_val in cond2_vals:
                trace = {
                    "type": "bar",
                    "x": cond1_vals,
                    "y": list(values.loc[cond1_val, cond2_val]['raw_value'] for cond1_val in cond1_vals),
                    "error_y": {
                        "type": "data",
                        "array": list(stats.loc[cond1_val, cond2_val]['raw_value'] for cond1_val in cond1_vals)
                    },
                    "name": cond2_val
                }
                data.append(trace)

            fig = {
                "data": data,
                "layout": {
                    "xaxis": {
                        "type": "category"
                    },
                    # "autosize": True,
                    # "height": 350,
                    "legend": {
                        "orientation": "h"
                    },
                    "margin": {
                        "l": 50,
                        "r": 50,
                        "b": 125,
                        "t": 25,
                        "pad": 4
                    }
                }
            }

        if len(groupby_cols) == 3:
            cond1_vals = tuple(values.index.get_level_values(groupby_cols[0]).unique())
            cond2_vals = tuple(values.index.get_level_values(groupby_cols[1]).unique())
            cond3_vals = tuple(values.index.get_level_values(groupby_cols[2]).unique())

            #Redirect stdout messages. tools.make_subplots() prints causing javascript-cgi to error
            sys.stdout = open(os.devnull, 'w')

            from plotly import tools

            # Get number of subplots
            fig_rows = len(cond1_vals)
            fig = tools.make_subplots(rows=fig_rows, cols=1,
                subplot_titles=cond1_vals)

            for cond1_val in cond1_vals:
                fig_row = cond1_vals.index(cond1_val) + 1
                for cond2_val in cond2_vals:
                    trace = {
                        "type": "bar",
                        "x": cond3_vals,
                        "y": values.loc[cond1_val, cond2_val, cond3_vals]['raw_value'],
                        "error_y": {
                            "type": "data",
                            "array": stats.loc[cond1_val, cond2_val, cond3_vals]['raw_value']
                        },
                        "name": cond2_val
                    }
                    fig.append_trace(trace, fig_row, 1)

            #Reset stdout
            sys.stdout = sys.__stdout__

            # Layout and Figure options
            # fig['layout']['autosize'] = False
            # fig['layout']['height'] = 350
            # fig['layout']['width'] = 450
            fig['layout']['legend'] = {'orientation': 'h'}
            fig['layout']['margin'] = {'l':50, 'r':50, 'b':75, 't':25, 'pad':4}

    if plot_type == 'line':
        #NOTE: Assumes plot is time course uses 'time_point' column

        #Get unique values for each grouped condition
        time_points = tuple(values.index.get_level_values('time_point').unique())
        groupby_cols.remove('time_point')

        if len(groupby_cols) == 1:
            cond1_vals = tuple(values.index.get_level_values(groupby_cols[0]).unique())
            data = []
            for cond1_val in cond1_vals:
                trace = {
                    "type": "scatter",
                    "x": time_points,
                    "y": [values.loc[cond1_val, time_point]['raw_value'] for time_point in time_points],
                    "error_y": {
                        "type": "data",
                        "array": [stats.loc[cond1_val, time_point]['raw_value'] for time_point in time_points]
                    },
                    "name": cond1_val
                }
                data.append(trace)

            fig = {
                "data": data,
                "layout": {
                    # "autosize": True,
                    # "height": 350,
                    "legend": {
                        "orientation": "h"
                    },
                    "margin": {
                        "l": 50,
                        "r": 50,
                        "b": 75,
                        "t": 25,
                        "pad": 4
                    }
                }
            }

        if len(groupby_cols) == 2:
            cond1_vals = tuple(values.index.get_level_values(groupby_cols[0]).unique())
            cond2_vals = tuple(values.index.get_level_values(groupby_cols[1]).unique())

            #Redirect stdout messages. tools.make_subplots() prints causing javascript-cgi to error
            sys.stdout = open(os.devnull, 'w')

            from plotly import tools

            # Get number of subplots
            fig_rows = len(cond1_vals)
            fig = tools.make_subplots(rows=fig_rows, cols=1,
                subplot_titles=cond1_vals)

            for cond1_val in cond1_vals:
                fig_row = cond1_vals.index(cond1_val) + 1
                for cond2_val in cond2_vals:
                    trace = {
                        "type": "scatter",
                        "x": time_points,
                        "y": values.loc[cond1_val, cond2_val, time_points]['raw_value'],
                        "error_y": {
                            "type": "data",
                            "array": stats.loc[cond1_val, cond2_val, time_points]['raw_value']
                        },
                        "name": cond2_val
                    }

                    fig.append_trace(trace, fig_row, 1)

            sys.stdout = sys.__stdout__

        # Layout and Figure options
        # fig['layout']['autosize'] = True
        # fig['layout']['height'] = 350
        fig['layout']['legend'] = {'orientation': 'h'}
        fig['layout']['margin'] = {'l':50, 'r':50, 'b':75, 't':25, 'pad':4}

    fig['layout']['height'] = 350
    fig['layout']['width'] = 400

    #Nevermind this. the plot div ends up being 700px.
    # https://github.com/plotly/plotly.py/blob/master/plotly/offline/offline.py#L472
    # fig['layout']['width'] = '100%'

    html = py.plot(fig, output_type='div', validate=False, include_plotlyjs=False, show_link=False)
    return html


class Plot:
    def __init__(self, dataset_id=None, gene_id=None, gene_symbol=None,
        gene_ensembl_id=None, expression=None, types=None, preferred_type=None):

        self.dataset_id = dataset_id
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol
        self.gene_ensembl_id = gene_ensembl_id
        self.expression = expression
        self.types = types
        self.preferred_type = preferred_type #Should be user_preference, else plot_default

        # Get gene_symbol and ensembl_id (if gene_id is given)
        if self.gene_id is not None:
            gene = geardb.get_gene_by_id(self.gene_id)
            if self.gene_symbol is None:
                self.gene_symbol = gene.gene_symbol
            if self.gene_ensembl_id is None:
                self.gene_ensembl_id = gene.ensembl_id

        # Get expression data (if dataset_id and gene_id are given)
        if self.expression is None:
            if self.dataset_id is not None and self.gene_id is not None:
                self.expression = get_expression(dataset_id=self.dataset_id, gene_id=self.gene_id)

        # Get list of plot types the given expression data can generate
        if self.types is None:
            if self.expression is not None:
                self.types = get_plot_types(self.expression)

    def get_preferred_plot(self, preferred_type=None):
        """
        This is just like get_plots() method below, but only generates the preferred
        plot type for a dataset and returns only the plotly generated HTML.
        """

        if preferred_type is None:
            if self.preferred_type is None:
                raise Exception("Unable to get plot. No preferred_type found.")
            else:
                preferred_type = self.preferred_type

        #Get plot HTML
        html = get_plot(plot_type=preferred_type, expression=self.expression)

        # Return plot HTML
        return html


    def get_plots(self, types=None):
        """
        Given expression dataframe and plot types, get a plot for each plot type.
        The idea is to supply each plot type so the user can later toggle between graphs.

        If no types are passed, the plot object's own supported types are checked instead.

        Returns a list of dicts where {plot_type: 'string of html of plot'}
            [{'bar': '<html>...</html>'}, {'line': '<html></html>'}]
        """
        if types is None:
            plot_types = self.types
        else:
            plot_types = types

        if len(plot_types) == 0:
            raise Exception("Unable to get plot. No plot types found.")

        plots = []

        # Get all possible plots for dataset
        for plot_type in plot_types:
            #Get plot HTML
            html = get_plot(plot_type=plot_type, expression=self.expression)

            # Add plot to
            plots.append({'dataset_id': self.dataset_id, 'plot_type': plot_type, 'plot_html': html})

        # Return a list of plots for the dataset
        return plots
