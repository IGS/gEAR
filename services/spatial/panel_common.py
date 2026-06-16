import traceback

import datashader as ds
import holoviews as hv
import hvplot
import hvplot.pandas  # noqa
import numpy as np
import panel as pn
import param
from common import (
    Settings,
    create_spatial_plot,
    create_umap_plot,
    create_violin_plot,
    has_selection,
    normalize_expression_name,
    retrieve_dataframe,
    retrieve_image_array,
)

# CRITICAL: Initialize the Bokeh backend for interactivity
hvplot.extension('bokeh', logo=False) # type: ignore
pn.extension(loading_indicator=True, defer_load=True, nthreads=4)  # type: ignore)

class BaseSpatialViewer(pn.viewable.Viewer):
    """
    Base Viewer component. Handles state and linking.
    """

    settings = Settings()

    # Any potential bound widgets should be declared here so they are accessible in the base class and subclasses
    use_clusters = param.Boolean(default=False, doc="Whether to show clusters or gene expression in the main plots")
    min_genes = param.Integer(default=0, doc="Minimum number of genes expressed to include a cell observation", bounds=(0, 500))

    def __init__(self, **params):
        super().__init__(**params)
        """
        DataFrame columns
        raw_value,spatial1,spatial2,n_genes_by_counts,UMAP1,UMAP2,clusters,clusters_cat_codes,colors
        """

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
                    "save": "save",
                    "display_name": "display_name",
                    "make_default": "make_default",
                },
            )

        self.orig_df = retrieve_dataframe(self.settings.dataset_id, self.settings.filename)
        self.orig_df['clusters'] = self.orig_df['clusters'].astype('category')

        # If min_genes is set, filter the dataframe to only include observations with at least that many genes
        self.min_genes = self.settings.min_genes
        self._filter_df()

        self.image_array = retrieve_image_array(self.settings.dataset_id)

        self.current_gene = normalize_expression_name(self.settings.filename)

        # If selection_x1/x2/y1/y2 are present save as a tuple in the form of (left, right, bottom, top)
        saved_bounds = None
        if has_selection(self.settings):
            saved_bounds = (
                self.settings.selection_x1,
                self.settings.selection_x2,
                self.settings.selection_y1,
                self.settings.selection_y2,
            )


        self.saved_bounds = saved_bounds

        # Initialize linking and streams
        self.linker = hv.link_selections.instance(unselected_alpha=0.5)
        self.bounds_stream_image = hv.streams.BoundsXY(bounds=self.saved_bounds)  # type: ignore
        self.bounds_stream_composite = hv.streams.BoundsXY(bounds=self.saved_bounds)  # type: ignore

        # Set up a callback to update the URL params whenever the user draws or clears a box
        self.bounds_stream_image.add_subscriber(self._update_bounds_callback)
        self.bounds_stream_composite.add_subscriber(self._update_bounds_callback)


        # Add some attributes that will be used in various places
        # This includes precomputing the datashader aggregations since they can be shared across multiple plots

        self.expression_agg = ds.max('raw_value')
        self.expression_cmap = 'YlOrRd'

        self.clusters_agg = ds.count_cat('clusters')
        self.cluster_cmap = dict(zip(self.df['clusters'], self.df['colors']))

        self.bg_image = None
        self.bg_image_dimmed = None
        self.img_height = None
        self.img_width = None
        if self.image_array is not None:
            # Ensure array contents are UInt8 (0-255) for proper display. If not, normalize to that range and convert.
            if self.image_array.dtype != np.uint8:
                self.image_array = (255 * (self.image_array - np.min(self.image_array)) / (np.ptp(self.image_array) + 1e-8)).astype(np.uint8)

            # use image_array shape to build image bounds. Image is 3-dimensional where shape is (y, x, c)
            img_bounds = (0, 0, self.image_array.shape[1], self.image_array.shape[0])
            self.img_height = self.image_array.shape[0]
            self.img_width = self.image_array.shape[1]

            # This is a fully opaque version
            self.bg_image = hv.RGB(self.image_array, bounds=img_bounds).opts(
                        xaxis=None, yaxis=None,
                        tools=["box_select"], default_tools=[]
                    )

            # Create a more opaque version of the image.
            alpha_channel = np.full((self.img_height, self.img_width, 1), int(255 * 0.5), dtype=np.uint8)

            # Check if image is RGB (3 channels). If so, append alpha. If already RGBA, just overwrite alpha.
            if self.image_array.shape[-1] == 3:
                dimmed_array = np.concatenate([self.image_array, alpha_channel], axis=-1)
            else:
                dimmed_array = self.image_array.copy()

                # SAdkins note - '...' means all rows in this case.
                dimmed_array[..., 3] = int(255 * 0.5)

            self.bg_image_dimmed = hv.RGB(dimmed_array, bounds=img_bounds).opts(
                        xaxis=None, yaxis=None, responsive=True,
                        tools=["box_select"], default_tools=[]
                    )

        self._init_widgets()

    def _filter_df(self):
        """Applies any necessary filtering to the dataframe based on the current settings."""
        df = self.orig_df
        if self.min_genes and self.min_genes > 0:
            self.df = df[df['n_genes_by_counts'] >= self.min_genes]
        else:
            self.df = df

    def _update_bounds_callback(self, bounds=None):
        """
        This callback fires automatically when the user draws or clears a box.
        It breaks the tuple into individual params, which location.sync pushes to the URL.
        """

        if bounds is None:
            # User clicked off/cleared the box
            self.settings.selection_x1 = None
            self.settings.selection_y1 = None
            self.settings.selection_x2 = None
            self.settings.selection_y2 = None
        else:
            # User drew a box
            (self.settings.selection_x1, self.settings.selection_y1, self.settings.selection_x2, self.settings.selection_y2) = bounds

        # Trigger the zoom update!
        self._update_zoom_panel(bounds)

    def _create_ghost_legend(self):
            """Creates a fake, invisible plot just to force Bokeh to draw a legend."""
            ghost_points = []

            # TODO: add some "click" interactivity to hide plot elements based on if legend item is enabled or disabled.

            # Hook to strip the invisible padding Bokeh reserves for axes
            def remove_bokeh_borders(plot, element):
                plot.state.min_border = 0
                plot.state.min_border_left = 0
                plot.state.min_border_right = 0
                plot.state.min_border_top = 0
                plot.state.min_border_bottom = 0

            # self.cluster_cmap should be a dict like {'Cluster 1': '#FF0000', ...}
            for cluster_name, hex_color in self.cluster_cmap.items():
                # Create a single point at (NaN, NaN) so it doesn't draw on the screen
                # But give it a label and a color so the legend picks it up
                pt = hv.Points(
                    [(np.nan, np.nan)],
                    label=str(cluster_name)
                ).opts(color=hex_color, size=10, tools=[], default_tools=[], hooks=[remove_bokeh_borders])

                ghost_points.append(pt)

            # Combine all the ghost points into a single overlay
            return hv.Overlay(ghost_points)

    def _init_widgets(self):
        """
        This is where you would initialize any Panel widgets (sliders, dropdowns, etc.) that you want to use in your app.
        You can then reference these widgets in your _build_layout method to include them in the layout and set up callbacks.
        """
        raise NotImplementedError("Subclasses must implement _init_widgets")

    def _build_layout(self):
        """
        This is where the main Panel layout is built.
        We leave it blank in the base class since the Condensed and Expanded viewers will have different layouts.
        """
        raise NotImplementedError("Subclasses must implement _build_layout")

    def _update_zoom_panel(self, bounds):
        """
        This is where you would implement the logic to update the zoomed-in plot based on the provided bounds.
        The bounds parameter will be a tuple in the form of (left, right, bottom, top) representing the coordinates of the box drawn by the user.
        You can use these bounds to set the xlim and ylim of the zoomed-in plot accordingly.
        """
        raise NotImplementedError("Subclasses must implement _update_zoom_panel")

    def __panel__(self):
        """
        Panel automatically looks for this method.
        It MUST return a Panel viewable object (Row, Column, Pane, etc.)
        """
        return self._build_layout()

class CondensedSpatialViewer(BaseSpatialViewer):
    """
    The specific component for your panel_app.py
    """

    def _build_layout(self):
        """Builds the 3-panel condensed row."""
        try:

            # One unfortunately annoyance is that datashader's default behavior is to flip the y-axis,
            # which is not what we want for spatial data. To fix this,
            # we can reverse the y-axis limits by setting ylim to (max, min) instead of (min, max).
            self.df["y_plot"] = self.df["spatial2"]
            if self.img_height is not None:
                self.df["y_plot"] = self.img_height - self.df["spatial2"]

            # Generate base plots
            image_panel = None
            if hasattr(self, 'bg_image'):
                image_panel = self.bg_image
                # Attach the stream to capture drawn boxes
                self.bounds_stream_image.source = image_panel

            linked_composite = self._add_center_plot
            self.zoom_pane = pn.pane.HoloViews(None, sizing_mode="stretch_both", linked_axes=False)

            # dedicated pane for the legend so plots can be resized without the legend affecting them
            self.legend_pane = pn.pane.HoloViews(None)

            self.legend_container = pn.Column(
                        self.legend_pane,
                        width=200,
                        height=275,
                        scroll=True,     # Activates the native CSS scrollbar
                        visible=False,
                        margin=(0,0,0,0)
                    )


            self.main_row = pn.Row(image_panel,
                            linked_composite,
                            self.zoom_pane,
                            self.legend_container,
                            sizing_mode='stretch_width'
                            )

            # Lay out the non-zoom panels side-by-side using HoloViews
            self.intro_markdown = pn.pane.Markdown(
                "### Click the Expand icon in the top right corner to see all plots", #width=markdown_width
            )

            self.pre_layout = pn.Row(
                self.intro_markdown,
                pn.Spacer(),
                self.switch_layout,
            )

            # Return final Panel layout
            return pn.Column(self.pre_layout,
                            self.main_row,
                            sizing_mode='stretch_both' # Fills the 100%x100% iframe
                            )
        except Exception as e:
            traceback.format_exc()
            return pn.pane.Alert(f"Error: {e}", alert_type="danger")

    @param.depends("use_clusters")
    def _add_center_plot(self):
        if self.use_clusters:
            plot =  create_spatial_plot(self.df, self.clusters_agg, y_col='y_plot', color_col='clusters', cmap=self.cluster_cmap, is_categorical=True) # type: ignore
        else:
            plot =  create_spatial_plot(self.df, self.expression_agg, y_col='y_plot', color_col='raw_value', cmap=self.expression_cmap) # type: ignore

        main_base = plot
        zoom_base = plot

        # Attach the stream to capture drawn boxes
        self.bounds_stream_composite.source = main_base

        # Apply the background image
        if hasattr(self, 'bg_image') and self.bg_image_dimmed is not None:
            main_base = self.bg_image_dimmed * plot # type: ignore
            zoom_base = self.bg_image_dimmed * plot # type: ignore

        if hasattr(self, 'legend_container'):
            self.legend_container.visible = self.use_clusters

        # Push the legend to the 4th pane instead of the center plot
        if hasattr(self, 'legend_pane'):
            if self.use_clusters:
                # Calculate the required width based on the longest cluster name
                needed_width = max(180, max(len(str(name)) for name in self.cluster_cmap) * 7)  # Adjust the multiplier as needed

                # Calculate the required canvas height: ~22px per cluster item + 50px for margins
                needed_height = max(300, len(self.cluster_cmap) * 22 + 50)

                ghost_legend = self._create_ghost_legend().opts(
                    show_legend=True,
                    legend_position="bottom_left",
                    xaxis=None, yaxis=None,
                    show_frame=False, toolbar=None, # Hide the empty plot canvas
                    margin=(0,0,0,0),
                    width=needed_width,
                    height=needed_height   # defined so things do not overlap elsewhere.
                )
                self.legend_pane.object = ghost_legend
            else:
                self.legend_pane.object = None # Clear it for Expression view
        # Overlay the background image over the other plots
        if hasattr(self, 'bg_image') and self.bg_image_dimmed is not None:
            # The image slides safely underneath the interactive linked points
            main_composite = self.bg_image_dimmed * main_base
            zoom_composite = self.bg_image_dimmed * zoom_base
        else:
            main_composite = main_base
            zoom_composite = zoom_base



        # Store the master composite objects so we can apply the zoom limits in the callback without having to rebuild the entire plot from scratch
        self.zoom_composite = zoom_composite

        # Update the zoom pane with the new composite plot (with or without background)
        if hasattr(self, 'zoom_pane'):
            self.zoom_pane.object = zoom_composite

        return main_composite

    def _init_widgets(self):
        """
        Initializes the widgets for the panel layout.

        This method sets up the following components:
        - A switch widget (`clusters_switch`) to toggle the visibility of clusters.
          It is created using the `pn.widgets.Switch.from_param` method and is styled
          with a margin.
        - A layout (`switch_layout`) that organizes the switch widget in a row.
          The layout has a fixed width (`switch_content_width`) and is intended to
          provide a structured arrangement for the widget.

        Note:
        - The layout includes commented-out HTML labels for potential future use
          to center labels using HTML, as referenced in a GitHub issue discussion.
        """
        self.clusters_switch = pn.widgets.Switch.from_param(self.param.use_clusters, name="Show Clusters", margin=(12, 0, 0, 10))

        self.switch_layout = pn.Row(
            self.clusters_switch,
        )

    def _update_zoom_panel(self, bounds):
        """
        This method updates the zoomed-in plot based on the provided bounds.
        The method applies the new limits to the master composite plot and updates the zoom pane accordingly.

        Parameters:
        - bounds: A tuple containing the new bounds in the format (left, bottom, right, top).
                  If bounds is None, it indicates that the user has cleared the selection.
        """

        # Safety check to ensure the master canvas exists
        if not hasattr(self, 'zoom_composite'):
            return

        # Apply the new limits to the master canvas
        if bounds is None:
            new_zoomed_obj = self.zoom_composite.opts(xlim=(None, None), ylim=(None, None), clone=True)
        else:
            x1, y1, x2, y2 = bounds
            new_zoomed_obj = self.zoom_composite.opts(
                xlim=(min(x1, x2), max(x1, x2)),
                ylim=(min(y1, y2), max(y1, y2)),
                clone=True
            )

        # Wrap it in a brand-new pane
        new_pane = pn.pane.HoloViews(new_zoomed_obj, sizing_mode='stretch_both', linked_axes=False)

        # Hot-swap it into the browser DOM (Index 2 is the 3rd panel in the row)
        self.zoom_pane = new_pane
        if hasattr(self, 'main_row'):
            self.main_row[2] = self.zoom_pane


class ExpandedSpatialViewer(BaseSpatialViewer):
    """
    The specific component for your panel_app_expanded.py
    """

    def _build_layout(self):
        # Build your spatial rows, UMAPs, and Violins here...

        try:

            main_plots = None
            zoom_plots = None
            umap_row = None
            violin_row = None

            return pn.Column(main_plots, zoom_plots, umap_row, violin_row)
        except Exception as e:
            traceback.format_exc()
            return pn.pane.Alert(f"Error: {e}", alert_type="danger")

    def _init_widgets(self):
        # ? This can be useful for filtering datasets even for projections, but how best to word it?
        min_slider_width = 300
        self.min_genes_slider = pn.widgets.IntSlider(
            name="Filter - Mininum genes per observation",
            start=0,
            end=500,
            step=25,
            width=min_slider_width,
            value=self.min_genes,
        )

        self.display_name = pn.widgets.TextInput(
            name="Display name",
            placeholder="Name this display to save...",
            width=250,
            visible=not self.nosave,
        )
        self.save_button = pn.widgets.Button(
            name="Save settings", button_type="primary", width=100, align="end"
            , visible=not self.nosave
        )
        self.make_default = pn.widgets.Checkbox(
            name="Make this the default display", value=False
            , visible=not self.nosave
        )

        markdown_width = 675
        spacer_width = markdown_width - min_slider_width    # Make default button should left-align with the above text input

        self.pre_layout = pn.Column(
            pn.Row(
                pn.pane.Markdown(
                    "## Select a region to modify zoomed in view in the bottom panel",
                    height=30,
                    width=markdown_width,
                ),
                self.display_name,
                self.save_button,
            ),
            pn.Row(self.min_genes_slider, pn.Spacer(width=spacer_width), self.make_default),
        )
