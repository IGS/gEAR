import traceback
import datashader as ds
import holoviews as hv
import hvplot
import hvplot.pandas  # noqa
import numpy as np
import panel as pn
import param
from common import (
    create_spatial_plot, create_umap_plot, create_violin_plot,
    retrieve_dataframe, retrieve_image_array, normalize_expression_name, has_selection, Settings
)

# CRITICAL: Initialize the Bokeh backend for interactivity
hvplot.extension('bokeh', logo=False) # type: ignore
pn.extension(loading_indicator=True, defer_load=True, nthreads=4)  # type: ignore)

class BaseSpatialViewer(pn.viewable.Viewer):
    """
    Base Viewer component. Handles state and linking.
    """

    settings = Settings()

    use_clusters = param.Boolean(default=False, doc="Whether to show clusters or gene expression in the main plots")

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
        # TODO: bind to a slider, which will again filter the original dataframe for updating the plots.
        if self.settings.min_genes and self.settings.min_genes > 0:
            self.df = self.orig_df[self.orig_df['n_genes_by_counts'] >= self.settings.min_genes]
        else:
            self.df = self.orig_df

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
        #self.bounds_stream_image.param.watch(self._update_bounds_callback, 'bounds')
        #self.bounds_stream_composite.param.watch(self._update_bounds_callback, 'bounds')
        self.bounds_stream_image.add_subscriber(self._update_bounds_callback)
        self.bounds_stream_composite.add_subscriber(self._update_bounds_callback)


        # Add some attributes that will be used in various places
        # This includes precomputing the datashader aggregations since they can be shared across multiple plots

        self.expression_agg = ds.max('raw_value')
        self.expression_cmap = 'YlOrRd'
        self.clusters_agg = ds.count_cat('clusters')
        self.cluster_cmap = dict(zip(self.df['clusters'], self.df['colors']))

        self.bg_image = None
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

            self.bg_image = hv.RGB(self.image_array, bounds=img_bounds).opts(
                        xaxis=None, yaxis=None, frame_width=300, frame_height=200,
                        tools=["box_select"], default_tools=[]
                    )

        self.layout_height = 312  # 360px - tile header height
        if self.settings.display_height and self.settings.display_height > 0: # type: ignore
            self.layout_height: int = self.settings.display_height # type: ignore

        self.layout_width = 1100  # Default width
        if self.settings.display_width and self.settings.display_width > 0: # type: ignore
            self.layout_width: int = self.settings.display_width # type: ignore

        self._init_widgets()

    def _update_bounds_callback(self, event):
        """
        This callback fires automatically when the user draws or clears a box.
        It breaks the tuple into individual params, which location.sync pushes to the URL.
        """
        new_bounds = event.new

        if new_bounds is None:
            # User clicked off/cleared the box
            self.settings.selection_x1 = None
            self.settings.selection_y1 = None
            self.settings.selection_x2 = None
            self.settings.selection_y2 = None
        else:
            # User drew a box
            self.settings.selection_x1, self.settings.selection_y1, self.settings.selection_x2, self.settings.selection_y2 = new_bounds

        # Instantly updates zoom panel with new bounds
        self._update_zoom_panel(new_bounds)

    def _create_ghost_legend(self):
            """Creates a fake, invisible plot just to force Bokeh to draw a legend."""
            ghost_points = []

            # self.cluster_cmap should be a dict like {'Cluster 1': '#FF0000', ...}
            for cluster_name, hex_color in self.cluster_cmap.items():
                # Create a single point at (NaN, NaN) so it doesn't draw on the screen
                # But give it a label and a color so the legend picks it up
                pt = hv.Points(
                    [(np.nan, np.nan)],
                    label=str(cluster_name)
                ).opts(color=hex_color, size=10, tools=[], default_tools=[])

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
            self.zoom_pane = pn.pane.HoloViews(None)

            main_row = pn.Row(image_panel, linked_composite, self.zoom_pane)

            # Lay out the non-zoom panels side-by-side using HoloView

            markdown_width = 400    # Manually measured.
            markdown_padding = 20  # Total left/right padding for the markdown pane
            self.intro_markdown = pn.pane.Markdown(
                "### Click the Expand icon in the top right corner to see all plots", width=markdown_width
            )

            # When width is too short, things break down.
            if self.layout_width < 1100:
                self.layout_width = 1100
                self.layout_height = 312

            spacer_width = self.layout_width - markdown_width - markdown_padding - self.switch_content_width
            if spacer_width < 0:
                spacer_width = 100

            self.pre_layout = pn.Row(
                self.intro_markdown,
                pn.Spacer(width=spacer_width),
                self.switch_layout,
                height=30
            )

            # Return final Panel layout
            return pn.Column(self.pre_layout, main_row)
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

        # Apply the background image
        image_panel = None
        if hasattr(self, 'bg_image') and self.bg_image is not None:
            image_panel = self.bg_image
            main_base = image_panel * plot # type: ignore
            zoom_base = image_panel * plot # type: ignore

        # Add a ghost legend overlay to properly display a cluster legend
        if self.use_clusters:
            ghost_legend = self._create_ghost_legend()
            # Overlay (multiply) the ghost legend on top, and explicitly tell it to show the legend
            main_base = (main_base * ghost_legend).opts(show_legend=True, legend_position='right')
            zoom_base = (zoom_base * ghost_legend).opts(show_legend=True, legend_position='right')

        # Overlay the background image over the other plots
        if hasattr(self, 'bg_image') and self.bg_image is not None:
            # The image slides safely underneath the interactive linked points
            main_composite = self.bg_image * main_base
            zoom_composite = self.bg_image * zoom_base
        else:
            main_composite = main_base
            zoom_composite = zoom_base

        # Attach the stream to capture drawn boxes
        self.bounds_stream_composite.source = main_composite

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
        self.clusters_switch = pn.widgets.Switch.from_param(self.param.use_clusters, name="Show Clusters", margin=(10, 10))

        # Using HTML to center the labels (https://github.com/holoviz/panel/issues/1313#issuecomment-1582731241)
        self.switch_content_width = 250
        self.switch_layout = pn.Row(
            self.clusters_switch,
            width=self.switch_content_width
        )

    def _update_zoom_panel(self, bounds):
        """
        Updates the zoom panel based on the provided bounds.

        This method is called when the user draws a box on one of the plots to zoom in on that area.
        It updates the `zoom_pane` dynamic map with the new bounds, which triggers a re-render of the zoomed-in view.

        Parameters:
        - bounds: A tuple containing the new bounds in the format (left, bottom, right, top).
                  If bounds is None, it indicates that the user has cleared the selection.
        """
        if not hasattr(self, 'zoom_pane') or self.zoom_pane.object is None:
            return

        if bounds is not None:
            l, b, r, t = bounds

            # 1. Protect against inverted Y-axes ignoring the zoom
            xlim = (min(l, r), max(l, r))
            ylim = (min(b, t), max(b, t))

            # 2. clone=True forces Panel to realize this is a new object to push to the frontend
            new_zoomed_obj = self.zoom_pane.object.opts(xlim=xlim, ylim=ylim, clone=True)
            self.zoom_pane.object = new_zoomed_obj
        else:
            # Clear the limits
            new_zoomed_obj = self.zoom_pane.object.opts(xlim=(None, None), ylim=(None, None), clone=True)
            self.zoom_pane.object = new_zoomed_obj

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