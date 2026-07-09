import asyncio
import gc
import json
import sys
import traceback

import datashader as ds
import holoviews as hv
import hvplot
import hvplot.pandas  # noqa
import numpy as np
import pandas as pd
import panel as pn
import param
from common import (
    Settings,
    autohide_toolbar,
    clip_expression_values,
    create_spatial_plot,
    create_umap_plot,
    create_violin_plot,
    has_selection,
    normalize_expression_name,
    retrieve_dataframe,
    retrieve_image_array,
    sort_clusters,
)

# CRITICAL: Initialize the Bokeh backend for interactivity
hvplot.extension('bokeh', logo=False) # type: ignore
pn.extension(loading_indicator=True, defer_load=True, nthreads=4)  # type: ignore)

class BaseSpatialViewer(pn.viewable.Viewer):
    """
    Base Viewer component. Handles state and linking.
    """

    def __init__(self, **params):
        # Override session_args if provided, otherwise use the default.
        # Useful for download_plugin where we want to pass in the request arguments to simulate a session.
        override_args = params.pop('session_args_override', None)

        # Instantiate an isolated settings object for this specific user session
        self.settings = Settings()

        super().__init__(**params)

        """
        DataFrame columns
        raw_value,spatial1,spatial2,n_genes_by_counts,UMAP1,UMAP2,clusters,clusters_cat_codes,colors
        """

        # Prioritize the injected override, fallback to the global WebSocket state
        args = override_args if override_args is not None else pn.state.session_args
        if args is not None:

            def get_arg(key, default=None):
                if key in args:
                    value = args[key][0].decode('utf-8')

                    # type coercion on some common keywords due to the stringification of the URL params.
                    if value.lower() in ('none', 'null', ''):
                        return None
                    if value.lower() == 'true':
                        return True
                    if value.lower() == 'false':
                        return False
                    return value

                return default

            # Decode the bytes back to standard Python strings
            # Format is {'filename': [b'my_data.parquet']}
            if 'dataset_id' in args:
                self.settings.dataset_id = get_arg('dataset_id', "")

            if 'filename' in args:
                self.settings.filename = get_arg('filename', "")

            if "hide_zeros" in args:
                hide_zeros = get_arg("hide_zeros", False)
                if hide_zeros is not None:
                    self.settings.hide_zeros = bool(int(hide_zeros))
                else:
                    self.settings.hide_zeros = False

            if "marker_shape" in args:
                marker_shape = get_arg("marker_shape", "square")
                if marker_shape is not None:
                    self.settings.marker_shape = str(marker_shape)
                else:
                    self.settings.marker_shape = "square"

            if 'selection_x1' in args:
                selection_x1 = get_arg('selection_x1')
                if selection_x1 is not None:
                    self.settings.selection_x1 = float(selection_x1)
            if 'selection_x2' in args:
                selection_x2 = get_arg('selection_x2')
                if selection_x2 is not None:
                    self.settings.selection_x2 = float(selection_x2)
            if 'selection_y1' in args:
                selection_y1 = get_arg('selection_y1')
                if selection_y1 is not None:
                    self.settings.selection_y1 = float(selection_y1)
            if 'selection_y2' in args:
                selection_y2 = get_arg('selection_y2')
                if selection_y2 is not None:
                    self.settings.selection_y2 = float(selection_y2)

            if 'expression_min_clip' in args:
                expression_min_clip = get_arg('expression_min_clip')
                if expression_min_clip is not None:
                    self.settings.expression_min_clip = float(expression_min_clip)

            if 'nosave' in args:
                nosave = get_arg('nosave', True)
                if nosave is not None:
                    self.settings.nosave = bool(int(nosave))
                else:
                    self.settings.nosave = False

            if 'display_name' in args:
                self.settings.display_name = get_arg('display_name', "")

            if 'make_default' in args:
                make_default = get_arg('make_default', False)
                if make_default is not None:
                    self.settings.make_default = bool(int(make_default))
                else:
                    self.settings.make_default = False

        # Add a hidden widget to sync url parameter state to the JS context.
        # This is so we can pass the params to the download widget so it can recreate the plot.
        # But first, give it a dummy string so the browser definitely recognizes a change later
        self.state_sync = pn.widgets.TextInput(value="INITIALIZING", visible=False)
        self.state_sync.jscallback(value="""
            if (cb_obj.value && cb_obj.value !== "INITIALIZING") {
                // Parse the updated state from Python
                const state = JSON.parse(cb_obj.value);

                // Construct a fresh URLSearchParams object
                const params = new URLSearchParams();
                for (const [key, val] of Object.entries(state)) {
                    params.append(key, val);
                }

                // Expose it globally to the frontend DOM
                window.gearSpatialUrlParams = params;
            }
        """)

        # Do an initial call to set state
        pn.state.onload(self._sync_state_to_js)

        self.orig_df = retrieve_dataframe(self.settings.dataset_id, self.settings.filename)
        self.orig_df['clusters'] = self.orig_df['clusters'].astype('category')
        self.orig_df = clip_expression_values(self.orig_df, self.settings.expression_min_clip)

        self.image_array = retrieve_image_array(self.settings.dataset_id)

        self.current_gene = normalize_expression_name(self.settings.filename)

        self.nosave = self.settings.nosave

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

        self.bounds_stream_image = hv.streams.BoundsXY(bounds=self.saved_bounds)  # type: ignore
        self.bounds_stream_composite = hv.streams.BoundsXY(bounds=self.saved_bounds)  # type: ignore

        # Track current viewport to retain zoom during widget interactions
        self.viewport_stream = hv.streams.RangeXY()

        # Set up a callback to update the URL params whenever the user draws or clears a box
        self.bounds_stream_image.add_subscriber(self._update_bounds_callback)
        self.bounds_stream_composite.add_subscriber(self._update_bounds_callback)

        # Set up the background image and its dimmed version for overlaying on the plots
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
                        tools=[], active_tools=[], default_tools=[], toolbar="below",
                        hooks=[autohide_toolbar]
                    )

            # Create a more opaque version of the image.
            opacity_value = int(255 * 0.5)  # Adjust the multiplier to set the desired opacity level
            alpha_channel = np.full((self.img_height, self.img_width, 1), opacity_value, dtype=np.uint8)

            # Check if image is RGB (3 channels). If so, append alpha. If already RGBA, just overwrite alpha.
            if self.image_array.shape[-1] == 3:
                dimmed_array = np.concatenate([self.image_array, alpha_channel], axis=-1)
            else:
                dimmed_array = self.image_array.copy()

                # SAdkins note - '...' means all rows in this case.
                dimmed_array[..., 3] = opacity_value

            self.bg_image_dimmed = hv.RGB(dimmed_array, bounds=img_bounds).opts(
                        xaxis=None, yaxis=None, responsive=True,
                        # No tools or toolbar as the composite plot takes care of this
                        tools=[], active_tools=[], default_tools=[],
                        hooks=[autohide_toolbar]
                    )

        self.df = self.orig_df.copy()

        # One unfortunately annoyance is that datashader's default behavior is to flip the y-axis,
        # which is not what we want for spatial data. To fix this,
        # we can reverse the y-axis limits by setting ylim to (max, min) instead of (min, max).
        self.df["y_plot"] = self.df["spatial2"]
        if self.img_height is not None:
            self.df["y_plot"] = self.img_height - self.df["spatial2"]

        self._init_widgets()


        ### Set up some attributes to use when plotting
        self.shape = "square"

        # Precompute the datashader aggregations since they can be shared across multiple plots
        self.expression_agg = ds.mean('raw_value')
        self.expression_cmap = 'YlOrRd'

        self.clusters_agg = ds.count_cat('clusters')
        self.cluster_cmap = dict(zip(self.df['clusters'], self.df['colors']))


    def _sync_state_to_js(self):
        """Serializes current active settings and pipes them to the frontend DOM."""

        hide_zeros = self.settings.hide_zeros
        if hasattr(self, "hide_zeros") and self.hide_zeros is not None:
            hide_zeros = self.hide_zeros

        marker_shape = self.settings.marker_shape
        if hasattr(self, "marker_shape") and self.marker_shape is not None:
            marker_shape = self.marker_shape

        state = {
            "dataset_id": self.settings.dataset_id,
            "filename": self.settings.filename,
            "hide_zeros": hide_zeros,  # Use the active class property
            "marker_shape": marker_shape,  # Use the active class property
        }

        # Include expression clip if it exists
        if self.settings.expression_min_clip is not None:
            state["expression_min_clip"] = self.settings.expression_min_clip

        # Include spatial bounds if an active selection exists
        if self.settings.selection_x1 is not None:
            state["selection_x1"] = self.settings.selection_x1
            state["selection_x2"] = self.settings.selection_x2
            state["selection_y1"] = self.settings.selection_y1
            state["selection_y2"] = self.settings.selection_y2

        # Writing to this value triggers the JS callback instantly
        self.state_sync.value = json.dumps(state)

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
            self.saved_bounds = None
        else:
            # User drew a box
            (self.settings.selection_x1, self.settings.selection_y1, self.settings.selection_x2, self.settings.selection_y2) = bounds
            self.saved_bounds = bounds

        # Trigger the zoom update!
        self._update_zoom_panel(bounds)

        # BROADCAST TO FRONTEND
        self._sync_state_to_js()

    def _create_ghost_legend(self):
            """Creates a fake, invisible plot just to force Bokeh to draw a legend."""

            # TODO: add some "click" interactivity to hide plot elements based on if legend item is enabled or disabled.

            # Hook to strip the invisible padding Bokeh reserves for axes
            def remove_bokeh_borders(plot, element):
                plot.state.min_border = 0
                plot.state.min_border_left = 0
                plot.state.min_border_right = 0
                plot.state.min_border_top = 0
                plot.state.min_border_bottom = 0

            ghost_points = []
            # self.cluster_cmap should be a dict like {'Cluster 1': '#FF0000', ...}

            sorted_clusters = sort_clusters(self.cluster_cmap.keys())
            for cluster_name in sorted_clusters:
                hex_color = self.cluster_cmap[cluster_name]

                # Create a single point at (NaN, NaN) so it doesn't draw on the screen
                # But give it a label and a color so the legend picks it up
                pt = hv.Points(
                    [(np.nan, np.nan)],
                    label=str(cluster_name)
                ).opts(color=hex_color, size=10, tools=[], default_tools=[], hooks=[remove_bokeh_borders])

                ghost_points.append(pt)

            # Combine all the ghost points into a single overlay
            return hv.Overlay(ghost_points)

    def _generate_spatial_grid(self):
        """
        This is where you would implement the logic to generate the spatial grid layout with the background image, expression plot, and cluster plot.
        Row 1 is the main view, and Row 2 (optional) is the zoomed-in view. The legend is also included in Row 1.
        """
        raise NotImplementedError("Subclasses must implement _generate_spatial_grid")

    def _init_widgets(self):
        """
        This is where you would initialize any Panel widgets (sliders, dropdowns, etc.) that you want to use in your app.
        You can then reference these widgets in your _build_layout method to include them in the layout and set up callbacks.
        """

        self.hide_zeros_toggle = pn.widgets.Checkbox(
                name='Hide Zero Expression Observations',
                value=False,
                margin=(10, 10)
            )

        self.marker_shape_toggle = pn.widgets.Select(
            options={'Square': 'square', 'Circle': 'circle'},
            value='square',
            width=100,
            align="center"
        )

        self.marker_shape_ui = pn.Row(
            pn.widgets.StaticText(value='Marker Shape:', align='start', margin=(10, 10)),
            self.marker_shape_toggle
        )

        # Bind the master update function to the exact dropdown widget, not the container
        self.hide_zeros_toggle.param.watch(self._update_plots, 'value')
        self.marker_shape_toggle.param.watch(self._update_plots, 'value')

    def _build_layout(self):
        """
        This is where the main Panel layout is built.
        We leave it blank in the base class since the Condensed and Expanded viewers will have different layouts.
        """
        raise NotImplementedError("Subclasses must implement _build_layout")

    def _update_plots(self, event):
        """
        This is where you would implement the logic to update the plots based on widget changes.
        The event parameter will contain information about which widget changed and its new value.
        You can use this information to filter your dataframe, update plot parameters, or trigger a redraw of the plots.
        """
        raise NotImplementedError("Subclasses must implement _update_plots")

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

            # Generate base plots
            self.main_row = self._generate_spatial_grid()

            # Return final Panel layout
            return pn.Column(
                self.state_sync, # Invisible DOM injector
                self.pre_layout,
                self.main_row,
                sizing_mode='stretch_both', # Fills the 100%x100% iframe
                margin=(0, 0, 0, 0)
                )
        except Exception as e:
            traceback.format_exc()
            return pn.pane.Alert(f"Error: {e}", alert_type="danger")

    def _generate_spatial_grid(self):
        """
        Generates the spatial grid layout with the background image, expression plot, and cluster plot.
        Row 1 is the main view, and Row 2 is the zoomed-in view. The legend is also included in Row 1.
        """

        if getattr(self, 'bg_image', None) is not None:
            self.bg_image = self.bg_image.opts(default_tools = ["box_zoom", "wheel_zoom", "pan", "reset"])

        # Drop the rows completely for maximum Datashader stability
        if hasattr(self, 'hide_zeros_toggle') and self.hide_zeros_toggle.value:
            spatial_df = self.df[self.df['raw_value'] > 0]
        else:
            spatial_df = self.df

        # Generate base plots
        expr_plot = create_spatial_plot(spatial_df, self.expression_agg, y_col="y_plot", color_col='raw_value', cmap=self.expression_cmap, title=f"Expression: {self.current_gene}", mode="standard", shape=self.shape)
        cluster_plot = create_spatial_plot(spatial_df, self.clusters_agg, y_col="y_plot", color_col='clusters', cmap=self.cluster_cmap, is_categorical=True, title="Clusters", mode="standard", shape=self.shape) # type: ignore

        # Capture the exact viewport limits right before the redraw
        current_x = self.viewport_stream.x_range
        current_y = self.viewport_stream.y_range

        # If the user is currently zoomed/panned, apply those limits to the new plots
        if current_x and current_y:
            expr_plot = expr_plot.opts(xlim=current_x, ylim=current_y)
            cluster_plot = cluster_plot.opts(xlim=current_x, ylim=current_y)
            if hasattr(self, 'bg_image') and getattr(self, 'bg_image', None) is not None:
                self.bg_image = self.bg_image.opts(xlim=current_x, ylim=current_y)

        # Attach the stream to the new expression plot to track future interactions
        self.viewport_stream.source = expr_plot

        master_image = self.bg_image if getattr(self, 'bg_image', None) is not None else pn.pane.HoloViews(None)

        # Create composites with the background image
        if hasattr(self, 'bg_image_dimmed') and self.bg_image_dimmed is not None:
            main_expr = self.bg_image_dimmed * expr_plot
            main_cluster = self.bg_image_dimmed * cluster_plot
        else:
            main_expr = expr_plot
            main_cluster = cluster_plot

        # Store master composites for the zoom callback
        self.master_image = master_image

        # Create the Ghost Legend and Container (Row 1 only)
        needed_width = max(180, max(len(str(name)) for name in self.cluster_cmap) * 7)
        needed_height = max(300, len(self.cluster_cmap) * 22 + 50)

        ghost_legend = self._create_ghost_legend().opts(
            show_legend=True, legend_position="top_left",
            xaxis=None, yaxis=None, show_frame=False, toolbar=None,
            width=needed_width, height=needed_height
        )

        legend_container = pn.Column(
            pn.pane.HoloViews(ghost_legend),
            width=200, height=300, scroll=True, margin=(0,0,0,0)
        )

        # This has a min_width to prevent plots from being squished in smaller layout tile configurations.
        return pn.Row(
            master_image, main_expr, main_cluster, legend_container,
            sizing_mode='stretch_both',
            min_width=1080,
            margin=(0, 0, 0, 0)
        )

    def _init_widgets(self):
        """Initializes widgets and builds the condensed control bar."""
        super()._init_widgets()

        # Build the top control bar
        self.pre_layout = pn.Row(
            self.hide_zeros_toggle,
            pn.Spacer(width=50),
            self.marker_shape_ui,
            sizing_mode='stretch_width',
            margin=(0, 0, 10, 0)
        )

    async def _update_plots(self, event):
        """Triggered when a widget change occurs, allowing for hot-swapping."""

        self.hide_zeros = self.hide_zeros_toggle.value
        self.marker_shape = self.marker_shape_toggle.value

        if self.marker_shape == "circle":
            self.shape="circle"
        else:
            self.shape="square"

        # Start event loading spinners on the affected containers
        self.main_row.loading = True

        # Yield control to the event loop to give the browser time to render spinners
        await asyncio.sleep(0.05)

        try:

            # Rebuild the Spatial Grid
            new_spatial_grid = self._generate_spatial_grid()
            self.main_row[:] = new_spatial_grid[:]
        except Exception as e:
            # Catch and print any hidden Datashader or Panel errors to your terminal
            print(f"Error updating plots: {e}")
            traceback.print_exc()
        finally:
            # Turn off the spinners and sweep memory
            self.main_row.loading = False

            # BROADCAST TO FRONTEND
            self._sync_state_to_js()

            gc.collect()


    def _update_zoom_panel(self, bounds):
        """
        Currently not needed.
        """
        pass


class ExpandedSpatialViewer(BaseSpatialViewer):
    """
    The specific component for your panel_app_expanded.py
    """

    def _build_layout(self):
        # Build your spatial rows, UMAPs, and Violins here...
        try:

            # Represents main row and zoom row
            self.spatial_grid_container = self._generate_spatial_grid()

            ### UMAP row
            expr_umap = create_umap_plot(
                self.df, self.expression_agg, color_col='raw_value', cmap="cividis_r", is_categorical=False, title=f"{self.current_gene} Expression"
            )
            cluster_umap = create_umap_plot(
                self.df, self.clusters_agg, color_col='clusters', cmap=self.cluster_cmap, is_categorical=True, title="Clusters"
            )

            # Create the UMAP legend
            needed_width = max(180, max(len(str(name)) for name in self.cluster_cmap) * 7)
            needed_height = max(300, len(self.cluster_cmap) * 22 + 50)

            umap_ghost_legend = self._create_ghost_legend().opts(
                show_legend=True,
                legend_position="top_left",
                xaxis=None, yaxis=None,
                show_frame=False, toolbar=None,
                width=needed_width,
                height=needed_height
            )

            umap_legend_container = pn.Column(
                umap_ghost_legend,
                width=200,
                height=300, # Match UMAP min_height
                scroll=True,
                margin=(0, 0, 0, 0)
            )

            # Layout the UMAP row
            self.umap_row_container = pn.Row(
                expr_umap,
                cluster_umap,
                umap_legend_container,
                sizing_mode='stretch_width'
            )

            ### Violin
            violin_base = create_violin_plot(
                self.df,
                y_col='raw_value',
                group_col='clusters',
                cmap=self.cluster_cmap,
                title=f"Expression Distribution: {self.current_gene} by Cluster"
            )

            self.violin_row_container = pn.Row(
                violin_base,
                sizing_mode='stretch_width'
            )

            return pn.Column(
                self.state_sync, # Invisible DOM injector
                self.pre_layout,
                self.spatial_grid_container,
                pn.layout.Divider(margin=(20, 0)), # Visual breathing room
                self.umap_row_container,
                pn.layout.Divider(margin=(20, 0)),
                self.violin_row_container,
                sizing_mode='stretch_both'
            )
        except Exception as e:
            traceback.format_exc()
            return pn.pane.Alert(f"Error: {e}", alert_type="danger")

    def _init_widgets(self):
        # ? This can be useful for filtering datasets even for projections, but how best to word it?

        super()._init_widgets()

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
        # Rebuild the pre_layout with the new toggles on the bottom row
        self.pre_layout = pn.Column(
            pn.Row(
                pn.pane.Markdown(
                    '## Use the "box select" tool on the top row to update the zoomed view.',
                    height=30,
                    width=markdown_width,
                ),
                self.display_name,
                self.save_button,
            ),
            pn.Row(
                self.hide_zeros_toggle,
                pn.Spacer(width=50),
                self.marker_shape_ui,
                pn.Spacer(), # Auto-fills space between toggles and the save checkbox
                self.make_default,
                sizing_mode='stretch_width'
            ),
        )

        # Emit a DOM event instead of changing a URL param
        # Pass the active states of the widgets down to the JS context (some things, like gene_symbol are already in the frontend)
        self.save_button.js_on_click(
            args={
                'name_input': self.display_name,
                'default_cb': self.make_default,
                'hide_zeros_toggle': self.hide_zeros_toggle,
                'marker_shape_toggle': self.marker_shape_toggle,
                'dataset_id': self.settings.dataset_id
            },
            code=f"""
            // Create a custom event containing all the widget values
            const evt = new CustomEvent(`save_spatial_display_${{dataset_id}}`, {{
                detail: {{
                    displayName: name_input.value,
                    makeDefault: default_cb.active,
                    hideZeros: hide_zeros_toggle.active,
                    marker_shape: marker_shape_toggle.active
                    // If you tracked bounds in JS, you could grab them, or grab them from Panel
                }}
            }});

            // Dispatch it to the host window so the vanilla JS listener catches it instantly
            window.dispatchEvent(evt);
            """)  # noqa: F541

    async def _update_plots(self, event):
        """Triggered when a widget change occurs, allowing for hot-swapping."""

        self.hide_zeros = self.hide_zeros_toggle.value
        self.marker_shape = self.marker_shape_toggle.value

        if self.marker_shape == "circle":
            self.shape="circle"
        else:
            self.shape="square"

        # Start event loading spinners on the affected containers
        self.spatial_grid_container.loading = True

        # Yield control to the event loop to give the browser time to render spinners
        await asyncio.sleep(0.05)

        try:
            # Rebuild the Spatial Grid
            new_spatial_grid = self._generate_spatial_grid()
            self.spatial_grid_container[:] = new_spatial_grid[:]

            # Re-trigger zoom if a box is currently drawn
            if self.saved_bounds is not None:
                self._update_zoom_panel(self.saved_bounds)
        except Exception as e:
            # Catch and print any hidden Datashader or Panel errors to your terminal
            print(f"Error updating plots: {e}")
            traceback.print_exc()
        finally:
            # Turn off the spinners and sweep memory
            self.spatial_grid_container.loading = False

            # BROADCAST TO FRONTEND
            self._sync_state_to_js()

            gc.collect()


    def _update_zoom_panel(self, bounds):
        """
        This method updates the row of zoomed-in plots based on the provided bounds.
        The method applies the new limits to the master composite plots and updates the zoom row accordingly.

        Parameters:
        - bounds: A tuple containing the new bounds in the format (left, bottom, right, top).
                  If bounds is None, it indicates that the user has cleared the selection.
        """

        if not hasattr(self, 'zoom_master_expr'):
            return

        # 1. Determine bounds
        if bounds is None:
            opts_dict = dict(xlim=(None, None), ylim=(None, None), clone=True)
        else:
            x1, y1, x2, y2 = bounds
            opts_dict = dict(
                xlim=(min(x1, x2), max(x1, x2)),
                ylim=(min(y1, y2), max(y1, y2)),
                clone=True,
                active_tools=['pan', 'wheel_zoom'] # Disable box select on zoomed plots
            )

        # 2. Apply limits to all three master canvases
        new_zoom_image = self.master_image.opts(**opts_dict) if self.master_image is not None else None
        new_zoom_expr = self.zoom_master_expr.opts(**opts_dict)
        new_zoom_cluster = self.zoom_master_cluster.opts(**opts_dict)

        # 3. Swing the Sledgehammer (Hot-swap all three panes)
        self.zoom_image_pane = pn.pane.HoloViews(new_zoom_image, sizing_mode='stretch_width', linked_axes=False)
        self.zoom_expr_pane = pn.pane.HoloViews(new_zoom_expr, sizing_mode='stretch_width', linked_axes=False)
        self.zoom_cluster_pane = pn.pane.HoloViews(new_zoom_cluster, sizing_mode='stretch_width', linked_axes=False)

        # 4. Inject back into Row 2 (preserving the spacer at index 3)
        if hasattr(self, 'zoom_row'):
            self.zoom_row[0] = self.zoom_image_pane
            self.zoom_row[1] = self.zoom_expr_pane
            self.zoom_row[2] = self.zoom_cluster_pane

    def _generate_spatial_grid(self):
        """
        Generates the spatial grid layout with the background image, expression plot, and cluster plot.
        Row 1 is the main view, and Row 2 is the zoomed-in view. The legend is also included in Row 1.
        """

        if getattr(self, 'bg_image', None) is not None:
            self.bg_image = self.bg_image.opts(default_tools = ["box_select", "reset"])

        # Drop the rows completely for maximum Datashader stability
        if hasattr(self, 'hide_zeros_toggle') and self.hide_zeros_toggle.value:
            spatial_df = self.df[self.df['raw_value'] > 0]
        else:
            spatial_df = self.df

        # Generate base plots
        expr_plot = create_spatial_plot(spatial_df, self.expression_agg, y_col="y_plot", color_col='raw_value', cmap=self.expression_cmap, title=f"Expression: {self.current_gene}", mode="expanded", shape=self.shape)
        cluster_plot = create_spatial_plot(spatial_df, self.clusters_agg, y_col="y_plot", color_col='clusters', cmap=self.cluster_cmap, is_categorical=True, title="Clusters", mode="expanded", shape=self.shape) # type: ignore

        # Capture the exact viewport limits right before the redraw
        current_x = self.viewport_stream.x_range
        current_y = self.viewport_stream.y_range

        # If the user is currently zoomed/panned, apply those limits to the new plots
        if current_x and current_y:
            expr_plot = expr_plot.opts(xlim=current_x, ylim=current_y)
            cluster_plot = cluster_plot.opts(xlim=current_x, ylim=current_y)
            if hasattr(self, 'bg_image') and getattr(self, 'bg_image', None) is not None:
                self.bg_image = self.bg_image.opts(xlim=current_x, ylim=current_y)

        # Attach the stream to the new expression plot to track future interactions
        self.viewport_stream.source = expr_plot

        master_image = self.bg_image if getattr(self, 'bg_image', None) is not None else pn.pane.HoloViews(None)

        # Update bounds streams for the image and expression plots.
        self.bounds_stream_image.source = master_image
        self.bounds_stream_composite.source = expr_plot

        # Create a 3rd stream for the cluster plot so it can also trigger the zoom
        if not hasattr(self, 'bounds_stream_cluster'):
            self.bounds_stream_cluster = hv.streams.BoundsXY(bounds=self.saved_bounds)
            self.bounds_stream_cluster.add_subscriber(self._update_bounds_callback)
            self.bounds_stream_cluster.source = cluster_plot

        # Create composites with the background image
        if hasattr(self, 'bg_image_dimmed') and self.bg_image_dimmed is not None:
            # Row 1
            main_expr = self.bg_image_dimmed
            main_cluster = self.bg_image_dimmed
            # Row 2 Zoom
            zoom_master_expr = self.bg_image_dimmed * expr_plot
            zoom_master_cluster = self.bg_image_dimmed * cluster_plot
        else:
            main_expr = expr_plot
            main_cluster = cluster_plot
            zoom_master_expr = expr_plot
            zoom_master_cluster = cluster_plot

        # Store master composites for the zoom callback
        self.master_image = master_image
        self.zoom_master_expr = zoom_master_expr
        self.zoom_master_cluster = zoom_master_cluster

        # Create the Ghost Legend and Container (Row 1 only)
        needed_width = max(180, max(len(str(name)) for name in self.cluster_cmap) * 7)
        needed_height = max(300, len(self.cluster_cmap) * 22 + 50)

        ghost_legend = self._create_ghost_legend().opts(
            show_legend=True, legend_position="top_left",
            xaxis=None, yaxis=None, show_frame=False, toolbar=None,
            width=needed_width, height=needed_height
        )

        legend_container = pn.Column(
            pn.pane.HoloViews(ghost_legend),
            width=200, height=300, scroll=True, margin=(0,0,0,0)
        )

        # Initialize Zoom Panes
        self.zoom_image_pane = pn.pane.HoloViews(master_image, sizing_mode='stretch_width', linked_axes=False)
        self.zoom_expr_pane = pn.pane.HoloViews(zoom_master_expr, sizing_mode='stretch_width', linked_axes=False)
        self.zoom_cluster_pane = pn.pane.HoloViews(zoom_master_cluster, sizing_mode='stretch_width', linked_axes=False)

        # Assemble Rows (4 slots each to maintain vertical alignment)
        self.main_row = pn.Row(master_image, main_expr, main_cluster, legend_container, sizing_mode='stretch_width')

        # Use a spacer in Row 2 to match the legend width perfectly
        zoom_spacer = pn.Spacer(width=200, margin=(0,0,0,0))
        self.zoom_row = pn.Row(self.zoom_image_pane, self.zoom_expr_pane, self.zoom_cluster_pane, zoom_spacer, sizing_mode='stretch_width')

        zoom_markdown = pn.pane.Markdown('### Zoomed View', height=30)

        return pn.Column(self.main_row, zoom_markdown, self.zoom_row, sizing_mode='stretch_width')
