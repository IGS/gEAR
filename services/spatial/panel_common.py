import asyncio
import gc
import json
import traceback

import datashader as ds
import holoviews as hv
import hvplot
import hvplot.pandas  # noqa
import numpy as np
import panel as pn
from bokeh.models import CustomJS
from common import (
    Settings,
    autohide_toolbar,
    clip_expression_values,
    compute_aggregation_params,
    create_datashader_agg,
    create_clusters_df,
    create_expression_df,
    create_spatial_plot,
    create_umap_plot,
    create_violin_plot,
    link_crosshairs,
    link_ranges,
    normalize_expression_name,
    retrieve_dataframe,
    retrieve_image_array,
    sort_clusters,
)

# CRITICAL: Initialize the Bokeh backend for interactivity
hvplot.extension('bokeh', logo=False) # type: ignore
pn.extension(loading_indicator=True, defer_load=True, nthreads=4)  # type: ignore)

MAX_CARD_WIDTH = 1800
EXPANDED_PLOT_MIN_HEIGHT = 380
EXPANDED_PLOT_MIN_WIDTH = 320

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

            # Preserve view range carried forward from a previous gene's
            # session. Switching genes tears down this
            # whole Panel/Bokeh server session and spins up a fresh one, so
            # we apply this as the initial range.
            def _get_restore_range(start_key, end_key):
                start_raw, end_raw = get_arg(start_key), get_arg(end_key)
                if start_raw is None or end_raw is None:
                    return None
                try:
                    return (float(start_raw), float(end_raw))
                except (TypeError, ValueError):
                    return None

            self.restore_xlim = _get_restore_range('x_range_start', 'x_range_end')
            self.restore_ylim = _get_restore_range('y_range_start', 'y_range_end')
        else:
            self.restore_xlim = None
            self.restore_ylim = None

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

                // Sync zoom range to frontend so switching genes
                // (which spins up a brand new Panel session) can carry
                // the last zoom/pan position forward into the next one.
                if (state.dataset_id &&
                    state.x_range_start !== undefined && state.x_range_end !== undefined &&
                    state.y_range_start !== undefined && state.y_range_end !== undefined) {
                    window.gearSpatialViewState = window.gearSpatialViewState || {};
                    window.gearSpatialViewState[state.dataset_id] = {
                        x_range_start: state.x_range_start,
                        x_range_end: state.x_range_end,
                        y_range_start: state.y_range_start,
                        y_range_end: state.y_range_end,
                    };
                }

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
            opacity_value = int(255 * 0.3)  # Adjust the multiplier to set the desired opacity level
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
                        tools=[], default_tools=[],
                        hooks=[autohide_toolbar],
                    )

        self.df = self.orig_df.copy()

        # One unfortunately annoyance is that datashader's default behavior is to flip the y-axis,
        # which is not what we want for spatial data. To fix this,
        # we can reverse the y-axis limits by setting ylim to (max, min) instead of (min, max).
        self.df["y_plot"] = self.df["spatial2"]
        if self.img_height is not None:
            self.df["y_plot"] = self.img_height - self.df["spatial2"]

        self.expr_df = self.df[self.df['raw_value'] > 0]

        # Create a mapping of cluster cat codes to names to re-add after datashader aggregation
        self.cluster_map = {
            code: self.df[self.df["clusters_cat_codes"] == code]["clusters"].to_numpy()[0]
            for code in self.df["clusters_cat_codes"].unique()
        }

        # Get 98th-percentile of valid expression for clipping the color scale. This is to avoid outliers dominating the color mapping.
        # Using all data runs into the issue of the value being 0.
        self.expression_98 = self.expr_df["raw_value"].quantile(0.98)
        if self.expression_98 == 0:
            self.expression_98 = None

        ### Set up some initial attributes to use when plotting
        self._init_widgets()

        # Precompute the datashader aggregations since they can be shared across multiple plots
        # For some reason, the x and y values are swapped for aggregation compared to what we eventually plot, so we need to swap them here.
        agg_width, agg_height, self.marker_radius = compute_aggregation_params(self.df, x_col="y_plot", y_col="spatial1")
        self.agg = create_datashader_agg(self.df, x="y_plot", y="spatial1", width=agg_width, height=agg_height)

        self.expression_agg = self.agg["expression"]
        self.umap_expression_agg = ds.max('raw_value')
        self.expression_df = create_expression_df(self.expression_agg)
        self.expression_df = self.expression_df[self.expression_df['raw_value'] > 0]
        self.expression_cmap = 'fire_r'

        self.clusters_agg = self.agg["clusters"]
        self.umap_clusters_agg = ds.count_cat('clusters')
        self.clusters_df = create_clusters_df(self.clusters_agg)
        self.clusters_df["clusters"] = self.clusters_df["clusters_cat_codes"].map(self.cluster_map)
        self.cluster_cmap = dict(zip(self.df['clusters'], self.df['colors']))


    def _sync_state_to_js(self):
        """Serializes current active settings and pipes them to the frontend DOM."""

        state = {
            "dataset_id": self.settings.dataset_id,
            "filename": self.settings.filename,
        }

        # Include expression clip if it exists
        if self.settings.expression_min_clip is not None:
            state["expression_min_clip"] = self.settings.expression_min_clip

        # Include the current view range, if zoomed in
        range_fig = getattr(self, '_range_figure', None)
        if range_fig is not None:
            try:
                if range_fig.x_range.start is not None \
                    and range_fig.x_range.end is not None \
                    and range_fig.y_range.start is not None \
                    and range_fig.y_range.end is not None:
                    state["x_range_start"] = float(range_fig.x_range.start)
                    state["x_range_end"] = float(range_fig.x_range.end)
                    state["y_range_start"] = float(range_fig.y_range.start)
                    state["y_range_end"] = float(range_fig.y_range.end)
            except (TypeError, ValueError):
                pass

        # Writing to this value triggers the JS callback instantly
        self.state_sync.value = json.dumps(state)

    # NOTE: We only need to watch one figure, as the range will propagate to the others via the link_ranges function
    def _apply_restored_range(self, figure):
        """
        Applies self.restore_xlim/self.restore_ylim (parsed in __init__ from
        the x_range_start/end, y_range_start/end URL args) as the initial
        view range on the provided Bokeh figure, if
        present. Call this before _wire_range_capture so we don't immediately
        "capture" our own restored range as a spurious change.
        """
        if figure is None:
            return
        if self.restore_xlim is not None and self.restore_ylim is not None:
            figure.x_range.start, figure.x_range.end = self.restore_xlim
            figure.y_range.start, figure.y_range.end = self.restore_ylim

    def _sync_cluster_legend_toggle(self, ghost_fig, real_renderers_by_cat):
        """
        Rewires the standalone (already-rendered) ghost legend's Bokeh Legend
        so clicking a swatch toggles visibility of the cluster-plot
        renderer, not just the invisible ghost one -- recreating the
        Plotly-style "click to hide" behavior
        """
        if not ghost_fig.legend:
            return

        # Map each legend item to the renderer for the Bokeh cluster spatial plot
        legend = ghost_fig.legend[0]
        for item in legend.items:
            label = item.label
            if isinstance(label, dict):
                label = label.get('value') or label.get('field')
            elif hasattr(label, 'value'):
                label = label.value

            real_renderer = real_renderers_by_cat.get(label)
            if real_renderer is not None:
                item.renderers = list(item.renderers) + [real_renderer]

        # toggles visibility
        legend.click_policy = 'hide'


    def _sync_range_updates(self, figure):
        """
        Watches the figure's x_range/y_range for live changes (pan, box-zoom,
        wheel-zoom) and pushes the current range up to the frontend via the
        state_sync -> window.gearSpatialViewState bridge, so switching genes
        can restore the same view. Bokeh Server automatically syncs client-
        side property changes back to these same server-side Range1d objects,
        so a plain `.on_change` here is standard Bokeh behavior.

        The spatial figures aren't linked to pan/zoom together today, so this
        just tracks whichever one was interacted with most recently as the
        "current" view to restore.
        """
        def make_callback(fig):
            def on_range_change(attr, old, new):
                self._range_figure = fig
                self._sync_state_to_js()
            return on_range_change

        if figure is None:
            return
        cb = make_callback(figure)
        figure.x_range.on_change('start', cb)
        figure.x_range.on_change('end', cb)
        figure.y_range.on_change('start', cb)
        figure.y_range.on_change('end', cb)


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

    def _force_legend_redraw_on_range_change(self, ghost_fig, range_source_fig):
        """
        Workaround for an observed Bokeh rendering quirk: after panning/
        zooming the (range-linked) spatial figures, the ghost legend's swatch
        icons can stop being painted, and only reliably come back once a
        legend item is actually clicked -- which apparently forces Bokeh to
        fully redraw the legend as a side effect of the click_policy='hide'
        interaction, regardless of which item was clicked or whether it was
        being hidden or shown.

        This forces that same redraw proactively, client-side only, whenever
        the range changes, instead of waiting for a click: briefly toggling
        the legend's own `visible` property off and back on forces Bokeh to
        recompute and repaint it, the same way toggling a renderer's
        `visible` via click_policy does.

        NOTE: SAdkins - Still don't fully understand the bug.
        """
        if not ghost_fig.legend:
            return

        legend = ghost_fig.legend[0]
        force_redraw = CustomJS(args=dict(legend=legend), code="""
            legend.visible = false;
            requestAnimationFrame(() => { legend.visible = true; });
        """)
        range_source_fig.x_range.js_on_change('start', force_redraw)
        range_source_fig.x_range.js_on_change('end', force_redraw)
        range_source_fig.y_range.js_on_change('start', force_redraw)
        range_source_fig.y_range.js_on_change('end', force_redraw)


    def _generate_spatial_plots(self):
        """
        This is where you would implement the logic to generate the plots to be returned for the spatial grid layout. Includes the background image, expression plot, and cluster plot.
        """
        raise NotImplementedError("Subclasses must implement _generate_spatial_plots")

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

    def _update_plots(self, event):
        """
        This is where you would implement the logic to update the plots based on widget changes.
        The event parameter will contain information about which widget changed and its new value.
        You can use this information to filter your dataframe, update plot parameters, or trigger a redraw of the plots.
        """
        raise NotImplementedError("Subclasses must implement _update_plots")

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
            # Generate raw plots
            master_image, main_expr, main_cluster, category_renderers = self._generate_spatial_plots()

            # Render to concrete Bokeh figures (instead of pn.pane.HoloViews), so we can synchronize some things across each
            bokeh_expr = hv.render(main_expr, backend='bokeh')
            bokeh_cluster = hv.render(main_cluster, backend='bokeh')
            figures_to_link = [bokeh_expr, bokeh_cluster]

            # If a master image exists, render it to Bokeh and add it to the list of figures to link
            if master_image is not None:
                bokeh_master = hv.render(master_image, backend='bokeh')
                figures_to_link.insert(0, bokeh_master)
                self.master_pane = pn.pane.Bokeh(bokeh_master, sizing_mode='stretch_width')
            else:
                self.master_pane = pn.pane.HoloViews(master_image, sizing_mode='stretch_width')

            link_ranges(figures_to_link)
            link_crosshairs(figures_to_link)

            # Restore a carried-over view range (from switching genes) before
            # wiring up capture, so we don't immediately re-report our own
            # restored range as a fresh "change".
            self._apply_restored_range(bokeh_expr)
            self._sync_range_updates(bokeh_expr)

            # Worth noting the tradeoff in using the Bokeh pane instead of Holoviews is that we lose
            # some reactivity, but gain access to the underlying Bokeh figure.
            # Hence, why some things like link_ranges has be implemented manually.
            self.expr_pane = pn.pane.Bokeh(bokeh_expr, sizing_mode='stretch_width')
            self.cluster_pane = pn.pane.Bokeh(bokeh_cluster, sizing_mode='stretch_width')

            # Create standalone legend
            needed_width = max(180, max(len(str(name)) for name in self.cluster_cmap) * 7)
            needed_height = max(300, len(self.cluster_cmap) * 22 + 50)
            ghost_legend = self._create_ghost_legend().opts(
                show_legend=True, legend_position="top_left", xaxis=None, yaxis=None, show_frame=False, toolbar=None, width=needed_width, height=needed_height
            )

            bokeh_ghost = hv.render(ghost_legend, backend='bokeh')
            self._sync_cluster_legend_toggle(bokeh_ghost, category_renderers)
            self._force_legend_redraw_on_range_change(bokeh_ghost, bokeh_expr)
            self.legend_pane = pn.Column(pn.pane.Bokeh(bokeh_ghost), width=200, height=300, scroll=True, margin=(0,0,0,0))

            # Build the permanent layout row
            self.main_row = pn.Row(
                self.master_pane, self.expr_pane, self.cluster_pane, self.legend_pane,
                sizing_mode='stretch_both', min_width=1080, margin=(0, 0, 0, 0)
            )

            return pn.Column(self.state_sync, self.main_row, sizing_mode='stretch_both', max_width=MAX_CARD_WIDTH, margin=(0, 0, 0, 0))

        except Exception as e:
            traceback.print_exc()
            return pn.pane.Alert(f"Error: {e}", alert_type="danger")

    def _generate_spatial_plots(self):
        """Generates the raw HoloViews objects without wrapping them in Panel layouts."""
        if getattr(self, 'bg_image', None) is not None:
            self.bg_image = self.bg_image.opts(default_tools=["box_zoom", "wheel_zoom", "pan", "reset"], active_tools=["wheel_zoom", "pan"], min_width=275, min_height=300)

        # This gets filled in during the cluster spatial plot creation
        # It is a registry of category to renderer, which can then be
        # used to wire up the ghost legend to toggle visibility of the real cluster plot.
        category_renderers = {}

        # Generate base plots
        expr_plot = create_spatial_plot(self.expression_df, y_col="y_plot", color_col='raw_value',
                                        cmap=self.expression_cmap, title=f"Expression: {self.current_gene}",
                                        cbar_max=self.expression_98, radius=self.marker_radius)
        cluster_plot = create_spatial_plot(self.clusters_df, y_col="y_plot", color_col='clusters',
                                        cmap=self.cluster_cmap, is_categorical=True, title="Clusters",
                                        category_renderers=category_renderers, radius=self.marker_radius)

        if hasattr(self, 'bg_image_dimmed') and self.bg_image_dimmed is not None:
            main_expr = self.bg_image_dimmed * expr_plot
            main_cluster = self.bg_image_dimmed * cluster_plot
        else:
            main_expr = expr_plot
            main_cluster = cluster_plot

        return self.bg_image, main_expr, main_cluster, category_renderers

    def _init_widgets(self):
        """Initializes widgets and builds the condensed control bar."""
        pass


class ExpandedSpatialViewer(BaseSpatialViewer):
    """
    The specific component for your panel_app_expanded.py
    """

    def _build_layout(self):
        try:
            # Generate Raw Plots
            master_image, main_expr, main_cluster, category_renderers = self._generate_spatial_plots()

            # Render to concrete Bokeh figures (instead of pn.pane.HoloViews), so we can synchronize some things across each
            bokeh_expr = hv.render(main_expr, backend='bokeh')
            bokeh_cluster = hv.render(main_cluster, backend='bokeh')
            figures_to_link = [bokeh_expr, bokeh_cluster]

            # If a master image exists, render it to Bokeh and add it to the list of figures to link
            if master_image is not None:
                bokeh_master = hv.render(master_image, backend='bokeh')
                figures_to_link.insert(0, bokeh_master)
                self.master_pane = pn.pane.Bokeh(bokeh_master, sizing_mode='stretch_width')
            else:
                self.master_pane = pn.pane.HoloViews(master_image, sizing_mode='stretch_width')

            link_ranges(figures_to_link)
            link_crosshairs(figures_to_link)

            # Restore a carried-over view range (from switching genes) before
            # wiring up capture, so we don't immediately re-report our own
            # restored range as a fresh "change".
            self._apply_restored_range(bokeh_expr)
            self._sync_range_updates(bokeh_expr)

            # Worth noting the tradeoff in using the Bokeh pane instead of Holoviews is that we lose
            # some reactivity, but gain access to the underlying Bokeh figure.
            # Hence, why some things like link_ranges has be implemented manually.
            self.expr_pane = pn.pane.Bokeh(bokeh_expr, sizing_mode='stretch_width')
            self.cluster_pane = pn.pane.Bokeh(bokeh_cluster, sizing_mode='stretch_width')

            needed_width = max(180, max(len(str(name)) for name in self.cluster_cmap) * 7)
            needed_height = max(300, len(self.cluster_cmap) * 22 + 50)
            ghost_legend = self._create_ghost_legend().opts(
                show_legend=True, legend_position="top_left", xaxis=None, yaxis=None, show_frame=False, toolbar=None, width=needed_width, height=needed_height
            )

            bokeh_ghost = hv.render(ghost_legend, backend='bokeh')
            self._sync_cluster_legend_toggle(bokeh_ghost, category_renderers)
            self._force_legend_redraw_on_range_change(bokeh_ghost, bokeh_expr)
            self.legend_pane = pn.Column(pn.pane.Bokeh(bokeh_ghost), width=200, height=300, scroll=True, margin=(0,0,0,0))

            self.main_row = pn.Row(self.master_pane, self.expr_pane, self.cluster_pane, self.legend_pane, sizing_mode='stretch_width')

            self.spatial_grid_container = pn.Column(self.main_row, sizing_mode='stretch_width')

            # Projections
            expr_umap = create_umap_plot(self.df, self.umap_expression_agg, color_col='raw_value', cmap="cividis_r", is_categorical=False, title=f"{self.current_gene} Expression", cbar_max=self.expression_98)
            cluster_umap = create_umap_plot(self.df, self.umap_clusters_agg, color_col='clusters', cmap=self.cluster_cmap, is_categorical=True, title="Clusters")

            umap_ghost_legend = self._create_ghost_legend().opts(show_legend=True, legend_position="top_left", xaxis=None, yaxis=None, show_frame=False, toolbar=None, width=needed_width, height=needed_height)
            umap_legend_container = pn.Column(umap_ghost_legend, width=200, height=300, scroll=True, margin=(0, 0, 0, 0))

            self.umap_row_container = pn.Row(expr_umap, cluster_umap, umap_legend_container, sizing_mode='stretch_width')

            violin_base = create_violin_plot(self.df, y_col='raw_value', group_col='clusters', cmap=self.cluster_cmap, title=f"Expression Distribution: {self.current_gene} by Cluster")
            self.violin_row_container = pn.Row(violin_base, sizing_mode='stretch_width')

            return pn.Column(
                self.state_sync, self.pre_layout, self.spatial_grid_container,
                pn.layout.Divider(margin=(20, 0)), self.umap_row_container,
                pn.layout.Divider(margin=(20, 0)), self.violin_row_container,
                sizing_mode='stretch_both', max_width=MAX_CARD_WIDTH
            )
        except Exception as e:
            traceback.print_exc()
            return pn.pane.Alert(f"Error: {e}", alert_type="danger")

    def _generate_spatial_plots(self):
        """Generates raw objects and prepares the background zoom attributes."""
        if getattr(self, 'bg_image', None) is not None:
            self.bg_image = self.bg_image.opts(default_tools=["box_zoom", "wheel_zoom", "pan", "reset"], active_tools=["wheel_zoom", "pan"], min_width=EXPANDED_PLOT_MIN_WIDTH, min_height=EXPANDED_PLOT_MIN_HEIGHT)

        # This gets filled in during the cluster spatial plot creation
        # It is a registry of category to renderer, which can then be
        # used to wire up the ghost legend to toggle visibility of the real cluster plot.
        category_renderers = {}

        # Main row
        expr_plot = create_spatial_plot(self.expression_df, y_col="y_plot", color_col='raw_value',
                                        cmap=self.expression_cmap, title=f"Expression: {self.current_gene}",
                                        cbar_max=self.expression_98, min_height=EXPANDED_PLOT_MIN_HEIGHT,
                                        min_width=EXPANDED_PLOT_MIN_WIDTH, radius=self.marker_radius)
        cluster_plot = create_spatial_plot(self.clusters_df, y_col="y_plot", color_col='clusters',
                                        cmap=self.cluster_cmap, is_categorical=True, title="Clusters",
                                        category_renderers=category_renderers, min_height=EXPANDED_PLOT_MIN_HEIGHT,
                                        min_width=EXPANDED_PLOT_MIN_WIDTH, radius=self.marker_radius)

        self.master_image = self.bg_image

        # Create the image + plot overlays.
        if hasattr(self, 'bg_image_dimmed') and self.bg_image_dimmed is not None:
            main_expr = self.bg_image_dimmed * expr_plot
            main_cluster = self.bg_image_dimmed * cluster_plot
        else:
            main_expr = expr_plot
            main_cluster = cluster_plot

        return self.bg_image, main_expr, main_cluster, category_renderers

    def _init_widgets(self):
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

        # Build the pre_layout
        self.left_pre = pn.Column(
            pn.Row(
                sizing_mode='stretch_width'
            ),
        )

        self.right_pre = pn.Column(
            pn.Row(

                self.display_name,
                self.save_button,
            ),
            self.make_default,
            align="end"
        )

        self.pre_layout = pn.Row(
            #self.left_pre,
            self.right_pre
        )

        # Emit a DOM event instead of changing a URL param
        # Pass the active states of the widgets down to the JS context (some things, like gene_symbol are already in the frontend)
        self.save_button.js_on_click(
            args={
                'name_input': self.display_name,
                'default_cb': self.make_default,
                'dataset_id': self.settings.dataset_id
            },
            code=f"""
            // Create a custom event containing all the widget values
            const evt = new CustomEvent(`save_spatial_display_${{dataset_id}}`, {{
                detail: {{
                    displayName: name_input.value,
                    makeDefault: default_cb.active,
                }}
            }});

            // Dispatch it to the host window so the vanilla JS listener catches it instantly
            window.dispatchEvent(evt);
            """)  # noqa: F541
