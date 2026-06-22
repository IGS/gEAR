import sys
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

    settings = Settings()   # Set default settings

    # Any potential bound widgets should be declared here so they are accessible in the base class and subclasses
    use_clusters = param.Boolean(default=False, doc="Whether to show clusters or gene expression in the main plots")
    min_genes = param.Integer(default=0, doc="Minimum number of genes expressed to include a cell observation", bounds=(0, 500))

    def __init__(self, **params):
        super().__init__(**params)
        """
        DataFrame columns
        raw_value,spatial1,spatial2,n_genes_by_counts,UMAP1,UMAP2,clusters,clusters_cat_codes,colors
        """

        if pn.state.session_args is not None:

            args = pn.state.session_args

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

            if 'min_genes' in args:
                self.settings.min_genes = int(get_arg('min_genes', 0))

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
                self.settings.nosave = bool(int(get_arg('nosave', True)))

            if 'display_name' in args:
                self.settings.display_name = get_arg('display_name', "")

            if 'make_default' in args:
                self.settings.make_default = bool(int(get_arg('make_default', False)))

        self.orig_df = retrieve_dataframe(self.settings.dataset_id, self.settings.filename)
        self.orig_df['clusters'] = self.orig_df['clusters'].astype('category')
        self.orig_df = clip_expression_values(self.orig_df, self.settings.expression_min_clip)

        # If min_genes is set, filter the dataframe to only include observations with at least that many genes
        self.min_genes = self.settings.min_genes
        self._filter_df()

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

        # Initialize linking (data filtering) and streams (coordinate reporting)
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
                        tools=["box_select"], active_tools=["box_select"], default_tools=[],
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
                        tools=["box_select"], active_tools=["box_select"], default_tools=[],
                        hooks=[autohide_toolbar]
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
            self.saved_bounds = None
        else:
            # User drew a box
            (self.settings.selection_x1, self.settings.selection_y1, self.settings.selection_x2, self.settings.selection_y2) = bounds
            self.saved_bounds = bounds

        # Trigger the zoom update!
        self._update_zoom_panel(bounds)

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
                # TODO: Figure out how to get self.linker and the bounds stream to work together

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
                "### Click the Expand icon in the top right corner for added functionality",
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

        # TODO: Figure out how to get self.linker and the bounds stream to work together
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

        # TODO: Trigger zoom update on initial load or if switch is toggled.
        #self._update_zoom_panel(self.saved_bounds)

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

            # One unfortunately annoyance is that datashader's default behavior is to flip the y-axis,
            # which is not what we want for spatial data. To fix this,
            # we can reverse the y-axis limits by setting ylim to (max, min) instead of (min, max).
            self.df["y_plot"] = self.df["spatial2"]
            if self.img_height is not None:
                self.df["y_plot"] = self.img_height - self.df["spatial2"]

            # Represents main row and zoom row
            spatial_grid = self._generate_spatial_grid()

            ### UMAP row
            expr_umap = create_umap_plot(
                self.df, self.expression_agg, color_col='raw_value', cmap="cividis_r", is_categorical=False, title=self.current_gene
            )
            cluster_umap = create_umap_plot(
                self.df, self.clusters_agg, color_col='clusters', cmap=self.cluster_cmap, is_categorical=True, title="Clusters"
            )

            # Wrap them in the cross-filtering linker
            linked_expr_umap = self.linker(expr_umap)
            linked_cluster_umap = self.linker(cluster_umap)

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
            umap_row = pn.Row(
                linked_expr_umap,
                linked_cluster_umap,
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

            linked_violin = self.linker(violin_base)

            violin_row = pn.Row(
                linked_violin,
                sizing_mode='stretch_width'
            )

            return pn.Column(
                self.pre_layout,
                spatial_grid,
                pn.layout.Divider(margin=(20, 0)), # Visual breathing room
                umap_row,
                pn.layout.Divider(margin=(20, 0)),
                violin_row,
                sizing_mode='stretch_both'
            )
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
                    "## Select a region to modify zoomed in view in the second row.",
                    height=30,
                    width=markdown_width,
                ),
                self.display_name,
                self.save_button,
            ),
            pn.Row(self.min_genes_slider, pn.Spacer(width=spacer_width), self.make_default),
        )

        # Emit a DOM event instead of changing a URL param
        # Pass the active states of the widgets down to the JS context (some things, like gene_symbol are already in the frontend)
        self.save_button.js_on_click(
            args={
                'name_input': self.display_name,
                'default_cb': self.make_default,
                'min_genes_slider': self.min_genes_slider,
                'dataset_id': self.settings.dataset_id
            },
            code=f"""
            // Create a custom event containing all the widget values
            const evt = new CustomEvent(`save_spatial_display_${{dataset_id}}`, {{
                detail: {{
                    displayName: name_input.value,
                    makeDefault: default_cb.active,
                    minGenes: min_genes_slider.value,
                    // If you tracked bounds in JS, you could grab them, or grab them from Panel
                }}
            }});

            // Dispatch it to the host window so the vanilla JS listener catches it instantly
            window.dispatchEvent(evt);
            """)

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
        # Generate base plots
        expr_plot = create_spatial_plot(self.df, self.expression_agg, y_col="y_plot", color_col='raw_value', cmap=self.expression_cmap, title=f"Expression: {self.current_gene}")
        cluster_plot = create_spatial_plot(self.df, self.clusters_agg, y_col="y_plot", color_col='clusters', cmap=self.cluster_cmap, is_categorical=True, title="Clusters") # type: ignore
        master_image = self.bg_image if hasattr(self, 'bg_image') else pn.pane.HoloViews(None)

        # Update bounds streams for the image and expression plots.
        self.bounds_stream_image.source = master_image
        self.bounds_stream_composite.source = expr_plot

        # Create a 3rd stream for the cluster plot so it can also trigger the zoom
        if not hasattr(self, 'bounds_stream_cluster'):
            self.bounds_stream_cluster = hv.streams.BoundsXY(bounds=self.saved_bounds)
            self.bounds_stream_cluster.add_subscriber(self._update_bounds_callback)
            self.bounds_stream_cluster.source = cluster_plot

        # Apply cross-filtering linker BEFORE adding images
        linked_expr = self.linker(expr_plot)
        linked_cluster = self.linker(cluster_plot)

        # Create composites with the background image
        if hasattr(self, 'bg_image_dimmed') and self.bg_image_dimmed is not None:
            # Row 1 (Linked)
            main_expr = self.bg_image_dimmed * linked_expr
            main_cluster = self.bg_image_dimmed * linked_cluster
            # Row 2 Zoom (Unlinked/Raw)
            zoom_master_expr = self.bg_image_dimmed * expr_plot
            zoom_master_cluster = self.bg_image_dimmed * cluster_plot
        else:
            main_expr = linked_expr
            main_cluster = linked_cluster
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
