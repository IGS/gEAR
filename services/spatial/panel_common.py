import traceback
import datashader as ds
import holoviews as hv
import hvplot
import hvplot.pandas  # noqa
import numpy as np
import panel as pn
from common import (
    create_spatial_plot, create_umap_plot, create_violin_plot,
    retrieve_dataframe, retrieve_image_array, normalize_expression_name, has_selection, ExpandedSettings
)

# CRITICAL: Initialize the Bokeh backend for interactivity
hvplot.extension('bokeh', logo=False) # type: ignore
pn.extension(loading_indicator=True, defer_load=True, nthreads=4)  # type: ignore)

class BaseSpatialViewer(pn.viewable.Viewer):
    """
    Base Viewer component. Handles state and linking.
    """

    settings = ExpandedSettings()

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
        self.linker = hv.link_selections.instance(unselected_alpha=1)
        self.bounds_stream = hv.streams.BoundsXY(bounds=self.saved_bounds)  # type: ignore

        # Set up a callback to update the URL params whenever the user draws or clears a box
        self.bounds_stream.param.watch(self._sync_stream_to_params, 'bounds')

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
                        xaxis=None, yaxis=None, frame_width=300, frame_height=200
                    )

        self.layout_height = 312  # 360px - tile header height
        if self.settings.display_height and self.settings.display_height > 0: # type: ignore
            self.layout_height: int = self.settings.display_height # type: ignore

        self.layout_width = 1100  # Default width
        if self.settings.display_width and self.settings.display_width > 0: # type: ignore
            self.layout_width: int = self.settings.display_width # type: ignore

        self.layout = pn.Column(pn.bind(self._build_layout), width=self.layout_width)

    def _sync_stream_to_params(self, event):
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

    def _loading_indicator(self, label):
        return pn.indicators.LoadingSpinner(
            value=True, name=label, align="center", color="info"
        )

    def _build_layout(self):
        """
        This is where the main Panel layout is built.
        We leave it blank in the base class since the Condensed and Expanded viewers will have different layouts.
        """
        raise NotImplementedError("Subclasses must implement _build_layout")

    def __panel__(self):
        """
        Panel automatically looks for this method.
        It MUST return a Panel viewable object (Row, Column, Pane, etc.)
        """
        return self.layout

class CondensedSpatialViewer(BaseSpatialViewer):
    """
    The specific component for your panel_app.py
    """

    def _build_layout(self):
        """Builds the 3-panel condensed row."""

        try:
            yield self._loading_indicator("Creating plots. Please wait...")

            # One unfortunately annoyance is that datashader's default behavior is to flip the y-axis,
            # which is not what we want for spatial data. To fix this,
            # we can reverse the y-axis limits by setting ylim to (max, min) instead of (min, max).
            self.df.loc[:, "y_plot"] = self.df["spatial2"]
            if self.img_height is not None:
                self.df.loc[:, "y_plot"] = self.img_height - self.df["spatial2"]

            # 1. Generate base plots
            image_panel = self.bg_image if hasattr(self, 'bg_image') else None
            plot_expr = create_spatial_plot(self.df, self.expression_agg, y_col='y_plot', color_col='raw_value', cmap=self.expression_cmap) # type: ignore
            plot_clust = create_spatial_plot(self.df, self.clusters_agg, y_col='y_plot', color_col='clusters', cmap=self.cluster_cmap) # type: ignore
            if image_panel:
                plot_expr = image_panel * plot_expr # type: ignore
                plot_clust = image_panel * plot_clust # type: ignore

            # If user has a saved box, draw it on the main plots so they see it
            if self.saved_bounds:
                saved_box = hv.Bounds(self.saved_bounds).opts(color='black', line_width=2)
                if image_panel is not None:
                    image_panel = image_panel * saved_box   # type: ignore
                plot_expr = plot_expr * saved_box   # type: ignore
                plot_clust = plot_clust * saved_box   # type: ignore

            # 2. Attach the stream to capture drawn boxes
            self.bounds_stream.source = plot_expr

            # 3. Apply the linker for cross-filtering
            linked_expr = self.linker(plot_expr)
            linked_clust = self.linker(plot_clust)

            # 4. Lay out the non-zoom panels side-by-side using HoloView
            main_row = (image_panel + linked_expr).opts(shared_axes=True)

            # 5. Define the dynamic Zoom Panel
            def zoomed_panel(bounds):
                zoom_expression_cmap = "YlGn"

                zoom_base = (create_spatial_plot(self.df, self.expression_agg, color_col='raw_value', cmap=zoom_expression_cmap)
                            .opts(title="Draw box in one of the other plots to zoom this one.")
                )
                if bounds is None:
                    return zoom_base.opts(title="Draw box to zoom")
                l, b, r, t = bounds
                return zoom_base.opts(xlim=(l, r), ylim=(b, t), title="Zoomed View")

            # 6. Wrap zoom panel in a DynamicMap to auto-update on box draw
            #plot_zoom = hv.DynamicMap(zoomed_panel, streams=[self.bounds_stream])
            plot_zoom = None

            markdown_width = 400    # Manually measured.
            markdown_padding = 20  # Total left/right padding for the markdown pane
            self.intro_markdown = pn.pane.Markdown(
                "### Click the Expand icon in the top right corner to see all plots", width=markdown_width
            )

            self.use_clusters = False
            self.use_clusters_switch = pn.widgets.Switch(
                value=self.use_clusters, sizing_mode="fixed", margin=(10, 10)
            )

            # Using HTML to center the labels (https://github.com/holoviz/panel/issues/1313#issuecomment-1582731241)
            switch_content_width = 250
            self.switch_layout = pn.Row(
                pn.pane.HTML("""<label><strong>Gene Expression</strong></label>""", sizing_mode="fixed", margin=(10, 10)),
                self.use_clusters_switch,
                pn.pane.HTML("""<label><strong>Clusters</strong></label>""", sizing_mode="fixed", margin=(10, 10)),
                width=switch_content_width
            )

            # When width is too short, things break down.
            if self.layout_width < 1100:
                self.layout_width = 1100
                self.layout_height = 312

            spacer_width = self.layout_width - markdown_width - markdown_padding - switch_content_width
            if spacer_width < 0:
                spacer_width = 100

            self.pre_layout = pn.Row(
                self.intro_markdown,
                pn.Spacer(width=spacer_width),
                self.switch_layout,
                height=30
            )

            # Return final Panel layout
            yield pn.Column(self.pre_layout, pn.pane.HoloViews(main_row), pn.pane.HoloViews(plot_zoom))
        except Exception as e:
            traceback.format_exc()
            yield pn.pane.Alert(f"Error: {e}", alert_type="danger")

class ExpandedSpatialViewer(BaseSpatialViewer):
    """
    The specific component for your panel_app_expanded.py
    """

    def _build_layout(self):
        # Build your spatial rows, UMAPs, and Violins here...

        try:
            yield self._loading_indicator("Creating plots. Please wait...")

            main_plots = None
            zoom_plots = None
            umap_row = None
            violin_row = None

            yield pn.Column(main_plots, zoom_plots, umap_row, violin_row)
        except Exception as e:
            traceback.format_exc()
            yield pn.pane.Alert(f"Error: {e}", alert_type="danger")