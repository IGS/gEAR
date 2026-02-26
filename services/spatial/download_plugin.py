from io import BytesIO

import panel as pn
from panel_app_expanded import ExpandedSettings, SpatialPanel
from tornado.web import HTTPError, RequestHandler

pn.extension("plotly", loading_indicator=True, defer_load=True, nthreads=4)

class DownloadHandler(RequestHandler):
    """
    throw‑away app that only exists to render the current spatial panel as a
    PNG.  It is mounted at a different prefix so every fetch creates a new
    session and `__panel__` runs.

    We cannot use the panel_app_expanded as the download endpoint because the download
    is triggered by a client-side fetch, which does not have access to the same session
    state as the main app (opened by a websocket).  By creating a separate app that reads settings from query
    parameters, we can ensure that the download endpoint has all the information it
    needs to render the panel correctly.
    """

    arg_mapping = {
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
    }

    def get(self):
        # self.request.arguments values are a list of a single-byte element.
        # this also decodes the arguments.
        query_args = {}
        for arg in self.request.arguments:
            # Ensure mappings work with the names used in "settings"
            mapped_key = self.arg_mapping.get(arg)
            if not mapped_key:
                continue
            query_args[mapped_key] = self.get_query_argument(arg)
            # if argument can be coaxed into an int or a float, do so
            # As the "settings" have specific types that are checked.
            try:
                query_args[mapped_key] = int(query_args[mapped_key])
            except ValueError:
                try:
                    query_args[mapped_key] = float(query_args[mapped_key])
                except ValueError:
                    pass

        # build settings from whatever query parameters were supplied on
        # this request; note that session_args is *not* used here.
        settings = ExpandedSettings(**query_args)

        # construct a transient SpatialPanel and force it to load its data
        try:
            panel = SpatialPanel(settings)
            # run the init_data generator to completion so df, maps, etc. exist
            for _ in panel.init_data():
                pass
        except Exception as e:
            raise HTTPError(400, f"Failed to create panel: {e}")

        # choose which parts of the layout we want in the image
        layout = pn.Column(
            pn.pane.Markdown("## Full view", height=30),
            panel.normal_pane,
            panel.fig_layout
            )

        buf = BytesIO()
        pn.io.save.save(layout, buf, resources="cdn", embed=True)   # Will eventually convert to PDF in the browser so do not need standalone HTML
        buf.seek(0)
        self.set_header("Content-Type", "text/html")
        self.write(buf.read())

# make it servable on a distinct prefix
ROUTES = [('/spatial_download', DownloadHandler, {})]