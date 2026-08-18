from io import BytesIO

import panel as pn
from panel_app_expanded import ExpandedSpatialViewer
from tornado.web import HTTPError, RequestHandler

pn.extension(loading_indicator=True, defer_load=True, nthreads=4)

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

    def get(self):
        try:
            # 1. Instantiate the refactored Viewer
            panel_app = ExpandedSpatialViewer(session_args_override=self.request.arguments)

            # 2. Clean up the UI for export
            # The exported HTML is offline. Datashader images will embed as static rasters,
            # and Python sliders won't work. We strip the controls (pre_layout)
            # so the user gets a clean, professional report.
            layout_components = [
                pn.pane.Markdown(f"## Spatial Display Export: {panel_app.current_gene}", height=30)
            ]

            # Extract all visualization rows, skipping the interactive controls at index 0
            built_ui = panel_app._build_layout()
            for item in built_ui.objects:
                if item is not panel_app.pre_layout:
                    layout_components.append(item)

            export_layout = pn.Column(*layout_components, sizing_mode='stretch_both')

            # 3. Generate the HTML file
            buf = BytesIO()
            # embed=True packs the current JS and Datashader rasters into a single file
            pn.io.save.save(export_layout, buf, resources="cdn", embed=True, title="gEAR Spatial Export")
            buf.seek(0)

            # 4. Send to browser
            self.set_header("Content-Type", "text/html")
            self.write(buf.read())

        except Exception as e:
            import traceback
            traceback.print_exc()
            raise HTTPError(400, f"Failed to create HTML export: {e}")

# make it servable on a distinct prefix
ROUTES = [('/spatial_download', DownloadHandler, {})]