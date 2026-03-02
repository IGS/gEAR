import logging
import sys

import panel as pn
from panel_common import CondensedSpatialViewer

# Reset logging level to "error" to suppress a bokeh "dropping patch" info message
# https://github.com/bokeh/bokeh/issues/13229
logging.getLogger().setLevel(logging.ERROR)

# If not params passed, just show OK as a way to test the app
if pn.state.location is None or not pn.state.location.query_params:
    pn.pane.Markdown("OK").servable()
else:
    # Instantiate the component
    spatial_dashboard = CondensedSpatialViewer()
    spatial_dashboard.servable()