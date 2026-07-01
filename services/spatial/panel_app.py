import logging

import panel as pn
from panel_common import CondensedSpatialViewer

# Reset logging level to "error" to suppress a bokeh "dropping patch" info message
# https://github.com/bokeh/bokeh/issues/13229
logging.getLogger().setLevel(logging.ERROR)

args = pn.state.session_args

# If no params passed, just show OK as a way to test the app
if not args:
    pn.pane.Markdown("OK").servable()
else:
    # Instantiate the component
    spatial_dashboard = CondensedSpatialViewer()
    spatial_dashboard.servable()