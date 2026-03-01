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

    # Drop it directly into the page template
    template = pn.template.GoldenTemplate(
        title='gEAR Spatial View',
        main=[spatial_dashboard] # Treats the class instance exactly like a pn.Row!
    )
    template.servable()