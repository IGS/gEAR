# Include our lib directory on system path
# so we have access to modules
import os
import sys
from pathlib import Path

# Prevent matplotlib backend rendering errors
# when imporing scanpy in resource modules
import matplotlib

matplotlib.use('Agg')

from flask import Flask
from flask_restful import Api

TWO_LEVELS_UP = 2
abs_path_gear = Path(__file__).resolve().parents[TWO_LEVELS_UP]
abs_path_lib = abs_path_gear.joinpath('lib')
# abs_path_lib is a Path object so we need to convert to string
sys.path.insert(0, str(abs_path_lib))

debug = os.environ.get('DEBUG', False)

from resources.aggregations import Aggregations  # noqa: E402
from resources.analyses import Analyses  # noqa: E402
from resources.available_display_types import (  # noqa: E402
    AvailableDisplayTypes,
    MGAvailableDisplayTypes,
)
from resources.dataset_display import DatasetDisplay  # noqa: E402
from resources.epiviz_data import EpivizData  # noqa: E402
from resources.gene_symbols import GeneSymbols  # noqa: E402
from resources.h5ad import H5ad  # noqa: E402
from resources.mg_plotly_data import MGPlotlyData  # noqa: E402
from resources.orthologs import Orthologs  # noqa: E402

# Import resources
from resources.plotly_data import PlotlyData  # noqa: E402
from resources.projectr import ProjectR, ProjectROutputFile  # noqa: E402
from resources.spatial_scanpy_data import SpatialScanpyData  # noqa: E402
from resources.svg_data import SvgData  # noqa: E402
from resources.top_pca_genes import TopPCAGenes  # noqa: E402
from resources.tsne_data import MGTSNEData, TSNEData  # noqa: E402

app = Flask(__name__)
api = Api(app)
# Add API endpoints to resources

# Set the maximum memory usage for the process
# to 95% of the system's maximum memory usage
import resource

soft, hard = resource.getrlimit(resource.RLIMIT_DATA)
soft95 = int(soft * 0.95)
resource.setrlimit(resource.RLIMIT_DATA, (soft95, hard))

api.add_resource(PlotlyData, '/plot/<dataset_id>'   # Default endpoint
                 , "/plot/<dataset_id>/plotly")     # add /plotly to this endpoint for name consistency with other endpoints
api.add_resource(MGPlotlyData,'/plot/<dataset_id>/mg_plotly')
api.add_resource(SvgData, '/plot/<dataset_id>/svg')
api.add_resource(TSNEData, '/plot/<dataset_id>/tsne')
api.add_resource(MGTSNEData, '/plot/<dataset_id>/mg_tsne')
api.add_resource(EpivizData, '/plot/<dataset_id>/epiviz')
api.add_resource(SpatialScanpyData, '/plot/<dataset_id>/spatial_scanpy')
api.add_resource(ProjectR, '/projectr/<dataset_id>')
api.add_resource(ProjectROutputFile, '/projectr/<dataset_id>/output_file')
api.add_resource(H5ad, '/h5ad/<dataset_id>')
api.add_resource(AvailableDisplayTypes, '/h5ad/<dataset_id>/availableDisplayTypes')
api.add_resource(MGAvailableDisplayTypes, '/h5ad/<dataset_id>/mg_availableDisplayTypes')
api.add_resource(Aggregations, '/h5ad/<dataset_id>/aggregations')
api.add_resource(Analyses, '/h5ad/<dataset_id>/analyses')
api.add_resource(Orthologs, '/h5ad/<dataset_id>/orthologs')
api.add_resource(GeneSymbols, '/h5ad/<dataset_id>/genes')
api.add_resource(TopPCAGenes, '/analysis/plotTopGenesPCA')
api.add_resource(DatasetDisplay, '/displays/<display_id>')

if __name__ == '__main__':
    # Problem:
    # Our client API makes requests to /api/... in order
    # to get data from our Flask Application. This prefix is
    # configured in Apache as an endpoint to communicate with Flask.
    # This is necessary in the production environment. However, in
    # development, we may not want Apache serving our Flask application
    # because it cannot detect changes and we will have to restart
    # Apache everytime.

    # Instead, we want to run Flask's development server, but this server
    # does not have the /api/ prefix. We would have to change all of our resource
    # routes to include the prefix, but then this same code would not work
    # in production.

    # Solution:
    # Add this middleware to our application that will be run if we are
    # running this code in development, e.g, starting the flask server
    # ourselves with python api.py

    # https://stackoverflow.com/a/36033627
    # class PrefixMiddleware(object):
    #     def __init__(self, app, prefix=''):
    #         self.app = app
    #         self.prefix = prefix

    #     def __call__(self, environ, start_response):

    #         if environ['PATH_INFO'].startswith(self.prefix):
    #             environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
    #             environ['SCRIPT_NAME'] = self.prefix
    #             return self.app(environ, start_response)
    #         else:
    #             start_response('404', [('Content-Type', 'text/plain')])
    #             return ["This url does not belong to the app.".encode()]
    # api.add_resource(PlotlyData, '/api/plot/<dataset_id>')
    # api.add_resource(H5ad, '/api/h5ad/<dataset_id>')
    # app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix='/api')
    app.run(debug=debug, threaded=True)

