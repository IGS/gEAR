from flask_restful import Resource
import geardb

class DatasetDisplay(Resource):
    """Dataset Display

    Returns
    -------
    dict

    """
    def get(self, display_id):
      display = geardb.get_display_by_id(display_id=display_id)
      return display, 200
