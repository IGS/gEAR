"""
dataset_display.py - Fetch a saved dataset display configuration.

Serves /displays/<display_id> in www/api/api.py.
"""

from flask_restful import Resource
import geardb

class DatasetDisplay(Resource):
    """Dataset Display

    Returns
    -------
    dict

    """
    def get(self, display_id):
      """
      Return the display with the given ID, or a 404 message if not found.
      """
      display = geardb.get_display_by_id(display_id=display_id)

      if not display:
          return {"message": "Display not found"}, 404

      return display, 200
