# Purpose

Personal notes for debugging how a specific analyses is shown on the front page in one of the data panels.

# The stack

In the database is a table dataset_display which lists all the configured displays for a dataset, and the user's preference is stored in the dataset_preference table (probably should have been called dataset_display_preference).

----------------
-= Python API =-

lib/geardb.py
    - Has two class methods which return the full row for dataset_displays:
      - get_default_display()
      - get_display_by_id()
      - get_displays_by_user_id()

www/cgi/get_dataset_display.cgi
    - calls geardb.get_display_by_id() and returns JSON of entire match

www/cgi/get_dataset_displays.cgi
    - calls geardb.get_displays_by_id() and returns JSON of entire match

---------------------
-= Javascript side =-

The entire panel of datasets is managed with ui-panel-dataset-collection.js (class DatasetCollectionPanel), which contains several DatasetPanel instances (from ui-panel-dataset.js).

DatasetCollectionPanel:
    - has load_frames() call which iterates through the datasets, instantiating several
      DatasetPanels
    - has update_by_search_result() call which iterates datasets, calling dataset.draw() on
      each with the current gene symbol
    
DatasetPanel extends Dataset:
    - has calls get get all dataset_displays as well as the current, default one.
    - has the draw() and draw_chart() methods which instantiate specific display
      types like PlotlyDisplay.
    - controls UI elements like loading spinner, no match display, etc.

Dataset:
    - Really just has attributes which correspond to the dataset table rows and shape() method

------------------------
-= Stack problem area =-

plotly_data.py calls:
   - geardb.Analysis.discover_type()
   - geardb.Analysis.dataset_path()
   - BOTH NEED type from dataset_display JSON!
   - called from api.add_resource(PlotlyData, '/plot/<dataset_id>')

js/classes/display.js
   - PlotlyDisplay extends Display and has API call to /api/plot/$dataset_id while
     passing this.analysis which doesn't have JSON?
