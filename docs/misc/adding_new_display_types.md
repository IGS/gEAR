# Purpose

The gene expression pages of gEAR instances support multiple types of dataset displays, implemented using several different technologies.  These include static images (Scanpy), colorized SVGs via D3.js/Snap.svg, interactive plots via Plotly.js, Gosling tracks, and the spatial Panel viewer.  This document describes the process of adding a NEW type to the system, and could eventually become the groundwork for end-user supplied display plug-ins.

## Process stack

1. `expression.html` (via `www/js/expression.js`) and `projection.html` (via `www/js/projection.js`) create a `TileGrid` (`www/js/classes/tilegrid.js`) for the selected dataset collection (layout) `shareId`.
2. `TileGrid.getLayout()` / `getDatasets()` load the collection members and dataset info through `apiCallsMixin` in `www/js/common.v2.js`; `addAllDisplays()` / `addDefaultDisplay()` fetch each dataset's saved displays and default display.
3. When genes are searched, `TileGrid.renderDisplays()` calls `DatasetTile.renderDisplay()` for each tile.
4. `DatasetTile.renderDisplay()` dispatches on `display.plot_type` to a renderer (`renderPlotlyDisplay`, `renderScanpyDisplay`, `renderSVG`, `renderGoslingDisplay`, `renderMultiGeneDisplay`, `renderSpatialPanelDisplay`). The expanded (zoomed) view is handled by `TileGrid.applySingleTileGrid(..., isZoomed=true)`, which builds a second `DatasetTile` and calls the same `renderDisplay()`.

Displays are created in the curators: `www/js/dataset_curator.js` (single gene) and `www/js/multigene_curator.js`, which share `www/js/curator_common.js`.

## Addition notes

Throughout the documentation, replace $TOOL with whatever the tool type name actually is.

### www/api/resources/$TOOL_data.py

In this directory, you need to create a module for the data representing this tool, showing how to access the data it requires for any given gene, and register it in `www/api/api.py` with `api.add_resource()` (existing examples: `plotly_data.py`, `tsne_data.py`, `svg_data.py`, `spatialpanel.py`).

### www/api/resources/available_display_types.py

Here you need to add the logic necessary to determine if the new display type is present for any given dataset.  `AvailableDisplayTypes` serves `/api/h5ad/<dataset_id>/availableDisplayTypes` (single gene) and `MGAvailableDisplayTypes` serves `/api/h5ad/<dataset_id>/mg_availableDisplayTypes` (multigene).  The curators call these via `fetchAvailablePlotTypes()`.

### www/js/classes/tilegrid.js

- Add the new `plot_type` to the dispatch in `DatasetTile.renderDisplay()` (or to the `plotlyPlots` / `scanpyPlots` lists at the top of the file if it reuses an existing renderer).
- Add a `render$TOOLDisplay()` method on `DatasetTile`, plus download handling if applicable.
- The zoomed view uses the same method, so check that it behaves with `this.isZoomed`.

### www/js/dataset_curator.js / multigene_curator.js / curator_common.js

Add the type to the curator's plot-type lists and implement its plot style (the curators register a `curatorSpecificPlotStyle()` function with `curator_common.js`), including pre/post plot option HTML in `www/include/plot_config/pre_plot/` and `post_plot/`.
