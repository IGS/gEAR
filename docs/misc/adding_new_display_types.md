# Purpose

The front page of gEAR instances support multiple types of dataset displays, implemented using several different technologies.  These include static images, colorized SVGs via D3.js, interactive plots via Plotly.js, etc.  This document describes the process of adding a NEW type to the system, and could eventually become the groundwork for end-user supplied display plug-ins.

## Process stack

TODO: update

1. When a profile is selected/changed dataset_collection_panel.set_layout() is called
2. set_layout() stores the selection, then calls dataset_collection_panel.load_frames()
3. load_frames() calls cgi/get_dataset_list.cgi

## Addition notes

Throughout the documentation, replace $TOOL with whatever the tool type name actually is.

### www/api/resources/$TOOL_data.py

In this directory, you need to create a module for the data representing this tool, showing how to access the data it requires for any given gene.

### www/api/resources/available_display_types.py

Here you need to add the logic necessary to determine if the new display type is present for any given dataset.

### www/js/index.js:zoom_on_dataset()

Add any code necessary to support the full-page zoomed view of this dataset/type.

### www/js/classes/display.js

Add a class for this which extends Display
