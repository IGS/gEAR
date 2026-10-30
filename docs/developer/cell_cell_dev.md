# Cell-cell communication

This page adds functionality to visualize the output of a [dominoSignal](https://github.com/FertigLab/dominoSignal) analysis run on a chosen data.

The page uses the underlying Python version of Plotly (Version 6.6.0 at the time of release) to render a series of dot and river plots. Plotting details can be found in `lib/dominosignal.py`.

The code for this page can be found in `cell_cell.js`, `cell_cell.css` and `cell_cell.html`, all in the `www/` folder. Tests for this page can be found in `tests/python/test_dominosignal.py` and in `tests/legacy/test/cell_cell.test.js`. The CGI scripts that link the underlying Python functions to this front end are `www/cgi/get_domino_options.cgi` and `www/cgi/get_domino_plots.cgi`.

Datasets that are loaded for this page are stored in `www/cache/domino_signal`, and are paired `.parquet` and `.JSON files`. The former contains the dominoSignal output (which traditionally is in a `.csv` file, but this file format is significantly smaller), and the latter contains metadata about the dataset. In future iterations, we plan for this page to accept all datasets already uploaded to gEAR that also contain dominoSignal output.

This page was added by @arun96, with significant help from Claude.