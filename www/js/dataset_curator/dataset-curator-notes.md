# Dataset curator notes

The dataset curator page is a Vue.JS app and here are a collection of notes that can aid in the development and troubleshooting of the codebase.

I wanted to make a chart of this, but LucidCharts has an entity limit :-(

## General notes

* The app uses a Vuex store to manage state for the app.
* Certain components will call /api scripts to retrieve some data.  In some cases it is mysql database data and in others it is h5ad data.
* Each plot's arguments are stored in the mysql database as JSON.  This is so the plot can be recreated outside of the Vue app, such as in the gene search display panel.
* Plotly plots are returned as JSON and drawn.  SVG and tSNE plots are returned as SVG format and base64-encoded PNG strings respectively.
* All Vue components use vue-bootstrap
* Plotly arguments uses vee-validate for validation purposes in a few areas

## Component Hierarchy

### Top-level app

* app
  * datasetCurator
    * datasetTitle
    * routes
      * /
        * datasetDisplays
          * userDisplays
            * 'edit' button routes to datasetDisplay
          * ownerDisplays
            * addDisplayBtn
      * /edit
        * datasetDisplay
      * /new
        * newDisplay

* Common components shared between newDisplay, datasetDisplay, userDisplays, ownerDisplays
  * plotlyChart
  * svgChart
  * tsneChart

### Route to create a brand new plot through arguments

* newDisplay
  * configurationPanel
    * primaryConfig
      * chooseDisplayType
      * plotlyDisplay
        * barDisplay
          * displayNameInput
          * displayColors
          * displayOrder
          * displayPalettes
          * geneSymbolInput (uses vue-bootstrap-typeahead)
          * saveDisplayBtn
          * plotlyArguments (uses vee-validate)
            * verticalLine
        * contourDisplay (extends barDisplay)
        * lineDisplay (extends barDisplay)
        * scatterDisplay (extends barDisplay)
        * tsnePlotlyDisplay(extends scatterDisplay)
        * violinDisplay (extends barDisplay)
      * svgDisplay
        * displayNameInput
        * geneSymbolInput (uses vue-bootstrap-typeahead)
        * saveDisplayBtn
      * tsneDisplay
        * displayNameInput
        * displayColors (uses vueDraggable)
        * geneSymbolInput (uses vue-bootstrap-typeahead)
        * saveDisplayBtn
        * tsneArguments
    * storedAnalysisConfig (also calls primaryConfig)
      * chooseStoredAnalysis

## Application hierarchy

* dataset_curator.html
* css/dataset_curator.css
* js/dataset_curator.js
  * Main Vue App,
  * routes
  * Vuex Store
* js/dataset_curator
  * components
    * Every individual vue component from "component hierarchy" has a unique file
