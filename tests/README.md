# Testing

Manual UI testing is time consuming and is one of the fastest ways to make people look for another job. To keep our great people, we've created this automated testing framework to exercise as many features of the site as possible. This should free us to do more development and identify issues faster.

Generally testing progresses in a few phases.

## UI Testing

### How to set up packages

Ensure npm is installed, and make sure you are in the `<gear_root>/tests` directory

```bash
npm install --save-dev mocha chai playwright
```

In the created package.json file, add within the outermost braces (you can append CLI options, like "--timeout" to this as well):

```json
"scripts": {
  "test": "mocha"
}
```

### About the UI testing packages

Mocha (https://mochajs.org/) is a widely-used testing framework that has a minimal setup, allows for organization of tests, and is flexible to let you choose your own libraries to extend its functionality. Mocha has interfaces for both Behavior-driven development (BDD) and Test-driven development (TDD).

Playwright (https://playwright.dev/) is a framework used for automated end-to-end testing. You can use it to interact and test on page selectors. It also provides headless browser testing by default and works with multiple browser types, including mobile ones.

### Current and pending UI tests

#### Account creation

* Navigate from home page to account creation
* Test error handling
  * Leaving out required field (e-mail)
  * Submitting with duplicate account info
* Saving of new account info

#### Main Page (index)

* Logging in
  * With incorrect credentials
  * With correct credentials

#### Comparison tool

* Select dataset
* Select conditions
* Ensure condition labels are reflected in the plot
* Ensure plot can be generated
  * Default options
  * Significance test
    * Filter
    * Color
* Gene highlighting
  * Found genes show in plot
  * Not found genes show in that div
* Select genes from plot
  * Ensure they show in table
  * Ensure highlighted genes are colored (Pou4f3)
* Visual regression testing
* Download gene selection table
* Name and save gene cart

#### Multigene Curator

* Create a heatmap
  * Must have 2+ genes
  * Alt heatmap with cluster observations checkbox
  * Alt heatmap with cluster genes checkbox
  * Alt heatmap with axes flipped
  * Distance metric for clustering observations/genes
  * Matrix plot
  * Sort by primary category
* Create a violin
  * Stacked violin plot
  * With jitter
* Create a volcano
  * REQUIRED - query/ref conditions
    * Window alert if not chosen or category is different
  * DE Algorithm
  * Annotate non-signficant p-values
  * Use adjusted p-vals
* Create a dotplot
* Create a quadrant plot
  * REQUIRED - query1/query2/ref conditions
    * Window alert if not chosen or category is different
  * DE Algorithm
  * Foldchange cutuff
  * FDR cutoff

#### Main Page - Display Panel mode

#### Contact Us

#### Analysis (single-cell) workbench

#### Colorblind mode

* Comparision tool
* Dataset curator
* Multigene curator
* Main page - display panel

#### Tests that use a private dataset instead of a public one (requires login)

#### Dataset (single-gene) curator

* Create a bar plot
* Create a scatter plot
  * Facet rows and cols
  * color
  * marker size
* Create a line plot
* Create a violin plot
* Create a tSNE/uMAP plot
  * split and color by a category
* Load an existing plot

#### Multigene curator

* Ensure plot is loaded when dataset is chosen (loaded from saved displays or default volcano plot)
* Volcanoes should be disabled when no categories have 2+ groups
* Quadrants should be disabled when no categories have 3+ groups
* Misc.
  * Primary category
  * Secondary category (may need to choose a new dataset)
  * Sort category (either primary or secondary)
* Load an existing plot
* Save a plot
* Save a new gene cart
  * From volcano or quadrant
* Use visual tool like needle (python package) for visual regression testing (plot doesn't differ)
* Ensure heatmap and matrixplot expression value for a single gene and observation/celltype is correct
  * This tests that sorting was fine
  * For heatmaps/violins/dotplots

#### Analysis (single-cell) Workbench

* Go through all steps in new analysis
* Resume unsaved analysis
* Resume saved analysis
* Load primary analysis

#### Manual Documentation

#### Dataset Uploader - Expression Data

#### Dataset Uploader - Epigenetic Data

#### Dataset Explorer

#### Gene Collection Manager

* New Gene Collection
  * Validation check on new label
  * Validation check on save (label and organism)
  * Cancel redirects correctly
* Searching
  * Clear search term
  * Search and return results
  * Sort By change
  * Facet filters - Just one
  * Facet filters - Select two
  * Facet filters - Select All
  * Facet filter - not logged in only show Public
* Existing Gene Collection
  * Expand + Collapse
  * Expand/Collapse All
  * Edit - Save change
    * action links visible
    * correct organism label
  * Edit - Cancel
  * Delete


#### Epiviz Panel Designer

## API Testing

TODO
