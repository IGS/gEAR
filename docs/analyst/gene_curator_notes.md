# Single- and Multi-gene curator notes

With the v2 rewrite of these two curators, I decided to take a shift in approach of how these are written

## Common curator code

Both the single gene curator and the multi gene curator have many shared components. These are all stored in `curator_common.js`, including the `window.onload()` function that runs when each page is opened. To make this work, many of the components have shared IDs and classes, so that the correct functions, actions, and events are called no matter the page.

Of course there are differences between pages, such as the number of the genes that can be selected for instance. To address this, many functions in `curator_common.js` have a `curatorSpecific<function>` function that is handled by each respective page. The downside of this approach is that this function must be declared on each page, even if empty.

## Plotting options

Each of the plot types available for both curators vary wildly in the relevant options available. When the user chooses a plot type, a new PlotStyle-class object will be created, which contains various properties and functions for handling that plot type.  For the single-gene curator, these are divided into PlotlyStyle, ScanpyStyle, and SvgStyle, and in the multi-gene curator we use the DashStyle class.

One of the class functions `loadPlotHtml` loads HTML with plot options for the backend plotting the request, and in some cases appending more HTML for plot-specific types. This HTML is included in `<root>/www/include/plot_config` and this is further broken down into HTML to render before the plot was made (`pre_plot`) and after the plot was made (`post_plot`)

## Error handling

There is a function in `common.js` called `logErrorInConsole` which checks various error scenarios and logs the error with `console.error`. Occasionally an error may trigger the `createToast` function from `curator_common.js` which creates a Toast-like notiification in the upper-right of the screen to inform the user.

## Element naming structures

Many of the elements responsible for the plotting options are duplicated... once for the pre-plot view and again for the post-plot view.  To keep things in sync, I assign each related element the same class under the name format `js-<backend>-<option>`.  A "change" event listener is also added to sync the values based on the event.target value. In addition, any element in the post-plot view will have "_post" added to their ID name.