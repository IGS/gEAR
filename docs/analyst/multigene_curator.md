# Multigene Curator Tool

> The multigene curator allows users to create custom displays for a dataset they own or other public datasets of interest. A user can create one or more displays for others to toggle or set as their default display when searching genes.

## Differences between the dataset curator and the multigene curator

There are two main differences between the two curator tools:

* Both tools were developed using different frameworks and different layouts. Most of the functionality remains the same though.
* The dataset curator can only display single gene plots, and the multigene curator, as the name would suggest, displays plots relevant to multigene visualization.

## Steps to create a curated plot

1. Choose a dataset. You can choose amongst datasets you have uploaded, datasets shared with you, or public datasets
2. Configure plot type. You can either:
  a. *Choose a plot type from the plot type dropdown.*
  b. *Load an existing saved display.* To do this, click "Load a Saved Display" and click "Edit" under the display you wish to load. This will pre-load the genes used to create that plot, as well as draw the plot (plot options are not currently preloaded). There is also an option to choose which display is your default when viewing profiles from the main page.
3. Choose genes to add. You can either:
  a. *Manually select genes.* The text box will filter possible genes as the name is typed. Either clicking the gene from the dropdown list or hitting Tab adds that gene to the selection.
  b. *Load from a gene cart.* You have the option to choose a public gene cart or your own saved gene cart. When a gene cart is chosen, all of those genes will be added to the gene selection box.
4. Configure plot options. Different plot types have different options available to them. Each option has tooltips to assist with the particular option
  a. In the case of the "Filter by Category" option, only the chosen groups in each category in each group will be included when creating the plots. Clicking "All" will use all groups in the category.
  b. Each plot type also has advanced options.
5. Click "Update plot" to generate the plot.

## Saving new displays

After a plot has been generated, you can assign the plot a name and click "Save new \<plot type\>" to save the display. Note that if you previously loaded a display, saving will create a new display rather than overwrite the existing one.

## Downloading images

After a plot has been generated, you can specify the dimensions and scale of the image to download (as a PNG).

## Volcano plot things

After a volcano plot has been generated, you can use the lasso tool or the rectangular selection tool to select some gene datapoints within the plot.  These selected genes will show in a new panel, where you can assign a name to this gene list and save as a new gene cart.