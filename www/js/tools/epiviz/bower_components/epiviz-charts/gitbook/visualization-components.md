### Visualization Components

_**epiviz-chart**_ components are a collection of reusable and extensible visualization components for genomic data. We use D3.js JavaScript library to render customizable and interactive charts. An epiviz-chart component requires two attributes to render a visualization on the page

* Data attribute - a JSON \(JavaScript Object Notation\) representation of genomic data.
* Dimensions \(or columns\) from the data attribute to visualize.

_**epiviz-chart**_ components are reactive components that render visualizations only after the data attribute is attached to the element. Any change to the data attribute triggers an event to revisualize the chart. Each visualization is extensible and easily customizable to define various settings and colors.

---

### Supported Visualizations

We support both track and feature \(plots\) based visualizations.

* Tracks

  * [epiviz-genes-track](/visualization-components/epiviz-genes-track.md)
  * [epiviz-blocks-track](/visualization-components/epiviz-blocks-track.md)
  * [epiviz-stacked-blocks-track](/visualization-components/epiviz-stacked-blocks-track.md)
  * [epiviz-line-track](/visualization-components/epiviz-line-track.md)
  * [epiviz-stacked-line-track](/visualization-components/epiviz-stacked-line-track.md)

* Plots
  * [epiviz-scatter-plot](/visualization-components/epiviz-scatter-plot.md)
  * [epiviz-heatmap-plot](/visualization-components/epiviz-heatmap-plot.md)
  * [epiviz-line-plot](/visualization-components/epiviz-line-plot.md)
  * [epiviz-stacked-line-plot](/visualization-components/epiviz-stacked-line-plot.md)



