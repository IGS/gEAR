### App Components

epiviz-chart components are simple user interface \(UI\) elements. They cannot make data requests or can directly interact with other epiviz elements on the page. Chart elements emit hover events that propagate up the DOM hierarchy. To build interactive web applications or to coordinate events, interactions and data requests across different chart elements, we encapsulate charts inside app components.

Another essential part of the epiviz design is that data and plots are separated: you can visualize multiple charts from the same data object without having to replicate the data. This way, data queries are made by the data object, not per chart, which leads to a more responsive design of the system.

epiviz-app components are abstract components that

1. Manage layouts
2. Coordinate interactions across charts by genomic position and 
3. Handle and cache data across charts.



---

### Supported App Components

* [epiviz-environment](/app-components/epiviz-environment.md)
* [epiviz-navigation](/app-components/epiviz-navigation.md)



