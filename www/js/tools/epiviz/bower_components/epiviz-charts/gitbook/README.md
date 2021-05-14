## **Introduction**

Epiviz Web Components are a collection of open source reusable and extensible data visualization and app elements for genomic data. Building upon the Web Component framework, we developed various elements to help both app developers and users to easily integrate data visualization into their workflow. Since these components can be integrated with any web-based application using minimal programming experience, the library reduces the effort to visualize and create applications for genomic data sets.

Epiviz web components are broadly classified into the following categories - [Visualization Components](/visualization-components.md), [App Components](/app-components.md), [Data Components](/data-components.md)

---

### Installation {#install}

Epiviz Web Components are built on the Web Component framework using the [Google Polymer](https://www.polymer-project.org/) library to support different browsers.

The components are available from three different GitHub repositories

To install using [bower](https://bower.io/)

```bash
bower install epiviz/epiviz-chart ## visualization & app components
bower install epiviz/epiviz-data-source ## data components
bower install epiviz/epiviz-workspace ## workspace components
```

---

### Run {#run}

To run the demo's, install the polymer-cli tool available at [Polymer Installation](https://www.polymer-project.org/1.0/start/). Simple run

```bash
polymer serve
```

and visit [http://localhost:8081/components/epiviz-charts.](http://localhost:8081/components/epiviz-charts)

**Note: port number might be different. please check the console after running polymer to browse to the right location.**

