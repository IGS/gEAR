### Data Components

epiviz-data-source component provides functionality for the epiviz view components to interact with an active web server or an R WebSocket connection. Data components require the location \(provider-url\) attribute that specifies where the WebSocket/web server is located. When the user interacts with the epiviz components, for example adding a new visualization or navigating to a genomic region, these interactions generate data requests that propagate to the data source element.

##### WebServer Data Provider

We developed a Python Flask based data provider \(epiviz-data-provider\) to respond to data requests. The Flask API queries a MySQL database to respond to data requests. The data provider enables scaling where we bin small regions together and average the measurements. Although this takes a bit longer to respond to data queries, we see a significant improvement in draw times of charts. We also implemented data import functions for commonly used Bioconductor datatypes like GenomicRanges, SummarizedExperiment, etc., in our [R/Bioconductor epivizrData](https://bioconductor.org/packages/release/bioc/html/epivizrData.html) package.

##### WebSocket Data Provider

Epiviz components are integrated with Bioconductor data types. This enables to have an interactive R-session and visualize Bioconductor data objects using the chart components. The [epivizrChart package](https://bioconductor.org/packages/release/bioc/html/epivizrChart.html) is an API to interactively visualize R/Bioconductor data objects.

```html
<epiviz-data-source 
        provider-type="epiviz.data.WebServerDataProvider" 
        provider-id="umd" 
        provider-url="http://epiviz-dev.cbcb.umd.edu/api/">
</epiviz-data-source>
```

---

### Workspace Component

epiviz-workspace component is built upon the Google Firebase infrastructure to manage user authentication, create shareable and reproducible workspaces. Workspace components are easily reconfigurable and allows developers to customize this component to their firebase instance.Since the suite of features as part of the web components framework are still being implemented across different browsers, we use the Google Polymer library to provide polyfills for missing implementations.

To setup workspace component for a new application, register for an app at [Google Firebase](https://firebase.google.com/)

```html
<epiviz-workspace 
        auth-domain="<APP_AUTH_DOMAIN>"
        database-url="<APP_DATABASE_URL>"
        api-key="<APP_API_KEY>"
        storage-bucket="<APP_STORAGE_BUCKET>"
        messaging-sender-id="<APP_MESSAGING_SENDER_ID>">
</epiviz-workspace>
```



