# Installation

`bower install epiviz/epiviz-chart`

# Documentation

run a local instance of polymer-server
`polymer serve`

Then navigate to http://localhost:8080/components/epiviz-chart/

# Demo

run a local instance of polymer-server
`polymer serve`

Then navigate to http://localhost:8080/components/epiviz-chart/demo/


# Update chartSettings:

```
# get chart
chart = document.querySelector("#chart1");
# get current chart settings
currentSettings = chart.ChartSettings;
# modify chart settings
...

# set settings back to chart
chart.setAttribute("chart-settings", JSON.stringify(currentSettings));
```

# Epiviz-environment

### must use polymer api to add charts to environment. Js dom api does not properly initialize elements

for example

if the page contains

```
<epiviz-environment id="env">
</epiviz-environment>
```

to add an epiviz chart for example line-track.

# create a new element
```
elem = document.createElement('epiviz-scatter-plot'); 
elem.dimS = ['affy1', 'affy2']; 
elem.className="charts"
```

# query dom for environment
`ot = document.querySelector('#env')`

# add chart
`Polymer.dom(ot).appendChild(elem)`


# Optimize elements for productions. 
```
npm install -g polymer-bundler

polymer-bundler --inline-scripts --inline-css --strip-comments epiviz-charts.html > dist/epiviz-charts.html
```



