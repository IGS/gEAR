# Installation

`bower install hcbravolab/epiviz-data-source`

# Documentation
run a local instance of polymer-server
`polymer serve`

Then navigate to http://localhost:8080/components/epiviz-data-source/

# Demo

run a local instance of polymer-server
`polymer serve`

Then navigate to http://localhost:8080/components/epiviz-data-source/demo/

# Access element DataManger

```
var elem = document.querySelector('epiviz-data-source');

elem.dataManager 

Use any available functions following the epiviz API
// get Measurements from dataManager
callback = function(data) {console.log(data);};

elem.dataManager.getMeasurements(callback);
```

