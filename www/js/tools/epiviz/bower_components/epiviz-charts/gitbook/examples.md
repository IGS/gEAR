### Demo

---

### Embedding App Components

```html
<epiviz-data-source provider-type="epiviz.data.WebServerDataProvider" provider-id="umd" provider-url="http://epiviz-dev.cbcb.umd.edu/api/"></epiviz-data-source>

<epiviz-environment chr="chr6" start=44076201 end=45076201 no-logo measurements='{"affy1":{
                    "id":"cancer",
                    "name":"Expression Colon Cancer",
                    "type":"feature",
                    "datasourceId":"gene_expression",
                    "datasourceGroup":"affymetrix_probeset",
                    "dataprovider":"umd",
                    "formula":null,
                    "defaultChartType":null,
                    "annotation":null,
                    "minValue":-3,
                    "maxValue":20,
                    "metadata":["probe"]
                },
                "affy2":{
                    "id":"normal",
                    "name":"Expression Colon Normal",
                    "type":"feature",
                    "datasourceId":"gene_expression",
                    "datasourceGroup":"affymetrix_probeset",
                    "dataprovider":"umd",
                    "formula":null,
                    "defaultChartType":null,
                    "annotation":null,
                    "minValue":-3,
                    "maxValue":20,
                    "metadata":["probe"]
                },
                "genes": {
                    "id": "genes",
                    "name": "Genes",
                    "type": "range",
                    "datasourceId": "genes",
                    "datasourceGroup": "genes",
                    "dataprovider": "umd",
                    "formula": null,
                    "defaultChartType": "Genes Track",
                    "annotation": null,
                    "minValue": null,
                    "maxValue": null,
                    "metadata": ["gene", "entrez", "exon_starts", "exon_ends"]
                }
              }'>
        <epiviz-genes-track id="test" class="charts" dim-s='["genes"]'></epiviz-genes-track>
        <epiviz-scatter-plot class="charts" dim-s='["affy1", "affy2"]'></epiviz-scatter-plot>
</epiviz-environment>
```

---

### Embedding Charts

#### Genes Track

```
<epiviz-genes-track use-default-data-provider="true"></epiviz-genes-track>
```

#### Line Track

```
<epiviz-line-track use-default-data-provider="true"></epiviz-line-track>
```

#### Scatter Plot

```
<epiviz-scatter-plot use-default-data-provider="true"></epiviz-scatter-plot>
```



