### epiviz-navigation

Navigation app component with genomic context linked to it. Navigation elements provide UI functionality to search for a gene/probe and update the location to a genomic region. In addition, `<epiviz-navigation>` has a collapse and expand feature. When collapsed, it hides its children and provides a smaller compact ideogram-view of the current genomic location of the element. When expanded, features such are navigating along the chromosome, or navigating to a gene location are available.

To create an navigation element on a HTML page, add

```html
<epiviz-navigation
  chr="chr11"
  start=1
  end=10000>
    <epiviz-genes-track></epiviz-genes-track>
    <epiviz-scatter-plot></epiviz-scatter-plot>
</epiviz-navigation>
```



