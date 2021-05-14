[![Published on webcomponents.org](https://img.shields.io/badge/webcomponents.org-published-blue.svg)](https://www.webcomponents.org/element/polymerEl/paper-tags)

# \<paper-tags\>

Series of elements for handling tags with Polymer and custom elements:

- `<paper-tags>` - just the tags 
- `<paper-tags-input>` - input element for entering tags
- `<paper-tags-dropdown>` - dropdown-elements for selecting tags in a dropdown-list

<div>
	<img src="https://raw.githubusercontent.com/PolymerEl/paper-tags/master/images/paper-tags.png" width="600"></img>
</div>

Example Usage:

<!--
```
<custom-element-demo>
  <template>
    <link rel="import" href="paper-tags.html">
  	<link rel="import" href="paper-tags-input.html">
  	<link rel="import" href="paper-tags-dropdown.html">
    <next-code-block></next-code-block>
  </template>
</custom-element-demo>
```
-->
```html
<paper-tags-input label="input label" show-counter="tags" items='["hello", "new"]'  maxLength="10"></paper-tags-input>
<paper-tags-dropdown items='[{"id":1,"desc":"description 1","label":"Title 1"},{"id":"12345","desc":"description 4","label":"new items"},{"id":"5","desc":"description 4","label":"Hello"},{"id":6,"desc":"description 4","label":"Bar"}]' noink label="label dropdown"  value-object='{"5": "true"}' ></paper-tags-dropdown>


```

