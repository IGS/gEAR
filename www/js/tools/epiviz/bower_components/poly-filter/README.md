
# Data Filter (`data-filter`) Polymer Component

A fast and customizable polymer 2.0 component for filtering arrays. 

### Features
- Automatically detect Attribute types from the `data` array.
- Client side filtering of large arrays
- Flattens `data` for easy viewing as table
- Token based filtering of attributes
- Display currently added filters
- Customizable filter components (see section below - types of filters)
- Supports Manual/Auto selection of rows
- If `Auto` selected, supports configurable selection Modes
  - Top - select top `x%` of rows
  - Random - selection random `x%` of rows

### Types of Filters

Attribute Type | Component | Component url
---------------|-----------|--------------
Categorical | `paper-tags` | [webcomponents.org](https://www.webcomponents.org/element/PolymerEl/paper-tags)
Numerical | `paper-rangle-slider` | [webcomponents.org](https://www.webcomponents.org/element/IftachSadeh/paper-range-slider)
Boolean | `paper-toggle-button` | [webcomponents.org](https://www.webcomponents.org/element/@polymer/paper-toggle-button)

### Adding a data filter component to an application
In typical use, just add `<data-filter>` to the HTML page:

     <data-filter data="[]"></data-filter>

## Installation

Install the component using `bower`

    bower install epiviz/table-filter --save

## Install the `Polymer-CLI`

First, make sure you have the [Polymer CLI](https://www.npmjs.com/package/polymer-cli) installed. Then run `polymer serve` to serve the element locally.

## License - MIT