
# gEAR plugin support

You can create your own plugins to add or even modify gEAR functionality.  These can operate on one or many pages across the site.

## Plugin architecture

### File structure and naming conventions

Each plugin is made up of a single directory within the www/plugins directory.  For any gEAR page your plugin will operate on, it is expected there will be three files with the same basename.  For this document we'll consider a plugin called *extrabuttons* which adds features to the index.html and expression.html pages.  The plugin directory should have this structure:

```text
extrabuttons/
├── expression.css
├── expression.html
├── expression.js
├── index.css
├── index.html
└── index.js
```

Each of these files needs to exist for each page your plugin operates, even if they're empty.

It's worth looking at the existing plugins to get a sense of how things work:

- `www/plugins/hrp_landing_tab/` and `www/plugins/aro_2026_tab/` - add content to `index.html`
- `www/plugins/deafness_gene_annotation/` - adds annotations to `expression.html`, with its own JSON data files

### How plugins are loaded

Loading is handled by `loadPlugins()` / `loadPlugin()` in `www/js/common.v2.js`, called during common page initialization. For the current page (`location.pathname`, defaulting to `index.html`), every plugin whose `enabled_plugins` entry lists that page is loaded: `plugins/<name>/<page>.html` is fetched into a new div appended to `<body>`, then `<page>.css` and `<page>.js` are appended to `<head>`.

#### The HTML file

It is expected that each of your plugin HTML files is an HTML snippet rather than a complete document.  The best practice is to make it a single div element with an ID you can reference as needed.

When a page is loaded the contents of that HTML file are added to a div at the end of the page and is accessible to your plugin JS file using the ID convention:

```text
  plugin_name + "_html_c"
```

So if I'm on the index page and my plugin is called extrabuttons I can add code to my index.js which can access the HTML elements via a selector like "#extrabuttons_html_c".  You can then move the content to whatever place in the DOM it belongs.

You should also take care to give unique, descriptive IDs and class attributes to your elements so they don't conflict with others. Prefixing these with your plugin name is a best practice.

#### The CSS file

Vanilla CSS file here just gets imported as a `<link>`, no special considerations.

#### The Javascript file

Add the Javascript here needed for your plugin. It is loaded as a classic (non-module) `<script>` after the page has loaded, so there is no need to wrap it in a document ready handler. The HTML snippet is inserted before the JS file is appended, so its elements are available when your script runs.

#### Other files

There are no restrictions to your adding other files (such as JSON data files) within your plugin directory as long as you handle their import/processing within your own plugin's javascript.

### Site preferences configuration

Each gEAR portal has a site preferences JSON file at www/site_domain_prefs.json.  In this file you'll need to create entries to inform gEAR about your plugin and on which pages it operates.  For our example, we'd add this:

```text
"enabled_plugins": {"extrabuttons": ["index.html", "expression.html"]}
```

## Submitting your plugin

If you have written a plugin you'd like to share with the community please submit a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests).

