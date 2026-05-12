# Webpage dependencies

This is just some scratch notes about page elements that need to load or be processed first before other things.  These can be addressed in a few different ways:

1. Using async/await or function () => {}.then () => {}
2. Creating event listeners (vanilla JS or jQuery) and then triggering these events

Many issues can relate to the speed of execution.  The `common.js` module calls `$(document).ready`, which executes when the DOM is ready and, more importantly, before `window.onload()` found in many of the `<page>.js` scripts.  But if the page has little to no content like images, styles, embedded stuff, etc. then there is a good chance that CGI or API calls executed from that context finish before those from `$(document).ready`. So it is important to determine what code within `window.onload()` is a dependency of something in `$(document).ready` and needs to be moved out of there.

Also, it's worth pointing out that `$(function() {})` and `$(() => {})` are both shorthand for `$(document).ready` and may appear on other pages.

## Common

These cases will appear on all or most pages:

### Login

1. The navigation bar loads
2. Login check occurs.  This changes the navigation bar from inputting user/pass to showing logged in user.
    1. A false logged-out display could indicate that the user logged in before the navigation bar loaded.

### JSTree objects

These pertain to retrieving the data to load the JSTree objects

1. Login check.  User's session ID will be returned after successful login or cookie retrieval
2. Session ID is applied in CGI script calls to retrieve dataset/profile/genecart data, so that user-specific data will also be retrieved in addition to public data
    1. JSTree objects with only public data populated could indicate the session ID was not retrieved yet.

## Index (front page)

### Setting the layout (AKA profile)

1. Profile JSTree object is retrieved, and `DatasetCollectionPanel.set_layout()` is run.  This sets the active layout ID.
2. Depending on if URL parameters were provided, the search button may be "clicked", which triggers `search_genes.py` and eventual drawing of plots for each display for the active layout ID. If there is no active layout ID, then the default layout ID is used.
    1. Seeing displays from the default layout ID when another profile is selected could mean that the profile tree (and active layout ID) was not loaded before the `windows.onload()` event started (which triggers the clicking of the search button if either gene symbols or a gene cart is in the URL).

TO BE CONTINUED