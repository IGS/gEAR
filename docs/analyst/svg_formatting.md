
### SVG formatting
SVGs can be used within the gEAR as a display type to color elements of the SVG by expression.

#### Element types colorized

SVGs can be made up of a wide variety of XML elements but only the following are considered to be elements of the image to be colorized:

 - path
 - circle
 - rect
 - ellipse

#### Linking elements to expression values

The class attribute of each element is the key to linking that portion of the image to its associated expression value so that it gets colorized.  So if you have this:

    <rect class="24h_treated" height="144.72" width="46.639999" x="264.34" y="280.26999" />

Then your dataset should have an observation (a sample, i.e. a row of the observations table) named "24h_treated".  These must match exactly.  Since this colorization is based on the class names, it means multiple elements in the SVG can have the same class name, colorizing them all together.

#### Required SVG elements/attributes

This namespace must be included in the header:

     xmlns="http://www.w3.org/2000/svg"

The root `<svg>` element must also have a `viewBox` attribute, and the viewBox must contain the whole drawing. gEAR scales the viewBox to fit the display panel, so anything drawn outside it spills past the panel edges or is cut off, in both the normal and expanded views. More on how `viewBox` controls scaling [here](https://css-tricks.com/scale-svg/).

#### Fitting the page to the drawing

Drawing programs write the `viewBox` from the document's page (Inkscape) or artboard (Illustrator), not from what you drew. If shapes or labels sit off the page, they're outside the viewBox. This happens easily in Inkscape, whose default page is A4: you draw a wide figure, it runs past the page edges, and the saved file still describes only the A4 page. For example:

    <svg width="210mm" height="297mm" viewBox="0 0 210 297">
      <rect class="alpha" x="441.4" y="61.3" width="32.6" height="32.2" />
      <rect class="Neurogen" x="-197.4" y="60.0" width="32.4" height="32.6" />

The page is 210 units wide, but one box starts at x = 441 (right of the page) and another at x = -197 (left of it). In gEAR this drawing appears too large, with its left and right columns cut off.

To fix it, resize the page to the drawing before you save:

- **Inkscape 1.2 and later:** open File > Document Properties and click "Resize to content".
- **Earlier Inkscape 1.x:** with nothing selected, choose Edit > Resize Page to Selection (Shift+Ctrl+R).
- **Adobe Illustrator:** choose Object > Artboards > Fit to Artwork Bounds, then export as SVG.

Resizing the page doesn't change any `class` attributes, so the coloring still works.

**Check before sending:** open the SVG file in a web browser. The browser shows only what's inside the viewBox, which is what gEAR will show. If any part of the drawing is missing or cut off at an edge, the page still doesn't cover the drawing.

#### Text

Use regular text objects for labels. In Inkscape, text made by dragging out a text box is "flowed text", which uses an SVG 2 feature (`shape-inside`) that browsers don't support, so line wrapping and alignment may not match what you see in Inkscape. To convert it, select the text and choose Text > Convert to Text.

#### Including elements describing the mean
Usually each element in the SVG is an experimental condition or single cell ID.  It is possible though to create a drawing with these types of values AND aggregate image elements.  For example, if you have 10 Deiter 3 cells but then want one portion of the image to represent the mean, this is accomplished using a class naming syntax like this:

 - Deiter3--1
 - Deiter3--2
 - ...
 - Deiter3--10
 - Deiter3--mean

For the `--mean` element to be colored, the dataset needs a `cell_type` observation column: gEAR computes the mean expression of all observations with each cell type and makes it available under the name `<cell_type>--mean` (here `Deiter3--mean`). The individual `Deiter3--1` ... elements are colored by the observations with those exact names.

The SVG file is attached to a dataset by the gEAR team (it is stored as `www/datasets_uploaded/<dataset_id>.svg`); there is no upload control for it. See also [SVG (anatomical image) displays](../wiki/gEARWiki.md#svg-anatomical-image-displays) in the user guide.
