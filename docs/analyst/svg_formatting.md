
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

    <rect class="24h_treated" height="144.72" width="46.639999" x="264.34" y="280.26999" \>

Then your dataset should have an observation column called "24h_treated".  These must match exactly.  Since this colorization is based on the class names, it means multiple elements in the SVG can have the same class name, colorizing them all together.

#### Required SVG elements/attributes

This namespace must be included in the header:

     xmlns="http://www.w3.org/2000/svg"

Also, to ensure the SVG scales in size with the rest of the interface the header must contain a viewbox attribute.  More on the importance of that [here](https://css-tricks.com/scale-svg/).

#### Including elements describing the mean
Usually each element in the SVG is an experimental condition or single cell ID.  It is possible though to create a drawing with these types of values AND aggregate image elements.  For example, if you have 10 Deiter 3 cells but then want on portion of the image to represet the mean, this is accomplished using a class naming syntax like this:

 - Deiter3--1
 - Deiter3--2
 - ...
 - Deiter3--10
 - Deiter3--mean

The interface sees the label before the -- symbol (Deiter3) and then calculates the mean of all those image elements which have that base name.
