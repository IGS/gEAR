var chromosomeFactory = cyto_chr.chromosome;

var x = chromosomeFactory()
  .segment("X")
  .target('#chrX')
  .resolution(550)
  .showAxis(true)
  .render();

var Y = chromosomeFactory()
  .segment("Y")
  .resolution(550)
  .target('#chrY')
  .showAxis(true)
  .selectionMode('multi')
  .render();

for (var i = 1; i < 23; i++) {
  var t = "#chr" + i;

  chromosomeFactory()
    .segment(i)
    .resolution(550)
    .target('#chr' + i)
    .showAxis(true)
    .render();
}
