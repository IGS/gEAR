(function(cyto_chr, d3){
  if (typeof angular === 'undefined') {
    return;
  }

  cyto_chr.modelLoader.setDataDir('./node_modules/cyto-chromosome-vis/data/');

  angular.module('cyto-chromosome-vis',[])
    .directive('cytochromosome',[function() {
      function link(scope, element, attr) {

        attr.resolution = cyto_chr.setOption(attr.resolution, "550");
        attr.width = cyto_chr.setOption(attr.width, 1000);
        attr.segment = cyto_chr.setOption(attr.segment, "1");
        attr.useRelative = cyto_chr.setOption(attr.useRelative, true);
        attr.showAxis = cyto_chr.setOption(attr.showAxis, false);

        cyto_chr.chromosome()
          .target(d3.select(element[0]))
          .width(attr.width)
          .segment(attr.segment)
          .resolution(attr.resolution)
          .useRelative(attr.useRelative == "true")
          .showAxis(attr.showAxis == "true")
          .render();
      }

      return {
        link: link,
        restrict: 'E'
      };
    }])
    .provider('cytochromosome', function(){
      this.build = function() {
        return cyto_chr.chromosome();
      };

      this.setDataDir = function(d) {
        cyto_chr.modelLoader.setDataDir(d);
      };

      this.$get = function() {
        return this;
      };
    });

})(cyto_chr || {}, d3);