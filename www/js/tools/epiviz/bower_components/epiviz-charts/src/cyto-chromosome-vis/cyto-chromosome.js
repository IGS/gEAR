(function() {
  var cyto_chr = {}, d3 = window.d3;

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = cyto_chr;
    d3 = require('d3');
  } else if(typeof angular !== 'undefined') {} else {
    window.cyto_chr = cyto_chr;
  }
  
(function(cyto_chr, d3) {
  'use strict'
  cyto_chr.margin = {
    top: 38,
    left: 14,
    right: 5
  };

  var CHR1_BP_END = 248956422;
  var CHR1_BP_MID = 121700000;

  var Chromosome = function() {
    this._segment = '1';
    this._selectionMode = 'single';
    this._domTarget = d3.select(document.documentElement);
    this._resolution = "550";
    this._width = null;
    this._height = 17;
    this._useRelative = true;
    this._showAxis = false;
    this.dispatch = d3.dispatch('bandclick', 'selectorchange', 'selectorend', 'selectordelete', "selectorhover", "selectorunhover");
    this.rendered = false;
    this.selectors = [];
    this.model = [];
    this.maxBasePair = 0;
    this.xscale = d3.scale.linear();
  };

  Chromosome.prototype.getMaxBasepair = function() {
    return this.maxBasePair;
  };

  Chromosome.prototype.segment = function (a) {
    if (typeof a === 'number') a = a.toString();
    return cyto_chr.InitGetterSetter.call(this, '_segment', a);
  };

  Chromosome.prototype.selectionMode = function (a) {
    return cyto_chr.InitGetterSetter.call(this, '_selectionMode', a);
  };

  Chromosome.prototype.target = function (a) {
    // if(typeof a === 'string') 
    a = d3.select(a);
    // if(a.empty()) {
    //   throw "Error: Invalid dom target";
    // }
    return cyto_chr.InitGetterSetter.call(this, '_domTarget', a);
  };

  Chromosome.prototype.resolution = function (a) {
    if (typeof a === 'number') a = a.toString();
    if (a === "400" || a === "550" || a ==="850" || a === "1200" || a === undefined) {
      return cyto_chr.InitGetterSetter.call(this, '_resolution', a);
    } else {
      throw "Error: Invalid resolution. Please enter 400, 550, 850, or 1200 only.";
    }
  };

  Chromosome.prototype.width = function (a) {
    return cyto_chr.InitGetterSetter.call(this, '_width', a);
  };

  Chromosome.prototype.height = function (a) {
    return cyto_chr.InitGetterSetter.call(this, '_height', a);
  };

  Chromosome.prototype.useRelative = function (a) {
    return cyto_chr.InitGetterSetter.call(this, '_useRelative', a);
  };

  Chromosome.prototype.showAxis = function (a) {
    return cyto_chr.InitGetterSetter.call(this, '_showAxis', a);
  };

  Chromosome.prototype.doHover = function() {
    // console.log(this.selectors[0]);
    this.svgTarget.selectAll(".extRect rect").style("stroke", "#ffc600").style("stroke-width", "5").style("stroke-opacity", "1.0");
  };

  Chromosome.prototype.doUnhover = function() {
    this.svgTarget.selectAll(".extRect rect").style("stroke-opacity", "0");
  };

  Chromosome.prototype.on = function(e, listener) {
    if (!this.dispatch.hasOwnProperty(e)) throw "Error: No event for " + e;
    this.dispatch.on(e, listener);
  };

  Chromosome.prototype.config = function(type, arg) {
    return this[type](arg);
  };

  Chromosome.prototype.renderAxis = function () {
    var bpAxis = d3.svg.axis()
      .scale(this.xscale)
      .tickFormat(d3.format('s'))
      .orient("bottom");

    if (this._useRelative && (this._segment === "Y" || this._segment === "22" || this._segment === "21" || this._segment === "20" || this._segment === "19")) {
      bpAxis.ticks(6);
    }

    var axisg = this.svgTarget.append('g')
      .classed('bp-axis', true)
      .attr('transform', 'translate('+ cyto_chr.margin.left + ',' + (this._height + cyto_chr.margin.top + 6) + ")");

      axisg.call(bpAxis);

    axisg.selectAll('text')
      // .attr('transform', 'rotate(-50)')
      .style('font', '8px sans-serif');

    axisg.selectAll('path, line')
      .style({
        "fill": "none",
        "stroke": "#666666",
        "shape-rendering": "crispEdges"
      });
  };

  Chromosome.prototype.remove = function() {
    this.svgTarget.remove();
  };

  // Chromosome.prototype.moveSelectorTo = function(start, stop) {
  //   if(arguments.length !== 2) {
  //     throw "Error moveSelectorTo: Invalid number of arguments. Both start and stop coordinates are required";
  //   }

  //   if (this.selectors.length === 0) {
  //     this.newSelector(0, "10000000");
  //   } else {
  //     this.selectors[0].move(start, stop);
  //   }

  // };

  Chromosome.prototype.newSelector = function(bp_start, bp_stop) {

    var self = this;

    function selectorRemoveCB(sel) {
      self.dispatch.selectordelete(sel);
      var index = self.selectors.indexOf(sel);
      self.selectors.splice(index, 1);
    }

    var ve = cyto_chr.selector(selectorRemoveCB)
      .x(cyto_chr.margin.left)
      .y(cyto_chr.margin.top - (this._height / 4))
      .height(this._height + (this._height / 2))
      .xscale(this.xscale)
      .extent([bp_start, bp_stop])
      .target(this.svgTarget)
      .render();

    // ve.dispatch.on('change', function(d) {
    //   self.dispatch.selectorchange(d);
    // });

    // ve.dispatch.on('changeend', function(d) {
    //   self.dispatch.selectorend(d);
    // });

    ve.dispatch.on('selectorhover', function(d) {
      self.dispatch.selectorhover(d);
    });

    ve.dispatch.on('selectorunhover', function(d) {
      self.dispatch.selectorunhover(d);
    });

    this.selectors.push(ve);
  };

  Chromosome.prototype.getSelections = function() {

    var ret = [];
    for(var i = 0; i < this.selectors.length; i++) {
      var sel = this.selectors[i].extent();
      ret.push({
        start: sel[0],
        stop: sel[1]
      })
    }
    return ret;
  };

  Chromosome.prototype.getSelectedBands = function(sensitivity) {

    var results = [];
    if (this.selectors.length > 0) {

      var se = this.selectors[0].extent();
      var selStart = +se[0];
      var selStop = +se[1];

      if (typeof sensitivity !== 'undefined') {
        selStart -= sensitivity;
        selStop += sensitivity;
      }

      results = this.model.slice().filter(function(e) {
        var bStart = +e.bp_start;
        var bStop = +e.bp_stop;

        if ((selStart >= bStart && selStart < bStop) ||
          (selStop > bStart && selStop <= bStop) ||
          (selStart <= bStart && selStop >= bStop)) {
          return true;
        } else {
          return false;
        }
      });
    }

    return results;
  };

  Chromosome.prototype.getSVGTarget = function() {
    return this.svgTarget;
  };

  Chromosome.prototype.render = function (zoomRange) {

    var self = this;

    if(self.rendered) {
      self.selectors = [];
      self.remove();
    }

    if(self._width === null) {
      var parentWidth = d3.select(self._domTarget[0][0].parentNode).node().getBoundingClientRect().width;
      self.width(parentWidth)
    }

    cyto_chr.modelLoader.load(this._segment, this._resolution, function(data) {
      self.model = data;
      self.maxBasePair = d3.max(data, function(d) {
        return +d.bp_stop;
      });

      self.segMid = 0;
      for(var j =0; j < data.length;j++) {
        if(data[j].stain ==="acen") {
          self.segMid = data[j].bp_stop;
          break;
        }
      }

      var rangeTo = self._useRelative ? (self.maxBasePair / CHR1_BP_END) * self._width : self._width;
      rangeTo -= (cyto_chr.margin.left + cyto_chr.margin.right);

      if (zoomRange) {
        self.xscale.domain([zoomRange[0], zoomRange[1]]);
      } else {
        self.xscale.domain([1, self.maxBasePair]);
      }

      self.xscale.range([0, rangeTo]);

      var h = self._height + 62;
      // var h = self._height;
      var w = self._width;

      self.svgTarget = self._domTarget
        .style('height', h + 'px')
        .style('width', (w+5) + 'px')
        .append('svg')
        .attr('width', w)
        .attr('height', h);

      var bands = self.svgTarget.selectAll('g')
        .data(data).enter();

      // self.svgTarget.append('text')
      //   .text("chr: " + self._segment)
      //   .attr('x', w - 77)
      //   .attr('y', cyto_chr.margin.top + (self._height/ 3) + 2)
      //   .attr('text-anchor','middle')
      //   .style('font', '12px sans-serif');

      function bpCoord(bp) {
        var xshift = 0;
        if(self.alignCentromere && self._segment !== "1") {
          xshift = self.xscale(CHR1_BP_MID) - self.xscale(self.segMid);
        }

        return self.xscale(bp) + cyto_chr.margin.left + xshift;
      }

      bands.append('g')
        .each(function(d, i) {

          var elem = d3.select(this);

          function applyBorder() {
            this
              .attr('stroke', '#000000')
              .attr('stroke-width', 0.2);
          }

          function drawRoundedRect(d, r, tl, tr, bl, br) {
            return this.append('path')
              .attr("d", cyto_chr.roundedRect(bpCoord(d.bp_start), cyto_chr.margin.top, bpCoord(d.bp_stop) - bpCoord(d.bp_start), self._height, r, tl, tr, bl, br))
              .style('fill', cyto_chr.getStainColour(d.stain, d.density));
          }

          var labelSkipFactor = self._resolution === '1200' ? 8 : 2;

          if (self._useRelative && self._resolution == "1200" && self._segment == 'Y') {
            labelSkipFactor = 12;
          }

          // if(i % labelSkipFactor === 0) {
          //   var bmid = (bpCoord(d.bp_stop) + bpCoord(d.bp_start)) / 2;
          //   elem.append('line')
          //     .attr('x1', bmid)
          //     .attr('y1', cyto_chr.margin.top)
          //     .attr('x2', bmid)
          //     .attr('y2', cyto_chr.margin.top - 4)
          //     .style('stroke', 'grey')
          //     .style('stroke-width',1);

          //   elem.append('text')
          //     .attr('transform', 'translate(' + bmid + ',' + (cyto_chr.margin.top - 6) + ')rotate(-50)')
          //     .style('font', '8px sans-serif')
          //     .text(d.arm + d.band);
          // }

          var rect;
          var w = bpCoord(d.bp_stop) - bpCoord(d.bp_start);
          var acenThreshold = (self._resolution == "1200") ? 7 : 6;
          if (i === 0 && w > 10) {
            rect = drawRoundedRect.call(elem, d, 4, true, false, true, false);
            applyBorder.call(rect);
          } else if (d.stain === "acen" && (w > acenThreshold)) {

            if (d.arm === "p") {
              rect = drawRoundedRect.call(elem, d, 5, false, true, false, true);

            } else if(d.arm === "q") {
              rect = drawRoundedRect.call(elem, d, 5, true, false, true, false);
            }
          } else if (i === data.length - 1) {

            rect = drawRoundedRect.call(elem, d, 5, false, true, false, true);
            applyBorder.call(rect);

          } else {

            var ys = d.stain === "stalk" ? cyto_chr.margin.top + (self._height / 4) : cyto_chr.margin.top;
            var hs = d.stain === "stalk" ? self._height / 2 : self._height;
            rect = elem.append('rect')
              .attr('x', bpCoord(d.bp_start))
              .attr('y', ys)
              .attr('height', hs)
              .attr('width', self.xscale(d.bp_stop) - self.xscale(d.bp_start))
              .style('fill', cyto_chr.getStainColour(d.stain, d.density));
            applyBorder.call(rect);
          }

          rect.append('title')
            .text(d.arm + d.band)

          rect.on('mouseover', function(d) {
            var e = d3.select(this)
              .style('opacity', "0.5")
              .style('cursor', 'pointer');

            if (d.stain === "gneg") {
              e.style('fill', cyto_chr.getStainColour("gpos", "25"));
            }

          });

          rect.on('mouseout', function(d) {
            var e = d3.select(this)
              .style('opacity', "1")
              .style('cursor', 'default');

            if (d.stain === "gneg") {
              e.style('fill', cyto_chr.getStainColour("gneg"));
            }
          });

          // rect.on('click', function(d) {
          //   // if (self.selectors.length === 0 || (self._selectionMode === 'multi' && d3.event.altKey)) {
          //   //   self.newSelector(d.bp_start, d.bp_stop);
          //   // }

          //   // if (self._selectionMode === 'single' && self.selectors.length > 0) {
          //   //   self.moveSelectorTo(d.bp_start, d.bp_stop);
          //   // }
          //   // self.dispatch.bandclick(d);
          // });
        });

      if (self._showAxis) {
        self.renderAxis();
      }

      self.rendered = true;
    });

    return self;
  };

  cyto_chr.chromosome = function() {
    return new Chromosome();
  };

})(cyto_chr || {}, d3);

(function (cyto_chr, d3) {

  var defaultDataURLs = {
    "400" : "ideogram_9606_GCF_000001305.14_400_V1",
    "550" : "ideogram_9606_GCF_000001305.14_550_V1",
    "850" : "ideogram_9606_GCF_000001305.14_850_V1",
    "1200" : "ideogram_9606_GCF_000001305.13_1200_v1"
  };

  var baseDir = 'data/';

  function CacheInstance() {
    this.status = "notloaded";
    this.cache = [];
  }

  var dataCache = {
    "400" : new CacheInstance,
    "550" : new CacheInstance,
    "850" : new CacheInstance,
    "1200" : new CacheInstance
  };

  var callQueue = [];

  function loadData(file, res, cb) {

    var c = dataCache[res];
    if (c.cache.length === 0) {
      if (c.status === "loading") {
        callQueue.push({
          res: res,
          cb: cb
        });

        return;
      } else if (c.status === "notloaded") {
        c.status = "loading";
        d3.tsv(file, function(d) {
          c.cache = d;
          c.status = "loaded";
          cb(d);

          while(callQueue.length > 0) {
            var cbq = callQueue.shift();
            cbq.cb(d);
          }

        });
      }
    } else {
      cb(c.cache);
    }
  }

  function getChromosomeData(chr, resolution, cb) {

    chr = chr || '1';
    resolution = resolution || "550";

    var fileName = defaultDataURLs[resolution];

    var data = [{"#chromosome":"1","arm":"p","band":"36.3","iscn_start":"0","iscn_stop":"451","bp_start":"1","bp_stop":"7100000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"36.2","iscn_start":"451","iscn_stop":"682","bp_start":"7100001","bp_stop":"15900000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"p","band":"36.1","iscn_start":"682","iscn_stop":"1259","bp_start":"15900001","bp_stop":"27600000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"35","iscn_start":"1259","iscn_stop":"1583","bp_start":"27600001","bp_stop":"34300000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"p","band":"34.3","iscn_start":"1583","iscn_stop":"1779","bp_start":"34300001","bp_stop":"39600000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"34.2","iscn_start":"1779","iscn_stop":"1999","bp_start":"39600001","bp_stop":"43700000","stain":"gpos","density":"25"},{"#chromosome":"1","arm":"p","band":"34.1","iscn_start":"1999","iscn_stop":"2184","bp_start":"43700001","bp_stop":"46300000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"33","iscn_start":"2184","iscn_stop":"2519","bp_start":"46300001","bp_stop":"50200000","stain":"gpos","density":"75"},{"#chromosome":"1","arm":"p","band":"32.3","iscn_start":"2519","iscn_stop":"2723","bp_start":"50200001","bp_stop":"55600000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"32.2","iscn_start":"2723","iscn_stop":"2825","bp_start":"55600001","bp_stop":"58500000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"p","band":"32.1","iscn_start":"2825","iscn_stop":"3050","bp_start":"58500001","bp_stop":"60800000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"31.3","iscn_start":"3050","iscn_stop":"3387","bp_start":"60800001","bp_stop":"68500000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"p","band":"31.2","iscn_start":"3387","iscn_stop":"3705","bp_start":"68500001","bp_stop":"69300000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"31.1","iscn_start":"3705","iscn_stop":"4598","bp_start":"69300001","bp_stop":"84400000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"p","band":"22.3","iscn_start":"4598","iscn_stop":"4933","bp_start":"84400001","bp_stop":"87900000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"22.2","iscn_start":"4933","iscn_stop":"5148","bp_start":"87900001","bp_stop":"91500000","stain":"gpos","density":"75"},{"#chromosome":"1","arm":"p","band":"22.1","iscn_start":"5148","iscn_stop":"5430","bp_start":"91500001","bp_stop":"94300000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"21","iscn_start":"5430","iscn_stop":"6274","bp_start":"94300001","bp_stop":"106700000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"p","band":"13.3","iscn_start":"6274","iscn_stop":"6560","bp_start":"106700001","bp_stop":"111200000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"13.2","iscn_start":"6560","iscn_stop":"6827","bp_start":"111200001","bp_stop":"115500000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"p","band":"13.1","iscn_start":"6827","iscn_stop":"7094","bp_start":"115500001","bp_stop":"117200000","stain":"gneg","density":""},{"#chromosome":"1","arm":"p","band":"12","iscn_start":"7094","iscn_stop":"7279","bp_start":"117200001","bp_stop":"120400000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"p","band":"11","iscn_start":"7279","iscn_stop":"7394","bp_start":"120400001","bp_stop":"123400000","stain":"acen","density":""},{"#chromosome":"1","arm":"q","band":"11","iscn_start":"7394","iscn_stop":"7554","bp_start":"123400001","bp_stop":"125100000","stain":"acen","density":""},{"#chromosome":"1","arm":"q","band":"12","iscn_start":"7554","iscn_stop":"8672","bp_start":"125100001","bp_stop":"143200000","stain":"gvar","density":""},{"#chromosome":"1","arm":"q","band":"21.1","iscn_start":"8672","iscn_stop":"8915","bp_start":"143200001","bp_stop":"147500000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"21.2","iscn_start":"8915","iscn_stop":"9094","bp_start":"147500001","bp_stop":"150600000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"q","band":"21.3","iscn_start":"9094","iscn_stop":"9350","bp_start":"150600001","bp_stop":"155100000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"22","iscn_start":"9350","iscn_stop":"9685","bp_start":"155100001","bp_stop":"156600000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"q","band":"23","iscn_start":"9685","iscn_stop":"10160","bp_start":"156600001","bp_stop":"165500000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"24","iscn_start":"10160","iscn_stop":"10565","bp_start":"165500001","bp_stop":"173000000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"q","band":"25","iscn_start":"10565","iscn_stop":"11005","bp_start":"173000001","bp_stop":"185800000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"31","iscn_start":"11005","iscn_stop":"12207","bp_start":"185800001","bp_stop":"198700000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"q","band":"32.1","iscn_start":"12207","iscn_stop":"12757","bp_start":"198700001","bp_stop":"207100000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"32.2","iscn_start":"12757","iscn_stop":"12988","bp_start":"207100001","bp_stop":"211300000","stain":"gpos","density":"25"},{"#chromosome":"1","arm":"q","band":"32.3","iscn_start":"12988","iscn_stop":"13272","bp_start":"211300001","bp_stop":"214400000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"41","iscn_start":"13272","iscn_stop":"13920","bp_start":"214400001","bp_stop":"223900000","stain":"gpos","density":"100"},{"#chromosome":"1","arm":"q","band":"42.1","iscn_start":"13920","iscn_stop":"14212","bp_start":"223900001","bp_stop":"230500000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"42.2","iscn_start":"14212","iscn_stop":"14296","bp_start":"230500001","bp_stop":"234600000","stain":"gpos","density":"50"},{"#chromosome":"1","arm":"q","band":"42.3","iscn_start":"14296","iscn_stop":"14406","bp_start":"234600001","bp_stop":"236400000","stain":"gneg","density":""},{"#chromosome":"1","arm":"q","band":"43","iscn_start":"14406","iscn_stop":"14776","bp_start":"236400001","bp_stop":"243500000","stain":"gpos","density":"75"},{"#chromosome":"1","arm":"q","band":"44","iscn_start":"14776","iscn_stop":"15100","bp_start":"243500001","bp_stop":"248956422","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"25.3","iscn_start":"0","iscn_stop":"193","bp_start":"1","bp_stop":"4400000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"25.2","iscn_start":"193","iscn_stop":"282","bp_start":"4400001","bp_stop":"6900000","stain":"gpos","density":"50"},{"#chromosome":"2","arm":"p","band":"25.1","iscn_start":"282","iscn_stop":"475","bp_start":"6900001","bp_stop":"12000000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"24","iscn_start":"475","iscn_stop":"931","bp_start":"12000001","bp_stop":"23800000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"p","band":"23","iscn_start":"931","iscn_stop":"1558","bp_start":"23800001","bp_stop":"31800000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"22","iscn_start":"1558","iscn_stop":"2145","bp_start":"31800001","bp_stop":"41500000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"p","band":"21","iscn_start":"2145","iscn_stop":"2650","bp_start":"41500001","bp_stop":"47500000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"16","iscn_start":"2650","iscn_stop":"3439","bp_start":"47500001","bp_stop":"61000000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"p","band":"15","iscn_start":"3439","iscn_stop":"3682","bp_start":"61000001","bp_stop":"63900000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"14","iscn_start":"3682","iscn_stop":"3976","bp_start":"63900001","bp_stop":"68400000","stain":"gpos","density":"50"},{"#chromosome":"2","arm":"p","band":"13","iscn_start":"3976","iscn_stop":"4572","bp_start":"68400001","bp_stop":"74800000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"12","iscn_start":"4572","iscn_stop":"5068","bp_start":"74800001","bp_stop":"83100000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"p","band":"11.2","iscn_start":"5068","iscn_stop":"5534","bp_start":"83100001","bp_stop":"91800000","stain":"gneg","density":""},{"#chromosome":"2","arm":"p","band":"11.1","iscn_start":"5534","iscn_stop":"5665","bp_start":"91800001","bp_stop":"93900000","stain":"acen","density":""},{"#chromosome":"2","arm":"q","band":"11.1","iscn_start":"5665","iscn_stop":"5828","bp_start":"93900001","bp_stop":"96000000","stain":"acen","density":""},{"#chromosome":"2","arm":"q","band":"11.2","iscn_start":"5828","iscn_stop":"6181","bp_start":"96000001","bp_stop":"102100000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"12","iscn_start":"6181","iscn_stop":"6602","bp_start":"102100001","bp_stop":"108700000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"q","band":"13","iscn_start":"6602","iscn_stop":"6928","bp_start":"108700001","bp_stop":"112200000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"14.1","iscn_start":"6928","iscn_stop":"7254","bp_start":"112200001","bp_stop":"118100000","stain":"gpos","density":"50"},{"#chromosome":"2","arm":"q","band":"14.2","iscn_start":"7254","iscn_stop":"7336","bp_start":"118100001","bp_stop":"121600000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"14.3","iscn_start":"7336","iscn_stop":"7662","bp_start":"121600001","bp_stop":"129100000","stain":"gpos","density":"50"},{"#chromosome":"2","arm":"q","band":"21.1","iscn_start":"7662","iscn_stop":"7920","bp_start":"129100001","bp_stop":"131700000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"21.2","iscn_start":"7920"
                ,"iscn_stop":"8177","bp_start":"131700001","bp_stop":"134300000","stain":"gpos","density":"25"},{"#chromosome":"2","arm":"q","band":"21.3","iscn_start":"8177","iscn_stop":"8314","bp_start":"134300001","bp_stop":"136100000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"22","iscn_start":"8314","iscn_stop":"9048","bp_start":"136100001","bp_stop":"147900000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"q","band":"23","iscn_start":"9048","iscn_stop":"9374","bp_start":"147900001","bp_stop":"154000000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"24.1","iscn_start":"9374","iscn_stop":"9650","bp_start":"154000001","bp_stop":"158900000","stain":"gpos","density":"75"},{"#chromosome":"2","arm":"q","band":"24.2","iscn_start":"9650","iscn_stop":"9868","bp_start":"158900001","bp_stop":"162900000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"24.3","iscn_start":"9868","iscn_stop":"10202","bp_start":"162900001","bp_stop":"168900000","stain":"gpos","density":"75"},{"#chromosome":"2","arm":"q","band":"31","iscn_start":"10202","iscn_stop":"10854","bp_start":"168900001","bp_stop":"182100000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"32.1","iscn_start":"10854","iscn_stop":"11271","bp_start":"182100001","bp_stop":"188500000","stain":"gpos","density":"75"},{"#chromosome":"2","arm":"q","band":"32.2","iscn_start":"11271","iscn_stop":"11427","bp_start":"188500001","bp_stop":"191100000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"32.3","iscn_start":"11427","iscn_stop":"11792","bp_start":"191100001","bp_stop":"196600000","stain":"gpos","density":"75"},{"#chromosome":"2","arm":"q","band":"33","iscn_start":"11792","iscn_stop":"12362","bp_start":"196600001","bp_stop":"208200000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"34","iscn_start":"12362","iscn_stop":"12811","bp_start":"208200001","bp_stop":"214500000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"q","band":"35","iscn_start":"12811","iscn_stop":"13164","bp_start":"214500001","bp_stop":"220700000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"36","iscn_start":"13164","iscn_stop":"13653","bp_start":"220700001","bp_stop":"230100000","stain":"gpos","density":"100"},{"#chromosome":"2","arm":"q","band":"37.1","iscn_start":"13653","iscn_stop":"13930","bp_start":"230100001","bp_stop":"234700000","stain":"gneg","density":""},{"#chromosome":"2","arm":"q","band":"37.2","iscn_start":"13930","iscn_stop":"14027","bp_start":"234700001","bp_stop":"236400000","stain":"gpos","density":"50"},{"#chromosome":"2","arm":"q","band":"37.3","iscn_start":"14027","iscn_stop":"14400","bp_start":"236400001","bp_stop":"242193529","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"26","iscn_start":"0","iscn_stop":"257","bp_start":"1","bp_stop":"8100000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"p","band":"25","iscn_start":"257","iscn_stop":"649","bp_start":"8100001","bp_stop":"16300000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"24.3","iscn_start":"649","iscn_stop":"1010","bp_start":"16300001","bp_stop":"23800000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"p","band":"24.2","iscn_start":"1010","iscn_stop":"1114","bp_start":"23800001","bp_stop":"26300000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"24.1","iscn_start":"1114","iscn_stop":"1355","bp_start":"26300001","bp_stop":"30800000","stain":"gpos","density":"75"},{"#chromosome":"3","arm":"p","band":"23","iscn_start":"1355","iscn_stop":"1574","bp_start":"30800001","bp_stop":"32000000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"22","iscn_start":"1574","iscn_stop":"1966","bp_start":"32000001","bp_stop":"43600000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"p","band":"21.3","iscn_start":"1966","iscn_stop":"2910","bp_start":"43600001","bp_stop":"50600000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"21.2","iscn_start":"2910","iscn_stop":"3041","bp_start":"50600001","bp_stop":"52300000","stain":"gpos","density":"25"},{"#chromosome":"3","arm":"p","band":"21.1","iscn_start":"3041","iscn_stop":"3302","bp_start":"52300001","bp_stop":"54400000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"14.3","iscn_start":"3302","iscn_stop":"3575","bp_start":"54400001","bp_stop":"58600000","stain":"gpos","density":"50"},{"#chromosome":"3","arm":"p","band":"14.2","iscn_start":"3575","iscn_stop":"3847","bp_start":"58600001","bp_stop":"63800000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"14.1","iscn_start":"3847","iscn_stop":"4103","bp_start":"63800001","bp_stop":"69700000","stain":"gpos","density":"50"},{"#chromosome":"3","arm":"p","band":"13","iscn_start":"4103","iscn_stop":"4456","bp_start":"69700001","bp_stop":"74100000","stain":"gneg","density":""},{"#chromosome":"3","arm":"p","band":"12","iscn_start":"4456","iscn_stop":"5439","bp_start":"74100001","bp_stop":"87100000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"p","band":"11","iscn_start":"5439","iscn_stop":"5573","bp_start":"87100001","bp_stop":"90900000","stain":"acen","density":""},{"#chromosome":"3","arm":"q","band":"11.1","iscn_start":"5573","iscn_stop":"5733","bp_start":"90900001","bp_stop":"94000000","stain":"acen","density":""},{"#chromosome":"3","arm":"q","band":"11.2","iscn_start":"5733","iscn_stop":"6034","bp_start":"94000001","bp_stop":"98600000","stain":"gvar","density":""},{"#chromosome":"3","arm":"q","band":"12","iscn_start":"6034","iscn_stop":"6175","bp_start":"98600001","bp_stop":"103100000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"13.1","iscn_start":"6175","iscn_stop":"6665","bp_start":"103100001","bp_stop":"111600000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"q","band":"13.2","iscn_start":"6665","iscn_stop":"6806","bp_start":"111600001","bp_stop":"113700000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"13.3","iscn_start":"6806","iscn_stop":"7389","bp_start":"113700001","bp_stop":"122200000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"q","band":"21","iscn_start":"7389","iscn_stop":"7926","bp_start":"122200001","bp_stop":"129500000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"22","iscn_start":"7926","iscn_stop":"8255","bp_start":"129500001","bp_stop":"139000000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"q","band":"23","iscn_start":"8255","iscn_stop":"8538","bp_start":"139000001","bp_stop":"143100000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"24","iscn_start":"8538","iscn_stop":"9206","bp_start":"143100001","bp_stop":"149200000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"q","band":"25.1","iscn_start":"9206","iscn_stop":"9334","bp_start":"149200001","bp_stop":"152300000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"25.2","iscn_start":"9334","iscn_stop":"9443","bp_start":"152300001","bp_stop":"155300000","stain":"gpos","density":"50"},{"#chromosome":"3","arm":"q","band":"25.3","iscn_start":"9443","iscn_stop":"9639","bp_start":"155300001","bp_stop":"161000000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"26.1","iscn_start":"9639","iscn_stop":"10147","bp_start":"161000001","bp_stop":"167900000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"q","band":"26.2","iscn_start":"10147","iscn_stop":"10269","bp_start":"167900001","bp_stop":"171200000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"26.3","iscn_start":"10269","iscn_stop":"10702","bp_start":"171200001","bp_stop":"183000000","stain":"gpos","density":"100"},{"#chromosome":"3","arm":"q","band":"27","iscn_start":"10702","iscn_stop":"11022","bp_start":"183000001","bp_stop":"188200000","stain":"gneg","density":""},{"#chromosome":"3","arm":"q","band":"28","iscn_start":"11022","iscn_stop":"11352","bp_start":"188200001","bp_stop":"192600000","stain":"gpos","density":"75"},{"#chromosome":"3","arm":"q","band":"29","iscn_start":"11352","iscn_stop":"11700","bp_start":"192600001","bp_stop":"198295559","stain":"gneg","density":""},{"#chromosome":"4","arm":"p","band":"16","iscn_start":"0","iscn_stop":"843","bp_start":"1","bp_stop":"11300000","stain":"gneg","density":""},{"#chromosome":"4","arm":"p","band":"15.3","iscn_start":"843","iscn_stop":"1330","bp_start":"11300001","bp_stop":"21300000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"p","band":"15.2","iscn_start":"1330","iscn_stop":"1449","bp_start":"21300001","bp_stop":"27700000","stain":"gneg","density":""},{"#chromosome":"4","arm":"p","band":"15.1","iscn_start":"1449","iscn_stop":"1857","bp_start":"27700001","bp_stop":"35800000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"p","band":"14","iscn_start":"1857","iscn_stop":"2292","bp_start":"35800001","bp_stop":"41200000","stain":"gneg","density":""},{"#chromosome":"4","arm":"p","band":"13","iscn_start":"2292","iscn_stop":"2647","bp_start":"41200001","bp_stop":"44600000","stain":"gpos","density":"50"},{"#chromosome":"4","arm":"p","band":"12","iscn_start":"2647","iscn_stop":"2950","bp_start":"44600001","bp_stop":"48200000","stain":"gneg","density":""},{"#chromosome":"4","arm":"p","band":"11","iscn_start":"2950","iscn_stop":"3095","bp_start":"48200001","bp_stop":"50000000","stain":"acen","density":""},{"#chromosome":"4","arm":"q","band":"11","iscn_start":"3095","iscn_stop":"3211","bp_start":"50000001","bp_stop":"51800000","stain":"acen","density":""},{"#chromosome":"4","arm":"q","band":"12","iscn_start":"3211","iscn_stop":"3649","bp_start":"51800001","bp_stop":"58500000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"13.1","iscn_start":"3649","iscn_stop":"3993","bp_start":"58500001","bp_stop":"65500000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"13.2","iscn_start":"3993","iscn_stop":"4147",
                "bp_start":"65500001","bp_stop":"69400000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"13.3","iscn_start":"4147","iscn_stop":"4422","bp_start":"69400001","bp_stop":"75300000","stain":"gpos","density":"75"},{"#chromosome":"4","arm":"q","band":"21.1","iscn_start":"4422","iscn_stop":"4577","bp_start":"75300001","bp_stop":"78000000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"21.2","iscn_start":"4577","iscn_stop":"4920","bp_start":"78000001","bp_stop":"86000000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"21.3","iscn_start":"4920","iscn_stop":"5169","bp_start":"86000001","bp_stop":"87100000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"22","iscn_start":"5169","iscn_stop":"5735","bp_start":"87100001","bp_stop":"97900000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"23","iscn_start":"5735","iscn_stop":"5838","bp_start":"97900001","bp_stop":"100100000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"24","iscn_start":"5838","iscn_stop":"6353","bp_start":"100100001","bp_stop":"106700000","stain":"gpos","density":"50"},{"#chromosome":"4","arm":"q","band":"25","iscn_start":"6353","iscn_stop":"6560","bp_start":"106700001","bp_stop":"113200000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"26","iscn_start":"6560","iscn_stop":"7294","bp_start":"113200001","bp_stop":"119900000","stain":"gpos","density":"75"},{"#chromosome":"4","arm":"q","band":"27","iscn_start":"7294","iscn_stop":"7577","bp_start":"119900001","bp_stop":"122800000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"28","iscn_start":"7577","iscn_stop":"8440","bp_start":"122800001","bp_stop":"138500000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"31.1","iscn_start":"8440","iscn_stop":"8813","bp_start":"138500001","bp_stop":"140600000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"31.2","iscn_start":"8813","iscn_stop":"9071","bp_start":"140600001","bp_stop":"150200000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"31.3","iscn_start":"9071","iscn_stop":"9431","bp_start":"150200001","bp_stop":"154600000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"32","iscn_start":"9431","iscn_stop":"10050","bp_start":"154600001","bp_stop":"169200000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"33","iscn_start":"10050","iscn_stop":"10281","bp_start":"169200001","bp_stop":"171000000","stain":"gneg","density":""},{"#chromosome":"4","arm":"q","band":"34","iscn_start":"10281","iscn_stop":"10732","bp_start":"171000001","bp_stop":"182300000","stain":"gpos","density":"100"},{"#chromosome":"4","arm":"q","band":"35","iscn_start":"10732","iscn_stop":"11170","bp_start":"182300001","bp_stop":"190214555","stain":"gneg","density":""},{"#chromosome":"5","arm":"p","band":"15.3","iscn_start":"0","iscn_stop":"574","bp_start":"1","bp_stop":"9900000","stain":"gneg","density":""},{"#chromosome":"5","arm":"p","band":"15.2","iscn_start":"574","iscn_stop":"815","bp_start":"9900001","bp_stop":"15000000","stain":"gpos","density":"50"},{"#chromosome":"5","arm":"p","band":"15.1","iscn_start":"815","iscn_stop":"1028","bp_start":"15000001","bp_stop":"18400000","stain":"gneg","density":""},{"#chromosome":"5","arm":"p","band":"14","iscn_start":"1028","iscn_stop":"1897","bp_start":"18400001","bp_stop":"28900000","stain":"gpos","density":"100"},{"#chromosome":"5","arm":"p","band":"13.3","iscn_start":"1897","iscn_stop":"2123","bp_start":"28900001","bp_stop":"33800000","stain":"gneg","density":""},{"#chromosome":"5","arm":"p","band":"13.2","iscn_start":"2123","iscn_stop":"2301","bp_start":"33800001","bp_stop":"38400000","stain":"gpos","density":"25"},{"#chromosome":"5","arm":"p","band":"13.1","iscn_start":"2301","iscn_stop":"2444","bp_start":"38400001","bp_stop":"42500000","stain":"gneg","density":""},{"#chromosome":"5","arm":"p","band":"12","iscn_start":"2444","iscn_stop":"2712","bp_start":"42500001","bp_stop":"46100000","stain":"gpos","density":"50"},{"#chromosome":"5","arm":"p","band":"11","iscn_start":"2712","iscn_stop":"2845","bp_start":"46100001","bp_stop":"48800000","stain":"acen","density":""},{"#chromosome":"5","arm":"q","band":"11.1","iscn_start":"2845","iscn_stop":"2990","bp_start":"48800001","bp_stop":"51400000","stain":"acen","density":""},{"#chromosome":"5","arm":"q","band":"11.2","iscn_start":"2990","iscn_stop":"3518","bp_start":"51400001","bp_stop":"59600000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"12","iscn_start":"3518","iscn_stop":"3966","bp_start":"59600001","bp_stop":"67400000","stain":"gpos","density":"100"},{"#chromosome":"5","arm":"q","band":"13.1","iscn_start":"3966","iscn_stop":"4249","bp_start":"67400001","bp_stop":"69100000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"13.2","iscn_start":"4249","iscn_stop":"4462","bp_start":"69100001","bp_stop":"74000000","stain":"gpos","density":"50"},{"#chromosome":"5","arm":"q","band":"13.3","iscn_start":"4462","iscn_stop":"4691","bp_start":"74000001","bp_stop":"77600000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"14","iscn_start":"4691","iscn_stop":"5496","bp_start":"77600001","bp_stop":"93000000","stain":"gpos","density":"100"},{"#chromosome":"5","arm":"q","band":"15","iscn_start":"5496","iscn_stop":"5747","bp_start":"93000001","bp_stop":"98900000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"21","iscn_start":"5747","iscn_stop":"6577","bp_start":"98900001","bp_stop":"110200000","stain":"gpos","density":"100"},{"#chromosome":"5","arm":"q","band":"22","iscn_start":"6577","iscn_stop":"6802","bp_start":"110200001","bp_stop":"115900000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"23.1","iscn_start":"6802","iscn_stop":"7130","bp_start":"115900001","bp_stop":"122100000","stain":"gpos","density":"75"},{"#chromosome":"5","arm":"q","band":"23.2","iscn_start":"7130","iscn_stop":"7513","bp_start":"122100001","bp_stop":"127900000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"23.3","iscn_start":"7513","iscn_stop":"7751","bp_start":"127900001","bp_stop":"131200000","stain":"gpos","density":"75"},{"#chromosome":"5","arm":"q","band":"31.1","iscn_start":"7751","iscn_stop":"8146","bp_start":"131200001","bp_stop":"136900000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"31.2","iscn_start":"8146","iscn_stop":"8381","bp_start":"136900001","bp_stop":"140100000","stain":"gpos","density":"25"},{"#chromosome":"5","arm":"q","band":"31.3","iscn_start":"8381","iscn_stop":"8740","bp_start":"140100001","bp_stop":"145100000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"32","iscn_start":"8740","iscn_stop":"9189","bp_start":"145100001","bp_stop":"150400000","stain":"gpos","density":"75"},{"#chromosome":"5","arm":"q","band":"33.1","iscn_start":"9189","iscn_stop":"9331","bp_start":"150400001","bp_stop":"153300000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"33.2","iscn_start":"9331","iscn_stop":"9432","bp_start":"153300001","bp_stop":"156300000","stain":"gpos","density":"50"},{"#chromosome":"5","arm":"q","band":"33.3","iscn_start":"9432","iscn_stop":"9558","bp_start":"156300001","bp_stop":"160500000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"34","iscn_start":"9558","iscn_stop":"10138","bp_start":"160500001","bp_stop":"169000000","stain":"gpos","density":"100"},{"#chromosome":"5","arm":"q","band":"35.1","iscn_start":"10138","iscn_stop":"10271","bp_start":"169000001","bp_stop":"173300000","stain":"gneg","density":""},{"#chromosome":"5","arm":"q","band":"35.2","iscn_start":"10271","iscn_stop":"10388","bp_start":"173300001","bp_stop":"177100000","stain":"gpos","density":"25"},{"#chromosome":"5","arm":"q","band":"35.3","iscn_start":"10388","iscn_stop":"10600","bp_start":"177100001","bp_stop":"181538259","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"25","iscn_start":"0","iscn_stop":"405","bp_start":"1","bp_stop":"7100000","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"24","iscn_start":"405","iscn_stop":"643","bp_start":"7100001","bp_stop":"13400000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"p","band":"23","iscn_start":"643","iscn_stop":"1047","bp_start":"13400001","bp_stop":"15200000","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"22.3","iscn_start":"1047","iscn_stop":"1381","bp_start":"15200001","bp_stop":"25200000","stain":"gpos","density":"75"},{"#chromosome":"6","arm":"p","band":"22.2","iscn_start":"1381","iscn_stop":"1541","bp_start":"25200001","bp_stop":"27100000","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"22.1","iscn_start":"1541","iscn_stop":"1773","bp_start":"27100001","bp_stop":"30500000","stain":"gpos","density":"50"},{"#chromosome":"6","arm":"p","band":"21.3","iscn_start":"1773","iscn_stop":"2428","bp_start":"30500001","bp_stop":"36600000","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"21.2","iscn_start":"2428","iscn_stop":"2571","bp_start":"36600001","bp_stop":"40500000","stain":"gpos","density":"25"},{"#chromosome":"6","arm":"p","band":"21.1","iscn_start":"2571","iscn_stop":"2999","bp_start":"40500001","bp_stop":"46200000","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"12","iscn_start":"2999","iscn_stop":"3618","bp_start":"46200001","bp_stop":"57200000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"p","band":"11.2","iscn_start":"3618","iscn_stop":"3714","bp_start":"57200001","bp_stop":"58500000","stain":"gneg","density":""},{"#chromosome":"6","arm":"p","band":"11.1","iscn_start":"3714","iscn_stop":"3856","bp_start":"58500001","bp_stop":"59800000","stain":"acen","density":""},{"#chromosome":"6","arm":"q","band":"11","iscn_start":"3856","iscn_stop":"3999","bp_start":"59800001","bp_stop":"62700000","stain":"acen",
                "density":""},{"#chromosome":"6","arm":"q","band":"12","iscn_start":"3999","iscn_stop":"4439","bp_start":"62700001","bp_stop":"69200000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"q","band":"13","iscn_start":"4439","iscn_stop":"4701","bp_start":"69200001","bp_stop":"75200000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"14","iscn_start":"4701","iscn_stop":"5129","bp_start":"75200001","bp_stop":"87300000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"q","band":"15","iscn_start":"5129","iscn_stop":"5402","bp_start":"87300001","bp_stop":"92500000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"16.1","iscn_start":"5402","iscn_stop":"5742","bp_start":"92500001","bp_stop":"98900000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"q","band":"16.2","iscn_start":"5742","iscn_stop":"5807","bp_start":"98900001","bp_stop":"100000000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"16.3","iscn_start":"5807","iscn_stop":"6068","bp_start":"100000001","bp_stop":"105000000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"q","band":"21","iscn_start":"6068","iscn_stop":"6746","bp_start":"105000001","bp_stop":"114200000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"22.1","iscn_start":"6746","iscn_stop":"6955","bp_start":"114200001","bp_stop":"117900000","stain":"gpos","density":"75"},{"#chromosome":"6","arm":"q","band":"22.2","iscn_start":"6955","iscn_stop":"7067","bp_start":"117900001","bp_stop":"118100000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"22.3","iscn_start":"7067","iscn_stop":"7793","bp_start":"118100001","bp_stop":"130000000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"q","band":"23.1","iscn_start":"7793","iscn_stop":"7962","bp_start":"130000001","bp_stop":"130900000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"23.2","iscn_start":"7962","iscn_stop":"8096","bp_start":"130900001","bp_stop":"134700000","stain":"gpos","density":"50"},{"#chromosome":"6","arm":"q","band":"23.3","iscn_start":"8096","iscn_stop":"8221","bp_start":"134700001","bp_stop":"138300000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"24","iscn_start":"8221","iscn_stop":"8851","bp_start":"138300001","bp_stop":"148500000","stain":"gpos","density":"100"},{"#chromosome":"6","arm":"q","band":"25.1","iscn_start":"8851","iscn_stop":"8996","bp_start":"148500001","bp_stop":"152100000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"25.2","iscn_start":"8996","iscn_stop":"9119","bp_start":"152100001","bp_stop":"155200000","stain":"gpos","density":"50"},{"#chromosome":"6","arm":"q","band":"25.3","iscn_start":"9119","iscn_stop":"9386","bp_start":"155200001","bp_stop":"160600000","stain":"gneg","density":""},{"#chromosome":"6","arm":"q","band":"26","iscn_start":"9386","iscn_stop":"9696","bp_start":"160600001","bp_stop":"164100000","stain":"gpos","density":"50"},{"#chromosome":"6","arm":"q","band":"27","iscn_start":"9696","iscn_stop":"10100","bp_start":"164100001","bp_stop":"170805979","stain":"gneg","density":""},{"#chromosome":"7","arm":"p","band":"22","iscn_start":"0","iscn_stop":"470","bp_start":"1","bp_stop":"7200000","stain":"gneg","density":""},{"#chromosome":"7","arm":"p","band":"21","iscn_start":"470","iscn_stop":"1300","bp_start":"7200001","bp_stop":"20900000","stain":"gpos","density":"100"},{"#chromosome":"7","arm":"p","band":"15.3","iscn_start":"1300","iscn_stop":"1555","bp_start":"20900001","bp_stop":"25500000","stain":"gneg","density":""},{"#chromosome":"7","arm":"p","band":"15.2","iscn_start":"1555","iscn_stop":"1700","bp_start":"25500001","bp_stop":"27900000","stain":"gpos","density":"50"},{"#chromosome":"7","arm":"p","band":"15.1","iscn_start":"1700","iscn_stop":"1894","bp_start":"27900001","bp_stop":"28800000","stain":"gneg","density":""},{"#chromosome":"7","arm":"p","band":"14","iscn_start":"1894","iscn_stop":"2339","bp_start":"28800001","bp_stop":"43300000","stain":"gpos","density":"100"},{"#chromosome":"7","arm":"p","band":"13","iscn_start":"2339","iscn_stop":"2698","bp_start":"43300001","bp_stop":"45400000","stain":"gneg","density":""},{"#chromosome":"7","arm":"p","band":"12","iscn_start":"2698","iscn_stop":"3069","bp_start":"45400001","bp_stop":"53900000","stain":"gpos","density":"100"},{"#chromosome":"7","arm":"p","band":"11.2","iscn_start":"3069","iscn_stop":"3267","bp_start":"53900001","bp_stop":"58100000","stain":"gneg","density":""},{"#chromosome":"7","arm":"p","band":"11.1","iscn_start":"3267","iscn_stop":"3391","bp_start":"58100001","bp_stop":"60100000","stain":"acen","density":""},{"#chromosome":"7","arm":"q","band":"11.1","iscn_start":"3391","iscn_stop":"3513","bp_start":"60100001","bp_stop":"62100000","stain":"acen","density":""},{"#chromosome":"7","arm":"q","band":"11.21","iscn_start":"3742","iscn_stop":"4012","bp_start":"62100001","bp_stop":"67500000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"11.22","iscn_start":"4012","iscn_stop":"4200","bp_start":"67500001","bp_stop":"72700000","stain":"gpos","density":"50"},{"#chromosome":"7","arm":"q","band":"11.23","iscn_start":"4200","iscn_stop":"4605","bp_start":"72700001","bp_stop":"77900000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"21.1","iscn_start":"4605","iscn_stop":"5263","bp_start":"77900001","bp_stop":"91500000","stain":"gpos","density":"100"},{"#chromosome":"7","arm":"q","band":"21.2","iscn_start":"5263","iscn_stop":"5371","bp_start":"91500001","bp_stop":"93300000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"21.3","iscn_start":"5371","iscn_stop":"5612","bp_start":"93300001","bp_stop":"98400000","stain":"gpos","density":"75"},{"#chromosome":"7","arm":"q","band":"22","iscn_start":"5612","iscn_stop":"6292","bp_start":"98400001","bp_stop":"107800000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"31.1","iscn_start":"6292","iscn_stop":"6584","bp_start":"107800001","bp_stop":"115000000","stain":"gpos","density":"75"},{"#chromosome":"7","arm":"q","band":"31.2","iscn_start":"6584","iscn_stop":"6836","bp_start":"115000001","bp_stop":"117700000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"31.3","iscn_start":"6836","iscn_stop":"7517","bp_start":"117700001","bp_stop":"127500000","stain":"gpos","density":"100"},{"#chromosome":"7","arm":"q","band":"32","iscn_start":"7517","iscn_stop":"7966","bp_start":"127500001","bp_stop":"132900000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"33","iscn_start":"7966","iscn_stop":"8221","bp_start":"132900001","bp_stop":"138500000","stain":"gpos","density":"50"},{"#chromosome":"7","arm":"q","band":"34","iscn_start":"8221","iscn_stop":"8440","bp_start":"138500001","bp_stop":"143400000","stain":"gneg","density":""},{"#chromosome":"7","arm":"q","band":"35","iscn_start":"8440","iscn_stop":"8719","bp_start":"143400001","bp_stop":"148200000","stain":"gpos","density":"75"},{"#chromosome":"7","arm":"q","band":"36","iscn_start":"8719","iscn_stop":"9350","bp_start":"148200001","bp_stop":"159345973","stain":"gneg","density":""},{"#chromosome":"8","arm":"p","band":"23.3","iscn_start":"0","iscn_stop":"102","bp_start":"1","bp_stop":"2300000","stain":"gneg","density":""},{"#chromosome":"8","arm":"p","band":"23.2","iscn_start":"102","iscn_stop":"293","bp_start":"2300001","bp_stop":"6300000","stain":"gpos","density":"75"},{"#chromosome":"8","arm":"p","band":"23.1","iscn_start":"293","iscn_stop":"610","bp_start":"6300001","bp_stop":"12800000","stain":"gneg","density":""},{"#chromosome":"8","arm":"p","band":"22","iscn_start":"610","iscn_stop":"1138","bp_start":"12800001","bp_stop":"19200000","stain":"gpos","density":"100"},{"#chromosome":"8","arm":"p","band":"21.3","iscn_start":"1138","iscn_stop":"1288","bp_start":"19200001","bp_stop":"23500000","stain":"gneg","density":""},{"#chromosome":"8","arm":"p","band":"21.2","iscn_start":"1288","iscn_stop":"1449","bp_start":"23500001","bp_stop":"27500000","stain":"gpos","density":"50"},{"#chromosome":"8","arm":"p","band":"21.1","iscn_start":"1449","iscn_stop":"1656","bp_start":"27500001","bp_stop":"29000000","stain":"gneg","density":""},{"#chromosome":"8","arm":"p","band":"12","iscn_start":"1656","iscn_stop":"2266","bp_start":"29000001","bp_stop":"36700000","stain":"gpos","density":"75"},{"#chromosome":"8","arm":"p","band":"11.2","iscn_start":"2266","iscn_stop":"2611","bp_start":"36700001","bp_stop":"43200000","stain":"gneg","density":""},{"#chromosome":"8","arm":"p","band":"11.1","iscn_start":"2611","iscn_stop":"2743","bp_start":"43200001","bp_stop":"45200000","stain":"acen","density":""},{"#chromosome":"8","arm":"q","band":"11.1","iscn_start":"2743","iscn_stop":"2891","bp_start":"45200001","bp_stop":"47200000","stain":"acen","density":""},{"#chromosome":"8","arm":"q","band":"11.21","iscn_start":"2946","iscn_stop":"3001","bp_start":"47200001","bp_stop":"51300000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"11.22","iscn_start":"3001","iscn_stop":"3070","bp_start":"51300001","bp_stop":"51700000","stain":"gpos","density":"75"},{"#chromosome":"8","arm":"q","band":"11.23","iscn_start":"3070","iscn_stop":"3148","bp_start":"51700001","bp_stop":"54600000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"12","iscn_start":"3148","iscn_stop":"3641","bp_start":"54600001","bp_stop":"65100000","stain":"gpos","density":"100"},{"#chromosome":"8","arm":"q","band":"13","iscn_start":"3641","iscn_stop":"4145","bp_start":"65100001","bp_stop":"72000000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"21.1","iscn_start":"4145","iscn_stop":"4756","bp_start":"72000001","bp_stop":"83500000","stain":"gpos","density":"100"},{"#chromosome":"8","arm":"q","band":"21.2","iscn_start":"4756","iscn_stop":"4885","bp_start":"83500001","bp_stop":"85900000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"21.3","iscn_start":"4885",
                "iscn_stop":"5506","bp_start":"85900001","bp_stop":"92300000","stain":"gpos","density":"100"},{"#chromosome":"8","arm":"q","band":"22.1","iscn_start":"5506","iscn_stop":"5856","bp_start":"92300001","bp_stop":"97900000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"22.2","iscn_start":"5856","iscn_stop":"5996","bp_start":"97900001","bp_stop":"100500000","stain":"gpos","density":"25"},{"#chromosome":"8","arm":"q","band":"22.3","iscn_start":"5996","iscn_stop":"6276","bp_start":"100500001","bp_stop":"105100000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"23","iscn_start":"6276","iscn_stop":"7026","bp_start":"105100001","bp_stop":"116700000","stain":"gpos","density":"100"},{"#chromosome":"8","arm":"q","band":"24.1","iscn_start":"7026","iscn_stop":"7628","bp_start":"116700001","bp_stop":"126300000","stain":"gneg","density":""},{"#chromosome":"8","arm":"q","band":"24.2","iscn_start":"7628","iscn_stop":"7885","bp_start":"126300001","bp_stop":"138900000","stain":"gpos","density":"100"},{"#chromosome":"8","arm":"q","band":"24.3","iscn_start":"7885","iscn_stop":"8250","bp_start":"138900001","bp_stop":"145138636","stain":"gneg","density":""},{"#chromosome":"9","arm":"p","band":"24","iscn_start":"0","iscn_stop":"355","bp_start":"1","bp_stop":"9000000","stain":"gneg","density":""},{"#chromosome":"9","arm":"p","band":"23","iscn_start":"355","iscn_stop":"841","bp_start":"9000001","bp_stop":"14200000","stain":"gpos","density":"75"},{"#chromosome":"9","arm":"p","band":"22","iscn_start":"841","iscn_stop":"942","bp_start":"14200001","bp_stop":"19900000","stain":"gneg","density":""},{"#chromosome":"9","arm":"p","band":"21","iscn_start":"942","iscn_stop":"1661","bp_start":"19900001","bp_stop":"33200000","stain":"gpos","density":"100"},{"#chromosome":"9","arm":"p","band":"13","iscn_start":"1661","iscn_stop":"2148","bp_start":"33200001","bp_stop":"39000000","stain":"gneg","density":""},{"#chromosome":"9","arm":"p","band":"12","iscn_start":"2148","iscn_stop":"2522","bp_start":"39000001","bp_stop":"40000000","stain":"gpos","density":"50"},{"#chromosome":"9","arm":"p","band":"11","iscn_start":"2522","iscn_stop":"2634","bp_start":"40000001","bp_stop":"43000000","stain":"acen","density":""},{"#chromosome":"9","arm":"q","band":"11","iscn_start":"2634","iscn_stop":"2773","bp_start":"43000001","bp_stop":"45500000","stain":"acen","density":""},{"#chromosome":"9","arm":"q","band":"12","iscn_start":"2773","iscn_stop":"3695","bp_start":"45500001","bp_stop":"61500000","stain":"gvar","density":""},{"#chromosome":"9","arm":"q","band":"13","iscn_start":"3695","iscn_stop":"3864","bp_start":"61500001","bp_stop":"65000000","stain":"gneg","density":""},{"#chromosome":"9","arm":"q","band":"21.1","iscn_start":"3864","iscn_stop":"4232","bp_start":"65000001","bp_stop":"76600000","stain":"gpos","density":"100"},{"#chromosome":"9","arm":"q","band":"21.2","iscn_start":"4232","iscn_stop":"4370","bp_start":"76600001","bp_stop":"78500000","stain":"gneg","density":""},{"#chromosome":"9","arm":"q","band":"21.3","iscn_start":"4370","iscn_stop":"4865","bp_start":"78500001","bp_stop":"87800000","stain":"gpos","density":"100"},{"#chromosome":"9","arm":"q","band":"22.1","iscn_start":"4865","iscn_stop":"5160","bp_start":"87800001","bp_stop":"89200000","stain":"gneg","density":""},{"#chromosome":"9","arm":"q","band":"22.2","iscn_start":"5160","iscn_stop":"5283","bp_start":"89200001","bp_stop":"91200000","stain":"gpos","density":"25"},{"#chromosome":"9","arm":"q","band":"22.3","iscn_start":"5283","iscn_stop":"5857","bp_start":"91200001","bp_stop":"99800000","stain":"gneg","density":""},{"#chromosome":"9","arm":"q","band":"31","iscn_start":"5857","iscn_stop":"6373","bp_start":"99800001","bp_stop":"112100000","stain":"gpos","density":"100"},{"#chromosome":"9","arm":"q","band":"32","iscn_start":"6373","iscn_stop":"6571","bp_start":"112100001","bp_stop":"114900000","stain":"gneg","density":""},{"#chromosome":"9","arm":"q","band":"33","iscn_start":"6571","iscn_stop":"6958","bp_start":"114900001","bp_stop":"127500000","stain":"gpos","density":"100"},{"#chromosome":"9","arm":"q","band":"34.1","iscn_start":"6958","iscn_stop":"7448","bp_start":"127500001","bp_stop":"133100000","stain":"gneg","density":""},{"#chromosome":"9","arm":"q","band":"34.2","iscn_start":"7448","iscn_stop":"7559","bp_start":"133100001","bp_stop":"134500000","stain":"gpos","density":"25"},{"#chromosome":"9","arm":"q","band":"34.3","iscn_start":"7559","iscn_stop":"7950","bp_start":"134500001","bp_stop":"138394717","stain":"gneg","density":""},{"#chromosome":"10","arm":"p","band":"15","iscn_start":"0","iscn_stop":"456","bp_start":"1","bp_stop":"6600000","stain":"gneg","density":""},{"#chromosome":"10","arm":"p","band":"14","iscn_start":"456","iscn_stop":"841","bp_start":"6600001","bp_stop":"12200000","stain":"gpos","density":"75"},{"#chromosome":"10","arm":"p","band":"13","iscn_start":"841","iscn_stop":"1261","bp_start":"12200001","bp_stop":"17300000","stain":"gneg","density":""},{"#chromosome":"10","arm":"p","band":"12.3","iscn_start":"1261","iscn_stop":"1658","bp_start":"17300001","bp_stop":"22300000","stain":"gpos","density":"100"},{"#chromosome":"10","arm":"p","band":"12.2","iscn_start":"1658","iscn_stop":"1711","bp_start":"22300001","bp_stop":"24300000","stain":"gneg","density":""},{"#chromosome":"10","arm":"p","band":"12.1","iscn_start":"1711","iscn_stop":"1923","bp_start":"24300001","bp_stop":"29300000","stain":"gpos","density":"50"},{"#chromosome":"10","arm":"p","band":"11.2","iscn_start":"1923","iscn_stop":"2477","bp_start":"29300001","bp_stop":"38000000","stain":"gneg","density":""},{"#chromosome":"10","arm":"p","band":"11.1","iscn_start":"2477","iscn_stop":"2620","bp_start":"38000001","bp_stop":"39800000","stain":"acen","density":""},{"#chromosome":"10","arm":"q","band":"11.1","iscn_start":"2620","iscn_stop":"2762","bp_start":"39800001","bp_stop":"41600000","stain":"acen","density":""},{"#chromosome":"10","arm":"q","band":"11.2","iscn_start":"2762","iscn_stop":"3222","bp_start":"41600001","bp_stop":"51100000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"21.1","iscn_start":"3222","iscn_stop":"3832","bp_start":"51100001","bp_stop":"59400000","stain":"gpos","density":"100"},{"#chromosome":"10","arm":"q","band":"21.2","iscn_start":"3832","iscn_stop":"3985","bp_start":"59400001","bp_stop":"62800000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"21.3","iscn_start":"3985","iscn_stop":"4442","bp_start":"62800001","bp_stop":"68800000","stain":"gpos","density":"100"},{"#chromosome":"10","arm":"q","band":"22.1","iscn_start":"4442","iscn_stop":"4845","bp_start":"68800001","bp_stop":"73100000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"22.2","iscn_start":"4845","iscn_stop":"5047","bp_start":"73100001","bp_stop":"75900000","stain":"gpos","density":"50"},{"#chromosome":"10","arm":"q","band":"22.3","iscn_start":"5047","iscn_stop":"5388","bp_start":"75900001","bp_stop":"80300000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"23.1","iscn_start":"5388","iscn_stop":"5566","bp_start":"80300001","bp_stop":"86100000","stain":"gpos","density":"100"},{"#chromosome":"10","arm":"q","band":"23.2","iscn_start":"5566","iscn_stop":"5667","bp_start":"86100001","bp_stop":"87700000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"23.3","iscn_start":"5667","iscn_stop":"6096","bp_start":"87700001","bp_stop":"95300000","stain":"gpos","density":"100"},{"#chromosome":"10","arm":"q","band":"24.1","iscn_start":"6096","iscn_stop":"6193","bp_start":"95300001","bp_stop":"97500000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"24.2","iscn_start":"6193","iscn_stop":"6371","bp_start":"97500001","bp_stop":"100100000","stain":"gpos","density":"50"},{"#chromosome":"10","arm":"q","band":"24.3","iscn_start":"6371","iscn_stop":"6644","bp_start":"100100001","bp_stop":"104000000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"25.1","iscn_start":"6644","iscn_stop":"7017","bp_start":"104000001","bp_stop":"110100000","stain":"gpos","density":"100"},{"#chromosome":"10","arm":"q","band":"25.2","iscn_start":"7017","iscn_stop":"7174","bp_start":"110100001","bp_stop":"113100000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"25.3","iscn_start":"7174","iscn_stop":"7351","bp_start":"113100001","bp_stop":"117300000","stain":"gpos","density":"75"},{"#chromosome":"10","arm":"q","band":"26.1","iscn_start":"7351","iscn_stop":"7722","bp_start":"117300001","bp_stop":"125700000","stain":"gneg","density":""},{"#chromosome":"10","arm":"q","band":"26.2","iscn_start":"7722","iscn_stop":"7852","bp_start":"125700001","bp_stop":"128800000","stain":"gpos","density":"50"},{"#chromosome":"10","arm":"q","band":"26.3","iscn_start":"7852","iscn_stop":"8050","bp_start":"128800001","bp_stop":"133797422","stain":"gneg","density":""},{"#chromosome":"11","arm":"p","band":"15.5","iscn_start":"0","iscn_stop":"212","bp_start":"1","bp_stop":"2800000","stain":"gneg","density":""},{"#chromosome":"11","arm":"p","band":"15.4","iscn_start":"212","iscn_stop":"425","bp_start":"2800001","bp_stop":"11700000","stain":"gpos","density":"50"},{"#chromosome":"11","arm":"p","band":"15.3","iscn_start":"425","iscn_stop":"687","bp_start":"11700001","bp_stop":"13800000","stain":"gneg","density":""},{"#chromosome":"11","arm":"p","band":"15.2","iscn_start":"687","iscn_stop":"862","bp_start":"13800001","bp_stop":"16900000","stain":"gpos","density":"50"},{"#chromosome":"11","arm":"p","band":"15.1","iscn_start":"862","iscn_stop":"1149","bp_start":"16900001","bp_stop":"22000000","stain":"gneg","density":""},{"#chromosome":"11","arm":"p","band":"14","iscn_start":"1149","iscn_stop":"1957","bp_start":"22000001","bp_stop":"31000000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"p","band":"13","iscn_start":"1957","iscn_stop":"2262","bp_start":"31000001","bp_stop":"36400000",
                "stain":"gneg","density":""},{"#chromosome":"11","arm":"p","band":"12","iscn_start":"2262","iscn_stop":"2639","bp_start":"36400001","bp_stop":"43400000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"p","band":"11.2","iscn_start":"2639","iscn_stop":"3186","bp_start":"43400001","bp_stop":"48800000","stain":"gneg","density":""},{"#chromosome":"11","arm":"p","band":"11.12","iscn_start":"3257","iscn_stop":"3309","bp_start":"48800001","bp_stop":"51000000","stain":"gpos","density":"75"},{"#chromosome":"11","arm":"p","band":"11.11","iscn_start":"2872","iscn_stop":"3035","bp_start":"51000001","bp_stop":"53400000","stain":"acen"},{"#chromosome":"11","arm":"q","band":"11","iscn_start":"3348","iscn_stop":"3478","bp_start":"53400001","bp_stop":"55800000","stain":"acen","density":""},{"#chromosome":"11","arm":"q","band":"12","iscn_start":"3478","iscn_stop":"3999","bp_start":"55800001","bp_stop":"63600000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"q","band":"13.1","iscn_start":"3999","iscn_stop":"4248","bp_start":"63600001","bp_stop":"66100000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"13.2","iscn_start":"4248","iscn_stop":"4353","bp_start":"66100001","bp_stop":"68700000","stain":"gpos","density":"25"},{"#chromosome":"11","arm":"q","band":"13.3","iscn_start":"4353","iscn_stop":"4584","bp_start":"68700001","bp_stop":"70500000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"13.4","iscn_start":"4584","iscn_stop":"4708","bp_start":"70500001","bp_stop":"75500000","stain":"gpos","density":"50"},{"#chromosome":"11","arm":"q","band":"13.5","iscn_start":"4708","iscn_stop":"4842","bp_start":"75500001","bp_stop":"77400000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"14.1","iscn_start":"4842","iscn_stop":"5186","bp_start":"77400001","bp_stop":"85900000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"q","band":"14.2","iscn_start":"5186","iscn_stop":"5324","bp_start":"85900001","bp_stop":"88600000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"14.3","iscn_start":"5324","iscn_stop":"5599","bp_start":"88600001","bp_stop":"93000000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"q","band":"21","iscn_start":"5599","iscn_stop":"5703","bp_start":"93000001","bp_stop":"97400000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"22.1","iscn_start":"5703","iscn_stop":"5967","bp_start":"97400001","bp_stop":"102300000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"q","band":"22.2","iscn_start":"5967","iscn_stop":"6114","bp_start":"102300001","bp_stop":"103000000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"22.3","iscn_start":"6114","iscn_stop":"6363","bp_start":"103000001","bp_stop":"110600000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"q","band":"23.1","iscn_start":"6363","iscn_stop":"6585","bp_start":"110600001","bp_stop":"112700000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"23.2","iscn_start":"6585","iscn_stop":"6793","bp_start":"112700001","bp_stop":"114600000","stain":"gpos","density":"50"},{"#chromosome":"11","arm":"q","band":"23.3","iscn_start":"6793","iscn_stop":"7311","bp_start":"114600001","bp_stop":"121300000","stain":"gneg","density":""},{"#chromosome":"11","arm":"q","band":"24","iscn_start":"7311","iscn_stop":"7658","bp_start":"121300001","bp_stop":"130900000","stain":"gpos","density":"100"},{"#chromosome":"11","arm":"q","band":"25","iscn_start":"7658","iscn_stop":"7980","bp_start":"130900001","bp_stop":"135086622","stain":"gneg","density":""},{"#chromosome":"12","arm":"p","band":"13.3","iscn_start":"0","iscn_stop":"522","bp_start":"1","bp_stop":"10000000","stain":"gneg","density":""},{"#chromosome":"12","arm":"p","band":"13.2","iscn_start":"522","iscn_stop":"665","bp_start":"10000001","bp_stop":"12600000","stain":"gpos","density":"75"},{"#chromosome":"12","arm":"p","band":"13.1","iscn_start":"665","iscn_stop":"760","bp_start":"12600001","bp_stop":"14600000","stain":"gneg","density":""},{"#chromosome":"12","arm":"p","band":"12.3","iscn_start":"760","iscn_stop":"1092","bp_start":"14600001","bp_stop":"19800000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"p","band":"12.2","iscn_start":"1092","iscn_stop":"1160","bp_start":"19800001","bp_stop":"21100000","stain":"gneg","density":""},{"#chromosome":"12","arm":"p","band":"12.1","iscn_start":"1160","iscn_stop":"1492","bp_start":"21100001","bp_stop":"26300000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"p","band":"11.2","iscn_start":"1492","iscn_stop":"1922","bp_start":"26300001","bp_stop":"33200000","stain":"gneg","density":""},{"#chromosome":"12","arm":"p","band":"11.1","iscn_start":"1922","iscn_stop":"2105","bp_start":"33200001","bp_stop":"35500000","stain":"acen","density":""},{"#chromosome":"12","arm":"q","band":"11","iscn_start":"2105","iscn_stop":"2202","bp_start":"35500001","bp_stop":"37800000","stain":"acen","density":""},{"#chromosome":"12","arm":"q","band":"12","iscn_start":"2202","iscn_stop":"2697","bp_start":"37800001","bp_stop":"46000000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"q","band":"13.1","iscn_start":"2697","iscn_stop":"3204","bp_start":"46000001","bp_stop":"54500000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"13.2","iscn_start":"3204","iscn_stop":"3340","bp_start":"54500001","bp_stop":"56200000","stain":"gpos","density":"25"},{"#chromosome":"12","arm":"q","band":"13.3","iscn_start":"3340","iscn_stop":"3430","bp_start":"56200001","bp_stop":"57700000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"14","iscn_start":"3430","iscn_stop":"3942","bp_start":"57700001","bp_stop":"67300000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"q","band":"15","iscn_start":"3942","iscn_stop":"4286","bp_start":"67300001","bp_stop":"71100000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"21.1","iscn_start":"4286","iscn_stop":"4460","bp_start":"71100001","bp_stop":"75300000","stain":"gpos","density":"75"},{"#chromosome":"12","arm":"q","band":"21.2","iscn_start":"4460","iscn_stop":"4664","bp_start":"75300001","bp_stop":"79900000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"21.3","iscn_start":"4664","iscn_stop":"5293","bp_start":"79900001","bp_stop":"92200000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"q","band":"22","iscn_start":"5293","iscn_stop":"5716","bp_start":"92200001","bp_stop":"95800000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"23","iscn_start":"5716","iscn_stop":"6229","bp_start":"95800001","bp_stop":"108600000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"q","band":"24.1","iscn_start":"6229","iscn_stop":"6714","bp_start":"108600001","bp_stop":"113900000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"24.2","iscn_start":"6714","iscn_stop":"6917","bp_start":"113900001","bp_stop":"120300000","stain":"gpos","density":"100"},{"#chromosome":"12","arm":"q","band":"24.31","iscn_start":"7227","iscn_stop":"7351","bp_start":"120300001","bp_stop":"125400000","stain":"gneg","density":""},{"#chromosome":"12","arm":"q","band":"24.32","iscn_start":"7351","iscn_stop":"7412","bp_start":"125400001","bp_stop":"128700000","stain":"gpos","density":"50"},{"#chromosome":"12","arm":"q","band":"24.33","iscn_start":"7412","iscn_stop":"7500","bp_start":"128700001","bp_stop":"133275309","stain":"gneg","density":""},{"#chromosome":"13","arm":"p","band":"13","iscn_start":"0","iscn_stop":"301","bp_start":"1","bp_stop":"4600000","stain":"gvar","density":""},{"#chromosome":"13","arm":"p","band":"12","iscn_start":"301","iscn_stop":"602","bp_start":"4600001","bp_stop":"10100000","stain":"stalk","density":""},{"#chromosome":"13","arm":"p","band":"11.2","iscn_start":"602","iscn_stop":"1004","bp_start":"10100001","bp_stop":"16500000","stain":"gvar","density":""},{"#chromosome":"13","arm":"p","band":"11.1","iscn_start":"1004","iscn_stop":"1204","bp_start":"16500001","bp_stop":"17700000","stain":"acen","density":""},{"#chromosome":"13","arm":"q","band":"11","iscn_start":"1204","iscn_stop":"1351","bp_start":"17700001","bp_stop":"18900000","stain":"acen","density":""},{"#chromosome":"13","arm":"q","band":"12.1","iscn_start":"1351","iscn_stop":"1735","bp_start":"18900001","bp_stop":"27200000","stain":"gneg","density":""},{"#chromosome":"13","arm":"q","band":"12.2","iscn_start":"1735","iscn_stop":"1821","bp_start":"27200001","bp_stop":"28300000","stain":"gpos","density":"25"},{"#chromosome":"13","arm":"q","band":"12.3","iscn_start":"1821","iscn_stop":"2019","bp_start":"28300001","bp_stop":"31600000","stain":"gneg","density":""},{"#chromosome":"13","arm":"q","band":"13","iscn_start":"2019","iscn_stop":"2295","bp_start":"31600001","bp_stop":"39500000","stain":"gpos","density":"100"},{"#chromosome":"13","arm":"q","band":"14.1","iscn_start":"2295","iscn_stop":"2690","bp_start":"39500001","bp_stop":"46700000","stain":"gneg","density":""},{"#chromosome":"13","arm":"q","band":"14.2","iscn_start":"2690","iscn_stop":"2841","bp_start":"46700001","bp_stop":"50300000","stain":"gpos","density":"50"},{"#chromosome":"13","arm":"q","band":"14.3","iscn_start":"2841","iscn_stop":"3028","bp_start":"50300001","bp_stop":"54700000","stain":"gneg","density":""},{"#chromosome":"13","arm":"q","band":"21.1","iscn_start":"3028","iscn_stop":"3294","bp_start":"54700001","bp_stop":"59000000","stain":"gpos","density":"100"},{"#chromosome":"13","arm":"q","band":"21.2","iscn_start":"3294","iscn_stop":"3444","bp_start":"59000001","bp_stop":"61800000","stain":"gneg","density":""},{"#chromosome":"13","arm":"q","band":"21.3","iscn_start":"3444","iscn_stop":"4093","bp_start":"61800001","bp_stop":"72800000","stain":"gpos","density":"100"},{"#chromosome":"13","arm":"q","band":"22","iscn_start":"4093","iscn_stop":"4625","bp_start":"72800001","bp_stop":"78500000","stain":"gneg","density":""},
                {"#chromosome":"13","arm":"q","band":"31","iscn_start":"4625","iscn_stop":"5280","bp_start":"78500001","bp_stop":"94400000","stain":"gpos","density":"100"},{"#chromosome":"13","arm":"q","band":"32","iscn_start":"5280","iscn_stop":"5742","bp_start":"94400001","bp_stop":"101100000","stain":"gneg","density":""},{"#chromosome":"13","arm":"q","band":"33","iscn_start":"5742","iscn_stop":"6065","bp_start":"101100001","bp_stop":"109600000","stain":"gpos","density":"100"},{"#chromosome":"13","arm":"q","band":"34","iscn_start":"6065","iscn_stop":"6510","bp_start":"109600001","bp_stop":"114364328","stain":"gneg","density":""},{"#chromosome":"14","arm":"p","band":"13","iscn_start":"0","iscn_stop":"305","bp_start":"1","bp_stop":"3600000","stain":"gvar","density":""},{"#chromosome":"14","arm":"p","band":"12","iscn_start":"305","iscn_stop":"611","bp_start":"3600001","bp_stop":"8000000","stain":"stalk","density":""},{"#chromosome":"14","arm":"p","band":"11.2","iscn_start":"611","iscn_stop":"1018","bp_start":"8000001","bp_stop":"16100000","stain":"gvar","density":""},{"#chromosome":"14","arm":"p","band":"11.1","iscn_start":"1018","iscn_stop":"1222","bp_start":"16100001","bp_stop":"17200000","stain":"acen","density":""},{"#chromosome":"14","arm":"q","band":"11.1","iscn_start":"1222","iscn_stop":"1421","bp_start":"17200001","bp_stop":"18200000","stain":"acen","density":""},{"#chromosome":"14","arm":"q","band":"11.2","iscn_start":"1421","iscn_stop":"1808","bp_start":"18200001","bp_stop":"24100000","stain":"gneg","density":""},{"#chromosome":"14","arm":"q","band":"12","iscn_start":"1808","iscn_stop":"2243","bp_start":"24100001","bp_stop":"32900000","stain":"gpos","density":"100"},{"#chromosome":"14","arm":"q","band":"13","iscn_start":"2243","iscn_stop":"2633","bp_start":"32900001","bp_stop":"37400000","stain":"gneg","density":""},{"#chromosome":"14","arm":"q","band":"21","iscn_start":"2633","iscn_stop":"3441","bp_start":"37400001","bp_stop":"50400000","stain":"gpos","density":"100"},{"#chromosome":"14","arm":"q","band":"22","iscn_start":"3441","iscn_stop":"3894","bp_start":"50400001","bp_stop":"57600000","stain":"gneg","density":""},{"#chromosome":"14","arm":"q","band":"23","iscn_start":"3894","iscn_stop":"4169","bp_start":"57600001","bp_stop":"67400000","stain":"gpos","density":"100"},{"#chromosome":"14","arm":"q","band":"24.1","iscn_start":"4169","iscn_stop":"4587","bp_start":"67400001","bp_stop":"69800000","stain":"gneg","density":""},{"#chromosome":"14","arm":"q","band":"24.2","iscn_start":"4587","iscn_stop":"4786","bp_start":"69800001","bp_stop":"73300000","stain":"gpos","density":"50"},{"#chromosome":"14","arm":"q","band":"24.3","iscn_start":"4786","iscn_stop":"5084","bp_start":"73300001","bp_stop":"78800000","stain":"gneg","density":""},{"#chromosome":"14","arm":"q","band":"31","iscn_start":"5084","iscn_stop":"5599","bp_start":"78800001","bp_stop":"89300000","stain":"gpos","density":"100"},{"#chromosome":"14","arm":"q","band":"32.1","iscn_start":"5599","iscn_stop":"5798","bp_start":"89300001","bp_stop":"95800000","stain":"gneg","density":""},{"#chromosome":"14","arm":"q","band":"32.2","iscn_start":"5798","iscn_stop":"5881","bp_start":"95800001","bp_stop":"100900000","stain":"gpos","density":"50"},{"#chromosome":"14","arm":"q","band":"32.3","iscn_start":"5881","iscn_stop":"6300","bp_start":"100900001","bp_stop":"107043718","stain":"gneg","density":""},{"#chromosome":"15","arm":"p","band":"13","iscn_start":"0","iscn_stop":"319","bp_start":"1","bp_stop":"4200000","stain":"gvar","density":""},{"#chromosome":"15","arm":"p","band":"12","iscn_start":"319","iscn_stop":"637","bp_start":"4200001","bp_stop":"9700000","stain":"stalk","density":""},{"#chromosome":"15","arm":"p","band":"11.2","iscn_start":"637","iscn_stop":"1062","bp_start":"9700001","bp_stop":"17500000","stain":"gvar","density":""},{"#chromosome":"15","arm":"p","band":"11.1","iscn_start":"1062","iscn_stop":"1275","bp_start":"17500001","bp_stop":"19000000","stain":"acen","density":""},{"#chromosome":"15","arm":"q","band":"11.1","iscn_start":"1275","iscn_stop":"1498","bp_start":"19000001","bp_stop":"20500000","stain":"acen","density":""},{"#chromosome":"15","arm":"q","band":"11.2","iscn_start":"1498","iscn_stop":"1657","bp_start":"20500001","bp_stop":"25500000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"12","iscn_start":"1657","iscn_stop":"1862","bp_start":"25500001","bp_stop":"27800000","stain":"gpos","density":"50"},{"#chromosome":"15","arm":"q","band":"13","iscn_start":"1862","iscn_stop":"2043","bp_start":"27800001","bp_stop":"33400000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"14","iscn_start":"2043","iscn_stop":"2411","bp_start":"33400001","bp_stop":"39800000","stain":"gpos","density":"75"},{"#chromosome":"15","arm":"q","band":"15","iscn_start":"2411","iscn_stop":"2832","bp_start":"39800001","bp_stop":"44500000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"21.1","iscn_start":"2832","iscn_stop":"3121","bp_start":"44500001","bp_stop":"49200000","stain":"gpos","density":"75"},{"#chromosome":"15","arm":"q","band":"21.2","iscn_start":"3121","iscn_stop":"3312","bp_start":"49200001","bp_stop":"52600000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"21.3","iscn_start":"3312","iscn_stop":"3600","bp_start":"52600001","bp_stop":"58800000","stain":"gpos","density":"75"},{"#chromosome":"15","arm":"q","band":"22.1","iscn_start":"3600","iscn_stop":"3724","bp_start":"58800001","bp_stop":"59000000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"22.2","iscn_start":"3724","iscn_stop":"3820","bp_start":"59000001","bp_stop":"63400000","stain":"gpos","density":"25"},{"#chromosome":"15","arm":"q","band":"22.3","iscn_start":"3820","iscn_stop":"4203","bp_start":"63400001","bp_stop":"67200000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"23","iscn_start":"4203","iscn_stop":"4471","bp_start":"67200001","bp_stop":"72400000","stain":"gpos","density":"25"},{"#chromosome":"15","arm":"q","band":"24","iscn_start":"4471","iscn_stop":"5002","bp_start":"72400001","bp_stop":"78000000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"25","iscn_start":"5002","iscn_stop":"5381","bp_start":"78000001","bp_stop":"88500000","stain":"gpos","density":"100"},{"#chromosome":"15","arm":"q","band":"26.1","iscn_start":"5381","iscn_stop":"5650","bp_start":"88500001","bp_stop":"93800000","stain":"gneg","density":""},{"#chromosome":"15","arm":"q","band":"26.2","iscn_start":"5650","iscn_stop":"5861","bp_start":"93800001","bp_stop":"98000000","stain":"gpos","density":"50"},{"#chromosome":"15","arm":"q","band":"26.3","iscn_start":"5861","iscn_stop":"6070","bp_start":"98000001","bp_stop":"101991189","stain":"gneg","density":""},{"#chromosome":"16","arm":"p","band":"13.3","iscn_start":"0","iscn_stop":"410","bp_start":"1","bp_stop":"7800000","stain":"gneg","density":""},{"#chromosome":"16","arm":"p","band":"13.2","iscn_start":"410","iscn_stop":"638","bp_start":"7800001","bp_stop":"10400000","stain":"gpos","density":"50"},{"#chromosome":"16","arm":"p","band":"13.1","iscn_start":"638","iscn_stop":"1012","bp_start":"10400001","bp_stop":"16700000","stain":"gneg","density":""},{"#chromosome":"16","arm":"p","band":"12","iscn_start":"1012","iscn_stop":"1522","bp_start":"16700001","bp_stop":"28500000","stain":"gpos","density":"100"},{"#chromosome":"16","arm":"p","band":"11.2","iscn_start":"1522","iscn_stop":"2060","bp_start":"28500001","bp_stop":"35300000","stain":"gneg","density":""},{"#chromosome":"16","arm":"p","band":"11.1","iscn_start":"2060","iscn_stop":"2187","bp_start":"35300001","bp_stop":"36800000","stain":"acen","density":""},{"#chromosome":"16","arm":"q","band":"11.1","iscn_start":"2187","iscn_stop":"2310","bp_start":"36800001","bp_stop":"38400000","stain":"acen","density":""},{"#chromosome":"16","arm":"q","band":"11.2","iscn_start":"2310","iscn_stop":"2811","bp_start":"38400001","bp_stop":"47000000","stain":"gvar","density":""},{"#chromosome":"16","arm":"q","band":"12.1","iscn_start":"2811","iscn_stop":"2969","bp_start":"47000001","bp_stop":"52600000","stain":"gneg","density":""},{"#chromosome":"16","arm":"q","band":"12.2","iscn_start":"2969","iscn_stop":"3118","bp_start":"52600001","bp_stop":"56000000","stain":"gpos","density":"50"},{"#chromosome":"16","arm":"q","band":"13","iscn_start":"3118","iscn_stop":"3364","bp_start":"56000001","bp_stop":"57300000","stain":"gneg","density":""},{"#chromosome":"16","arm":"q","band":"21","iscn_start":"3364","iscn_stop":"3750","bp_start":"57300001","bp_stop":"66600000","stain":"gpos","density":"100"},{"#chromosome":"16","arm":"q","band":"22","iscn_start":"3750","iscn_stop":"4189","bp_start":"66600001","bp_stop":"74100000","stain":"gneg","density":""},{"#chromosome":"16","arm":"q","band":"23","iscn_start":"4189","iscn_stop":"4655","bp_start":"74100001","bp_stop":"84100000","stain":"gpos","density":"100"},{"#chromosome":"16","arm":"q","band":"24","iscn_start":"4655","iscn_stop":"5120","bp_start":"84100001","bp_stop":"90338345","stain":"gneg","density":""},{"#chromosome":"17","arm":"p","band":"13","iscn_start":"0","iscn_stop":"602","bp_start":"1","bp_stop":"10800000","stain":"gneg","density":""},{"#chromosome":"17","arm":"p","band":"12","iscn_start":"602","iscn_stop":"1051","bp_start":"10800001","bp_stop":"16100000","stain":"gpos","density":"75"},{"#chromosome":"17","arm":"p","band":"11.2","iscn_start":"1051","iscn_stop":"1492","bp_start":"16100001","bp_stop":"22700000","stain":"gneg","density":""},{"#chromosome":"17","arm":"p","band":"11.1","iscn_start":"1492","iscn_stop":"1618","bp_start":"22700001","bp_stop":"25100000","stain":"acen","density":""},{"#chromosome":"17","arm":"q","band":"11.1","iscn_start":"1618","iscn_stop":"1762","bp_start":"25100001","bp_stop":"27400000","stain":"acen","density":""},{"#chromosome":"17","arm":"q","band":"11.2","iscn_start":"1762","iscn_stop":"2050","bp_start":"27400001",
                "bp_stop":"33500000","stain":"gneg","density":""},{"#chromosome":"17","arm":"q","band":"12","iscn_start":"2050","iscn_stop":"2383","bp_start":"33500001","bp_stop":"39800000","stain":"gpos","density":"50"},{"#chromosome":"17","arm":"q","band":"21.1","iscn_start":"2383","iscn_stop":"2554","bp_start":"39800001","bp_stop":"40200000","stain":"gneg","density":""},{"#chromosome":"17","arm":"q","band":"21.2","iscn_start":"2554","iscn_stop":"2669","bp_start":"40200001","bp_stop":"42800000","stain":"gpos","density":"25"},{"#chromosome":"17","arm":"q","band":"21.3","iscn_start":"2669","iscn_stop":"3149","bp_start":"42800001","bp_stop":"52100000","stain":"gneg","density":""},{"#chromosome":"17","arm":"q","band":"22","iscn_start":"3149","iscn_stop":"3680","bp_start":"52100001","bp_stop":"59500000","stain":"gpos","density":"75"},{"#chromosome":"17","arm":"q","band":"23","iscn_start":"3680","iscn_stop":"3950","bp_start":"59500001","bp_stop":"64600000","stain":"gneg","density":""},{"#chromosome":"17","arm":"q","band":"24","iscn_start":"3950","iscn_stop":"4441","bp_start":"64600001","bp_stop":"72900000","stain":"gpos","density":"100"},{"#chromosome":"17","arm":"q","band":"25","iscn_start":"4441","iscn_stop":"4950","bp_start":"72900001","bp_stop":"83257441","stain":"gneg","density":""},{"#chromosome":"18","arm":"p","band":"11.32","iscn_start":"0","iscn_stop":"77","bp_start":"1","bp_stop":"2900000","stain":"gneg","density":""},{"#chromosome":"18","arm":"p","band":"11.31","iscn_start":"77","iscn_stop":"208","bp_start":"2900001","bp_stop":"7200000","stain":"gpos","density":"50"},{"#chromosome":"18","arm":"p","band":"11.2","iscn_start":"623","iscn_stop":"1190","bp_start":"7200001","bp_stop":"15400000","stain":"gneg","density":""},{"#chromosome":"18","arm":"p","band":"11.1","iscn_start":"1190","iscn_stop":"1301","bp_start":"15400001","bp_stop":"18500000","stain":"acen","density":""},{"#chromosome":"18","arm":"q","band":"11.1","iscn_start":"1301","iscn_stop":"1433","bp_start":"18500001","bp_stop":"21500000","stain":"acen","density":""},{"#chromosome":"18","arm":"q","band":"11.2","iscn_start":"1433","iscn_stop":"1862","bp_start":"21500001","bp_stop":"27500000","stain":"gneg","density":""},{"#chromosome":"18","arm":"q","band":"12.1","iscn_start":"1862","iscn_stop":"2223","bp_start":"27500001","bp_stop":"35100000","stain":"gpos","density":"100"},{"#chromosome":"18","arm":"q","band":"12.2","iscn_start":"2223","iscn_stop":"2419","bp_start":"35100001","bp_stop":"39500000","stain":"gneg","density":""},{"#chromosome":"18","arm":"q","band":"12.3","iscn_start":"2419","iscn_stop":"2721","bp_start":"39500001","bp_stop":"45900000","stain":"gpos","density":"75"},{"#chromosome":"18","arm":"q","band":"21.1","iscn_start":"2721","iscn_stop":"3010","bp_start":"45900001","bp_stop":"50700000","stain":"gneg","density":""},{"#chromosome":"18","arm":"q","band":"21.2","iscn_start":"3010","iscn_stop":"3183","bp_start":"50700001","bp_stop":"56200000","stain":"gpos","density":"75"},{"#chromosome":"18","arm":"q","band":"21.3","iscn_start":"3183","iscn_stop":"3449","bp_start":"56200001","bp_stop":"63900000","stain":"gneg","density":""},{"#chromosome":"18","arm":"q","band":"22","iscn_start":"3449","iscn_stop":"4229","bp_start":"63900001","bp_stop":"75400000","stain":"gpos","density":"100"},{"#chromosome":"18","arm":"q","band":"23","iscn_start":"4229","iscn_stop":"4650","bp_start":"75400001","bp_stop":"80373285","stain":"gneg","density":""},{"#chromosome":"19","arm":"p","band":"13.3","iscn_start":"0","iscn_stop":"668","bp_start":"1","bp_stop":"6900000","stain":"gneg","density":""},{"#chromosome":"19","arm":"p","band":"13.2","iscn_start":"668","iscn_stop":"1069","bp_start":"6900001","bp_stop":"12600000","stain":"gpos","density":"25"},{"#chromosome":"19","arm":"p","band":"13.1","iscn_start":"1069","iscn_stop":"1461","bp_start":"12600001","bp_stop":"19900000","stain":"gneg","density":""},{"#chromosome":"19","arm":"p","band":"12","iscn_start":"1461","iscn_stop":"1746","bp_start":"19900001","bp_stop":"24200000","stain":"gvar","density":""},{"#chromosome":"19","arm":"p","band":"11","iscn_start":"1746","iscn_stop":"1880","bp_start":"24200001","bp_stop":"26200000","stain":"acen","density":""},{"#chromosome":"19","arm":"q","band":"11","iscn_start":"1880","iscn_stop":"2019","bp_start":"26200001","bp_stop":"28100000","stain":"acen","density":""},{"#chromosome":"19","arm":"q","band":"12","iscn_start":"2019","iscn_stop":"2342","bp_start":"28100001","bp_stop":"31900000","stain":"gvar","density":""},{"#chromosome":"19","arm":"q","band":"13.1","iscn_start":"2342","iscn_stop":"2943","bp_start":"31900001","bp_stop":"38200000","stain":"gneg","density":""},{"#chromosome":"19","arm":"q","band":"13.2","iscn_start":"2943","iscn_stop":"3344","bp_start":"38200001","bp_stop":"42900000","stain":"gpos","density":"100"},{"#chromosome":"19","arm":"q","band":"13.3","iscn_start":"3344","iscn_stop":"3710","bp_start":"42900001","bp_stop":"50900000","stain":"gneg","density":""},{"#chromosome":"19","arm":"q","band":"13.4","iscn_start":"3710","iscn_stop":"4120","bp_start":"50900001","bp_stop":"58617616","stain":"gpos","density":"100"},{"#chromosome":"20","arm":"p","band":"13","iscn_start":"0","iscn_stop":"471","bp_start":"1","bp_stop":"5100000","stain":"gneg","density":""},{"#chromosome":"20","arm":"p","band":"12","iscn_start":"471","iscn_stop":"1005","bp_start":"5100001","bp_stop":"17900000","stain":"gpos","density":"100"},{"#chromosome":"20","arm":"p","band":"11.2","iscn_start":"1005","iscn_stop":"1593","bp_start":"17900001","bp_stop":"25700000","stain":"gneg","density":""},{"#chromosome":"20","arm":"p","band":"11.1","iscn_start":"1593","iscn_stop":"1729","bp_start":"25700001","bp_stop":"28100000","stain":"acen","density":""},{"#chromosome":"20","arm":"q","band":"11.1","iscn_start":"1729","iscn_stop":"1870","bp_start":"28100001","bp_stop":"30400000","stain":"acen","density":""},{"#chromosome":"20","arm":"q","band":"11.2","iscn_start":"1870","iscn_stop":"2362","bp_start":"30400001","bp_stop":"39000000","stain":"gneg","density":""},{"#chromosome":"20","arm":"q","band":"12","iscn_start":"2362","iscn_stop":"2732","bp_start":"39000001","bp_stop":"43100000","stain":"gpos","density":"75"},{"#chromosome":"20","arm":"q","band":"13.1","iscn_start":"2732","iscn_stop":"3075","bp_start":"43100001","bp_stop":"51200000","stain":"gneg","density":""},{"#chromosome":"20","arm":"q","band":"13.2","iscn_start":"3075","iscn_stop":"3409","bp_start":"51200001","bp_stop":"56400000","stain":"gpos","density":"75"},{"#chromosome":"20","arm":"q","band":"13.3","iscn_start":"3409","iscn_stop":"3770","bp_start":"56400001","bp_stop":"64444167","stain":"gneg","density":""},{"#chromosome":"21","arm":"p","band":"13","iscn_start":"0","iscn_stop":"301","bp_start":"1","bp_stop":"3100000","stain":"gvar","density":""},{"#chromosome":"21","arm":"p","band":"12","iscn_start":"301","iscn_stop":"601","bp_start":"3100001","bp_stop":"7000000","stain":"stalk","density":""},{"#chromosome":"21","arm":"p","band":"11.2","iscn_start":"601","iscn_stop":"1002","bp_start":"7000001","bp_stop":"10900000","stain":"gvar","density":""},{"#chromosome":"21","arm":"p","band":"11.1","iscn_start":"1002","iscn_stop":"1203","bp_start":"10900001","bp_stop":"12000000","stain":"acen","density":""},{"#chromosome":"21","arm":"q","band":"11.1","iscn_start":"1203","iscn_stop":"1341","bp_start":"12000001","bp_stop":"13000000","stain":"acen","density":""},{"#chromosome":"21","arm":"q","band":"11.2","iscn_start":"1341","iscn_stop":"1514","bp_start":"13000001","bp_stop":"15000000","stain":"gneg","density":""},{"#chromosome":"21","arm":"q","band":"21","iscn_start":"1514","iscn_stop":"2204","bp_start":"15000001","bp_stop":"30200000","stain":"gpos","density":"100"},{"#chromosome":"21","arm":"q","band":"22.1","iscn_start":"2204","iscn_stop":"2631","bp_start":"30200001","bp_stop":"38300000","stain":"gneg","density":""},{"#chromosome":"21","arm":"q","band":"22.2","iscn_start":"2631","iscn_stop":"2808","bp_start":"38300001","bp_stop":"41200000","stain":"gpos","density":"50"},{"#chromosome":"21","arm":"q","band":"22.3","iscn_start":"2808","iscn_stop":"3200","bp_start":"41200001","bp_stop":"46709983","stain":"gneg","density":""},{"#chromosome":"22","arm":"p","band":"13","iscn_start":"0","iscn_stop":"308","bp_start":"1","bp_stop":"4300000","stain":"gvar","density":""},{"#chromosome":"22","arm":"p","band":"12","iscn_start":"308","iscn_stop":"604","bp_start":"4300001","bp_stop":"9400000","stain":"stalk","density":""},{"#chromosome":"22","arm":"p","band":"11.2","iscn_start":"604","iscn_stop":"1093","bp_start":"9400001","bp_stop":"13700000","stain":"gvar","density":""},{"#chromosome":"22","arm":"p","band":"11.1","iscn_start":"1093","iscn_stop":"1232","bp_start":"13700001","bp_stop":"15000000","stain":"acen","density":""},{"#chromosome":"22","arm":"q","band":"11.1","iscn_start":"1232","iscn_stop":"1371","bp_start":"15000001","bp_stop":"17400000","stain":"acen","density":""},{"#chromosome":"22","arm":"q","band":"11.2","iscn_start":"1371","iscn_stop":"1887","bp_start":"17400001","bp_stop":"25500000","stain":"gneg","density":""},{"#chromosome":"22","arm":"q","band":"12.1","iscn_start":"1887","iscn_stop":"2000","bp_start":"25500001","bp_stop":"29200000","stain":"gpos","density":"50"},{"#chromosome":"22","arm":"q","band":"12.2","iscn_start":"2000","iscn_stop":"2123","bp_start":"29200001","bp_stop":"31800000","stain":"gneg","density":""},{"#chromosome":"22","arm":"q","band":"12.3","iscn_start":"2123","iscn_stop":"2287","bp_start":"31800001","bp_stop":"37200000","stain":"gpos","density":"50"},{"#chromosome":"22","arm":"q","band":"13.1","iscn_start":"2287","iscn_stop":"2596","bp_start":"37200001","bp_stop":"40600000","stain":"gneg","density":""},{"#chromosome":"22","arm":"q","band":"13.2","iscn_start":"2596","iscn_stop":"2782","bp_start":"40600001","bp_stop":"43800000","stain":"gpos","density":"50"},{"#chromosome":"22","arm":"q","band":"13.3","iscn_start":"2782","iscn_stop":"3400",
                "bp_start":"43800001","bp_stop":"50818468","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"22.3","iscn_start":"0","iscn_stop":"428","bp_start":"1","bp_stop":"9600000","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"22.2","iscn_start":"428","iscn_stop":"642","bp_start":"9600001","bp_stop":"17400000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"p","band":"22.1","iscn_start":"642","iscn_stop":"1131","bp_start":"17400001","bp_stop":"24900000","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"21.3","iscn_start":"1131","iscn_stop":"1454","bp_start":"24900001","bp_stop":"29300000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"p","band":"21.2","iscn_start":"1454","iscn_stop":"1575","bp_start":"29300001","bp_stop":"31500000","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"21.1","iscn_start":"1575","iscn_stop":"1977","bp_start":"31500001","bp_stop":"37800000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"p","band":"11.4","iscn_start":"1977","iscn_stop":"2252","bp_start":"37800001","bp_stop":"42500000","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"11.3","iscn_start":"2252","iscn_stop":"2446","bp_start":"42500001","bp_stop":"47600000","stain":"gpos","density":"75"},{"#chromosome":"X","arm":"p","band":"11.23","iscn_start":"2729","iscn_stop":"2912","bp_start":"47600001","bp_stop":"50100000","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"11.22","iscn_start":"2912","iscn_stop":"3014","bp_start":"50100001","bp_stop":"54800000","stain":"gpos","density":"25"},{"#chromosome":"X","arm":"p","band":"11.21","iscn_start":"3014","iscn_stop":"3057","bp_start":"54800001","bp_stop":"58100000","stain":"gneg","density":""},{"#chromosome":"X","arm":"p","band":"11.1","iscn_start":"3108","iscn_stop":"3241","bp_start":"58100001","bp_stop":"61000000","stain":"acen","density":""},{"#chromosome":"X","arm":"q","band":"11","iscn_start":"3241","iscn_stop":"3341","bp_start":"61000001","bp_stop":"65400000","stain":"acen","density":""},{"#chromosome":"X","arm":"q","band":"12","iscn_start":"3341","iscn_stop":"3672","bp_start":"65400001","bp_stop":"68500000","stain":"gpos","density":"50"},{"#chromosome":"X","arm":"q","band":"13","iscn_start":"3672","iscn_stop":"4405","bp_start":"68500001","bp_stop":"76800000","stain":"gneg","density":""},{"#chromosome":"X","arm":"q","band":"21.1","iscn_start":"4405","iscn_stop":"4740","bp_start":"76800001","bp_stop":"85400000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"q","band":"21.2","iscn_start":"4740","iscn_stop":"4830","bp_start":"85400001","bp_stop":"87000000","stain":"gneg","density":""},{"#chromosome":"X","arm":"q","band":"21.3","iscn_start":"4830","iscn_stop":"5559","bp_start":"87000001","bp_stop":"99100000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"q","band":"22.1","iscn_start":"5559","iscn_stop":"5730","bp_start":"99100001","bp_stop":"103300000","stain":"gneg","density":""},{"#chromosome":"X","arm":"q","band":"22.2","iscn_start":"5730","iscn_stop":"5820","bp_start":"103300001","bp_stop":"104500000","stain":"gpos","density":"50"},{"#chromosome":"X","arm":"q","band":"22.3","iscn_start":"5820","iscn_stop":"5951","bp_start":"104500001","bp_stop":"109400000","stain":"gneg","density":""},{"#chromosome":"X","arm":"q","band":"23","iscn_start":"5951","iscn_stop":"6352","bp_start":"109400001","bp_stop":"117400000","stain":"gpos","density":"75"},{"#chromosome":"X","arm":"q","band":"24","iscn_start":"6352","iscn_stop":"6653","bp_start":"117400001","bp_stop":"121800000","stain":"gneg","density":""},{"#chromosome":"X","arm":"q","band":"25","iscn_start":"6653","iscn_stop":"7195","bp_start":"121800001","bp_stop":"129500000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"q","band":"26","iscn_start":"7195","iscn_stop":"7647","bp_start":"129500001","bp_stop":"138900000","stain":"gneg","density":""},{"#chromosome":"X","arm":"q","band":"27","iscn_start":"7647","iscn_stop":"8168","bp_start":"138900001","bp_stop":"148000000","stain":"gpos","density":"100"},{"#chromosome":"X","arm":"q","band":"28","iscn_start":"8168","iscn_stop":"8610","bp_start":"148000001","bp_stop":"156040895","stain":"gneg","density":""},{"#chromosome":"Y","arm":"p","band":"11.3","iscn_start":"0","iscn_stop":"286","bp_start":"1","bp_stop":"600000","stain":"gpos","density":"100"},{"#chromosome":"Y","arm":"p","band":"11.2","iscn_start":"286","iscn_stop":"751","bp_start":"600001","bp_stop":"10300000","stain":"gneg","density":""},{"#chromosome":"Y","arm":"p","band":"11.1","iscn_start":"751","iscn_stop":"883","bp_start":"10300001","bp_stop":"10400000","stain":"acen","density":""},{"#chromosome":"Y","arm":"q","band":"11.1","iscn_start":"883","iscn_stop":"1002","bp_start":"10400001","bp_stop":"10600000","stain":"acen","density":""},{"#chromosome":"Y","arm":"q","band":"11.21","iscn_start":"1143","iscn_stop":"1268","bp_start":"10600001","bp_stop":"12400000","stain":"gneg","density":""},{"#chromosome":"Y","arm":"q","band":"11.221","iscn_start":"1268","iscn_stop":"1568","bp_start":"12400001","bp_stop":"17100000","stain":"gpos","density":"50"},{"#chromosome":"Y","arm":"q","band":"11.222","iscn_start":"1568","iscn_stop":"1727","bp_start":"17100001","bp_stop":"19600000","stain":"gneg","density":""},{"#chromosome":"Y","arm":"q","band":"11.223","iscn_start":"1727","iscn_stop":"1992","bp_start":"19600001","bp_stop":"23800000","stain":"gpos","density":"50"},{"#chromosome":"Y","arm":"q","band":"11.23","iscn_start":"1992","iscn_stop":"2169","bp_start":"23800001","bp_stop":"26600000","stain":"gneg","density":""},{"#chromosome":"Y","arm":"q","band":"12","iscn_start":"2169","iscn_stop":"3650","bp_start":"26600001","bp_stop":"57227415","stain":"gvar","density":""}];

    var filteredResults = filterByChromosome(data, chr);
    cb(filteredResults);

    // loadData(baseDir + fileName, resolution, function (d) {

    //   if(d) {
    //     var filteredResults = filterByChromosome(d, chr);
    //     cb(filteredResults);
    //   }

    // });
  }

  function filterByChromosome(data, chr) {
    var newAry = [];
    for(var i = 0; i < data.length; i++) {
      if (data[i]['#chromosome'] === chr) {
        newAry.push(data[i]);
      }
    }
    return newAry;
  }

  cyto_chr.modelLoader = {
    load: getChromosomeData,
    setDataDir: function(d) {baseDir = d;},
    getDataDir: function() {return baseDir;}
  };

})(cyto_chr || {}, d3);

(function(cyto_chr, d3) {

  var Selector = function(closecb) {
    this._brush = d3.svg.brush();
    this.dispatch = d3.dispatch('change', 'changeend', 'selectordelete', "selectorhover", "selectorunhover");
    this._x = 0;
    this._y = 0;
    this._extent = [0,0];
    this._height = 0;
    this._closecb = closecb;
  };

  Selector.prototype.test = function(e) {
    var self = this;
    return cyto_chr.InitGetterSetter.call(this, "_test", e, function(){
      self._another = "_that";
    });
  };

  Selector.prototype.extent = function (a) {
    var self = this;
    if(typeof a === "undefined") {
      return self._brush.extent();
    }

    return cyto_chr.InitGetterSetter.call(this, "_extent", a, function(){
      self._brush.extent(a);
    });

  };

  Selector.prototype.xscale = function(a) {
    var self = this;
    return cyto_chr.InitGetterSetter.call(this, "_xscale", a, function(){
      self._brush.x(a);
    });
  };

  Selector.prototype.target = function(a) {
    return cyto_chr.InitGetterSetter.call(this, '_target', a);
  };

  Selector.prototype.height = function(a) {
    return cyto_chr.InitGetterSetter.call(this, '_height', a);
  };

  Selector.prototype.x = function(a) {
    return cyto_chr.InitGetterSetter.call(this, '_x', a);
  };

  Selector.prototype.y = function(a) {
    return cyto_chr.InitGetterSetter.call(this, '_y', a);
  };

  Selector.prototype.render = function() {
    var self = this;

    this.selector = this._target.append('g')
      .classed("extRect", true)
      .attr('transform', 'translate(' + this._x + ',' + this._y + ')')
      .call(this._brush);
    
    this.selector.selectAll(".resize").remove();

    this.selector.selectAll('rect')
      .attr('height', this._height);

    this.selector.select('.background').remove();

    var e = this.selector.select('.extent')
      .style('fill', 'steelblue')
      .style('opacity', '0.5');

    this.selector.selectAll("rect").classed("extent", false).style("cursor", "default");

    this.selector.selectAll("rect")
      .on('mouseover', function() {
        var ext = self._brush.extent();
        self.dispatch.selectorhover(ext);
      })
      .on("mouseout", function() {
        var ext = self._brush.extent();
        self.dispatch.selectorunhover(ext);
      });

    // var cbg_xpos = this._xscale(this._extent[1]) + cyto_chr.margin.left;
    // var cbg_ypos = cyto_chr.margin.top - 3;
    // var cbg = this._target.append('g');
    // cbg.append('title').text('remove');

    // this.deleteButton = cbg.append('circle')
    //   .attr('cx', cbg_xpos)
    //   .attr('cy', cbg_ypos)
    //   .attr('r', 5)
    //   .attr('fill', 'red')
    //   .on('mouseover', function() {
    //     d3.select(this)
    //       .style('cursor', 'pointer');
    //   })
    //   .on('mouseout', function(){
    //     d3.select(this)
    //       .style('cursor', 'default');
    //   })
    //   .on('click', function() {
    //     self.remove();
    //   });

    // this._brush.on('brush', function() {
    //   self.updateXButton();
    //   var ext = self._brush.extent();
    //   // self.dispatch.change(ext);
    // });

    // this._brush.on('brushend', function(d) {
    //   var ext = self._brush.extent();
    //   self.dispatch.changeend(ext);
    // });
    return this;
  };

  Selector.prototype.remove = function() {
    this.selector.remove();
    // this.deleteButton.remove();
    this._closecb(this);
    return this;
  };

  // Selector.prototype.updateXButton = function() {
  //   var e = this._brush.extent();
  //   var new_xpos = this._xscale(e[1]) + cyto_chr.margin.left;
  //   // this.deleteButton.attr('cx', new_xpos);
  // };

  // Selector.prototype.move = function(start, stop) {
  //   this._brush.extent([start, stop]);
  //   this.selector.call(this._brush);
  //   this.updateXButton();
  // };

  cyto_chr.selector = function(cb){
    return new Selector(cb);
  };

})(cyto_chr || {}, d3);

(function(cyto_chr, d3) {

  cyto_chr.initPattern = function () {
    var pg = this.append('pattern')
      .attr('id', 'acen-fill')
      .attr('patternUnits', 'userSpaceOnUse')
      .attr('x', '0')
      .attr('y', '0')
      .attr('width', '10')
      .attr('height', '10')
      .append('g')
      .style({
        "fill": "none",
        "stroke": "#708090",
        "stroke-width": "2"
      });

    pg.append('path')
      .attr('d', "M0,0 l10,10");
    pg.append('path')
      .attr('d','M10,0 l-10,10');
  };

  cyto_chr.roundedRect = function (x, y, w, h, r, tl, tr, bl, br) {
      var retval;
      retval = "M" + (x + r) + "," + y;
      retval += "h" + (w - 2 * r);
      if (tr) {
        retval += "a" + r + "," + r + " 0 0 1 " + r + "," + r;
      } else {
        retval += "h" + r;
        retval += "v" + r;
      }
      retval += "v" + (h - 2 * r);
      if (br) {
        retval += "a" + r + "," + r + " 0 0 1 " + -r + "," + r;
      } else {
        retval += "v" + r;
        retval += "h" + -r;
      }
      retval += "h" + (2 * r - w);
      if (bl) {
        retval += "a" + r + "," + r + " 0 0 1 " + -r + "," + -r;
      } else {
        retval += "h" + -r;
        retval += "v" + -r;
      }
      retval += "v" + (2 * r - h);
      if (tl) {
        retval += "a" + r + "," + r + " 0 0 1 " + r + "," + -r;
      } else {
        retval += "v" + -r;
        retval += "h" + r;
      }
      retval += "z";
      return retval;
    };

  cyto_chr.getStainColour = function (bandtype, density) {

    if(bandtype === "gpos") {
      if(density === "" || density === null) { return "#000000"; }

      switch(density) {
        case "100":
          return "#000000";
        case "75":
          return "#666666";
        case "50":
          return "#999999";
        case "25":
          return "#d9d9d9";
      }
    }

    if (bandtype === "gneg") {
      return "#ffffff";
    }

    if (bandtype === "acen") {
      //return "url(#acen-fill)";
      return "#708090";
    }

    if (bandtype === "gvar") {
      return "#e0e0e0";
    }

    if(bandtype === "stalk") {
      return "#708090";
    }

    return "green";
  };

  cyto_chr.setOption = function (userOption, def) {
      if(typeof userOption !== "undefined") {
        return userOption;
      } else {
        return def;
      }
    };

  cyto_chr.InitGetterSetter = function(prop, arg, cb) {
    if(typeof arg !== 'undefined') {
      this[prop] =  arg;
      if(typeof cb === 'function') {
        cb();
      }
      return this;
    } else {
      return this[prop];
    }
  };

})(cyto_chr || {}, d3);
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
}());