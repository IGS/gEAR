
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
    this.dispatch = d3.dispatch('bandclick', 'selectorchange', 'selectorend', 'selectordelete');
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
    if(typeof a === 'string') a = d3.select(a);
    if(a.empty()) {
      throw "Error: Invalid dom target";
    }
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
      .style('font', '10px sans-serif');

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

  Chromosome.prototype.moveSelectorTo = function(start, stop) {
    if(arguments.length !== 2) {
      throw "Error moveSelectorTo: Invalid number of arguments. Both start and stop coordinates are required";
    }

    if (this.selectors.length === 0) {
      this.newSelector(0, "10000000");
    } else {
      this.selectors[0].move(start, stop);
    }

  };

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

    ve.dispatch.on('change', function(d) {
      self.dispatch.selectorchange(d);
    });

    ve.dispatch.on('changeend', function(d) {
      self.dispatch.selectorend(d);
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
      var w = self._width;

      self.svgTarget = self._domTarget
        .style('height', h + 'px')
        .style('width', (w+5) + 'px')
        .append('svg')
        .attr('width', w)
        .attr('height', h);

      var bands = self.svgTarget.selectAll('g')
        .data(data).enter();

      self.svgTarget.append('text')
        .text(self._segment)
        .attr('x', 5)
        .attr('y', cyto_chr.margin.top + (self._height/ 2) + 2)
        .attr('text-anchor','middle')
        .style('font', '10px sans-serif');

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

          if(i % labelSkipFactor === 0) {
            var bmid = (bpCoord(d.bp_stop) + bpCoord(d.bp_start)) / 2;
            elem.append('line')
              .attr('x1', bmid)
              .attr('y1', cyto_chr.margin.top)
              .attr('x2', bmid)
              .attr('y2', cyto_chr.margin.top - 4)
              .style('stroke', 'grey')
              .style('stroke-width',1);

            elem.append('text')
              .attr('transform', 'translate(' + bmid + ',' + (cyto_chr.margin.top - 6) + ')rotate(-50)')
              .style('font', '10px sans-serif')
              .text(d.arm + d.band);
          }

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

          rect.on('click', function(d) {
            if (self.selectors.length === 0 || (self._selectionMode === 'multi' && d3.event.altKey)) {
              self.newSelector(d.bp_start, d.bp_stop);
            }

            if (self._selectionMode === 'single' && self.selectors.length > 0) {
              self.moveSelectorTo(d.bp_start, d.bp_stop);
            }
            self.dispatch.bandclick(d);
          });
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
  };1

})(cyto_chr || {}, d3);