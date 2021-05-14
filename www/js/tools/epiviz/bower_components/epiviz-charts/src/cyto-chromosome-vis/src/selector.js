
(function(cyto_chr, d3) {

  var Selector = function(closecb) {
    this._brush = d3.svg.brush();
    this.dispatch = d3.dispatch('change', 'changeend', 'selectordelete');
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
      .attr('transform', 'translate(' + this._x + ',' + this._y + ')')
      .call(this._brush);

    this.selector.selectAll('rect')
      .attr('height', this._height);

    this.selector.select('.background').remove();

    var e = this.selector.select('.extent')
      .style('fill', 'steelblue')
      .style('opacity', '0.5');

    var cbg_xpos = this._xscale(this._extent[1]) + cyto_chr.margin.left;
    var cbg_ypos = cyto_chr.margin.top - 3;
    var cbg = this._target.append('g');
    cbg.append('title').text('remove');

    this.deleteButton = cbg.append('circle')
      .attr('cx', cbg_xpos)
      .attr('cy', cbg_ypos)
      .attr('r', 5)
      .attr('fill', 'red')
      .on('mouseover', function() {
        d3.select(this)
          .style('cursor', 'pointer');
      })
      .on('mouseout', function(){
        d3.select(this)
          .style('cursor', 'default');
      })
      .on('click', function() {
        self.remove();
      });

    this._brush.on('brush', function() {
      self.updateXButton();
      var ext = self._brush.extent();
      self.dispatch.change(ext);
    });

    this._brush.on('brushend', function(d) {
      var ext = self._brush.extent();
      self.dispatch.changeend(ext);
    });
    return this;
  };

  Selector.prototype.remove = function() {
    this.selector.remove();
    this.deleteButton.remove();
    this._closecb(this);
    return this;
  };

  Selector.prototype.updateXButton = function() {
    var e = this._brush.extent();
    var new_xpos = this._xscale(e[1]) + cyto_chr.margin.left;
    this.deleteButton.attr('cx', new_xpos);
  };

  Selector.prototype.move = function(start, stop) {
    this._brush.extent([start, stop]);
    this.selector.call(this._brush);
    this.updateXButton();
  };

  cyto_chr.selector = function(cb){
    return new Selector(cb);
  };

})(cyto_chr || {}, d3);