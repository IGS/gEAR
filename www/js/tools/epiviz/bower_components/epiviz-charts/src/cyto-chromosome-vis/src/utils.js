
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