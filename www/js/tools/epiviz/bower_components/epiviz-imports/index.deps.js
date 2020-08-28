'use strict';

var jquery = require('./src/js/lib/jquery/jquery-1.8.2.js');
var d3 = require('d3');
var sprintf = require('sprintf');
var epiviz = require('./node_modules/epiviz/index.min.js');

module.exports = {
    sprintf: sprintf,
    epiviz : epiviz
};

// window.sprintf = sprintf;
// window.epiviz = epiviz;
