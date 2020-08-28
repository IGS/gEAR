
var COLORS = {};
COLORS['netherlands'] = ['rgb(0,0,255)','rgb(2,2,255)','rgb(4,4,255)','rgb(6,6,255)','rgb(8,8,255)','rgb(10,10,255)','rgb(12,12,255)','rgb(14,14,255)','rgb(16,16,255)','rgb(18,18,255)','rgb(20,20,255)','rgb(22,22,255)','rgb(24,24,255)','rgb(26,26,255)','rgb(28,28,255)','rgb(30,30,255)','rgb(32,32,255)','rgb(34,34,255)','rgb(36,36,255)','rgb(38,38,255)','rgb(40,40,255)','rgb(42,42,255)','rgb(44,44,255)','rgb(46,46,255)','rgb(48,48,255)','rgb(50,50,255)','rgb(52,52,255)','rgb(54,54,255)','rgb(56,56,255)','rgb(58,58,255)','rgb(60,60,255)','rgb(62,62,255)','rgb(64,64,255)','rgb(66,66,255)','rgb(68,68,255)','rgb(70,70,255)','rgb(72,72,255)','rgb(74,74,255)','rgb(76,76,255)','rgb(78,78,255)','rgb(80,80,255)','rgb(82,82,255)','rgb(84,84,255)','rgb(86,86,255)','rgb(88,88,255)','rgb(90,90,255)','rgb(92,92,255)','rgb(94,94,255)','rgb(96,96,255)','rgb(98,98,255)','rgb(100,100,255)','rgb(102,102,255)','rgb(104,104,255)','rgb(106,106,255)','rgb(108,108,255)','rgb(110,110,255)','rgb(112,112,255)','rgb(114,114,255)','rgb(116,116,255)','rgb(118,118,255)','rgb(120,120,255)','rgb(122,122,255)','rgb(124,124,255)','rgb(126,126,255)','rgb(128,128,255)','rgb(130,130,255)','rgb(132,132,255)','rgb(134,134,255)','rgb(136,136,255)','rgb(138,138,255)','rgb(140,140,255)','rgb(142,142,255)','rgb(144,144,255)','rgb(146,146,255)','rgb(148,148,255)','rgb(150,150,255)','rgb(152,152,255)','rgb(154,154,255)','rgb(156,156,255)','rgb(158,158,255)','rgb(160,160,255)','rgb(162,162,255)','rgb(164,164,255)','rgb(166,166,255)','rgb(168,168,255)','rgb(170,170,255)','rgb(172,172,255)','rgb(174,174,255)','rgb(176,176,255)','rgb(178,178,255)','rgb(180,180,255)','rgb(182,182,255)','rgb(184,184,255)','rgb(186,186,255)','rgb(188,188,255)','rgb(190,190,255)','rgb(192,192,255)','rgb(194,194,255)','rgb(196,196,255)','rgb(198,198,255)','rgb(200,200,255)','rgb(202,202,255)','rgb(204,204,255)','rgb(206,206,255)','rgb(208,208,255)','rgb(210,210,255)','rgb(212,212,255)','rgb(214,214,255)','rgb(216,216,255)','rgb(218,218,255)','rgb(220,220,255)','rgb(222,222,255)','rgb(224,224,255)','rgb(226,226,255)','rgb(228,228,255)','rgb(230,230,255)','rgb(232,232,255)','rgb(234,234,255)','rgb(236,236,255)','rgb(238,238,255)','rgb(240,240,255)','rgb(242,242,255)','rgb(244,244,255)','rgb(246,246,255)','rgb(248,248,255)','rgb(250,250,255)','rgb(252,252,255)','rgb(255,253,253)','rgb(255,251,251)','rgb(255,249,249)','rgb(255,247,247)','rgb(255,245,245)','rgb(255,243,243)','rgb(255,241,241)','rgb(255,239,239)','rgb(255,237,237)','rgb(255,235,235)','rgb(255,233,233)','rgb(255,231,231)','rgb(255,229,229)','rgb(255,227,227)','rgb(255,225,225)','rgb(255,223,223)','rgb(255,221,221)','rgb(255,219,219)','rgb(255,217,217)','rgb(255,215,215)','rgb(255,213,213)','rgb(255,211,211)','rgb(255,209,209)','rgb(255,207,207)','rgb(255,205,205)','rgb(255,203,203)','rgb(255,201,201)','rgb(255,199,199)','rgb(255,197,197)','rgb(255,195,195)','rgb(255,193,193)','rgb(255,191,191)','rgb(255,189,189)','rgb(255,187,187)','rgb(255,185,185)','rgb(255,183,183)','rgb(255,181,181)','rgb(255,179,179)','rgb(255,177,177)','rgb(255,175,175)','rgb(255,173,173)','rgb(255,171,171)','rgb(255,169,169)','rgb(255,167,167)','rgb(255,165,165)','rgb(255,163,163)','rgb(255,161,161)','rgb(255,159,159)','rgb(255,157,157)','rgb(255,155,155)','rgb(255,153,153)','rgb(255,151,151)','rgb(255,149,149)','rgb(255,147,147)','rgb(255,145,145)','rgb(255,143,143)','rgb(255,141,141)','rgb(255,139,139)','rgb(255,137,137)','rgb(255,135,135)','rgb(255,133,133)','rgb(255,131,131)','rgb(255,129,129)','rgb(255,127,127)','rgb(255,125,125)','rgb(255,123,123)','rgb(255,121,121)','rgb(255,119,119)','rgb(255,117,117)','rgb(255,115,115)','rgb(255,113,113)','rgb(255,111,111)','rgb(255,109,109)','rgb(255,107,107)','rgb(255,105,105)','rgb(255,103,103)','rgb(255,101,101)','rgb(255,99,99)','rgb(255,97,97)','rgb(255,95,95)','rgb(255,93,93)','rgb(255,91,91)','rgb(255,89,89)','rgb(255,87,87)','rgb(255,85,85)','rgb(255,83,83)','rgb(255,81,81)','rgb(255,79,79)','rgb(255,77,77)','rgb(255,75,75)','rgb(255,73,73)','rgb(255,71,71)','rgb(255,69,69)','rgb(255,67,67)','rgb(255,65,65)','rgb(255,63,63)','rgb(255,61,61)','rgb(255,59,59)','rgb(255,57,57)','rgb(255,55,55)','rgb(255,53,53)','rgb(255,51,51)','rgb(255,49,49)','rgb(255,47,47)','rgb(255,45,45)','rgb(255,43,43)','rgb(255,41,41)','rgb(255,39,39)','rgb(255,37,37)','rgb(255,35,35)','rgb(255,33,33)','rgb(255,31,31)','rgb(255,29,29)','rgb(255,27,27)','rgb(255,25,25)','rgb(255,23,23)','rgb(255,21,21)','rgb(255,19,19)','rgb(255,17,17)','rgb(255,15,15)','rgb(255,13,13)','rgb(255,11,11)','rgb(255,9,9)','rgb(255,7,7)','rgb(255,5,5)','rgb(255,3,3)','rgb(255,1,1)','rgb(255,0,0)'];
var P2_snap_paths = null;

// each element is a time point
var data = [
    // time point label
    // per tissue:
    //   color (from mean)
    //   sampling values
    { 'l': 'Baseline',
      'v': [
          {'t':'Non-sensory_cells', 'c':200, 's': [42881,36382,47346]},
          {'t':'Hair_cells', 'c':150, 's': [37354,33641,24343]}
      ]
    },
    { 'l': '6 hrs',
      'v': [
          {'t':'Non-sensory_cells', 'c':225, 's': [34120,37059,60594]},
          {'t':'Hair_cells', 'c':125, 's': [19315,35393,27311]}
      ]
    },
    { 'l': '24 hrs',
      'v' : [
          {'t':'Non-sensory_cells', 'c':150, 's': [28525,32715,34991]},
          {'t':'Hair_cells', 'c':100, 's': [21512,27178,27763]}
      ]
    }
];

var myBarChart = null;
var chart_data = null;

window.onload=function() {
    var s = Snap("#demo_svg_container");
    img_src = "demo_timecourse.svg";

    $('#demo_svg_container').empty();

    Snap.load(img_src, function (svg) {
        P2_snap_paths = svg.selectAll("path");
        colors_by_tissue = new Object();
        colors_by_tissue['Non-sensory_cells'] = "rgba(220,220,220,0.8)";
        colors_by_tissue['Hair_cells'] = "rgba(151,187,205,0.8)";

        P2_snap_paths.forEach(
            function(elem,i) {
                var elem_class = elem.attr("class");
                if (elem_class.length > 0) {
                    var pieces = elem_class.split(/[\s,]+/);
                    elem_class = pieces[pieces.length-1]

                    if (typeof colors_by_tissue[elem_class] !== 'undefined') {
                        elem.attr("fill", colors_by_tissue[elem_class]);
                    }
                }
            }
        );

        s.append(svg);
    });

    $('#time_point_count').text(data.length);

    // All things chart.
    chart_data = { labels:[], datasets:[] };
    var chart_tissue_values = {};
    var chart_sites = [];

    $.each(data, function(i, item) {
        chart_data['labels'].push(item['l']);
        //chart_data['labels'].push('');
        $.each(item['v'], function(ii, item2) {
            body_site = item2['t'];
            console.log("Setting site = " + body_site);
            if (chart_tissue_values.hasOwnProperty(body_site)) {
                chart_tissue_values[body_site].push(item2['c']);
            } else {
                chart_tissue_values[body_site] = [ item2['c'] ];
                chart_sites.push(item2['t']);
            }
        });
    });

    for (var site in chart_tissue_values) {
        console.log("Processing site " + site);
        chart_data['datasets'].push(
            {
                label: site,
                fillColor: "rgba(220,220,220,0.5)",
                strokeColor: "rgba(220,220,220,0.8)",
                highlightFill: "rgba(220,220,220,0.75)",
                highlightStroke: "rgba(220,220,220,1)",
                data: chart_tissue_values[site]
            }
        );
    }

    var ctx = document.getElementById("myChart").getContext("2d");
    // second argument is options
    myBarChart = new Chart(ctx).Bar(chart_data, {});
    myBarChart.datasets[0].bars[0].fillColor = "rgba(220,220,220,0.8)";
    myBarChart.datasets[0].bars[1].fillColor = "rgba(120,120,120,0.8)";
    myBarChart.datasets[0].bars[2].fillColor = "rgba(20,20,20,0.8)";
    myBarChart.datasets[1].bars[0].fillColor = "rgba(151,187,205,0.8)";
    myBarChart.datasets[1].bars[1].fillColor = "rgba(151,187,205,0.8)";
    myBarChart.datasets[1].bars[2].fillColor = "rgba(151,187,205,0.8)";
    myBarChart.update();
    
    $(function() {
        $( "#slider" ).slider({
            value:0,
            min: 0,
            max: 100,
            step: 50,
            slide: function( event, ui ) {
                colors_by_tissue = new Object();
                
                if (ui.value == 0) {
                    $("#current_time_point").text("Baseline");
                    colors_by_tissue['Non-sensory_cells'] = "rgba(220,220,220,0.8)";
                    myBarChart.datasets[0].bars[0].fillColor = "rgba(220,220,220,0.8)";
                    myBarChart.datasets[0].bars[1].fillColor = "rgba(120,120,120,0.8)";
                    myBarChart.datasets[0].bars[2].fillColor = "rgba(20,20,20,0.8)";
                    colors_by_tissue['Hair_cells'] = "rgba(151,187,205,0.8)";
                    myBarChart.datasets[1].bars[0].fillColor = "rgba(151,187,205,0.8)";
                    myBarChart.datasets[1].bars[1].fillColor = "rgba(151,187,205,0.8)";
                    myBarChart.datasets[1].bars[2].fillColor = "rgba(151,187,205,0.8)";
                } else if (ui.value == 50) {
                    $("#current_time_point").text("6 hours after noise exposure");
                    colors_by_tissue['Non-sensory_cells'] = "rgba(120,120,120,0.8)";
                    myBarChart.datasets[0].bars[0].fillColor = "rgba(220,220,220,0.8)";
                    myBarChart.datasets[0].bars[1].fillColor = "rgba(120,120,120,0.8)";
                    myBarChart.datasets[0].bars[2].fillColor = "rgba(20,20,20,0.8)";
                    colors_by_tissue['Hair_cells'] = "rgba(121,157,175,0.8)";
                    myBarChart.datasets[1].bars[0].fillColor = "rgba(151,187,205,0.8)";
                    myBarChart.datasets[1].bars[1].fillColor = "rgba(121,157,175,0.8)";
                    myBarChart.datasets[1].bars[2].fillColor = "rgba(101,137,155,0.8)";
                } else {
                    $("#current_time_point").text("24 hours after noise exposure");
                    colors_by_tissue['Non-sensory_cells'] = "rgba(20,20,20,0.8)";
                    myBarChart.datasets[0].bars[0].fillColor = "rgba(220,220,220,0.8)";
                    myBarChart.datasets[0].bars[1].fillColor = "rgba(120,120,120,0.8)";
                    myBarChart.datasets[0].bars[2].fillColor = "rgba(20,20,20,0.8)";
                    colors_by_tissue['Hair_cells'] = "rgba(101,137,155,0.8)";
                    myBarChart.datasets[1].bars[0].fillColor = "rgba(151,187,205,0.8)";
                    myBarChart.datasets[1].bars[1].fillColor = "rgba(121,157,175,0.8)";
                    myBarChart.datasets[1].bars[2].fillColor = "rgba(101,137,155,0.8)";
                }

                myBarChart.update();
                
                var s = Snap("#demo_svg_container");

                Snap.load(img_src, function (svg) {
                    P2_snap_paths = svg.selectAll("path");
                    P2_snap_paths.forEach(
                        function(elem,i) {
                            var elem_class = elem.attr("class");
                            if (elem_class.length > 0) {
                                var pieces = elem_class.split(/[\s,]+/);
                                elem_class = pieces[pieces.length-1]

                                if (typeof colors_by_tissue[elem_class] !== 'undefined') {
                                    elem.attr("fill", colors_by_tissue[elem_class]);
                                }
                            }
                        }
                    );

                    $('#demo_svg_container').empty();
                    s.append(svg);
                });
            }
        });
    });
/*
    var data = {
        labels: ["January", "February", "March", "April", "May", "June", "July"],
        datasets: [
            {
                label: "My First dataset",
                fillColor: "rgba(220,220,220,0.5)",
                strokeColor: "rgba(220,220,220,0.8)",
                highlightFill: "rgba(220,220,220,0.75)",
                highlightStroke: "rgba(220,220,220,1)",
                data: [65, 59, 80, 81, 56, 55, 40]
            },
            {
                label: "My Second dataset",
                fillColor: "rgba(151,187,205,0.5)",
                strokeColor: "rgba(151,187,205,0.8)",
                highlightFill: "rgba(151,187,205,0.75)",
                highlightStroke: "rgba(151,187,205,1)",
                data: [28, 48, 40, 19, 86, 27, 90]
            }
        ]
    };
*/

};
