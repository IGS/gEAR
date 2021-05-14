window.onload=function() {
    // check if the user is already logged in

    window.addEventListener('WebComponentsReady', function() {
        // show body now that everything is ready
        files = {};
        files["display_id"] = getUrlParameter("dataset_id");
    
        $.ajax({
            url: './cgi/get_dataset_display.cgi',
            type: "POST",
            data: files,
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                if (data) {
                    if (data.id == files["display_id"] && data.plot_type == "epiviz") {
                        add_tracks(data.plotly_config);
                    }
                }
            },
            error: function (jqXHR, textStatus, errorThrown) {
                console.log('jqXHR: ', jqXHR);
                console.log('textStatus: ', textStatus);
                console.log('errorThrown: ', errorThrown);
    
                var msg = "Unable to query for datasets. Please try again.";
                alert_user(msg);
            }
        }); //end ajax
      });
};

/*jslint unparam: true */
/*global window, $ */
$(function () {
    'use strict';
});

// from common.js
function getUrlParameter(sParam) {
    var sPageURL = decodeURIComponent(window.location.search.substring(1)),
        sURLVariables = sPageURL.split('&'),
        sParameterName,
        i;

    for (i = 0; i < sURLVariables.length; i++) {
        sParameterName = sURLVariables[i].split('=');

        if (sParameterName[0] === sParam) {
            return sParameterName[1] === undefined ? true : sParameterName[1];
        }
    }
};

function add_tracks(config) {

    var extendRangeRatio = 10;
    // var epivizNav = document.querySelector("epiviz-navigation");

    var epivizNav = document.createElement("epiviz-navigation");
    epivizNav.chr = getUrlParameter("chr");
    // epivizNav.range = epivizNav.getGenomicRange(epivizNav.chr, epivizNav.start, epivizNav.end);

    var pstart = parseInt(getUrlParameter("start"));
    var pend = parseInt(getUrlParameter("end"));

    const nstart = pstart - Math.round((pend - pstart) * extendRangeRatio);
    const nend = pend + Math.round((pend - pstart) * extendRangeRatio);

    epivizNav.start = nstart;
    epivizNav.end = nend;

    epivizNav.hideSearch = true;

//     <epiviz-navigation
//     chr="chr10"
//     start=31218495
//     end=31429814
//     hide-search>
// </epiviz-navigation>

    document.querySelector("#epiviz_container").appendChild(epivizNav);

    if (epivizNav) {
        for (const track in config.tracks) {
            const track_config = config.tracks[track];
            track_config.forEach(function(tc) {
                var ttrack = document.createElement(track);
                ttrack.slot = "charts";
                ttrack.setAttribute("measurements", JSON.stringify(tc.measurements));

                if (tc.colors != null) {
                    ttrack.setAttribute("chart-colors", JSON.stringify(tc.colors));
                }
            
                if (tc.settings != null) {
                    ttrack.setAttribute("chart-settings", JSON.stringify(tc.settings));
                }
            
                ttrack.style = "min-height:200px;";
            
                epivizNav.appendChild(ttrack);
            });
        }
    }
  }