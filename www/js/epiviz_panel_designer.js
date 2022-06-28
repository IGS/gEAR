/*
 This script relies on the source having also included the
 common.js within this project
*/

window.onload=function() {
    // Set the curator link href for when the user gets that far
    var home_url = `./index.html`;
    var data_manager_url = `./dataset_explorer.html`;
    var self = this;
    $("#goHome").click(function() {
        window.location.replace(home_url);
    });
    $("#goDataManager").click(function() {
        window.location.replace(data_manager_url);
    });
};

window.addEventListener('WebComponentsReady', function() {
    var self = this;
    // // $("#continue_container").hide();
    self.setTimeout(function() {
        files = {};
        files["session_id"] = CURRENT_USER.session_id;

        $.ajax({
            url: './cgi/get_epigenetic_dataset_list.cgi',
            type: "POST",
            data: files,
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                if (data.datasets.length > 0) {
                    var measurements = {};
                    var mAll = new epiviz.measurements.MeasurementSet();

                    var epivizNavElem = document.querySelector("epiviz-navigation");

                    data.datasets.forEach(function(m) {
                        type = "range";
                        chart = "blocksTrack";
                        if (m.type == "bigwig" || m.type == "bw") {
                            type = "feature";
                            chart = "lineTrack";
                        }
                        measurements[m.url] = {
                            id: m.id + "." + m.type,
                            name: m.title,
                            type: type,
                            datasourceId: m.url,
                            datasourceGroup: m.id + "." + m.type,
                            dataprovider: "fileapi",
                            formula: null,
                            defaultChartType: chart,
                            annotation: null,
                            minValue: "default",
                            maxValue: "default",
                            metadata: []
                        };

                        mAll.add(new epiviz.measurements.Measurement(
                            m.id + "." + m.type,
                            m.title,
                            type,
                            m.url,
                            m.id + "." + m.type,
                            "fileapi",
                            null,
                            chart,
                            null,
                            "default",
                            "default",
                            []
                        ))
                    });

                    epivizNavElem.measurementSet = measurements;
                    epivizNavElem.measurementsSet = mAll;
                    epivizNavElem.measurementAll = mAll;

                    $("#epiviz_container").removeClass("disable_epiviz");

                } else {
                    var msg = "Unable to query for datasets. Please try again.";
                    alert_user(msg);
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

        changeGenome();
    }, 2000);
});

/*jslint unparam: true */
/*global window, $ */
$(function () {
    'use strict';
});

function changeGenome() {
    var genome = $("#genome").val();
    var nav_elem = document.querySelector("epiviz-navigation");
    var genes_track = nav_elem.querySelector("epiviz-genes-track");

    if (genes_track) {
        genes_track.remove();
    }

    var gtrack = document.createElement("epiviz-genes-track");
    gtrack.slot = "charts";
    gtrack.measurements = [{
        id: genome,
        name: genome,
        type: "range",
        datasourceId: "https://obj.umiacs.umd.edu/genomes/" + genome + "/" + genome + ".txt.gz.tbi",
        datasourceGroup: genome,
        dataprovider: "fileapi",
        formula: null,
        defaultChartType: "GenesTrack",
        annotation: null,
        minValue: "default",
        maxValue: "default",
        metadata: ["geneid", "exons_start", "exons_end", "gene"]
    }];
    nav_elem.appendChild(gtrack);

}

function get_epiviz_profile() {

    var econfig = {
        "genome":"mm10",
        "fullviewer":"http://epiviz.umgear.org/",
        "dataserver":"/epivizfile",
        "tracks":{}
    };

    var tracks = {};
    var navElem = document.querySelector("epiviz-navigation");
    let navChildren = Polymer.FlattenedNodesObserver.getFlattenedNodes(
        navElem
    ).filter(
      (n) =>
        n.nodeType === Node.ELEMENT_NODE && n.nodeName.indexOf("EPIVIZ") !== -1
    );

    var numChildren = navChildren.length;

    if (numChildren == 0) {
        return null;
    }

    for (var index = 0; index < numChildren; index++) {
      var currentChild = navChildren[index];
      if (!Object.keys(tracks).includes(currentChild.nodeName)) {
        tracks[currentChild.nodeName] = [];
      }

      tracks[currentChild.nodeName].push({
          "measurements": currentChild.measurements,
          "settings": currentChild.chartSettings,
          "colors": currentChild.chartColors
      });
    }

    econfig.tracks = tracks;

    return econfig;
}

$("#submit").click(function(e) {
	// remove any warnings, if present
	$('.label-warning').remove();

    // disable the submit button and put the UID into the form
    $("#submit").text("Processing");
    $("#submit").attr("disabled", true);
    $("#epiviz_container").addClass("disable_epiviz");

    e.preventDefault();

    var fdata = {};

    var formData = $("#upload_form").serializeArray();

    formData.forEach(function(m) {
        if (m.name == "dataset_uid") {
            fdata[m.name] = guid('long');
        } else if (m.name == "share_uid") {
            fdata[m.name] = guid('short');
        } else {
            fdata[m.name] = m.value;
        }
    });

     var econfig = get_epiviz_profile();

     if(!econfig) {
        var msg = "Unable to add configuration. No Charts in Epiviz workspace. Please try again.";
        alert_user(msg);

        $("#submit").text("Add Epiviz Panel to gEAR");
        $("#submit").attr("disabled", false);
        $("#epiviz_container").removeClass("disable_epiviz");

     } else {
        fdata["epiviz_config"] = JSON.stringify(econfig);

        $.ajax({
            url: './cgi/add_epiviz_dataset_profile.cgi',
            type: "POST",
            data : fdata,
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                if (data['success'] == 1) {
                    console.log("upload success! check db entries!")
                    $("#upload_form").hide();
                    $("#continue_container").show();
                } else {
                    var msg = "An upload error has occurred. Please try again.";
                    alert_user(msg);

                    //reset submit button area
                    $("#submit").text("Submit");
                    $("#submit").attr("disabled", false);
                }
            },
            error: function (jqXHR, textStatus, errorThrown) {
                console.log('jqXHR: ', jqXHR);
                console.log('textStatus: ', textStatus);
                console.log('errorThrown: ', errorThrown);

                var msg = "Unable to finish dataset upload. Please try again.";
                alert_user(msg);
            }
        }); //end ajax
     }
});

/**
 * Generates a GUID string.
 * @returns {String} The generated GUID.
 * @example af8a8416-6e18-a307-bd9c-f2c947bbb3aa
 * @author Slavik Meltser (slavik@meltser.info).
 * @link http://slavik.meltser.info/?p=142
 */
function guid(uid_length) {
    function _p8(s) {
        var p = (Math.random().toString(16)+"000000000").substr(2,8);
        return s ? "-" + p.substr(0,4) + "-" + p.substr(4,4) : p ;
    }
    if (uid_length == 'long') {
        return _p8() + _p8(true) + _p8(true) + _p8();
    }
    if (uid_length == 'short') {
        return _p8();
    }
}
