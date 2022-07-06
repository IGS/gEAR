/*
plot_display_config.js - This script can be passed to various scripts that generate plots
    (emphasis on Plotly) to pass in as post-plotting options.
    Examples include trying to adjust layout parameters so that a plot looks better on a particular page.
*/

// These assume a grid_width of 4 and an mg_grid_width of 6

const post_plotly_config = {
    index: [
        {
            plot_type: "all"
            , config: {
                showLink: false
                , displaylogo: false
                , responsive: false
                , modeBarButtonsToRemove: [
                    "zoom2d",
                    "autoScale2d",
                    "hoverClosestCartesian",
                    "hoverCompareCartesian",
                    "zoom3d",
                    "pan3d",
                    "resetCameraDefault3d",
                    "resetCameraLastSave3d",
                    "hoverClosest3d",
                    "orbitRotation",
                    "tableRotation",
                    "zoomInGeo",
                    "zoomOutGeo",
                    "resetGeo",
                    "hoverClosestGeo",
                    "sendDataToCloud",
                    "hoverClosestGl2d",
                    "hoverClosestPie",
                    "toggleHover",
                    "resetViews",
                    "toggleSpikelines",
                    "resetViewMapbox"
                ]
            }, layout: {
                modebar: {
                    orientation: "h"
                }
            }
        }, {
            plot_type: "volcano"
            , layout:{
                legend: {
                    x: 1.05,
                    y: -0.05,
                    xanchor: "center",
                    yanchor: "top",
                    font: {
                        size: 8
                    }
                }
                , title : {
                    x: 0.5,
                    xref: "paper",
                    y: 0.8,
                }
            }
        }, {
            plot_type: "heatmap"
            , layout:{
                margin: {
                    b: 10
                }
            }
        }, {
            plot_type: "violin"
            , layout:{
                margin: {
                    t: 30
                    , r: 10
                }
            }
        }, {
            plot_type: "bar"
            , layout:{
                margin: {
                    t: 30
                    , r: 10
                }
                /* , legend: {
                    orientation: "h"
                }
                */
            }
        }
    ]
    , curator: [
        {
            plot_type:"all"
            , config: {
                showLink: false
                , displaylogo: false
                , responsive: false
                , modeBarButtonsToRemove: [
                    "zoom2d",
                    "autoScale2d",
                    "hoverClosestCartesian",
                    "hoverCompareCartesian",
                    "zoom3d",
                    "pan3d",
                    "resetCameraDefault3d",
                    "resetCameraLastSave3d",
                    "hoverClosest3d",
                    "orbitRotation",
                    "tableRotation",
                    "zoomInGeo",
                    "zoomOutGeo",
                    "resetGeo",
                    "hoverClosestGeo",
                    "sendDataToCloud",
                    "hoverClosestGl2d",
                    "hoverClosestPie",
                    "toggleHover",
                    "resetViews",
                    "toggleSpikelines",
                    "resetViewMapbox"
                ]
            }, layout:{
                modebar: {
                    orientation: "v"
                }
                , dragMode: "select"
            }
        }, {
            plot_type: "volcano"
            , layout:{
                "height": 800
            }
        }, {
            plot_type: "quadrant"
            , layout:{
                "height": 800
            }
        }
    ]
}

// Functions that cannot be encapsulated by a general config change
// ! These will modify an object reference in place if the param argument is an object property

function adjustExpressionColorbar(plotData) {
    // The colorbar is outside of the graph div.  Need to adjust to bring back in.
    for (const element of plotData) {
        if ("colorbar" in element && element.name === "expression") {
            element.colorbar.len = 1.5;
            element.colorbar.xpad = 10;
            element.colorbar.x = 1.15;
            element.colorbar.title = {text: "Expression"};
        }
    }
}

function adjustClusterColorbars(plotData) {
    const plot_min_domain = 0;
    const plot_max_domain = 1;
    const new_min_domain = -0.5;
    const new_max_domain = 1.5;

    const num_clusterbars = plotData.filter(element => "colorbar" in element && element.name === "clusterbar").length;

    // The colorbar is outside of the graph div.  Need to adjust to bring back in.
    for (const element of plotData) {
        if ("colorbar" in element && element.name === "clusterbar") {
            element.colorbar.xpad = 10;
            element.colorbar.x = -0.3;
            element.colorbar.xanchor = "left";
            element.colorbar.len = 2 / num_clusterbars;
            // Scale the colorbar y positions to match the plot container instead of the plot
            element.colorbar.y = scaleBetween(element.colorbar.y, new_min_domain, new_max_domain, plot_min_domain, plot_max_domain);
        }
    }
}

function setHeatmapHeightBasedOnGenes(plotLayout, genesFilter) {
    if (genesFilter.length > 20) {
        plotLayout.height = genesFilter.length * 20;
      }
}

function adjustStackedViolinHeight(plotLayout) {
    plotLayout.height = 800;
}

function adjustStackedViolinSharedAxes(plotLayout) {
    return;
}

// Scaling a number from one range to another range (Source: https://stackoverflow.com/a/31687097)
function scaleBetween(unscaledNum, minAllowed, maxAllowed, min, max) {
    return (maxAllowed - minAllowed) * (unscaledNum - min) / (max - min) + minAllowed;
}