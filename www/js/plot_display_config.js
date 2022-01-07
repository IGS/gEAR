/*
plot_display_config.js - This script can be passed to various scripts that generate plots
    (emphasis on Plotly) to pass in as post-plotting options.
    Examples include trying to adjust layout parameters so that a plot looks better on a particular page.
*/

// These assume a grid_width of 4 and an mg_grid_width of 6

const post_plotly_config = {
    "index":[
        {
            "plot_type":"volcano"
            , "layout":{
                "legend": {
                    x: 1.05,
                    y: -0.05,
                    xanchor: "center",
                    yanchor:"top",
                    font: {size:8}
                }
                , "title" : {
                    x: 0.5,
                    xref: "paper",
                    y: 0.8,
                }
            }
        }, {
            "plot_type":"heatmap"
            , "layout":{
                "margin": {
                    "b":10
                }
                /*, "title": {
                    "pad": {
                        "b":10
                    }
                }*/
            }
        }
    ]
    , "curator":[
        {
            "plot_type": "volcano"
            , "layout":{
                "height": 800
            }
        }
        , {
            "plot_type": "quadrant"
            , "layout":{
                "height": 800
            }
        }
    ]
}

// Functions that cannot be encapsulated by a general config change

function adjustExpressionColorbar(plotData) {
    // The colorbar is outside of the graph div.  Need to adjust to bring back in.
    for (let i = 0; i < plotData.length; i++) {
        if ("colorbar" in plotData[i] && plotData[i].name === "expression") {
            plotData[i].colorbar.len = 0.7;
            plotData[i].colorbar.xpad = 10;
            plotData[i].colorbar.x = 1.15;
            plotData[i].colorbar.title = {text: "Expression"};
        }
    }
}

function adjustClusterColorbars(plotData) {
    // The colorbar is outside of the graph div.  Need to adjust to bring back in.
    for (let i = 0; i < plotData.length; i++) {
        if ("colorbar" in plotData[i] && plotData[i].name === "clusterbar") {
            plotData[i].colorbar.xpad = 10;
            plotData[i].colorbar.x = -0.3;
            plotData[i].colorbar.xanchor = "left";
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
