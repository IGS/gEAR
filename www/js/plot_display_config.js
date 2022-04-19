/*
plot_display_config.js - This script can be passed to various scripts that generate plots
    (emphasis on Plotly) to pass in as post-plotting options.
    Examples include trying to adjust layout parameters so that a plot looks better on a particular page.
*/

// These assume a grid_width of 4 and an mg_grid_width of 6

const TICK_LABEL_MAX_LEN_ALLOWED=10  // Any text over this limit will be truncated
const TICK_LABEL_TRUNCATION_LEN=7    // How much of the original text to use (followed by ellipses)


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

// Truncate axis labels to a fixed length. Add a hover property to the label to show the full label. Edits inplace.
function truncateAxisLabels() {
    $('.xaxislayer-above text').each(function() {
        const sublabel = $(this).text().length > TICK_LABEL_MAX_LEN_ALLOWED
            ? `${$(this).text().substring(0, TICK_LABEL_TRUNCATION_LEN)}...`
            : $(this).text();

        $(this).html(`<a style="fill:inherit;">${sublabel}</a>`);
    });
    /*
    const axis_ticktexts = {}

    for (const element in plotLayout) {
        if (element.includes("axis")) {
            const axis = plotLayout[element];
            if ("ticktext" in axis) {
                // If tick label exceeds max alloned length, truncate it and add hover property to show full label
                for (let i = 0; i < axis.ticktext.length; i++) {
                    const fulllabel = axis.ticktext[i];
                    if (axis.ticktext[i].length > TICK_LABEL_MAX_LEN_ALLOWED) {
                        const sublabel = `${fulllabel.substring(0, TICK_LABEL_TRUNCATION_LEN)}...`;
                        axis.ticktext[i] = sublabel;
                    }
                    axis_ticktexts[fulllabel] = axis.ticktext[i];
                }
            }
        }
    }
    */

    /*
    Issues with modifying before plot generation - Cannot store both full and shortened name.  Plotly is very restrictive about what goes in the tick text
    Issues with modifying after plot generation - Plot margin and layout settings are based on the original labels, and shortening the labels does not auto-adjust these.
    */
}

// If tick label is hovered over, give full label in designated hoverspace
$(document).on('mouseover', '.xaxislayer-above a', function () {
    // Mouse enter

    const h5ad_container = $(this).closest('.h5ad-container');
    const hoverarea = h5ad_container.children(".hoverarea");

    console.log($(this));

    // Show the full label when hovering over the truncated label (assuming we are on a page where the element exists)
    if (hoverarea.length) {
        hoverarea.text($(this).data("unformatted") );
    }
})
$(document).on('mouseleave', '.xaxislayer-above a', function () {
    // Mouse out

    const h5ad_container = $(this).closest('.h5ad-container');
    const hoverarea = h5ad_container.children(".hoverarea");

    if (hoverarea.length) {
        hoverarea.text('');
    }
});
