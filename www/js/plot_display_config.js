/*
plot_display_config.js - This script can be passed to various scripts that generate plots
    (emphasis on Plotly) to pass in as post-plotting options.
    Examples include trying to adjust layout parameters so that a plot looks better on a particular page.
*/

// These assume a grid_width of 4 and an mg_grid_width of 6

const TICK_LABEL_MAX_LEN_ALLOWED=10  // Any text over this limit will be truncated
const TICK_LABEL_TRUNCATION_LEN=7    // How much of the original text to use (followed by ellipses)

// color palettes
// Obtained from https://plotly.com/python/builtin-colorscales/
// and https://plotly.com/python/discrete-color/
// The "dotplot" property means these are possible for dotplots
const availablePalettes = [
    {
        label: "Qualitative Scales",
        continuous: false,
        options: [
            // These need to be kept up-to-date with the "color_swatch_map" in lib/mg_plotting.py
            { value: "alphabet", text: "Alphabet (26 colors)" },
            { value: "bold", text: "Bold (11 colors)" },
            { value: "d3", text: "D3 (10 colors)" },    // default for categorical plots
            { value: "dark24", text: "Dark (24 colors)" },
            { value: "light24", text: "Light (24 colors)" },
            { value: "safe", text: "Safe (11 colors)" },
            { value: "vivid", text: "Vivid (11 colors)" },
        ]
    },
    {
        label: "Sequential Scales",
        continuous: true,
        options: [
            { value: "greys", text: "Greys" },
            { value: "blues", text: "Blues" },
            { value: "purp", text: "Purples" }, // default for single-gene continouous plotly
            { value: "reds", text: "Reds" },
            { value: "bluered", text: "Blue-Red" }, // default for multigene continouous plotly
            { value: "ylgnbu", text: "Yellow-Green-Blue" },
            { value: "ylorrd", text: "Yellow-Orange-Red" }, // default for scanpy plots
        ],
    },
    {
        label: "Diverging Scales",
        continuous: true,
        options: [
            { value: "brbg", text: "Brown-Blue-Green" },
            { value: "piyg", text: "Pink-Green" },
            { value: "prgn", text: "Purple-Green" },
            { value: "rdbu", text: "Red-Blue" },
            { value: "rdylbu", text: "Red-Yellow-Blue" },
            { value: "bublrd", text: "Blue-Black-Red"},
            { value: "multicolor_diverging", text: "Purple-Blue-Black-Red-Yellow" },
        ],
    },
    {
        label: "Color Vision Accessibility Scales",
        continuous: true,
        options: [
            { value: "cividis", text: "Cividis" },
            { value: "inferno", text: "Inferno" },
            { value: "viridis", text: "Viridis" },
        ]
    },
];

// Convert the plotly colorscale names from availablePalettes to their matplotlib equivalent names
const plotly2MatplotlibNames = {
    "cividis": "cividis",
    "inferno": "inferno",
    "viridis": "viridis",
    "greys": "Greys",
    "blues": "Blues",
    "purp": "Purples",
    "reds": "Reds",
    "bluered": "bluered",   // Custom colorscale
    "ylgnbu": "YlGnBu",
    "ylorrd": "YlOrRd",
    "brbg": "BrBG",
    "piyg": "PiYG",
    "prgn": "PRGn",
    "rdbu": "RdBu",
    "rdylbu": "RdYlBu",
    // anything below is custom colorscales
    "bublrd": "bublrd",
    "multicolor_diverging": "multicolor_diverging"
}

// invert the plotly2MatplotlibNames object
const matplotlib2PlotlyNames = Object.fromEntries(Object.entries(plotly2MatplotlibNames).map(([key, value]) => [value, key]));

const postPlotlyConfig = {
    expression: [
        {
            plot_type: "all"
            , config: {
                showLink: false
                , displaylogo: false
                , responsive: true
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
                autosize: true,
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
        }
    ]
    , curator: [
        {
            plot_type:"all"
            , config: {
                showLink: false
                , displaylogo: false
                , responsive: true
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
                dragMode: "select"
                , height: window.innerHeight * 0.8
                , autosize: true
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
    const plotMinDomain = 0;
    const plotMaxDomain = 1;
    const newMinDomain = -0.5;
    const newMaxDomain = 1.5;

    const numClusterbars = plotData.filter(element => "colorbar" in element && element.name === "clusterbar").length;

    // The colorbar is outside of the graph div.  Need to adjust to bring back in.
    for (const element of plotData) {
        if ("colorbar" in element && element.name === "clusterbar") {
            element.colorbar.xpad = 10;
            element.colorbar.x = -0.3;
            element.colorbar.xanchor = "left";
            element.colorbar.len = 2 / numClusterbars;
            // Scale the colorbar y positions to match the plot container instead of the plot
            element.colorbar.y = scaleBetween(element.colorbar.y, newMinDomain, newMaxDomain, plotMinDomain, plotMaxDomain);
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

// Truncate axis labels to a fixed length. Add a hover property to the label to show the full label. Edits inplace.
function truncateAxisLabels(plotLayout) {
    const selector = document.querySelectorAll('.xaxislayer-above text');
    for (const el of selector) {
        const elText = el.textContent;
        const sublabel = elText.length > TICK_LABEL_MAX_LEN_ALLOWED
            ? `${elText.substring(0, TICK_LABEL_TRUNCATION_LEN)}...`
            : elText;


        const subLabelElt = document.createElement('a');
        subLabelElt.setAttribute('style', 'fill:inherit;');
        subLabelElt.textContent = sublabel;
        el.textContent = "";
        el.appendChild(subLabelElt);
    }

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

    /*
    Issues with modifying before plot generation - Cannot store both full and shortened name.  Plotly is very restrictive about what goes in the tick text
    Issues with modifying after plot generation - Plot margin and layout settings are based on the original labels, and shortening the labels does not auto-adjust these.
    */
}

// If tick label is hovered over, give full label in designated hoverspace
/*  !!! Currently event.target.closest does not work
document.addEventListener("mouseover", (event) => {
    // mouse enter
    if (!event.target.closest('.xaxislayer-above a')) {
        return;
    }
    const h5adContainer = event.target.closest('.h5ad-container');
    const hoverarea = h5adContainer.children(".hoverarea");

    // Show the full label when hovering over the truncated label (assuming we are on a page where the element exists)
    if (hoverarea.length) {
        hoverarea.textContent = event.target.dataset.unformatted;
    }
});

document.addEventListener('mouseleave', (event) => {
    // mouse out
    if (!event.target.closest('.xaxislayer-above a')) {
        return;
    }
    const h5adContainer = event.target.closest('.h5ad-container');
    const hoverarea = h5adContainer.children(".hoverarea");

    if (hoverarea.length) {
        hoverarea.textContent = "";
    }
});
*/
