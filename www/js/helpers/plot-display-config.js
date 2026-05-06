/*
plot-display-config.js - This script can be passed to various scripts that generate plots
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
export const availablePalettes = [
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
export const plotly2MatplotlibNames = {
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

export const postPlotlyConfig = {
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
                },
                // modern tweaks
                font: {
                        family: 'Inter, Roboto, Arial, sans-serif',
                        size: 12,
                        color: '#333333'
                },
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

export const adjustClusterColorbars = (plotData) => {
    return;   // Disable for now since we are no longer doing clusterbars as separate traces.  Will need to re-evalulate if we want to go back to that approach.
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

export const setHeatmapHeightBasedOnGenes = (plotLayout, genesFilter) => {
    if (genesFilter.length > 20) {
        plotLayout.height = genesFilter.length * 20;
    }
}

export const adjustStackedViolinHeight = (plotLayout) => {
    plotLayout.height = 800;
}

// Scaling a number from one range to another range (Source: https://stackoverflow.com/a/31687097)
const scaleBetween = (unscaledNum, minAllowed, maxAllowed, min, max) => {
    return (maxAllowed - minAllowed) * (unscaledNum - min) / (max - min) + minAllowed;
}

/**
 * Attaches hover tooltips to truncated Plotly axis tick labels using the
 * axis_label_mapping stored in the plot's layout metadata.
 *
 * @param {string} plotDivId - The ID of the Plotly plot div.
 */
export const attachAxisLabelTooltips = (plotDivId) => {
    const plotDiv = document.getElementById(plotDivId);
    if (!plotDiv) return;

    const layout = plotDiv.layout;
    const axisLabelMapping = layout?.meta?.axis_label_mapping;
    if (!axisLabelMapping || !Object.keys(axisLabelMapping).length) return;

    // Build a map of tick label text -> full name for quick lookup
    const tickLabelMap = new Map();
    const tickLabels = plotDiv.querySelectorAll('.xaxislayer-above .xtick text, .xaxislayer .xtick text');
    for (const tickLabel of tickLabels) {
        const labelText = tickLabel.textContent?.trim();
        if (labelText && axisLabelMapping[labelText]) {
            // map element's text property to full name for tooltip lookup later
            tickLabelMap.set(tickLabel, axisLabelMapping[labelText]);
        }
    }

    if (!tickLabelMap.size) return;

    // Create a shared tooltip element
    const tooltip = document.createElement('div');
    tooltip.classList.add('tooltip');   // Bulma's style
    tooltip.style.position = 'absolute';
    tooltip.style.pointerEvents = 'none'; // Ensure the tooltip doesn't "intercept" mouse events
    tooltip.style.fontSize = "12px";
    tooltip.style.bottom = "0px";
    tooltip.style.left = "0px";
    tooltip.style.backgroundColor = 'white';
    tooltip.style.opacity = 0.8;
    tooltip.style.color = 'black';
    tooltip.style.padding = '5px';
    tooltip.style.border = '1px solid gray';
    tooltip.style.zIndex = 3;
    tooltip.style.display = 'none'; // Initial state
    plotDiv.appendChild(tooltip);

    /**
     * Check if a mouse event's coordinates fall within a tick label's bounding rect.
     * Using getBoundingClientRect() works for both HTML and SVG elements and returns
     * viewport-relative coordinates, matching clientX/clientY directly.
     */
    const findHoveredTick = (clientX, clientY) => {
        for (const [tickEl, fullName] of tickLabelMap) {
            const rect = tickEl.getBoundingClientRect();
            if (
                clientX >= rect.left &&
                clientX <= rect.right &&
                clientY >= rect.top &&
                clientY <= rect.bottom
            ) {
                return fullName;
            }
        }
        return null;
    };

    // Plotly renders a transparent drag layer on top that blocks mouse events to SVG elements.
    // Attach to the plotDiv itself and use bounding rect comparison to detect tick proximity.
    plotDiv.addEventListener('mousemove', (event) => {
        const fullName = findHoveredTick(event.clientX, event.clientY);
        if (fullName) {
            tooltip.innerHTML = `<strong>${fullName}</strong>`;
            tooltip.style.display = 'block';
        } else {
            tooltip.style.display = 'none';
        }
    });

    plotDiv.addEventListener('mouseleave', () => {
        tooltip.style.display = 'none';
    });

 // Clean up tooltip from body when plot div is removed from DOM
    const observer = new MutationObserver(() => {
        if (!document.body.contains(plotDiv)) {
            tooltip.remove();
            observer.disconnect();
        }
    });
    observer.observe(document.body, { childList: true, subtree: true });
};

