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
            { value: "d3", text: "D3 (10 colors)" },
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
            { value: "purp", text: "Purples" }, // Cannot use in dotplot
            { value: "reds", text: "Reds" },
            { value: "bluered", text: "Blue-Red" },
            { value: "dense", text: "Dense" },
            { value: "electric", text: "Electric" },
            { value: "ylgnbu", text: "Yellow-Green-Blue" },
            { value: "ylorrd", text: "Yellow-Orange-Red" },
        ],
    },
    {
        label: "Diverging Scales",
        continuous: true,
        options: [
            { value: "earth", text: "Earth" },
            { value: "piyg", text: "Pink-Green" },
            { value: "prgn", text: "Purple-Green" },
            { value: "rdbu", text: "Red-Blue" },
            { value: "rdylbu", text: "Red-Yellow-Blue" },
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

// Store information about the colorscales
// See lib/gear/mg_plotting.get_colorscale for how I got the colorscales
const paletteInformation = {
    "alphabet":[[0,"#AA0DFE"],[0.04,"#3283FE"],[0.08,"#85660D"],[0.12,"#782AB6"],[0.16,"#565656"],[0.2,"#1C8356"],[0.24,"#16FF32"],[0.28,"#F7E1A0"],[0.32,"#E2E2E2"],[0.36,"#1CBE4F"],[0.4,"#C4451C"],[0.44,"#DEA0FD"],[0.48,"#FE00FA"],[0.52,"#325A9B"],[0.56,"#FEAF16"],[0.6,"#F8A19F"],[0.64,"#90AD1C"],[0.68,"#F6222E"],[0.72,"#1CFFCE"],[0.76,"#2ED9FF"],[0.8,"#B10DA1"],[0.84,"#C075A6"],[0.88,"#FC1CBF"],[0.92,"#B00068"],[0.96,"#FBE426"],[1,"#FA0087"]],
    "bold":[[0,"rgb(127, 60, 141)"],[0.1,"rgb(17, 165, 121)"],[0.2,"rgb(57, 105, 172)"],[0.3,"rgb(242, 183, 1)"],[0.4,"rgb(231, 63, 116)"],[0.5,"rgb(128, 186, 90)"],[0.6,"rgb(230, 131, 16)"],[0.7,"rgb(0, 134, 149)"],[0.8,"rgb(207, 28, 144)"],[0.9,"rgb(249, 123, 114)"],[1,"rgb(165, 170, 153)"]]
    ,"d3":[[0,"#1F77B4"],[0.1111111111111111,"#FF7F0E"],[0.2222222222222222,"#2CA02C"],[0.3333333333333333,"#D62728"],[0.4444444444444444,"#9467BD"],[0.5555555555555556,"#8C564B"],[0.6666666666666666,"#E377C2"],[0.7777777777777778,"#7F7F7F"],[0.8888888888888888,"#BCBD22"],[1,"#17BECF"]]
    ,"dark24":[[0,"#2E91E5"],[0.043478260869565216,"#E15F99"],[0.08695652173913043,"#1CA71C"],[0.13043478260869565,"#FB0D0D"],[0.17391304347826086,"#DA16FF"],[0.21739130434782608,"#222A2A"],[0.2608695652173913,"#B68100"],[0.30434782608695654,"#750D86"],[0.34782608695652173,"#EB663B"],[0.391304347826087,"#511CFB"],[0.43478260869565216,"#00A08B"],[0.4782608695652174,"#FB00D1"],[0.5217391304347826,"#FC0080"],[0.5652173913043478,"#B2828D"],[0.6086956521739131,"#6C7C32"],[0.6521739130434783,"#778AAE"],[0.6956521739130435,"#862A16"],[0.7391304347826086,"#A777F1"],[0.782608695652174,"#620042"],[0.8260869565217391,"#1616A7"],[0.8695652173913043,"#DA60CA"],[0.9130434782608695,"#6C4516"],[0.9565217391304348,"#0D2A63"],[1,"#AF0038"]]
    ,"light24":[[0,"#FD3216"],[0.043478260869565216,"#00FE35"],[0.08695652173913043,"#6A76FC"],[0.13043478260869565,"#FED4C4"],[0.17391304347826086,"#FE00CE"],[0.21739130434782608,"#0DF9FF"],[0.2608695652173913,"#F6F926"],[0.30434782608695654,"#FF9616"],[0.34782608695652173,"#479B55"],[0.391304347826087,"#EEA6FB"],[0.43478260869565216,"#DC587D"],[0.4782608695652174,"#D626FF"],[0.5217391304347826,"#6E899C"],[0.5652173913043478,"#00B5F7"],[0.6086956521739131,"#B68E00"],[0.6521739130434783,"#C9FBE5"],[0.6956521739130435,"#FF0092"],[0.7391304347826086,"#22FFA7"],[0.782608695652174,"#E3EE9E"],[0.8260869565217391,"#86CE00"],[0.8695652173913043,"#BC7196"],[0.9130434782608695,"#7E7DCD"],[0.9565217391304348,"#FC6955"],[1,"#E48F72"]]
    ,"safe":[[0,"rgb(136, 204, 238)"],[0.1,"rgb(204, 102, 119)"],[0.2,"rgb(221, 204, 119)"],[0.3,"rgb(17, 119, 51)"],[0.4,"rgb(51, 34, 136)"],[0.5,"rgb(170, 68, 153)"],[0.6,"rgb(68, 170, 153)"],[0.7,"rgb(153, 153, 51)"],[0.8,"rgb(136, 34, 85)"],[0.9,"rgb(102, 17, 0)"],[1,"rgb(136, 136, 136)"]]
    ,"vivid":[[0,"rgb(229, 134, 6)"],[0.1,"rgb(93, 105, 177)"],[0.2,"rgb(82, 188, 163)"],[0.3,"rgb(153, 201, 69)"],[0.4,"rgb(204, 97, 176)"],[0.5,"rgb(36, 121, 108)"],[0.6,"rgb(218, 165, 27)"],[0.7,"rgb(47, 138, 196)"],[0.8,"rgb(118, 78, 159)"],[0.9,"rgb(237, 100, 90)"],[1,"rgb(165, 170, 153)"]]
    ,"greys":[[0,"rgb(255,255,255)"],[0.125,"rgb(240,240,240)"],[0.25,"rgb(217,217,217)"],[0.375,"rgb(189,189,189)"],[0.5,"rgb(150,150,150)"],[0.625,"rgb(115,115,115)"],[0.75,"rgb(82,82,82)"],[0.875,"rgb(37,37,37)"],[1,"rgb(0,0,0)"]]
    ,"blues":[[0,"rgb(247,251,255)"],[0.125,"rgb(222,235,247)"],[0.25,"rgb(198,219,239)"],[0.375,"rgb(158,202,225)"],[0.5,"rgb(107,174,214)"],[0.625,"rgb(66,146,198)"],[0.75,"rgb(33,113,181)"],[0.875,"rgb(8,81,156)"],[1,"rgb(8,48,107)"]]
    ,"purp":[[0,"rgb(243, 224, 247)"],[0.16666666666666666,"rgb(228, 199, 241)"],[0.3333333333333333,"rgb(209, 175, 232)"],[0.5,"rgb(185, 152, 221)"],[0.6666666666666666,"rgb(159, 130, 206)"],[0.8333333333333333,"rgb(130, 109, 186)"],[1,"rgb(99, 88, 159)"]]
    ,"reds":[[0,"rgb(255,245,240)"],[0.125,"rgb(254,224,210)"],[0.25,"rgb(252,187,161)"],[0.375,"rgb(252,146,114)"],[0.5,"rgb(251,106,74)"],[0.625,"rgb(239,59,44)"],[0.75,"rgb(203,24,29)"],[0.875,"rgb(165,15,21)"],[1,"rgb(103,0,13)"]]
    ,"bluered":[[0,"rgb(0,0,255)"],[1,"rgb(255,0,0)"]]
    ,"dense":[[0,"rgb(230, 240, 240)"],[0.09090909090909091,"rgb(191, 221, 229)"],[0.18181818181818182,"rgb(156, 201, 226)"],[0.2727272727272727,"rgb(129, 180, 227)"],[0.36363636363636365,"rgb(115, 154, 228)"],[0.4545454545454546,"rgb(117, 127, 221)"],[0.5454545454545454,"rgb(120, 100, 202)"],[0.6363636363636364,"rgb(119, 74, 175)"],[0.7272727272727273,"rgb(113, 50, 141)"],[0.8181818181818182,"rgb(100, 31, 104)"],[0.9090909090909092,"rgb(80, 20, 66)"],[1,"rgb(54, 14, 36)"]]
    ,"electric":[[0,"rgb(0,0,0)"],[0.2,"rgb(30,0,100)"],[0.4,"rgb(120,0,100)"],[0.6000000000000001,"rgb(160,90,0)"],[0.8,"rgb(230,200,0)"],[1,"rgb(255,250,220)"]]
    ,"ylgnbu":[[0,"rgb(255,255,217)"],[0.125,"rgb(237,248,177)"],[0.25,"rgb(199,233,180)"],[0.375,"rgb(127,205,187)"],[0.5,"rgb(65,182,196)"],[0.625,"rgb(29,145,192)"],[0.75,"rgb(34,94,168)"],[0.875,"rgb(37,52,148)"],[1,"rgb(8,29,88)"]]
    ,"ylorrd":[[0,"rgb(255,255,204)"],[0.125,"rgb(255,237,160)"],[0.25,"rgb(254,217,118)"],[0.375,"rgb(254,178,76)"],[0.5,"rgb(253,141,60)"],[0.625,"rgb(252,78,42)"],[0.75,"rgb(227,26,28)"],[0.875,"rgb(189,0,38)"],[1,"rgb(128,0,38)"]]
    ,"earth":[[0,"rgb(161, 105, 40)"],[0.16666666666666666,"rgb(189, 146, 90)"],[0.3333333333333333,"rgb(214, 189, 141)"],[0.5,"rgb(237, 234, 194)"],[0.6666666666666666,"rgb(181, 200, 184)"],[0.8333333333333333,"rgb(121, 167, 172)"],[1,"rgb(40, 135, 161)"]]
    ,"piyg":[[0,"rgb(142,1,82)"],[0.1,"rgb(197,27,125)"],[0.2,"rgb(222,119,174)"],[0.30000000000000004,"rgb(241,182,218)"],[0.4,"rgb(253,224,239)"],[0.5,"rgb(247,247,247)"],[0.6000000000000001,"rgb(230,245,208)"],[0.7000000000000001,"rgb(184,225,134)"],[0.8,"rgb(127,188,65)"],[0.9,"rgb(77,146,33)"],[1,"rgb(39,100,25)"]]
    ,"prgn":[[0,"rgb(64,0,75)"],[0.1,"rgb(118,42,131)"],[0.2,"rgb(153,112,171)"],[0.30000000000000004,"rgb(194,165,207)"],[0.4,"rgb(231,212,232)"],[0.5,"rgb(247,247,247)"],[0.6000000000000001,"rgb(217,240,211)"],[0.7000000000000001,"rgb(166,219,160)"],[0.8,"rgb(90,174,97)"],[0.9,"rgb(27,120,55)"],[1,"rgb(0,68,27)"]]
    ,"rdbu":[[0,"rgb(103,0,31)"],[0.1,"rgb(178,24,43)"],[0.2,"rgb(214,96,77)"],[0.30000000000000004,"rgb(244,165,130)"],[0.4,"rgb(253,219,199)"],[0.5,"rgb(247,247,247)"],[0.6000000000000001,"rgb(209,229,240)"],[0.7000000000000001,"rgb(146,197,222)"],[0.8,"rgb(67,147,195)"],[0.9,"rgb(33,102,172)"],[1,"rgb(5,48,97)"]]
    ,"rdylbu":[[0,"rgb(165,0,38)"],[0.1,"rgb(215,48,39)"],[0.2,"rgb(244,109,67)"],[0.30000000000000004,"rgb(253,174,97)"],[0.4,"rgb(254,224,144)"],[0.5,"rgb(255,255,191)"],[0.6000000000000001,"rgb(224,243,248)"],[0.7000000000000001,"rgb(171,217,233)"],[0.8,"rgb(116,173,209)"],[0.9,"rgb(69,117,180)"],[1,"rgb(49,54,149)"]]
    ,"cividis":[[0,"#00224e"],[0.1111111111111111,"#123570"],[0.2222222222222222,"#3b496c"],[0.3333333333333333,"#575d6d"],[0.4444444444444444,"#707173"],[0.5555555555555556,"#8a8678"],[0.6666666666666666,"#a59c74"],[0.7777777777777777,"#c3b369"],[0.8888888888888888,"#e1cc55"],[1,"#fee838"]]
    ,"inferno":[[0,"#000004"],[0.1111111111111111,"#1b0c41"],[0.2222222222222222,"#4a0c6b"],[0.3333333333333333,"#781c6d"],[0.4444444444444444,"#a52c60"],[0.5555555555555556,"#cf4446"],[0.6666666666666666,"#ed6925"],[0.7777777777777777,"#fb9b06"],[0.8888888888888888,"#f7d13d"],[1,"#fcffa4"]]
    ,"viridis":[[0,"#440154"],[0.1111111111111111,"#482878"],[0.2222222222222222,"#3e4989"],[0.3333333333333333,"#31688e"],[0.4444444444444444,"#26828e"],[0.5555555555555556,"#1f9e89"],[0.6666666666666666,"#35b779"],[0.7777777777777777,"#6ece58"],[0.8888888888888888,"#b5de2b"],[1,"#fde725"]]
}


const postPlotlyConfig = {
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
        }, {
            plot_type: "heatmap"
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
function truncateAxisLabels() {
    const selector = document.querySelectorAll('.xaxislayer-above text');
    for (const el of selector) {
        const elText = el.textContent;
        const sublabel = elText.length > TICK_LABEL_MAX_LEN_ALLOWED
            ? `${elText.substring(0, TICK_LABEL_TRUNCATION_LEN)}...`
            : elText;

        el.innerHTML(`<a style="fill:inherit;">${sublabel}</a>`);
    }
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