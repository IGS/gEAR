
var gene_chr = "chr17";
var gene_start = 63921044;
var gene_end = 63951044;
// var id = "ENCFF995ODT";
// var datasource = "ENCFF995ODT";
// configuration is an array of tracks
// chartType can be "LineTrack", "ScatterPlot", "HeatmapPlot" etc
//  {"id": "", "chartType": "LineTrack", "name": ""}
var configuration = [{ "id": "ENCFF995ODT", "chartType": "LineTrack", "name": "ENCFF995ODT" }];

// If the gEAR frame has no epiviz-environment tag,
// This will add the environment and then add tracks
function renderEpivizChart(target, configuration) {
    var eNav = document.querySelector("#" + target).querySelector("epiviz-environment");
    if (!eNav) {
        eNav = document.createElement("epiviz-environment");
        elem.setAttribute("chr", gene_chr);
        elem.setAttribute("start", gene_start);
        elem.setAttribute("end", gene_end);
        elem.setAttribute("no-logo", true);
        document.querySelector("#" + target).appendChild(eNav);
    }

    configuration.forEach(function (track) {
        var cTag = get_chart_by_type(track.chartType);
        var dataManagerElem = document.querySelector('epiviz-data-source');
        var data = dataManagerElem.measurementSet.subset(function (m) { return m.datasourceGroup() == track.id; });
        var measurements = data.subset(function (m) { return m.datasourceGroup() === track.name });
        var dataTrack = document.createElement(cTag);
        dataTrack.slot = "charts";
        dataTrack.setAttribute("measurements", JSON.stringify(measurements.raw()));
        Polymer.dom(eNav).appendChild(dataTrack);
    });
}

// all possible charts
function get_chart_by_type(type) {
    var charts = {
        "GenesTrack": "epiviz-genes-track",
        "BlocksTrack": 'epiviz-blocks-track',
        "StackedBlocksTrack": 'epiviz-stacked-blocks-track',
        "HeatmapPlot": "epiviz-heatmap-plot",
        "ScatterPlot": "epiviz-scatter-plot",
        "LineTrack": 'epiviz-line-track',
        "StackedLineTrack": 'epiviz-stacked-line-track',
        "LinePlot": 'epiviz-line-plot',
        "StackedLinePlot": 'epiviz-stacked-line-plot'
    }

    return (charts[type]);
}