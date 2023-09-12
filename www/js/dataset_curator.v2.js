// I use camelCase for my variable/function names to adhere to JS style standards
// Exception being functions that do fetch calls, so we can use JS destructuring on the payload

'use strict';

//TODO  - Create "error" reporting in each step

let plotStyle;  // Plot style object

//let plotConfig = {};  // Plot config that is passed to API or stored in DB
let allColumns = [];
let catColumns = [];
let levels = {};    // categorical columns + groups

let datasetId = null;
let organismId = null;
let analysisObj = null;

let analysisSelect = null;
let plotTypeSelect = null;
let geneSelect = null;

const userId = 622;     // ! It's me
const sessionId = "ee95e48d-c512-4083-bf05-ca9f65e2c12a"    // ! It's me
const session_id = sessionId;
const colorblindMode = false;

Cookies.set('gear_session_id', sessionId, { expires: 7 });

const plotlyPlots = ["bar", "line", "scatter", "tsne_dyna", "violin"];
const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];

/*
! Quick note -
This code uses a lot of "for...in..." and "for...of..." syntax
"for...in..." is generally used to iterate over object keys or array indexes
"for...of..." is generally used to iterate over array (or collection) values, and will not work on objects

Using "for...of..." over array.forEach means we don't have to transform HTMLCollections into arrays
*/

// Create PlotStyle abstract class and Plotly, Scanpy, SVG subclasses
class PlotHandler {
    constructor() {
        // Check if this is an abstract class
        // (new PlotStyle() will fail)
        if (new.target === PlotHandler) throw new TypeError("Cannot construct PlotHandler instances directly");
    }

    // Certain groups like "display order" and "colors" and "vlines" will be custom-handled
    classElt2Prop = {}; // This will be overridden by subclasses
    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    cloneDisplay() {
        throw new Error("You have to implement the method cloneDisplay!");
    }

    async createPlot() {
        throw new Error("You have to implement the method createPlot!");
    }

    async loadPlotHtml() {
        throw new Error("You have to implement the method loadPlotHtml!");
    }

    populatePlotConfig() {
        throw new Error("You have to implement the method populatePlotConfig!");
    }

    setupPlotSpecificEvents() {
        throw new Error("You have to implement the method setupPlotSpecificEvents!");
    }

}

// Class for methods that are common to all Plotly plot types
class PlotlyHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
    }

    // TODO: Add color palette for continuous colors

    classElt2Prop = {
        "js-plotly-x-axis":"x_axis"
        , "js-plotly-y-axis":"y_axis"
        , "js-plotly-label":"point_label"
        , "js-plotly-hide-x-ticks":"hide_x_labels"
        , "js-plotly-hide-y-ticks":"hide_y_labels"
        , "js-plotly-color":"color_name"
        , "js-plotly-size":"size_by_group"
        , "js-plotly-facet-row":"facet_row"
        , "js-plotly-facet-col":"facet_col"
        , "js-plotly-x-title":"x_title"
        , "js-plotly-y-title":"y_title"
        , "js-plotly-x-min":"x_min"
        , "js-plotly-y-min":"y_min"
        , "js-plotly-x-max":"x_max"
        , "js-plotly-y-max":"y_max"
        , "js-plotly-hide-legend":"hide_legend"
        , "js-plotly-add-jitter":"jitter"
        , "js-plotly-marker-size":"marker_size"
    }

    configProp2ClassElt = super.configProp2ClassElt;

    plotConfig = {};  // Plot config that is passed to API

    cloneDisplay(config) {
        // plotly plots
        for (const prop in this.configProp2ClassElt) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            const orderContainer = document.getElementById("order_container");
            for (const series of config["order"]) {
                renderOrderSortableSeries(series);
            }

            document.getElementById("order_section").style.display = "";
        }

        // Handle colors
        if (config["colors"]) {
            // do nothing if color_name is not set
            if (!config["color_name"]) return;
            const series = config["color_name"];
            renderColorPicker(series);
            for (const group in config["colors"]) {
                const color = config["colors"][group];
                const colorField = document.getElementById(`${group}_color`);
                colorField.value = color;
            }
        }

        // Handle vlines
        if (config["vlines"]) {
            const vLinesBody = document.getElementById("vlines_body");
            const vlineField = document.querySelector(".js-plotly-vline-field");
            // For each vline object create and populate a new vline field
            for (const vlineObj in config["vlines"]) {
                const newVlineField = vlineField.cloneNode(true);
                newVlineField.querySelector(":scope .js-vline-pos").value = vlineObj["vl_pos"];
                newVlineField.querySelector(":scope .js-vline-style-select select").value = vlineObj["vl_style"];
                vLinesBody.prepend(newVlineField);
            }
        }
    }

    async createPlot(datasetId, analysisObj, userId, colorblindMode) {
        // Get data and set up the image area
        let plotJson;
        try {
            const data = await fetchPlotlyData(this.plotConfig, datasetId, this.plotType, analysisObj, userId, colorblindMode);
            ({plot_json: plotJson} = data);
        } catch (error) {
            // TODO: Make failure indicator
            return;
        }

        const plotContainer = document.getElementById("plot_container");
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const divElt = generateElements('<div class="container is-max-desktop" id="plotly_preview"></div>');
        plotContainer.append(divElt);
        Plotly.purge("plotly_preview"); // clear old Plotly plots

        // TODO: Throw error if "plotJson" is null
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "config");
        Plotly.newPlot("plotly_preview", plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "layout")
        Plotly.relayout("plotly_preview", custonLayout)

    }

    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot_options_collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post_plot_adjustments");
        postPlotOptionsElt.replaceChildren();

        const prePlot = "../include/plot_config/pre_plot/single_gene_plotly.html";
        const preResponse = await fetch(prePlot, {cache: "reload"});
        const preBody = await preResponse.text();
        prePlotOptionsElt.innerHTML = preBody;

        const postPlot = "../include/plot_config/post_plot/single_gene_plotly.html"
        const postResponse = await fetch(postPlot, {cache: "reload"});
        const postBody = await postResponse.text();
        postPlotOptionsElt.innerHTML = postBody;

        // populate advanced options for specific plot types
        const prePlotSpecificOptionsElt = document.getElementById("plot_specific_options");
        const postPlotSpecificOptionselt = document.getElementById("post_plot_specific_options");

        if (["scatter", "tsne_dyna"].includes(this.plotType)) {
            const prePlot = "../include/plot_config/pre_plot/advanced_scatter.html";
            const preResponse = await fetch(prePlot, {cache: "reload"});
            const preBody = await preResponse.text();
            prePlotSpecificOptionsElt.innerHTML = preBody;

            const postPlot = "../include/plot_config/post_plot/advanced_scatter.html";
            const postResponse = await fetch(postPlot, {cache: "reload"});
            const postBody = await postResponse.text();
            postPlotSpecificOptionselt.innerHTML = postBody;
            return;
        }
        if (this.plotType === "violin") {
            const prePlot = "../include/plot_config/pre_plot/advanced_violin.html";
            const preResponse = await fetch(prePlot, {cache: "reload"});
            const preBody = await preResponse.text();
            prePlotSpecificOptionsElt.innerHTML = preBody;

            const postPlot = "../include/plot_config/post_plot/advanced_violin.html";
            const postResponse = await fetch(postPlot, {cache: "reload"});
            const postBody = await postResponse.text();
            postPlotSpecificOptionselt.innerHTML = postBody;
            return;
        }
    }

    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = getPlotConfigValueFromClassName(classElt)
        }

        // Small fix for tsne/umap dynamic plots
        if (this.plotType.toLowerCase() === "tsne_dyna") {
            this.plotType = "tsne/umap_dynamic";
        }

        // Get order
        this.plotConfig["order"] = getPlotOrderFromSortable();

        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("color_series_post").value;
        if (colorSeries && colorElts.length) {
            // Input is either color mapping or just the series
            this.plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("_color", "");
                this.plotConfig["colors"][group] = field.value;
            })
        }

        // Get vlines
        const vlineFields = document.getElementsByClassName("js-plotly-vline-field");
        this.plotConfig["vlines"] = [...vlineFields].map((field) => {
            const vlinePos = field.querySelector(":scope .js-plotly-vline-pos").value;
            const vlineStyle = field.querySelector(":scope .js-plotly-vline-style-select select").value;
            // Return either objects or nothing (which will be filtered out)
            return vlinePos ?  {"vl_pos":vlinePos, "vl_style":vlineStyle} : null;
        }).filter(x => x !== null);
    }

    setupPlotSpecificEvents() {
        setupPlotlyOptions();
    }

}

// Class for methods that are common to all Scanpy plot types
class ScanpyHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
    }

    classElt2Prop = {
        "js-tsne-x-axis":"x_axis"
        , "js-tsne-y-axis":"y_axis"
        , "js-tsne-colorize-legend-by":"colorize_legend_by"
        , "js-tsne-plot-by-series":"plot_by_group"
        , "js-tsne-max-columns":"max_columns"
        , "js-tsne-skip-gene-plot":"skip_gene_plot"
        , "js-tsne-horizontal-legend":"horizontal_legend"
    }

    configProp2ClassElt = super.configProp2ClassElt;

    plotConfig = {};  // Plot config that is passed to API

    cloneDisplay(config) {
        for (const prop in this.configProp2ClassElt) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            const orderContainer = document.getElementById("order_container");
            for (const series of config["order"]) {
                renderOrderSortableSeries(series);
            }

            document.getElementById("order_section").style.display = "";
        }

        // Restoring some disabled/checked elements in UI
        const plotBySeries = document.getElementsByClassName("js-tsne-plot-by-series");
        const maxColumns = document.getElementsByClassName('js-tsne-max-columns');
        const skipGenePlot = document.getElementsByClassName("js-tsne-skip-gene-plot");
        const horizontalLegend = document.getElementsByClassName("js-tsne-horizontal-legend");

        if (config["colorize_legend_by"]) {
            const series = config["colorize_legend_by"];
            for (const targetElt of [...plotBySeries, ...maxColumns, ...horizontalLegend]) {
                targetElt.disabled = (catColumns.includes(series)) ? false : true;
            }

            // Handle colors
            if (config["colors"]) {
                renderColorPicker(series);
                for (const group in config["colors"]) {
                    const color = config["colors"][group];
                    const colorField = document.getElementById(`${group}_color`);
                    colorField.value = color;
                }
            }
        }

        if (config["plot_by_group"]) {
            for (const targetElt of skipGenePlot) {
                targetElt.disabled = true;
                targetElt.checked = false;
            }
            for (const targetElt of maxColumns) {
                targetElt.disabled = false;
            }
        }
    }

    async createPlot(datasetId, analysisObj, userId, colorblindMode) {
        let image;
        try {
            const data = await fetchTsneImageData(this.plotConfig, datasetId, this.plotType, analysisObj, userId, colorblindMode);
            ({image} = data);
        } catch (error) {
            // TODO: Make failure indicator
            return;
        }

        const plotContainer = document.getElementById("plot_container");
        plotContainer.replaceChildren();    // erase plot

        const imgElt = generateElements('<img class="image" id="tsne_preview"></img>');
        plotContainer.append(imgElt);

        if (image) {
            document.getElementById("tsne_preview").setAttribute("src", `data:image/png;base64,${image}`);
        }
    }

    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot_options_collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post_plot_adjustments");
        postPlotOptionsElt.replaceChildren();

        const prePlot = "../include/plot_config/pre_plot/tsne_static.html";
        const preResponse = await fetch(prePlot, {cache: "reload"});
        const preBody = await preResponse.text();
        prePlotOptionsElt.innerHTML = preBody;

        const postPlot = "../include/plot_config/post_plot/tsne_static.html"
        const postResponse = await fetch(postPlot, {cache: "reload"});
        const postBody = await postResponse.text();
        postPlotOptionsElt.innerHTML = postBody;
    }

    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = getPlotConfigValueFromClassName(classElt)
        }

        // Get order
        this.plotConfig["order"] = getPlotOrderFromSortable();

        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("colorize_legend_by_post").textContent;
        if (colorSeries && colorElts.length) {
            this.plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("_color", "");
                this.plotConfig["colors"][group] = field.value;
            })
        }

        // If user did not want to have a colorized annotation, ensure it does not get passed to the scanpy code
        if (!(colorSeries)) {
            this.plotConfig["plot_by_group"] = null;
            this.plotConfig["max_columns"] = null;
            this.plotConfig["skip_gene_plot"] = false;
            this.plotConfig["horizontal_legend"] = false;
        }
    }

    setupPlotSpecificEvents() {
        setupScanpyOptions();
    }

}

// Class for methods that are common to all SVG images
class SvgHandler extends PlotHandler {
    constructor() {
        super();
    }

    // These do not get passed into the API call, but want to keep the same data structure for cloning display
    classElt2Prop = {
        "js-svg-low-color":"low_color"
        , "js-svg-mid-color":"mid_color"
        , "js-svg-high-color":"high_color"
    }

    configProp2ClassElt = super.configProp2ClassElt;

    plotConfig = {colors: {}};  // Plot config to color SVG

    cloneDisplay(config) {
        // Props are in a "colors" dict
        for (const prop in this.configProp2ClassElt) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config.colors[prop]);
        }

        // If a mid-level color was provided, ensure checkbox to enable it is checked (for aesthetics)
        if (config.colors["mid_color"]) {
            for (const elt of document.getElementsByClassName("js-svg-enable-mid")) {
                elt.checked = true;
            }
        }

    }

    async createPlot(datasetId, geneSymbol) {
        let data;
        try {
            data = await fetchSvgData(datasetId, geneSymbol)
        } catch (error) {
            // TODO: Make failure indicator
            return;
        }
        colorSVG(data, this.plotConfig["colors"]);
    }

    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot_options_collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post_plot_adjustments");
        postPlotOptionsElt.replaceChildren();

        const prePlot = "../include/plot_config/pre_plot/svg.html";
        const preResponse = await fetch(prePlot, {cache: "reload"});
        const preBody = await preResponse.text();
        prePlotOptionsElt.innerHTML = preBody;

        const postPlot = "../include/plot_config/post_plot/svg.html"
        const postResponse = await fetch(postPlot, {cache: "reload"});
        const postBody = await postResponse.text();
        postPlotOptionsElt.innerHTML = postBody;
    }

    populatePlotConfig() {
        this.plotConfig["colors"] = {};   // Reset plot config

        this.plotConfig["colors"]["low_color"] = document.getElementById("low_color").value;
        this.plotConfig["colors"]["mid_color"] = document.getElementById("mid_color").value;
        this.plotConfig["colors"]["high_color"] = document.getElementById("high_color").value;

        // If user did not choose a mid-color, set it as null instead of to black
        if (!(document.getElementById("enable_mid_color").checked)) {
            this.plotConfig["colors"]["mid_color"] = null;
        }
    }

    setupPlotSpecificEvents() {
        setupSVGOptions();
    }

}

const datasetTree = new DatasetTree({
    element: document.getElementById("dataset_tree")
    , searchElement: document.getElementById("dataset_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current_dataset_c").style.display = "";
        document.getElementById("current_dataset").textContent = e.node.title;
        document.getElementById("current_dataset_post").textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        // Click to get to next step
        document.getElementById("load_plot_s").click();
        document.getElementById('new_display').classList.add("is-loading");

        // Clear "current <whatever>" text
        document.getElementById("current_gene_c").style.display = "none";
        document.getElementById("current_analysis_c").style.display = "none";
        document.getElementById("current_plot_type_c").style.display = "none";

        // Clear "success/failure" icons
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.style.display = "none";
        }

        datasetId = newDatasetId;

        document.getElementById("dataset_s_success").style.display = "";

        // Clear (and update) options within nice-select2 structure.
        analysisObj = null;
        analysisSelect.clear(); // BUG: Figure out why this is triggering gene-population twice (check function in Github)
        geneSelect.clear();
        plotTypeSelect.clear();

        // Fetch dataset information
        let ownerId;
        try {
            // Must wrap in parentheses - https://stackoverflow.com/a/48714713
            ({owner_id: ownerId} = await fetchDatasetInfo(datasetId));
        } catch (error) {
            ownerId = -1;   // Owner displays wont be fetched regardless
        }

        // displays
        const userDisplays = await fetchDatasetDisplays(userId, datasetId);
        const ownerDisplays = userId === ownerId ? [] : await fetchDatasetDisplays(ownerId, datasetId);
        let defaultDisplayId;
        try {
            defaultDisplayId = await fetchDefaultDisplay(userId, datasetId);
        } catch (error) {
            defaultDisplayId = -1;  // Cannot make any display a default.
        }
        renderDisplayCards(userDisplays, ownerDisplays, defaultDisplayId);
        document.getElementById('new_display').classList.remove("is-loading");
        document.getElementById('new_display').disabled = false;

    })
});

const chooseAnalysis = async (event) => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const analysisText = analysisId || "Primary Analysis";

    // Display current selected analysis
    document.getElementById("current_analysis_c").style.display = "";
    document.getElementById("current_analysis").textContent = analysisText;
    document.getElementById("current_analysis_post").textContent = analysisText;


    // NOTE: For now, we can just pass analysis id only to tSNE and be fine
    // Any private dataset will belong to our user. Any public datasets can be found by the API "get_analysis" code.
    if (analysisId) {
        analysisObj = {id: analysisId};
    }

    if (analysisSelect.data.length < 5) return; // Have not retrieved analyses from API yet

    // Populate gene select element
    try {
        const geneSymbols = await fetchGeneSymbols(datasetId, analysisId);
        updateGeneSymbolOptions(geneSymbols);
    } catch (error) {
        document.getElementById("gene_s_failed").style.display = "";
    }

}

/* Display has been chosen, so display analysis and plot type options */
const chooseDisplay = async (event) => {
    document.getElementById('new_display').classList.add("is-loading");
    document.getElementById("analysis_select").disabled = false;

    document.getElementById("plot_type_select").disabled = false;


    // analyses
    try {
        const {publicAnalyses, privateAnalyses} = await fetchAnalyses(datasetId);
        updateAnalysesOptions(privateAnalyses, publicAnalyses);
        analysisSelect.update();
        document.getElementById("analysis_type_select_c_success").style.display = "";   // Default analysis is good
    } catch (error) {
        console.error(error);
        // Show failure state things.
        document.getElementById("plot_type_s_failed").style.display = "";
        document.getElementById("analysis_type_select_c_failed").style.display = "";
        document.getElementById('new_display').classList.remove("is-loading"); // Don't give impression display is still loading
        return;
    } finally {
        document.getElementById("load_plot_s_success").style.display = "";
    }

    // plot types
    // NOTE: Believe updating "disabled" properties triggers the plotSelect "change" element
    try {
        const availablePlotTypes = await fetchAvailablePlotTypes(userId, sessionId, datasetId, undefined);
        for (const plotType in availablePlotTypes) {
            const isAllowed = availablePlotTypes[plotType];
            setPlotTypeDisabledState(plotType, isAllowed);

        }
        plotTypeSelect.update();
    } catch (error) {
        console.error(error);
        document.getElementById("plot_type_s_failed").style.display = "";
        document.getElementById("plot_type_select_c_failed").style.display = "";
        document.getElementById("plot_type_s_success").style.display = "none";
        document.getElementById("plot_type_select_c_success").style.display = "none";

        return;
    } finally {
        document.getElementById('new_display').classList.remove("is-loading");
    }

    // Populate gene select element
    try {
        const geneSymbols = await fetchGeneSymbols(datasetId, null);
        updateGeneSymbolOptions(geneSymbols);
    } catch (error) {
        document.getElementById("gene_s_failed").style.display = "";
    }

    document.getElementById("plot_type_s").click();
}

const chooseGene = (event) => {
    if (!geneSelect.selectedOptions.length) return;   // Do not trigger after initial population

    document.getElementById("gene_s_success").style.display = "";

    // Display current selected gene
    document.getElementById("current_gene_c").style.display = "";
    document.getElementById("current_gene").textContent = getSelect2Value(geneSelect);
    document.getElementById("current_gene_post").textContent = getSelect2Value(geneSelect);

    document.getElementById("plot_options_s").click();
}

const choosePlotType = (event) => {
    if (!plotTypeSelect.selectedOptions.length) return;   // Do not trigger after setting disable/enable on options

    // Do not display if default opt is chosen
    const plotType = getSelect2Value(plotTypeSelect)
    if (plotType === "nope") {
        document.getElementById("plot_type_select_c_success").style.display = "none";
        document.getElementById("plot_type_s_success").style.display = "none";
        return;
    }

    document.getElementById("plot_type_select_c_failed").style.display = "none";
    document.getElementById("plot_type_select_c_success").style.display = "";

    document.getElementById("plot_type_s_failed").style.display = "none";
    document.getElementById("plot_type_s_success").style.display = "";

    // Display current selected plot type
    document.getElementById("current_plot_type_c").style.display = "";
    document.getElementById("current_plot_type").textContent = plotType;

    includePlotParamOptions();
    document.getElementById("gene_s").click();
}

const cloneDisplay = async (event, display) => {

    const cloneId = event.target.id;
    document.getElementById(cloneId).classList.add("is-loading");

    document.getElementById("analysis_select").disabled = false;
    document.getElementById("plot_type_select").disabled = false;


    // analyses
    try {
        const {publicAnalyses, privateAnalyses} = await fetchAnalyses(datasetId);
        updateAnalysesOptions(privateAnalyses, publicAnalyses);
        analysisSelect.update();
        document.getElementById("analysis_type_select_c_success").style.display = "";   // Default analysis is good
    } catch (error) {
        console.error(error);
        // Show failure state things.
        document.getElementById("plot_type_s_failed").style.display = "";
        document.getElementById("analysis_type_select_c_failed").style.display = "";
        document.getElementById('new_display').classList.remove("is-loading"); // Don't give impression display is still loading
        return;
    } finally {
        document.getElementById("load_plot_s_success").style.display = "";
    }


    // plot types

    // Read clone config to populate analysis, plot type, gnee and plot-specific options
    let plotType = display.plot_type;

    // ? Move this to class constructor to handle
    if (plotType.toLowerCase() === "tsne") {
        // Handle legacy plots
        plotType = "tsne_static";
    } else if (["tsne/umap_dynamic", "tsne_dynamic"].includes(plotType.toLowerCase())) {
        plotType = "tsne_dyna";
    }

    try {
        const availablePlotTypes = await fetchAvailablePlotTypes(userId, sessionId, datasetId, undefined);
        for (const plotType in availablePlotTypes) {
            const isAllowed = availablePlotTypes[plotType];
            setPlotTypeDisabledState(plotType, isAllowed);
        }

        setSelectBoxByValue("plot_type_select", plotType);
        plotTypeSelect.update();    // This should trigger fetching the right "plot" html files
        trigger(document.getElementById("plot_type_select"), "change");
    } catch (error) {
        console.error(error);
        document.getElementById("plot_type_s_failed").style.display = "";
        document.getElementById("plot_type_select_c_failed").style.display = "";
        document.getElementById("plot_type_s_success").style.display = "none";
        document.getElementById("plot_type_select_c_success").style.display = "none";

        return;
    } finally {
        document.getElementById(cloneId).classList.remove("is-loading");
    }

    // Populate gene select element
    try {
        const geneSymbols = await fetchGeneSymbols(datasetId, null);
        updateGeneSymbolOptions(geneSymbols);
    } catch (error) {
        document.getElementById("gene_s_failed").style.display = "";
    }

    const config = display.plotly_config;
    setSelectBoxByValue("gene_select", config.gene_symbol);
    geneSelect.update();
    trigger(document.getElementById("gene_select"), "change");


    if (plotlyPlots.includes(plotType)) {
        plotStyle = new PlotlyHandler(plotType);
    } else if (scanpyPlots.includes(plotType)) {
        plotStyle = new ScanpyHandler(plotType);
    } else if (plotType === "svg") {
        plotStyle = new SvgHandler();
    } else {
        console.warn(`Plot type ${plotType} not recognized.`)
        return;
    }

    plotStyle.cloneDisplay(config);
}

const colorSVG = (chartData, plotConfig) => {
    // I found adding the mid color for the colorblind mode  skews the whole scheme towards the high color
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : plotConfig["low_color"];
    const midColor = colorblindMode ? null : plotConfig["mid_color"];
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : plotConfig["high_color"];

    // for those fields which have no reading, a specific value is sometimes put in instead
    // These are colored a neutral color
    const NA_FIELD_PLACEHOLDER = -0.012345679104328156;
    const NA_FIELD_COLOR = '#808080';

    //const scoreMethod = document.getElementById("scoring_method").value;
    const score = chartData.scores["gene"]
    const { min, max } = score;
    let color = null;
    // are we doing a three- or two-color gradient?
    if (midColor) {
        if (min >= 0) {
            // All values greater than 0, do right side of three-color
            color = d3
                .scaleLinear()
                .domain([min, max])
                .range([midColor, highColor]);
        } else if (max <= 0) {
            // All values under 0, do left side of three-color
            color = d3
                .scaleLinear()
                .domain([min, max])
                .range([lowColor, midColor]);
        } else {
            // We have a good value range, do the three-color
            color = d3
                .scaleLinear()
                .domain([min, 0, max])
                .range([lowColor, midColor, highColor]);
        }
    } else {
        color = d3
            .scaleLinear()
            .domain([min, max])
            .range([lowColor, highColor]);
    }


    // Load SVG file and set up the window
    const svg = document.getElementById("plot_container");
    const snap = Snap(svg);
    const svg_path = `datasets_uploaded/${datasetId}.svg`;
    Snap.load(svg_path, async (path) => {
        await snap.append(path)

        snap.select("svg").attr({
            width: "100%",
        });

        // Fill in tissue classes with the expression colors
        const {data: expression} = chartData;
        const tissues = Object.keys(chartData.data);   // dataframe
        const paths = Snap.selectAll("path, circle");

        // NOTE: This must use the SnapSVG API Set.forEach function to iterate
        paths.forEach(path => {
            const tissue = path.node.className.baseVal;
            if (tissues.includes(tissue)) {
                if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                    path.attr('fill', NA_FIELD_COLOR);
                } else {
                    path.attr('fill', color(expression[tissue]));
                }
            }
        });

        // TODO: Potentially replicate some of the features in display.js like log-transforms and tooltips
    });

}

const createAnalysisSelectInstance = () => {
    analysisSelect = NiceSelect.bind(document.getElementById("analysis_select"), {
        placeholder: 'Select an analysis.',
        allowClear: true,
    });
}

const createGeneSelectInstance = () => {
    geneSelect = NiceSelect.bind(document.getElementById("gene_select"), {
        placeholder: 'To search, start typing a gene name',
        searchtext: 'To search, start typing a gene name',
        searchable: true,
        allowClear: true,
    });

    // BUG: Figure out "scrollbar"
}

const createPlot = async (event) => {

    const plotType = getSelect2Value(plotTypeSelect);

    const plotBtns = document.getElementsByClassName("js-plot-btn");

    // Set loading
    for (const plotBtn of plotBtns) {
        plotBtn.classList.add("is-loading");
    }

    // Delete unsupplied parameter configs
    /*
    for (const param in plotConfig) {
        if (!plotConfig[param]) {
            delete plotConfig[param];
            continue
        };
        if (["object", "array"].includes(typeof(plotConfig[param]))) {
            if (typeof(plotConfig[param]) === "object" && !Object.keys(plotConfig[param]).length) {
                delete plotConfig[param];
                continue
            }
            if (typeof(plotConfig[param]) === "array" && !(plotConfig[param].length)) {
                delete plotConfig[param];
                continue
            };
        }
    }
    */

    plotStyle.populatePlotConfig();

    const geneSymbol = getSelect2Value(geneSelect);
    plotStyle.plotConfig["gene_symbol"] = geneSymbol;

    // Call API route by plot type
    if (plotlyPlots.includes(plotType)) {
        plotStyle.createPlot(datasetId, analysisObj, userId, colorblindMode);

    } else if (scanpyPlots.includes(plotType)) {
        plotStyle.createPlot(datasetId, analysisObj, userId, colorblindMode);

    } else if (plotType === "svg") {
        plotStyle.createPlot(datasetId, geneSymbol);
    } else {
        console.warn(`Plot type ${plotType} selected for plotting is not a valid type.`)
        return;
    }

    // Stop loader
    for (const plotBtn of plotBtns) {
        plotBtn.classList.remove("is-loading");
    }

    // Hide this view
    document.getElementById("content_c").style.display = "none";
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").style.display = "";

    // TODO: Add alert for non-success w/ message

}

const createPlotSelectInstance = () => {
    // Initialize plot types
    plotTypeSelect = NiceSelect.bind(document.getElementById("plot_type_select"), {
        placeholder: 'Choose how to plot',
        minimumResultsForSearch: -1
    });
}

/* Creates a Toast-style message in the upper-corner of the screen. */
const createToast = (msg) => {
    const template = `
    <div class="notification is-danger animate__animated animate__fadeInUp animate__faster">
        <button class="delete"></button>
        ${msg}
    </div>
    `
    const html = generateElements(template);
    document.getElementById("main_c").prepend(html);

    // This should get the newly added notification since it is now the first
    document.querySelector(".notification .delete").addEventListener("click", (event) => {
        const notification = event.target.parentNode;
        notification.parentNode.removeChild(notification);
    });
}

const deleteDisplay = async(user_id, display_id) => {
    const payload = {user_id, display_id};
    try {
        await axios.post("/cgi/delete_dataset_display.cgi", convertToFormData(payload));
        // Remove display card
        const displayCard = document.getElementById(`${display_id}_display`);
        displayCard.remove();
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not delete this display. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchAnalyses = async (datasetId) => {
    try {
        const { data } = await axios.get(`./api/h5ad/${datasetId}/analyses`);

        const { public: publicAnalyses, private: privateAnalyses } = data;
        return {publicAnalyses, privateAnalyses};

    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch saved analyses for this dataset. Please contact the gEAR team."
        createToast(msg);
        //throw new Error(msg, {cause: error}); // ECMAscript 2022 nice-to-have
        throw new Error(msg);
    }
}

const fetchAvailablePlotTypes = async (user_id, session_id, dataset_id, analysis_id) => {
    const payload = {user_id, session_id, dataset_id, analysis_id};
    try {
        const {data} = await axios.post(`/api/h5ad/${dataset_id}/availableDisplayTypes`, payload);
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch compatible plot types for this dataset. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchDatasetDisplay = async (display_id) => {
    const payload = {display_id};
    try {
        // POST due to payload variables being sensitive
        const {data} = await axios.post("/cgi/get_dataset_display.cgi", convertToFormData(payload));
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch this saved display for this dataset. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchDatasetDisplayImage = async (dataset_id, display_id) => {
    const payload = {dataset_id, display_id};
    try {
        // POST due to payload variables being sensitive
        const {data} = await axios.post("/cgi/get_dataset_display_image.cgi", convertToFormData(payload));
        return data
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch the image preview for this dataset display. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchDatasetDisplays = async (user_id, dataset_id) => {
    const payload = {user_id, dataset_id};
    try {
        // POST due to payload variables being sensitive
        const {data} = await axios.post("/cgi/get_dataset_displays.cgi", convertToFormData(payload));
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch the saved displays for this dataset. Please contact the gEAR team."
        createToast(msg);
        return [];  // Send an empty list of displays
    }
}

const fetchDatasetInfo = async (dataset_id) => {
    const payload = {dataset_id};
    try {
        const {data} = await axios.post("/cgi/get_dataset_info.cgi", convertToFormData(payload));
        const {title, is_public, owner_id} = data;
        return {title, is_public, owner_id};
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch metadata for this dataset. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchDatasets = async (session_id) => {
    const payload = {session_id}
    try {
        const {data} = await axios.post("cgi/get_h5ad_dataset_list.cgi", convertToFormData(payload));
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch datasets. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchDefaultDisplay = async (user_id, dataset_id) => {
    const payload = {user_id, dataset_id};
    try {
        // POST due to payload variables being sensitive
        const {data} =  await axios.post("/cgi/get_default_display.cgi", convertToFormData(payload));
        const {default_display_id} = data;
        return default_display_id;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch default display for this dataset. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchGeneSymbols = async (datasetId, analysisId) => {
    let url = `./api/h5ad/${datasetId}/genes`;
    if (analysisId) url += `?analysis_id=${analysisId}`;

    try {
        const { data } = await axios.get(url);
        return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene symbols for this dataset. Please contact the gEAR team."
        createToast(msg);
        return [];
    }
}

const fetchH5adInfo = async (datasetId, analysisId) => {
    let url = `/api/h5ad/${datasetId}`
    if (analysisId) url += `?analysis_id=${analysisId}`;
    try {
        const {data} = await axios.get(url);
        const { obs_columns, obs_levels } = data;
        return { obs_columns, obs_levels };
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch H5AD observation data for this dataset. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchPlotlyData = async (plotConfig, datasetId, plot_type, analysis, analysis_owner_id, colorblind_mode)  => {
    // NOTE: gene_symbol already passed to plotConfig
    const payload = { ...plotConfig, plot_type, analysis, analysis_owner_id, colorblind_mode };
    try {
        const { data } = await axios.post(`/api/plot/${datasetId}`, payload);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        return data
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not create Plotly plot for this dataset and parameters. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchSvgData = async (datasetId, geneSymbol) => {
    try {
        const { data } = await axios.get(`/api/plot/${datasetId}/svg?gene=${geneSymbol}`);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        return data
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch SVG data for this dataset and parameters. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
};

const fetchTsneImageData = async (plotConfig, datasetId, plot_type, analysis, analysis_owner_id, colorblind_mode) => {
    // NOTE: gene_symbol already passed to plotConfig
    const payload = { ...plotConfig, plot_type, analysis, analysis_owner_id, colorblind_mode };
    try {
        const { data } = await axios.post(`/api/plot/${datasetId}/tsne`, payload);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not create plot image for this dataset and parameters. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

/* Get HTML element value to save into plot config */
const getPlotConfigValueFromClassName = (className) => {
    // NOTE: Some elements are only present in certain plot type configurations

    const classElts = document.getElementsByClassName(className);

    if (classElts) {
        const elt = classElts[0];   // All elements for this class should have their checks and values synced up
        if (elt?.type == "checkbox") {
            return elt.checked;
        }
        return elt.value || undefined;
    }
    return undefined;
}

/* Get order of series from sortable lists. Return object */
const getPlotOrderFromSortable = () => {
    const order = {};
    for (const elt of document.getElementById("order_container").children) {
        const series = elt.querySelector("p").textContent;
        const serialized = sortable(`#${series}_order_list`, 'serialize')[0].items;
        // Sort by "sortable" index position
        order[series] = serialized.map((val) => val.label);
    }
    return order;
}

const getPlotlyDisplayUpdates = (plotConfObj, plotType, category) => {
    // Get updates and additions to plot from the plot_display_config JS object
    let updates = {};
    for (const idx in plotConfObj) {
        const conf = plotConfObj[idx];
        // Get config (data and/or layout info) for the plot type chosen, if it exists
        if (conf.plot_type == "all" || conf.plot_type == plotType) {
            const update = category in conf ? conf[category] : {};
            updates = {...updates, ...update};    // Merge updates
        }
    }
    return updates;
}

const getSelect2Value = (select) => {
    // Get value from select2 element
    return select.selectedOptions[0].data.value;
}

const includePlotParamOptions = async () => {
    const plotType = getSelect2Value(plotTypeSelect);

    // include plotting backend options
    if (plotlyPlots.includes(plotType)) {
        plotStyle = new PlotlyHandler(plotType);
    } else if (scanpyPlots.includes(plotType)) {
        plotStyle = new ScanpyHandler(plotType);
    } else if (plotType === "svg") {
        plotStyle = new SvgHandler();
    } else {
        console.warn(`Plot type ${plotType} not recognized.`)
        return;
    }

    await plotStyle.loadPlotHtml();

    // NOTE: Events are triggered in the order they are regstered.
    // We want to trigger a param sync event before the plot requirement validation
    // (since it checks every element in the js-plot-req class)
    setupParamValueCopyEvents(Object.keys(plotStyle.classElt2Prop));    // Ensure pre- and post- plot view params are synced up
    setupValidationEvents();        // Set up validation events required to plot at a minimum
    plotStyle.setupPlotSpecificEvents()       // Set up plot-specific events

}

const logErrorInConsole = (error) => {
    if (error.response) {
        // The request was made and the server responded with a status code
        // that falls out of the range of 2xx
        console.error(error.response.data);
        console.error(error.response.status);
        console.error(error.response.headers);
    } else if (error.request) {
        // The request was made but no response was received
        // `error.request` is an instance of XMLHttpRequest in the browser and an instance of
        // http.ClientRequest in node.js
        console.error(error.request);
    } else {
        // Something happened in setting up the request that triggered an Error
        console.error('Error', error.message);
    }
    console.error(error.config);
}


/* Transform and load dataset data into a "tree" format */
const loadDatasetTree = async () => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const datasetData = await fetchDatasets(sessionId);

        let counter = 0;

        // Populate select box with dataset information owned by the user
        if (datasetData.user.datasets.length > 0) {
            // User has some profiles
            for (const item of datasetData.user.datasets) {
                if (item) {
                    userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        // Next, add datasets shared with the user
        if (datasetData.shared_with_user.datasets.length > 0) {
            for (const item of datasetData.shared_with_user.datasets) {
                if (item) {
                    sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        // Now, add public datasets
        if (datasetData.public.datasets.length > 0) {
            for (const item of datasetData.public.datasets) {
                if (item) {
                    domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
    } catch (error) {
        document.getElementById("dataset_s_failed").style.display = "";
    }
    datasetTree.userDatasets = userDatasets;
    datasetTree.sharedDatasets = sharedDatasets;
    datasetTree.domainDatasets = domainDatasets;
    datasetTree.generateTree();
}

const renderColorPicker = (seriesName) => {
    const colorsContainer = document.getElementById("colors_container");
    const colorsSection = document.getElementById("colors_section");

    colorsSection.style.display = "none";
    colorsContainer.replaceChildren();
    if (!seriesName) {
        return;
    }

    if (!(catColumns.includes(seriesName))) {
        // ? Continuous series colorbar picker
        return;
    }

    const seriesNameHtml = generateElements(`<p class="has-text-weight-bold is-underlined">${seriesName}</p>`);
    colorsContainer.append(seriesNameHtml);

    // Otherwise d3 category10 colors
    const swatchColors = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"];

    let counter = 0;
    for (const group of levels[seriesName]) {
        const darkerLevel = Math.floor(counter / 10);
        const baseColor = swatchColors[counter%10];
        const groupColor = darkerLevel > 0
            ? d3.color(baseColor).darker(darkerLevel).formatHex()
            : baseColor;    // Cycle through swatch but make darker if exceeding 10 groups
        counter++;
        const groupHtml = generateElements(`
            <p class="is-flex is-justify-content-space-between">
                <span class="has-text-weight-medium">${group}</span>
                <input class="js-plot-color" id="${group}_color" type="color" value="${groupColor}" aria-label="Select a color" />
            </p>
        `);
        colorsContainer.append(groupHtml)
    }

    colorsSection.style.display = "";
}

const renderDisplayCards = async (userDisplays, ownerDisplays, defaultDisplayId) => {
    const userDisplaysElt = document.getElementById("user_displays");
    const ownerDisplaysElt = document.getElementById("owner_displays");

    // Empty existing displays
    userDisplaysElt.replaceChildren();
    ownerDisplaysElt.replaceChildren();

    for (const display of userDisplays) {
        renderUserDisplayCard(display, defaultDisplayId)
    }

    for (const display of ownerDisplays) {
        renderOwnerDisplayCard(display, defaultDisplayId);
    }
}

// Render the specific series as a sortable list, if it is not already
const renderOrderSortableSeries = (series) => {
    const orderContainer = document.getElementById("order_container");

    // If continouous series, cannot sort.
    if (!catColumns.includes(series)) return;

    // If series is used in another param, also return
    // NOTE: if the original series is removed, the param will not simply be switch over to the new one.
    if (document.getElementById(`${series}_order`)) return;

    // Create parent template
    // Designed so the title is a full row and the draggables are 50% width
    const parentList = `<ul id="${series}_order_list" class="content column is-two-thirds js-plot-order-sortable"></ul>`;
    const template = `
        <div id="${series}_order" class="columns is-multiline">
        <p id="${series}_order_title" class="has-text-weight-bold column is-full">${series}</p>
        ${parentList}
        </div
    `;

    const htmlCollection = generateElements(template);
    orderContainer.append(htmlCollection);

    // Add in list elements
    for (const group of levels[series]) {
        const listElt = `<li class="has-background-grey-lighter has-text-dark">${group}</li>`;
        const listCollection = generateElements(listElt);
        document.getElementById(`${series}_order_list`).append(listCollection);
    }

    // Create sortable for this series
    sortable(`#${series}_order_list`, {
        hoverClass: "has-text-weight-bold"
        , itemSerializer(item, container) {
            item.label = item.node.textContent
            return item
        },
    });

}

const renderOwnerDisplayCard = async (display, defaultDisplayId) => {
    let displayUrl = "";
    try {
        displayUrl = await fetchDatasetDisplayImage(datasetId, display.id);
    } catch (error) {
        displayUrl = "/img/dataset_previews/missing.png";
    }
    const geneSymbol = display.plotly_config.gene_symbol;

    const label = display.label || `Unnamed ${display.plot_type} display`;

    const template = `
                <div id="${display.id}_display" class="column is-one-quarter">
                    <div class="box card has-background-primary-light has-text-primary">
                        <header class="card-header">
                            <p class="card-header-title has-text-black">${label}</p>
                        </header>
                        <div class="card-image">
                            <figure class="image is-4by3">
                                <img src="${displayUrl} " alt="Saved display">
                            </figure>
                        </div>
                        <div class="card-content">
                            <p class="subtitle is-6">Gene: ${geneSymbol}</p>
                        </div>
                        <footer class="card-footer buttons">
                            <button class="js-display-default card-footer-item button is-primary" id="${display.id}_default">Set as Default</button>
                            <button class="card-footer-item button is-primary" id="${display.id}_clone">Clone</button>
                        </footer>
                    </div>
                </div>`;

    const htmlCollection = generateElements(template);
    const ownerDisplaysElt = document.getElementById("owner_displays");
    ownerDisplaysElt.append(htmlCollection);

    // Edit default properties if this is the default display
    const defaultElt = document.getElementById(`${display.id}_default`);
    if (display.id === defaultDisplayId) {
        defaultElt.textContent = "Default";
        defaultElt.disabled = true;
    }

    defaultElt.addEventListener("click", (event) => saveDefaultDisplay(display.id));
    document.getElementById(`${display.id}_clone`).addEventListener("click", (event) => cloneDisplay(event, display));
}

const renderUserDisplayCard = async (display, defaultDisplayId) => {

    let displayUrl = "";
    try {
        displayUrl = await fetchDatasetDisplayImage(datasetId, display.id);
    } catch (error) {
        displayUrl = "/img/dataset_previews/missing.png";
    }
    const geneSymbol = display.plotly_config.gene_symbol;

    const label = display.label || "";

    // TODO - Get footer button styles correct
    const template = `
                <div id="${display.id}_display" class="column is-one-quarter">
                    <div class="box card has-background-primary-light has-text-primary">
                        <header class="card-header">
                            <p class="card-header-title has-text-black">${label}</p>
                        </header>
                        <div class="card-image">
                            <figure class="image is-4by3">
                                <img src="${displayUrl} " alt="Saved display">
                            </figure>
                        </div>
                        <div class="card-content">
                            <p class="subtitle">Gene: ${geneSymbol}</p>
                        </div>
                        <footer class="card-footer buttons">
                            <button class="js-display-default card-footer-item button is-primary" id="${display.id}_default">Set as Default</button>
                            <button class="card-footer-item button is-primary" id="${display.id}_clone">Clone</button>
                            <button class="card-footer-item button is-danger" id="${display.id}_delete">Delete</button>
                        </footer>
                    </div>
                </div>`;

    const htmlCollection = generateElements(template);
    const userDisplaysElt = document.getElementById("user_displays");
    userDisplaysElt.append(htmlCollection);

    // Edit default properties if this is the default display
    const defaultElt = document.getElementById(`${display.id}_default`);
    if (display.id === defaultDisplayId) {
        defaultElt.textContent = "Default";
        defaultElt.disabled = true;
    }

    defaultElt.addEventListener("click", (event) => saveDefaultDisplay(display.id));
    document.getElementById(`${display.id}_clone`).addEventListener("click", (event) => cloneDisplay(event, display));
    document.getElementById(`${display.id}_delete`).addEventListener("click", (event) => deleteDisplay(userId, display.id));
}

const saveDatasetDisplay = async(displayId, dataset_id, user_id, label, plot_type, plotConfig) => {
    // NOTE: Saving all displays as new displays (clone) instead of overwriting. User can always delete excess displays
    const payload = {
        id: displayId,
        dataset_id,
        user_id,
        label,
        plot_type,
        plotly_config: JSON.stringify({
            ...plotConfig,  // depending on display type, this object will have different properties
        }),
    };
    if (!displayId) delete payload.id;  // Prevent passing in "null" as a string.

    try {
        const data = await axios.post("/cgi/save_dataset_display.cgi", convertToFormData(payload));
        const {display_id, success} = data;
        if (!success) {
            throw new Error("Could not save this new display. Please contact the gEAR team.");
        }

        // Make new display card and make it the default display
        renderUserDisplayCard({id: display_id, label, plot_type, plotly_config: payload.plotly_config}, display_id);

        return display_id;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not save this new display. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const saveDefaultDisplay = async (displayId) => {
    const payload = {display_id: displayId, user_id: userId};
    try {
        const {success} = await axios.post("/cgi/save_default_display.cgi", convertToFormData(payload))
        if (!success) {
            throw new Error("Could not save this display as your default. Please contact the gEAR team.");
        }
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not save this display as your default. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    };

    //Update labels of displays... this becomes "Default", others become "Make Default"
    for (const elt of document.getElementsByClassName("js-display-default")) {
        elt.disabled = false;
        elt.textContent = "Set as Default";
    }

    const currentDefaultElt = document.getElementById(`${displayId}_default`);
    currentDefaultElt.disabled = true;
    currentDefaultElt.textContent = "Default";
}

/* Set HTML element value from the plot config value */
const setPlotEltValueFromConfig = (classSelector, confVal) => {
    for (const elt of document.getElementsByClassName(classSelector)) {
        if (elt.type == "checkbox") {
            elt.checked = confVal;
            return;
        }
        elt.value = confVal;
    }
}

/* Set disabled state for the given plot type. Also normalize plot type labels */
const setPlotTypeDisabledState = (plotType, isAllowed) => {
    if (plotType === "tsne/umap_dynamic") {
        document.getElementById("tsne_dyna_opt").disabled = !isAllowed;
    } else {
        document.getElementById(`${plotType}_opt`).disabled = !isAllowed;
    }
}

/**
 * Set Select Box Selection By Value
 * Modified to set value and "selected" so nice-select2 update() will catch it
 * Taken from https://stackoverflow.com/a/20662180
 * @param eid Element ID
 * @param eval Element value
 */
const setSelectBoxByValue = (eid, val) => {
    const elt = document.getElementById(eid);
    for (const i in elt.options) {
        if (elt.options[i].value === val) {
            elt.value = val;
            elt.options[i].selected = true;
            return;
        }
    }
}

/* Ensure all elements in this class have the same value */
const setupParamValueCopyEvents = (plotSpecificParams) => {
    for (const classSelector of plotSpecificParams) {
        const classElts = document.getElementsByClassName(classSelector)
        for (const elt of classElts) {
            elt.addEventListener("change", (event) => {
                for (const classElt of classElts) {
                    classElt.value = event.target.value;
                    // Believe that programmatically changing the value does not trigger "change" (AKA no cascading)
                    classElt.disabled = event.target.disabled;
                    if (event.target.type == "checkbox") classElt.checked = event.target.checked;
                }
            });
        }
    }
}

/* Set up any Plotly-based plot options and events for the pre- and post- plot views */
const setupPlotlyOptions = async () => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const plotType = getSelect2Value(plotTypeSelect);
    try {
        ({obs_columns: allColumns, obs_levels: levels} = await fetchH5adInfo(datasetId, analysisId));
    } catch (error) {
        document.getElementById("plot_options_s_failed").style.display = "";
        return;
    }
    catColumns = Object.keys(levels);

    // class name, list of columns, add expression, default category
    updateSeriesOptions("js-plotly-x-axis", allColumns, true);
    updateSeriesOptions("js-plotly-y-axis", allColumns, true, "raw_value");
    updateSeriesOptions("js-plotly-color", allColumns, true);
    updateSeriesOptions("js-plotly-label", allColumns, true);
    updateSeriesOptions("js-plotly-facet-row", catColumns, false);
    updateSeriesOptions("js-plotly-facet-col", catColumns, false);

    if (["scatter", "tsne_dyna"].includes(plotType)) {
        const difference = (arr1, arr2) => arr1.filter(x => !arr2.includes(x))
        const continuousColumns = difference(allColumns, catColumns);

        updateSeriesOptions("js-plotly-size", continuousColumns, true);
    }

    if (["scatter", "tsne_dyna"].includes(plotType)) {
        const difference = (arr1, arr2) => arr1.filter(x => !arr2.includes(x))
        const continuousColumns = difference(allColumns, catColumns);

        updateSeriesOptions("js-plotly-size", continuousColumns, true);

        // If x-axis is continuous show vline stuff, otherwise hide
        // If x-axis is categorical, enable jitter plots
        const xAxisSeriesElts = document.getElementsByClassName("js-plotly-x-axis");
        for (const elt of xAxisSeriesElts) {
            elt.addEventListener("change", (event) => {

                const jitterElts = document.getElementsByClassName("js-plotly-add-jitter");
                const vLinesContainer = document.getElementById("vlines_container")

                if ((catColumns.includes(event.target.value))) {
                    // categorical x-axis
                    for (const jitterElt of jitterElts) {
                        jitterElt.style.display = "";
                    }
                    vLinesContainer.style.display = "none";
                    // Remove all but first existing vline
                    const toRemove = document.querySelectorAll(".js-plotly-vline-field:not(:first-of-type)");
                    for (const elt of toRemove) {
                        elt.remove();
                    }
                    // Clear the first vline
                    document.querySelector(".js-plotly-vline-pos").value = "";
                    document.querySelector(".js-plotly-vline-style-select").value = "solid";
                } else {
                    for (const jitterElt of jitterElts) {
                        jitterElt.style.display = "none";
                        jitterElt.checked = false;
                    }

                    vLinesContainer.style.display = "";

                }

            });
        }


        // Vertical line add and remove events
        const vLinesBody = document.getElementById("vlines_body");
        const vLineField = document.querySelector(".js-plotly-vline-field");
        document.getElementById("vline_add_btn").addEventListener("click", (event) => {
            vLinesBody.append(vLineField.cloneNode(true));
            // clear the values of the clone
            document.querySelector(".js-plotly-vline-field:last-of-type .js-plotly-vline-pos").value = "";
            document.querySelector(".js-plotly-vline-field:last-of-type .js-plotly-vline-style-select").value = "solid";
            // NOTE: Currently if original is set before cloning, values are copied to clone
            document.getElementById("vline_remove_btn").disabled = false;
        })
        document.getElementById("vline_remove_btn").addEventListener("click", (event) => {
            // Remove last vline
            const lastVLine = document.querySelector(".js-plotly-vline-field:last-of-type");
            lastVLine.remove();
            if (vLineField.length < 2) document.getElementById("vline_remove_btn").disabled = true;
        })
    }

    // If color series is selected, let user choose colors.
    const colorSeriesElts = document.getElementsByClassName("js-plotly-color");
    for (const elt of colorSeriesElts) {
        elt.addEventListener("change", (event) => {
            if ((catColumns.includes(event.target.value))) {
                renderColorPicker(event.target.value);
            }
        })
    }


    // Certain elements trigger plot order
    const plotOrderElts = document.getElementsByClassName("js-plot-order");
    for (const elt of plotOrderElts) {
        elt.addEventListener("change", (event) => {
            const paramId = event.target.id;
            const param = paramId.replace("_series", "").replace("_post", "");
            // NOTE: continuous series will be handled in the function
            updateOrderSortable();
        });
    }

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    if (document.getElementById("x_axis_series").value) {
        trigger(document.getElementById("x_axis_series"), "change");
    }
    if (document.getElementById("y_axis_series").value) {
        trigger(document.getElementById("y_axis_series"), "change");
    }

    // Setup the dropdown menu on the post-plot view
    const plotlyDropdown = document.getElementById("plotly_param_dropdown");
    plotlyDropdown.addEventListener("click", (event) => {
        event.stopPropagation();    // This prevents the document from being clicked as well.
        plotlyDropdown.classList.toggle("is-active");
    })

    // Close dropdown if it is clicked off of, or ESC is pressed
    // https://siongui.github.io/2018/01/19/bulma-dropdown-with-javascript/#footnote-1
    document.addEventListener('click', () => {
        plotlyDropdown.classList.remove("is-active");
    });
    document.addEventListener('keydown', (event) => {
        if (event.key === "Escape") {
            plotlyDropdown.classList.remove("is-active");
        }
    });

    const plotlyDropdownMenuItems = document.querySelectorAll("#plotly_param_dropdown .dropdown-item");
    for (const item of plotlyDropdownMenuItems) {
        item.addEventListener("click", showPostPlotlyParamSubsection);
    }

}

/* Set up any Scanpy-based plot options and events for the pre- and post- plot views */
const setupScanpyOptions = async () => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const plotType = getSelect2Value(plotTypeSelect);
    try {
        ({obs_columns: allColumns, obs_levels: levels} = await fetchH5adInfo(datasetId, analysisId));
    } catch (error) {
        document.getElementById("plot_options_s_failed").style.display = "";
        return;
    }
    catColumns = Object.keys(levels);

    let xDefaultOption = null;
    let yDefaultOption = null;

    // If these exist, make the default option
    switch (plotType) {
        case "pca_static":
            xDefaultOption = "X_pca_1";
            yDefaultOption = "X_pca_2";
            break;
        case "tsne_static":
            xDefaultOption = "X_tsne_1";
            yDefaultOption = "X_tsne_2";
            break;
        case "umap_static":
            xDefaultOption = "X_umap_1";
            yDefaultOption = "X_umap_2";
            break;
    }

    updateSeriesOptions("js-tsne-x-axis", allColumns, true, xDefaultOption);
    updateSeriesOptions("js-tsne-y-axis", allColumns, true, yDefaultOption);
    updateSeriesOptions("js-tsne-colorize-legend-by", allColumns, false);
    updateSeriesOptions("js-tsne-plot-by-series", catColumns, false);

    const colorizeLegendBy = document.getElementsByClassName("js-tsne-colorize-legend-by");
    const plotBySeries = document.getElementsByClassName("js-tsne-plot-by-series");
    const maxColumns = document.getElementsByClassName('js-tsne-max-columns');
    const skipGenePlot = document.getElementsByClassName("js-tsne-skip-gene-plot");
    const horizontalLegend = document.getElementsByClassName("js-tsne-horizontal-legend");

    // Do certain things if the chosen annotation series is categorical or continuous
    for (const elt of colorizeLegendBy) {
        elt.addEventListener("change", (event) => {
            for (const targetElt of [...plotBySeries, ...maxColumns, ...horizontalLegend]) {
                targetElt.disabled = true;
                // If colorized legend is continuous, we cannot plot by group
                // So all dependencies need to be disabled.
                if ((catColumns.includes(event.target.value))) {
                    renderColorPicker(event.target.value);
                    targetElt.disabled = false;
                }
            }
        });

        elt.addEventListener("change", (event) => {
            renderColorPicker(event.target.value);
            return;
        })
    }

    // Plotting by group plots gene expression, so cannot skip gene plots.
    for (const elt of plotBySeries) {
        elt.addEventListener("change", (event) => {
            // Must plot gene expression if series value selected
            for (const targetElt of skipGenePlot) {
                targetElt.disabled = event.target.value ? true : false;
                if (event.target.value) targetElt.checked = false;
            }
            // Must be allowed to specify max columns if series value selected
            for (const targetElt of maxColumns) {
                targetElt.disabled = event.target.value ? false : true;
            }
            updateOrderSortable();

        });
    }

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    if (document.getElementById("x_axis_series").value) {
        trigger(document.getElementById("x_axis_series"), "change");
    }
    if (document.getElementById("y_axis_series").value) {
        trigger(document.getElementById("y_axis_series"), "change");
    }

}

/* Set up any SVG options and events for the pre- and post- plot views */
const setupSVGOptions = () => {
    const enableMidColorElts = document.getElementsByClassName("js-svg-enable-mid");
    const midColorElts = document.getElementsByClassName("js-svg-mid-color");
    const midColorFields = document.getElementsByClassName("js-mid-color-field");
    for (const elt of enableMidColorElts) {
        elt.addEventListener("change", (event) => {
            for (const field of midColorFields) {
                field.style.display = (event.target.checked) ? "" : "none";
            }
            for (const midColor of midColorElts) {
                midColor.disabled = !(event.target.checked);
            }
        });
    }


    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    if (document.getElementById("low_color").value) {
        trigger(document.getElementById("low_color"), "change");
    }
    if (document.getElementById("high_color").value) {
        trigger(document.getElementById("high_color"), "change");
    }
}

/* Setup a fail-fast validation trigger. */
const setupValidationEvents = () => {
    const validationElts = document.getElementsByClassName("js-plot-req");
    for (const elt of validationElts ) {
        elt.addEventListener("change", () => {
            // Reset "status" classes
            elt.classList.remove("is-success", "is-danger");
            if (elt.value) {
                elt.parentElement.classList.remove("is-danger");
                elt.parentElement.classList.add("is-success");

                // If every validation param has been filled out, it's OK to plot
                // NOTE: need to ensure pre- and post- param elements are filled before this function is called
                if ([...validationElts].every(element => element.value)) {
                    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
                        plotBtn.disabled = false;
                    }
                    document.getElementById("plot_options_s_success").style.display = "";
                    document.getElementById("plot_options_s_failed").style.display = "none";
                }

                return;
            }

            // Required paramater has no value. Indicate it and disable plot buttons
            elt.parentElement.classList.add("is-danger");
            elt.parentElement.classList.remove("is-success");

            for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
                plotBtn.disabled = true;
            }
            document.getElementById("plot_options_s_success").style.display = "none";
            document.getElementById("plot_options_s_failed").style.display = "";
        })
    }

}

const showPostPlotlyParamSubsection = (event) => {
    for (const subsection of document.getElementsByClassName("js-plot-config-section")) {
        subsection.style.display = "none";
    }

    switch (event.target.textContent.trim()) {
        case "X-axis":
            document.getElementById("x_axis_section_post").style.display = "";
            break;
        case "Y-axis":
            document.getElementById("y_axis_section_post").style.display = "";
            break;
        case "Color":
            document.getElementById("color_section_post").style.display = "";
            break;
        case "Marker Size":
            document.getElementById("size_section_post").style.display = "";
            break;
        case "Subplots":
            document.getElementById("subplots_section_post").style.display = "";
            break;
        default:
            document.getElementById("misc_section_post").style.display = "";
            break;
    }
    event.preventDefault(); // Prevent "link" clicking from "a" elements
}

const updateAnalysesOptions = (privateAnalyses, publicAnalyses) => {
    const privateAnalysesElt = document.getElementById("private_analyses");
    const publicAnalysesElt = document.getElementById("public_analyses");

    // Empty the old optgroups
    privateAnalysesElt.replaceChildren();
    publicAnalysesElt.replaceChildren();

    // Show message that no analyses are present if none exist
    if (!privateAnalyses?.length) {
        const option = document.createElement("option");
        option.text = "You have no saved analyses for this dataset.";
        option.disabled = true;
        privateAnalysesElt.append(option);
    }

    if (!publicAnalyses?.length) {
        const option = document.createElement("option");
        option.text = "There are no public analyses for this dataset.";
        option.disabled = true;
        publicAnalysesElt.append(option);
    }

    // Load each analysis as an option
    // NOTE: For now, we can just pass analysis id only to tSNE and be fine
    for (const analysis of privateAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.value = analysis.id;
        //option.dataset.type = analysis.type;
        //option.dataset.owner_id = analysis.user_id;
        privateAnalysesElt.append(option);
    }
    for (const analysis of publicAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.value = analysis.id;
        //option.dataset.type = analysis.type;
        //option.dataset.owner_id = analysis.user_id;
        publicAnalysesElt.append(option);
    }

    // Update seel
    analysisSelect.update();
}

const updateGeneSymbolOptions = (geneSymbols) => {
    const geneSelectElt = document.getElementById("gene_select");

    geneSelectElt.replaceChildren();

    // Append empty placeholder element
    const firstOption = document.createElement("option");
    firstOption.textContent = "Please select a gene";
    geneSelectElt.append(firstOption);

    for (const gene of geneSymbols.sort()) {
        const option = document.createElement("option");
        option.textContent = gene;
        option.value = gene;
        geneSelectElt.append(option);
    }

    // Update the nice-select2 element to reflect this.
    geneSelect.update();
}

const updateOrderSortable = () => {
    // TODO: Need to come up with a way to a) make a set of series, and remove a series if no params use it

    // Get all current plot param series for plotting order and save as a set
    const plotOrderElts = document.getElementsByClassName("js-plot-order");
    const seriesSet = new Set();
    for (const elt of plotOrderElts) {
        const series = elt.value;
        // Only include categorical series
        if (series && (catColumns.includes(series))) {
            seriesSet.add(series);
        }
    }

    // Get all current plotting order series and save as a set
    const sortableElts = document.querySelectorAll(".js-plot-order-sortable p");
    const sortableSet = new Set();
    for (const elt of sortableElts) {
        const series = elt.value;
        // These series already are categorical
        if (series) {
            sortableSet.add(series);
        }
    }

    // Three scenarios:
    for (const series of seriesSet) {
        // 1. Series is in both seriesSet and sortableSet, do nothing
        if (sortableSet.has(series)) {
            continue;
        }
        // 2. Series is in seriesSet but not sortableSet, add <series>_order element
        renderOrderSortableSeries(series);
    }

    for (const series of sortableSet) {
        // 3. Series is in sortableSet but not seriesSet, remove <series>_order element
        if (!seriesSet.has(series)) {
            const orderElt = document.getElementById(`${series}_order`);
            orderElt.remove();
        }
    }


    const orderContainer = document.getElementById("order_container");
    const orderSection = document.getElementById("order_section");

    // Pre-emptively hide the container but show it ass
    if (!orderContainer.children.length) {
        orderSection.style.display = "none";
        return;
    }

    orderSection.style.display = "";

}

// For plotting options, populate select menus with category groups
const updateSeriesOptions = (classSelector, seriesArray, addExpression, defaultOption) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add an expression option (since expression is not in the categories)
        if (addExpression) {
            const expression = document.createElement("option");
            elt.append(expression);
            expression.textContent = "expression";
            expression.value = "raw_value";
            if ("raw_value" === defaultOption) {
                expression.selected = true;
            }
        }

        // Add categories
        for (const group of seriesArray.sort()) {
            // Skip columns listed as "_colors" as they just provide colors for another series
            if (group.includes("_colors")) continue;

            const option = document.createElement("option");
            option.textContent = group;
            // Change X_pca/X_tsne/X_umap text content to be more user_friendly
            if (group.includes("X_") && (
                group.includes("pca")
                || group.includes("tsne")
                || group.includes("umap")
            )) {
                option.textContent = `${group} (from selected analysis)`;
            }
            option.value = group;
            elt.append(option);
            // NOTE: It is possible for a default option to not be in the list of groups.
            if (group === defaultOption) {
                option.selected = true;
            }
        }
    }

}

window.onload = () => {

    loadDatasetTree();

    // If brought here by the "gene search results" page, curate on the dataset ID that referred us
    const urlParams = new URLSearchParams(window.location.search);
    if (urlParams.has("dataset_id")) {
        const linkedDatasetId = URLSearchParams.get("dataset_id");
        try {
            // TODO: test this
            // find DatasetTree node and trigger "activate"
            const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
            foundNode.setActive(true);
            datasetId = linkedDatasetId;
        } catch {
            console.error(`Dataset id ${linkedDatasetId} was not returned as a public/private/shared dataset`);
        }
    }

    createAnalysisSelectInstance();
    createPlotSelectInstance();
    createGeneSelectInstance();
};

document.getElementById("new_display").addEventListener("click", chooseDisplay);
document.getElementById("analysis_select").addEventListener("change", chooseAnalysis);
document.getElementById("plot_type_select").addEventListener("change", choosePlotType);
document.getElementById("gene_select").addEventListener("change", chooseGene);
const plotBtns = document.getElementsByClassName("js-plot-btn");
for (const plotBtn of plotBtns) {
    plotBtn.addEventListener("click", createPlot);
}

document.getElementById("save_json_config").addEventListener("click", () => {
    // Config plot configuration to JSON for sharing (or passing to API by hand)
    const blob = new Blob([JSON.stringify({...this.plotConfig, plot_type:this.plotType})]);
    const link = document.createElement("a");
    link.download = "gEAR_plot_configuration.json";
    link.href = window.URL.createObjectURL(blob);
    link.click()
    // Give confirmation
    document.getElementById("save_json_config").classList.add("is-success");
    setTimeout(() => {
        document.getElementById("save_json_config").classList.remove("is-success");
    }, 1000);
});
document.getElementById("save_display_btn").addEventListener("click", async () => {
    // Save new plot display.
    const label = document.getElementById("new_display_label").value;
    const displayId = await saveDatasetDisplay(null, datasetId, userId, label, this.plotType, this.plotConfig);
    saveDefaultDisplay(displayId);
    // Give confirmation
    document.getElementById("save_display_btn").classList.add("is-success");
    setTimeout(() => {
        document.getElementById("save_display_btn").classList.remove("is-success");
    }, 1000);
    // TODO: Reload new displays if user goes back to parameter editing
});

document.getElementById("edit_params").addEventListener("click", () => {
    // Hide this view
    document.getElementById("content_c").style.display = "";
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").style.display = "none";
})
