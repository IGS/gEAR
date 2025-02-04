// I use camelCase for my variable/function names to adhere to JS style standards
// Exception being functions that do fetch calls, so we can use JS destructuring on the payload

'use strict';

const isMultigene = 0;

let selectedGene = null;

const plotlyPlots = ["bar", "line", "scatter", "tsne_dyna", "violin"];
const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];

/**
 * Represents a PlotlyHandler, a class that handles Plotly plots.
 * @class
 * @extends PlotHandler
 */
class PlotlyHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
        this.apiPlotType = plotType;
    }

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
        , "js-plotly-color-palette":"color_palette"
        , "js-plotly-reverse-palette":"reverse_palette"
    }

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    /**
     * Clones the display based on the provided configuration.
     *
     * @param {Object} config - The configuration object.
     */
    cloneDisplay(config) {
        // plotly plots
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            for (const series in config["order"]) {
                const order = config["order"][series];
                // sort "levels" series by order
                levels[series].sort((a, b) => order.indexOf(a) - order.indexOf(b));
                renderOrderSortableSeries(series);
            }

            document.getElementById("order-section").classList.remove("is-hidden");
        }

        // Handle filters
        if (config["obs_filters"]) {
            facetWidget.filters = config["obs_filters"];
        }

        // Handle colors
        if (config["colors"]) {
            // do nothing if color_name is not set
            if (!config["color_name"]) return;

            try {
                const series = config["color_name"];
                renderColorPicker(series);
                for (const group in config["colors"]) {
                    const color = config["colors"][group];
                    // Sometimes we need to escape the group name
                    // Found a case where the group in config was truncated compared to the (older) dataset's actual group
                    // But there are cases where the group name should not be escaped (such as if slashes are in the name)
                    try {
                        const colorField = document.getElementById(`${group}-color`);
                        colorField.value = color;
                    } catch (error) {
                        const colorField = document.getElementById(`${CSS.escape(group)}-color`);
                        colorField.value = color;
                    }
                }
            } catch (error) {
                console.error(error);
                // pass
            }
        }

        if (config["color_palette"]) {
            setSelectBoxByValue("color-palette-post", config["color_palette"]);
        }

        // Handle vlines
        if (config["vlines"]) {
            const vLinesBody = document.getElementById("vlines-body");
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

    /**
     * Creates a plot using the provided dataset ID and analysis object.
     * @param {string} datasetId - The ID of the dataset.
     * @param {Object} analysisObj - The analysis object.
     * @returns {void}
     */
    async createPlot(datasetId, analysisObj) {
        // Get data and set up the image area
        let plotJson;
        try {
            const data = await fetchPlotlyData(datasetId, analysisObj, this.apiPlotType, this.plotConfig);
            ({plot_json: plotJson} = data);
        } catch (error) {
            return;
        }

        const plotContainer = document.getElementById("plot-container");
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const plotlyPreview = document.createElement("div");
        plotlyPreview.classList.add("container", "is-max-desktop");
        plotlyPreview.id = "plotly-preview";
        plotContainer.append(plotlyPreview);
        Plotly.purge("plotly-preview"); // clear old Plotly plots

        if (!plotJson) {
            createToast("Could not retrieve plot information. Cannot make plot.");
            return;
        }
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "config");
        Plotly.newPlot("plotly-preview", plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "layout")
        Plotly.relayout("plotly-preview", custonLayout)

        addOvercrowdedSeriesWarning(plotContainer);

    }

    /**
     * Loads the plot HTML by replacing the content of prePlotOptionsElt and postPlotOptionsElt elements.
     * Populates advanced options for specific plot types.
     * @returns {Promise<void>} A promise that resolves when the plot HTML is loaded.
     */
    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot-options-collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post-plot-adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/single_gene_plotly.html");
        postPlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/post_plot/single_gene_plotly.html");

        // populate advanced options for specific plot types
        const prePlotSpecificOptionsElt = document.getElementById("plot-specific-options");
        const postPlotSpecificOptionselt = document.getElementById("post-plot-specific-options");

        // Load color palette select options
        if (["violin"].includes(this.plotType)) {
            // TODO: Discrete scale should go to color mapping
            loadColorscaleSelect(false);
        } else {
            loadColorscaleSelect(true);
        }

        if (["scatter", "tsne_dyna"].includes(this.plotType)) {
            prePlotSpecificOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/advanced_scatter.html");
            postPlotSpecificOptionselt.innerHTML = await includeHtml("../include/plot_config/post_plot/advanced_scatter.html");
            return;
        }
        if (this.plotType === "violin") {
            prePlotSpecificOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/advanced_violin.html");
            postPlotSpecificOptionselt.innerHTML = await includeHtml("../include/plot_config/post_plot/advanced_violin.html");
            return;
        }
    }

    /**
     * Populates the plot configuration based on the current state of the dataset curator.
     */
    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = getPlotConfigValueFromClassName(classElt)
        }

        // Small fix for tsne/umap dynamic plots
        if (this.plotType.toLowerCase() === "tsne_dyna") {
            this.apiPlotType = "tsne/umap_dynamic";
        }

        // Violin plots will error (from API) if color_palette is provided
        if (this.plotType.toLowerCase() === "violin") {
            this.plotConfig["color_palette"] = null;
        }

        // Filtered observation groups
        this.plotConfig["obs_filters"] = facetWidget?.filters || {};

        // Get order
        this.plotConfig["order"] = getPlotOrderFromSortable();

        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("color-series-post").value;
        if (colorSeries && colorElts.length) {
            // Input is either color mapping or just the series
            this.plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("-color", "");
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

    /**
     * Sets up the event for copying parameter values.
     * @async
     * @function setupParamValueCopyEvent
     * @returns {Promise<void>}
     */
    async setupParamValueCopyEvent() {
        //pass
    }

    /**
     * Sets up plot-specific events.
     * @async
     * @function setupPlotSpecificEvents
     * @returns {Promise<void>}
     */
    async setupPlotSpecificEvents() {
        await setupPlotlyOptions();
    }

}

/**
 * Represents a ScanpyHandler class that extends PlotHandler.
 * This class is responsible for creating and manipulating plots for a given dataset using the Scanpy analysis object.
 */
class ScanpyHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
        this.apiPlotType = plotType;
    }

    classElt2Prop = {
        "js-tsne-x-axis":"x_axis"
        , "js-tsne-y-axis":"y_axis"
        , "js-tsne-flip-x":"flip_x"
        , "js-tsne-flip-y":"flip_y"
        , "js-tsne-colorize-legend-by":"colorize_legend_by"
        , "js-tsne-plot-by-series":"plot_by_group"
        , "js-tsne-max-columns":"max_columns"
        , "js-tsne-skip-gene-plot":"skip_gene_plot"
        , "js-tsne-horizontal-legend":"horizontal_legend"
        , "js-tsne-marker-size":"marker_size"
        , "js-tsne-color-palette":"expression_palette"
        , "js-tsne-reverse-palette":"reverse_palette"
        , "js-tsne-two-way-palette":"two_way_palette"
        , "js-tsne-center-around-median":"center_around_median"
    }

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    /**
     * Clones the display based on the given configuration.
     * @param {Object} config - The configuration object.
     */
    cloneDisplay(config) {
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            for (const series in config["order"]) {
                const order = config["order"][series];
                // sort "levels" series by order
                levels[series].sort((a, b) => order.indexOf(a) - order.indexOf(b));
                renderOrderSortableSeries(series);
            }

            document.getElementById("order-section").classList.remove("is-hidden");
        }

        // Handle filters
        if (config["obs_filters"]) {
            facetWidget.filters = config["obs_filters"];
        }

        // Restoring some disabled/checked elements in UI
        const plotBySeries = document.getElementsByClassName("js-tsne-plot-by-series");
        const maxColumns = document.getElementsByClassName('js-tsne-max-columns');
        const skipGenePlot = document.getElementsByClassName("js-tsne-skip-gene-plot");
        const horizontalLegend = document.getElementsByClassName("js-tsne-horizontal-legend");

        if (config["colorize_legend_by"]) {
            const series = config["colorize_legend_by"];
            for (const targetElt of [...plotBySeries, ...horizontalLegend]) {
                targetElt.disabled = true;
                if (catColumns.includes(series)) {
                    targetElt.disabled = false;
                }

                // Applies to horizontal legend
                disableCheckboxLabel(targetElt, targetElt.disabled);
            }

            // The "max columns" parameter is only available for categorical series
            if (!(catColumns.includes(series))) {
                for (const targetElt of maxColumns) {
                    targetElt.disabled = true;
                }
            }

            // Handle colors
            if (config["colors"]) {
                renderColorPicker(series);
                for (const group in config["colors"]) {
                    const color = config["colors"][group];
                    // Sometimes we need to escape the group name
                    // Found a case where the group in config was truncated compared to the (older) dataset's actual group
                    // But there are cases where the group name should not be escaped (such as if slashes are in the name)
                    try {
                        const colorField = document.getElementById(`${group}-color`);
                        colorField.value = color;
                    } catch (error) {
                        const colorField = document.getElementById(`${CSS.escape(group)}-color`);
                        colorField.value = color;
                    }
                }
            }
        }

        if (config["expression_palette"]) {
            setSelectBoxByValue("color-palette-post", config["expression_palette"]);
        }

        if (config["plot_by_group"]) {
            for (const targetElt of skipGenePlot) {
                targetElt.disabled = true;
                targetElt.checked = false;
                disableCheckboxLabel(targetElt, targetElt.disabled);
            }
            for (const targetElt of maxColumns) {
                targetElt.disabled = false;
            }
        }

        // If marker size is present, enable the override option
        if (config["marker_size"]) {
            for (const classElt of document.getElementsByClassName("js-tsne-marker-size")) {
                classElt.disabled = false;
            }
            for (const classElt of document.getElementsByClassName("js-tsne-override-marker-size")) {
                classElt.checked = true;
            }
        }

        // If skip gene plot is checked, disable the expression palette, reverse palette, and two-way palette
        if (config["skip_gene_plot"]) {
            for (const classElt of document.getElementsByClassName("js-tsne-expression-palette")) {
                classElt.disabled = true;
            }
            for (const classElt of document.getElementsByClassName("js-tsne-reverse-palette")) {
                classElt.disabled = true;
                classElt.checked = false;
            }
            for (const classElt of document.getElementsByClassName("js-tsne-two-way-palette")) {
                classElt.disabled = true;
                classElt.checked = false;
            }
        }
    }

    /**
     * Creates a plot for a given dataset using the provided analysis object.
     * @param {string} datasetId - The ID of the dataset.
     * @param {Object} analysisObj - The analysis object.
     * @returns {void}
     */
    async createPlot(datasetId, analysisObj) {
        let image;
        try {
            const data = await fetchTsneImage(datasetId, analysisObj, this.apiPlotType, this.plotConfig);
            ({image} = data);
        } catch (error) {
            return;
        }

        const plotContainer = document.getElementById("plot-container");
        plotContainer.replaceChildren();    // erase plot

        const tsnePreview = document.createElement("img");
        tsnePreview.classList.add("image");
        tsnePreview.id = "tsne-preview";
        plotContainer.append(tsnePreview);

        if (!image) {
            createToast("Could not retrieve plot image. Cannot make plot.");
            return;
        }
        const blob = await fetch(`data:image/webp;base64,${image}`).then(r => r.blob());
        tsnePreview.src = URL.createObjectURL(blob);
        tsnePreview.onload = () => {
            // Revoke the object URL to free up memory
            URL.revokeObjectURL(tsnePreview.src);
        }
        return;

    }

    /**
     * Loads the plot HTML by replacing the content of prePlotOptionsElt and postPlotOptionsElt elements.
     * @returns {Promise<void>} A promise that resolves when the plot HTML is loaded.
     */
    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot-options-collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post-plot-adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/tsne_static.html");
        postPlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/post_plot/tsne_static.html");

        loadColorscaleSelect(true, true);
    }

    /**
     * Populates the plot configuration based on various elements and values.
     */
    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = getPlotConfigValueFromClassName(classElt)
        }

        // Get order
        this.plotConfig["order"] = getPlotOrderFromSortable();

        // Filtered observation groups
        this.plotConfig["obs_filters"] = facetWidget?.filters || {};

        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("colorize-legend-by-post").textContent;
        if (colorSeries && colorElts.length) {
            this.plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("-color", "");
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

        // If override marker size is not checked, ensure it does not get passed to the scanpy code
        if (!(document.getElementById("override-marker-size-post").checked)) {
            this.plotConfig["marker_size"] = null;
        }

        // If skip gene plot is checked, ensure expression palette, reverse palette, and two-way palette do not get passed to the scanpy code
        if (document.getElementById("skip-gene-plot-post").checked) {
            this.plotConfig["expression_palette"] = null;
            this.plotConfig["reverse_palette"] = false;
            this.plotConfig["two_way_palette"] = false;
        }
    }

    /**
     * Sets up the event for copying parameter values.
     * @returns {Promise<void>} A promise that resolves when the event setup is complete.
     */
    async setupParamValueCopyEvent() {
        //pass
    }

    /**
     * Sets up plot-specific events.
     * @returns {Promise<void>} A promise that resolves when the setup is complete.
     */
    async setupPlotSpecificEvents() {
        await setupScanpyOptions();
    }

}

/**
 * Represents a SvgHandler, a class that handles SVG plots.
 * @extends PlotHandler
 */
class SvgHandler extends PlotHandler {
    constructor() {
        super();
        this.plotType = "svg";
        this.apiPlotType = "svg";
    }

    // These do not get passed into the API call, but want to keep the same data structure for cloning display
    classElt2Prop = {
        "js-svg-low-color":"low_color"
        , "js-svg-mid-color":"mid_color"
        , "js-svg-high-color":"high_color"
    }

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {colors: {}};  // Plot config to color SVG

    /**
     * Clones the display based on the provided configuration.
     * @param {Object} config - The configuration object.
     */
    cloneDisplay(config) {
        // Props are in a "colors" dict
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config.colors[prop]);
        }

        // If a mid-level color was provided, ensure checkbox to enable it is checked (for aesthetics)
        if (config.colors["mid_color"]) {
            for (const elt of document.getElementsByClassName("js-svg-enable-mid")) {
                elt.checked = true;
            }
        }

    }

    /**
     * Creates a plot for a given dataset and gene symbol.
     * @param {string} datasetId - The ID of the dataset.
     * @returns {void}
     */
    async createPlot(datasetId) {
        let data;
        try {
            data = await fetchSvgData(datasetId, this.plotConfig)
        } catch (error) {
            return;
        }
        const plotContainer = document.getElementById("plot-container");
        plotContainer.replaceChildren();    // erase plot

        colorSVG(data, this.plotConfig["colors"]);
    }

    /**
     * Loads the plot HTML and updates the DOM elements accordingly.
     * @returns {Promise<void>} A promise that resolves once the plot HTML is loaded and the DOM elements are updated.
     */
    async loadPlotHtml() {
        document.getElementById("facet-content").classList.add("is-hidden");
        document.getElementById("selected-facets").classList.add("is-hidden");

        const prePlotOptionsElt = document.getElementById("plot-options-collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post-plot-adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/svg.html");
        postPlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/post_plot/svg.html");
    }

    /**
     * Populates the plot configuration with color values based on user input.
     */
    populatePlotConfig() {
        this.plotConfig["colors"] = {};   // Reset plot config

        this.plotConfig["colors"]["low_color"] = document.getElementById("low-color").value;
        this.plotConfig["colors"]["mid_color"] = document.getElementById("mid-color").value;
        this.plotConfig["colors"]["high_color"] = document.getElementById("high-color").value;

        // If user did not choose a mid-color, set it as null instead of to black
        if (!(document.getElementById("enable-mid-color").checked)) {
            this.plotConfig["colors"]["mid_color"] = null;
        }

    }

    /**
     * Sets up an event listener for copying parameter values.
     * @returns {Promise<void>} A promise that resolves when the event listener is set up.
     */
    async setupParamValueCopyEvent() {
        setupParamValueCopyEvent("js-svg-enable-mid");
    }

    /**
     * Sets up plot-specific events.
     */
    setupPlotSpecificEvents() {
        setupSVGOptions();
    }
}

/**
 * Adds a warning message to the plot container if any categorical series value from a series in ".js-plot-req" has more than 20 groups.
 * This warning message alerts the user that the plot may be difficult to read or render properly.
 *
 * @param {HTMLElement} plotContainer - The container element where the warning message will be added.
 */
const addOvercrowdedSeriesWarning = (plotContainer) => {
    const plotlyReqSeries = document.getElementsByClassName("js-plot-req");
    const overcrowdedSeries = [...plotlyReqSeries].filter((series) => {
        const seriesValue = series.value;
        if (!levels[seriesValue]) {
            return false;
        }
        const seriesGroups = levels[seriesValue];
        return seriesGroups.length > 20;
    });
    if (!overcrowdedSeries.length) {
        return;
    }
    const overcrowdedSeriesWarning = document.createElement("article");
    overcrowdedSeriesWarning.classList.add("message", "is-warning");
    overcrowdedSeriesWarning.id = "overcrowded-series-warning";
    overcrowdedSeriesWarning.innerHTML = `
            <div class="message-header">
                <p class="mb-0">Overcrowding Warning</p>
                <button class="delete" aria-label="delete"></button>
            </div>
            <div class="message-body">
                <strong>WARNING:</strong> One or more of the selected categorical series has more than 20 groups. This may cause the plot to be more difficult to read or render properly.
            </div>
        `;
    plotContainer.prepend(overcrowdedSeriesWarning);

    // Add event listener to delete button
    const deleteButton = document.getElementById("overcrowded-series-warning").querySelector(".delete");
    deleteButton.addEventListener("click", (event) => {
        event.target.parentElement.parentElement.remove();
    });
}

/**
 * Function to handle the selection of a gene.
 */
const chooseGene = () => {

    // Cannot plot if no gene is selected
    if (!validateGeneSelected()){
        document.getElementById("gene-s-failed").classList.remove("is-hidden");
        document.getElementById("gene-s-success").classList.add("is-hidden");
        document.getElementById("current-gene").textContent = "";
        document.getElementById("current-gene-post").textContent = "";

        for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
            plotBtn.disabled = true;
        }
        return;
    }

    document.getElementById("gene-s-failed").classList.add("is-hidden");
    document.getElementById("gene-s-success").classList.remove("is-hidden");
    // Display current selected gene
    document.getElementById("current-gene-c").classList.remove("is-hidden");
    document.getElementById("current-gene").textContent = selectedGene;
    document.getElementById("current-gene-post").textContent = selectedGene;
    // Force validationcheck to see if plot button should be enabled
    trigger(document.querySelector(".js-plot-req"), "change");
    document.getElementById("plot-options-s").click();
}

/**
 * Applies color to an SVG chart based on the provided data and plot configuration.
 * @param {Object} chartData - The data used to color the chart.
 * @param {Object} plotConfig - The configuration settings for the chart.
 */
const colorSVG = (chartData, plotConfig) => {
    // I found adding the mid color for the colorblind mode  skews the whole scheme towards the high color
    const colorblindMode = CURRENT_USER.colorblind_mode;
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
    const svg = document.getElementById("plot-container");
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

/**
 * Creates an autocomplete instance.
 * @param {string} selector - The selector for the input element.
 * @param {Array} dataSource - The data source for autocomplete suggestions.
 * @param {Object} otherAutocomplete - The linked autocomplete instance.
 * @returns {Object} - The created autocomplete instance.
 */
const createAutocomplete = (selector, dataSource, otherAutocomplete) => {
    const autoCompleteJS =  new autoComplete({
        linkedAutocomplete: otherAutocomplete,
        selector,
        placeHolder: "Enter a gene",
        data: {
            src: dataSource,
            cache: true,
        },
        resultItem: {
            class: "dropdown-item",
            highlight: true,
        },
        resultsList: {
            class: "dropdown-content",
            element: (list, data) => {
                if (data.results.length) {
                    return;
                }
                // Create "No Results" message element
                const message = document.createElement("div");
                // Add class to the created element
                message.setAttribute("class", "no_result");
                // Add message text content
                message.innerHTML = `<span>Found No Results for "${data.query}"</span>`;
                // Append message element to the results list
                list.prepend(message);
            },
            noResults: true,
        },
        events: {
            input: {
                selection: (event) => {
                    // change the input value to the selected value
                    const selection = event.detail.selection.value;
                    autoCompleteJS.input.value = selection;
                    if (otherAutocomplete) {
                        otherAutocomplete.input.value = selection;
                    }
                    if (autoCompleteJS.linkedAutocomplete) {
                        autoCompleteJS.linkedAutocomplete.input.value = selection;
                    }
                    selectedGene = selection;
                    chooseGene();
                }
            }
        }
    });
    return autoCompleteJS
}

/**
 * Creates a plot based on the specified plot type.
 * @param {string} plotType - The type of plot to create.
 * @returns {Promise<void>} - A promise that resolves when the plot is created.
 */
const curatorSpecifcCreatePlot = async (plotType) => {
    // Call API route by plot type
    if (plotlyPlots.includes(plotType)) {
        await plotStyle.createPlot(datasetId, analysisObj);

    } else if (scanpyPlots.includes(plotType)) {
        await plotStyle.createPlot(datasetId, analysisObj);

    } else if (plotType === "svg") {
        await plotStyle.createPlot(datasetId);
    } else {
        console.warn(`Plot type ${plotType} selected for plotting is not a valid type.`)
        return;
    }

}

/**
 * Callback function for curator specific dataset tree.
 * @returns {void}
 */
const curatorSpecifcDatasetTreeCallback = () => {
    document.getElementById("current-gene").textContent = "";
    document.getElementById("current-gene-post").textContent = "";
}

/**
 * Callback function for selecting a specific facet item in the curator.
 * @param {string} seriesName - The name of the series.
 */
const curatorSpecifcFacetItemSelectCallback = (seriesName) => {
    // Update the color picker in case some elements of the color series were filtered out
    if(plotStyle.plotConfig?.color_name) {
        renderColorPicker(plotStyle.plotConfig.color_name);
    }
}

/**
 * Updates the curator-specific navbar with the current page information.
 */
const curatorSpecificNavbarUpdates = () => {
	document.getElementById("page-header-label").textContent = "Single-gene Displays";
}

const curatorSpecificOnLoad = async () => {
    // pass
}

/**
 * Returns a specific plot style handler based on the given plot type.
 * @param {string} plotType - The type of plot.
 * @returns {PlotlyHandler|ScanpyHandler|SvgHandler|null} - The plot style handler.
 */
const curatorSpecificPlotStyle = (plotType) => {
    // include plotting backend options
    if (plotlyPlots.includes(plotType)) {
        return new PlotlyHandler(plotType);
    } else if (scanpyPlots.includes(plotType)) {
        return new ScanpyHandler(plotType);
    } else if (plotType === "svg") {
        return new SvgHandler();
    } else {
        return null;
    }
}

/**
 * Adjusts the plot type for the dataset curator.
 * @param {string} plotType - The original plot type.
 * @returns {string} - The adjusted plot type.
 */
const curatorSpecificPlotTypeAdjustments = (plotType) => {
    // ? Move this to class constructor to handle
    if (plotType.toLowerCase() === "tsne") {
        // Handle legacy plots
        plotType = "tsne_static";
    } else if (["tsne/umap_dynamic", "tsne_dynamic"].includes(plotType.toLowerCase())) {
        plotType = "tsne_dyna";
    }
    return plotType
}

/**
 * Updates the dataset genes for the curator.
 *
 * @param {Array<string>} geneSymbols - The gene symbols to update.
 */
const curatorSpecificUpdateDatasetGenes = (geneSymbols) => {
    const geneAutocomplete = createAutocomplete("#gene-autocomplete", geneSymbols);
    const genePostAutocomplete = createAutocomplete("#gene-autocomplete-post", geneSymbols, geneAutocomplete);

    // Set the otherAutocomplete reference for geneAutocomplete after genePostAutocomplete has been created
    geneAutocomplete.linkedAutocomplete = genePostAutocomplete;
}

/**
 * Performs specific validation checks for the curator.
 * @returns {boolean} Returns true if all validation checks pass, otherwise false.
 */
const curatorSpecificValidationChecks = () => {
    if (!validateGeneSelected()) {
        return false;
    }
    return true;
}

/**
 * Fetches Plotly data for a given dataset, analysis, plot type, and plot configuration.
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} analysis - The analysis to perform.
 * @param {string} plotType - The type of plot to create.
 * @param {object} plotConfig - The configuration options for the plot.
 * @returns {Promise<object>} - The fetched Plotly data.
 * @throws {Error} - If the data fetch fails or an error occurs.
 */
const fetchPlotlyData = async (datasetId, analysis, plotType, plotConfig)  => {
    // NOTE: gene_symbol already passed to plotConfig
    try {
        const data = await apiCallsMixin.fetchPlotlyData(datasetId, analysis, plotType, plotConfig);
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

/**
 * Fetches SVG data for a given dataset and gene symbol.
 * @param {string} datasetId - The ID of the dataset.
 * @param {object} plotConfig - The configuration options for the plot.
 * @returns {Promise<Object>} - The fetched SVG data.
 * @throws {Error} - If there is an error fetching the SVG data.
 */
const fetchSvgData = async (datasetId, plotConfig) => {
    try {
        const {gene_symbol: geneSymbol} = plotConfig;
        const data = await apiCallsMixin.fetchSvgData(datasetId, geneSymbol);
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

/**
 * Fetches the TSNE image for a given dataset, analysis, plot type, and plot configuration.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} analysis - The analysis type.
 * @param {string} plotType - The type of plot.
 * @param {object} plotConfig - The configuration for the plot.
 * @returns {Promise<object>} - The fetched data.
 * @throws {Error} - If there is an error fetching the data or creating the plot image.
 */
const fetchTsneImage = async (datasetId, analysis, plotType, plotConfig) => {
    // NOTE: gene_symbol already passed to plotConfig
    try {
        const data = await apiCallsMixin.fetchTsneImage(datasetId, analysis, plotType, plotConfig);
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


/**
 * Renders the color picker for a given series name.
 *
 * @param {string} seriesName - The name of the series.
 */
const renderColorPicker = (seriesName) => {
    const colorsContainer = document.getElementById("colors-container");
    const colorsSection = document.getElementById("colors-section");

    colorsSection.classList.add("is-hidden");
    colorsContainer.replaceChildren();
    if (!seriesName) {
        return;
    }

    if (!(catColumns.includes(seriesName))) {
        // ? Continuous series colorbar picker
        return;
    }

    const seriesNameElt = document.createElement("p");
    seriesNameElt.classList.add("has-text-weight-bold", "is-underlined");
    seriesNameElt.textContent = seriesName;
    colorsContainer.append(seriesNameElt);

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

        const groupElt = document.createElement("p");
        groupElt.classList.add("is-flex", "is-justify-content-space-between", "pr-3");

        const groupText = document.createElement("span");
        groupText.classList.add("has-text-weight-medium");
        groupText.textContent = group;

        const colorInput = document.createElement("input");
        colorInput.classList.add("js-plot-color");
        colorInput.id = `${group}-color`;
        colorInput.type = "color";
        colorInput.value = groupColor;
        colorInput.setAttribute("aria-label", "Select a color");

        groupElt.append(groupText, colorInput);
        colorsContainer.append(groupElt);
    }

    colorsSection.classList.remove("is-hidden");
}

/**
 * Sets up the options for Plotly.
 * @returns {Promise<void>} A promise that resolves when the options are set up.
 */
const setupPlotlyOptions = async () => {
    const analysisId = getAnalysisId();
    const plotType = getSelect2Value(plotTypeSelect);
    try {
        ({obs_columns: allColumns, obs_levels: levels} = await curatorApiCallsMixin.fetchH5adInfo(datasetId, analysisId));
    } catch (error) {
        document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
        return;
    }
    // Filter out values we don't want of "levels", like "colors"
    allColumns = allColumns.filter((col) => !col.includes("_colors"));
    for (const key in levels) {
        if (key.includes("_colors")) {
            delete levels[key];
        }
    }

    if (!allColumns.length) {
        document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
        createToast("No metadata columns found in dataset. Cannot create a plot. Please choose another analysis or choose another dataset.");
        return;
    }

    catColumns = Object.keys(levels);


    const difference = (arr1, arr2) => arr1.filter(x => !arr2.includes(x))
    const continuousColumns = difference(allColumns, catColumns);

    const xColumns = ["bar", "violin"].includes(plotType) ? catColumns : allColumns;
    const xUseRaw = ["bar", "violin"].includes(plotType) ? false : true;
    const yColumns = ["bar", "violin"].includes(plotType) ? continuousColumns : allColumns;
    const colorColumns = ["bar", "line", "violin"].includes(plotType) ? catColumns : allColumns;
    const colorUseRaw = ["bar", "line", "violin"].includes(plotType) ? false : true;

    // Arguments - class name, list of columns, add expression, default category
    updateSeriesOptions("js-plotly-x-axis", xColumns, xUseRaw);
    updateSeriesOptions("js-plotly-y-axis", yColumns, true, "raw_value");
    updateSeriesOptions("js-plotly-color", colorColumns, colorUseRaw);
    updateSeriesOptions("js-plotly-label", allColumns, true);
    updateSeriesOptions("js-plotly-facet-row", catColumns, false);
    updateSeriesOptions("js-plotly-facet-col", catColumns, false);

    // If plot_type is bar or line, disable the marker size options
    if (["bar", "line", "violin"].includes(plotType)) {
        for (const elt of document.getElementsByClassName("js-plotly-size")) {
            elt.disabled = true;
            elt.value = "";
        }

        for (const elt of document.getElementsByClassName("js-plotly-marker-size")) {
            elt.disabled = true;
            elt.value = "";
        }
    }

    const xAxisSeriesElts = document.getElementsByClassName("js-plotly-x-axis");
    // If x-axis is categorical, enable jitter plots
    if (["violin", "scatter", "tsne_dyna"].includes(plotType)) {
        for (const elt of xAxisSeriesElts) {
            elt.addEventListener("change", (event) => {
                const jitterElts = document.getElementsByClassName("js-plotly-add-jitter");
                if ((catColumns.includes(event.target.value))) {
                    // categorical x-axis
                    for (const jitterElt of jitterElts) {
                        jitterElt.disabled = false;
                        disableCheckboxLabel(jitterElt, false);
                    }
                } else {
                    for (const jitterElt of jitterElts) {
                        jitterElt.disabled = true;
                        jitterElt.checked = false;
                        disableCheckboxLabel(jitterElt, true);
                    }

                }

            });
        }
    }

    if (["scatter", "tsne_dyna"].includes(plotType)) {
        updateSeriesOptions("js-plotly-size", continuousColumns, true);
        // If x-axis is continuous show vline stuff, otherwise hide
        // If x-axis is categorical, enable jitter plots
        for (const elt of xAxisSeriesElts) {
            elt.addEventListener("change", (event) => {
                const vLinesContainer = document.getElementById("vlines-container")
                if ((catColumns.includes(event.target.value))) {
                    vLinesContainer.classList.add("is-hidden");
                    // Remove all but first existing vline
                    const toRemove = document.querySelectorAll(".js-plotly-vline-field:not(:first-of-type)");
                    for (const elt of toRemove) {
                        elt.remove();
                    }
                    // Clear the first vline
                    document.querySelector(".js-plotly-vline-pos").value = "";
                    document.querySelector(".js-plotly-vline-style-select").value = "solid";
                    return;
                }
                vLinesContainer.classList.remove("is-hidden");

            });
        }


        // Vertical line add and remove events
        const vLinesBody = document.getElementById("vlines-body");
        const vLineField = document.querySelector(".js-plotly-vline-field");
        document.getElementById("vline-add-btn").addEventListener("click", (event) => {
            vLinesBody.append(vLineField.cloneNode(true));
            // clear the values of the clone
            document.querySelector(".js-plotly-vline-field:last-of-type .js-plotly-vline-pos").value = "";
            document.querySelector(".js-plotly-vline-field:last-of-type .js-plotly-vline-style-select").value = "solid";
            // NOTE: Currently if original is set before cloning, values are copied to clone
            document.getElementById("vline-remove-btn").disabled = false;
        })
        document.getElementById("vline-remove-btn").addEventListener("click", (event) => {
            // Remove last vline
            const lastVLine = document.querySelector(".js-plotly-vline-field:last-of-type");
            lastVLine.remove();
            if (vLineField.length < 2) document.getElementById("vline-remove-btn").disabled = true;
        })
    }

    // If color series is selected, let user choose colors.
    const colorSeriesElts = document.getElementsByClassName("js-plotly-color");
    const colorPaletteElts = document.getElementsByClassName("js-plotly-color-palette");
    const reversePaletteElts = document.getElementsByClassName("js-plotly-reverse-palette");
    const hideLegend = document.getElementsByClassName("js-plotly-hide-legend");
    for (const elt of colorSeriesElts) {
        elt.addEventListener("change", (event) => {
            if ((catColumns.includes(event.target.value))) {
                renderColorPicker(event.target.value);
                for (const paletteElt of [...colorPaletteElts, ...reversePaletteElts]) {
                    paletteElt.disabled = true;
                }
                for (const legendElt of hideLegend) {
                    legendElt.disabled = false;
                    disableCheckboxLabel(legendElt, false);
                }
                document.getElementById("color-palette-post").disabled = true;
            } else {
                // remove color picker
                const colorsContainer = document.getElementById("colors-container");
                const colorsSection = document.getElementById("colors-section");
                colorsSection.classList.add("is-hidden");
                colorsContainer.replaceChildren();

                // Enable the color palette select
                for (const paletteElt of [...colorPaletteElts, ...reversePaletteElts]) {
                    paletteElt.disabled = false;
                }
                for (const legendElt of hideLegend) {
                    legendElt.disabled = true;
                    legendElt.checked = false;
                    disableCheckboxLabel(legendElt, true);
                }
                document.getElementById("color-palette-post").disabled = false;
            }
        })
    }


    // Certain elements trigger plot order
    const plotOrderElts = document.getElementsByClassName("js-plot-order");
    for (const elt of plotOrderElts) {
        elt.addEventListener("change", (event) => {
            const paramId = event.target.id;
            const param = paramId.replace("-series", "").replace("-post", "");
            // NOTE: continuous series will be handled in the function
            updateOrderSortable();
        });
    }

    // Ensure that the same series is not selected for both x-axis and y-axis (plot would be meaningless)
    for (const classElt of [...document.getElementsByClassName("js-plotly-x-axis"), ...document.getElementsByClassName("js-plotly-y-axis")]) {
        // if x-axis series changes, disable that option in y-axis series, and vice versa
        classElt.addEventListener("change", (event) => {
            const series = event.target.value;

            // if no series selected, do nothing
            if (!series) {
                return;
            }

            const otherClass = event.target.classList.contains("js-plotly-x-axis") ? "js-plotly-y-axis" : "js-plotly-x-axis";

            for (const otherClassElt of [...document.getElementsByClassName(otherClass)]) {
                // enable all series
                for (const opt of otherClassElt.options) {
                    opt.removeAttribute("disabled");
                }
                // disable selected series in other series
                const opt = otherClassElt.querySelector(`option[value="${series}"]`);

                // If bar/violin plot, continuous y-axis but categorical x-axis
                if (!opt) {
                    continue;
                }

                opt.setAttribute("disabled", "disabled");
                // If this option was selected, unselect it
                if (otherClassElt.value === series) {
                    otherClassElt.value = "";
                }
            }
        });
    }

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    const xSeries = document.getElementById("x-axis-series");
    const ySeries = document.getElementById("y-axis-series");
    if (xSeries.value) {
        // If value is categorical, disable min and max boundaries
        for (const elt of [...document.getElementsByClassName("js-plotly-x-min"), ...document.getElementsByClassName("js-plotly-x-max")]) {
            elt.disabled = catColumns.includes(xSeries.value)
        }
        trigger(xSeries, "change");
    }
    if (ySeries.value) {
        // If value is categorical, disable min and max boundaries
        for (const elt of [...document.getElementsByClassName("js-plotly-y-min"), ...document.getElementsByClassName("js-plotly-y-max")]) {
            elt.disabled = catColumns.includes(ySeries.value)
        }
        trigger(ySeries, "change");
    }

    // Setup the dropdown menu on the post-plot view
    const plotlyDropdown = document.getElementById("plotly-param-dropdown");
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

    const plotlyDropdownMenuItems = document.querySelectorAll("#plotly-param-dropdown .dropdown-item");
    for (const item of plotlyDropdownMenuItems) {
        item.addEventListener("click", showPostPlotlyParamSubsection);
    }

}

/**
 * Sets up the options for Scanpy analysis.
 * @returns {Promise<void>} A promise that resolves when the setup is complete.
 */
const setupScanpyOptions = async () => {
    const analysisId = getAnalysisId();
    const plotType = getSelect2Value(plotTypeSelect);
    try {
        ({obs_columns: allColumns, obs_levels: levels} = await curatorApiCallsMixin.fetchH5adInfo(datasetId, analysisId));
    } catch (error) {
        document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
        return;
    }

    // Filter out values we don't want of "levels", like "colors"
    allColumns = allColumns.filter((col) => !col.includes("_colors"));
    for (const key in levels) {
        if (key.includes("_colors")) {
            delete levels[key];
        }
    }

    if (!allColumns.length) {
        document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
        createToast("No metadata columns found in dataset. Cannot create a plot. Please choose another analysis or choose another dataset.");
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
            for (const targetElt of [...plotBySeries, ...horizontalLegend]) {
                targetElt.disabled = true;
                // If colorized legend is continuous, we cannot plot by group
                // So all dependencies need to be disabled.
                if ((catColumns.includes(event.target.value))) {
                    targetElt.disabled = false;
                    disableCheckboxLabel(targetElt, false);
                }
            }

            // The "max columns" parameter should only be disabled if the colorized legend is continuous
            if (!(catColumns.includes(event.target.value))) {
                for (const targetElt of maxColumns) {
                    targetElt.disabled = true;
                    disableCheckboxLabel(targetElt, true);
                }
            }
        });

        //
        elt.addEventListener("change", (event) => {
            // if series is empty or not categorical, remove color picker
            if ((catColumns.includes(event.target.value))) {
                renderColorPicker(event.target.value);
                return;
            }
            const colorsContainer = document.getElementById("colors-container");
            const colorsSection = document.getElementById("colors-section");
            colorsSection.classList.add("is-hidden");
            colorsContainer.replaceChildren();
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
                disableCheckboxLabel(targetElt, targetElt.disabled);
            }
            // Must be allowed to specify max columns if series value selected
            for (const targetElt of maxColumns) {
                targetElt.disabled = event.target.value ? false : true;
            }
            updateOrderSortable();

        });
    }

    // Ensure that the same series is not selected for both x-axis and y-axis (plot would be meaningless)
    for (const classElt of [...document.getElementsByClassName("js-tsne-x-axis"), ...document.getElementsByClassName("js-tsne-y-axis")]) {
        // if x-axis series changes, disable that option in y-axis series, and vice versa
        classElt.addEventListener("change", (event) => {
            const series = event.target.value;
            if (!series) {
                return;
            }

            const otherClass = event.target.classList.contains("js-tsne-x-axis") ? "js-tsne-y-axis" : "js-tsne-x-axis";

            for (const otherClassElt of [...document.getElementsByClassName(otherClass)]) {
                // enable all series
                for (const opt of otherClassElt.options) {
                    opt.removeAttribute("disabled");
                }
                // disable selected series in other series
                const opt = otherClassElt.querySelector(`option[value="${series}"]`);
                opt.setAttribute("disabled", "disabled");
                // If this option was selected, unselect it
                if (otherClassElt.value === series) {
                    otherClassElt.value = "";
                }
            }
        });
    }

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    if (document.getElementById("x-axis-series").value) {
        trigger(document.getElementById("x-axis-series"), "change");
    }
    if (document.getElementById("y-axis-series").value) {
        trigger(document.getElementById("y-axis-series"), "change");
    }

    // If override marker size is checked, enable the marker size field
    const overrideMarkerSize = document.getElementsByClassName("js-tsne-override-marker-size");
    const markerSize = document.getElementsByClassName("js-tsne-marker-size");
    for (const elt of overrideMarkerSize) {
        elt.addEventListener("change", (event) => {
            for (const targetElt of markerSize) {
                targetElt.disabled = event.target.checked ? false : true;
            }
        });
    }

    // If skip gene plot is checked, disable the expression palette, reverse palette, and two-way palette options
    for (const elt of skipGenePlot) {
        elt.addEventListener("change", (event) => {
            const colorPaletteElts = document.getElementsByClassName("js-tsne-color-palette");
            const reversePaletteElts = document.getElementsByClassName("js-tsne-reverse-palette");
            const twoWayPaletteElts = document.getElementsByClassName("js-tsne-two-way-palette");

            for (const targetElt of [...colorPaletteElts, ...reversePaletteElts, ...twoWayPaletteElts]) {
                targetElt.disabled = event.target.checked ? true : false;
            }

            for (const targetElt of [...reversePaletteElts, ...twoWayPaletteElts]) {
                targetElt.checked = false;
                disableCheckboxLabel(targetElt, targetElt.disabled);
            }
        });
    }

}

/**
 * Sets up SVG options for dataset curator.
 */
const setupSVGOptions = () => {
    const enableMidColorElts = document.getElementsByClassName("js-svg-enable-mid");
    const midColorElts = document.getElementsByClassName("js-svg-mid-color");
    const midColorFields = document.getElementsByClassName("js-mid-color-field");
    for (const elt of enableMidColorElts) {
        elt.addEventListener("change", (event) => {
            for (const field of midColorFields) {
                field.classList.add("is-hidden");
                if (event.target.checked) {
                    field.classList.remove("is-hidden");
                }
            }
            for (const midColor of midColorElts) {
                midColor.disabled = !(event.target.checked);
            }
        });
    }


    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    if (document.getElementById("low-color").value) {
        trigger(document.getElementById("low-color"), "change");
    }
    if (document.getElementById("high-color").value) {
        trigger(document.getElementById("high-color"), "change");
    }
}

/**
 * Shows the corresponding subsection based on the selected option in the plot configuration menu.
 * @param {Event} event - The event triggered by the user's selection.
 */
const showPostPlotlyParamSubsection = (event) => {
    for (const subsection of document.getElementsByClassName("js-plot-config-section")) {
        subsection.classList.add("is-hidden");
    }

    switch (event.target.textContent.trim()) {
        case "X-axis":
            document.getElementById("x-axis-section-post").classList.remove("is-hidden");
            break;
        case "Y-axis":
            document.getElementById("y-axis-section-post").classList.remove("is-hidden");
            break;
        case "Color":
            document.getElementById("color-section-post").classList.remove("is-hidden");
            break;
        case "Marker Size":
            document.getElementById("size-section-post").classList.remove("is-hidden");
            break;
        case "Subplots":
            document.getElementById("subplots-section-post").classList.remove("is-hidden");
            break;
        default:
            document.getElementById("misc-section-post").classList.remove("is-hidden");
            break;
    }
    event.preventDefault(); // Prevent "link" clicking from "a" elements
}

/**
 * Updates the series options in a select element based on the provided parameters.
 *
 * @param {string} classSelector - The class selector for the select elements to update.
 * @param {Array<string>} seriesArray - An array of series names.
 * @param {boolean} addExpression - Indicates whether to add an expression option.
 * @param {string} defaultOption - The default option to select.
 */
const updateSeriesOptions = (classSelector, seriesArray, addExpression, defaultOption) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Create continuous and categorical optgroups
        const contOptgroup = document.createElement("optgroup");
        contOptgroup.setAttribute("label", "Continuous data");
        const catOptgroup = document.createElement("optgroup");
        catOptgroup.setAttribute("label", "Categorical data");

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add an expression option (since expression is not in the categories)
        if (addExpression) {
            const expression = document.createElement("option");
            contOptgroup.append(expression);
            expression.textContent = "expression";
            expression.value = "raw_value";
            if ("raw_value" === defaultOption) {
                expression.selected = true;
            }
        }

        // Add categories
        for (const group of seriesArray.sort()) {

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
            if (catColumns.includes(group)) {
                catOptgroup.append(option);
            } else {
                contOptgroup.append(option);
            }
            // NOTE: It is possible for a default option to not be in the list of groups.
            if (group === defaultOption) {
                option.selected = true;
            }
        }

        // Only append optgroup if it has children
        if (contOptgroup.children.length) elt.append(contOptgroup);
        if (catOptgroup.children.length) elt.append(catOptgroup);

    }
}

/**
 * Validates the selected gene.
 *
 * @returns {boolean} Returns true if a gene is selected, false otherwise.
 */
const validateGeneSelected = () => {
    return (selectedGene ? true : false);
}
