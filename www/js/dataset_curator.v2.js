// I use camelCase for my variable/function names to adhere to JS style standards
// Exception being functions that do fetch calls, so we can use JS destructuring on the payload

'use strict';

const isMultigene = 0;

let geneSelect = null;

const plotlyPlots = ["bar", "line", "scatter", "tsne_dyna", "violin"];
const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];

let geneSelectPost = null;

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
        , "js-plotly-color-palette":"color_palette"
        , "js-plotly-reverse-palette":"reverse_palette"
    }

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    cloneDisplay(config) {
        // plotly plots
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            for (const series in config["order"]) {
                renderOrderSortableSeries(series);
            }

            document.getElementById("order_section").style.display = "";
        }

        // Handle colors
        if (config["colors"]) {
            // do nothing if color_name is not set
            if (!config["color_name"]) return;

            // BUG: Occasionally the color_name element is not populated yet, which prevents the series from being rendered
            try {
                const series = config["color_name"];
                renderColorPicker(series);
                for (const group in config["colors"]) {
                    const color = config["colors"][group];
                    const colorField = document.getElementById(`${group}_color`);
                    colorField.value = color;
                }
            } catch (error) {
                console.error(error);
                // pass
            }
        }

        if (config["color_palette"]) {
            setSelectBoxByValue("color_palette_post", config["color_palette"]);
            colorscaleSelect.update();
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
            return;
        }

        const plotContainer = document.getElementById("plot_container");
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const divElt = generateElements('<div class="container is-max-desktop" id="plotly_preview"></div>');
        plotContainer.append(divElt);
        Plotly.purge("plotly_preview"); // clear old Plotly plots

        if (!plotJson) {
            createToast("Could not retrieve plot information. Cannot make plot.");
            return;
        }
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "config");
        Plotly.newPlot("plotly_preview", plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "layout")
        Plotly.relayout("plotly_preview", custonLayout)

        // If any categorical series in ".js_plot_req", and the series has more then 20 groups, display a warning about overcrowding
        const plotlyReqSeries = document.getElementsByClassName("js_plot_req");
        const overcrowdedSeries = [...plotlyReqSeries].filter((series) => {
            const seriesName = series.id.replace("_color", "");
            const seriesGroups = levels[seriesName];
            return seriesGroups.length > 20;
        });
        if (overcrowdedSeries.length) {
            const template = `
                <article class="message is-warning" id="overcrowded_series_warning">
                    <div class="message-body">
                        <strong>WARNING:</strong> One or more of the selected categorical series has more than 20 groups. This may cause the plot to be more difficult to read or render properly.
                    </div>
                </article>
            `
            const elt = generateElements(template);
            plotContainer.prepend(elt);

            // Add event listener to delete button
            const deleteButton = document.getElementById("overcrowded_series_warning").querySelector(".delete");
            deleteButton.addEventListener("click", (event) => {
                event.target.parentElement.parentElement.remove();
            });
        }

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

        // Load color palette select options
        loadColorscaleSelect(true);

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

    async setupPlotSpecificEvents() {
        await setupPlotlyOptions();
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

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    cloneDisplay(config) {
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            for (const series in config["order"]) {
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
                    const colorField = document.getElementById(`${group}_color`);
                    colorField.value = color;
                }
            }
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
    }

    async createPlot(datasetId, analysisObj, userId, colorblindMode) {
        let image;
        try {
            const data = await fetchTsneImageData(this.plotConfig, datasetId, this.plotType, analysisObj, userId, colorblindMode);
            ({image} = data);
        } catch (error) {
            return;
        }

        const plotContainer = document.getElementById("plot_container");
        plotContainer.replaceChildren();    // erase plot

        const imgElt = generateElements('<img class="image" id="tsne_preview"></img>');
        plotContainer.append(imgElt);

        if (image) {
            document.getElementById("tsne_preview").setAttribute("src", `data:image/png;base64,${image}`);
        } else {
            createToast("Could not retrieve plot image. Cannot make plot.");
            return;
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

    async setupPlotSpecificEvents() {
        await setupScanpyOptions();
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

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {colors: {}};  // Plot config to color SVG

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

    async createPlot(datasetId, geneSymbol) {
        let data;
        try {
            data = await fetchSvgData(datasetId, geneSymbol)
        } catch (error) {
            return;
        }
        const plotContainer = document.getElementById("plot_container");
        plotContainer.replaceChildren();    // erase plot

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

const cloneDisplay = async (event, display) => {

    const cloneId = event.target.id;
    document.getElementById(cloneId).classList.add("is-loading");

    // Populate gene select element
    // Will be overwritten if an analysis was in config
    try {
        const geneSymbols = await fetchGeneSymbols(datasetId, null);
        updateGeneOptions(geneSymbols);
    } catch (error) {
        document.getElementById("gene_s_failed").style.display = "";
    }

    document.getElementById("analysis_select").disabled = false;
    document.getElementById("plot_type_select").disabled = false;

    // analyses
    await analysisSelectUpdate();
    // TODO update analysis select with config analysis


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
        plotTypeSelect.update();
        await choosePlotType();
        // In this step, a PlotStyle object is instantiated onto "plotStyle", and we will use that
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

    // Choose gene from config
    const config = display.plotly_config;
    setSelectBoxByValue("gene_select", config.gene_symbol);
    geneSelect.update();
    trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to set the other select2

    plotStyle.cloneDisplay(config);

    // Mark plot params as success
    document.getElementById("plot_options_s_success").style.display = "";

    // Click "submit" button to load plot
    document.getElementById("plot_btn").click();    // updates geneSelectPost by triggered "click" event

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

const curatorSpecifcChooseGene = (event) => {
    // If one select element was updated ensure the other is updated as well
    const select2 = event.target.id === "gene_select" ? geneSelect : geneSelectPost;
    const oppSelect2 = event.target.id === "gene_select" ? geneSelectPost : geneSelect;
    const oppEltId = event.target.id === "gene_select" ? "gene_select_post" : "gene_select";

    if (!select2.selectedOptions.length) return;   // Do not trigger after initial population

    const val = getSelect2Value(select2);

    // NOTE: I thought about updating the select2 element directly with updateSelectValue()
    // but this triggers the "change" event for the regular "select" element, which causes a max stack call error
    setSelectBoxByValue(oppEltId, val);

    // copy data from one select2 to the other
    // Render the dropdown for the other select2
    oppSelect2.data = select2.data;
    oppSelect2.options = select2.options;
    oppSelect2.selectedOptions = select2.selectedOptions;

    // Recreate update() function without the extractData() call, which is causing noticeable slowdown/hanging
    if (oppSelect2.dropdown) {
        const open = oppSelect2.dropdown.classList.contains("open");
        oppSelect2.dropdown.parentNode.removeChild(oppSelect2.dropdown);
        oppSelect2.create();

        if (open) {
        triggerClick(oppSelect2.dropdown);
        }
    }

    // Cannot plot if no gene is selected
    if (document.getElementById("gene_select").value === "Please select a gene") {
        document.getElementById("gene_s_failed").style.display = "";
        document.getElementById("gene_s_success").style.display = "none";
        document.getElementById("current_gene").textContent = "";
        for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
            plotBtn.disabled = true;
        }
        return;
    }

    document.getElementById("gene_s_failed").style.display = "none";
    document.getElementById("gene_s_success").style.display = "";
    // Display current selected gene
    document.getElementById("current_gene_c").style.display = "";
    document.getElementById("current_gene").textContent = val;
    // Force validationcheck to see if plot button should be enabled
    trigger(document.querySelector(".js-plot-req"), "change");
    document.getElementById("plot_options_s").click();
}

const curatorSpecifcCreatePlot = async (plotType) => {
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

}

const curatorSpecifcDatasetTreeCallback = () => {

    // Create gene select2 elements for both views
    geneSelect = createGeneSelectInstance(geneSelect);
    geneSelectPost = createGeneSelectInstance(geneSelectPost);

}

const curatorSpecificOnLoad = async () => {
    // pass
}

const curatorSpecificUpdateGeneOptions = (geneSymbols) => {
    // copy to "#gene_select_post"
    const geneSelectEltPost = document.getElementById("gene_select_post");
    geneSelectEltPost.replaceChildren();
    for (const gene of geneSymbols.sort()) {
        const option = document.createElement("option");
        option.textContent = gene;
        option.value = gene;
        geneSelectEltPost.append(option);
    }

}

const deleteDisplay = async(user_id, displayId) => {
    const payload = {user_id, id: displayId};
    try {
        await axios.post("/cgi/delete_dataset_display.cgi", convertToFormData(payload));
        // Remove display card
        const displayCard = document.getElementById(`${displayId}_display`);
        displayCard.remove();
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not delete this display. Please contact the gEAR team."
        createToast(msg);
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

    if (classElts.length) {
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

    // NOTE: Changing plots within the same plot style will clear the plot config as fresh templates are loaded
    await plotStyle.loadPlotHtml();

    // NOTE: Events are triggered in the order they are regstered.
    // We want to trigger a param sync event before the plot requirement validation
    // (since it checks every element in the js-plot-req class)
    setupParamValueCopyEvents(Object.keys(plotStyle.classElt2Prop));    // Ensure pre- and post- plot view params are synced up
    setupValidationEvents();        // Set up validation events required to plot at a minimum
    await plotStyle.setupPlotSpecificEvents()       // Set up plot-specific events

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
            <p class="is-flex is-justify-content-space-between pr-3">
                <span class="has-text-weight-medium">${group}</span>
                <input class="js-plot-color" id="${group}_color" type="color" value="${groupColor}" aria-label="Select a color" />
            </p>
        `);
        colorsContainer.append(groupHtml)
    }

    colorsSection.style.display = "";
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

/* Set HTML element value from the plot config value */
const setPlotEltValueFromConfig = (classSelector, confVal) => {
    for (const elt of document.getElementsByClassName(classSelector)) {
        if (elt.type === "checkbox") {
            elt.checked = confVal;
            continue;
        }
        elt.value = confVal;
        trigger(elt, "change");
    }
}

/**
 * Set Select Box Selection By Value
 * Modified to set value and "selected" so nice-select2 extractData() will catch it
 * Taken from https://stackoverflow.com/a/20662180
 * @param eid Element ID
 * @param eval Element value
 */
const setSelectBoxByValue = (eid, val) => {
    const elt = document.getElementById(eid);
    for (const i in elt.options) {
        if (elt.options[i].value === val) {
            elt.value = val;
            elt.options[i].setAttribute("selected", true);
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
                    disableCheckboxLabel(classElt, classElt.disabled);
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

    const difference = (arr1, arr2) => arr1.filter(x => !arr2.includes(x))
    const continuousColumns = difference(allColumns, catColumns);

    // class name, list of columns, add expression, default category

    const xColumns = ["bar", "violin"].includes(plotType) ? catColumns : allColumns;
    const xUseRaw = ["bar", "violin"].includes(plotType) ? false : true;
    const yColumns = ["bar", "violin"].includes(plotType) ? continuousColumns : allColumns;

    updateSeriesOptions("js-plotly-x-axis", xColumns, xUseRaw);
    updateSeriesOptions("js-plotly-y-axis", yColumns, true, "raw_value");
    updateSeriesOptions("js-plotly-color", allColumns, true);
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
                const vLinesContainer = document.getElementById("vlines_container")
                if ((catColumns.includes(event.target.value))) {
                    vLinesContainer.style.display = "none";
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
                vLinesContainer.style.display = "";

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
                colorscaleSelect.disable()
            } else {
                // Enable the color palette select
                for (const paletteElt of [...colorPaletteElts, ...reversePaletteElts]) {
                    paletteElt.disabled = false;
                }
                for (const legendElt of hideLegend) {
                    legendElt.disabled = true;
                    legendElt.checked = false;
                    disableCheckboxLabel(legendElt, true);
                }
                colorscaleSelect.enable()
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
    const xSeries = document.getElementById("x_axis_series");
    const ySeries = document.getElementById("y_axis_series");
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
            for (const targetElt of [...plotBySeries, ...horizontalLegend]) {
                targetElt.disabled = true;
                // If colorized legend is continuous, we cannot plot by group
                // So all dependencies need to be disabled.
                if ((catColumns.includes(event.target.value))) {
                    renderColorPicker(event.target.value);
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
                disableCheckboxLabel(targetElt, targetElt.disabled);
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

// Update the params that will comprise the "order" section in post-plot view
const updateOrderSortable = () => {
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