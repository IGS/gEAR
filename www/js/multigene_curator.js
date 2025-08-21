'use strict';

import { apiCallsMixin, getCurrentUser, initCommonUI, logErrorInConsole, trigger } from './common.v2.js';
import { curatorCommon } from './curator_common.js';
import { Gene, WeightedGene } from "./classes/gene.js";
import { GeneCart, WeightedGeneCart } from './classes/genecart.v2.js';
import { adjustStackedViolinHeight, postPlotlyConfig, setHeatmapHeightBasedOnGenes } from './plot_display_config.js';
import { fetchGeneCartData, geneCollectionState } from '../include/gene-collection-selector/gene-collection-selector.js';

// Pre-initialize some stuff
initCommonUI();

curatorCommon.setIsMultigene(1);

let allColumns = [];

const notFoundGenes = new Set(); // Store genes that are not found in datasetGenes

let manuallyEnteredGenes = new Set();
let datasetGenes = new Set();

let plotSelectedGenes = []; // genes selected from plot "select" utility

const genesAsAxisPlots = ["dotplot", "heatmap", "mg_violin"];
const genesAsDataPlots = ["quadrant", "volcano"];
const scanpyPlots = ["mg_pca_static", "mg_tsne_static", "mg_umap_static"];

class GenesAsAxisHandler extends curatorCommon.PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
        this.apiPlotType = plotType;
    }

    classElt2Prop = {
        "js-dash-primary": "primary_col"
        , "js-dash-secondary": "secondary_col"
        , "js-dash-color-palette": "colorscale"
        , "js-dash-reverse-palette": "reverse_colorscale"
        , "js-dash-distance-metric": "distance_metric"
        , "js-dash-matrixplot": "matrixplot"
        , "js-dash-center-mean": "center_around_zero"
        , "js-dash-cluster-obs": "cluster_obs"
        , "js-dash-cluster-genes": "cluster_genes"
        , "js-dash-flip-axes": "flip_axes"
        , "js-dash-hide-obs-labels": "hide_obs_labels"
        , "js-dash-hide-gene-labels": "hide_gene_labels"
        , "js-dash-add-jitter": "violin_add_points"
        , "js-dash-subsampling-limit": "subsample_limit"
        , "js-dash-stacked-violin": "stacked_violin"
        , "js-dash-plot-title": "plot_title"
        , "js-dash-legend-title": "legend_title"
    };

    // deal with clusterbar separately since it is an array of selected options

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    plotJson = null;  // Plotly plot JSON

    cloneDisplay(config) {
        // load plot values
        for (const prop in config) {
            curatorCommon.setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["sort_order"]) {
            for (const series in config["sort_order"]) {
                const order = config["sort_order"][series];
                const levels = curatorCommon.getLevels();
                // sort "levels" series by order
                levels[series].sort((a, b) => order.indexOf(a) - order.indexOf(b));
                curatorCommon.renderOrderSortableSeries(series);
            }

            document.getElementById("order-section").classList.remove("is-hidden");
        }

        // Handle filters
        if (config["obs_filters"]) {
            curatorCommon.getFacetWidget().filters = config["obs_filters"];
        }

        // handle colors
        if (config["color_palette"]) {
            curatorCommon.setSelectBoxByValue("color-palette-post", config["color_palette"]);
        }

        // handle clusterbar values
        if (config["clusterbar_fields"]) {
            for (const field of config["clusterbar_fields"]) {
                const elt = document.querySelector(`#clusterbar-c .js-dash-clusterbar-checkbox[value="${field}"]`);
                elt.checked = true;
            }
        }

        // restore subsample limit
        if (config["subsample_limit"]) {
            for (const classElt of document.getElementsByClassName("js-dash-enable-subsampling")) {
                classElt.checked = true;
            }
        }

    }

    async createPlot(datasetId, analysisObj) {
        // Get data and set up the image area
        try {
            const data = await fetchMgPlotlyData(datasetId, analysisObj, this.apiPlotType, this.plotConfig);
            ({ plot_json: this.plotJson } = data);
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

        if (!this.plotJson) {
            createToast("Could not retrieve plot information. Cannot make plot.");
            return;
        }

        if (this.plotType === 'heatmap') {
            setHeatmapHeightBasedOnGenes(this.plotJson.layout, this.plotConfig.gene_symbols);
        } else if (this.plotType === "mg_violin" && this.plotConfig.stacked_violin) {
            adjustStackedViolinHeight(this.plotJson.layout);
        }

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = curatorCommon.getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "config");
        Plotly.newPlot("plotly-preview", this.plotJson.data, this.plotJson.layout, custonConfig);   // HIGH MEM/CPU with heatmap no matrixplot
        const custonLayout = curatorCommon.getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "layout");
        Plotly.relayout("plotly-preview", custonLayout);

        document.getElementById("legend-title-container").classList.remove("is-hidden");
        if (this.plotType === "dotplot") {
            document.getElementById("legend-title-container").classList.add("is-hidden");
        }

    }

    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot-options-collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post-plot-adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/multi_gene_as_axis.html");
        postPlotOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/multi_gene_as_axis.html");

        // populate advanced options for specific plot types
        const prePlotSpecificOptionsElt = document.getElementById("plot-specific-options");
        const postPlotSpecificOptionselt = document.getElementById("post-plot-specific-options");

        // Load color palette select options
        const isContinuous = ["dotplot", "heatmap"].includes(this.plotType) ? true : false;
        curatorCommon.loadColorscaleSelect(isContinuous);


        if (this.plotType === "heatmap") {
            prePlotSpecificOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/advanced_heatmap.html");
            postPlotSpecificOptionselt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/advanced_heatmap.html");
            return;
        }
        if (this.plotType === "mg_violin") {
            prePlotSpecificOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/advanced_mg_violin.html");
            postPlotSpecificOptionselt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/advanced_mg_violin.html");
            return;
        }
    }

    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = curatorCommon.getPlotConfigValueFromClassName(classElt);
        }

        // Get checked clusterbar values
        if (this.plotType === "heatmap") {
            const clusterbarValues = [];
            // They should be synced so just grab the first set of clusterbar values
            for (const elt of document.querySelectorAll("#clusterbar-c .js-dash-clusterbar-checkbox")) {
                if (elt.checked) {
                    clusterbarValues.push(elt.value);
                }
            }
            this.plotConfig["clusterbar_fields"] = clusterbarValues;


            // if subsampling is checked, add subsampling limit to plot config
            this.plotConfig["subsample_limit"] = 0;
            if (document.querySelector(".js-dash-enable-subsampling").checked) {
                this.plotConfig["subsample_limit"] = Number(document.querySelector(".js-dash-subsampling-limit").value);
            }
        }

        // Filtered observation groups
        this.plotConfig["obs_filters"] = curatorCommon.getFacetWidget()?.filters || {};
        if (Object.keys(this.plotConfig["obs_filters"]).length === 0) {
            this.plotConfig["obs_filters"] = null; // don't send empty filters
        }

        // Get order
        this.plotConfig["sort_order"] = curatorCommon.getPlotOrderFromSortable();
        if (!curatorCommon.getSortOrderChanged()) {
            // If no order was changed, set to null so API does not try to sort by the default order
            this.plotConfig["sort_order"] = null;
        }

    }

    async setupParamValueCopyEvent() {
        // These plot parameters do not directly correlate to a plot config property
        //setupParamValueCopyEvent("js-dash-enable-subsampling")
    }

    async setupPlotSpecificEvents(datasetId) {
        const catColumns = await getCategoryColumns(datasetId);
        curatorCommon.setCatColumns(catColumns);

        if (!catColumns.length) {
            document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
            createToast("No metadata columns found in dataset. Cannot create a plot. Please choose another analysis or choose another dataset.");
            return;
        }

        updateSeriesOptions("js-dash-primary", catColumns);
        updateSeriesOptions("js-dash-secondary", catColumns);

        // If primary series changes, disable chosen option in secondary series
        for (const classElt of document.getElementsByClassName("js-dash-primary")) {
            classElt.addEventListener("change", (event) => {
                const primarySeries = event.target.value;
                for (const secondaryClassElt of document.getElementsByClassName("js-dash-secondary")) {
                    // enable all series
                    for (const opt of secondaryClassElt.options) {
                        opt.removeAttribute("disabled");
                    }
                    // disable primary series in secondary series
                    const opt = secondaryClassElt.querySelector(`option[value="${primarySeries}"]`);
                    if (opt) {
                        opt.setAttribute("disabled", "disabled");
                        // If this option was selected, unselect it
                        if (secondaryClassElt.value === primarySeries) {
                            secondaryClassElt.value = "";
                        }
                    }
                }
            });
        }

        // if subsampling is checked, enable subsampling limit
        for (const classElt of document.getElementsByClassName("js-dash-enable-subsampling")) {
            classElt.addEventListener("change", (event) => {
                const checked = event.target.checked;
                for (const innerClassElt of document.getElementsByClassName("js-dash-subsampling-limit")) {
                    innerClassElt.disabled = !(checked);
                }
            });
        }

        // For clusterbar options, create checkboxes for all catColumns
        for (const classElt of document.getElementsByClassName("js-dash-clusterbar")) {
            for (const catColumn of catColumns) {

                const clusterbarElt = document.createElement("div");
                clusterbarElt.classList.add("control");

                const label = document.createElement("label");
                label.classList.add("checkbox");

                const input = document.createElement("input");
                input.type = "checkbox";
                input.value = catColumn;
                input.disabled = true;
                input.classList.add("js-dash-clusterbar-checkbox");
                label.appendChild(input);
                label.innerHTML += ` ${catColumn}`;
                clusterbarElt.appendChild(label);

                classElt.appendChild(clusterbarElt);

            }
        }

        // Add event listener to sync checkboxes with same value
        for (const inputElt of document.getElementsByClassName("js-dash-clusterbar-checkbox")) {
            inputElt.addEventListener("change", (event) => {
                const checked = event.target.checked;
                const value = event.target.value;
                for (const innerClassElt of document.getElementsByClassName("js-dash-clusterbar-checkbox")) {
                    if (innerClassElt.value === value) {
                        innerClassElt.checked = checked;
                    }
                }
            });
        }

        // If stacked violin is checked, disable legend title (since there is no legend)
        for (const classElt of document.getElementsByClassName("js-dash-stacked-violin")) {
            classElt.addEventListener("change", (event) => {
                const checked = event.target.checked;
                for (const legendTitleElt of document.getElementsByClassName("js-dash-legend-title")) {
                    if (checked) {
                        legendTitleElt.setAttribute("disabled", "disabled");
                    } else {
                        legendTitleElt.removeAttribute("disabled");
                    }
                }
            });
        }

        // Certain elements trigger plot order
        const plotOrderElts = document.getElementsByClassName("js-plot-order");
        for (const elt of plotOrderElts) {
            elt.addEventListener("change", (event) => {
                const paramId = event.target.id;
                const param = paramId.replace("-series", "").replace("-post", "");
                // NOTE: continuous series will be handled in the function
                curatorCommon.updateOrderSortable();
            });
        }

        // Setup the dropdown menu on for heatmap options
        if (this.plotType === "heatmap") {
            const heatmapDropdown = document.getElementById("heatmap-param-dropdown");
            heatmapDropdown.addEventListener("click", (event) => {
                event.stopPropagation();    // This prevents the document from being clicked as well.
                heatmapDropdown.classList.toggle("is-active");
            });

            // Close dropdown if it is clicked off of, or ESC is pressed
            // https://siongui.github.io/2018/01/19/bulma-dropdown-with-javascript/#footnote-1
            document.addEventListener('click', () => {
                heatmapDropdown.classList.remove("is-active");
            });
            document.addEventListener('keydown', (event) => {
                if (event.key === "Escape") {
                    heatmapDropdown.classList.remove("is-active");
                }
            });

            const heatmapDropdownMenuItems = document.querySelectorAll("#heatmap-param-dropdown .dropdown-item");
            for (const item of heatmapDropdownMenuItems) {
                item.addEventListener("click", showPostHeatmapParamSubsection);
            }

            // If matrixplot is unchecked enable all clusterbars.
            for (const classElt of document.getElementsByClassName("js-dash-matrixplot")) {
                classElt.addEventListener("change", (event) => {
                    const checked = event.target.checked;
                    if (checked) {
                        // Get primary and secondary series
                        const primarySeries = document.querySelector(".js-dash-primary").value;
                        const secondarySeries = document.querySelector(".js-dash-secondary").value;
                        const values = [];
                        if (primarySeries) {
                            values.push(primarySeries);
                        }
                        if (secondarySeries) {
                            values.push(secondarySeries);
                        }
                        // If either series are selected, enable those checkboxes
                        for (const innerClassElt of document.getElementsByClassName("js-dash-clusterbar-checkbox")) {
                            if (values.includes(innerClassElt.value)) {
                                innerClassElt.disabled = false;
                            } else {
                                innerClassElt.checked = false;
                                innerClassElt.disabled = true;
                            }
                        }
                    } else {
                        for (const innerClassElt of document.getElementsByClassName("js-dash-clusterbar-checkbox")) {
                            innerClassElt.disabled = false;
                        }
                    }
                });
            }


            // Ensure only clusterbar checkboxes that match the primary and secondary series are enabled to be checked
            // If a checkbox is checked and the primary or secondary series is changed, uncheck the checkbox
            for (const classElt of document.getElementsByClassName("js-dash-primary")) {
                classElt.addEventListener("change", (event) => {
                    const primarySeries = event.target.value;
                    const secondarySeries = document.querySelector(".js-dash-secondary").value;
                    const values = [];
                    if (primarySeries) {
                        values.push(primarySeries);
                    }
                    if (secondarySeries) {
                        values.push(secondarySeries);
                    }

                    // Get matrixplot value
                    const matrixplot = document.querySelector(".js-dash-matrixplot").checked;

                    if (matrixplot) {
                        for (const innerClassElt of document.getElementsByClassName("js-dash-clusterbar-checkbox")) {
                            if (values.includes(innerClassElt.value)) {
                                innerClassElt.disabled = false;
                            } else {
                                innerClassElt.checked = false;
                                innerClassElt.disabled = true;
                            }
                        }
                    }
                });
            }
            for (const classElt of document.getElementsByClassName("js-dash-secondary")) {
                classElt.addEventListener("change", (event) => {
                    const secondarySeries = event.target.value;
                    const primarySeries = document.querySelector(".js-dash-primary").value;
                    const values = [];
                    if (primarySeries) {
                        values.push(primarySeries);
                    }
                    if (secondarySeries) {
                        values.push(secondarySeries);
                    }

                    // Get matrixplot value
                    const matrixplot = document.querySelector(".js-dash-matrixplot").checked;

                    if (matrixplot) {
                        for (const innerClassElt of document.getElementsByClassName("js-dash-clusterbar-checkbox")) {
                            if (values.includes(innerClassElt.value)) {
                                innerClassElt.disabled = false;
                            } else {
                                innerClassElt.checked = false;
                                innerClassElt.disabled = true;
                            }
                        }
                    }
                });
            }


        }

    }

}

class GenesAsDataHandler extends curatorCommon.PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
        this.apiPlotType = plotType;
    }

    compareSeparator = ";-;";

    classElt2Prop = {
        "js-dash-de-test": "de_test_algo"
        , "js-dash-fold-change-cutoff": "fold_change_cutoff"
        , "js-dash-fdu-cutoff": "fdr_cutoff"
        , "js-dash-include-zero-fc": "include_zero_fc"
        , "js-dash-annot-nonsig": "annot_nonsignificant"
        , "js-dash-pvalue-threshold": "pvalue_threshold"
        , "js-dash-use-adj-pvalues": "adj_pvals"
        , "js-dash-lower-logfc-threshold": "lower_logfc_threshold"
        , "js-dash-upper-logfc-threshold": "upper_logfc_threshold"
        , "js-dash-plot-title": "plot_title"
        , "js-dash-legend-title": "legend_title"
    };
    // Deal with js-dash-query/reference/compare1/compare2 separately since they combine with js-dash-compare

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    plotJson = null;  // Plotly plot JSON

    cloneDisplay(config) {
        // load plot values
        for (const prop in config) {
            curatorCommon.setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle filters
        if (config["obs_filters"]) {
            curatorCommon.getFacetWidget().filters = config["obs_filters"];
        }

        // Split compare series and groups
        const refCondition = config["ref_condition"];
        const [compareSeries, refGroup] = refCondition.split(this.compareSeparator);
        for (const classElt of document.getElementsByClassName("js-dash-compare")) {
            classElt.value = compareSeries;
        }

        const levels = curatorCommon.getLevels();
        // populate group options
        updateGroupOptions("js-dash-reference", levels[compareSeries]);
        for (const classElt of document.getElementsByClassName("js-dash-reference")) {
            classElt.value = refGroup;
        }

        if (this.plotType === "volcano") {
            updateGroupOptions("js-dash-query", levels[compareSeries]);
            const queryCondition = config["query_condition"];
            const queryGroup = queryCondition.split(this.compareSeparator)[1];
            for (const classElt of document.getElementsByClassName("js-dash-query")) {
                classElt.value = queryGroup;
            }
            trigger(document.querySelector(".js-dash-query"), "change"); // trigger change event to start validation

        }
        if (this.plotType === "quadrant") {
            updateGroupOptions("js-dash-compare1", levels[compareSeries]);
            updateGroupOptions("js-dash-compare2", levels[compareSeries]);

            const compare1Condition = config["compare1_condition"];
            const compare1Group = compare1Condition.split(this.compareSeparator)[1];
            const compare2Condition = config["compare2_condition"];
            const compare2Group = compare2Condition.split(this.compareSeparator)[1];
            for (const classElt of document.getElementsByClassName("js-dash-compare1")) {
                classElt.value = compare1Group;
            }
            for (const classElt of document.getElementsByClassName("js-dash-compare2")) {
                classElt.value = compare2Group;
            }

            trigger(document.querySelector(".js-dash-compare1"), "change"); // trigger change event to start validation
        }

        // for some reason triggering .js-dash-compare did not populate the compare groups into the plot config
    }

    async createPlot(datasetId, analysisObj) {
        // Get data and set up the image area
        try {
            const data = await fetchMgPlotlyData(datasetId, analysisObj, this.apiPlotType, this.plotConfig);
            ({ plot_json: this.plotJson } = data);
        } catch (error) {
            return;
        }

        const plotContainer = document.getElementById("plot-container");
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const plotlyPreviewElt = document.createElement("div");
        plotlyPreviewElt.classList.add("container", "is-max-desktop");
        plotlyPreviewElt.id = "plotly-preview";
        plotContainer.append(plotlyPreviewElt);
        Plotly.purge("plotly-preview"); // clear old Plotly plots

        if (!this.plotJson) {
            createToast("Could not retrieve plot information. Cannot make plot.");
            return;
        }
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = curatorCommon.getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "config");
        Plotly.newPlot("plotly-preview", this.plotJson.data, this.plotJson.layout, custonConfig);
        const custonLayout = curatorCommon.getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "layout");
        Plotly.relayout("plotly-preview", custonLayout);

        // Show button to add genes to gene cart
        document.getElementById("gene-cart-btn-c").classList.remove("is-hidden");

        const plotlyPreview = document.getElementById("plotly-preview");

        // Append small note about using the Plotly selection utilities
        const plotlyNote = document.createElement("div");
        plotlyNote.classList.add("notification", "is-info", "is-light");
        plotlyNote.innerHTML = `<p><strong>Tip:</strong> Use the Plotly box and lasso select tools (upper-right) to select genes to view as a table.</p>`;
        plotlyPreview.append(plotlyNote);

        // If plot data is selected, create the right-column table and do other misc things
        plotlyPreview.on("plotly_selected", async (eventData) => {

            // Hide selected genes table and disable unweighted radio button if no genes are selected
            document.getElementById("tbl-selected-genes").classList.add("is-hidden");
            document.getElementById("download-selected-genes-btn").classList.add("is-hidden");
            document.querySelector("input[name='genecart_type'][value='unweighted']").disabled = true;
            document.querySelector("input[name='genecart_type'][value='unweighted']").parentElement.removeAttribute("disabled");
            if (eventData?.points.length) {
                document.getElementById("tbl-selected-genes").classList.remove("is-hidden");
                document.getElementById("download-selected-genes-btn").classList.remove("is-hidden");
                document.querySelector("input[name='genecart_type'][value='unweighted']").disabled = false;
                document.querySelector("input[name='genecart_type'][value='unweighted']").parentElement.setAttribute("disabled", "disabled");

                adjustGeneTableLabels(this.plotType);
                populateGeneTable(eventData, this.plotType);
            }

            // Highlight table rows that match searched genes
            const searchedGenes = this.plotConfig.gene_symbols;
            if (searchedGenes) {
                const geneTableBody = document.getElementById("gene-table-body");
                // Select the first column (gene_symbols) in each row
                for (const row of geneTableBody.children) {
                    const tableGene = row.children[0].textContent;
                    for (const gene of searchedGenes) {
                        if (gene.toLowerCase() === tableGene.toLowerCase()) {
                            row.classList.add("has-background-success-light");
                        }
                    };
                }
            }
        });

    }

    async loadPlotHtml() {

        const prePlotOptionsElt = document.getElementById("plot-options-collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post-plot-adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/multi_gene_as_data.html");
        postPlotOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/multi_gene_as_data.html");

        // populate advanced options for specific plot types
        const prePlotSpecificOptionsElt = document.getElementById("plot-specific-options");
        const postPlotSpecificOptionselt = document.getElementById("post-plot-specific-options");

        // For quadrants and volcanos we load the "series" options in the plot-specific HTML, so that should come first
        if (this.plotType === "quadrant") {
            prePlotSpecificOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/advanced_quadrant.html");
            postPlotSpecificOptionselt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/advanced_quadrant.html");
            return;
        }
        if (this.plotType === "volcano") {
            prePlotSpecificOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/advanced_volcano.html");
            postPlotSpecificOptionselt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/advanced_volcano.html");
            return;
        }
    }

    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = curatorCommon.getPlotConfigValueFromClassName(classElt);
        }

        // convert numerical inputs from plotConfig into Number type
        for (const prop of ["fold_change_cutoff", "fdr_cutoff", "pvalue_threshold", "lower_logfc_threshold", "upper_logfc_threshold"]) {
            if (this.plotConfig[prop]) {
                this.plotConfig[prop] = Number(this.plotConfig[prop]);
            }
        }

        // Get compare series and groups and combine
        const combineSeries = document.querySelector(".js-dash-compare").value;
        curatorCommon.getFacetWidget().filters[combineSeries] = [];
        const refGroup = document.querySelector(".js-dash-reference").value;
        this.plotConfig["ref_condition"] = combineSeries + this.compareSeparator + refGroup;
        curatorCommon.getFacetWidget().filters[combineSeries].push(refGroup);

        if (this.plotType === "volcano") {
            const queryGroup = document.querySelector(".js-dash-query").value;
            this.plotConfig["query_condition"] = combineSeries + this.compareSeparator + queryGroup;
            curatorCommon.getFacetWidget().filters[combineSeries].push(queryGroup);
        }
        if (this.plotType === "quadrant") {
            const compare1Group = document.querySelector(".js-dash-compare1").value;
            const compare2Group = document.querySelector(".js-dash-compare2").value;
            this.plotConfig["compare1_condition"] = combineSeries + this.compareSeparator + compare1Group;
            this.plotConfig["compare2_condition"] = combineSeries + this.compareSeparator + compare2Group;
            curatorCommon.getFacetWidget().filters[combineSeries].push(compare1Group);
            curatorCommon.getFacetWidget().filters[combineSeries].push(compare2Group);
        }

        // Filtered observation groups
        this.plotConfig["obs_filters"] = curatorCommon.getFacetWidget()?.filters || {};
        if (Object.keys(this.plotConfig["obs_filters"]).length === 0) {
            this.plotConfig["obs_filters"] = null; // don't send empty filters
        }

    }

    async setupParamValueCopyEvent() {
        // These plot parameters do not directly correlate to a plot config property
        for (const classElt of ["js-dash-compare", "js-dash-reference", "js-dash-query", "js-dash-compare1", "js-dash-compare2"]) {
            curatorCommon.setupParamValueCopyEvent(classElt);
        }
    }

    async setupPlotSpecificEvents(datasetId) {

        const catColumns = await getCategoryColumns(datasetId);
        curatorCommon.setCatColumns(catColumns);

        if (!catColumns.length) {
            document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
            createToast("No metadata columns found in dataset. Cannot create a plot. Please choose another analysis or choose another dataset.");
            return;
        }

        updateSeriesOptions("js-dash-compare", catColumns);

        const levels = curatorCommon.getLevels();
        // When compare series changes, update the compare groups
        for (const classElt of document.getElementsByClassName("js-dash-compare")) {
            classElt.addEventListener("change", async (event) => {
                const compareSeries = event.target.value;
                updateGroupOptions("js-dash-reference", levels[compareSeries]);
                if (this.plotType === "quadrant") {
                    updateGroupOptions("js-dash-compare1", levels[compareSeries]);
                    updateGroupOptions("js-dash-compare2", levels[compareSeries]);
                }
                if (this.plotType === "volcano") {
                    updateGroupOptions("js-dash-query", levels[compareSeries]);
                }
            });
        }

        // When compare groups change, prevent the same group from being selected in the other compare groups
        for (const classElt of document.getElementsByClassName("js-compare-groups")) {
            classElt.addEventListener("change", (event) => {
                const compareGroups = [...document.getElementsByClassName("js-compare-groups")].map((elt) => elt.value);
                // Filter out empty values and duplicates
                const uniqueGroups = [...new Set(compareGroups)].filter(x => x);
                // Get all unselected groups
                const series = document.querySelector(".js-dash-compare").value;
                const unselectedGroups = levels[series].filter((group) => !uniqueGroups.includes(group));

                for (const innerClassElt of document.getElementsByClassName("js-compare-groups")) {
                    // enable all unselected groups
                    for (const group of unselectedGroups) {
                        const opt = innerClassElt.querySelector(`option[value="${group}"]`);
                        opt.removeAttribute("disabled");
                    }
                    // disable unique groups in other compare groups
                    for (const group of uniqueGroups) {
                        if (innerClassElt.id !== event.target.id) {
                            const opt = innerClassElt.querySelector(`option[value="${group}"]`);
                            opt.setAttribute("disabled", "disabled");
                        }
                    }
                }
            });
        }
    }
}

/**
 * Represents a ScanpyHandler class that extends curatorCommon.PlotHandler.
 * This class is responsible for creating and manipulating plots for a given dataset using the Scanpy analysis object.
 */
class ScanpyHandler extends curatorCommon.PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
        this.apiPlotType = plotType;
    }

    classElt2Prop = {
        "js-tsne-x-axis": "x_axis"
        , "js-tsne-y-axis": "y_axis"
        , "js-tsne-flip-x": "flip_x"
        , "js-tsne-flip-y": "flip_y"
        , "js-tsne-colorize-legend-by": "colorize_legend_by"
        , "js-tsne-max-columns": "max_columns"
        , "js-tsne-horizontal-legend": "horizontal_legend"
        , "js-tsne-marker-size": "marker_size"
        , "js-tsne-color-palette": "expression_palette"
        , "js-tsne-reverse-palette": "reverse_palette"
        , "js-tsne-center-around-median": "center_around_median"
    };

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    /**
     * Clones the display based on the given configuration.
     * @param {Object} config - The configuration object.
     */
    cloneDisplay(config) {
        for (const prop in config) {
            curatorCommon.setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["order"]) {
            for (const series in config["order"]) {
                const order = config["order"][series];
                const levels = curatorCommon.getLevels();
                // sort "levels" series by order
                levels[series].sort((a, b) => order.indexOf(a) - order.indexOf(b));
                curatorCommon.renderOrderSortableSeries(series);
            }

            document.getElementById("order-section").classList.remove("is-hidden");
        }

        // Handle filters
        if (config["obs_filters"]) {
            curatorCommon.getFacetWidget().filters = config["obs_filters"];
        }

        // Restoring some disabled/checked elements in UI
        const horizontalLegend = document.getElementsByClassName("js-tsne-horizontal-legend");

        const catColumns = curatorCommon.getCatColumns();

        if (config["colorize_legend_by"]) {
            const series = config["colorize_legend_by"];
            for (const targetElt of [...horizontalLegend]) {
                targetElt.disabled = true;
                if (catColumns.includes(series)) {
                    targetElt.disabled = false;
                }

                // Applies to horizontal legend
                curatorCommon.disableCheckboxLabel(targetElt, targetElt.disabled);
            }

            // Handle colors
            if (config["colors"]) {
                curatorCommon.renderColorPicker(series);
                for (const group in config["colors"]) {
                    const color = config["colors"][group];
                    const colorField = document.getElementById(`${CSS.escape(group)}-color`);
                    if (colorField) {
                        colorField.value = color;
                    } else {
                        // Optionally log or handle missing field
                        console.warn(`Color field for group "${group}" not found.`);
                    }
                }
            }
        }

        if (config["expression_palette"]) {
            curatorCommon.setSelectBoxByValue("color-palette-post", config["expression_palette"]);
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
            const data = await fetchMgTsneImage(datasetId, analysisObj, this.apiPlotType, this.plotConfig);
            ({ image } = data);
        } catch (error) {
            return;
        }

        const plotContainer = document.getElementById("plot-container");
        plotContainer.replaceChildren();    // erase plot

        // Add a message saying the image dimensions will change when viewed in the gene expression page
        const dimensionMessage = document.createElement("div");
        dimensionMessage.innerHTML = "<strong>Note:</strong> The plot image may render differently when viewed in the gene expression page based on the dimensions of the display tile.";
        dimensionMessage.classList.add("notification", "is-info", "is-light");
        plotContainer.append(dimensionMessage);

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
        };
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

        prePlotOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/pre_plot/tsne_static.html");
        postPlotOptionsElt.innerHTML = await curatorCommon.includeHtml("../include/plot_config/post_plot/tsne_static.html");

        // Remove some single-gene options from the post-plot adjustments
        const plotBySeries = document.querySelector(".js-tsne-plot-by-series");
        const skipGenePlot = document.querySelector(".js-tsne-skip-gene-plot");
        const twoWayPalette = document.querySelector(".js-tsne-two-way-palette");
        for (const targetElt of [plotBySeries, skipGenePlot, twoWayPalette]) {
            targetElt.closest(".is-justify-content-space-between").remove();
        }

        document.querySelector(".js-tsne-max-columns").disabled = false;

        curatorCommon.loadColorscaleSelect(true, true);
    }

    /**
     * Populates the plot configuration based on various elements and values.
     */
    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = curatorCommon.getPlotConfigValueFromClassName(classElt);
        }

        // Get order
        this.plotConfig["order"] = curatorCommon.getPlotOrderFromSortable();

        // Filtered observation groups
        this.plotConfig["obs_filters"] = curatorCommon.getFacetWidget()?.filters || {};
        if (Object.keys(this.plotConfig["obs_filters"]).length === 0) {
            this.plotConfig["obs_filters"] = null;
        }

        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("colorize-legend-by-post").textContent;
        if (colorSeries && colorElts.length) {
            this.plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("-color", "");
                this.plotConfig["colors"][group] = field.value;
            });
        }

        // If user did not want to have a colorized annotation, ensure it does not get passed to the scanpy code
        if (!(colorSeries)) {
            this.plotConfig["max_columns"] = null;
            this.plotConfig["horizontal_legend"] = false;
        }

        // If override marker size is not checked, ensure it does not get passed to the scanpy code
        if (!(document.getElementById("override-marker-size-post").checked)) {
            this.plotConfig["marker_size"] = null;
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
     * @param {string} datasetId - The ID of the dataset.
     * @returns {Promise<void>} A promise that resolves when the setup is complete.
     */
    async setupPlotSpecificEvents(datasetId) {
        await setupScanpyOptions();
    }

}

const adjustGeneTableLabels = (plotType) => {
    const geneX = document.getElementById("tbl-gene-x");
    const geneY = document.getElementById("tbl-gene-y");

    // Adjust headers to the plot type
    if (plotType === "quadrant") {
        geneX.innerHTML = 'X Log2 FC <span class="icon"><i class="mdi mdi-sort-numeric-ascending" aria-hidden="true"></i></span>';
        geneY.innerHTML = 'Y Log2 FC <span class="icon"><i class="mdi mdi-sort-numeric-ascending" aria-hidden="true"></i></span>';
    } else {
        // volcano
        geneX.innerHTML = 'Log2 FC <span class="icon"><i class="mdi mdi-sort-numeric-ascending" aria-hidden="true"></i></span>';
        geneY.innerHTML = 'P-value <span class="icon"><i class="mdi mdi-sort-numeric-ascending" aria-hidden="true"></i></span>';
    }
};


const appendGeneTagButton = (geneTagElt) => {
    // Add delete button
    const deleteBtnElt = document.createElement("button");
    deleteBtnElt.classList.add("delete", "is-small");
    geneTagElt.appendChild(deleteBtnElt);
    deleteBtnElt.addEventListener("click", (event) => {
        // Remove gene from selectedGenes
        const gene = event.target.parentNode.textContent;
        geneCollectionState.selectedGenes.delete(gene);
        event.target.parentNode.remove();

        // Remove checkmark from gene lists dropdown
        const geneListLabel = document.querySelector(`#dropdown-content-genes .gene-item-label[text="${gene}"]`);
        if (!geneListLabel) {
            return;
        }
        const geneListElt = geneListLabel.parentElement;
        const geneListI = geneListElt.querySelector("i.toggler");
        if (!geneListI.classList.contains("mdi-check")) {
            return;
        }
        geneListI.classList.replace("mdi-check", "mdi-plus");
        geneListI.classList.replace("gene-list-item-remove", "gene-list-item-add");
        geneListElt.classList.remove("is-selected");
    });

};

const clearGenes = (event) => {
    document.getElementById("clear-genes-btn").classList.add("is-loading");
    document.getElementById("gene-tags").replaceChildren();
    geneCollectionState.selectedGenes.clear();
    document.getElementById("dropdown-gene-list-cancel").click();	// clear the dropdown
    document.getElementById("clear-genes-btn").classList.remove("is-loading");
    document.getElementById("genes-manually-entered").value = "";
};

/**
 * Triggered when a gene is selected.
 *
 * @param {Array} genes - The selected genes.
 */
const chooseGenes = (genes=[]) => {

    manuallyEnteredGenes = new Set(genes);
    geneCollectionState.selectedGenes = manuallyEnteredGenes;
    const geneSymbolString = genes.join(" ");
    document.getElementById('genes-manually-entered').value = geneSymbolString;

    // Delete existing tags
    const geneTagsElt = document.getElementById("gene-tags");
    geneTagsElt.replaceChildren();

    const genesNotFoundTagsElt = document.getElementById("genes-not-found-tags");
    genesNotFoundTagsElt.replaceChildren();

    // Remove selected genes that are not in datasetGenes (case-insensitive)
    const datasetGenesLower = new Set(Array.from(datasetGenes, g => g.toLowerCase()));

    for (const gene of geneCollectionState.selectedGenes) {
        if (gene.trim() === "") continue;  // Skip empty genes

        if (!datasetGenesLower.has(gene.toLowerCase())) {
            geneCollectionState.selectedGenes.delete(gene);
            notFoundGenes.add(gene);
        }
    }

    // Create tags for each gene not found.
    document.getElementById("genes-not-found-c").classList.remove("is-hidden");
    for (const gene of notFoundGenes) {
        const geneTagElt = document.createElement("span");
        geneTagElt.classList.add("tag", "is-primary", "is-light", "mx-1");
        geneTagElt.textContent = gene;
        genesNotFoundTagsElt.appendChild(geneTagElt);
    }

    document.getElementById("num-selected-genes-c").classList.remove("is-hidden");
    document.getElementById("num-selected-genes").textContent = geneCollectionState.selectedGenes.size;
    document.getElementById("num-selected-genes-post").textContent = geneCollectionState.selectedGenes.size;

    if (geneCollectionState.selectedGenes.size == 0) return;  // Do not trigger after initial population

    // Update list of gene tags
    const sortedGenes = Array.from(geneCollectionState.selectedGenes).sort();
    for (const opt in sortedGenes) {
        const geneTagElt = document.createElement("span");
        geneTagElt.classList.add("tag", "is-primary", "mx-1");
        geneTagElt.textContent = sortedGenes[opt];
        appendGeneTagButton(geneTagElt);
        geneTagsElt.appendChild(geneTagElt);
    }

    document.getElementById("gene-tags-c").classList.remove("is-hidden");

    if (geneCollectionState.selectedGenes.size < 2) {
        document.getElementById("gene-s-failed").classList.remove("is-hidden");
        document.getElementById("gene-s-success").classList.add("is-hidden");
        for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
            plotBtn.disabled = true;
        }
        document.getElementById("continue-to-plot-options").classList.add("is-hidden");
        return;
    }

    // If more than 10 tags, hide the rest and add a "show more" button
    if (geneCollectionState.selectedGenes.size > 10) {
        const geneTags = geneTagsElt.querySelectorAll("span.tag");
        for (let i = 10; i < geneTags.length; i++) {
            geneTags[i].classList.add("is-hidden");
        }
        // Add show more button
        const showMoreBtnElt = document.createElement("button");
        showMoreBtnElt.classList.add("tag", "button", "is-small", "is-primary", "is-light");
        const numToDisplay = geneCollectionState.selectedGenes.size - 10;
        showMoreBtnElt.textContent = `+${numToDisplay} more`;
        showMoreBtnElt.addEventListener("click", (event) => {
            const geneTags = geneTagsElt.querySelectorAll("span.tag");
            for (let i = 10; i < geneTags.length; i++) {
                geneTags[i].classList.remove("is-hidden");
            }
            event.target.remove();
        });
        geneTagsElt.appendChild(showMoreBtnElt);
    }

    document.getElementById("gene-s-failed").classList.add("is-hidden");
    document.getElementById("gene-s-success").classList.remove("is-hidden");

    document.getElementById("continue-to-plot-options").classList.remove("is-hidden");

};
curatorCommon.registerChooseGenes(chooseGenes);


const curatorSpecifcCreatePlot = async (plotType, datasetId, analysisObj) => {

    const plotStyle = curatorCommon.getPlotStyle();
    // Add gene or genes to plot config
    plotStyle.plotConfig["gene_symbols"] = Array.from(geneCollectionState.selectedGenes);

    // Call API route by plot type
    await plotStyle.createPlot(datasetId, analysisObj);
};
curatorCommon.registerCuratorSpecifcCreatePlot(curatorSpecifcCreatePlot);

const curatorSpecifcDatasetTreeCallback = async () => {

    document.getElementById("num-selected-genes").textContent = 0;
    document.getElementById("num-selected-genes-post").textContent = 0;
};
curatorCommon.registerCuratorSpecifcDatasetTreeCallback(curatorSpecifcDatasetTreeCallback);


const curatorSpecificNavbarUpdates = () => {
    // Update with current page info
    document.getElementById("page-header-label").textContent = "Multi-gene Displays";
};
curatorCommon.registerCuratorSpecificNavbarUpdates(curatorSpecificNavbarUpdates);


const curatorSpecificOnLoad = async () => {
    await fetchGeneCartData();
    // Should help with lining things up on index page
    document.getElementById("dropdown-gene-lists").classList.remove("is-right");
};
curatorCommon.registerCuratorSpecificOnLoad(curatorSpecificOnLoad);

const curatorSpecificPlotStyle = (plotType) => {
    // include plotting backend options
    if (genesAsAxisPlots.includes(plotType)) {
        return new GenesAsAxisHandler(plotType);
    } else if (genesAsDataPlots.includes(plotType)) {
        return new GenesAsDataHandler(plotType);
    } else if (scanpyPlots.includes(plotType)) {
        return new ScanpyHandler(plotType);
    } else {
        return null;
    }
};
curatorCommon.registerCuratorSpecificPlotStyle(curatorSpecificPlotStyle);

const curatorSpecificPlotTypeAdjustments = (plotType) => {
    return plotType;
};
curatorCommon.registerCuratorSpecificPlotTypeAdjustments(curatorSpecificPlotTypeAdjustments);

const curatorSpecificUpdateDatasetGenes = async (geneSymbols) => {
    // Convert geneSymbols to a set
    // This will be used to check manually entered genes and those from gene lists
    datasetGenes = new Set(geneSymbols);
};
registerCuratorSpecificUpdateDatasetGenes(curatorSpecificUpdateDatasetGenes);

const curatorSpecificValidationChecks = () => {
    return true;
};
curatorCommon.registerCuratorSpecificValidationChecks(curatorSpecificValidationChecks);

const downloadSelectedGenes = (event) => {
    event.preventDefault();

    // Builds a file in memory for the user to download.  Completely client-side.
    // plot_data contains three keys: x, y and symbols
    // build the file string from this

    const plotStyle = curatorCommon.getPlotStyle();

    const plotType = plotStyle.plotType;
    const plotConfig = plotStyle.plotConfig;

    // Adjust headers to the plot type
    let xLabel;
    let yLabel;

    if (plotType === "quadrant") {
        const query1 = plotConfig.compare1_condition.split(';-;')[1];
        const query2 = plotConfig.compare2_condition.split(';-;')[1];
        const ref = plotConfig.ref_condition.split(';-;')[1];
        xLabel = `${query1} vs ${ref} Log2FC`;
        yLabel = `${query2} vs ${ref} Log2FC`;
    } else {
        const query = plotConfig.query_condition.split(';-;')[1];
        let ref = plotConfig.ref_condition.split(';-;')[1];

        ref = ref === "Union of the rest of the groups" ? "rest" : ref;

        // volcano
        xLabel = `${query} vs ${ref} Log2FC`;
        yLabel = `${query} vs ${ref} p-value`;
    }
    let fileContents = `gene_symbol\t${xLabel}\t${yLabel}\n`;

    // Entering genes and info now.
    plotSelectedGenes.forEach((gene) => {
        fileContents +=
            `${gene.gene_symbol}\t`
            + `${gene.x}\t`
            + `${gene.y}\n`;
    });

    const element = document.createElement("a");
    element.setAttribute(
        "href",
        `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(fileContents)}`
    );
    element.setAttribute("download", "geneCollectionState.selectedGenes.tsv");
    element.style.display = "none";
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
};

const fetchMgPlotlyData = async (datasetId, analysis, plotType, plotConfig) => {
    try {
        const data = await apiCallsMixin.fetchMgPlotlyData(datasetId, analysis, plotType, plotConfig);
        if (data?.success < 1) {
            throw new Error(data?.message || "Unknown error.");
        }
        return data;
    } catch (error) {
        const data = error?.response?.data;
        if (data?.success < 1) {
            const msg = "Exceeded memory limit. Please try to filter or subsample the dataset in the post-plotting view, or contact the gEAR team.";
            createToast(msg);
            throw error;
        }

        logErrorInConsole(error);
        const msg = "Could not create plot for this dataset and parameters. Please contact the gEAR team.";
        createToast(msg);
        throw new Error(msg);
    }
};

/**
 * Fetches the multigene TSNE image for a given dataset, analysis, plot type, and plot configuration.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} analysis - The analysis type.
 * @param {string} plotType - The type of plot.
 * @param {object} plotConfig - The configuration for the plot.
 * @returns {Promise<object>} - The fetched data.
 * @throws {Error} - If there is an error fetching the data or creating the plot image.
 */
const fetchMgTsneImage = async (datasetId, analysis, plotType, plotConfig) => {
    // NOTE: gene_symbol already passed to plotConfig
    try {
        const data = await apiCallsMixin.fetchMgTsneImage(datasetId, analysis, plotType, plotConfig);
        if (data?.success < 1) {
            throw new Error(data?.message ? data.message : "Unknown error.");
        }
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not create plot image for this dataset and parameters. Please contact the gEAR team.";
        createToast(msg);
        throw new Error(msg);
    }
};

const getCategoryColumns = async (datasetId) => {
    const analysisId = curatorCommon.getAnalysisId();
    let levels;
    try {
        ({ obs_columns: allColumns, obs_levels: levels } = await curatorCommon.curatorApiCallsMixin.fetchH5adInfo(datasetId, analysisId));
    } catch (error) {
        document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
        return;
    }
    // Filter out values we don't want of "levels", like "colors"
    for (const key in levels) {
        if (key.includes("_colors")) {
            delete levels[key];
        }
    }
    curatorCommon.setLevels(levels);
    return Object.keys(levels);
};

// Invert a log function
const invertLogFunction = (value, base = 10) => {
    return base ** value;
};


const populateGeneTable = (data, plotType) => {
    plotSelectedGenes = [];

    for (const pt of data.points) {
        plotSelectedGenes.push({
            gene_symbol: pt.data.text[pt.pointNumber],
            ensembl_id: pt.data.customdata[pt.pointNumber], // Ensembl ID stored in "customdata" property
            x: pt.data.x[pt.pointNumber].toFixed(1),
            y: plotType === "volcano" ? invertLogFunction(-pt.data.y[pt.pointNumber]).toExponential(2) : pt.data.y[pt.pointNumber].toFixed(2),
        });
    };

    // Sort in alphabetical order
    plotSelectedGenes.sort();


    const geneTableBody = document.getElementById("gene-table-body");
    geneTableBody.replaceChildren();

    for (const gene of plotSelectedGenes) {
        const row = document.createElement("tr");
        row.innerHTML = `<td>${gene.gene_symbol}</td><td>${gene.x}</td><td>${gene.y}</td>`;
        geneTableBody.appendChild(row);
    }
};

const saveGeneCart = () => {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: sessionId
        , label: document.getElementById("new-genecart-label").value
        , gctype: "unweighted-list"
        , organism_id: getOrganismId()
        , is_public: 0
    });

    for (const sg of plotSelectedGenes) {
        const gene = new Gene({
            id: sg.ensembl_id, // Ensembl ID stored in "customdata" property
            gene_symbol: sg.gene_symbol,
        });
        gc.addGene(gene);
    }

    gc.save(() => { }, updateUIAfterGeneCartSaveFailure);
};

const saveWeightedGeneCart = () => {

    const plotStyle = curatorCommon.getPlotStyle();

    // must have access to USER_SESSION_ID
    const plotType = plotStyle.plotType;
    const plotConfig = plotStyle.plotConfig;

    // Saving raw FC by default so it is easy to transform weight as needed
    const foldchangeLabel = "FC";
    let weightLabels = [foldchangeLabel];


    if (plotType === "quadrant") {
        const query1 = plotConfig.compare1_condition.split(';-;')[1];
        const query2 = plotConfig.compare2_condition.split(';-;')[1];
        const ref = plotConfig.ref_condition.split(';-;')[1];
        const xLabel = `${query1}-vs-${ref}`;
        const yLabel = `${query2}-vs-${ref}`;

        const fcl1 = `${xLabel}-${foldchangeLabel}`;
        const fcl2 = `${yLabel}-${foldchangeLabel}`;

        weightLabels = [fcl1, fcl2];
    }

    const gc = new WeightedGeneCart({
        session_id: sessionId
        , label: document.getElementById("new-genecart-label").value
        , gctype: 'weighted-list'
        , organism_id: getOrganismId()
        , is_public: 0
    }, weightLabels);

    // Volcano and Quadrant plots have multiple traces of genes, broken into groups.
    // Loop through these to get the info we need.
    for (const trace of plotStyle.plotJson.data) {
        for (const pt in trace.x) {

            const foldchange = Number(trace.x[pt].toFixed(1));
            const weights = [foldchange];

            // If quadrant plot was specified, there are fold changes in x and y axis
            if (plotType === "quadrant") {
                const foldchange2 = Number(trace.y[pt].toFixed(1));
                weights.push(foldchange2);
            }

            const gene = new WeightedGene({
                id: trace.customdata[pt], // Ensembl ID stored in "customdata" property
                gene_symbol: trace.text[pt]
            }, weights);
            gc.addGene(gene);
        }
    };

    gc.save(() => { }, updateUIAfterGeneCartSaveFailure);
};

/**
 * Sets up the options for Scanpy analysis.
 * @param {string} datasetId - The ID of the dataset.
 * @returns {Promise<void>} A promise that resolves when the setup is complete.
 */
const setupScanpyOptions = async (datasetId) => {
    const analysisId = curatorCommon.getAnalysisId();
    const plotType = curatorCommon.getSelect2Value(curatorCommon.getPlotTypeSelect());
    let levels;
    try {
        ({ obs_columns: allColumns, obs_levels: levels } = await curatorCommon.curatorApiCallsMixin.fetchH5adInfo(datasetId, analysisId));
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

    curatorCommon.setLevels(levels);

    if (!allColumns.length) {
        document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
        createToast("No metadata columns found in dataset. Cannot create a plot. Please choose another analysis or choose another dataset.");
        return;
    }

    curatorCommon.setCatColumns(Object.keys(levels));
    const catColumns = curatorCommon.getCatColumns();

    let xDefaultOption = null;
    let yDefaultOption = null;

    // If these exist, make the default option
    switch (plotType) {
        case "mg_pca_static":
            xDefaultOption = "X_pca_1";
            yDefaultOption = "X_pca_2";
            break;
        case "mg_tsne_static":
            xDefaultOption = "X_tsne_1";
            yDefaultOption = "X_tsne_2";
            break;
        case "mg_umap_static":
            xDefaultOption = "X_umap_1";
            yDefaultOption = "X_umap_2";
            break;
    }

    updateSeriesOptions("js-tsne-x-axis", allColumns, true, xDefaultOption);
    updateSeriesOptions("js-tsne-y-axis", allColumns, true, yDefaultOption);
    updateSeriesOptions("js-tsne-colorize-legend-by", allColumns, false);

    const colorizeLegendBy = document.getElementsByClassName("js-tsne-colorize-legend-by");
    const horizontalLegend = document.getElementsByClassName("js-tsne-horizontal-legend");

    // Do certain things if the chosen annotation series is categorical or continuous
    for (const elt of colorizeLegendBy) {
        elt.addEventListener("change", (event) => {
            for (const targetElt of [...horizontalLegend]) {
                targetElt.disabled = true;
                // If colorized legend is continuous, we cannot plot by group
                // So all dependencies need to be disabled.
                if ((catColumns.includes(event.target.value))) {
                    targetElt.disabled = false;
                    curatorCommon.disableCheckboxLabel(targetElt, false);
                }
            }
        });

        //
        elt.addEventListener("change", (event) => {
            // if series is empty or not categorical, remove color picker
            if ((catColumns.includes(event.target.value))) {
                curatorCommon.renderColorPicker(event.target.value);
                return;
            }
            const colorsContainer = document.getElementById("colors-container");
            const colorsSection = document.getElementById("colors-section");
            colorsSection.classList.add("is-hidden");
            colorsContainer.replaceChildren();
            return;
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

};

/**
 * Shows the corresponding subsection based on the selected option in the plot configuration menu.
 * @param {Event} event - The event triggered by the user's selection.
 */
const showPostHeatmapParamSubsection = (event) => {
    for (const subsection of document.getElementsByClassName("js-plot-config-section")) {
        subsection.classList.add("is-hidden");
    }

    switch (event.target.textContent.trim()) {
        case "Clustering":
            document.getElementById("clustering-section-post").classList.remove("is-hidden");
            break;
        default:
            document.getElementById("misc-section-post").classList.remove("is-hidden");
            break;
    }
    event.preventDefault(); // Prevent "link" clicking from "a" elements
};

// Taken from https://www.w3schools.com/howto/howto_js_sort_table.asp
const sortGeneTable = (mode) => {
    let table;
    let rows;
    let switching;
    let i;
    let x;
    let y;
    let shouldSwitch;
    let dir;
    let switchcount = 0;
    table = document.getElementById("tbl-selected-genes");

    switching = true;
    // Set the sorting direction to ascending:
    dir = "asc";
    /* Make a loop that will continue until
        no switching has been done: */
    while (switching) {
        // Start by saying: no switching is done:
        switching = false;
        rows = table.rows;
        /* Loop through all table rows (except the
        first, which contains table headers): */
        for (i = 1; i < rows.length - 1; i++) {
            // Start by saying there should be no switching:
            shouldSwitch = false;
            /* Get the two elements you want to compare,
                one from current row and one from the next: */
            x = rows[i].getElementsByTagName("td")[mode];
            y = rows[i + 1].getElementsByTagName("td")[mode];
            /* Check if the two rows should switch place,
                based on the direction, asc or desc: */
            if (dir == "asc") {
                // First column is gene_symbol... rest are numbers
                if (mode === 0 && x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                    // If so, mark as a switch and break the loop:
                    shouldSwitch = true;
                    break;
                }
                if (Number(x.innerHTML) > Number(y.innerHTML)) {
                    shouldSwitch = true;
                    break;
                }
            } else if (dir == "desc") {
                if (mode === 0 && x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
                    // If so, mark as a switch and break the loop:
                    shouldSwitch = true;
                    break;
                }
                if (Number(x.innerHTML) < Number(y.innerHTML)) {
                    shouldSwitch = true;
                    break;
                }
            }
        }
        if (shouldSwitch) {
            /* If a switch has been marked, make the switch
                and mark that a switch has been done: */
            rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
            switching = true;
            // Each time a switch is done, increase this count by 1:
            switchcount++;

        } else {
            /* If no switching has been done AND the direction is "asc",
                set the direction to "desc" and run the while loop again. */
            if (switchcount == 0 && dir == "asc") {
                dir = "desc";
                switching = true;
            }
        }
    }

    // Reset other sort icons to "ascending" state, to show what direction they will sort when clicked
    const otherTblHeaders = document.querySelectorAll(`.js-tbl-gene-header:not(:nth-child(${mode + 1}))`);
    for (const tblHeader of otherTblHeaders) {
        const currIcon = tblHeader.querySelector("i");
        if (mode == 0) {
            currIcon.classList.remove("mdi-sort-alphabetical-descending");
            currIcon.classList.add("mdi-sort-alphabetical-ascending");
        } else {
            currIcon.classList.remove("mdi-sort-numeric-descending");
            currIcon.classList.add("mdi-sort-numeric-ascending");
        }
    }

    // toggle the mdi icons between ascending / descending
    // icon needs to reflect the current state of the sort
    const selectedTblHeader = document.querySelector(`.js-tbl-gene-header:nth-child(${mode + 1})`);
    const currIcon = selectedTblHeader.querySelector("i");
    if (dir == "asc") {
        if (mode == 0) {
            currIcon.classList.remove("mdi-sort-alphabetical-descending");
            currIcon.classList.add("mdi-sort-alphabetical-ascending");
        } else {
            currIcon.classList.remove("mdi-sort-numeric-descending");
            currIcon.classList.add("mdi-sort-numeric-ascending");
        }
    } else {
        if (mode == 0) {
            currIcon.classList.remove("mdi-sort-alphabetical-ascending");
            currIcon.classList.add("mdi-sort-alphabetical-descending");
        } else {
            currIcon.classList.remove("mdi-sort-numeric-ascending");
            currIcon.classList.add("mdi-sort-numeric-descending");
        }
    }
};


// For plotting options, populate select menus with category groups
const updateSeriesOptions = (classSelector, seriesArray) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add categories
        for (const group of seriesArray.sort()) {

            const option = document.createElement("option");
            option.textContent = group;
            option.value = group;
            elt.append(option);
        }
    }
};


// For a given categorical series (e.g. "celltype"), update the "group" options
const updateGroupOptions = (classSelector, groupsArray) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add categories
        for (const group of groupsArray.sort()) {
            const option = document.createElement("option");
            option.textContent = group;
            option.value = group;
            elt.append(option);
        }
    }

};

const updateUIAfterGeneCartSaveFailure = (gc, message) => {
    createToast(message);
};

document.getElementById("clear-genes-btn").addEventListener("click", clearGenes);

// code from Bulma documentation to handle modals
document.getElementById("gene-cart-btn").addEventListener("click", ($trigger) => {
    const closestButton = $trigger.target.closest(".button");
    const modal = closestButton.dataset.target;
    const $target = document.getElementById(modal);
    openModal($target);

});

document.getElementById("new-genecart-label").addEventListener("input", (event) => {
    const saveBtn = document.getElementById("save-genecart-btn");
    saveBtn.disabled = event.target.value ? false : true;
});

document.getElementById("save-genecart-btn").addEventListener("click", (event) => {
    event.preventDefault();
    event.target.classList.add("is-loading");
    // get value of genecart radio button group
    const geneCartName = document.querySelector("input[name='genecart_type']:checked").value;
    if (getCurrentUser()) {
        if (geneCartName === "unweighted") {
            saveGeneCart();
        } else {
            saveWeightedGeneCart();
        }
    }
    event.target.classList.remove("is-loading");
});

document.getElementById("download-selected-genes-btn").addEventListener("click", downloadSelectedGenes);

// handle when the dropdown-gene-list-search-input input box is changed
document.getElementById('genes-manually-entered').addEventListener('change', (event) => {
    const searchTermString = event.target.value;
    const newManuallyEnteredGenes = searchTermString.length > 0 ? new Set(searchTermString.split(/[ ,]+/)) : new Set();

    // Remove genes that have been deleted from the geneCollectionState.selectedGenes set
    for (const gene of manuallyEnteredGenes) {
        if (!newManuallyEnteredGenes.has(gene)) {
            geneCollectionState.selectedGenes.delete(gene);
            notFoundGenes.delete(gene);
        }
    }

    // Add new genes to the selectedGenes set
    for (const gene of newManuallyEnteredGenes) {
        geneCollectionState.selectedGenes.add(gene);
    }

    manuallyEnteredGenes = newManuallyEnteredGenes;
    chooseGenes([]);
});

document.getElementById('dropdown-gene-list-proceed').addEventListener('click', (event) => {chooseGenes([])});
