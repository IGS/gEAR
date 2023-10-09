const isMultigene = 1;

let geneSelect = null;

const genesAsAxisPlots = ["dotplot", "heatmap", "mg_violin"];
const genesAsDataPlots = ["quadrant", "volcano"];

class GenesAsAxisHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
    }

    classElt2Prop = {
        "js-dash-primary":"primary_col"
        , "js-dash-secondary":"secondary_col"
        , "js-dash-color-palette":"colorscale"
        , "js-dash-reverse-palette":"reverse_colorscale"
        , "js-dash-distance-metric":"distance_metric"
        , "js-dash-matrixplot":"matrixplot"
        , "js-dash-center-mean": "center_around_zero"
        , "js-dash-cluster-obs": "cluster_obs"
        , "js-dash-cluster-genes": "cluster_genes"
        , "js-dash-flip-axes": "flip_axes"
        , "js-dash-hide-obs-labels": "hide_obs_labels"
        , "js-dash-hide-gene-labels": "hide_gene_labels"
        , "js-dash-add-jitter": "violin_add_points"
        , "js-dash-stacked-violin": "stacked_violin"
        , "js-dash-plot-title": "plot_title"
        , "js-dash-legend-title": "legend_title"
    }

    // deal with clusterbar separately since it is an array of selected options

    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    cloneDisplay(config) {
        // load plot values
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle order
        if (config["sort_order"]) {
            for (const series in config["sort_order"]) {
                const order = config["sort_order"][series];
                // sort "levels" series by order
                levels[series].sort((a, b) => order.indexOf(a) - order.indexOf(b));
                renderOrderSortableSeries(series);
            }

            document.getElementById("order_section").classList.remove("is-hidden");
        }

        // Handle filters
        if (config["obs_filters"]) {
            facetWidget.filters = config["obs_filters"];
        }

        // handle colors
        if (config["color_palette"]) {
            setSelectBoxByValue("color_palette_post", config["color_palette"]);
            //colorscaleSelect.update();
        }

        // handle clusterbar values
        if (config["clusterbar_fields"]) {
            for (const field of config["clusterbar_fields"]) {
                const elt = document.querySelector(`#clusterbar_c .js-dash-clusterbar-checkbox[value="${field}"]`);
                elt.checked = true;
            }
        }

    }

    async createPlot(datasetId, analysisObj, userId, colorblindMode) {
        // Get data and set up the image area
        let plotJson;
        try {
            const data = await fetchDashData(this.plotConfig, datasetId, this.plotType, analysisObj, userId, colorblindMode);
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

        if (this.plotType === 'heatmap') {
            setHeatmapHeightBasedOnGenes(plotJson.layout, this.plotConfig.gene_symbols);
        } else if (this.plotType === "mg_violin" && this.plotConfig.stacked_violin){
            adjustStackedViolinHeight(plotJson.layout);
        }

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "config");
        Plotly.newPlot("plotly_preview", plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(curatorDisplayConf, this.plotType, "layout")
        Plotly.relayout("plotly_preview", custonLayout)

        document.getElementById("legend_title_container").classList.remove("is-hidden");
        if (this.plotType === "dotplot") {
            document.getElementById("legend_title_container").classList.add("is-hidden");
        }

    }

    async loadPlotHtml() {
        const prePlotOptionsElt = document.getElementById("plot_options_collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post_plot_adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/multi_gene_as_axis.html");
        postPlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/post_plot/multi_gene_as_axis.html");

        // populate advanced options for specific plot types
        const prePlotSpecificOptionsElt = document.getElementById("plot_specific_options");
        const postPlotSpecificOptionselt = document.getElementById("post_plot_specific_options");

        // Load color palette select options
        const isContinuous = ["dotplot", "heatmap"].includes(this.plotType) ? true : false;
        loadColorscaleSelect(isContinuous);


        if (this.plotType === "heatmap") {
            prePlotSpecificOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/advanced_heatmap.html");
            postPlotSpecificOptionselt.innerHTML = await includeHtml("../include/plot_config/post_plot/advanced_heatmap.html");
            return;
        }
        if (this.plotType === "mg_violin") {
            prePlotSpecificOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/advanced_mg_violin.html");
            postPlotSpecificOptionselt.innerHTML = await includeHtml("../include/plot_config/post_plot/advanced_mg_violin.html");
            return;
        }
    }

    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = getPlotConfigValueFromClassName(classElt)
        }

        // Get checked clusterbar values
        if (this.plotType === "heatmap") {
            const clusterbarValues = [];
            // They should be synced so just grab the first set of clusterbar values
            for (const elt of document.querySelectorAll("#clusterbar_c .js-dash-clusterbar-checkbox")) {
                if (elt.checked) {
                    clusterbarValues.push(elt.value);
                }
            }
            this.plotConfig["clusterbar_fields"] = clusterbarValues;
        }

        // Filtered observation groups
        this.plotConfig["obs_filters"] = facetWidget?.filters || {};

        // Get order
        this.plotConfig["sort_order"] = getPlotOrderFromSortable();

    }

    async setupParamValueCopyEvent() {
        //pass
    }

    async setupPlotSpecificEvents() {
        catColumns = await getCategoryColumns();

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
                    opt.setAttribute("disabled", "disabled");
                    // If this option was selected, unselect it
                    if (secondaryClassElt.value === primarySeries) {
                        secondaryClassElt.value = "";
                    }
                }
            })
        }

        // For clusterbar options, create checkboxes for all catColumns
        for (const classElt of document.getElementsByClassName("js-dash-clusterbar")) {
            for (const catColumn of catColumns) {

                const template = `<div class="control">
                    <label class="checkbox">
                        <input type="checkbox" value="${catColumn}" class="js-dash-clusterbar-checkbox" />
                        ${catColumn}
                    </label>
                <div>`

                const html = generateElements(template);
                classElt.append(html);
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
            })
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
    }

}

class GenesAsDataHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
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

    cloneDisplay(config) {
        // load plot values
        for (const prop in config) {
            setPlotEltValueFromConfig(this.configProp2ClassElt[prop], config[prop]);
        }

        // Handle filters
        if (config.hasOwnProperty("obs_filters")) {
            facetWidget.filters = config["obs_filters"];
        }

        // Split compare series and groups
        const refCondition = this.plotConfig["ref_condition"];
        const [combineSeries, refGroup] = refCondition.split(this.compareSeparator);
        for (const classElt of document.getElementsByClassName("js-dash-compare")) {
            classElt.value = combineSeries;
        }
        for (const classElt of document.getElementsByClassName("js-dash-reference")) {
            classElt.value = refGroup;
        }

        if (this.plotType === "volcano") {
            const queryCondition = this.plotConfig["query_condition"];
            const queryGroup = queryCondition.split(this.compareSeparator)[1];
            for (const classElt of document.getElementsByClassName("js-dash-query")) {
                classElt.value = queryGroup;
            }
        }
        if (this.plotType === "quadrant") {
            const compare1Condition = this.plotConfig["compare1_condition"];
            const compare1Group = compare1Condition.split(this.compareSeparator)[1];
            const compare2Condition = this.plotConfig["compare2_condition"];
            const compare2Group = compare2Condition.split(this.compareSeparator)[1];
            for (const classElt of document.getElementsByClassName("js-dash-compare1")) {
                classElt.value = compare1Group;
            }
            for (const classElt of document.getElementsByClassName("js-dash-compare2")) {
                classElt.value = compare2Group;
            }
        }
    }

    async createPlot(datasetId, analysisObj, userId, colorblindMode) {
        // Get data and set up the image area
        let plotJson;
        try {
            const data = await fetchDashData(this.plotConfig, datasetId, this.plotType, analysisObj, userId, colorblindMode);
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

    }

    async loadPlotHtml() {

        const prePlotOptionsElt = document.getElementById("plot_options_collapsable");
        prePlotOptionsElt.replaceChildren();

        const postPlotOptionsElt = document.getElementById("post_plot_adjustments");
        postPlotOptionsElt.replaceChildren();

        prePlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/multi_gene_as_data.html");
        postPlotOptionsElt.innerHTML = await includeHtml("../include/plot_config/post_plot/multi_gene_as_data.html");

        // populate advanced options for specific plot types
        const prePlotSpecificOptionsElt = document.getElementById("plot_specific_options");
        const postPlotSpecificOptionselt = document.getElementById("post_plot_specific_options");

        // For quadrants and volcanos we load the "series" options in the plot-specific HTML, so that should come first
        if (this.plotType === "quadrant") {
            prePlotSpecificOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/advanced_quadrant.html");
            postPlotSpecificOptionselt.innerHTML = await includeHtml("../include/plot_config/post_plot/advanced_quadrant.html");
            return;
        }
        if (this.plotType === "volcano") {
            prePlotSpecificOptionsElt.innerHTML = await includeHtml("../include/plot_config/pre_plot/advanced_volcano.html");
            postPlotSpecificOptionselt.innerHTML = await includeHtml("../include/plot_config/post_plot/advanced_volcano.html");
            return;
        }
    }

    populatePlotConfig() {
        this.plotConfig = {};   // Reset plot config

        for (const classElt in this.classElt2Prop) {
            this.plotConfig[this.classElt2Prop[classElt]] = getPlotConfigValueFromClassName(classElt)
        }

        // Get compare series and groups and combine
        const combineSeries = document.querySelector(".js-dash-compare").value;
        this.plotConfig["ref_condition"] = combineSeries + this.compareSeparator + document.querySelector(".js-dash-reference").value;
        if (this.plotType === "volcano") {
            const queryGroup = document.querySelector(".js-dash-query").value;
            this.plotConfig["query_condition"] = combineSeries + this.compareSeparator + queryGroup;
        }
        if (this.plotType === "quadrant") {
            const compare1Group = document.querySelector(".js-dash-compare1").value;
            const compare2Group = document.querySelector(".js-dash-compare2").value;
            this.plotConfig["compare1_condition"] = combineSeries + this.compareSeparator + compare1Group;
            this.plotConfig["compare2_condition"] = combineSeries + this.compareSeparator + compare2Group;
        }

        // Filtered observation groups
        this.plotConfig["obs_filters"] = facetWidget?.filters || {};

    }

    async setupParamValueCopyEvent() {
        // These plot parameters do not directly correlate to a plot config property
        for (const classElt of ["js-dash-compare", "js-dash-reference", "js-dash-query", "js-dash-compare1", "js-dash-compare2"]) {
            setupParamValueCopyEvent(classElt)
        }
    }

    async setupPlotSpecificEvents() {

        catColumns = await getCategoryColumns();
        updateSeriesOptions("js-dash-compare", catColumns);

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
            })
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
            })
        }
    }
}

const geneCartTree = new GeneCartTree({
    element: document.getElementById("genecart_tree")
    , searchElement: document.getElementById("genecart_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "genecart") {
            return;
        }

        // Get gene symbols from gene cart
        const geneCartId = e.node.data.orig_id;
        const geneCartMembers = await fetchGeneCartMwembers(sessionId, geneCartId);
        const geneCartSymbols = geneCartMembers.map((item) => item.label);

        // Normalize gene symbols to lowercase
        const geneSelectSymbols = geneSelect.data.map((opt) => opt.value);
        const geneCartSymbolsLowerCase = geneCartSymbols.map((x) => x.toLowerCase());

        const geneSelectedOptions = geneSelect.selectedOptions.map((opt) => opt.data.value);

        // Get genes from gene cart that are present in dataset's genes.  Preserve casing of dataset's genes.
        const geneCartIntersection = geneSelectSymbols.filter((x) => geneCartSymbolsLowerCase.includes(x.toLowerCase()));
        // Add in already selected genes (union)
        const geneSelectIntersection = [...new Set(geneCartIntersection.concat(geneSelectedOptions))];

        // change all options to be unselected
        const origSelect = document.getElementById("gene_select");
        for (const opt of origSelect.options) {
            opt.removeAttribute("selected");
        }

        // Assign intersection genes to geneSelect "selected" options
        for (const gene of geneSelectIntersection) {
            const opt = origSelect.querySelector(`option[value="${gene}"]`);
            try {
                opt.setAttribute("selected", "selected");
            } catch (error) {
                // sanity check
                const msg = `Could not add gene ${gene} to gene select.`;
                console.warn(msg);
            }
        }

        geneSelect.update();
        trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to load tags
    })
});

const appendGeneTagButton = (geneTagElt) => {
    // Add delete button
    const deleteBtnElt = document.createElement("button");
    deleteBtnElt.classList.add("delete", "is-small");
    geneTagElt.appendChild(deleteBtnElt);
    deleteBtnElt.addEventListener("click", (event) => {
        // Remove gene from geneSelect
        const gene = event.target.parentNode.textContent;
        const geneSelectElt = document.getElementById("gene_select");
        geneSelectElt.querySelector(`option[value="${gene}"]`).removeAttribute("selected");

        geneSelect.update();
        trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to load tags
    });

    // ? Should i add ellipses for too many genes? Should I make the box collapsable?
}

const clearGenes = (event) => {
    document.getElementById("clear_genes_btn").classList.add("is-loading");
    geneSelect.clear();
    document.getElementById("clear_genes_btn").classList.remove("is-loading");
}

const curatorSpecifcChooseGene = (event) => {
    // Triggered when a gene is selected

    // Delete existing tags
    const geneTagsElt = document.getElementById("gene_tags");
    geneTagsElt.replaceChildren();

    if (!geneSelect.selectedOptions.length) return;   // Do not trigger after initial population

    // Update list of gene tags
    const sortedGenes = geneSelect.selectedOptions.map((opt) => opt.data.value).sort();
    for (const opt in sortedGenes) {
        const geneTagElt = document.createElement("span");
        geneTagElt.classList.add("tag", "is-primary");
        geneTagElt.textContent = sortedGenes[opt];
        appendGeneTagButton(geneTagElt);
        geneTagsElt.appendChild(geneTagElt);
    }

    document.getElementById("gene_tags_c").classList.remove("is-hidden");
    if (!geneSelect.selectedOptions.length) {
        document.getElementById("gene_tags_c").classList.add("is-hidden");
    }

    // Cannot plot if 2+ genes are not selected
    if (geneSelect.selectedOptions.length < 2) {
        document.getElementById("gene_s_failed").classList.remove("is-hidden");
        document.getElementById("gene_s_success").classList.add("is-hidden");
        for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
            plotBtn.disabled = true;
        }
        document.getElementById("continue_to_plot_options").classList.add("is-hidden");
        return;
    }

    // If more than 10 tags, hide the rest and add a "show more" button
    if (geneSelect.selectedOptions.length > 10) {
        const geneTags = geneTagsElt.querySelectorAll("span.tag");
        for (let i = 10; i < geneTags.length; i++) {
            geneTags[i].classList.add("is-hidden");
        }
        // Add show more button
        const showMoreBtnElt = document.createElement("button");
        showMoreBtnElt.classList.add("tag", "button", "is-small", "is-primary", "is-light");
        const numToDisplay = geneSelect.selectedOptions.length - 10;
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


    document.getElementById("gene_s_failed").classList.add("is-hidden");
    document.getElementById("gene_s_success").classList.remove("is-hidden");

    // Force validation check to see if plot button should be enabled
    //trigger(document.querySelector(".js-plot-req"), "change");

    document.getElementById("continue_to_plot_options").classList.remove("is-hidden");

}

const curatorSpecifcCreatePlot = async (plotType) => {
    // Call API route by plot type
    await plotStyle.createPlot(datasetId, analysisObj, userId, colorblindMode);
}

const curatorSpecifcDatasetTreeCallback = async () => {
    // Creates gene select instance that allows for multiple selection
    geneSelect = createGeneSelectInstance("gene_select", geneSelect);
}

const curatorSpecificOnLoad = async () => {
    // Load gene carts
    await loadGeneCarts();
}

const curatorSpecificPlotStyle = (plotType) => {
    // include plotting backend options
    if (genesAsAxisPlots.includes(plotType)) {
        return new GenesAsAxisHandler(plotType);
    } else if (genesAsDataPlots.includes(plotType)) {
        return new GenesAsDataHandler(plotType);
    } else {
        return null;
    }
}

const curatorSpecificPlotTypeAdjustments = (plotType) => {
    return plotType;
}

const curatorSpecificUpdateGeneOptions = async (geneSymbols) => {
    //pass
}

const fetchAvailablePlotTypes = async (session_id, dataset_id, analysis_id) => {
    // Plot types will depend on the number of comparabie categorical conditions
    // Volcano plots must have at least two conditions
    // Quadrant plots must have at least three conditions

    const payload = {session_id, dataset_id, analysis_id};
    try {
        const {data} = await axios.post(`/api/h5ad/${dataset_id}/mg_availableDisplayTypes`, payload);
        if (data.hasOwnProperty("success") && data.success < 1) {
            throw new Error(data?.message || "Could not fetch compatible plot types for this dataset. Please contact the gEAR team.");
        }
        return data;
    } catch (error) {
        logErrorInConsole(error);
        createToast(error.message);
        throw new Error(msg);
    }
}

const fetchDashData = async (plotConfig, datasetId, plot_type, analysis, analysis_owner_id, colorblind_mode)  => {
    // NOTE: gene_symbol already passed to plotConfig
    const payload = { ...plotConfig, plot_type, analysis, analysis_owner_id, colorblind_mode };
    try {
        const { data } = await axios.post(`/api/plot/${datasetId}/mg_dash`, payload);
        if (data?.success < 1) {
            throw new Error (data?.message || "Unknown error.")
        }
        return data
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not create plot for this dataset and parameters. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

/* Fetch gene collections */
const fetchGeneCarts = async (session_id) => {
    const payload = {session_id};
    try {
        const {data} = await axios.post(`/cgi/get_user_gene_carts.cgi`, convertToFormData(payload));
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene collections. You can still enter genes manually.";
        createToast(msg);
        throw new Error(msg);
    }
}

/* Fetch gene collection members */
const fetchGeneCartMwembers = async (session_id, geneCartId) => {
    const payload = { session_id, gene_cart_id: geneCartId };
    try {
        const {data} = await axios.post(`/cgi/get_gene_cart_members.cgi`, convertToFormData(payload));
        const {gene_symbols, success} = data;
        if (!success) {
            throw new Error("Could not fetch gene collection members. You can still enter genes manually.");
        }
        return gene_symbols;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene collection members. You can still enter genes manually.";
        createToast(msg);
        throw new Error(msg);
    }
}

const getCategoryColumns = async () => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    try {
        ({obs_columns: allColumns, obs_levels: levels} = await fetchH5adInfo(datasetId, analysisId));
    } catch (error) {
        document.getElementById("plot_options_s_failed").classList.remove("is-hidden");
        return;
    }
    // Filter out values we don't want of "levels", like "colors"
    for (const key in levels) {
        if (key.includes("_colors")) {
            delete levels[key];
        }
    }
    return Object.keys(levels);
}

// Invert a log function
function invertLogFunction(value, base=10) {
    return base ** value;
}

/* Transform and load gene collection data into a "tree" format */
const loadGeneCarts = async () => {
    try {
        const geneCartData = await fetchGeneCarts(sessionId);
        const carts = {};
        const cartTypes = ['domain', 'user', 'group', 'shared', 'public'];
        let cartsFound = false;

        // Loop through the different types of gene collections and add them to the carts object
        for (const ctype of cartTypes) {
            carts[ctype] = [];

            if (geneCartData[`${ctype}_carts`].length > 0) {
                cartsFound = true;

                for (const item of geneCartData[`${ctype}_carts`]) {
                    carts[ctype].push({value: item.id, text: item.label });
                };
            }
        }

        geneCartTree.domainGeneCarts = carts.domain;
        geneCartTree.userGeneCarts = carts.user;
        geneCartTree.groupGeneCarts = carts.group;
        geneCartTree.sharedGeneCarts = carts.shared;
        geneCartTree.publicGeneCarts = carts.public;
        geneCartTree.generateTree();
        /*if (!cartsFound ) {
            // ? Put some warning if carts not found
            $('#gene_cart_container').show();
        }*/

    } catch (error) {
        document.getElementById("gene_s_failed").classList.remove("is-hidden");
    }

}

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
}


// For a given categorical series (e.g. "celltype"), update the "group" options
const updateGroupOptions = (classSelector, levels) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add categories
        for (const group of levels.sort()) {
            const option = document.createElement("option");
            option.textContent = group;
            option.value = group;
            elt.append(option);
        }
    }

}
document.getElementById("clear_genes_btn").addEventListener("click", clearGenes);