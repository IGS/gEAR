'use strict';

/* These are functions that are common to the "curator" pages, (i.e. single-gene, multi-gene) */

let plotStyle;  // Plot style object

let facetWidget = null;

//let plotConfig = {};  // Plot config that is passed to API or stored in DB
let allColumns = [];
let catColumns = [];
let levels = {};    // categorical columns + groups

let datasetId = null;
let organismId = null;
let analysisObj = null;

let analysisSelect = null;
let plotTypeSelect = null;
let colorscaleSelect = null;

/*
! Quick note -
This code uses a lot of "for...in..." and "for...of..." syntax
"for...in..." is generally used to iterate over object keys or array indexes
"for...of..." is generally used to iterate over array (or collection) values, and will not work on objects

Using "for...of..." over array.forEach means we don't have to transform HTMLCollections into arrays

-- Nice-select2 --
This code uses the nice-select2 library to style select boxes. The "select" is hidden and replaced with a new "div" element.
If you select an option from the "select" element, or if you disable the "select" element,
you must update the nice-select2 instance for it to reflect the change.
To make it a "multiple" select object, add "multiple" to the original "select" element.

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
    configProp2ClassElt = {};   // This will be overridden by subclasses (note: cannot "super" instance properties)

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

    async setupParamValueCopyEvent() {
        throw new Error("You have to implement the method setupParamValueCopyEvent!");
    }

    setupPlotSpecificEvents() {
        throw new Error("You have to implement the method setupPlotSpecificEvents!");
    }

}

/* API Mixin */

const curatorApiCallsMixin = {

    async deleteDisplay(displayId) {
        try {
            await super.deleteDisplay(displayId);
            // Remove display card
            const displayCard = document.getElementById(`${displayId}_display`);
            displayCard.remove();
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not delete this display. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },
    async fetchAggregations(datasetId, analysisId, filters){
        try {
            const data = await super.fetchAggregations(datasetId, analysisId, filters);
            if (data.hasOwnProperty("success") && data.success < 1) {
                throw new Error(data?.message || "Could not fetch number of observations for this dataset. Please contact the gEAR team.");
            }
            const {aggregations, total_count} = data;
            return {aggregations, total_count};
        } catch (error) {
            logErrorInConsole(error);
        }
    },

    async fetchAnalyses (datasetId) {
        try {
            const { public: publicAnalyses, private: privateAnalyses } = await super.fetchAnalyses(datasetId);
            return {publicAnalyses, privateAnalyses};

        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch saved analyses for this dataset. You can still create a plot but it will be based on the original dataset."
            createToast(msg);
            throw new Error(msg);
        }
    },

    async fetchAvailablePlotTypes(datasetId, analysisId, isMultigene=false){
        try {
            const data = await super.fetchAvailablePlotTypes(datasetId, analysisId, isMultigene);
            if (data.hasOwnProperty("success") && data.success < 1) {
                throw new Error(data?.message || "Could not fetch compatible plot types for this dataset. Please contact the gEAR team.");
            }

            // Multigene plot types will depend on the number of comparabie categorical conditions
            // Volcano plots must have at least two conditions
            // Quadrant plots must have at least three conditions

            return data;
        } catch (error) {
            logErrorInConsole(error);
            createToast(error.message);
            throw new Error(msg);
        }
    },

    async fetchDatasets() {
        try {
            return await super.fetchDatasets();
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch datasets. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    async fetchDatasetDisplayImage(datasetId, displayId) {
        try {
            return await super.fetchDatasetDisplayImage(datasetId, displayId);
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch the image preview for this dataset display. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    async fetchDatasetDisplays(datasetId) {
        try {
            // POST due to payload variables being sensitive
            const {user, owner} = await super.fetchDatasetDisplays(datasetId);
            // Filter only the single-gene displays
            if (isMultigene) {
                const userDisplays = user.filter( display => display.plotly_config.hasOwnProperty('gene_symbols'));
                const ownerDisplays = owner.filter( display => display.plotly_config.hasOwnProperty('gene_symbols'));
                return {userDisplays, ownerDisplays};
            }
            const userDisplays = user.filter( display => display.plotly_config.hasOwnProperty('gene_symbol'));
            const ownerDisplays = owner.filter( display => display.plotly_config.hasOwnProperty('gene_symbol'));
            return {userDisplays, ownerDisplays};
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch the saved displays for this dataset. Please contact the gEAR team."
            createToast(msg);
            return [];  // Send an empty list of displays
        }
    },

    async fetchDefaultDisplay(datasetId) {
        try {
            // POST due to payload variables being sensitive
            const {default_display_id} =  await super.fetchDefaultDisplay(datasetId);
            return default_display_id;
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch default display for this dataset. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    async fetchGeneSymbols(datasetId, analysisId) {
        try {
            const data = await super.fetchGeneSymbols(datasetId, analysisId);
            return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch gene symbols for this dataset. Please contact the gEAR team."
            createToast(msg);
            return [];
        }
    },

    async fetchH5adInfo(datasetId, analysisId) {
        try {
            const {obs_columns, obs_levels} = await super.fetchH5adInfo(datasetId, analysisId);
            return { obs_columns, obs_levels };
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch H5AD observation data for this dataset. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    async saveDatasetDisplay(datasetId, displayId, label, plotType, plotConfig){
        // NOTE: Saving all displays as new displays (clone) instead of overwriting. User can always delete excess displays
        try {
            const {display_id, success} = await super.saveDatasetDisplay(datasetId, displayId, label, plotType, plotConfig);
            if (!success) {
                throw new Error("Could not save this new display. Please contact the gEAR team.");
            }

            // Make new display card and make it the default display
            renderUserDisplayCard({id: display_id, label, plot_type, plotly_config: plotConfig}, display_id);

            return display_id;
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not save this new display. Please contact the gEAR team."
            createToast(error);
            throw new Error(msg);
        }
    },

    async saveDefaultDisplay(displayId) {
        try {
            const {success} = await super.saveDefaultDisplay(datasetId, displayId, isMultigene);
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

}
Object.setPrototypeOf(curatorApiCallsMixin, apiCallsMixin);


const datasetTree = new DatasetTree({
    element: document.getElementById("dataset_tree")
    , searchElement: document.getElementById("dataset_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current_dataset_c").classList.remove("is-hidden");
        document.getElementById("current_dataset").textContent = e.node.title;
        document.getElementById("current_dataset_post").textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        datasetId = newDatasetId;

        // Click to get to next step
        document.getElementById("load_plot_s").click();
        document.getElementById('new_display').classList.add("is-loading");

        // Clear "success/failure" icons
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.classList.add("is-hidden");
        }
        for (const elt of document.getElementsByClassName("js-step-failure")) {
            elt.classList.add("is-hidden");
        }

        document.getElementById("dataset_s_success").classList.remove("is-hidden");

        // displays
        const {userDisplays, ownerDisplays} = await curatorApiCallsMixin.fetchDatasetDisplays(datasetId);
        let defaultDisplayId;
        try {
            defaultDisplayId = await fetchDefaultDisplay(datasetId);
        } catch (error) {
            defaultDisplayId = -1;  // Cannot make any display a default.
        }
        renderDisplayCards(userDisplays, ownerDisplays, defaultDisplayId);
        document.getElementById('new_display').classList.remove("is-loading");
        document.getElementById('new_display').disabled = false;

        // Clear (and update) options within nice-select2 structure.
        // Not providing the object will duplicate the nice-select2 structure
        analysisObj = null;
        analysisSelect = createAnalysisSelectInstance("analysis_select", analysisSelect);
        plotTypeSelect = createPlotTypeSelectInstance("plot_type_select", plotTypeSelect);

        // Call any curator-specific callbacks
        await curatorSpecifcDatasetTreeCallback();
    })
});

const analysisSelectUpdate = async () => {
    try {
        const {publicAnalyses, privateAnalyses} = await curatorApiCallsMixin.fetchAnalyses(datasetId);
        updateAnalysesOptions(privateAnalyses, publicAnalyses);
        document.getElementById("analysis_type_select_c_success").classList.remove("is-hidden");   // Default analysis is good
    } catch (error) {
        // Show failure state things.
        document.getElementById("plot_type_s_failed").classList.remove("is-hidden");
        document.getElementById("analysis_type_select_c_failed").classList.remove("is-hidden");
        document.getElementById('new_display').classList.remove("is-loading"); // Don't give impression display is still loading
    } finally {
        document.getElementById("load_plot_s_success").classList.remove("is-hidden");
    }

}

const chooseAnalysis = async (event) => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const analysisText = (analysisId.length) ? analysisId : "Primary Analysis";

    // Display current selected analysis
    document.getElementById("current_analysis").textContent = analysisText;
    document.getElementById("current_analysis_post").textContent = analysisText;


    // NOTE: For now, we can just pass analysis id only to tSNE and be fine
    // Any private dataset will belong to our user. Any public datasets can be found by the API "get_analysis" code.
    if (analysisId) {
        analysisObj = {id: analysisId};
    }

    if (analysisSelect.data.length < 5) return; // Have not retrieved analyses from API yet

    if (analysisId) {
        await Promise.all([
            plotTypeSelectUpdate(analysisId)
            , geneSelectUpdate(analysisId)
        ]);

        // Create facet widget
        facetWidget = await createFacetWidget(datasetId, analysisId, {});
    }
}

const chooseGene = async (event) => {
    // Each page will deal with this separately
    curatorSpecifcChooseGene(event);
}

/* New display has been chosen, so display analysis and plot type options */
const chooseNewDisplay = async (event) => {
    document.getElementById('new_display').classList.add("is-loading");
    document.getElementById("analysis_select").disabled = false;

    document.getElementById("plot_type_select").disabled = false;

    // update genes, analysis, and plot type selects in parallel
    await Promise.all([
        geneSelectUpdate(),
        analysisSelectUpdate(),
        plotTypeSelectUpdate()      // NOTE: Believe updating "disabled" properties triggers the plotTypeSelect "change" element

    ]);

    document.getElementById('new_display').classList.remove("is-loading");
    document.getElementById("plot_type_s").click();
}

const choosePlotType = async (event) => {
    if (!plotTypeSelect.selectedOptions.length) return;   // Do not trigger after setting disable/enable on options

    // Do not display if default opt is chosen
    const plotType = getSelect2Value(plotTypeSelect)
    if (plotType === "nope") {
        document.getElementById("plot_type_select_c_success").classList.add("is-hidden");
        document.getElementById("plot_type_s_success").classList.add("is-hidden");
        return;
    }

    document.getElementById("plot_type_select_c_failed").classList.add("is-hidden");
    document.getElementById("plot_type_select_c_success").classList.remove("is-hidden");

    document.getElementById("plot_type_s_failed").classList.add("is-hidden");
    document.getElementById("plot_type_s_success").classList.remove("is-hidden");

    document.getElementById("plot_options_s_failed").classList.add("is-hidden");
    document.getElementById("plot_options_s_success").classList.add("is-hidden");

    // Display current selected plot type
    document.getElementById("current_plot_type_c").classList.remove("is-hidden");
    document.getElementById("current_plot_type").textContent = plotType;

    // Create facet widget, which will refresh filters
    facetWidget = await createFacetWidget(datasetId, null, {});
    document.getElementById("facet_content").classList.remove("is-hidden");
    document.getElementById("selected_facets").classList.remove("is-hidden");

    // Reset sortable lists
    document.getElementById("order_section").classList.add("is-hidden");
    document.getElementById("order_container").replaceChildren();


    await includePlotParamOptions();
    document.getElementById("gene_s").click();
}


const cloneDisplay = async (event, display) => {

    const cloneId = event.target.id;
    document.getElementById(cloneId).classList.add("is-loading");

    // Populate gene select element
    // Will be overwritten if an analysis was in config
    try {
        const geneSymbols = await curatorApiCallsMixin.fetchGeneSymbols(datasetId, null);
        updateGeneOptions(geneSymbols);
    } catch (error) {
        document.getElementById("gene_s_failed").classList.remove("is-hidden");
    }

    document.getElementById("analysis_select").disabled = false;
    document.getElementById("plot_type_select").disabled = false;

    // analyses
    await analysisSelectUpdate();
    // TODO update analysis select with config analysis


    // plot types

    // Read clone config to populate analysis, plot type, gnee and plot-specific options
    let plotType = display.plot_type;
    plotType = curatorSpecificPlotTypeAdjustments(plotType);

    try {
        const availablePlotTypes = await curatorApiCallsMixin.fetchAvailablePlotTypes(datasetId, undefined, isMultigene);
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
        document.getElementById("plot_type_s_failed").classList.remove("is-hidden");
        document.getElementById("plot_type_select_c_failed").classList.remove("is-hidden");
        document.getElementById("plot_type_s_success").classList.add("is-hidden");
        document.getElementById("plot_type_select_c_success").classList.add("is-hidden");

        return;
    } finally {
        document.getElementById(cloneId).classList.remove("is-loading");
    }

    // Choose gene from config
    const config = display.plotly_config;
    if (isMultigene) {
        // Update gene_select with genes from config
        const geneSymbols = config.gene_symbols;
        for (const geneSymbol of geneSymbols) {
            setSelectBoxByValue("gene_select", geneSymbol);
        }
    } else {
        setSelectBoxByValue("gene_select", config.gene_symbol);
    }
    geneSelect.update();
    trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to set the other select2 (single-gene only)

    try {
        plotStyle.cloneDisplay(config);
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not clone this display. Please contact the gEAR team."
        createToast(msg);
        return;
    }

    // Mark plot params as success
    document.getElementById("plot_options_s_success").classList.remove("is-hidden");

    // Click "submit" button to load plot
    document.getElementById("plot_btn").click();    // updates geneSelectPost by triggered "click" event

}

const createAnalysisSelectInstance = (idSelector, analysisSelect=null) => {
    // If object exists, just update it with the revised data and return
    if (analysisSelect) {
        analysisSelect.update();
        return analysisSelect;
    }

    return NiceSelect.bind(document.getElementById(idSelector), {
        placeholder: 'Select an analysis.',
        allowClear: true,
    });
}

// Create the gradient for the canvas element using a given colorscale's information and the element HTML object
const createCanvasGradient = (elem) => {
    // Get ID of canvas element and remove "gradient_" from the name
    const id = elem.id.replace("gradient_", "");
    // Get the colorscale info for the given element (object is in plot_display_config.js)
    const data = paletteInformation[id];

    const ctx = elem.getContext("2d");  // canvas element
    const grid = ctx.createLinearGradient(0, 0, elem.width, 0);    // Fill across but not down
    // Add the colors to the gradient
    for (const color of data) {
        grid.addColorStop(color[0], color[1]);
    }
    // Fill the canvas with the gradient
    ctx.fillStyle = grid;
    ctx.fillRect(0, 0, elem.width, 20);
}

const createCanvasScale = (elem) => {
    // Get ID of canvas element and remove "gradient_" from the name
    const id = elem.id.replace("gradient_", "");
    // Get the colorscale info for the given element (object is in plot_display_config.js)
    const data = paletteInformation[id];

    const elemWidth = elem.width;
    const ctx = elem.getContext("2d");  // canvas element
    // Add the colors to the scale
    const { length } = data;
    const width = elemWidth/length;   // 150 is length of canvas
    for (const color of data) {
        ctx.fillStyle = color[1];
        // The length/length+1 is to account for the fact that the last color has a value of 1.0
        // Otherwise the last color would be cut off
        const x = color[0] * (length/(length+1)) * elemWidth;
        ctx.fillRect(x, 0, width, 20);
    }
}

const createColorscaleSelectInstance = (idSelector, colorscaleSelect=null) => {
    // If object exists, just update it with the revised data and return
    if (colorscaleSelect) {
        colorscaleSelect.update();
        return colorscaleSelect;
    }

    return NiceSelect.bind(document.getElementById(idSelector), {
        placeholder: 'Choose a color palette',
        width: '50%',
        minimumResultsForSearch: -1
    });
}

const createFacetWidget = async (datasetId, analysisId, filters) => {
    document.getElementById("selected_facets_loader").classList.remove("is-hidden")

    const {aggregations, total_count:totalCount} = await curatorApiCallsMixin.fetchAggregations(datasetId, analysisId, filters);
    document.getElementById("num_selected").textContent = totalCount;


    const facetWidget = new FacetWidget({
        aggregations,
        filters,
        onFilterChange: async (filters) => {
            if (filters) {
                try {
                    const {aggregations, total_count:totalCount} = await curatorApiCallsMixin.fetchAggregations(datasetId, analysisId, filters);
                    facetWidget.updateAggregations(aggregations);
                    document.getElementById("num_selected").textContent = totalCount;
                } catch (error) {
                    logErrorInConsole(error);
                }
            } else {
                // Save an extra API call
                facetWidget.updateAggregations(facetWidget.aggregations);
            }
            // Sortable lists need to reflect groups filtered out or unfiltered
            updateOrderSortable();
        }
    });
    document.getElementById("selected_facets_loader").classList.add("is-hidden")
    return facetWidget;
}

const createGeneSelectInstance = (idSelector, geneSelect=null) => {
    // NOTE: Updating the list of genes can be memory-intensive if there are a lot of genes
    // and (I've noticed) if multiple select2 elements for genes are present.

    // If object exists, just update it with the revised data and return
    if (geneSelect) {
        geneSelect.update();
        return geneSelect;
    }

    return NiceSelect.bind(document.getElementById(idSelector), {
        placeholder: 'To search, start typing a gene name',
        searchtext: 'To search, start typing a gene name',
        searchable: true,
        allowClear: true,
    });
}

const createPlotTypeSelectInstance = (idSelector, plotTypeSelect=null) => {
    // If object exists, just update it with the revised data and return
    if (plotTypeSelect) {
        plotTypeSelect.update();
        return geneSelect;
    }

    // Initialize fixed plot types
    return NiceSelect.bind(document.getElementById(idSelector), {
        placeholder: 'Choose how to plot',
        minimumResultsForSearch: -1
    });
}

const createPlot = async (event) => {

    const plotType = getSelect2Value(plotTypeSelect);

    const plotBtns = document.getElementsByClassName("js-plot-btn");

    // Set loading
    for (const plotBtn of plotBtns) {
        plotBtn.classList.add("is-loading");
    }

    plotStyle.populatePlotConfig();

    // Add gene or genes to plot config
    if (isMultigene) {
        plotStyle.plotConfig["gene_symbols"] = geneSelect.selectedOptions.map(e => e.data.value);
    } else {
        plotStyle.plotConfig["gene_symbol"] = getSelect2Value(geneSelect);
    }

    await curatorSpecifcCreatePlot(plotType);


    // Stop loader
    for (const plotBtn of plotBtns) {
        plotBtn.classList.remove("is-loading");
    }

    // Hide this view
    document.getElementById("content_c").classList.add("is-hidden");
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").classList.remove("is-hidden");

}

/* Creates a Toast-style message in the upper-corner of the screen. */
const createToast = (msg, levelClass="is-danger") => {
    const template = `
    <div class="notification js-toast ${levelClass} animate__animated animate__fadeInUp animate__faster">
        <button class="delete"></button>
        ${msg}
    </div>
    `
    const html = generateElements(template);

    const numToasts = document.querySelectorAll(".js-toast.notification").length;

    if (document.querySelector(".js-toast.notification")) {
        // If .js-toast notifications are present, append under final notification
        // This is to prevent overlapping toast notifications
        document.querySelector(".js-toast.notification:last-of-type").insertAdjacentElement("afterend", html);
        // Position new toast under previous toast with CSS
        html.style.setProperty("top", `${(numToasts * 70) + 30}px`);
    } else {
        // Otherwise prepend to top of main content
        document.getElementById("main_c").prepend(html);
    }

    // This should get the newly added notification since it is now the first
    html.querySelector(".js-toast.notification .delete").addEventListener("click", (event) => {
        console.log(event.target);
        const notification = event.target.closest(".js-toast.notification");
        notification.remove(notification);
    });

    // For a success message, remove it after 3 seconds
    if (levelClass === "is-success") {
        const notification = document.querySelector(".js-toast.notification:last-of-type");
        notification.classList.remove("animate__fadeInUp");
        notification.classList.remove("animate__faster");
        notification.classList.add("animate__fadeOutDown");
        notification.classList.add("animate__slower");
    }
}

const disableCheckboxLabel = (checkboxElt, state) => {
    // if parent element is a .checkbox class, disable it too (uses Bulma CSS styling)
    // Meant for checkboxes where the label is also a clickable element
    // NOTE: ".disable" attribute only applies to certain elements (https://www.w3schools.com/tags/att_disabled.asp)
    if (checkboxElt.parentElement.classList.contains("checkbox")) {
        if (state) {
            checkboxElt.parentElement.setAttribute("disabled", "");
        } else {
            checkboxElt.parentElement.removeAttribute("disabled");
        }
    }
}

// Create the template for the colorscale select2 option dropdown
const formatColorscaleOptionText = (option, text, isContinuous=false) => {

    const fragment = document.createDocumentFragment();
    const canvas = document.createElement("canvas");
    canvas.id = `gradient_${option.value}`;
    canvas.width = 100;
    canvas.height = 20;
    canvas.classList.add("js-palette-canvas");

    // BUG: Gradient does not show in the nice-select2 rendered option
    if (isContinuous) {
        createCanvasGradient(canvas);
    } else {
        createCanvasScale(canvas);
    }
    fragment.append(canvas);

    const text_span = document.createElement("span");
    text_span.classList.add("pl-1");
    text_span.textContent = text;
    fragment.append(text_span);
    return fragment;
}

const geneSelectUpdate = async (analysisId=null) => {
    // Populate gene select element
    try {
        const geneSymbols = await curatorApiCallsMixin.fetchGeneSymbols(datasetId, analysisId);
        updateGeneOptions(geneSymbols); // Come from curator specific code
    } catch (error) {
        document.getElementById("gene_s_failed").classList.remove("is-hidden");
    }
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
        const serialized = sortable(`#${CSS.escape(series)}_order_list`, 'serialize')[0].items;
        // Sort by "sortable" index position
        order[series] = serialized.map((val) => val.label);
    }
    return order;
}

const getSelect2Value = (select) => {
    // Get value from select2 element
    return select.selectedOptions[0].data.value;
}

const includeHtml = async (url) => {
    const preResponse = await fetch(url, {cache: "reload"});
    return await preResponse.text();
}

/* Load custom plot options */
const includePlotParamOptions = async () => {
    const plotType = getSelect2Value(plotTypeSelect);

    // New plot... so disable plot button
    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
        plotBtn.disabled = true;
    }

    plotStyle = curatorSpecificPlotStyle(plotType);
    if (!plotStyle) {
        console.warn(`Plot type ${plotType} not recognized.`)
        document.getElementById("plot_type_s_failed").classList.remove("is-hidden");
        document.getElementById("plot_type_s_success").classList.add("is-hidden");
        return;
    }
    document.getElementById("plot_type_s_failed").classList.add("is-hidden");


    // NOTE: Changing plots within the same plot style will clear the plot config as fresh templates are loaded
    await plotStyle.loadPlotHtml();

    // NOTE: Events are triggered in the order they are regstered.
    // We want to trigger a param sync event before the plot requirement validation
    // (since it checks every element in the js-plot-req class)
    for (const classSelector of Object.keys(plotStyle.classElt2Prop)) {
        // Ensure pre- and post- plot view params are synced up
        setupParamValueCopyEvent(classSelector);
    }
    plotStyle.setupParamValueCopyEvent();   // handle some copy events that could not be handled in the loop above
    setupValidationEvents();        // Set up validation events required to plot at a minimum
    await plotStyle.setupPlotSpecificEvents()       // Set up plot-specific events

}

// Load colorscale select2 object and populate with data
const loadColorscaleSelect = (isContinuous=false) => {

    let filteredPalettes = availablePalettes;

    // If plot that uses continuous colorscales is chosen, then filter availablePalettes to only those for continuous plots
    // If not continuous, then filter to only those for discrete plots
    filteredPalettes = isContinuous ?
        availablePalettes.filter((type) => type.continuous) :
        availablePalettes.filter((type) => !type.continuous);

    for (const palette of filteredPalettes) {
        const label = palette.label;
        const optgroup = document.createElement("optgroup");
        optgroup.setAttribute("label", label);
        for (const option of palette.options) {
            const optionElt = document.createElement("option");
            optionElt.value = option.value;
            // Add canvas element information to option, which is converted to innerHTML by nice-select2._renderItem
            optionElt.append(formatColorscaleOptionText(optionElt, option.text, isContinuous));
            optgroup.append(optionElt);
        }
        document.getElementById("color_palette_post").append(optgroup);
    }


    // set default to purples
    setSelectBoxByValue("color_palette_post", "purp");

    return;

    //colorscaleSelect = createColorscaleSelectInstance("color_palette_post", colorscaleSelect);

}

/* Transform and load dataset data into a "tree" format */
const loadDatasetTree = async () => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const datasetData = await curatorApiCallsMixin.fetchDatasets();

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
        datasetTree.userDatasets = userDatasets;
        datasetTree.sharedDatasets = sharedDatasets;
        datasetTree.domainDatasets = domainDatasets;
        datasetTree.generateTree();
    } catch (error) {
        document.getElementById("dataset_s_failed").classList.remove("is-hidden");
    }

}

// Update the options within the plot type select element, specifically the disabled state
const plotTypeSelectUpdate = async (analysisId=null) => {
    // NOTE: Believe updating "disabled" properties triggers the plotTypeSelect "change" element
    try {
        plotTypeSelect.clear();
        const availablePlotTypes = await curatorApiCallsMixin.fetchAvailablePlotTypes(datasetId, analysisId, isMultigene);
        for (const plotType in availablePlotTypes) {
            const isAllowed = availablePlotTypes[plotType];
            setPlotTypeDisabledState(plotType, isAllowed);
        }
        plotTypeSelect.update();
    } catch (error) {
        document.getElementById("plot_type_s_failed").classList.remove("is-hidden");
        document.getElementById("plot_type_select_c_failed").classList.remove("is-hidden");
        document.getElementById("plot_type_s_success").classList.add("is-hidden");
        document.getElementById("plot_type_select_c_success").classList.add("is-hidden");;
    }
}

const renderDisplayCards = async (userDisplays, ownerDisplays, defaultDisplayId) => {
    const userDisplaysElt = document.getElementById("user_displays");
    const ownerDisplaysElt = document.getElementById("owner_displays");

    // Empty existing displays
    userDisplaysElt.replaceChildren();
    ownerDisplaysElt.replaceChildren();

    // Add titles to each section if there are displays
    if (userDisplays.length) {
        const userTitle = generateElements(`<p class="has-text-weight-bold is-underlined column is-full">Your Displays</p>`);
        userDisplaysElt.append(userTitle);
    }

    if (ownerDisplays.length) {
    const ownerTitle = generateElements(`<p class="has-text-weight-bold is-underlined column is-full">Displays by Dataset Owner</p>`);
    ownerDisplaysElt.append(ownerTitle);
    }

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

    // Start with a fresh template
    const orderElt = document.getElementById(`${series}_order`);
    if (orderElt) {
        orderElt.remove();
    }

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
        // If filters are present and group is not in filters, skip
        if (facetWidget.filters.hasOwnProperty(series) && !facetWidget.filters[series].includes(group)) continue;

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
        displayUrl = await curatorApiCallsMixin.fetchDatasetDisplayImage(datasetId, display.id);
    } catch (error) {
        displayUrl = "/img/dataset_previews/missing.png";
    }
    const geneSymbol = display.plotly_config.gene_symbol;

    const label = display.label || `Unnamed ${display.plot_type} display`;

    let geneCardContent;
    if (isMultigene) {
        // Card content should be number of genes
        const numGenes = display.plotly_config.gene_symbols.length;
        geneCardContent = `<div class="card-content">
            <p class="subtitle">Number of genes: ${numGenes}</p>
        </div>`
    } else {
        // Card content should be gene symbol
        geneCardContent = `<div class="card-content">
            <p class="subtitle">Gene: ${geneSymbol}</p>
        </div>`
    }

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
                        ${geneCardContent}
                        <footer class="card-footer buttons">
                            <p class="card-footer-item is-paddingless">
                                <button class="js-display-default button is-responsive is-fullwidth is-primary" id="${display.id}_default">Set as Default</button>
                            </p>
                            <p class="card-footer-item is-paddingless">
                                <button class="button is-fullwidth is-responsive is-primary" id="${display.id}_clone">Clone</button>
                            </p>
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
        displayUrl = await curatorApiCallsMixin.fetchDatasetDisplayImage(datasetId, display.id);
    } catch (error) {
        displayUrl = "/img/dataset_previews/missing.png";
    }

    const geneSymbol = display.plotly_config.gene_symbol;

    const label = display.label || "Unnamed display"; // Added text to keep "p" tag from collapsing

    let geneCardContent;
    if (isMultigene) {
        // Card content should be number of genes
        const numGenes = display.plotly_config.gene_symbols.length;
        geneCardContent = `<div class="card-content">
            <p class="subtitle">Number of genes: ${numGenes}</p>
        </div>`
    } else {
        // Card content should be gene symbol
        geneCardContent = `<div class="card-content">
            <p class="subtitle">Gene: ${geneSymbol}</p>
        </div>`
    }

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
                        ${geneCardContent}
                        <footer class="card-footer ">
                            <p class="card-footer-item is-paddingless">
                                <button class="js-display-default button is-responsive is-fullwidth is-primary" id="${display.id}_default">Set as Default</button>
                            </p>
                            <p class="card-footer-item is-paddingless">
                                <button class="button is-fullwidth is-responsive is-primary" id="${display.id}_clone">Clone</button>
                            </p>
                            <p class="card-footer-item is-paddingless">
                                <button class="button is-fullwidth is-responsive is-danger" id="${display.id}_delete">Delete</button>
                            </p>
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
    document.getElementById(`${display.id}_delete`).addEventListener("click", (event) => curatorApiCallsMixin.deleteDisplay(display.id));
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
 * Modified to set value and "selected" so nice-select2 extractData() will catch it
 * Taken from https://stackoverflow.com/a/20662180
 * @param eid Element ID
 * @param eval Element value
 */
const setSelectBoxByValue = (eid, val) => {
    const elt = document.getElementById(eid);
    for (const i in elt.options) {
        if (elt.options[i].value === val) {
            // By using "selected" instead of "value", we can account for the "multiple" property
            elt.options[i].setAttribute("selected", true);
            return;
        }
    }
}

/* Ensure all elements in this class have the same value */
const setupParamValueCopyEvent = (classSelector) => {
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

/* Setup a fail-fast validation trigger. */
const setupValidationEvents = () => {
    const validationElts = document.getElementsByClassName("js-plot-req");
    for (const elt of validationElts ) {
        elt.addEventListener("change", validateRequirements);
    }
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

    // Update select2
    analysisSelect.update();
}

const updateGeneOptions = (geneSymbols) => {

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

    curatorSpecificUpdateGeneOptions(geneSymbols);

    // Update the nice-select2 element to reflect this.
    // This function is always called in the 1st view, so only update that
    geneSelect.update();

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

    for (const series of seriesSet) {
        // If series already exists, it will be removed before rendered
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
        orderSection.classList.add("is-hidden");
        return;
    }

    orderSection.classList.remove("is-hidden");

}

const validateRequirements = (event) => {
    const elt = event.target;
    // Reset "status" classes
    elt.classList.remove("is-success", "is-danger");
    if (elt.value) {
        elt.parentElement.classList.remove("is-danger");
        elt.parentElement.classList.add("is-success");

        const validationElts = document.getElementsByClassName("js-plot-req");

        // If every validation param has been filled out, it's OK to plot
        // NOTE: need to ensure pre- and post- param elements are filled before this function is called
        if ([...validationElts].every(element => element.value)) {
            for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
                plotBtn.disabled = false;
            }
            document.getElementById("plot_options_s_success").classList.remove("is-hidden");
            document.getElementById("plot_options_s_failed").classList.add("is-hidden");
        }
        return;
    }

    // Required paramater has no value. Indicate it and disable plot buttons
    elt.parentElement.classList.add("is-danger");
    elt.parentElement.classList.remove("is-success");

    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
        plotBtn.disabled = true;
    }
    document.getElementById("plot_options_s_success").classList.add("is-hidden");
    document.getElementById("plot_options_s_failed").classList.remove("is-hidden");
}

document.getElementById("new_display").addEventListener("click", chooseNewDisplay);
document.getElementById("analysis_select").addEventListener("change", chooseAnalysis);
document.getElementById("plot_type_select").addEventListener("change", choosePlotType);

const geneSelectElts = document.querySelectorAll("select.js-gene-select");
for (const geneSelectElt of geneSelectElts) {
    geneSelectElt.addEventListener("change", chooseGene);
}

const plotBtns = document.getElementsByClassName("js-plot-btn");
for (const plotBtn of plotBtns) {
    plotBtn.addEventListener("click", createPlot);
}

document.getElementById("save_json_config").addEventListener("click", () => {
    // Config plot configuration to JSON for sharing (or passing to API by hand)
    const blob = new Blob([JSON.stringify({...plotStyle.plotConfig, plot_type:plotStyle.plotType})]);
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

document.getElementById("save_display_btn").addEventListener("click", async (event) => {
    // Save new plot display.
    const label = document.getElementById("new_display_label").value;
    event.target.classList.add("is-loading");
    try {
        const displayId = await curatorApiCallsMixin.saveDatasetDisplay(datasetId, null, label, plotStyle.plotType, plotStyle.plotConfig);
        createToast("Display saved.", "is-success");

        if (document.getElementById("make_default_display_check").checked) {
            saveDefaultDisplay(displayId);
        }
    } catch (error) {
        //pass - handled in functions
    } finally {
        event.target.classList.remove("is-loading");
    }
});

document.getElementById("edit_params").addEventListener("click", (event) => {
    event.target.classList.add("is-loading");
    // Hide this view
    document.getElementById("content_c").classList.remove("is-hidden");
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").classList.add("is-hidden");

    event.target.classList.remove("is-loading");
})

// Set up .js-post-plot-collapsable to collapse and uncollapse the content element
const collapsableElts = document.getElementsByClassName("js-collapsable-trigger");
for (const classElt of collapsableElts) {
    classElt.addEventListener("click", (event) => {
        // find the sibling .is-collapsable element and toggle its "is-hidden" class
        // Since the user could click the icon, we get closest parent element with class .js-collapsable-trigger
        const elt = event.target.closest(".js-collapsable-trigger")
        const contentElt = elt.parentElement.querySelector(".js-collapsable-content");
        contentElt.classList.toggle("is-hidden");

        // switch toggle icon to up or down
        const iconElt = elt.querySelector("span.icon.is-pulled-right");
        if (iconElt.innerHTML.trim() === '<i class="mdi mdi-chevron-down"></i>') {
            iconElt.innerHTML = '<i class="mdi mdi-chevron-up"></i>';
        } else {
            iconElt.innerHTML = '<i class="mdi mdi-chevron-down"></i>';
        }
    });
}

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

    curatorSpecificNavbarUpdates();

    const sessionId = CURRENT_USER.session_id;
    if (! sessionId ) {
        createToast("Not logged in so saving displays is disabled.");
        document.getElementById("save_display_btn").disabled = true;
    }

	try {
		await loadDatasetTree()
        // If brought here by the "gene search results" page, curate on the dataset ID that referred us
        const urlParams = new URLSearchParams(window.location.search);
        if (urlParams.has("dataset_id")) {
            const linkedDatasetId = urlParams.get("dataset_id");
            try {
                // find DatasetTree node and trigger "activate"
                const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
                foundNode.setActive(true);
                datasetId = linkedDatasetId;
            } catch (error) {
                createToast(`Dataset id ${linkedDatasetId} was not found as a public/private/shared dataset`);
                throw new Error(error);
            }
        }
	} catch (error) {
		logErrorInConsole(error);
	}

    // Load any script-specific code
    curatorSpecificOnLoad();

};
