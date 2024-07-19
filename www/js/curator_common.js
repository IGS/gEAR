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

/**
 * Represents a base class for handling plots.
 * @class
 */
class PlotHandler {
    constructor() {
        // Check if this is an abstract class
        // (new PlotStyle() will fail)
        if (new.target === PlotHandler) throw new TypeError("Cannot construct PlotHandler instances directly");
    }

    // Certain groups like "display order" and "colors" and "vlines" will be custom-handled
    classElt2Prop = {}; // This will be overridden by subclasses
    configProp2ClassElt = {};   // This will be overridden by subclasses (note: cannot "super" instance properties)

    /**
     * Creates a clone of the plot display.
     * @abstract
     * @throws {Error} You have to implement the method cloneDisplay!
     */
    cloneDisplay() {
        throw new Error("You have to implement the method cloneDisplay!");
    }

    /**
     * Asynchronously creates the plot.
     * @abstract
     * @throws {Error} You have to implement the method createPlot!
     */
    async createPlot() {
        throw new Error("You have to implement the method createPlot!");
    }

    /**
     * Asynchronously loads the plot HTML.
     * @abstract
     * @throws {Error} You have to implement the method loadPlotHtml!
     */
    async loadPlotHtml() {
        throw new Error("You have to implement the method loadPlotHtml!");
    }

    /**
     * Populates the plot configuration.
     * @abstract
     * @throws {Error} You have to implement the method populatePlotConfig!
     */
    populatePlotConfig() {
        throw new Error("You have to implement the method populatePlotConfig!");
    }

    /**
     * Sets up the event for copying parameter values.
     * @abstract
     * @throws {Error} You have to implement the method setupParamValueCopyEvent!
     */
    async setupParamValueCopyEvent() {
        throw new Error("You have to implement the method setupParamValueCopyEvent!");
    }

    /**
     * Sets up plot-specific events.
     * @abstract
     * @throws {Error} You have to implement the method setupPlotSpecificEvents!
     */
    setupPlotSpecificEvents() {
        throw new Error("You have to implement the method setupPlotSpecificEvents!");
    }
}

/* API Mixin */

/**
 * Mixin containing various API call methods for the curator module.
 * @mixin
 */
const curatorApiCallsMixin = {

    /**
     * Deletes a display by its ID.
     *
     * @param {string} displayId - The ID of the display to delete.
     * @throws {Error} If the display cannot be deleted.
     */
    async deleteDisplay(displayId) {
        try {
            await super.deleteDisplay(displayId);
            // Remove display card
            const displayCard = document.getElementById(`${displayId}-display`);
            displayCard.remove();
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not delete this display. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    /**
     * Fetches aggregations for a dataset and analysis.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis.
     * @param {object} filters - The filters to apply to the aggregations.
     * @returns {Promise<{aggregations: object, total_count: number}>} The fetched aggregations and total count.
     */
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

    /**
     * Fetches the analyses for a given dataset.
     * @param {string} datasetId - The ID of the dataset.
     * @returns {Promise<{publicAnalyses: Array, privateAnalyses: Array}>} - The fetched public and private analyses.
     * @throws {Error} - If the analyses cannot be fetched.
     */
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

    /**
     * Fetches the available plot types for a given dataset and analysis.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis.
     * @param {boolean} [isMultigene=false] - Flag indicating if the plot types are for multigene analysis.
     * @returns {Promise<Object>} - A promise that resolves to the available plot types data.
     * @throws {Error} - If there is an error fetching the plot types.
     */
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
            throw error;
        }
    },

    /**
     * Fetches datasets.
     * @returns {Promise<any>} A promise that resolves with the fetched datasets.
     * @throws {Error} If the datasets cannot be fetched.
     */
    async fetchAllDatasets() {
        try {
            return await super.fetchAllDatasets();
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch datasets. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    /**
     * Fetches the display image for a dataset.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} displayId - The ID of the display.
     * @returns {Promise} - A promise that resolves with the fetched display image.
     * @throws {Error} - If the image preview cannot be fetched.
     */
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

    /**
     * Fetches the dataset displays for a given dataset ID.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @returns {Promise<{userDisplays: Array, ownerDisplays: Array}>} - The user displays and owner displays.
     * @throws {Error} - If there is an error fetching the displays.
     */
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

    /**
     * Fetches the default display for a dataset.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @returns {Promise<string>} The ID of the default display.
     * @throws {Error} If the default display cannot be fetched.
     */
    async fetchDefaultDisplay(datasetId) {
        try {
            // POST due to payload variables being sensitive
            const {default_display_id} =  await super.fetchDefaultDisplay(datasetId, isMultigene);
            return default_display_id;
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch default display for this dataset. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },


    /**
     * Fetches gene symbols for a given dataset and analysis.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis.
     * @returns {Promise<string[]>} - A promise that resolves to an array of unique gene symbols.
     * @throws {Error} - If the gene symbols cannot be fetched, an error is thrown.
     */
    async fetchGeneSymbols(datasetId, analysisId) {
        try {
            const data = await super.fetchGeneSymbols(datasetId, analysisId);
            return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not fetch gene symbols for this dataset. Please contact the gEAR team."
            createToast(msg);
            throw new Error(msg);
        }
    },

    /**
     * Fetches H5AD information for a dataset and analysis.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis.
     * @returns {Promise<{obs_columns: any, obs_levels: any}>} The observation columns and levels.
     * @throws {Error} If the H5AD observation data cannot be fetched.
     */
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

    /**
     * Saves a dataset display as a new display.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} displayId - The ID of the display.
     * @param {string} label - The label of the display.
     * @param {string} plotType - The type of plot for the display.
     * @param {object} plotConfig - The configuration for the plot.
     * @returns {string} The ID of the saved display.
     * @throws {Error} If the new display could not be saved.
     */
    async saveDatasetDisplay(datasetId, displayId, label, plotType, plotConfig){
        // NOTE: Saving all displays as new displays (clone) instead of overwriting. User can always delete excess displays
        if (analysisObj) {
            plotConfig["analysis"] = analysisObj;
        }

        try {
            const {display_id, success} = await super.saveDatasetDisplay(datasetId, displayId, label, plotType, plotConfig);
            if (!success) {
                throw new Error("Could not save this new display. Please contact the gEAR team.");
            }

            // Make new display card and make it the default display
            renderUserDisplayCard({id: display_id, label, plot_type: plotType, plotly_config: plotConfig}, display_id);

            return display_id;
        } catch (error) {
            logErrorInConsole(error);
            const msg = "Could not save this new display. Please contact the gEAR team."
            createToast(error);
            throw new Error(msg);
        }
    },

    /**
     * Saves the default display with the specified displayId.
     *
     * @param {string} displayId - The ID of the display to be saved as default.
     * @returns {Promise<void>} - A promise that resolves when the default display is saved successfully.
     * @throws {Error} - If the default display cannot be saved.
     */
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

        const currentDefaultElt = document.getElementById(`${displayId}-default`);
        currentDefaultElt.disabled = true;
        currentDefaultElt.textContent = "Default";
    }

}
Object.setPrototypeOf(curatorApiCallsMixin, apiCallsMixin);


/**
 * Represents a dataset tree.
 *
 * @class
 * @constructor
 * @param {Object} options - The options for the dataset tree.
 * @param {HTMLElement} options.element - The element to render the dataset tree.
 * @param {HTMLElement} options.searchElement - The element for searching the dataset tree.
 * @param {Function} options.selectCallback - The callback function to be called when a dataset is selected.
 */
const datasetTree = new DatasetTree({
    element: document.getElementById("dataset-tree")
    , searchElement: document.getElementById("dataset-query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current-dataset-c").classList.remove("is-hidden");
        document.getElementById("current-dataset").textContent = e.node.title;
        document.getElementById("current-dataset-post").textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        datasetId = newDatasetId;

        // Click to get to next step
        document.getElementById("load-plot-s").click();
        document.getElementById('new-display').classList.add("is-loading");

        // Clear "success/failure" icons
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.classList.add("is-hidden");
        }
        for (const elt of document.getElementsByClassName("js-step-failure")) {
            elt.classList.add("is-hidden");
        }

        // collapse tree
        e.node.tree.expandAll(false);

        document.getElementById("dataset-s-success").classList.remove("is-hidden");

        // displays
        const {userDisplays, ownerDisplays} = await curatorApiCallsMixin.fetchDatasetDisplays(datasetId);
        let defaultDisplayId;
        try {
            defaultDisplayId = await curatorApiCallsMixin.fetchDefaultDisplay(datasetId);
        } catch (error) {
            defaultDisplayId = -1;  // Cannot make any display a default.
        }
        renderDisplayCards(userDisplays, ownerDisplays, defaultDisplayId);
        document.getElementById('new-display').classList.remove("is-loading");
        document.getElementById('new-display').disabled = false;

        // Clear (and update) options within nice-select2 structure.
        // Not providing the object will duplicate the nice-select2 structure
        analysisObj = null;
        analysisSelect = createAnalysisSelectInstance("analysis-select", analysisSelect);
        plotTypeSelect = createPlotTypeSelectInstance("plot-type-select", plotTypeSelect);

        // Call any curator-specific callbacks
        await curatorSpecifcDatasetTreeCallback();
    })
});

/**
 * Updates the analysis select options based on the fetched public and private analyses.
 * @returns {Promise<void>} A promise that resolves when the analysis select options are updated.
 */
const analysisSelectUpdate = async () => {
    try {
        const {publicAnalyses, privateAnalyses} = await curatorApiCallsMixin.fetchAnalyses(datasetId);
        updateAnalysesOptions(privateAnalyses, publicAnalyses);
        document.getElementById("analysis-type-select-c-success").classList.remove("is-hidden");   // Default analysis is good
    } catch (error) {
        // Show failure state things.
        document.getElementById("plot-type-s-failed").classList.remove("is-hidden");
        document.getElementById("analysis-type-select-c-failed").classList.remove("is-hidden");
        document.getElementById('new-display').classList.remove("is-loading"); // Don't give impression display is still loading
    } finally {
        document.getElementById("load-plot-s-success").classList.remove("is-hidden");
    }

}

/**
 * Chooses an analysis based on the selected options and updates the UI accordingly.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the analysis is chosen and the UI is updated.
 */
const chooseAnalysis = async (event) => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const analysisText = (analysisId?.length) ? analysisId : "Primary Analysis";

    // Display current selected analysis
    document.getElementById("current-analysis").textContent = analysisText;
    document.getElementById("current-analysis-post").textContent = analysisText;


    // NOTE: For now, we can just pass analysis id only to tSNE and be fine
    // Any private dataset will belong to our user. Any public datasets can be found by the API "get_analysis" code.
    if (analysisId) {
        analysisObj = {id: analysisId};
    }

    if (analysisSelect.data.length < 5) return; // Have not retrieved analyses from API yet

    if (analysisId) {
        await Promise.all([
            plotTypeSelectUpdate(analysisId)
            , updateDatasetGenes(analysisId)
        ]);

        // Create facet widget
        facetWidget = await createFacetWidget(datasetId, analysisId, {});
    }
}

/* New display has been chosen, so display analysis and plot type options */
/**
 * Updates the display by enabling the analysis and plot type selects, and updating the gene, analysis, and plot type selects in parallel.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the display is updated.
 */
const chooseNewDisplay = async (event) => {
    document.getElementById('new-display').classList.add("is-loading");
    document.getElementById("analysis-select").disabled = false;

    document.getElementById("plot-type-select").disabled = false;

    document.getElementById("current-plot-type-c").classList.add("is-hidden");
    document.getElementById("plot-type-s-success").classList.add("is-hidden");
    document.getElementById("plot-type-select-c-success").classList.add("is-hidden");
    document.getElementById("analysis-type-select-c-failed").classList.add("is-hidden");
    document.getElementById("plot-type-select-c-failed").classList.remove("is-hidden");
    document.getElementById("plot-options-s-success").classList.add("is-hidden");

    // update genes, analysis, and plot type selects in parallel
    await Promise.all([
        updateDatasetGenes(),
        analysisSelectUpdate(),
        plotTypeSelectUpdate()      // NOTE: Believe updating "disabled" properties triggers the plotTypeSelect "change" element

    ]);

    document.getElementById('new-display').classList.remove("is-loading");
    document.getElementById("plot-type-s").click();
}

/**
 * Handles the selection of a plot type.
 *
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves once the plot type is chosen.
 */
const choosePlotType = async (event) => {
    if (!plotTypeSelect.selectedOptions.length) return;   // Do not trigger after setting disable/enable on options

    // Do not display if default opt is chosen
    const plotType = getSelect2Value(plotTypeSelect)
    if (plotType === "nope") {
        document.getElementById("plot-type-select-c-success").classList.add("is-hidden");
        document.getElementById("plot-type-s-success").classList.add("is-hidden");
        return;
    }

    document.getElementById("plot-type-select-c-failed").classList.add("is-hidden");
    document.getElementById("plot-type-select-c-success").classList.remove("is-hidden");

    document.getElementById("plot-type-s-failed").classList.add("is-hidden");
    document.getElementById("plot-type-s-success").classList.remove("is-hidden");

    document.getElementById("plot-options-s-failed").classList.add("is-hidden");
    document.getElementById("plot-options-s-success").classList.add("is-hidden");

    // Display current selected plot type
    document.getElementById("current-plot-type-c").classList.remove("is-hidden");
    document.getElementById("current-plot-type").textContent = plotType;

    // Create facet widget, which will refresh filters
    facetWidget = await createFacetWidget(datasetId, null, {});
    document.getElementById("facet-content").classList.remove("is-hidden");
    document.getElementById("selected-facets").classList.remove("is-hidden");

    // Reset sortable lists
    document.getElementById("order-section").classList.add("is-hidden");
    document.getElementById("order-container").replaceChildren();

    await includePlotParamOptions();
    document.getElementById("gene-s").click();
}


/**
 * Clones a display based on the event and display object.
 * @param {Event} event - The event object.
 * @param {Object} display - The display object to be cloned.
 * @returns {Promise<void>} - A promise that resolves when the display is cloned.
 */
const cloneDisplay = async (event, display) => {

    const cloneElt = event.currentTarget;
    cloneElt.classList.add("is-loading");

    updateDatasetGenes(),

    document.getElementById("analysis-select").disabled = false;
    document.getElementById("plot-type-select").disabled = false;

    const config = display.plotly_config;

    // analyses
    await analysisSelectUpdate();
    if (config.analysis_id) {
        analysisObj = {id: config.analysis_id};
    } else if (config.analysis) {
        analysisObj = config.analysis;
    }
    // select analysis
    if (analysisObj) {
        setSelectBoxByValue("analysis-select", analysisObj.id);
        analysisSelect.update();
        await chooseAnalysis();
        // NOTE: Analysis will not be chosen if the user cannot access it (i.e. owner curation, private analysis)
    }

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

        setSelectBoxByValue("plot-type-select", plotType);
        plotTypeSelect.update();
        await choosePlotType();
        // In this step, a PlotStyle object is instantiated onto "plotStyle", and we will use that
    } catch (error) {
        console.error(error);
        document.getElementById("plot-type-s-failed").classList.remove("is-hidden");
        document.getElementById("plot-type-select-c-failed").classList.remove("is-hidden");
        document.getElementById("plot-type-s-success").classList.add("is-hidden");
        document.getElementById("plot-type-select-c-success").classList.add("is-hidden");

        return;
    } finally {
        cloneElt.classList.remove("is-loading");
    }

    // Choose gene from config
    if (isMultigene) {
        manuallyEnteredGenes = new Set(config.gene_symbols);
        selected_genes = manuallyEnteredGenes;
        const geneSymbolString = config.gene_symbols.join(" ");
        document.getElementById('genes-manually-entered').value = geneSymbolString;
    } else {
        selectedGene = config.gene_symbol;
    }

    try {
        plotStyle.cloneDisplay(config);
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not clone this display. Please contact the gEAR team."
        createToast(msg);
        return;
    }

    // Mark plot params as success
    document.getElementById("plot-options-s-success").classList.remove("is-hidden");

    // Click "submit" button to load plot
    document.getElementById("plot-btn").click();

}

/**
 * Creates an instance of the analysis select.
 * @param {string} idSelector - The ID of the HTML element to bind the select to.
 * @param {object} analysisSelect - Optional. The existing analysis select object to update.
 * @returns {object} - The analysis select instance.
 */
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

/**
 * Creates a color scale select instance.
 * If the colorscaleSelect object exists, it updates it with the revised data and returns it.
 * Otherwise, it creates a new NiceSelect instance with the provided options.
 * @param {string} idSelector - The ID of the HTML element to bind the NiceSelect instance to.
 * @param {object} colorscaleSelect - The existing colorscaleSelect object to update (optional).
 * @returns {object} - The NiceSelect instance.
 */
const createColorscaleSelectInstance = (idSelector, colorscaleSelect=null) => {
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

/**
 * Creates a facet widget for filtering data based on the provided dataset ID, analysis ID, and filters.
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} analysisId - The ID of the analysis.
 * @param {object} filters - The filters to apply to the data.
 * @returns {FacetWidget} The created facet widget.
 */
const createFacetWidget = async (datasetId, analysisId, filters) => {
    document.getElementById("selected-facets-loader").classList.remove("is-hidden")

    const {aggregations, total_count:totalCount} = await curatorApiCallsMixin.fetchAggregations(datasetId, analysisId, filters);
    document.getElementById("num-selected").textContent = totalCount;

    const facetWidget = new FacetWidget({
        aggregations,
        filters,
        onFilterChange: async (filters, seriesName) => {
            if (filters) {
                try {
                    const {aggregations, total_count:totalCount} = await curatorApiCallsMixin.fetchAggregations(datasetId, analysisId, filters);
                    facetWidget.updateAggregations(aggregations);
                    document.getElementById("num-selected").textContent = totalCount;
                } catch (error) {
                    logErrorInConsole(error);
                }
            } else {
                // Save an extra API call
                facetWidget.updateAggregations(facetWidget.aggregations);
            }
            // Sortable lists need to reflect groups filtered out or unfiltered
            updateOrderSortable();

            // Update levels for based on chosen filter groups
            for (const filterGroup of Object.keys(facetWidget.filters)) {
                levels[filterGroup] = facetWidget.filters[filterGroup];
            }

            curatorSpecifcFacetItemSelectCallback(seriesName);
        }
    });
    document.getElementById("selected-facets-loader").classList.add("is-hidden")
    return facetWidget;
}

/**
 * Creates a plot type select instance.
 * @param {string} idSelector - The ID of the selector element.
 * @param {object} plotTypeSelect - The plot type select object (optional).
 * @returns {object} - The plot type select instance.
 */
const createPlotTypeSelectInstance = (idSelector, plotTypeSelect=null) => {
    // If object exists, just update it with the revised data and return
    if (plotTypeSelect) {
        plotTypeSelect.update();
        return plotTypeSelect;
    }

    // Initialize fixed plot types
    return NiceSelect.bind(document.getElementById(idSelector), {
        placeholder: 'Choose how to plot',
        minimumResultsForSearch: -1
    });
}

/**
 * Creates a plot based on the selected plot type and gene(s).
 * @param {Event} event - The event that triggered the plot creation.
 * @returns {Promise<void>} - A promise that resolves when the plot creation is complete.
 */
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
        plotStyle.plotConfig["gene_symbols"] = Array.from(selected_genes);
    } else {
        plotStyle.plotConfig["gene_symbol"] = selectedGene
    }

    await curatorSpecifcCreatePlot(plotType);


    // Stop loader
    for (const plotBtn of plotBtns) {
        plotBtn.classList.remove("is-loading");
    }

    // Hide this view
    document.getElementById("content-c").classList.add("is-hidden");
    // Generate and display "post-plotting" view/container
    document.getElementById("post-plot-content-c").classList.remove("is-hidden");

}

/**
 * Disables or enables a checkbox and its associated label.
 * If the parent element is a .checkbox class, it will also be disabled or enabled.
 * @param {HTMLElement} checkboxElt - The checkbox element.
 * @param {boolean} state - The state to set for the checkbox and its label. True to disable, false to enable.
 */
const disableCheckboxLabel = (checkboxElt, state) => {
    // NOTE: ".disable" attribute only applies to certain elements (https://www.w3schools.com/tags/att_disabled.asp)
    if (checkboxElt.parentElement.classList.contains("checkbox")) {
        if (state) {
            checkboxElt.parentElement.setAttribute("disabled", "");
        } else {
            checkboxElt.parentElement.removeAttribute("disabled");
        }
    }
}

/**
 * Retrieves updates and additions to the plot from the plot_display_config JS object.
 *
 * @param {Object[]} plotConfObj - The plot configuration object.
 * @param {string} plotType - The type of plot.
 * @param {string} category - The category of updates to retrieve.
 * @returns {Object} - The updates and additions to the plot.
 */
const getPlotlyDisplayUpdates = (plotConfObj, plotType, category) => {
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

/**
 * Retrieves the value of a plot configuration element based on its class name.
 *
 * @param {string} className - The class name of the plot configuration element.
 * @returns {boolean|string|undefined} - The value of the plot configuration element, or undefined if not found.
 */
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

/**
 * Retrieves the plot order from the sortable container.
 * @returns {Object} The plot order object.
 */
const getPlotOrderFromSortable = () => {
    const order = {};
    for (const elt of document.getElementById("order-container").children) {
        const series = elt.querySelector("p").textContent;
        const serialized = sortable(`#${CSS.escape(series)}-order-list`, 'serialize')[0].items;
        // Sort by "sortable" index position
        order[series] = serialized.map((val) => val.label);
    }
    return order;
}

/**
 * Retrieves the value from a select2 element.
 * @param {HTMLSelectElement} select - The select2 element.
 * @returns {string} - The value of the selected option.
 */
const getSelect2Value = (select) => {
    return select.selectedOptions[0].data.value;
}

/**
 * Fetches the content of an HTML file from the specified URL.
 *
 * @param {string} url - The URL of the HTML file to fetch.
 * @returns {Promise<string>} - A promise that resolves with the content of the HTML file as a string.
 */
const includeHtml = async (url) => {
    const preResponse = await fetch(url, {cache: "reload"});
    return await preResponse.text();
}

/**
 * Includes plot parameter options and performs necessary setup for plotting.
 * @returns {Promise<void>} A promise that resolves once the plot parameter options are included and setup is complete.
 */
const includePlotParamOptions = async () => {
    const plotType = getSelect2Value(plotTypeSelect);

    // New plot... so disable plot button
    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
        plotBtn.disabled = true;
    }

    plotStyle = curatorSpecificPlotStyle(plotType);
    if (!plotStyle) {
        console.warn(`Plot type ${plotType} not recognized.`)
        document.getElementById("plot-type-s-failed").classList.remove("is-hidden");
        document.getElementById("plot-type-s-success").classList.add("is-hidden");
        return;
    }
    document.getElementById("plot-type-s-failed").classList.add("is-hidden");


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

/**
 * Loads the colorscale select options based on the given parameters.
 *
 * @param {boolean} [isContinuous=false] - Indicates whether the plot uses continuous colorscales.
 * @param {boolean} [isScanpy=false] - Indicates whether the plot is a scanpy plot.
 * @returns {void}
 */
const loadColorscaleSelect = (isContinuous=false, isScanpy=false) => {

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
            // if the plot is a scanpy plot, then the colorscales are in plotly2MatplotlibNames
            if (isScanpy) {
                optionElt.value  = plotly2MatplotlibNames[option.value]
            }
            optionElt.textContent = option.text;
            optgroup.append(optionElt);
        }
        document.getElementById("color-palette-post").append(optgroup);
    }


    // set default color
    let defaultColor = "d3";
    if (isContinuous) {
        defaultColor = "purp";
        if (isMultigene) {
            // I personally don't like purp for multigene plots
            defaultColor = "bluered";
        }
    }
    if (isScanpy) {
        defaultColor = "YlOrRd";
    }

    setSelectBoxByValue("color-palette-post", defaultColor);

    return;

}

/**
 * Loads the dataset tree by fetching dataset information from the curator API.
 * Populates the userDatasets, sharedDatasets, and domainDatasets arrays with dataset information.
 * Generates the dataset tree using the generateTree method of the datasetTree object.
 * @throws {Error} If there is an error fetching the dataset information.
 */
const loadDatasetTree = async () => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const datasetData = await curatorApiCallsMixin.fetchAllDatasets();

        let counter = 0;

        // Create data structure with  dataset information owned by the user
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
        logErrorInConsole(error);
        const msg = "Could not fetch datasets. Please contact the gEAR team."
        createToast(msg);
        document.getElementById("dataset-s-failed").classList.remove("is-hidden");
    }

}

/**
 * Updates the plot type select element based on the available plot types.
 * @param {string|null} analysisId - The ID of the analysis (optional).
 * @returns {Promise<void>} - A promise that resolves when the plot type select is updated.
 */
const plotTypeSelectUpdate = async (analysisId=null) => {
    // NOTE: Believe updating "disabled" properties triggers the plotTypeSelect "change" element
    try {
        // Clear selected options  Need this incase a plot type was selected previously (i.e. cloned display)
        //plotTypeSelect.clear();

        const availablePlotTypes = await curatorApiCallsMixin.fetchAvailablePlotTypes(datasetId, analysisId, isMultigene);
        for (const plotType in availablePlotTypes) {
            const isAllowed = availablePlotTypes[plotType];
            setPlotTypeDisabledState(plotType, isAllowed);
        }

        // set plot type to first option
        setSelectBoxByValue("plot-type-select", "nope");
        plotTypeSelect.update();
    } catch (error) {
        document.getElementById("plot-type-s-failed").classList.remove("is-hidden");
        document.getElementById("plot-type-select-c-failed").classList.remove("is-hidden");
        document.getElementById("plot-type-s-success").classList.add("is-hidden");
        document.getElementById("plot-type-select-c-success").classList.add("is-hidden");;
    }
}

/**
 * Renders display cards for user and owner displays.
 *
 * @param {Array} userDisplays - The array of user displays.
 * @param {Array} ownerDisplays - The array of owner displays.
 * @param {string} defaultDisplayId - The default display ID.
 * @returns {Promise<void>} - A promise that resolves when the display cards are rendered.
 */
const renderDisplayCards = async (userDisplays, ownerDisplays, defaultDisplayId) => {
    const userDisplaysElt = document.getElementById("user-displays");
    const ownerDisplaysElt = document.getElementById("owner-displays");

    // Empty existing displays
    userDisplaysElt.replaceChildren();
    ownerDisplaysElt.replaceChildren();

    // Add titles to each section if there are displays
    if (userDisplays.length) {
        const userTitle = document.createElement("p");
        userTitle.classList.add("has-text-weight-bold", "is-underlined", "column", "is-full");
        userTitle.textContent = "Your Displays";
        userDisplaysElt.append(userTitle);

    }

    if (ownerDisplays.length) {
        const ownerTitle = document.createElement("p");
        ownerTitle.classList.add("has-text-weight-bold", "is-underlined", "column", "is-full");
        ownerTitle.textContent = "Displays by Dataset Owner";
        ownerDisplaysElt.append(ownerTitle);

    }

    for (const display of userDisplays) {
        renderUserDisplayCard(display, defaultDisplayId)
    }

    for (const display of ownerDisplays) {
        renderOwnerDisplayCard(display, defaultDisplayId);
    }
}

/**
 * Renders the order sortable series.
 *
 * @param {string} series - The series to render.
 */
const renderOrderSortableSeries = (series) => {
    const orderContainer = document.getElementById("order-container");

    // If continouous series, cannot sort.
    if (!catColumns.includes(series)) return;

    // Start with a fresh template
    const orderElt = document.getElementById(`${series}-order`);
    if (orderElt) {
        orderElt.remove();
    }

    // Create parent template
    // Designed so the title is a full row and the draggables are 50% width
    const parentList = `<ul id="${series}-order-list" class="content column is-two-thirds js-plot-order-sortable"></ul>`;

    const orderDiv = document.createElement("div");
    orderDiv.id = `${series}-order`;
    orderDiv.classList.add("columns", "is-multiline");
    orderDiv.innerHTML = `
        <p id="${series}-order-title" class="has-text-weight-bold column is-full">${series}</p>
        ${parentList}
    `;
    orderContainer.append(orderDiv);

    // Add in list elements
    for (const group of levels[series]) {
        // If filters are present and group is not in filters, skip
        if (facetWidget.filters.hasOwnProperty(series) && !facetWidget.filters[series].includes(group)) continue;

        const listElt = document.createElement("li");
        listElt.classList.add("has-background-grey-lighter", "has-text-dark");
        listElt.textContent = group;
        document.getElementById(`${series}-order-list`).append(listElt);
    }

    // Create sortable for this series
    sortable(`#${series}-order-list`, {
        hoverClass: "has-text-weight-bold"
        , itemSerializer(item, container) {
            item.label = item.node.textContent
            return item
        },
    });

}

/**
 * Renders the owner display card.
 *
 * @param {Object} display - The display object.
 * @param {string} defaultDisplayId - The ID of the default display.
 * @returns {Promise<void>} - A promise that resolves when the display card is rendered.
 */
const renderOwnerDisplayCard = async (display, defaultDisplayId) => {
    let displayUrl = "";
    try {
        displayUrl = await curatorApiCallsMixin.fetchDatasetDisplayImage(datasetId, display.id);
    } catch (error) {
        displayUrl = "/img/dataset_previews/missing.png";
    }
    const geneSymbol = display.plotly_config.gene_symbol;

    const label = display.label || `Unnamed ${display.plot_type} display`;

    const template = document.getElementById("owner-display-card");
    const displayCard = template.content.cloneNode(true);
    const displayCardElt = displayCard.querySelector(".column");
    displayCardElt.id = `${display.id}-display`;

    const displayCardHeader = displayCard.querySelector(".card-header-title");
    displayCardHeader.textContent = label;

    const displayCardImage = displayCard.querySelector(".card-image img");
    displayCardImage.src = displayUrl;

    const displayCardSubtitle = displayCard.querySelector(".card-content .subtitle");
    if (isMultigene) {
        // Card content should be number of genes
        const numGenes = display.plotly_config.gene_symbols.length;
        displayCardSubtitle.textContent = `Number of genes: ${numGenes}`;
    } else {
        // Card content should be gene symbol
        displayCardSubtitle.textContent = `Gene: ${geneSymbol}`;
    }

    // Edit default properties if this is the default display
    const displayCardDefaultBtn = displayCard.querySelector(".js-display-default");
    displayCardDefaultBtn.id = `${display.id}-default`;
    if (display.id === defaultDisplayId) {
        displayCardDefaultBtn.textContent = "Default";
        displayCardDefaultBtn.disabled = true;
    }
    // Add event listeners
    displayCardDefaultBtn.addEventListener("click", (event) => curatorApiCallsMixin.saveDefaultDisplay(display.id));
    const displayCardCloneBtn = displayCard.querySelector(".js-display-clone");
    displayCardCloneBtn.addEventListener("click", (event) => cloneDisplay(event, display));

    const ownerDisplaysElt = document.getElementById("owner-displays");
    ownerDisplaysElt.append(displayCard);
}

/**
 * Renders a user display card with the given display and default display ID.
 * @param {Object} display - The display object.
 * @param {string} defaultDisplayId - The ID of the default display.
 * @returns {Promise<void>} - A promise that resolves when the display card is rendered.
 */
const renderUserDisplayCard = async (display, defaultDisplayId) => {

    let displayUrl = "";
    try {
        displayUrl = await curatorApiCallsMixin.fetchDatasetDisplayImage(datasetId, display.id);
    } catch (error) {
        displayUrl = "/img/dataset_previews/missing.png";
    }

    const geneSymbol = display.plotly_config.gene_symbol;

    const label = display.label || "Unnamed display"; // Added text to keep "p" tag from collapsing

    const template = document.getElementById("user-display-card");
    const displayCard = template.content.cloneNode(true);
    const displayCardElt = displayCard.querySelector(".column");
    displayCardElt.id = `${display.id}-display`;

    const displayCardHeader = displayCard.querySelector(".card-header-title");
    displayCardHeader.textContent = label;

    const displayCardImage = displayCard.querySelector(".card-image img");
    displayCardImage.src = displayUrl;

    const displayCardSubtitle = displayCard.querySelector(".card-content .subtitle");
    if (isMultigene) {
        // Card content should be number of genes
        const numGenes = display.plotly_config.gene_symbols.length;
        displayCardSubtitle.textContent = `Number of genes: ${numGenes}`;
    } else {
        // Card content should be gene symbol
        displayCardSubtitle.textContent = `Gene: ${geneSymbol}`;
    }

    // Edit default properties if this is the default display
    const displayCardDefaultBtn = displayCard.querySelector(".js-display-default");
    displayCardDefaultBtn.id = `${display.id}-default`;
    if (display.id === defaultDisplayId) {
        displayCardDefaultBtn.textContent = "Default";
        displayCardDefaultBtn.disabled = true;
    }
    // Add event listeners
    displayCardDefaultBtn.addEventListener("click", (event) => curatorApiCallsMixin.saveDefaultDisplay(display.id));
    const displayCardCloneBtn = displayCard.querySelector(".js-display-clone");
    displayCardCloneBtn.addEventListener("click", (event) => cloneDisplay(event, display));
    const displayCardDeleteBtn = displayCard.querySelector(".js-display-delete");
    displayCardDeleteBtn.addEventListener("click", (event) => curatorApiCallsMixin.deleteDisplay(display.id));

    const userDisplaysElt = document.getElementById("user-displays");
    userDisplaysElt.append(displayCard);
}

/**
 * Sets the value of plot elements based on the provided configuration value.
 * @param {string} classSelector - The class selector for the plot elements.
 * @param {boolean|string|number} confVal - The configuration value to set.
 */
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
 * Sets the disabled state of a plot type option.
 * @param {string} plotType - The plot type.
 * @param {boolean} isAllowed - Whether the plot type is allowed or not.
 */
const setPlotTypeDisabledState = (plotType, isAllowed) => {
    if (plotType === "tsne/umap_dynamic") {
        document.getElementById("tsne-dyna-opt").disabled = !isAllowed;
    } else {
        // replace _ with - for id
        const fixedPlotType = plotType.replace("_", "-");
        document.getElementById(`${fixedPlotType}-opt`).disabled = !isAllowed;
    }
}

/**
 * Sets the selected option in a select box by its value. Assumes single selection
 * @param {string} eid - The ID of the select box element.
 * @param {string} val - The value of the option to be selected.
 */
const setSelectBoxByValue = (eid, val) => {
    // Modified to set value and "selected" so nice-select2 extractData() will catch it
    // Taken from https://stackoverflow.com/a/20662180
    const elt = document.getElementById(eid);

    // Clear selected attribute from the selected value (if multiple, only the first value)
    elt.options[elt.selectedIndex].removeAttribute("selected");

    for (const i in elt.options) {
        if (elt.options[i].value === val) {
            // By using "selected" instead of "value", we can account for the "multiple" property
            elt.options[i].setAttribute("selected", true);
            return;
        }
    }
}

/**
 * Sets up an event listener on elements with the specified class selector.
 * When the value of any element changes, it copies the new value to all other elements with the same class selector.
 * Additionally, it updates the disabled state and checked state of the elements based on the triggering element.
 * @param {string} classSelector - The class selector for the elements to attach the event listener to.
 */
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

/**
 * Sets up validation events for elements with the class "js-plot-req".
 * @returns {void}
 */
const setupValidationEvents = () => {
    const validationElts = document.getElementsByClassName("js-plot-req");
    for (const elt of validationElts ) {
        elt.addEventListener("change", validateRequirements);
    }
}

/**
 * Updates the options for the private and public analyses select elements.
 *
 * @param {Array} privateAnalyses - The array of private analyses.
 * @param {Array} publicAnalyses - The array of public analyses.
 */
const updateAnalysesOptions = (privateAnalyses, publicAnalyses) => {
    const privateAnalysesElt = document.getElementById("private-analyses");
    const publicAnalysesElt = document.getElementById("public-analyses");

    // Empty the old optgroups
    privateAnalysesElt.replaceChildren();
    publicAnalysesElt.replaceChildren();

    const analysisElt = document.getElementById("analysis-select");
    analysisElt.parentElement.classList.add("is-loading");

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
        privateAnalysesElt.append(option);
    }
    for (const analysis of publicAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.value = analysis.id;
        publicAnalysesElt.append(option);
    }

    // Update select2
    analysisSelect.update();

    analysisElt.parentElement.classList.remove("is-loading");

}

/**
 * Updates the gene select element with gene symbols.
 * @param {string|null} analysisId - The analysis ID (optional).
 * @returns {Promise<void>} - A promise that resolves when the gene select element is updated.
 */
const updateDatasetGenes = async (analysisId=null) => {
    try {
        const geneSymbols = await curatorApiCallsMixin.fetchGeneSymbols(datasetId, analysisId);
        curatorSpecificUpdateDatasetGenes(geneSymbols);
    } catch (error) {
        document.getElementById("gene-s-failed").classList.remove("is-hidden");
    }
}

/**
 * Updates the sortable order of plot param series based on the current selection.
 */
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
        // 3. Series is in sortableSet but not seriesSet, remove <series>-order element
        if (!seriesSet.has(series)) {
            const orderElt = document.getElementById(`${series}-order`);
            orderElt.remove();
        }
    }


    const orderContainer = document.getElementById("order-container");
    const orderSection = document.getElementById("order-section");

    // Pre-emptively hide the container but show it ass
    if (!orderContainer.children.length) {
        orderSection.classList.add("is-hidden");
        return;
    }

    orderSection.classList.remove("is-hidden");

}

/**
 * Validates the requirements for plotting and updates the UI accordingly.
 *
 * @param {Event} event - The event object triggered by the input element.
 * @returns {void}
 */
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
            document.getElementById("plot-options-s-success").classList.remove("is-hidden");
            document.getElementById("plot-options-s-failed").classList.add("is-hidden");
        }

    // Perform misc. validation checks
    if (!curatorSpecificValidationChecks()){
        for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
            plotBtn.disabled = true;
        }
    }

        return;
    }

    // Required paramater has no value. Indicate it and disable plot buttons
    elt.parentElement.classList.add("is-danger");
    elt.parentElement.classList.remove("is-success");

    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
        plotBtn.disabled = true;
    }
    document.getElementById("plot-options-s-success").classList.add("is-hidden");
    document.getElementById("plot-options-s-failed").classList.remove("is-hidden");
}

document.getElementById("new-display").addEventListener("click", chooseNewDisplay);
document.getElementById("analysis-select").addEventListener("change", chooseAnalysis);
document.getElementById("plot-type-select").addEventListener("change", choosePlotType);

const plotBtns = document.getElementsByClassName("js-plot-btn");
for (const plotBtn of plotBtns) {
    plotBtn.addEventListener("click", createPlot);
}

document.getElementById("save-json-config").addEventListener("click", () => {
    // Config plot configuration to JSON for sharing (or passing to API by hand)
    const blob = new Blob([JSON.stringify({...plotStyle.plotConfig, plot_type:plotStyle.plotType})]);
    const link = document.createElement("a");
    link.download = "gEAR_plot_configuration.json";
    link.href = window.URL.createObjectURL(blob);
    link.click()
    // Give confirmation
    document.getElementById("save-json-config").classList.add("is-success");
    setTimeout(() => {
        document.getElementById("save-json-config").classList.remove("is-success");
    }, 1000);
});

document.getElementById("save-display-btn").addEventListener("click", async (event) => {
    // Save new plot display.
    const label = document.getElementById("new-display-label").value;
    event.target.classList.add("is-loading");
    try {
        const displayId = await curatorApiCallsMixin.saveDatasetDisplay(datasetId, null, label, plotStyle.plotType, plotStyle.plotConfig);
        createToast("Display saved.", "is-success");

        if (document.getElementById("make-default-display-check").checked) {
            apiCallsMixin.saveDefaultDisplay(displayId);
        }
    } catch (error) {
        //pass - handled in functions
    } finally {
        event.target.classList.remove("is-loading");
    }
});

document.getElementById("edit-params").addEventListener("click", (event) => {
    event.target.classList.add("is-loading");
    // Hide this view
    document.getElementById("content-c").classList.remove("is-hidden");
    // Generate and display "post-plotting" view/container
    document.getElementById("post-plot-content-c").classList.add("is-hidden");

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
/**
 * Handles page-specific login UI updates.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are completed.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {

    curatorSpecificNavbarUpdates();

    const sessionId = CURRENT_USER.session_id;
    if (! sessionId ) {
        createToast("Not logged in so saving displays is disabled.", "is-warning");
        document.getElementById("save-display-btn").disabled = true;
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
                datasetTree.tree.setActiveNode(foundNode);
                datasetTree.selectCallback({node: foundNode});  // manually trigger the "activate" event.
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
