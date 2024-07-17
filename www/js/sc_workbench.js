"use strict";

let currentAnalysis;
let clickedMarkerGenes = new Set();
let typedMarkerGenes = new Set();
let currentLabel = null;
let datasetId = null;

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
        //document.querySelector(UI.currentDatasetContainer).classList.remove("is-hidden");
        document.querySelector(UI.currentDatasetElt).textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        const organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        createToast("Loading dataset", "is-info");

        datasetId = newDatasetId;

        // Clear "success/failure" icons
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.classList.add("is-hidden");
        }
        for (const elt of document.getElementsByClassName("js-step-failure")) {
            elt.classList.add("is-hidden");
        }

        // Hide Steppers
        for (const elt of document.querySelectorAll(UI.stepsElts)) {
            elt.classList.add("is-hidden");
        }

        // collapse tree
        e.node.tree.expandAll(false);

        // We only want to reset these plots only if the dataset changes.
        document.querySelector(UI.primaryInitialScatterContainer).replaceChildren();
        document.querySelector(UI.primaryInitialViolinContainer).replaceChildren();

        if (currentAnalysis) {
            resetWorkbench();
        }

        // Reset to the default state, so that users are forced to choose an analysis.
        document.querySelector(UI.analysisSelect).querySelector("option[data-analysis-id='-1']").selected = true;
        document.querySelector(UI.analysisWorkflowElt).classList.add("is-hidden");
        document.querySelector(UI.btnProgressGuideElt).classList.add("is-hidden");

        document.querySelector(UI.currentAnalysisElt).textContent = "None selected";

        // This is a placeholder to retrieve preliminary figures which are stored in the "primary" directory
        currentAnalysis = new Analysis({id: datasetId, type: "primary", datasetIsRaw: true});

        // Technically these could load asynchronously, but logically the progress logs make more sense sequentially
        try {
            await getDatasetInfo(datasetId);
            await currentAnalysis.loadPreliminaryFigures();
        } catch (error) {
            logErrorInConsole(error);
            // pass
        }

    })
});

    /**
     * Counts and highlights duplicate values in the clustering group labels step.
     * @returns {number} The number of duplicate values found.
     */
    const countAndHighlightDuplicates = () => {
        const allValues = [];
        const dupValues = [];

        // first remove any duplicate-labeled ones
        for (const elt of document.querySelectorAll(UI.clusterGroupLabelsInputElts)) {
            elt.classList.remove('has-background-danger');

            const clusterLabel = elt.value.trim();

            // this means it WAS found
            if (allValues.includes(clusterLabel)) {
                dupValues.push(clusterLabel);
                elt.classList.add('has-background-danger');
            } else {
                allValues.push(clusterLabel);
            }
        }

        return dupValues.length;
    }

/**
 * Downloads the table data as an Excel file.
 *
 * @param {string} tableId - The ID of the table element.
 * @param {string} filename - The name of the downloaded file.
 */
const downloadTableAsExcel = (tableId, filename) => {
    let tableStr = '';

    // Loop through the table header and add each column to the table string
    for (const elt of document.querySelectorAll(`${tableId} thead tr th`)){
        tableStr += `${elt.textContent}\t`;
    }
    tableStr = `${tableStr.trim()}\n`;

    // Loop through the table body and add each row to the table string
    for (const row of document.querySelectorAll(`${tableId} tbody tr`)){
        for (const cell of row.querySelectorAll('td')){
            tableStr += `${cell.textContent}\t`;
        }
        tableStr = `${tableStr.trim()}\n`;
    }

    // Create a link element, set the href to the table string, and download it
    const element = document.createElement('a');
    element.setAttribute('href', `data:text/plain;charset=utf-8,${encodeURIComponent(tableStr)}`);
    element.setAttribute('download', filename);
    element.classList.add('is-hidden');
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

/**
 * Retrieves dataset information and performs necessary UI updates.
 *
 * @param {string} datasetId - The ID of the dataset to retrieve information for.
 * @returns {Promise<void>} - A promise that resolves when the dataset information is retrieved and UI updates are complete.
 */
const getDatasetInfo = async (datasetId) => {
    document.querySelector(UI.analysisSelect).disabled = true;

    try {
        const data = await apiCallsMixin.fetchDatasetInfo(datasetId);

        const ds = new Dataset(data);

        currentAnalysis.dataset = ds;

        document.querySelector(UI.primaryFilterSection).classList.remove("is-hidden");
        analysisLabels = currentAnalysis.getSavedAnalysesList(ds.id, -1, 'sc_workbench');   // select first "selct an analysis" option

        document.querySelector(UI.primaryInitialInfoSection).classList.remove("is-hidden");
        document.querySelector(UI.selectedDatasetShapeInitialElt).textContent = currentAnalysis.dataset.shape();

        document.querySelector(UI.analysisSelect).disabled = false;
        createToast("Dataset loaded", "is-success");
    } catch (error) {
        createToast("Failed to access dataset");
        logErrorInConsole(`Failed ID was: ${datasetId} because msg: ${error.message}`);
        document.querySelector(UI.analysisSelect).disabled = true;

    }
}

const getEmbeddedTsneDisplay = async (datasetId) => {
    const {data} = await axios.post("./cgi/get_embedded_tsne_display.cgi", convertToFormData({ dataset_id: datasetId }));
    return data;
}

/**
 * Retrieves an array of genes from a collection of cells.
 *
 * @param {NodeList} cells - The collection of cells.
 * @returns {string[]} An array of genes extracted from the cells.
 */
const getGenesFromCells = (cells) => {
    return [...cells].map(cell => cell.textContent.trim());
}

/**
 * Retrieves the t-SNE image data for a given gene symbol and configuration.
 *
 * @param {string} geneSymbol - The gene symbol to retrieve t-SNE image data for.
 * @param {object} config - The configuration object.
 * @returns {Promise<string>} - The t-SNE image data.
 */
const getTsneImageData = async (geneSymbol, config) => {
    config.colorblind_mode = CURRENT_USER.colorblind_mode;
    config.gene_symbol = geneSymbol;

    // in order to avoid circular references (since analysis is referenced in the individual step objects),
    //  we need to create a smaller analysis object to pass to the API

    const analysis = {
        "id": currentAnalysis.id,
        "type": currentAnalysis.type,
    }

    const data = await apiCallsMixin.fetchTsneImage(currentAnalysis.dataset.id, analysis, "tsne_static", config);

    if (!data.success || data.success < 1) {
        const message = data.message || "Unknown error";
        throw new Error(message);
    }
    return data.image;
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
        const datasetData = await apiCallsMixin.fetchAllDatasets();

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
 * Resets the manual marker gene entries.
 *
 * @returns {void}
 */
const resetManualMarkerGeneEntries = () => {
    clickedMarkerGenes = new Set();

    // remember which GOI are from the table
    for (const elt of document.querySelectorAll(UI.markerGenesTableHighlightedElts)) {
        clickedMarkerGenes.add(elt.textContent);
    }
}

const resetWorkbench = () => {
    // Performs all the steps needed to reset the workbench if the input dataset changes

    // handle any data structures
    // ? Is this necessary - new instance should clear this. Maybe we only do the UI reset?
    currentAnalysis.reset();

    // re-enable all buttons in the blocked step
    document.querySelectorAll(".js-step-collapsable button").forEach((button) => {
        button.disabled = false;
    });

    /*
    for (const elt of document.querySelectorAll('.reset-on-change')) {
        // TODO - replace
        elt.classList.add("is-hidden");
    }

    for (const elt of document.querySelectorAll('.empty-on-change')) {
        // TODO - replace
        elt.replaceChildren();
    }
    */
    document.querySelector(UI.newAnalysisLabelElt).value = '';
}

/**
 * Saves the marker gene list.
 * @returns {void}
 */
const saveMarkerGeneList = () => {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: CURRENT_USER.session_id,
        label: document.querySelector(UI.markerGenesListNameElt).value,
        gctype: 'unweighted-list',
        organism_id: currentAnalysis.dataset.organism_id,
        is_public: 0
    });

    for (const geneId of currentAnalysis.markerGenes.genesOfInterest) {
        const gene = new Gene({
            //id: geneId,    // TODO: figure out how to get ensembl ID for this
            gene_symbol: geneId,
        });
        gc.addGene(gene);
    };

    gc.save(updateUiAfterMarkerGeneListSaveSuccess, updateUiAfterMarkerGeneListSaveFailure);
}

/**
 * Saves the PCA gene list.
 *
 * @async
 * @function savePcaGeneList
 * @returns {Promise<void>} A promise that resolves when the PCA gene list is saved successfully.
 * @throws {Error} If there is an error saving the PCA gene list.
 */
const savePcaGeneList = async () => {
    // ? Move to PCA class?

    try {
        const {data} = await axios.post("./cgi/get_PCs_from_anndata.cgi", convertToFormData({
            'dataset_id': currentAnalysis.dataset.id,
            'analysis_id': currentAnalysis.id,
            'analysis_type': currentAnalysis.type,
            'session_id': currentAnalysis.userSessionId,
        }));

        if (!data.success || data.success < 1) {
            const message = data.msg || "Unknown error";
            throw new Error(message);
        }

        const weightLabels = data.pc_data.columns;

        const geneList = new WeightedGeneCart({
                session_id: CURRENT_USER.session_id,
                label: document.querySelector(UI.pcaGeneListNameElt).value,
                gctype: 'weighted-list',
                organism_id: currentAnalysis.dataset.organism_id,
                is_public: 0
            }, weightLabels
        );

        data.pc_data.index.forEach((geneId, i) => {
            const weights = data.pc_data.data[i];
            const gene = new WeightedGene({
                id: geneId,
                gene_symbol: data.gene_symbols[i]
            }, weights
            );
            geneList.addGene(gene);
        });

        geneList.save(updateUiAfterPcaGeneListSaveSuccess, updateUiAfterPcaGeneListSaveFailure);

    } catch (error) {
        createToast(`Error saving PCs as gene list: ${error.message}`);
        document.querySelector(UI.btnSavePcaGeneListElt).disabled = false;
    }
}

/**
 * Updates the highlight status of an element or a collection of elements.
 * @param {HTMLElement|HTMLElement[]} element - The element or collection of elements to update.
 * @param {string[]} genes - The genes associated with the element(s).
 * @param {boolean} addHighlight - Indicates whether to add or remove the highlight.
 */
const updateHighlightStatus = (element, genes, addHighlight) => {
    const classlistMethod = addHighlight ? 'add' : 'remove';

    // If element is an array, it is a collection of elements
    if (Array.isArray(element)) {
        for (const elt of element) {
            elt.classList[classlistMethod]('js-highlighted', "has-background-info", "has-text-white");
        }
    } else {
        element.classList[classlistMethod]('js-highlighted', "has-background-info", "has-text-white");
    }
    for (const gene of genes) {
        const setMethod = addHighlight ? 'add' : 'delete';
        clickedMarkerGenes[setMethod](gene);
    };
    currentAnalysis.markerGenes.addRemoveGenesOfInterest(genes, addHighlight);
}

/**
 * Updates the manual marker gene entries based on the provided gene string.
 *
 * @param {string} geneString - The gene string to update the marker gene entries with.
 */
const updateManualMarkerGeneEntries = (geneString) => {
    typedMarkerGenes = new Set();
    const geneSyms = geneString.split(',');

    for (let geneSym of geneSyms) {
        geneSym = geneSym.trim();
        if (geneSym) {
            typedMarkerGenes.add(geneSym);
        }
    }

    const counterSet = new Set([...typedMarkerGenes, ...clickedMarkerGenes]);
    document.querySelector(UI.markerGenesUniqueCountElt).textContent = counterSet.size;
    document.querySelector(UI.markerGenesEnteredCountElt).textContent = typedMarkerGenes.size;
}

/**
 * Updates the UI after a marker gene list is successfully saved.
 *
 * @param {Object} geneCart - The saved marker gene list object.
 */
const updateUiAfterMarkerGeneListSaveSuccess = (geneCart) => {
    createToast("Saved marker gene list", "is-success");
}

/**
 * Updates the UI after a failure to save the marker gene list.
 *
 * @param {Object} geneCart - The marker gene list object.
 * @param {string} message - The error message.
 */
const updateUiAfterMarkerGeneListSaveFailure = (geneCart, message) => {
    createToast(`Error saving gene list: ${geneCart.label}`);
    logErrorInConsole(message);
    document.querySelector(UI.btnSaveMarkerGeneListElt).disabled = false;
}

/**
 * Updates the UI after successfully saving the PCA gene list.
 *
 * @param {Object} geneCart - The saved gene list object.
 */
const updateUiAfterPcaGeneListSaveSuccess = (geneCart) => {
    createToast("Saved weighted gene list", "is-success");
}

/**
 * Updates the UI after a failure to save the weighted gene list.
 *
 * @param {Object} geneCart - The gene list object.
 * @param {string} message - The error message.
 */
const updateUiAfterPcaGeneListSaveFailure = (geneCart, message) => {
    createToast(`Error saving gene list: ${geneCart.label}`);
    logErrorInConsole(message);
    document.querySelector(UI.btnSavePcaGeneListElt).disabled = false;
}

/**
 * Validates the marker gene selection and updates the UI accordingly.
 *
 * @returns {boolean} - Returns true if there are marker genes selected, false otherwise.
 */
const validateMarkerGeneSelection = () => {
    currentAnalysis.markerGenes.genesOfInterest = new Set([...typedMarkerGenes, ...clickedMarkerGenes]);
    // Only allow saving of gene list if genes are selected
    document.querySelector(UI.btnSaveMarkerGeneListElt).disabled = !currentAnalysis.markerGenes.genesOfInterest.size;
    return Boolean(currentAnalysis.markerGenes.genesOfInterest.size)
}

/* -- page entrypoint -- */

/**
 * Handles page-specific login UI updates (after the login event is triggered).
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are complete.
 */

const handlePageSpecificLoginUIUpdates = async (event) => {
	document.getElementById("page-header-label").textContent = "Single Cell Workbench";

    const sessionId = CURRENT_USER.session_id;
    if (! sessionId ) {
        createToast("Not logged in so saving analyses is disabled.", "is-warning");
        document.querySelector(UI.btnSaveAnalysisElt).disabled = true;

        // TODO: Other actions
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

        // ? This could be used to pre-select an analysis
        //currentAnalysis = new Analysis();

	} catch (error) {
		logErrorInConsole(error);
	}

}

/* Event listeners for elements already loaded */

// General

// Show/Hide #progress-guide
document.querySelector(UI.btnProgressGuideElt).addEventListener("click", (event) => {
    document.querySelector(UI.progressGuideElt).classList.toggle("is-hidden");
});

// if scrolling makes #summary-s go above top of screen, make it sticky
window.addEventListener("scroll", (event) => {
    const summary = document.querySelector(UI.summarySection);
    summary.classList.remove("stick-to-top");
    if (summary.getBoundingClientRect().top < 0) {
        summary.classList.add("stick-to-top");
    }
});

// Handle the "click" event for the steps
for (const step of document.querySelectorAll(UI.stepSegmentElts)) {
    step.addEventListener("click", (event) => {
        event.preventDefault(); // stop normal click behavior
        const summary = document.querySelector(UI.summarySection);

        // if the summary is not already sticky, make it so
        // Then scroll to under the sticky header
        if (!summary.classList.contains("stick-to-top")) {
            summary.classList.add("stick-to-top");
        }

        // get height of anchor href and scroll to it (minus the height of the sticky summary section)
        const href = step.querySelector(".steps-marker").getAttribute("href");
        const hrefHeight = document.querySelector(href).offsetTop - summary.offsetHeight;
        window.scrollTo(0, hrefHeight);

    });
}

// TODO: Fix this to scroll to the correct offset height
for (const step of document.querySelectorAll(".js-step h5")) {
    step.addEventListener("click", (event) => {
        event.preventDefault(); // stop normal click behavior

        const summary = document.querySelector(UI.summarySection);

        // get height of this step container and scroll to it (minus the height of the sticky summary section)
        const stepHeight = event.target.parentNode.offsetTop - summary.offsetHeight;
        window.scrollTo(0, stepHeight);

    });
}

document.querySelector(UI.btnDeleteSavedAnalysisElt).addEventListener("click", async (event) => {
    // Delete the current analysis
    await currentAnalysis.delete();
});

document.querySelector(UI.btnMakePublicCopyElt).addEventListener("click", async (event) => {
    // Make a public copy of the current analysis
    await currentAnalysis.makePublicCopy();
});

document.querySelector(UI.btnSaveAnalysisElt).addEventListener("click", async (event) => {
    // Save the current analysis to the user area
    event.target.classList.add("is-loading");
    await currentAnalysis.saveToUserArea();
    event.target.classList.remove("is-loading");

});

for (const button of document.querySelectorAll(UI.analysisDeleteElts)) {
    button.addEventListener("click", (event) => {
        // Delete the current analysis
        currentAnalysis.delete();
    })
}

// Show the "rename" analysis label input when the button is clicked
for (const button of document.querySelectorAll(UI.analysisRenameElts)) {
    button.addEventListener("click", (event) => {
        document.querySelector(UI.newAnalysisLabelElt).value = currentAnalysis.label;
        document.querySelector(UI.newAnalysisLabelContainer).classList.remove("is-hidden");
        document.querySelector(UI.currentAnalysisElt).textContent = currentAnalysis.label;
    });
}

document.querySelector(UI.btnNewAnalysisLabelSaveElt).addEventListener("click", async (event) => {
    // Save the new label to the current analysis
    currentAnalysis.label = document.querySelector(UI.newAnalysisLabelElt).value;
    await currentAnalysis.save();
    document.querySelector(UI.newAnalysisLabelContainer).classList.add("is-hidden");
});

document.querySelector(UI.btnNewAnalysisLabelCancelElt).addEventListener("click", async (event) => {
    // Reset the label to the current analysis label
    document.querySelector(UI.newAnalysisLabelElt).value = currentAnalysis.label;
    document.querySelector(UI.newAnalysisLabelContainer).classList.add("is-hidden");
});

// Handle the analysis selection change
document.querySelector(UI.analysisSelect).addEventListener("change", async (event) => {

    // Grab the dataset ID from the current analysis to reuse it
    const datasetObj = currentAnalysis.dataset;


    document.querySelector(UI.currentAnalysisElt).textContent = event.target.selectedOptions[0].textContent;

    document.querySelector(UI.deNovoStepsElt).classList.add("is-hidden");
    document.querySelector(UI.primaryStepsElt).classList.add("is-hidden");
    document.querySelector(UI.btnProgressGuideElt).classList.add("is-hidden");
    document.querySelector(UI.progressGuideElt).classList.add("is-hidden");

    // Analysis ID -1 is "select an analysis"
    if (event.target.value === "-1") {
        document.querySelector(UI.analysisWorkflowElt).classList.add("is-hidden");
        return;
    }

    document.querySelector(UI.analysisWorkflowElt).classList.remove("is-hidden");

    // Analysis ID 0 is the blank 'New' one
    if (event.target.value === "0") {
        resetWorkbench();

        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.add("is-hidden");

        currentAnalysis = new Analysis({
            "datasetObj": datasetObj,
            'type': 'primary',  // This is to start with the original dataset, but will eventually change to 'user_unsaved'
            'datasetIsRaw': true}
        );

        document.querySelector(UI.deNovoStepsElt).classList.remove("is-hidden");
        document.querySelector(UI.btnProgressGuideElt).classList.remove("is-hidden");
        // Reset the stepper
        resetStepperWithHrefs("#primary-filter-s");

        // Jump to the primary filter step
        document.querySelector(UI.primaryFilterSection).click();
        document.querySelector(`a[href='${UI.primaryFilterSection}']`).click();
        return;
    }
    createToast("Loading stored analysis", "is-info");

    document.querySelector(UI.newAnalysisLabelContainer).classList.add("is-hidden");

    resetWorkbench();


    const selectedOption = event.target.selectedOptions[0];
    currentAnalysis.type = selectedOption.dataset.analysisType;
    currentAnalysis.id = selectedOption.dataset.analysisId;

    await currentAnalysis.getStoredAnalysis();    // await-able

    if (currentAnalysis.type === 'primary') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.remove("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
    }

    if (currentAnalysis.type === 'user_saved') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
        createToast("This analysis is stored in your profile.", "is-info", true);
        document.querySelector(UI.btnMakePublicCopyElt).classList.remove("is-hidden");
    }

    if (currentAnalysis.type === 'user_unsaved') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
    }

    if (currentAnalysis.type === 'public') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        createToast("Changes made to this public analysis will spawn a local copy within your profile.", "is-info", true);
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
    }

});

document.querySelector(UI.newAnalysisLabelElt).addEventListener("keyup", (event) => {
    // Update the new analysis label if it is not a duplicate
    if (analysisLabels.has(event.target.value.trim())) {
        if (event.target.value.trim() !== currentLabel) {
            event.target.classList.add("duplicate");
            document.querySelector(UI.btnNewAnalysisLabelSaveElt).disabled = true;
            document.querySelector(UI.duplicateLabelWarningElt).classList.remove("is-hidden");
        }
        return;
    }

    document.querySelector(UI.duplicateLabelWarningElt).classList.add("is-hidden");
    event.target.classList.remove("duplicate");
    document.querySelector(UI.btnNewAnalysisLabelSaveElt).disabled = false;
});

document.querySelector(UI.newAnalysisLabelElt).addEventListener("focus", (event) => {
    // Reset the new analysis label
    currentLabel = event.target.value.trim();
});

// Labeled tSNE

document.querySelector(UI.btnLabeledTsneRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    document.querySelector(UI.labeledTsnePlotContainer).replaceChildren();
    document.querySelector(UI.labeledTsnePlotContainer).classList.remove("is-hidden");

    createToast("Generating labeled tSNE plot", "is-info");
    await currentAnalysis.labeledTsne.runAnalysis();
    event.target.classList.remove("is-loading");

});

// Primary Filter

document.querySelector(UI.btnApplyPrimaryFilterElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Apply the primary filter to the dataset
    await currentAnalysis.primaryFilter.applyPrimaryFilter();
    event.target.classList.remove("is-loading");
});

// QC by Mito
document.querySelector(UI.btnDoAnalysisQcByMitoElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the QC by Mito analysis
    await currentAnalysis.qcByMito.runAnalysis(0);
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.btnQbmSaveElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    disableAndHideElement(document.querySelector(UI.btnDoAnalysisQcByMitoElt));
    // Save the QC by Mito analysis
    await currentAnalysis.qcByMito.runAnalysis(1);
    event.target.classList.remove("is-loading");
});

// Select Variable Genes

document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the select variable genes analysis
    await currentAnalysis.selectVariableGenes.runAnalysis(0);
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.btnAsvgSaveElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    disableAndHideElement(document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt));

    // Save the select variable genes analysis
    await currentAnalysis.selectVariableGenes.runAnalysis(1);
    event.target.classList.remove("is-loading");
});

// PCA

document.querySelector(UI.btnPcaRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the PCA analysis
    await currentAnalysis.checkDependenciesAndRun(currentAnalysis.pca.runAnalysis.bind(currentAnalysis.pca));
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.btnPcaTopGenesElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the PCA top genes
    await currentAnalysis.checkDependenciesAndRun(currentAnalysis.pca.runAnalysisTopGenes.bind(currentAnalysis.pca));
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.pcaGeneListNameElt).addEventListener("input", (event) => {
    // enable the save button if a name is entered
    document.querySelector(UI.btnSavePcaGeneListElt).disabled = true;
    if (event.target.value) {
        document.querySelector(UI.btnSavePcaGeneListElt).disabled = false;
    }
});

document.querySelector(UI.btnSavePcaGeneListElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    if (CURRENT_USER) {
        await savePcaGeneList();
    } else {
        createToast("You must be signed in to save a PCA gene list.");
    }
    event.target.classList.remove("is-loading");
});

// tSNE

document.querySelector(UI.btnTsneRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the tSNE analysis (and/or UMAP)
    await currentAnalysis.checkDependenciesAndRun(currentAnalysis.tsne.runAnalysis.bind(currentAnalysis.tsne));
    event.target.classList.remove("is-loading");
});

// Clustering

document.querySelector(UI.btnClusteringRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the clustering analysis
    await currentAnalysis.checkDependenciesAndRun(currentAnalysis.clustering.runAnalysis.bind(currentAnalysis.clustering));
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.btnClusteringEditRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");

    // Check for duplicate labels and run the clustering edit analysis
    //  if there are none (or if the user has selected to merge clusters)

    // Remove any duplicate labels
    for (const elt of document.querySelectorAll(UI.clusterGroupLabelsInputElts)) {
        elt.classList.remove("has-background-danger");
    }

    // If the user has selected to merge clusters, run the clustering edit analysis
    if (document.querySelector(UI.clusteringMergeClustersElt).checked) {
        await currentAnalysis.checkDependenciesAndRun(currentAnalysis.clusteringEdit.runAnalysis.bind(currentAnalysis.clusteringEdit));
        event.target.classList.remove("is-loading");
        return;
    }

    // Otherwise, check for duplicates and run the clustering analysis if there are none
    const duplicateCount = countAndHighlightDuplicates();
    if (duplicateCount === 0) {
        await currentAnalysis.checkDependenciesAndRun(currentAnalysis.clusteringEdit.runAnalysis.bind(currentAnalysis.clusteringEdit));
    }

    event.target.classList.remove("is-loading");
});

// Marker Genes

document.querySelector(UI.btnMarkerGenesRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the marker genes analysis
    await currentAnalysis.checkDependenciesAndRun(currentAnalysis.markerGenes.runAnalysis.bind(currentAnalysis.markerGenes));
    event.target.classList.remove("is-loading");

});

document.querySelector(UI.btnVisualizeMarkerGenesElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Visualize the marker genes
    await currentAnalysis.markerGenes.performMarkerGeneVisualization();
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.markerGenesManuallyEnteredElt).addEventListener("keyup", (event) => {
    // Update the manual marker gene entries based on the entered gene string
    const geneString = event.target.value;
    updateManualMarkerGeneEntries(geneString);
});

document.querySelector(UI.btnDownloadMarkerGenesElt).addEventListener("click", (event) => {
    // Download the marker genes table
    currentAnalysis.markerGenes.downloadMarkerGenesTable();
});

document.querySelector(UI.markerGenesTableElt).addEventListener("click", (event) => {
    // Handle the marker genes table cell clicks

    const clickedCell = event.target.classList.contains("js-col-idx") || event.target.classList.contains("js-row-idx") ? event.target.closest("th") : event.target.closest("td");

    const isHighlighted = clickedCell.classList.contains("js-highlighted");
    let genesOfInterest;

    switch (true) {
        // Check if the clicked cell is a row index cell
        case clickedCell.classList.contains("js-row-idx"):
            const rowCells = [...clickedCell.parentNode.children].filter(child => child !== clickedCell);
            genesOfInterest = getGenesFromCells(rowCells);

            // highlight all cells in the row
            updateHighlightStatus(rowCells, genesOfInterest, !isHighlighted);

            // highlight clicked cell so toggle works
            updateHighlightStatus(clickedCell, [], !isHighlighted);

            break;

        // Check if the clicked cell is a column index cell
        case clickedCell.classList.contains("js-col-idx"):
            // if this cell is the first one, it overlaps with the row indexes and should be ignored
            if (clickedCell.cellIndex === 0) {
                break;
            }
            const clickedHeaderIndex = clickedCell.cellIndex + 1;
            const colCells = document.querySelectorAll(`table tr td:nth-child(${clickedHeaderIndex})`);
            genesOfInterest = getGenesFromCells(colCells);

            // highlight all cells in the column
            updateHighlightStatus([...colCells], genesOfInterest, !isHighlighted);

            // highlight clicked cell so toggle works
            updateHighlightStatus(clickedCell, [], !isHighlighted);

            break;

        // Check if the clicked cell is a data cell
        default:
            const geneOfInterest = clickedCell.textContent.trim();
            updateHighlightStatus(clickedCell, [geneOfInterest], !isHighlighted);
            break;
    }

    clickedMarkerGenes.delete("");

    document.querySelector(UI.markerGenesSelectedCountElt).textContent = clickedMarkerGenes.size;
    const counterSet = new Set([...typedMarkerGenes, ...clickedMarkerGenes]);
    document.querySelector(UI.markerGenesUniqueCountElt).textContent = counterSet.size;
});

document.querySelector(UI.markerGenesManuallyEnteredElt).addEventListener("focus", (event) => {
    // Reset the manual marker gene entries
    resetManualMarkerGeneEntries();
});

document.querySelector(UI.markerGenesManuallyEnteredElt).addEventListener("change", (event) => {
    // Process the manual marker gene entries
    validateMarkerGeneSelection();
});

document.querySelector(UI.markerGenesListNameElt).addEventListener("input", (event) => {
    // enable the save button if there are genes to save and a name is entered
    document.querySelector(UI.btnSaveMarkerGeneListElt).disabled = true;
    if (event.target.value && validateMarkerGeneSelection()) {
        document.querySelector(UI.btnSaveMarkerGeneListElt).disabled = false;
    }
});

document.querySelector(UI.btnSaveMarkerGeneListElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    if (CURRENT_USER) {
        await saveMarkerGeneList();
    } else {
        createToast("You must be signed in to save a marker gene list.");
    }
    event.target.classList.remove("is-loading");
});

// Compare Genes

document.querySelector(UI.queryClusterSelectElt).addEventListener("change", (event) => {
    // If it has value, enable the compare genes button
    document.querySelector(UI.btnCompareGenesRunElt).disabled = !event.target.value;
});

document.querySelector(UI.btnCompareGenesRunElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");
    // Run the compare genes analysis
    await currentAnalysis.compareGenes.runAnalysis();
    event.target.classList.remove("is-loading");
});

document.querySelector(UI.btnCompareGenesDownloadTableFElt).addEventListener("click", (event) => {
    // Download the compare genes table for the forward comparison
    const queryId = document.querySelector(UI.queryClusterSelectElt).value;
    const referenceId = document.querySelector(UI.referenceClusterSelectElt).value;

    downloadTableAsExcel(UI.compareGenesTableFElt, `cluster_comparison_${queryId}_vs_${referenceId}.xls`);
});

document.querySelector(UI.btnCompareGenesDownloadTableRElt).addEventListener("click", (event) => {
    // Download the compare genes table for the reverse comparison
    const queryId = document.querySelector(UI.queryClusterSelectElt).value;
    const referenceId = document.querySelector(UI.referenceClusterSelectElt).value;

    downloadTableAsExcel(UI.compareGenesTableRElt, `cluster_comparison_${referenceId}_vs_${queryId}.xls`);
});

document.querySelector(UI.btnCompareGenesShowTableFElt).addEventListener("click", (event) => {
    // Show or hide the compare genes table for the forward comparison
    if (document.querySelector(UI.compareGenesTableFElt).classList.contains("is-hidden")) {
        document.querySelector(UI.compareGenesTableFElt).classList.remove("is-hidden");
        event.target.textContent = "Hide table";
        // TODO: Change mdi eye to mdi eye-off
    } else {
        document.querySelector(UI.compareGenesTableFElt).classList.add("is-hidden");
        event.target.textContent = "Show table";
        // TODO: Change mdi eye-off to mdi eye
    }
});

document.querySelector(UI.btnCompareGenesShowTableRElt).addEventListener("click", (event) => {
    // Show or hide the compare genes table for the reverse comparison
    if (document.querySelector(UI.compareGenesTableRElt).classList.contains("is-hidden")) {
        document.querySelector(UI.compareGenesTableRElt).classList.remove("is-hidden");
        event.target.textContent = "Hide table";
    } else {
        document.querySelector(UI.compareGenesTableRElt).classList.add("is-hidden");
        event.target.textContent = "Show table";
    }
});

document.querySelector(UI.compareGenesMethodSelectElt).addEventListener("change", (event) => {
    // Show or hide the p-value correction method select element based on the selected method
    /*
        p-value correction method. Used only for ‘t-test’, ‘t-test_overestim_var’,
        and ‘wilcoxon’ methods.

        The only other option is 'logreg'
    */
    if (event.target.value === "logreg") {
        document.querySelector(UI.compareGenesCorrMethodSelectElt).classList.add("is-hidden");
        document.querySelector(UI.compareGenesCorrMethodSelectElt).value = "";
    } else {
        document.querySelector(UI.compareGenesCorrMethodSelectElt).classList.remove("is-hidden");
    }

});
