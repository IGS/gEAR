"use strict";

import { Analysis, getAnalysisLabels, setAnalysisLabels } from "./classes/analysis.js";
import { UI } from "./classes/analysis-ui.js";
import { Dataset } from "./classes/dataset.js";
import { Gene, WeightedGene } from "./classes/gene.js";
import { GeneCart, WeightedGeneCart } from "./classes/genecart.v2.js";
import { DatasetTree } from "./classes/tree.js";
import { resetStepperWithHrefs } from "./helpers/stepper-fxns.js";
import { apiCallsMixin, convertToFormData, createToast, disableAndHideElement, getCurrentUser, initCommonUI, logErrorInConsole, registerPageSpecificLoginUIUpdates } from "./common.v2.js";

let currentAnalysis;
let clickedMarkerGenes = new Set();
let typedMarkerGenes = new Set();
let currentLabel = null;
let datasetId = null;

// Floating UI function alias. See https://floating-ui.com/docs/getting-started#umd
// Subject to change if we ever need these common names for other things.
const computePosition = window.FloatingUIDOM.computePosition;
const flip = window.FloatingUIDOM.flip;
const shift = window.FloatingUIDOM.shift;
const offset = window.FloatingUIDOM.offset;
const arrow = window.FloatingUIDOM.arrow;

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
        for (const elt of document.querySelectorAll(UI.currentDatasetElts)) {
            elt.textContent = e.node.title;
        }

        const newDatasetId = e.node.data.dataset_id;
        const organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        createToast("Loading dataset", "is-info");

        datasetId = newDatasetId;

        // Clear "success/failure" icons
        // these are cleared in the "resetWorkbench" function
        /*
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.classList.add("is-hidden");
        }
        for (const elt of document.getElementsByClassName("js-step-failure")) {
            elt.classList.add("is-hidden");
        }
        */

        // Hide Steppers
        for (const elt of document.querySelectorAll(UI.stepsElts)) {
            elt.classList.add("is-hidden");
        }

        // Click the Show/Hide button to collapse the tree
        document.querySelector(UI.btnToggleDatasetTreeElt).click();

        // collapse tree
        e.node.tree.expandAll(false);

        document.querySelector(UI.datasetSectionSuccessElt).classList.remove("is-hidden");
        document.querySelector(UI.datasetSectionFailedElt).classList.add("is-hidden");

        // We only want to reset these plots only if the dataset changes.
        document.querySelector(UI.primaryInitialScatterContainer).replaceChildren();
        document.querySelector(UI.primaryInitialViolinContainer).replaceChildren();

        if (currentAnalysis) {
            document.querySelector(UI.analysisSelectFailedElt).classList.add("is-hidden");
            document.querySelector(UI.analysisSelectSuccessElt).classList.add("is-hidden");
            resetWorkbench();
        }

        // Reset to the default state, so that users are forced to choose an analysis.
        document.querySelector(UI.analysisSelect).querySelector("option[data-analysis-id='-1']").selected = true;
        document.querySelector(UI.analysisWorkflowElt).classList.add("is-hidden");
        document.querySelector(UI.btnProgressGuideElt).classList.add("is-hidden");

        for (const elt of document.querySelectorAll(UI.currentAnalysisElts)) {
            elt.textContent = "None selected";
        }

        // This is a placeholder to retrieve preliminary figures which are stored in the "primary" directory
        currentAnalysis = new Analysis({id: datasetId, type: "primary", datasetIsRaw: true});

        try {
            document.querySelector(UI.analysisSelect).disabled = true;
            const labels = await currentAnalysis.getSavedAnalysesList(datasetId, -1, 'sc_workbench');
            setAnalysisLabels(labels);
        } catch (error) {
            createToast("Failed to access analyses for this dataset");
            logErrorInConsole(error);
        }

        document.querySelector(UI.primaryInitialInfoSection).classList.remove("is-hidden");
        document.querySelector(UI.primaryInitialLoadingElt).classList.remove("is-hidden");
        try {
            await getDatasetInfo(datasetId);
            await currentAnalysis.loadPreliminaryFigures(); // depends on dataset.id from getDatasetInfo
            document.querySelector(UI.analysisSelect).disabled = false;
        } catch (error) {
            logErrorInConsole(error);

            // Cannot run analyses without a dataset
            document.querySelector(UI.analysisSelect).disabled = true;
            // pass
        } finally {
            document.querySelector(UI.primaryInitialLoadingElt).classList.add("is-hidden");
        }

    })
});

/**
 * Activates a dataset in the dataset tree based on a URL parameter.
 *
 * This function checks if the specified URL parameter exists, optionally fetches additional
 * dataset information using a provided function, and then activates the corresponding dataset
 * node in the dataset tree. If the dataset cannot be found or accessed, a toast notification
 * is displayed and an error is thrown.
 *
 * @async
 * @param {URLSearchParams} urlParams - The URLSearchParams object containing the URL parameters.
 * @param {string} paramName - The name of the URL parameter to look for.
 * @param {function} [fetchInfoFn] - Optional async function to fetch dataset info using the parameter value.
 *        Should return a Promise that resolves to an array of objects containing a `dataset_id` property.
 * @throws {Error} If the dataset cannot be accessed or found in the dataset tree.
 */
const activateDatasetFromParam = async (urlParams, paramName, fetchInfoFn) => {
    if (!urlParams.has(paramName)) {
        return;
    }
    const paramValue = urlParams.get(paramName);
    let linkedDatasetId;
    try {
        if (fetchInfoFn) {
            const data = await fetchInfoFn(paramValue);
            linkedDatasetId = data.datasets[0].id;
            if (!linkedDatasetId) {
                throw new Error(`Accessible dataset for ${paramName} ${paramValue} was not found`);
            }
        } else {
            linkedDatasetId = paramValue;
        }
    } catch (error) {
        createToast(error.message);
        throw new Error(error);
    }

    try {
        // find DatasetTree node and trigger "activate"
        const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
        foundNode.setActive(true, {focusTree:true});
        datasetTree.tree.setActiveNode(foundNode);
        datasetTree.selectCallback({node: foundNode});  // manually trigger the "activate" event.
        datasetId = linkedDatasetId;
    } catch (error) {
        createToast(`Dataset id ${linkedDatasetId} was not found as a public/private/shared dataset`);
        throw new Error(error);
    }
}

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
    try {
        const data = await apiCallsMixin.fetchDatasetInfo(datasetId);

        const ds = new Dataset(data);

        currentAnalysis.dataset = ds;

        document.querySelector(UI.primaryFilterSection).classList.remove("is-hidden");
        document.querySelector(UI.selectedDatasetShapeInitialElt).textContent = currentAnalysis.dataset.shape();

        createToast("Dataset loaded", "is-success");
    } catch (error) {
        createToast("Failed to access dataset");
        logErrorInConsole(`Failed ID was: ${datasetId} because msg: ${error.message}`);
    }
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
 * Loads the dataset tree by fetching dataset information from the curator API.
 * Populates the userDatasets, sharedDatasets, and domainDatasets arrays with dataset information.
 * Generates the dataset tree using the generateTree method of the datasetTree object.
 * @throws {Error} If there is an error fetching the dataset information.
 */
const loadDatasetTree = async (shareId) => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const datasetData = await apiCallsMixin.fetchAllDatasets(shareId);

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
        document.querySelector(UI.datasetSectionSuccessElt).classList.add("is-hidden");
        document.querySelector(UI.datasetSectionFailedElt).classList.remove("is-hidden");
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
}

/**
 * Saves the marker gene list.
 * @returns {void}
 */
const saveMarkerGeneList = async () => {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: getCurrentUser()?.session_id,
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

    await gc.save(updateUiAfterMarkerGeneListSaveSuccess, updateUiAfterMarkerGeneListSaveFailure);
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
            'session_id': currentAnalysis.analysisSessionId,
        }));

        if (!data.success || data.success < 1) {
            const message = data.msg || "Unknown error";
            throw new Error(message);
        }

        const weightLabels = data.pc_data.columns;

        const geneList = new WeightedGeneCart({
                session_id: getCurrentUser()?.session_id,
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

        await geneList.save(updateUiAfterPcaGeneListSaveSuccess, updateUiAfterPcaGeneListSaveFailure);

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

    // Enable the button if there are genes entered
    document.querySelector(UI.btnVisualizeMarkerGenesElt).disabled = !counterSet.size;
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

    const sessionId = getCurrentUser()?.session_id;
    if (! sessionId ) {
        createToast("Not logged in so saving analyses is disabled.", "is-warning");
        document.querySelector(UI.btnSaveAnalysisElt).disabled = true;

        // TODO: Other actions
    }

	try {

        // If brought here by the "gene search results" page, curate on the dataset ID that referred us
        const urlParams = new URLSearchParams(window.location.search);
        const shareId = urlParams.get("share_id");
		await loadDatasetTree(shareId);

        // Usage inside handlePageSpecificLoginUIUpdates
        if (urlParams.has("share_id")) {
            return await activateDatasetFromParam(urlParams, "share_id", async (shareId) =>
                await apiCallsMixin.fetchDatasetListInfo({permalink_share_id: shareId})
            );
        } else if (urlParams.has("dataset_id")) {
    		// Legacy support for dataset_id

            await activateDatasetFromParam(urlParams, "dataset_id");
        }

        // ? This could be used to pre-select an analysis
        //currentAnalysis = new Analysis();

	} catch (error) {
		logErrorInConsole(error);
	}
}
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

// Pre-initialize some stuff
await initCommonUI();

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

for (const button of document.querySelectorAll(UI.analysisDownloadElts)) {
    button.addEventListener("click", (event) => {
        // Download the current analysis
        currentAnalysis.download();
    });
}

// Show the "rename" analysis label input when the button is clicked
for (const button of document.querySelectorAll(UI.analysisRenameElts)) {
    button.addEventListener("click", (event) => {

        // remove existing popovers
        const existingPopover = document.getElementById('rename-analysis-popover');
        if (existingPopover) {
            existingPopover.remove();
        }

        // Create popover content
        const popoverContent = document.createElement('article');
        popoverContent.id = 'rename-analysis-popover';
        popoverContent.classList.add("message", "is-dark");
        popoverContent.setAttribute("role", "tooltip");
        popoverContent.style.width = "500px";
        popoverContent.innerHTML = `
            <div class='message-header'>
                <p>Rename collection</p>
            </div>
            <div class='message-body'>
                <p>Please provide a new name for this analysis</p>
                <div class='field'>
                    <div class='control'>
                        <input id='new-analysis-label' class='input' type='text' placeholder='Analysis name'>
                    </div>
                </div>
                <div id="duplicate-label-warning" class='is-hidden help has-text-danger-dark'>
                    This analysis label has already been used.
                </div>
                <div class='field is-grouped' style='width:250px'>
                    <p class="control">
                        <button id='confirm-analysis-rename' class='button is-dark' disabled>Update</button>
                    </p>
                    <p class="control">
                        <button id='cancel-analysis-rename' class='button' value='cancel_rename'>Cancel</button>
                    </p>
                </div>
            </div>
            <div id="arrow"></div>
        `;

        // append element to DOM to get its dimensions
        document.body.appendChild(popoverContent);

        document.getElementById("new-analysis-label").value = currentAnalysis.label;
        for (const elt of document.querySelectorAll(UI.currentAnalysisElts)) {
            elt.textContent = currentAnalysis.label;
        }

        const arrowElement = document.getElementById('arrow');

        // Create popover (help from https://floating-ui.com/docs/tutorial)
        computePosition(button, popoverContent, {
            placement: 'bottom',
            middleware: [
                flip(), // flip to top if there is not enough space on butotm
                shift(), // shift the popover to the right if there is not enough space on the left
                offset(5), // offset relative to the button
                arrow({ element: arrowElement }) // add an arrow pointing to the button
            ],
        }).then(({ x, y, placement, middlewareData }) => {
            // Position the popover
            Object.assign(popoverContent.style, {
                left: `${x}px`,
                top: `${y}px`,
            });
            // Accessing the data
            const { x: arrowX, y: arrowY } = middlewareData.arrow;

            // Position the arrow relative to the popover
            const staticSide = {
                top: 'bottom',
                right: 'left',
                bottom: 'top',
                left: 'right',
            }[placement.split('-')[0]];

            // Set the arrow position
            Object.assign(arrowElement.style, {
                left: arrowX != null ? `${arrowX}px` : '',
                top: arrowY != null ? `${arrowY}px` : '',
                right: '',
                bottom: '',
                [staticSide]: '-4px',
            });
        });

        document.getElementById("new-analysis-label").addEventListener("keyup", (event) => {
            // Update the new analysis label if it is not a duplicate
            if (getAnalysisLabels().has(event.target.value.trim())) {
                if (event.target.value.trim() !== currentLabel) {
                    event.target.classList.add("duplicate");
                    document.getElementById("confirm-analysis-rename").disabled = true;
                    document.getElementById("duplicate-label-warning").classList.remove("is-hidden");
                }
                return;
            }

            document.getElementById("duplicate-label-warning").classList.add("is-hidden");
            event.target.classList.remove("duplicate");
            document.getElementById("confirm-analysis-rename").disabled = false;
        });

        document.getElementById("new-analysis-label").addEventListener("focus", (event) => {
            // Reset the new analysis label
            currentLabel = event.target.value.trim();
        });

        document.getElementById("cancel-analysis-rename").addEventListener("click", async (event) => {
            // Reset the label to the current analysis label
            document.getElementById("new-analysis-label").value = currentAnalysis.label;
            popoverContent.remove();
        });


        document.getElementById("confirm-analysis-rename").addEventListener("click", async (event) => {
            event.target.classList.add("is-loading");
            currentAnalysis.label = document.getElementById("new-analysis-label").value;

            try {
                // Save the new label to the current analysis
                await currentAnalysis.save();
            } catch (error) {
                // pass - handled in the save method
            } finally {
                event.target.classList.remove("is-loading");
                popoverContent.remove();
            }
        });

    });
}





/**
 * Creates a confirmation popover for deleting an analysis.
 */
const deleteButtons = document.querySelectorAll(UI.analysisDeleteElts);
for (const button of deleteButtons) {
    button.addEventListener('click', (e) => {
        // remove existing popovers
        const existingPopover = document.getElementById('delete-analysis-popover');
        if (existingPopover) {
            existingPopover.remove();
        }

        // Create popover content
        const popoverContent = document.createElement('article');
        popoverContent.id = 'delete-analysis-popover';
        popoverContent.classList.add("message", "is-danger");
        popoverContent.setAttribute("role", "tooltip");
        popoverContent.style.width = "500px";
        popoverContent.innerHTML = `
            <div class='message-header'>
                <p>Delete this analysis</p>
            </div>
            <div class='message-body'>
                <p>Are you sure you want to delete this analysis? This will affect any saved dataset displays using this analysis for yourself and for other users.</p>
                <div class='field is-grouped' style='width:250px'>
                    <p class="control">
                        <button id='confirm-analysis-delete' class='button is-danger'>Delete</button>
                    </p>
                    <p class="control">
                        <button id='cancel-analysis-delete' class='button' value='cancel_delete'>Cancel</button>
                    </p>
                </div>
            </div>
            <div id="arrow"></div>
        `;

        // append element to DOM to get its dimensions
        document.body.appendChild(popoverContent);

        const arrowElement = document.getElementById('arrow');

        // Create popover (help from https://floating-ui.com/docs/tutorial)
        computePosition(button, popoverContent, {
            placement: 'bottom',
            middleware: [
                flip(), // flip to top if there is not enough space on butotm
                shift(), // shift the popover to the right if there is not enough space on the left
                offset(5), // offset relative to the button
                arrow({ element: arrowElement }) // add an arrow pointing to the button
            ],
        }).then(({ x, y, placement, middlewareData }) => {
            // Position the popover
            Object.assign(popoverContent.style, {
                left: `${x}px`,
                top: `${y}px`,
            });
            // Accessing the data
            const { x: arrowX, y: arrowY } = middlewareData.arrow;

            // Position the arrow relative to the popover
            const staticSide = {
                top: 'bottom',
                right: 'left',
                bottom: 'top',
                left: 'right',
            }[placement.split('-')[0]];

            // Set the arrow position
            Object.assign(arrowElement.style, {
                left: arrowX != null ? `${arrowX}px` : '',
                top: arrowY != null ? `${arrowY}px` : '',
                right: '',
                bottom: '',
                [staticSide]: '-4px',
            });
        });

        // Add event listener to cancel button
        document.getElementById('cancel-analysis-delete').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-analysis-delete').addEventListener('click', async (event) => {
            event.target.classList.add("is-loading");

            try {
                // Delete the current analysis
                await currentAnalysis.delete();
            } catch (error) {
                // pass - handled in the delete method
            } finally {
                event.target.classList.remove("is-loading");
                popoverContent.remove();
            }
        });
    });
}


// Handle the analysis selection change
document.querySelector(UI.analysisSelect).addEventListener("change", async (event) => {

    // Grab the dataset ID from the current analysis to reuse it
    const datasetObj = currentAnalysis.dataset;

    for (const elt of document.querySelectorAll(UI.currentAnalysisElts)) {
        elt.textContent = event.target.selectedOptions[0].textContent;
    }

    document.querySelector(UI.deNovoStepsElt).classList.add("is-hidden");
    document.querySelector(UI.primaryStepsElt).classList.add("is-hidden");
    document.querySelector(UI.btnProgressGuideElt).classList.add("is-hidden");
    document.querySelector(UI.progressGuideElt).classList.add("is-hidden");

    // Analysis ID -1 is "select an analysis"
    if (event.target.value === "-1") {
        document.querySelector(UI.analysisWorkflowElt).classList.add("is-hidden");
        document.querySelector(UI.analysisSelectFailedElt).classList.remove("is-hidden");
        document.querySelector(UI.analysisSelectSuccessElt).classList.add("is-hidden");
        return;
    }

    document.querySelector(UI.analysisWorkflowElt).classList.remove("is-hidden");
    document.querySelector(UI.analysisSelectFailedElt).classList.add("is-hidden");
    document.querySelector(UI.analysisSelectSuccessElt).classList.remove("is-hidden");

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
        resetStepperWithHrefs(UI.primaryFilterSection);

        // Icon counter adjustment from if we were previously in primary analysis (and tSNE was calculated)
        document.querySelector(`${UI.markerGenesSection} i`).classList.replace("mdi-numeric-3-circle", "mdi-numeric-9-circle")
        document.querySelector(`${UI.markerGenesSection} i`).classList.replace("mdi-numeric-4-circle", "mdi-numeric-9-circle")

        // Jump to the primary filter step
        document.querySelector(`a[href='${UI.primaryFilterSection}']`).click();
        return;
    }
    createToast("Loading stored analysis", "is-info");

    resetWorkbench();

    const selectedOption = event.target.selectedOptions[0];
    currentAnalysis.type = selectedOption.dataset.analysisType;
    currentAnalysis.id = selectedOption.dataset.analysisId;
    currentAnalysis.analysisSessionId = selectedOption.dataset.analysisSessionId;

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
        createToast("This analysis is stored in your private user profile.", "is-info", true);
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
        createToast("Changes made to this public analysis will create a local private copy within your profile.", "is-info", true);
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
    }

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

document.querySelector(UI.topPcaGenesElt).addEventListener("input", (event) => {
    // enable the 'plot top genes' button if at least one number is entered
    document.querySelector(UI.btnPcaTopGenesElt).disabled = true;
    if (!event.target.value) {
        return;
    }
    const splitValues = event.target.value.split(",");
    // check if first value is an int
    if (parseInt(splitValues[0])) {
        document.querySelector(UI.btnPcaTopGenesElt).disabled = false;
    };
});

document.querySelector(UI.btnPcaTopGenesElt).addEventListener("click", async (event) => {
    event.target.classList.add("is-loading");

    // validate that all values are integers
    const topGenes = document.querySelector(UI.topPcaGenesElt).value.split(",").map(x => parseInt(x));
    if (topGenes.some(isNaN)) {
        createToast("Please enter a comma-separated list of integers.");
        event.target.classList.remove("is-loading");
        return;
    }

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
    if (getCurrentUser()) {
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
    if (getCurrentUser()) {
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

document.querySelector(UI.btnToggleDatasetTreeElt).addEventListener("click", (event) => {
    // Toggle the dataset selection div
    const selectionDiv = document.querySelector(UI.datasetSelectionContainer);
    if (selectionDiv.classList.contains("is-hidden")) {
        selectionDiv.classList.remove("is-hidden");
        event.target.textContent = "Collapse dataset selection tool";
    } else {
        selectionDiv.classList.add("is-hidden");
        event.target.textContent = "Expand dataset selection tool";
    }
});
