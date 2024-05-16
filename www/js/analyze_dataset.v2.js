"use strict";

import UI from './classes/analysis.v2.js';  // ! Take out import when we add to CDN

let currentAnalysis = new Analysis();
let clickedMarkerGenes = new Set();
let enteredMarkerGenes = new Set();
let currentLabel = null;

// TODO:  Check font sizes on all instruction blocks
// TODO:  Check if mitochrondrial QC actually returned anything
// TODO:  Complete work on limiting the gene count
// TODO:  Louvain options are escaping their box
// TODO:  Make sure all plotting buttons either disable or show something else while the compute runs

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
        const organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        createToast("Loading dataset", "is-info");

        datasetId = newDatasetId;

        resetWorkbench();

        currentAnalysis = new Analysis({datasetId, type: "primary", datasetIsRaw: true});

        document.querySelector(UI.initialInstructionsElt).classList.add("is-hidden");

        // Technically these could load asynchronously, but logically the progress logs make more sense sequentially
        await getDatasetInfo(datasetId);
        await loadPreliminaryFigures(datasetId);
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
            elt.classList.remove('duplicate');

            const clusterLabel = elt.value.trim();

            // this means it WAS found
            if (allValues.includes(clusterLabel)) {
                dupValues.push(clusterLabel);
                elt.classList.add('duplicate');
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
    const tableStr = '';

    // Loop through the table header and add each column to the table string
    for (const elt of document.querySelector(`#${tableId} thead tr th`)){
        tableStr += `${elt.textContent}\t`;
    }
    tableStr = `${tableStr.trim()}\n`;

    // Loop through the table body and add each row to the table string
    for (const row of document.querySelector(`#${tableId} tbody tr`)){
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
    document.querySelector("#stored_analyses_c").classList.add("is-hidden");

    try {
        const data = await apiCallsMixin.getDatasetInfo(datasetId);

        const ds = new Dataset(data);

        // TODO: make Analysis.dataset the replacement for Analysis.dataset_id
        currentAnalysis.dataset = ds;

        // TODO: set dataset title and shape in the UI

        document.querySelector("#dataset_info").classList.remove("is-hidden");
        document.querySelector("#analysis_list_c").classList.remove("is-hidden");
        analysisLabels = currentAnalysis.getSavedAnalysesList(ds.id, 0, 'sc_workbench');
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
function getGenesFromCells(cells) {
    return [...cells].map(cell => cell.textContent.trim());
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
        const {data: datasetData} = await apiCallsMixin.fetchAllDatasets();

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
        throw new Error(msg);
    }

}

/**
 * Resets the manual marker gene entries.
 *
 * @returns {void}
 */
function resetManualMarkerGeneEntries() {
    clickedMarkerGenes = new Set();

    // remember which GOI are from the table
    for (const elt of document.querySelectorAll(UI.markerGenesTableHighlightedElts)) {
        clickedMarkerGenes.add(elt.textContent);
    }
}

function resetWorkbench() {
    // Performs all the steps needed to reset the workbench if the input dataset changes

    // handle any data structures
    currentAnalysis.reset();

    // Reset the analysis tool fields and some display labels
    //$("form.reset_on_change").trigger("reset");
    for (const elt of document.querySelectorAll('.reset-on-change')) {
        // TODO - replace
        elt.classList.add("is-hidden");
    }

    for (const elt of document.querySelectorAll('.empty-on-change')) {
        // TODO - replace
        elt.replaceChildren();
    }
    document.querySelector(UI.newAnalysisLabelElt).textContent = '';
    document.querySelector(UI.asvgTopGenesListElt).replaceChildren();

    // Hide any non-analysis-flow steps
    document.querySelector(UI.labeledTsneElt).classList.add("is-hidden");
    document.querySelector(UI.groupLabelsContainer).classList.add("is-hidden");

    // Toggle the tool buttons to hide them in the UI
    //$('.tooltoggle').bootstrapToggle('off');

    // Disable those analysis steps which have previous requirements
    //$('.tooltoggle').bootstrapToggle('disable');
}

/**
 * Updates the highlight status of an element and performs additional actions on the specified genes.
 *
 * @param {HTMLElement} element - The element to update the highlight status for.
 * @param {Array<string>} genes - The genes to perform additional actions on.
 * @param {boolean} addHighlight - A flag indicating whether to add or remove the highlight.
 */
function updateHighlightStatus(element, genes, addHighlight) {
    const method = addHighlight ? 'add' : 'remove';
    element.classList[method]('highlighted');
    for (const gene of genes) {
        clickedMarkerGenes[method](gene);
        currentAnalysis[`${method}GeneOfInterest`](gene);
    };
}

/**
 * Updates the manual marker gene entries based on the provided gene string.
 *
 * @param {string} geneString - The gene string to update the marker gene entries with.
 */
const updateManualMarkerGeneEntries = (geneString) => {
    enteredMarkerGenes = new Set();
    const geneSyms = geneString.split(',');

    for (let geneSym of geneSyms) {
        geneSym = geneSym.trim();
        if (geneSym) {
            enteredMarkerGenes.add(geneSym);
        }
    }

    const counterSet = new Set([...enteredMarkerGenes, ...clickedMarkerGenes]);
    document.querySelector(UI.markerGenesUniqueCountElt).textContent = counterSet.size;
    document.querySelector(UI.markerGenesEnteredCountElt).textContent = enteredMarkerGenes.size;
}

/**
 * Validates the gene selection and updates the gene cart save button accordingly.
 */
function validateMarkerGeneSelection() {
    currentAnalysis.genes_of_interest = new Set([...enteredMarkerGenes, ...clickedMarkerGenes]);
    // Only allow saving of gene cart if genes are selected
    document.querySelector(UI.btnSaveMarkerGeneCartElt).disabled = !currentAnalysis.genes_of_interest.size;
}

/* -- page entrypoint -- */

/**
 * Handles page-specific login UI updates (after the login event is triggered).
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are complete.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {
	document.getElementById("page-header-label").textContent = "Single Cell Workbench";

    for (const elt of document.querySelectorAll("#primary-nav .menu-list a.is-active")) {
        elt.classList.remove("is-active");
    }

    document.querySelector("a[tool='sc_workbench'").classList.add("is-active");


    const sessionId = CURRENT_USER.session_id;
    if (! sessionId ) {
        createToast("Not logged in so saving analyses is disabled.", "is-warning");
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
	} catch (error) {
		logErrorInConsole(error);
	}

}

/* Event listenters for elements already loaded */

// General

document.querySelector(UI.btnDeleteSavedAnalysisElt).addEventListener("click", async (event) => {
    // Delete the current analysis
    await currentAnalysis.delete();
});

document.querySelector(UI.btnMakePublicCopy).addEventListener("click", async (event) => {
    // Make a public copy of the current analysis
    await currentAnalysis.makePublicCopy();
});

document.querySelector(UI.btnSaveAnalysisElt).addEventListener("click", async (event) => {
    // Save the current analysis to the user area
    event.target.textContent = "Saving ...";
    event.target.disabled = true;
    await currentAnalysis.saveToUserArea();

});

document.querySelector(UI.btnNewAnalysisLabelSaveElt).addEventListener("click", async (event) => {
    // Save the new label to the current analysis
    currentAnalysis.label = document.querySelector(UI.btnNewAnalysisLabelElt).textContent;
    await currentAnalysis.save();
    document.querySelector(UI.btnNewAnalysisLabelContainer).classList.add("is-hidden");
});

document.querySelector(UI.btnNewAnalysisLabelCancelElt).addEventListener("click", async (event) => {
    // Reset the label to the current analysis label
    document.querySelector(UI.btnNewAnalysisLabelElt).textContent = currentAnalysis.label;
    document.querySelector(UI.btnNewAnalysisLabelContainer).classList.add("is-hidden");
});

document.querySelector(UI.analysisSelectElt).addEventListener("change", async (event) => {
    // Handle the analysis selection change

    // Grab the dataset ID from the current analysis to reuse it
    const datasetId = currentAnalysis.dataset.id;

    resetWorkbench();
    // The first analysis ID is the blank 'New' one
    if (event.target.dataset.analysisId === "0") {

        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoElt).textContent = "";
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.add("is-hidden");

        currentAnalysis = new Analysis({'datasetId': datasetId,
            'type': 'primary',
            'datasetIsRaw': true});
        // Need to reload prelim step so qc_by_mito toggle will work
        document.querySelector(UI.datasetInfoElt).classList.remove("is-hidden");
        loadPreliminaryFigures($("#dataset_id").val());

        // TODO: Update UI to match select state
        document.querySelector(UI.storedAnalysisContainer).classList.remove("is-hidden");
        return;
    }
    createToast("Loading stored analysis", "is-info");

    document.querySelector(UI.btnNewAnalysisLabelContainer).classList.add("is-hidden");

    currentAnalysis.type = event.target.dataset.analysisType;
    currentAnalysis.id = event.target.dataset.analysisId;
    currentAnalysis.datasetId = event.target.dataset.datasetId;
    currentAnalysis.getStoredAnalysis();

    if (currentAnalysis.type == 'primary') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.remove("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.add("is-hidden");
    }
    if (currentAnalysis.type == 'user_saved') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
        document.querySelector(UI.analysisStatusInfoElt).textContent = "This analysis is stored in your profile.";
        document.querySelector(UI.btnMakePublicCopyElt).classList.remove("is-hidden");
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.remove("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.add("is-hidden");
    }
    if (currentAnalysis.type == 'user_unsaved') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.remove("is-hidden");
    }
    if (currentAnalysis.type == 'public') {
        document.querySelector(UI.analysisPrimaryNotificationElt).classList.add("is-hidden");
        document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
        document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
        document.querySelector(UI.analysisStatusInfoElt).textContent = "Changes made to this public analysis will spawn a local copy within your profile.";
        document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.add("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.add("is-hidden");
    }
});

// Primary Filter

document.querySelector(UI.btnApplyPrimaryFilterElt).addEventListener("click", async (event) => {
    // Apply the primary filter to the dataset
    await currentAnalysis.primaryFilter.applyPrimaryFilter();
});

// QC by Mito

document.querySelector(UI.btnDoAnalysisQcByMitoElt).addEventListener("click", async (event) => {
    // Run the QC by Mito analysis
    await currentAnalysis.qcByMito.runAnalysis(0);
});

document.querySelector(UI.btnQbmSaveElt).addEventListener("click", async (event) => {
    // Save the QC by Mito analysis
    await currentAnalysis.qcByMito.runAnalysis(1);
    disableAndHideElement(document.querySelector(UI.btnDoAnalysisQcByMitoElt));
    event.target.textContent = "Saved";
    event.target.disabled = true;
});

document.querySelector(UI.newAnalysisLabelElt).addEventListener("focus", (event) => {
    // Reset the new analysis label
    currentLabel = event.target.textContent;
});

document.querySelector(UI.newAnalysisLabelElt).addEventListener("keyup", (event) => {
    // Update the new analysis label if it is not a duplicate
    if (analysisLabels.has(event.target.textContent)) {
        if (event.target.textContent !== currentLabel) {
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

// Select Variable Genes

document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).addEventListener("click", async (event) => {
    // Run the select variable genes analysis
    await currentAnalysis.selectMarkerGenes.runAnalysis(0);
});

document.querySelector(UI.btnAsvgSaveElt).addEventListener("click", async (event) => {
    // Save the select variable genes analysis
    await currentAnalysis.selectVariableGenes.runAnalysis(1);
    disableAndHideElement(document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt));
    event.target.textContent = "Saved";
    event.target.disabled = true;
});

// PCA

document.querySelector(UI.btnPcaRunElt).addEventListener("click", (event) => {
    // Run the PCA analysis
    currentAnalysis.checkDependenciesAndRun(currentAnalysis.pca.runAnalysis);
});

document.querySelector(UI.btnPcaTopGenesElt).addEventListener("click", (event) => {
    // Run the PCA top genes
    currentAnalysis.checkDependenciesAndRun(currentAnalysis.pca.runAnalysisTopGenes);
});

// tSNE

document.querySelector(UI.btnTsneRunElt).addEventListener("click", (event) => {
    // Run the tSNE analysis (and/or UMAP)
    currentAnalysis.checkDependenciesAndRun(currentAnalysis.tsne.runAnalysis);
});

// Clustering

document.querySelector(UI.btnClusteringRunElt).addEventListener("click", (event) => {
    // Run the clustering analysis
    currentAnalysis.checkDependenciesAndRun(currentAnalysis.clustering.runAnalysis);
});

document.querySelector(UI.btnClusteringEditRunElt).addEventListener("click", (event) => {
    // Check for duplicate labels and run the clustering edit analysis
    //  if there are none (or if the user has selected to merge clusters)

    // Remove any duplicate labels
    for (const elt of document.querySelectorAll(UI.clusterGroupLabelsInputElts)) {
        elt.classList.remove("duplicate");
    }

    // If the user has selected to merge clusters, run the clustering edit analysis
    if (document.querySelector(UI.clusteringMergeClustersElt).checked) {
        currentAnalysis.checkDependenciesAndRun(currentAnalysis.clusteringEdit.runAnalysis);
        return;
    }

    // Otherwise, check for duplicates and run the clustering analysis if there are none
    const duplicateCount = countAndHighlightDuplicates();
    if (duplicateCount === 0) {
        currentAnalysis.checkDependenciesAndRun(currentAnalysis.clusteringEdit.runAnalysis);
    }
});

// Marker Genes

document.querySelector(UI.btnMarkerGenesRunElt).addEventListener("click", (event) => {
    // Run the marker genes analysis
    currentAnalysis.checkDependenciesAndRun(currentAnalysis.markerGenes.runAnalysis);
});

document.querySelector(UI.btnVisualizeMarkerGenesElt).addEventListener("click", async (event) => {
    // Visualize the marker genes
    await performMarkerGeneVisualization();
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

document.querySelector(UI.markerGenesTableElt).addEventListener("onmouseover", (event) => {
    // Add the clickable class to the table cells
    event.target.classList.add("is-clickable");
    // ? could just put this directly in html template
});

document.querySelector(UI.markerGenesTableElt).addEventListener("click", (event) => {
    // Handle the marker genes table cell clicks

    let clickedCell = event.target.closest("td");
    const geneOfInterest = clickedCell.textContent.trim();

    if (event.target.classList.contains("js-col-idx")) {
        clickedCell = event.target.closest("th");
    }

    const isHighlighted = clickedCell.classList.contains("highlighted");
    let genesOfInterest;

    switch (true) {
        // Check if the clicked cell is a row index cell
        case clickedCell.classList.contains("js-row-idx"):
            const rowCells = [...clickedCell.parentNode.children].filter(child => child !== el);
            genesOfInterest = getGenesFromCells(rowCells);
            updateHighlightStatus(clickedCell, genesOfInterest, !isHighlighted);
            break;

        // Check if the clicked cell is a column index cell
        case clickedCell.classList.contains("js-col-idx"):
            const clickedHeaderIndex = clickedCell.cellIndex + 1;
            const colCells = document.querySelectorAll(`table tr td:nth-child(${clickedHeaderIndex})`);
            genesOfInterest = getGenesFromCells(colCells);
            updateHighlightStatus(clickedCell, genesOfInterest, !isHighlighted);
            break;

        // Check if the clicked cell is a data cell
        default:
            updateHighlightStatus(clickedCell, [geneOfInterest], !isHighlighted);
            break;
    }

    clickedMarkerGenes.delete("");

    document.querySelector(UI.markerGenesSelectedCountElt).textContent = clickedMarkerGenes.size;
    const counterSet = new Set([...enteredMarkerGenes, ...clickedMarkerGenes]);
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

// Compare Genes

document.querySelector(UI.btnCompareGenesRunElt).addEventListener("click", async (event) => {
    // Run the compare genes analysis
    await currentAnalysis.compareGenes.runAnalysis();
});

document.querySelector(UI.btnCompareGenesDownloadTableFElt).addEventListener("click", (event) => {
    // Download the compare genes table for the forward comparison
    const queryId = document.querySelector(UI.queryClusterSelectElt).value;
    const referenceId = document.querySelector(UI.referenceClusterSelectElt).value;

    downloadTableAsExcel("compare_genes_table_f", `cluster_comparison_${queryId}_vs_${referenceId}.xls`);
});

document.querySelector(UI.btnCompareGenesDownloadTableRElt).addEventListener("click", (event) => {
    // Download the compare genes table for the reverse comparison
    const queryId = document.querySelector(UI.queryClusterSelectElt).value;
    const referenceId = document.querySelector(UI.referenceClusterSelectElt).value;

    downloadTableAsExcel("compare_genes_table_r", `cluster_comparison_${referenceId}_vs_${queryId}.xls`);
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


/* -------------------------------------------------------- */

window.onload=() => {
    $('[data-toggle="tooltip"]').tooltip()

    $('.tooltoggle').bootstrapToggle('disable');

    $('#asvg_flavor_tooltip').tooltip({
        placement: 'bottom',
        html: true,
        trigger: 'hover',
        // I cannot get this to work, even though I can in this pen:
        //  https://jsbin.com/gebiju/edit?html,js,output
        delay: { "show": 100, "hide": 2000 }
    });


    $('.tooltoggle').change( function() {
        const analysis_block_id = `#analysis_${$(this).data('analysis-name')}`;

        if ( $(this).prop('checked') ) {
            $(analysis_block_id).show();
        } else {
            $(analysis_block_id).hide();
        }
    });

    $("#btn_sbs_tsne_run").click(async function(e) {
        /*
          When the user clicks to search a gene for tSNE display,
          generate a tSNE image based on the default layout.
        */
        e.preventDefault();
        $("#sbs_tsne_plot_c").empty();
        $("#sbs_tsne_gene_not_found").hide();
        $('#sbs_tsne_plot_c').show();

        show_working("Generating tSNE display");
        const ds = currentAnalysis.dataset;
        const dsp = new DatasetPanel({...ds});

        const data = await dsp.get_embedded_tsne_display(dsp.id);
        dsp.config = data.plotly_config;

        const img = document.createElement('img');
        img.className = 'img-fluid';

        get_tsne_image_data($("#sbs_tsne_gene_symbol").val(), dsp.config).then(
            data => {
                if (typeof data === 'object' || typeof data === "undefined") {
                    // If this is true, there was an error
                    $("#sbs_tsne_gene_not_found").show();
                } else {
                    // place image
                    img.src = `data:image/png;base64,${data}`
                    document.getElementById('sbs_tsne_plot_c').appendChild(img);
                }
            }
        );

        $("#btn_sbs_tsne_run").prop('disabled', false);

        done_working("tSNE visualized", false);
    });

    $(".show_analysis_renamer").click(function(e) {
        $("#new_analysis_label").val(currentAnalysis.label);
        $("#new_analysis_label_c").show(500);
    });

    // Handling all buttons manually
    $("button").click(function(e) {
        e.preventDefault();
    });


	$("#save_weighted_gene_cart").on("click", () => {
		$("#save_weighted_gene_cart").prop("disabled", true);
		if (CURRENT_USER) {
			save_weighted_gene_cart();
		} else {
			alert("You must be signed in to do that.");
		}
	});

    $("#save_marker_gene_cart").on("click", () => {
		$("#save_marker_gene_cart").prop("disabled", true);
		if (CURRENT_USER) {
			save_marker_gene_cart();
		} else {
			alert("You must be signed in to do that.");
		}
	});

    // Create observer to watch if user changes (ie. successful login does not refresh page)
    // See: https://developer.mozilla.org/en-US/docs/Web/API/MutationObserver

    // But we need to wait for navigation_bar to load first (in common.js) so do some polling
    // See: https://stackoverflow.com/q/38881301

    // Select the node that will be observed for mutations
    const target_node = document.getElementById('loggedin_controls');
    const safer_node = document.getElementById("navigation_bar");   // Empty div until loaded
    // Create an observer instance linked to the callback function
    const observer = new MutationObserver(function(mutationList, observer) {
        if (target_node) {
            populate_dataset_selection();
            this.disconnect();  // Don't need to reload once the trees are updated
        }
    });
    // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
    // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
    observer.observe(target_node || safer_node , { attributes: true });
}

function save_marker_gene_cart() {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: CURRENT_USER.session_id,
        label: $("#marker_gene_cart_name").val(),
        gctype: 'unweighted-list',
        organism_id: $("#dataset_id").data('organism-id'),
        is_public: 0
    });

    currentAnalysis.genes_of_interest.forEach((gene_id) => {
        const gene = new Gene({
            //id: gene_id,    // TODO: figure out how to get ensembl ID for this
            gene_symbol: gene_id,
        });
        gc.add_gene(gene);
    });

    gc.save(update_ui_after_marker_gene_cart_save_success, update_ui_after_marker_gene_cart_save_failure);

}

function update_ui_after_marker_gene_cart_save_success(gc) {
	$("#saved_marker_gene_cart_info_c > p").html(`Cart "${gc.label}" successfully saved.`);
	$("#saved_marker_gene_cart_info_c > p").removeClass("text-danger").addClass("text-success");
	$("#saved_marker_gene_cart_info_c").show();
    done_working("Saved marker gene cart", false);
}

function update_ui_after_marker_gene_cart_save_failure(gc, message) {
	$("#saved_marker_gene_cart_info_c > p").html("There was an issue saving the marker gene cart.");
	$("#saved_marker_gene_cart_info_c > p").removeClass("text-success").addClass("text-danger");
	$("#saved_marker_gene_cart_info_c").show();
    report_error(`Error saving gene cart: ${gc.label}`);
    report_error(message);
    $('#save_marker_gene_cart').attr("disabled", false);
}

function save_weighted_gene_cart() {
    // Return PC data from Anndata object
    $.ajax({
        type: "POST",
        url: "./cgi/get_PCs_from_anndata.cgi",
        data: {
            'dataset_id': currentAnalysis.dataset_id,
            'analysis_id': currentAnalysis.id,
            'analysis_type': currentAnalysis.type,
            'session_id': currentAnalysis.user_session_id,
        },
        datatype: "json",
    }).done((data) => {
        if (! data.success) {
            report_error(`Error getting PCs for saving: ${data.msg}`);
            $('#save_weighted_gene_cart').attr("disabled", false);
            return;
        }

        const weight_labels = data.pc_data.columns;

        // must have access to USER_SESSION_ID
        const gc = new WeightedGeneCart({
            session_id: CURRENT_USER.session_id,
            label: $("#weighted_gene_cart_name").val(),
            gctype: 'weighted-list',
            organism_id: $("#dataset_id").data('organism-id'),
            is_public: 0
        }, weight_labels
        );

        data.pc_data.index.forEach((gene_id, i) => {
            const weights = data.pc_data.data[i];
            const gene = new WeightedGene({
                id: gene_id,
                gene_symbol: data.gene_symbols[i]
            }, weights
            );
            gc.add_gene(gene);
        });

        gc.save(update_ui_after_weighted_gene_cart_save_success, update_ui_after_weighted_gene_cart_save_failure);
    }).fail((xhr, status, msg) => {
        report_error(`Error saving gene cart: ${msg}`);
        $('#save_weighted_gene_cart').attr("disabled", false);
    });

}

function update_ui_after_weighted_gene_cart_save_success(gc) {
	$("#saved_weighted_gene_cart_info_c > p").html(`Cart "${gc.label}" successfully saved.`);
	$("#saved_weighted_gene_cart_info_c > p").removeClass("text-danger").addClass("text-success");
	$("#saved_weighted_gene_cart_info_c").show();
    done_working("Saved weighted gene cart", false);
}

function update_ui_after_weighted_gene_cart_save_failure(gc, message) {
	$("#saved_weighted_gene_cart_info_c > p").html("There was an issue saving the weighted gene cart.");
	$("#saved_weighted_gene_cart_info_c > p").removeClass("text-success").addClass("text-danger");
	$("#saved_weighted_gene_cart_info_c").show();
    report_error(`Error saving gene cart: ${gc.label}`);
    report_error(message);
    $('#save_weighted_gene_cart').attr("disabled", false);
}

function get_tsne_image_data(gene_symbol, config) {
    config.colorblind_mode = CURRENT_USER.colorblind_mode;
    // then craziness: https://stackoverflow.com/a/48980526
    return axios.post(`/api/plot/${currentAnalysis.dataset_id}/tsne`, {
        gene_symbol,
        analysis: currentAnalysis,
        colorize_legend_by: config.colorize_legend_by,
        plot_type: 'tsne_static',
        plot_by_group: config.plot_by_group,
        max_columns: config.max_columns,
        horizontal_legend: config.horizontal_legend,
        x_axis: config.x_axis,
        y_axis: config.y_axis,
        analysis_owner_id: currentAnalysis.user_id,
        colors: config.colors,
        colorblind_mode: config.colorblind_mode,
        // helps stop caching issues
        timestamp: new Date().getTime()
    }).then(response => {
        return response.data.image
    }); // end axios
}


function loadPreliminaryFigures(dataset_id) {
    $("#stored_analyses_c").hide();

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_preview_primary_filter.cgi",
        data: {'dataset_id': currentAnalysis.dataset_id, 'analysis_id': currentAnalysis.id,
               'analysis_type': currentAnalysis.type, 'session_id': currentAnalysis.user_session_id
              },
        dataType: "json"
    }).done((data) => {
        $("#primary_initial_plot_loading_c").hide();

        if (data['success'] == 1) {
            $('#primary_initial_violin_c').html(`<a target="_blank" href="./datasets/${dataset_id}.prelim_violin.png"><img src="./datasets/` + dataset_id + '.prelim_violin.png" class="img-fluid img-zoomed" /></a>');
            $('#primary_initial_scatter_c').html(`<a target="_blank" href="./datasets/${dataset_id}.prelim_n_genes.png"><img src="./datasets/` + dataset_id + '.prelim_n_genes.png" class="img-fluid img-zoomed" /></a>');
            done_working("Prelim plots displayed");
        } else {
            $('#primary_initial_violin_c').html('Preliminary plots not yet generated. Continue your analysis.');
            done_working("Prelim figures missing.");
        }

        $("#primary_initial_plot_c").show(500);
    }).fail((xhr, status, msg) => {
        report_error("Failed to access dataset");
        report_error(`Failed ID was: ${dataset_id} because msg: ${msg}`);
    });
}
