'use strict';

import { apiCallsMixin, createToast, disableAndHideElement, enableAndShowElement, getCurrentUser, getUrlParameter, initCommonUI, logErrorInConsole, rebindUrlParam, registerPageSpecificLoginUIUpdates } from "./common.v2.js?v=cbfcd86";
import { datasetCollectionState, fetchDatasetCollections, registerEventListeners as registerDatasetCollectionEventListeners, selectDatasetCollection } from "../include/dataset-collection-selector/dataset-collection-selector.js?v=cbfcd86";
import { fetchPatternsData, getFlatPatternCartData, getSelectedPattern, populatePatternWeights, registerEventListeners as registerPatternEventListeners, selectPatternWeights, setSelectedPattern } from "../include/pattern-collection-selector/pattern-collection-selector.js?v=cbfcd86";
import { TileGrid } from "./classes/tilegrid.js?v=cbfcd86";

let selectedPattern;
let urlParamsPassed = false;
let isMulti = false;
let tilegrid = null;
let svgScoringMethod = 'gene';
let projectionOpts = { patternSource: null, algorithm: null, gctype: null };
let weightedGeneData = null;
let datasetShareId = null;
let layoutShareId = null;

// if certain legacy or shorthand URL parameters are passed, change the parameter to the new ones
const urlParams = new URLSearchParams(window.location.search);

// imported from pattern-collection-selector.js
// selectedPattern = {shareId: null, label: null, gctype: null, selectedWeights: []};

/**
 * Creates a proxy object for the selected pattern.
 *
 * @param {Object} selectedPattern - The selected pattern object.
 * @returns {Proxy} - The proxy object for the selected pattern.
 */
const createSelectedPatternProxy = (selectedPattern) => {
    return new Proxy(selectedPattern, {
        set: (target, key, value) => {
            target[key] = value;
            const algorithmElt = document.getElementById('algorithm');

            // NOTE: When checking keys, if multiple keys are set at once, the order of the if statements matters
            // The Proxy keys are triggered in order they were set in the object.

            if (key === "selectedWeights") {
                enableAndShowElement(document.getElementById("single-multi-multi"), true);
                if (!value.length) {
                    // Reset button was hit
                    enableAndShowElement(document.getElementById("single-multi-multi"), true);
                } else if (value.length < 2) {
                    disableAndHideElement(document.getElementById("single-multi-multi"), true);
                    document.getElementById("single-multi-single").checked = true;
                    isMulti = false;
                }
                // Enable or disable the "binary" option if all selectedWeights have "binary" property set to True
                const binary = value.every((w) => w.binary);
                algorithmElt.querySelector('option[value="binary"]').disabled = !binary;

                // if "binary" is selected, reset the algorithm "select" box to "pca"
                if (algorithmElt.value === 'binary') {
                    algorithmElt.value = 'pca';
                }

            } else if (key === "gctype") {
                // Adjust algorithm options based on gctype
                algorithmElt.querySelector('option[value="nmf"]').disabled = false;
                algorithmElt.querySelector('option[value="fixednmf"]').disabled = false;

                if (value === "unweighted-list") {
                    algorithmElt.querySelector('option[value="nmf"]').disabled = true;
                    algorithmElt.querySelector('option[value="fixednmf"]').disabled = true;
                    // Reset the algorithm "select" box
                    if (['nmf', 'fixednmf'].includes(algorithmElt.value)) {
                        algorithmElt.value = 'pca';
                    }
                }
            }
            return true;
        }
    });
};

/**
 * Builds the state URL with the selected parameters.
 * @returns {string} The state URL.
 */
const buildStateUrl = () => {

    // Create a new URL object (with no search params)
    const url = new URL('/projection.html', window.location.origin);

    // Add the projection algorithm to the URL
    const algorithm = document.getElementById('algorithm').value;
    url.searchParams.set('projection_algorithm', algorithm);

    const zscore = document.getElementById('zscore').checked ? 1 : 0;
    url.searchParams.set('zscore', zscore);

    // Add the minclip value to the URL if checked.
    // minclip can be a number
    const minclip = document.getElementById('minclip').checked;
    if (minclip) {
        url.searchParams.set('minclip', 0);
    }

    // Add the multipattern_plots value to the URL
    const multipatternPlots = document.getElementById('single-multi-multi').checked ? 1 : 0;
    url.searchParams.set('multipattern_plots', multipatternPlots);

    // Add the pattern source to the URL
    url.searchParams.set('projection_source', selectedPattern.shareId);

    // Add the dataset collection to the URL
    if (datasetShareId) {
        url.searchParams.append('share_id', datasetShareId);
    } else if (datasetCollectionState.selectedShareId) {
        url.searchParams.append('layout_id', datasetCollectionState.selectedShareId);
    }

    // Add the selected pattern weights to the URL
    const weights = selectedPattern.selectedWeights.map((w) => w.label);
    url.searchParams.set('projection_patterns', weights.join(','));

    return url.toString();
};

/**
 * Sorts an array of strings in ascending order based on the numeric value at the end of each string. (i.e. PC1, PC2, etc.)
 *
 * @param {string} a - The first string to compare.
 * @param {string} b - The second string to compare.
 * @returns {number} The difference between the numeric values at the end of the strings.
 */
const customNumericSort = (a, b) => {
    // NOTE: Ignores the leading string altogether, so this still applies even if that is not consistent.
    return (Number(a.match(/(\d+)$/g)[0]) - Number((b.match(/(\d+)$/g)[0])));
};

/**
 * Populates the pattern results list with weights.
 */
const populatePatternResultsList = () => {
    const template = document.getElementById('tmpl-pattern-result-item');
    document.getElementById('pattern-result-list').innerHTML = '';

    // sort the selectedWeights array based on the numeric value at the end of each string
    // If this does not work, just sort alphabetically
    let sortedLabels;
    try {
        sortedLabels = selectedPattern.selectedWeights.map((weight) => weight.label).sort(customNumericSort);
    } catch (error) {
        sortedLabels = selectedPattern.selectedWeights.map((weight) => weight.label).sort();
    }

    for (const label of sortedLabels) {
        const row = template.content.cloneNode(true);
        row.querySelector('li').innerHTML = label;
        row.querySelector('li').dataset.weight = label;
        document.getElementById('pattern-result-list').appendChild(row);

        const thisRow = document.querySelector(`.pattern-result-list-item[data-weight="${label}"]`);
        thisRow.addEventListener('click', (event) => {

            // remove is-selected from all the existing rows, then add it to this one
            const rows = document.getElementsByClassName('pattern-result-list-item');
            for (const row of rows) {
                row.classList.remove('is-selected');
            }

            event.currentTarget.classList.add('is-selected');
            selectPatternWeightResult(label);
        });
    }
};

/**
 * Parses the URL parameters and updates the UI based on the values.
 */
const parsePatternCartURLParams = async () => {

    // if projection algorithm is passed, set it in #algorithm
    const projectionAlgorithm = getUrlParameter('projection_algorithm', urlParams);
    if (projectionAlgorithm) {
        document.getElementById('algorithm').value = projectionAlgorithm;
    }

    // if zscore is passed, set it in #zscore
    const zscore = getUrlParameter('zscore', urlParams);
    if (zscore === '1') {
        document.getElementById('zscore').checked = true;
    }

    // if minclip is passed, set it in #minclip
    const minclip = getUrlParameter('minclip', urlParams);
    if (minclip) {
        document.getElementById('minclip').checked = true;
    }

    // single or multiple pattern view (convert to boolean)?
    // NOTE: This will be adjusted if the pattern only has one weight
    const isMultiParam = getUrlParameter('multipattern_plots', urlParams);
    isMulti = isMultiParam === '1';
    if (isMulti) {
        document.getElementById('single-multi-multi').checked = true;
    } else {
        document.getElementById('single-multi-single').checked = true;
    }

    // handle passed pattern lists
    const pattern = getUrlParameter('projection_source', urlParams);
    if (!pattern) {
        return;
    }
    urlParamsPassed = true;
    const foundPattern = getFlatPatternCartData().find((p) => p.share_id === pattern);
    if (!foundPattern) {
        console.warn(`Pattern ${pattern} not found in pattern cart data. Perhaps the user does not have access to it.`);
        return;
    }
    setSelectedPattern({ shareId: foundPattern.share_id, label: foundPattern.label, gctype: foundPattern.gctype, selectedWeights: [] });

    // Update proxy so that multi-gene radio button can be enabled/disabled
    selectedPattern = createSelectedPatternProxy(getSelectedPattern());

    // we cannot the click event, since the pattern list items only render when an intiial category is selected
    // so we need to manually populate the pattern weights
    await populatePatternWeights();

    // If no weights were passed, select the first weight for the pattern
    const rows = document.getElementsByClassName('dropdown-weight-item');
    const labels = Array.from(rows).map((row) => row.dataset.label);

    // handle manually-entered pattern symbols
    const urlWeights = getUrlParameter('projection_patterns', urlParams);
    if (urlWeights) {
        // Cannot have weights without a source pattern
        const urlLabels = urlWeights.split(',');

        const labelsToDeselect = labels.filter((label) => !urlLabels.includes(label));
        // deselect all weights that are not in the urlWeights
        selectPatternWeights(labelsToDeselect);
    }

    // click "proceed" button in pattern selector to update the UI
    document.getElementById('dropdown-pattern-list-proceed').click();
};

/**
 * Parses the URL parameters to extract the dataset collection information.
 * @returns {Promise<void>} A promise that resolves once the dataset collection information is parsed.
 */
const parseDatasetCollectionURLParams = () => {
    // handle passed dataset collection
    const layoutShareId = getUrlParameter('layout_id');

    if (!layoutShareId) {
        return;
    }

    datasetCollectionState.selectedShareId = layoutShareId;
    datasetCollectionState.selectedLabel = datasetCollectionState.labelIndex[layoutShareId];
    document.querySelector('#dropdown-dc-selector-label').innerHTML = datasetCollectionState.selectedLabel;
};

/**
 * Selects a pattern weight and performs various actions based on the selected weight.
 * For unweighted patterns, "unweighted" itself is the pattern weight which is all 1 values.
 * @param {string} label - The selected weight label.
 */
const selectPatternWeightResult = async (label) => {

    // get the selected pattern object
    const obj = selectedPattern.selectedWeights.find((w) => w.label === label);

    // if projection algorithm is "nmf", then hide the top_down genes
    const projectionAlgorithm = document.getElementById('algorithm').value;

    const topDownGenesElement = document.getElementById('top-down-genes');

    // Hide the unweighted gene list button
    document.getElementById("btn-view-unweighted-genes").classList.add('is-hidden');

    // Hide the element initially
    topDownGenesElement.classList.add('is-hidden');

    // Reset weighted gene data
    weightedGeneData = null;

    if (projectionAlgorithm === 'nmf') {
        // No additional action needed as the element is already hidden
    } else if (obj.top_down) {
        // Show the element if 'top_down' property is true
        topDownGenesElement.classList.remove('is-hidden');
    }

    // clear the top-up and top-down genes
    document.querySelector("#top-up-genes p").textContent = '';
    document.querySelector("#top-down-genes p").textContent = '';

    // populate top-up and top-down with the array of genes for that weight
    if (obj.top_up) {
        // format to add a space after each comma
        const topUp = obj.top_up.split(',').join(', ');
        document.querySelector("#top-up-genes p").textContent = topUp;
    }
    if (obj.top_down) {
        const topDown = obj.top_down.split(',').join(', ');
        document.querySelector("#top-down-genes p").textContent = topDown;
    }

    // if isMulti=false, show top_up and top_down genes
    document.getElementById("top-genes-c").classList.remove('is-hidden');
    if (isMulti || selectedPattern.gctype === "unweighted-list") {
        document.getElementById("top-genes-c").classList.add('is-hidden');
    }

    if (selectedPattern.gctype === "weighted-list") {
        document.getElementById("btn-view-weighted-genes").classList.remove("is-hidden");
        try {
            const data = await apiCallsMixin.fetchPatternWeightedGenes(selectedPattern.shareId, label);
            weightedGeneData = data;
        } catch (error) {
            logErrorInConsole(error);
            document.getElementById("btn-view-weighted-genes").classList.add("is-hidden");
        }
    } else if (selectedPattern.gctype === "unweighted-list") {
        // populate the gene list from this cart
        const geneListMemberData = await apiCallsMixin.fetchGeneCartMembers(selectedPattern.shareId,);
        const geneData = [];
        for (const member of geneListMemberData.gene_symbols) {
            geneData.push({ gene: member.label, weight: 1 });
        }
        weightedGeneData = geneData;

        document.getElementById("btn-view-unweighted-genes").classList.remove('is-hidden');
    }

    // Other things can be called next, such as plotting calls
    if (tilegrid) {
        // Revert back to "#result-panel-grid" display before rendering the new gene displays
        document.getElementById("result-panel-grid").classList.remove("is-hidden");
        document.getElementById("zoomed-panel-grid").classList.add("is-hidden");

        const minclip = document.getElementById('minclip').checked ? 0 : null;
        await tilegrid.renderDisplays(label, isMulti, svgScoringMethod, minclip, projectionOpts);
    }
};

/**
 * Sets up the tile grid for a given shareId and type.
 *
 * @param {string} shareId - The shareId to set up the tile grid for.
 * @param {string} [type="layout"] - The type of the tile grid. Defaults to "layout".
 * @returns {Promise<TileGrid>} - A promise that resolves to the initialized TileGrid object.
 */
const setupTileGrid = async (shareId, type = "layout") => {

    // Cannot proceed without a shareId
    if (!shareId) {
        return;
    }

    const tilegrid = new TileGrid(shareId, type, "#result-panel-grid");
    try {
        tilegrid.datasets = await tilegrid.getDatasets();
        tilegrid.layout = await tilegrid.getLayout();
        await tilegrid.addAllDisplays();

        tilegrid.applyTileGrid(isMulti);

        const algorithm = document.getElementById('algorithm').value;
        const zscore = document.getElementById('zscore').checked;
        const minclip = document.getElementById('minclip').checked ? 0 : null;

        // create projectionOpts object out of selectedPattern.shareId, algorithm, and selectedPattern.gctype
        projectionOpts = {
            patternSource: selectedPattern.shareId,
            algorithm,
            gctype: selectedPattern.gctype,
            zscore
        };

        for (const tile of tilegrid.tiles) {
            tile.enableProjectR();
        }

        // NOTE - the tilegrid.renderDisplays() call below can check and use the first array element of the geneCollectionState.selectedGenes array if single_pattern
        // We do not render for single-gene searches because the first pattern result is "clicked" and the tilegrid is rendered in the event listener.

        if (isMulti && selectedPattern.selectedWeights.length) {
            // create array of selected weight labels
            const selectedWeights = Array.from(selectedPattern.selectedWeights).map((w) => w.label);
            await tilegrid.renderDisplays(selectedWeights, isMulti, svgScoringMethod, minclip, projectionOpts);
        }
    } catch (error) {
        logErrorInConsole(error);
    } finally {
        return tilegrid;
    }
};

/**
 * Validates the projection search form.
 *
 * @returns {boolean} Returns true if the form is valid, otherwise false.
 */
const validateProjectionSearchForm = () => {

    // User passed in a single dataset share ID.
    if (datasetShareId) {
        return true;
    }

    // User must have either selected a pattern list or entered patterns manually. Either of these
    // will populate the selected_patterns array
    document.querySelector('#dropdown-pattern-lists button').classList.remove('is-danger');
    if (selectedPattern.shareId === null) {
        createToast('Please enter at least one pattern source to proceed');
        document.querySelector('#dropdown-pattern-lists button').classList.add('is-danger');
        return false;
    }

    // Check if the user has selected any dataset collections
    document.querySelector('#dropdown-dc button').classList.remove('is-danger');
    if (!datasetCollectionState.selectedShareId) {
        createToast('Please select at least one dataset to proceed');
        document.querySelector('#dropdown-dc button').classList.add('is-danger');
        return false;
    }

    // If multi, check that at least two weights are selected
    if (document.getElementById('single-multi-multi').checked && selectedPattern.selectedWeights.length < 2) {
        createToast('Please select at least two patterns to proceed');
        return false;
    }

    return true;
};

/**
 * Handles the UI updates specific to the page login.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are completed.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {

    // Set the page header title
    document.getElementById('page-header-label').textContent = 'Projection Search';
    datasetShareId = getUrlParameter('share_id');
    layoutShareId = getUrlParameter('layout_id');

    // There are some shorthand URL parameters (not on the shorthand URL) that need to be converted to the longform
    rebindUrlParam(urlParams, "multi", "multipattern_plots");
    rebindUrlParam(urlParams, "c", "projection_source");
    rebindUrlParam(urlParams, "ptrns", "projection_patterns");
    rebindUrlParam(urlParams, "algo", "projection_algorithm");


    // Register event listeners for pattern and dataset collection selectors
    registerPatternEventListeners(apiCallsMixin);
    registerDatasetCollectionEventListeners(apiCallsMixin, getCurrentUser());

    selectedPattern = createSelectedPatternProxy(getSelectedPattern());

    // add event listener for when the submit-projection-search button is clicked
    document.getElementById('submit-projection-search').addEventListener('click', async (event) => {

        const currentTarget = event.currentTarget;
        currentTarget.classList.add("is-loading");

        const status = validateProjectionSearchForm();

        if (!status) {
            console.info("Aborting search");
            event.currentTarget.classList.remove("is-loading");
            return;
        }

        document.getElementById("result-panel-initial-notification").classList.add('is-hidden');
        document.getElementById("result-panel-loader").classList.remove('is-hidden');

        // update multi/single pattern
        isMulti = document.getElementById('single-multi-multi').checked;

        populatePatternResultsList();

        // if multi, clear the selected pattern symbol and hide the pattern-result-list container
        document.getElementById("pattern-result-list-c").classList.remove('is-hidden');
        document.getElementById("scoring-method-div").classList.remove('is-hidden');
        if (isMulti) {
            document.getElementById("pattern-result-list-c").classList.add('is-hidden');
            document.getElementById("scoring-method-div").classList.add('is-hidden');
        }

        try {
            const setupTileGridFn = (datasetShareId) ? setupTileGrid(datasetShareId, "dataset") : setupTileGrid(datasetCollectionState.selectedShareId);
            tilegrid = await setupTileGridFn;

            // auto-select the first pattern in the list
            const firstPattern = document.querySelector('.pattern-result-list-item');
            if (!isMulti && firstPattern) {
                firstPattern.click();
            }

        } catch (error) {
            logErrorInConsole(error);
            return;
        } finally {
            currentTarget.classList.remove("is-loading");
            document.getElementById("result-panel-loader").classList.add('is-hidden');

        }

        const url = buildStateUrl();
        // add to state history
        history.pushState(null, '', url);
    });

    // Wait until all pending API calls have completed before checking if we need to search
    document.getElementById("submit-projection-search").classList.add("is-loading");
    try {
        const pattern = getUrlParameter('projection_source', urlParams);

        // SAdkins note - Promise.all fails fast,
        // but Promise.allSettled waits until all resolve/reject and lets you know which ones failed
        const [cartResult, dcResult,] = await Promise.all([
            fetchPatternsData(pattern),
            fetchDatasetCollections(layoutShareId),
        ]);

        parseDatasetCollectionURLParams();
        await parsePatternCartURLParams();

        // Should help with lining things up on index page
        document.getElementById("dropdown-dc").classList.remove("is-right");

    } catch (error) {
        logErrorInConsole(error);
    } finally {
        document.getElementById("submit-projection-search").classList.remove("is-loading");
    }

    // Trigger the default dataset collection to be selected in the
    if (datasetShareId) {
        selectDatasetCollection(null);  // Clear the label
        urlParamsPassed = true;
    } else if (layoutShareId) {
        datasetCollectionState.selectedShareId = layoutShareId;
        selectDatasetCollection(layoutShareId);
        urlParamsPassed = true;
    } else if (getCurrentUser()?.layout_share_id) {
        selectDatasetCollection(getCurrentUser().layout_share_id);
    }

    // Now, if URL params were passed and we have both patterns and a dataset collection,
    //  run the search
    if (urlParamsPassed) {

        if ((datasetShareId || datasetCollectionState.selectedShareId) && selectedPattern.shareId !== null && selectedPattern.selectedWeights.length > 0) {
            document.getElementById('submit-projection-search').click();
        }
    }

    // Add mutation observer to watch if #dropdown-dc-selector-label changes
    const observer = new MutationObserver((mutationsList, observer) => {
        for (const mutation of mutationsList) {
            if (mutation.type === 'childList') {
                // If the user selects a collection, clear the datasetShareId as scope has changed
                datasetShareId = null;
            }
        }
    });

    observer.observe(document.getElementById("dropdown-dc-selector-label"), { childList: true });
};
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

// Pre-initialize some stuff
await initCommonUI();

// If one of the "view genes" buttons is clicked, show the genes in a new window
for (const btn of document.getElementsByClassName("js-view-genes")) {
    btn.addEventListener('click', (event) => {
        let htmlStream = "<table>";
        // Gather genes and weights and show in new page
        for (const row of weightedGeneData) {
            htmlStream += `<tr><td>${row["gene"]}</td><td>${row["weight"]}</td></tr>`;
        }
        htmlStream += "</table>";
        const tab = window.open('about:blank', '_blank');
        tab.document.write(htmlStream);
        tab.document.close();
    });
}

// Change the svg scoring method when select element is changed
document.getElementById('svg-scoring-method').addEventListener('change', (event) => {
    if (isMulti) return;   // multi does not use this

    svgScoringMethod = event.target.value;

    // Loop through all tiles with svgData and update the display based on the selected method
    for (const tile of tilegrid.tiles) {
        if (tile.svg) {
            tile.updateSVGDisplay(svgScoringMethod);
        }
    }

});
