'use strict';

let urlParamsPassed = false;
let isMulti = false;
let tilegrid = null;
let svgScoringMethod = 'gene';
let projectionOpts = {patternSource: null, algorithm: null, gctype: null};
let weightedGeneData = null;
let datasetShareId = null;
let layoutShareId = null;

// imported from pattern-collection-selector.js
// selectedPattern = {shareId: null, label: null, gctype: null, selectedWeights: []};

// imported from dataset-collection-selector.js
// selected_dc_share_id = null;
// selected_dc_label = null;

// Proxy the selectedPattern object to watch for changes to the selectedWeights array
selectedPattern = new Proxy(selectedPattern, {
    set: (target, key, value) => {
        target[key] = value;
        const algorithmElt = document.getElementById('algorithm');

        if (key === "selectedWeights") {
            document.getElementById("single-multi-multi").disabled = false;
            if(value.length < 2) {
                document.getElementById("single-multi-multi").disabled = true;
                document.getElementById("single-multi-single").checked = true;
                isMulti = false;
            }
            // Enable or disable the "binary" option if all selectedWeights have "binary" property set to True
            const binary = value.every((w) => w.binary);
            algorithmElt.querySelector('option[value="binary"]').disabled  = !binary;

        } else if (key === "gctype") {
            // Adjust algorithm options based on gctype
            algorithmElt.querySelector('option[value="nmf"]').disabled = false;
            algorithmElt.querySelector('option[value="fixednmf"]').disabled = false;

            if (value === "unweighted-list") {
                algorithmElt.querySelector('option[value="nmf"]').disabled = true;
                algorithmElt.querySelector('option[value="fixednmf"]').disabled = true;
            }
        }
        return true;
    }
});

/**
 * Handles the UI updates specific to the page login.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are completed.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {

    // Set the page header title
    document.getElementById('page-header-label').textContent = 'Projection Search';

    // Set current sidebar menu item to active
	for (const elt of document.querySelectorAll("#primary-nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}

	document.querySelector("a[tool='projection'").classList.add("is-active");

    datasetShareId = getUrlParameter('share_id');
    layoutShareId = getUrlParameter('layout_id');

    // add event listener for when the submit-projection-search button is clicked
    document.querySelector('#submit-projection-search').addEventListener('click', async (event) => {

        const status = validateProjectionSearchForm();

        if (! status) {
            console.info("Aborting search");
            return;
        }

        // update multi/single pattern
        isMulti = document.querySelector('#single-multi-multi').checked;

        populatePatternResultsList();

        // if multi, clear the selected pattern symbol and hide the pattern-result-list container
        document.getElementById("pattern-result-list-c").classList.remove('is-hidden');
        document.getElementById("scoring-method-div").classList.remove('is-hidden');
        if (isMulti) {
            document.getElementById("pattern-result-list-c").classList.add('is-hidden');
            document.getElementById("scoring-method-div").classList.add('is-hidden');
        }

        try {
            const setupTileGridFn = (datasetShareId) ? setupTileGrid(datasetShareId, "dataset") : setupTileGrid(selected_dc_share_id);
            tilegrid =  await setupTileGridFn;

            // auto-select the first pattern in the list
            const first_pattern = document.querySelector('.pattern-result-list-item');
            if (!isMulti && first_pattern) {

                first_pattern.click();
            }

        } catch (error) {
            logErrorInConsole(error);
            return;
        }

        const url = buildStateUrl();
        // add to state history
        history.pushState(null, '', url);

    });

    // Change the svg scoring method when select element is changed
    document.getElementById('svg-scoring-method').addEventListener('change', (event) => {
        if (isMulti) return;   // multi does not use this

        svgScoringMethod = event.target.value;
        // Get pattern symbol from currently selected list item
        let listItem = document.querySelector('.pattern-result-list-item.is-selected');
        if (!listItem) {
            listItem = document.querySelector('.pattern-result-list-item');
        }

        const pattern = listItem.textContent;
        selectPatternWeightResult(pattern);
    });

    // Wait until all pending API calls have completed before checking if we need to search
    try {
        // SAdkins note - Promise.all fails fast,
        // but Promise.allSettled waits until all resolve/reject and lets you know which ones failed
        const [cartResult, dcResult,] = await Promise.all([
            fetchPatternsData(parsepatternCartURLParams),
            fetchDatasetCollections(parseDatasetCollectionURLParams),
        ]);

        // Should help with lining things up on index page
        document.getElementById("dropdown-dc").classList.remove("is-right");

    } catch (error) {
        logErrorInConsole(error);
    }

    // Trigger the default dataset collection to be selected in the
    if (datasetShareId) {
        selectDatasetCollection(null);  // Clear the label
        urlParamsPassed = true;
    } else if (layoutShareId) {
        selected_dc_share_id = layoutShareId;
        selectDatasetCollection(layoutShareId);
        urlParamsPassed = true;
    } else if (!layoutShareId && CURRENT_USER.default_profile_share_id) {
        selectDatasetCollection(CURRENT_USER.default_profile_share_id);
    }

    // Now, if URL params were passed and we have both patterns and a dataset collection,
    //  run the search
    if (urlParamsPassed) {
        if ((datasetShareId || selected_dc_share_id) && selectedPattern.shareId !== null && selectedPattern.selectedWeights.length > 0) {
            document.querySelector('#submit-projection-search').click();
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
}

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

    // Add the multipattern_plots value to the URL
    const multipatternPlots = document.querySelector('#single-multi-multi').checked ? 1 : 0;
    url.searchParams.set('multipattern_plots', multipatternPlots);

    // Add the pattern source to the URL
    url.searchParams.set('projection_source', selectedPattern.shareId);

    // Add the dataset collection to the URL
    if (datasetShareId) {
        url.searchParams.append('share_id', datasetShareId);
    } else if (selected_dc_share_id) {
        url.searchParams.append('layout_id', selected_dc_share_id);
    }

    // Add the selected pattern weights to the URL
    const weights = selectedPattern.selectedWeights.map((w) => w.label);
    url.searchParams.set('projection_patterns', weights.join(','));

    return url.toString();
}

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
}

/**
 * Populates the pattern results list with weights.
 */
const populatePatternResultsList = () => {
    const template = document.querySelector('#tmpl-pattern-result-item');
    document.querySelector('#pattern-result-list').innerHTML = '';

    // sort the selectedWeights array based on the numeric value at the end of each string
    const sortedLabels = selectedPattern.selectedWeights.map((weight) => weight.label).sort(customNumericSort);

    for (const label of sortedLabels) {
        const row = template.content.cloneNode(true);
        row.querySelector('li').innerHTML = label;
        row.querySelector('li').dataset.weight = label;
        document.querySelector('#pattern-result-list').appendChild(row);

        const thisRow = document.querySelector(`.pattern-result-list-item[data-weight="${label}"]`);
        thisRow.addEventListener('click', (event) => {

            // remove is-selected from all the existing rows, then add it to this one
            const rows = document.querySelectorAll('.pattern-result-list-item');
            for (const row of rows) {
                row.classList.remove('is-selected');
            }

            event.currentTarget.classList.add('is-selected');
            selectPatternWeightResult(label);
        });
    }
}

/**
 * Parses the URL parameters and updates the UI based on the values.
 */
const parsepatternCartURLParams = () => {
    // if projection algorithm is passed, set it in #algorithm
    const projectionAlgorithm = getUrlParameter('projection_algorithm');
    if (projectionAlgorithm) {
        document.getElementById('algorithm').value = projectionAlgorithm;
    }

    // single or multiple pattern view (convert to boolean)?
    // NOTE: This will be adjusted if the pattern only has one weight
    const isMultiParam = getUrlParameter('multipattern_plots');
    isMulti = isMultiParam === '1';
    if (isMulti) {
        document.querySelector('#single-multi-multi').checked = true;
    } else {
        document.querySelector('#single-multi-single').checked = true;
    }

    // handle passed pattern lists
    const pattern = getUrlParameter('projection_source')
    if (pattern) {
        urlParamsPassed = true;
        const foundPattern = flatPatternsCartData.find((p) => p.share_id === pattern);
        selectedPattern = {shareId: foundPattern.share_id, label: foundPattern.label, gctype: foundPattern.gctype, selectedWeights: []};
        updatePatternListSelectorLabel()
    }

    // handle manually-entered pattern symbols
    const urlWeights = getUrlParameter('projection_patterns');
    if (pattern && urlWeights) {
        // Cannot have weights without a source pattern
        const labels = urlWeights.split(',');
        selectedPattern.selectedWeights = labels.map((label) => ({label, top_up: null, top_down: null}));
    }
}

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

    selected_dc_share_id = layoutShareId;
    selected_dc_label = dataset_collection_label_index[layoutShareId];
    document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;
}

/**
 * Selects a pattern weight and performs various actions based on the selected weight.
 * @param {string} label - The selected weight label.
 */
const selectPatternWeightResult = async (label) => {

    // get the selected pattern object
    const obj = selectedPattern.selectedWeights.find((w) => w.label === label);

    // if projection algorithm is "nmf", then hide the top_down genes
    const projectionAlgorithm = document.getElementById('algorithm').value;
    if (projectionAlgorithm === 'nmf') {
        document.getElementById('svg-scoring-method').value = 'top_up';
        document.getElementById('top-down-genes').classList.add('is-hidden');
    } else {
        document.getElementById('top-down-genes').classList.remove('is-hidden');
    }

    // if isMulti=false, show top_up and top_down genes
    document.getElementById("top-genes-c").classList.remove('is-hidden');
    if (isMulti || selectedPattern.gctype === "unweighted-list") {
        document.getElementById("top-genes-c").classList.add('is-hidden');
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

    document.getElementById("btn-view-weighted-genes").classList.remove("is-hidden");
    try {
        const data = await apiCallsMixin.fetchPatternWeightedGenes(selectedPattern.shareId, label);
        weightedGeneData = data;
    } catch (error) {
        logErrorInConsole(error);
        document.getElementById("btn-view-weighted-genes").classList.add("is-hidden");
    }

    // Other things can be called next, such as plotting calls
    if (tilegrid) {
        tilegrid.renderDisplays(label, isMulti, svgScoringMethod, projectionOpts);
    }
}

/**
 * Sets up the tile grid for projection.
 *
 * @param {string} shareId - The share ID of the layout.
 * @returns {Promise<TileGrid>} - A promise that resolves to the initialized TileGrid object.
 */
const setupTileGrid = async (shareId, type="layout") => {

    // Cannot proceed without a shareId
    if (!shareId) {
        return;
    }

    const tilegrid = new TileGrid(shareId, type, "#result-panel-grid");
    try {
        tilegrid.layout = await tilegrid.getLayout();
        await tilegrid.addAllDisplays();

        tilegrid.applyTileGrid(isMulti);

        const algorithm = document.getElementById('algorithm').value;

        // create projectionOpts object out of selectedPattern.shareId, algorithm, and selectedPattern.gctype
        projectionOpts = {
            patternSource: selectedPattern.shareId,
            algorithm,
            gctype: selectedPattern.gctype
        };

        for (const tile of tilegrid.tiles) {
            tile.enableProjectR();
        }

        await tilegrid.addDefaultDisplays();

        // NOTE - the tilegrid.renderDisplays() call below can check and use the first array element of the selected_genes array if single_pattern
        // We do not render for single-gene searches because the first pattern result is "clicked" and the tilegrid is rendered in the event listener.

        if (isMulti && selectedPattern.selectedWeights.length) {
            // create array of selected weight labels
            const selectedWeights = Array.from(selectedPattern.selectedWeights).map((w) => w.label);
            await tilegrid.renderDisplays(selectedWeights, isMulti, svgScoringMethod, projectionOpts);
        }
    } catch (error) {
        logErrorInConsole(error);
    } finally {
        return tilegrid;
    }
}

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
    if (!selected_dc_share_id) {
        createToast('Please select at least one dataset to proceed');
        document.querySelector('#dropdown-dc button').classList.add('is-danger');
        return false;
    }

    return true;
}

document.getElementById('btn-view-weighted-genes').addEventListener('click', (event) => {
    let htmlStream = "<table>";
    // Gather genes and weights and show in new page
    for (const row of weightedGeneData) {
        htmlStream += `<tr><td>${row["gene"]}</td><td>${row["weight"]}</td></tr>`;
    }
    htmlStream += "</table>"
    const tab = window.open('about:blank', '_blank');
    tab.document.write(htmlStream);
    tab.document.close();
});
