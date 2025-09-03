'use strict';

import { apiCallsMixin, createToast, getCurrentUser, getUrlParameter, initCommonUI, logErrorInConsole, registerPageSpecificLoginUIUpdates } from "./common.v2.js?v=2860b88";
import { datasetCollectionState, fetchDatasetCollections, selectDatasetCollection } from "../include/dataset-collection-selector/dataset-collection-selector.js?v=2860b88";
import { fetchGeneCartData, geneCollectionState, selectGeneLists } from "../include/gene-collection-selector/gene-collection-selector.js?v=2860b88";
import { TileGrid } from "./classes/tilegrid.js?v=2860b88";

let urlParamsPassed = false;
let isMultigene = false;
let tilegrid = null;
let svgScoringMethod = 'gene';

let datasetShareId = null;
let layoutShareId = null;
let shareUsed = false;

let currentlySelectedOrgId = "";
let annotationData = null;
let currentlySelectedGeneSymbol = null;  // goes to the plugin expression.js

/*
TODOs:
- Scrolling of datasets in collection should still show gene list
- When changing genes the tiles need to show 'loading' states before redrawing
- When changing scopes the tiles need to show 'loading' states before redrawing
*/

const organismSelector = document.getElementById('organism-selector');

// handle when the dropdown-gene-list-search-input input box is changed
document.getElementById('genes-manually-entered').addEventListener('change', (event) => {
    const searchTermString = event.target.value;
    updateGenesSelected(searchTermString);
});

document.getElementById('functional-annotation-toggle').addEventListener('click', (event) => {
    const annotationPanel = document.getElementById('extended-annotation-panel');
    const toggleIcon = document.querySelector('#functional-annotation-toggle i');
    const organismSelector = document.getElementById('annotation-panel-organism-selector-c');

    if (annotationPanel.classList.contains('is-hidden')) {
        annotationPanel.classList.remove('is-hidden');
        organismSelector.classList.remove('is-hidden');
        toggleIcon.classList.remove('mdi-chevron-down');
        toggleIcon.classList.add('mdi-chevron-up');
        return;
    }
    annotationPanel.classList.add('is-hidden');
    organismSelector.classList.add('is-hidden');
    toggleIcon.classList.remove('mdi-chevron-up');
    toggleIcon.classList.add('mdi-chevron-down');
});

// add event listener for when the submit-expression-search button is clicked
document.getElementById('submit-expression-search').addEventListener('click', async (event) => {
    const currentTarget = event.currentTarget;
    currentTarget.classList.add('is-loading');
    const isExactMatch = document.getElementById('gene-search-exact-match').checked;

    // Validate here. However for multigene fuzzy searches, we need the gene list to be populated
    if (!(isMultigene && !isExactMatch)) {
        const status = validateExpressionSearchForm();

        if (!status) {
            console.info("Aborting search");
            event.currentTarget.classList.remove('is-loading');
            return;
        }
    }

    document.getElementById("result-panel-initial-notification").classList.add('is-hidden');
    document.getElementById("result-panel-loader").classList.remove('is-hidden');

    // update multigene/single gene
    isMultigene = document.getElementById('single-multi-multi').checked;

    // if multigene, clear the selected gene symbol and hide the gene-result-list container
    document.getElementById("gene-result-list-c").classList.remove('is-hidden');
    document.getElementById("currently-selected-gene-header").classList.remove('is-hidden');
    document.getElementById("annotation-panel").classList.remove('is-hidden');
    document.getElementById("scoring-method-div").classList.remove('is-hidden');
    if (isMultigene) {
        currentlySelectedGeneSymbol = null;
        document.getElementById("gene-result-list-c").classList.add('is-hidden');
        document.getElementById("currently-selected-gene-header").classList.add('is-hidden');
        document.getElementById("annotation-panel").classList.add('is-hidden');
        document.getElementById("scoring-method-div").classList.add('is-hidden');
    }
    try {

        const setupTileGridFn = (datasetShareId) ? setupTileGrid(datasetShareId, "dataset") : setupTileGrid(datasetCollectionState.selectedShareId);

        const [annotRes, tilegridRes] = await Promise.allSettled([fetchGeneAnnotations(), setupTileGridFn]);
        tilegrid = tilegridRes.value;

        // auto-select the first gene in the list
        const firstGene = document.querySelector('.gene-result-list-item');
        if (!isMultigene && firstGene) {
            // Will call renderDisplays on the selected gene downstream
            firstGene.click();

        } else if (isMultigene) {
            let genes = Array.from(geneCollectionState.selectedGenes);
            if (!isExactMatch) {
                const status = validateExpressionSearchForm();

                if (!status) {
                    console.info("Aborting search");
                    event.currentTarget.classList.remove('is-loading');
                    return;
                }

                // If multigene and not exact match, render all the genes
                const allGenesElts = document.querySelectorAll('.gene-result-list-item');
                genes = Array.from(allGenesElts).map(elt => elt.textContent);
            }
            await tilegrid.renderDisplays(genes, isMultigene);

        }

        // If the user isn't logged in, set the first organism's annotation as the default
        if (!getCurrentUser().session_id && tilegrid?.datasets.length > 0) {
            const first_organism_id = tilegrid.datasets[0].organism_id;
            currentlySelectedOrgId = parseInt(first_organism_id);
            organismSelector.value = currentlySelectedOrgId;
            updateAnnotationDisplay();
        }

    } catch (error) {
        logErrorInConsole(error);
        return;
    } finally {
        currentTarget.classList.remove('is-loading');
        document.getElementById("result-panel-loader").classList.add('is-hidden');
    }

    const url = buildStateURL();

    // add to state history
    history.pushState(null, '', url);

});

const setDefaultOrganism = document.getElementById('set-default-organism');

// handle when the organism-selector select box is changed
organismSelector.addEventListener('change', (event) => {
    if (organismSelector?.value) {
        currentlySelectedOrgId = parseInt(organismSelector.value);
    } else {
        currentlySelectedOrgId = "";
    }

    if (currentlySelectedOrgId === "") {
        showOrganismSelectorTooltip();
        setDefaultOrganism.classList.add('is-hidden');
        return;
    }
    hideOrganismSelectorTooltip();
    updateAnnotationDisplay();

    // If the user is logged in and doesn't have a default org ID or it's different from their current,
    //  show the control
    if (!getCurrentUser().session_id) {
        return;
    }
    const shouldHide = getCurrentUser().default_org_id === currentlySelectedOrgId;
    setDefaultOrganism.classList.toggle('is-hidden', shouldHide);
});

setDefaultOrganism.addEventListener('click', (event) => {
    // we don't want to set to null, and the UI should have prevented this, but check just in case
    if (currentlySelectedOrgId === "") {
        return;
    }
    getCurrentUser().default_org_id = currentlySelectedOrgId;
    apiCallsMixin.saveUserDefaultOrgId(getCurrentUser());
    setDefaultOrganism.classList.add('is-hidden');
});

// Change the svg scoring method when select element is changed
document.getElementById('svg-scoring-method').addEventListener('change', (event) => {
    if (isMultigene) return;   // multigene does not use this

    svgScoringMethod = event.target.value;

    // Loop through all tiles with svgData and update the display based on the selected method
    for (const tile of tilegrid.tiles) {
        if (tile.svg) {
            tile.updateSVGDisplay(svgScoringMethod);
        }
    }
});

/**
 * Builds the URL for a GET request based on the current state.
 * @returns {string} The URL string.
 */
const buildStateURL = () => {

    // Create a new URL object (with no search params)
    const url = new URL('/expression.html', window.location.origin);

    // add the manually-entered genes

    const manuallyEnteredGenes = Array.from(geneCollectionState.manuallyEnteredGenes);
    if (manuallyEnteredGenes.length > 0) {
        url.searchParams.append('gene_symbol', manuallyEnteredGenes.join(','));
    }

    // are we doing exact matches?
    if (document.getElementById('gene-search-exact-match').checked) {
        url.searchParams.append('gene_symbol_exact_match', '1');
    }

    // get the value of the single-multi radio box
    const singleMulti = document.querySelector('input[name="single-multi"]:checked').value;
    url.searchParams.append('is_multigene', singleMulti === 'single' ? '0' : '1');

    // add the gene lists
    //  NOTE: This will only be for labeling purposes, since individual genes could have been
    //    deselected within
    if (geneCollectionState.selectedGeneLists.size > 0) {
        const geneCartShareIds = Array.from(geneCollectionState.selectedGeneLists);
        url.searchParams.append('gene_lists', geneCartShareIds.join(','));
    }

    // add the dataset collections
    if (datasetShareId) {
        url.searchParams.append('share_id', datasetShareId);
    } else if (datasetCollectionState.selectedShareId) {
        url.searchParams.append('layout_id', datasetCollectionState.selectedShareId);
    }

    return url.toString();
};

/**
 * Fetches gene annotations.
 * @param {Function} callback - The callback function to be executed after fetching gene annotations.
 * @returns {Promise<void>} - A promise that resolves when the gene annotations are fetched.
 */
const fetchGeneAnnotations = async (callback) => {
    try {
        annotationData = await apiCallsMixin.fetchGeneAnnotations(
            Array.from(geneCollectionState.selectedGenes).join(','),
            document.getElementById('gene-search-exact-match').checked,
            false,  // case-insensitive
            datasetCollectionState.selectedShareId,
            isMultigene
        );

        const geneResultCountElt = document.getElementById("gene-result-count");
        geneResultCountElt.innerHTML = Object.keys(annotationData).length;
        geneResultCountElt.parentElement.classList.remove('is-hidden');

        // Render template based on the number of annotations
        if (Object.keys(annotationData).length === 0) {
            document.getElementById('gene-result-list').innerHTML = '';

            const noHistoryTemplate = document.getElementById('tmpl-gene-result-none-found');
            document.getElementById('gene-result-list').appendChild(noHistoryTemplate.content.cloneNode(true));
        } else {
            const template = document.getElementById('tmpl-gene-result-item');
            document.getElementById('gene-result-list').innerHTML = '';

            for (const geneSymbol in annotationData) {
                const row = template.content.cloneNode(true);
                row.querySelector('li').innerHTML = geneSymbol;
                document.getElementById('gene-result-list').appendChild(row);

                // due to a python issue, at some point in depth the data becomes a string. Parse it.
                for (const organismId in annotationData[geneSymbol]['by_organism']) {
                    const annot = JSON.parse(annotationData[geneSymbol]['by_organism'][organismId][0]);
                    annotationData[geneSymbol]['by_organism'][organismId] = annot;
                }
            }

            // add event listeners to the gene result list items
            for (const geneResult of document.getElementsByClassName('gene-result-list-item')) {
                geneResult.addEventListener('click', (event) => {
                    const geneSymbol = event.target.textContent;
                    document.getElementById('currently-selected-gene').textContent = geneSymbol;

                    // remove is-selected from all the existing rows, then add it to this one
                    const rows = document.getElementsByClassName('gene-result-list-item');
                    for (const row of rows) {
                        row.classList.remove('is-selected');
                    }

                    event.currentTarget.classList.add('is-selected');
                    selectGeneResult(geneSymbol);
                });
            }
        }
    } catch (error) {
        console.error(error);
    }
};

/**
 * Fetches organisms and populates the organism selector dropdown.
 *
\ * @returns {Promise<void>} - A promise that resolves when the organisms are fetched and the dropdown is populated.
 */
const fetchOrganisms = async () => {
    try {
        const orgs = await apiCallsMixin.fetchOrganismList();
        const template = document.getElementById('tmpl-organism-option');

        for (const organism of orgs['organisms']) {
            const row = template.content.cloneNode(true);
            row.querySelector('option').innerHTML = organism.label;
            row.querySelector('option').value = organism.id;

            // if this matches the user's default organism, select it
            if (getCurrentUser().default_org_id === organism.id) {
                row.querySelector('option').selected = true;
                currentlySelectedOrgId = organism.id;
            }

            organismSelector.appendChild(row);
        }

    } catch (error) {
        createToast(`There was an error fetching the organism list (${error})`);
    }
};


/**
 * We need to maintain a few data structures with the selected and manually entered genes,
 * but these can be updated when the user types or via a share URL, so let's shift this
 * to a function.
 */
const updateGenesSelected = (searchTermString) => {
    const newManuallyEnteredGenes = searchTermString?.length > 0 ? new Set(searchTermString.split(/[ ,]+/)) : new Set();

    // Remove genes that have been deleted from the selectedGenes set
    for (const gene of geneCollectionState.manuallyEnteredGenes) {
        if (!newManuallyEnteredGenes.has(gene)) {
            geneCollectionState.selectedGenes.delete(gene);
        }
    }

    // Add new genes to the selectedGenes set
    for (const gene of newManuallyEnteredGenes) {
        geneCollectionState.selectedGenes.add(gene);
    }

    geneCollectionState.manuallyEnteredGenes = newManuallyEnteredGenes;

    //console.log("Selected genes updated:", Array.from(geneCollectionState.selectedGenes));
    //console.log("Manually entered genes updated:", Array.from(geneCollectionState.manuallyEnteredGenes));
};

/**
 * Hides the organism selector tooltip.
 */
const hideOrganismSelectorTooltip = () => {
    document.getElementById('organism-selector-control').classList.remove('has-tooltip-top', 'has-tooltip-arrow', 'has-tooltip-active');
    document.getElementById('organism-selector-control').removeAttribute('data-tooltip');
};

/**
 * Shows the organism selector tooltip.
 */
const showOrganismSelectorTooltip = () => {
    document.getElementById('organism-selector-control').setAttribute('data-tooltip', 'Select an organism to view annotation');
    document.getElementById('organism-selector-control').classList.add('has-tooltip-top', 'has-tooltip-arrow', 'has-tooltip-active');
};

/**
 * Parses the URL parameters related to gene cart and updates the UI accordingly.
 */
const parseGeneListURLParams = async () => {
    // handle manually-entered gene symbols
    const geneSymbols = getUrlParameter('gene_symbol');
    if (geneSymbols) {
        document.getElementById('genes-manually-entered').value = geneSymbols.replaceAll(',', ' ');
        // Update the geneCollectionState with the gene symbols from the URL
        geneCollectionState.selectedGenes.clear();
        const splitGeneSymbols = geneSymbols.split(',');
        for (const gs of splitGeneSymbols) {
            geneCollectionState.selectedGenes.add(gs);
        }
        geneCollectionState.manuallyEnteredGenes = geneCollectionState.selectedGenes;
        urlParamsPassed = true;
    }

    // handle passed gene lists
    let gene_lists = [];

    // This is a legacy option
    if (getUrlParameter('gene_cart_share_id')) {
        gene_lists.push(getUrlParameter('gene_cart_share_id'));
        urlParamsPassed = true;
    }

    await fetchGeneCartData();

    if (getUrlParameter('gene_lists')) {
        const geneLists = getUrlParameter('gene_lists').split(',');
        for (const geneList of geneLists) {

            const geneListMembers = await apiCallsMixin.fetchGeneCartMembers(geneList);
            const genes = [];

            for (const member of geneListMembers.gene_symbols) {
                genes.push(member.label);
            }

            geneCollectionState.shareIdsToGenes[geneList] = genes;
            geneCollectionState.selectedGenes = new Set([...geneCollectionState.selectedGenes, ...genes]);

        }
        selectGeneLists(geneLists); // declared in gene-collection-selector.js
        urlParamsPassed = true;
    }

    // are we doing exact matches?
    const exact_match = getUrlParameter('gene_symbol_exact_match');

    if (exact_match !== null && exact_match !== undefined) {
        document.getElementById('gene-search-exact-match').checked = exact_match === '1';
    } else {
        document.getElementById('gene-search-exact-match').checked = false;
    }

    // single or multiple gene view (convert to boolean)?
    const isMultigeneParam = getUrlParameter('is_multigene');
    isMultigene = isMultigeneParam === '1';
    if (isMultigene) {
        document.getElementById('single-multi-multi').checked = true;
    } else {
        document.getElementById('single-multi-single').checked = true;
    }
};

/**
 * Parses the URL parameters to handle the passed dataset collection.
 */
const parseDatasetCollectionURLParams = () => {
    // handle passed dataset collection
    const layoutShareId = getUrlParameter('layout_id');

    if (!layoutShareId) {
        return;
    }

    datasetCollectionState.selectedShareId = layoutShareId;
    datasetCollectionState.selectedLabel = datasetCollectionState.labelIndex[layoutShareId];
    document.getElementById('dropdown-dc-selector-label').innerHTML = datasetCollectionState.selectedLabel;
};

/**
 * Selects a gene result and performs necessary actions based on the selected gene symbol.
 *
 * @param {string} geneSymbol - The symbol of the gene to be selected.
 */
const selectGeneResult = async (geneSymbol) => {
    const selectedOrganismId = organismSelector.value;
    currentlySelectedGeneSymbol = geneSymbol;

    // if no organism is selected, display a tooltip to choose one
    if (selectedOrganismId === "") {
        showOrganismSelectorTooltip();
    } else {
        updateAnnotationDisplay();
    }

    // Other things can be called next, such as plotting calls
    if (tilegrid) {
        // Revert back to "#result-panel-grid" display before rendering the new gene displays
        document.getElementById("result-panel-grid").classList.remove("is-hidden");
        document.getElementById("zoomed-panel-grid").classList.add("is-hidden");

        await tilegrid.renderDisplays(currentlySelectedGeneSymbol, isMultigene, svgScoringMethod);
    }

    // call any callbacks that have been added (usually by plugins)
    if (window.geneChangeCallbacks) {
        for (const callback of window.geneChangeCallbacks) {
            callback(currentlySelectedGeneSymbol);
        }
    }

};

/**
 * Sets up the tile grid with the provided shareId and type.
 *
 * @param {string} shareId - The shareId for the tile grid.
 * @param {string} [type="layout"] - The type of the tile grid. Default is "layout".
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

        tilegrid.applyTileGrid(isMultigene);

    } catch (error) {
        logErrorInConsole(error);
    } finally {
        return tilegrid;
    }
};

/**
 * Updates the annotation display based on the currently selected gene symbol and organism ID.
 */
const updateAnnotationDisplay = () => {
    // these make some of the syntax below shorter
    const gs = currentlySelectedGeneSymbol;
    const oid = currentlySelectedOrgId;

    // clear the external resource links and GO terms
    document.getElementById('external-resource-links').innerHTML = '';
    document.getElementById('go-terms').innerHTML = '';
    document.getElementById('go-term-count').innerHTML = '';

    // Did we find annotation for this gene symbol at all?
    if (!annotationData.hasOwnProperty(gs)) {
        document.getElementById('currently-selected-gene').innerHTML = "";
        document.getElementById('currently-selected-gene-product').innerHTML = "(annotation not available for this gene)";
        document.getElementById('currently-selected-gene-product').classList.remove('is-hidden');

        const dbxref_template = document.querySelector('#tmpl-external-resource-link-none-found');
        const dbxref_template_row = dbxref_template.content.cloneNode(true);
        document.getElementById('external-resource-links').appendChild(dbxref_template_row);
        return;
    }

    // if the selected organism is not in the annotation data, show a message
    if (!annotationData[gs]['by_organism'].hasOwnProperty(oid)) {
        document.getElementById('currently-selected-gene-product').innerHTML = " - (annotation not available for this organism)";
        document.getElementById('currently-selected-gene-product').classList.remove('is-hidden');

        const dbxref_template = document.querySelector('#tmpl-external-resource-link-none-found');
        const dbxref_template_row = dbxref_template.content.cloneNode(true);
        document.getElementById('external-resource-links').appendChild(dbxref_template_row);
        return;
    }

    // if we got this far, we have annotation for this one. let's display it
    const annotation = annotationData[gs]['by_organism'][oid];
    document.getElementById('annotation-panel-gene-symbol').innerHTML = gs;

    // Gene product
    document.getElementById('currently-selected-gene-product').innerHTML = " - " + annotation['product'];
    document.getElementById('currently-selected-gene-product').classList.remove('is-hidden');
    document.getElementById('annotation-panel-gene-product').innerHTML = annotation['product'];

    // aliases and Ensembl ID
    document.getElementById('annotation-panel-gene-aliases').innerHTML = annotation['aliases'].length > 1 ? annotation['aliases'].join(', ') : "None found";

    document.getElementById('annotation-panel-gene-ensembl-release').innerHTML = annotation['ensembl_release'];
    document.getElementById('annotation-panel-gene-ensembl-id').innerHTML = annotation['ensembl_id'];

    const ensembl_url = "https://www.ensembl.org/Multi/Search/Results?q=" + annotation['ensembl_id'] + ";site=ensembl";
    document.getElementById('annotation-panel-gene-ensembl-id').setAttribute('href', ensembl_url);

    // External database references
    let good_dbxref_count = 0;

    for (const dbxref of annotation['dbxrefs']) {
        if (dbxref['url'] !== null) {
            const dbxref_template = document.getElementById('tmpl-external-resource-link');
            const row = dbxref_template.content.cloneNode(true);
            row.querySelector('a').innerHTML = dbxref['source'];
            row.querySelector('a').href = dbxref['url'];
            document.getElementById('external-resource-links').appendChild(row);
            good_dbxref_count++;
        }
    }

    if (good_dbxref_count === 0) {
        const dbxref_template = document.getElementById('tmpl-external-resource-link-none-found');
        const row = dbxref_template.content.cloneNode(true);
        document.getElementById('external-resource-links').appendChild(row);
    }

    // GO terms
    document.getElementById('go-term-count').innerHTML = `(${annotation.go_terms.length})`;
    if (annotation.go_terms.length === 0) {
        const go_term_template = document.getElementById('tmpl-go-term-none-found');
        const row = go_term_template.content.cloneNode(true);
        document.getElementById('go-terms').appendChild(row);
        return;
    }
    for (const go_term of annotation['go_terms']) {
        const go_term_template = document.getElementById('tmpl-go-term');
        const go_term_url = "https://amigo.geneontology.org/amigo/search/ontology?q=" + go_term['go_id'];

        const row = go_term_template.content.cloneNode(true);
        row.querySelector('.go-term-id').innerHTML = go_term['go_id'];
        row.querySelector('.go-term-id').href = go_term_url;
        row.querySelector('.go-term-label').innerHTML = go_term['name'];
        document.querySelector('#go-terms').appendChild(row);
    }
};

/**
 * Validates the expression search form.
 *
 * @returns {boolean} Returns true if the form is valid, false otherwise.
 */
const validateExpressionSearchForm = () => {

    // User passed in a single dataset share ID.
    if (datasetShareId) {
        return true;
    }

    // User must have either selected a gene list or entered genes manually. Either of these
    // will populate the geneCollectionState.selectedGenes array
    if (geneCollectionState.selectedGenes.size + geneCollectionState.manuallyEnteredGenes.size === 0) {
        createToast('Please enter at least one gene to proceed');
        return false;
    }

    // Check if the user has selected any dataset collections
    if (!datasetCollectionState.selectedShareId) {
        createToast('Please select at least one dataset to proceed');
        return false;
    }

    const isMulti = document.getElementById('single-multi-multi').checked;
    const isExactMatch = document.getElementById('gene-search-exact-match').checked;

    // If multi, check that at least two genes are selected
    if (isMulti) {
        // if not exact match, test if there are at least two genes in results
        if (!isExactMatch) {
            const allGenesElts = document.querySelectorAll('.gene-result-list-item');
            const allGenes = Array.from(allGenesElts).map(elt => elt.textContent);
            // Add in selected genes
            const comboGenes = allGenes.concat(Array.from(geneCollectionState.selectedGenes));
            if (comboGenes.length < 2) {
                createToast('Need at least two genes to proceed');
                return false;
            }

        } else if (geneCollectionState.selectedGenes.size < 2) {
            createToast('Please select at least two genes to proceed');
            return false;
        }
    }

    return true;
};

/**
 * Handles the UI updates specific to the page after login.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when all UI updates are completed.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {
    // Set the page header title
    document.getElementById('page-header-label').textContent = 'Gene Expression Search';

    datasetShareId = getUrlParameter('share_id');
    layoutShareId = getUrlParameter('layout_id');
    const cartShareId = getUrlParameter('gene_lists');
    const geneSymbolId = getUrlParameter('gene_symbol');
    shareUsed = getUrlParameter('share_used') === '1';

    // Wait until all pending API calls have completed before checking if we need to search
    document.getElementById("submit-expression-search").classList.add("is-loading");
    try {
        // SAdkins note - Promise.all fails fast,
        // but Promise.allSettled waits until all resolve/reject and lets you know which ones failed
        const [dc_result, org_result] = await Promise.all([
            parseGeneListURLParams(),
            fetchDatasetCollections(layoutShareId),
            fetchOrganisms()
        ]);
    } catch (error) {
        logErrorInConsole(error);
    } finally {
        document.getElementById("submit-expression-search").classList.remove("is-loading");
    }

    parseDatasetCollectionURLParams();

    // Trigger the default dataset collection to be selected in the
    if (datasetShareId) {
        selectDatasetCollection(null);  // Clear the label
        urlParamsPassed = true;
    } else if (layoutShareId) {
        datasetCollectionState.selectedShareId = layoutShareId;
        selectDatasetCollection(layoutShareId);
        urlParamsPassed = true;
    } else if (getCurrentUser().layout_share_id) {
        selectDatasetCollection(getCurrentUser().layout_share_id);
    }

    // Now, if URL params were passed and we have both genes and a dataset collection,
    //  run the search
    if (urlParamsPassed) {
        if ((datasetShareId || datasetCollectionState.selectedShareId) && geneCollectionState.selectedGenes.size > 0) {
            document.getElementById('submit-expression-search').click();
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

    // If we entered via a share link, conditionally display the notification
    if (shareUsed) {
        const shareData = await apiCallsMixin.fetchShareData(datasetShareId, layoutShareId);
        document.getElementById("share-entrance-notification").classList.remove('is-hidden');

        // if an individual dataset was shared, show its info
        if (datasetShareId) {
            document.getElementById("share-entrance-dataset-label").innerHTML = shareData['dataset_label'];
            document.getElementById("share-entrance-dataset").classList.remove('is-hidden');

            if (shareData['owner_name']) {
                document.getElementById("share-entrance-dataset-owner-label").innerHTML = shareData['owner_name'];
            } else {
                document.getElementById("share-entrance-dataset-owner-label").innerHTML = 'Unknown';
            }

            // if a dataset collection was shared, show its info
        } else if (layoutShareId) {
            document.getElementById("share-entrance-layout-label").innerHTML = shareData['layout_label'];
            document.getElementById("share-entrance-layout").classList.remove('is-hidden');

            if (shareData['owner_name']) {
                document.getElementById("share-entrance-layout-owner-label").innerHTML = shareData['owner_name'];
            } else {
                document.getElementById("share-entrance-layout-owner-label").innerHTML = 'Unknown';
            }
        }

        if (cartShareId || geneSymbolId) {
            document.getElementById("share-entrance-genes-preselected").classList.remove('is-hidden');
        } else {
            if (shareData['gene_symbol']) {
                document.getElementById("share-entrance-genes-autoselected").classList.remove('is-hidden');
                document.getElementById("share-entrance-genes-label").innerHTML = shareData['gene_symbol'];
                document.getElementById("share-entrance-genes").classList.remove('is-hidden');

                // insert the gene symbol and trigger a search
                document.getElementById('genes-manually-entered').value = shareData['gene_symbol'];
                updateGenesSelected(shareData['gene_symbol']);
                document.getElementById('submit-expression-search').click();
            } else {
                document.getElementById("share-entrance-genes-noneselected").classList.remove('is-hidden');
            }
        }
    }
};
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

// Pre-initialize some stuff
await initCommonUI();
