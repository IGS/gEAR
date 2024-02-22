'use strict';

let urlParamsPassed = false;
let currentlySelectedPattern = null;
let currentlySelectedOrgId = "";
let isMulti = false;
let tilegrid = null;
let svgScoringMethod = 'gene';

const handlePageSpecificLoginUIUpdates = async (event) => {

    // Set the page header title
    document.getElementById('page-header-label').textContent = 'Projection Search';

    // Set current sidebar menu item to active
	for (const elt of document.querySelectorAll("#primary_nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}


    // add event listener for when the submit-projection-search button is clicked
    document.querySelector('#submit-projection-search').addEventListener('click', async (event) => {
        const status = validateProjectionSearchForm();

        if (! status) {
            console.log("Aborting search");
            return;
        }

        // update multi/single pattern
        isMulti = document.querySelector('#single-multi-multi').checked;

        // if multi, clear the selected pattern symbol and hide the pattern-result-list container
        document.getElementById("pattern-result-list-c").classList.remove('is-hidden');
        document.getElementById("currently-selected-pattern-header").classList.remove('is-hidden');
        document.getElementById("scoring-method-div").classList.remove('is-hidden');
        if (isMulti) {
            currentlySelectedPattern = null;
            document.getElementById("pattern-result-list-c").classList.add('is-hidden');
            document.getElementById("currently-selected-pattern-header").classList.add('is-hidden');
            document.getElementById("scoring-method-div").classList.add('is-hidden');
        }


        try {
            const [tilegridRes] = await setupTileGrid(selected_dc_share_id);
            tilegrid = tilegridRes.value;

            // auto-select the first pattern in the list
            const first_pattern = document.querySelector('.pattern-result-list-item');
            if (!isMulti && first_pattern) {
                first_pattern.click();
            }

        } catch (error) {
            logErrorInConsole(error);
        }
    });

    // Add event listeners to the pattern result list items even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (!event.target.classList.contains('pattern-result-list-item')) {
            return;
        }

        const pattern_symbol = event.target.textContent;
        document.querySelector('#currently-selected-pattern').innerHTML = pattern_symbol;

        // remove is-selected from all the existing rows, then add it to this one
        const rows = document.querySelectorAll('.pattern-result-list-item');
        rows.forEach((row) => {
            row.classList.remove('is-selected');
        });

        event.target.classList.add('is-selected');
        selectPatternResult(pattern_symbol);
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
        selectPatternResult(pattern);
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

    // Now, if URL params were passed and we have both patterns and a dataset collection,
    //  run the search
    if (urlParamsPassed) {
        if (selected_dc_share_id && selectedPatterns.size > 0) {
            document.querySelector('#submit-projection-search').click();
        }
    }
}


const parsepatternCartURLParams = () => {
    // handle manually-entered pattern symbols
    const urlPatterns = getUrlParameter('patterns');
    if (urlPatterns) {
        selectedPatterns = new Set(urlPatterns.split(','));
        urlParamsPassed = true;
    }

    // handle passed pattern lists
    const pattern = getUrlParameter('projection_source')
    if (pattern) {
        urlParamsPassed = true;
        selectPatternLists(pattern); // declared in pattern-collection-selector.js
    }

    // single or multiple pattern view (convert to boolean)?
    const isMultiParam = getUrlParameter('multipattern_plots');
    isMulti = isMultiParam === '1';
    if (isMulti) {
        document.querySelector('#single-multi-multi').checked = true;
    } else {
        document.querySelector('#single-multi-single').checked = true;
    }
}

const parseDatasetCollectionURLParams = async () => {
    // handle passed dataset collection
    const layoutShareId = getUrlParameter('layout_id');

    if (!layoutShareId) {
        return;
    }

    selected_dc_share_id = layoutShareId;
    selected_dc_label = dataset_collection_label_index[layoutShareId];
    document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;
}

const selectPatternResult = (pattern_symbol) => {
    currentlySelectedPattern = pattern_symbol;

    // Other things can be called next, such as plotting calls
    if (tilegrid) {
        tilegrid.renderDisplays(currentlySelectedPattern, isMulti, svgScoringMethod);
    }
}

const setupTileGrid = async (layout_share_id) => {
    const tilegrid = new TileGrid(layout_share_id, "#result-panel-grid");
    try {
        tilegrid.layout = await tilegrid.getLayout();
        await tilegrid.addAllDisplays();

        tilegrid.patternrateTileGrid(isMulti);
        tilegrid.applyTileGrid(isMulti);

        await tilegrid.performProjection(); // TODO:

        await tilegrid.addDefaultDisplays();

        // NOTE - the tilegrid.renderDisplays() call below can check and use the first array element of the selected_patterns array if single_pattern
        // But we are using a string for clarity.
        if (isMulti) {
            // Don't render yet if a pattern is not selected
            if (selected_patterns.size) {
                await tilegrid.renderDisplays(Array.from(selectedPatterns), isMulti);
            }
        } else {
            // Don't render yet if a pattern is not selected
            if (currentlySelectedPattern) {
                await tilegrid.renderDisplays(currentlySelectedPattern, isMulti, svgScoringMethod);
            }
        }
    } catch (error) {
        logErrorInConsole(error);
    } finally {
        return tilegrid;
    }
}

const validateProjectionSearchForm = () => {
    // User must have either selected a pattern list or entered patterns manually. Either of these
    // will populate the selected_patterns array
    if (selectedPatterns.size) {
        createToast('Please enter at least one pattern to proceed');
        return false;
    }

    // Check if the user has selected any dataset collections
    if (!selected_dc_share_id) {
        createToast('Please select at least one dataset to proceed');
        return false;
    }

    return true;
}
