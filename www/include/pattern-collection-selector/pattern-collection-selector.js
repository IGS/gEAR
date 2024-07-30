"use strict";

// NOTE: This component depends on common.js and on Bulma CSS being imported in the parent HTML file

// Terminology:
// - patterns - The pattern source data
// - weights - The individual patterns

let patternsCartData = null;
let flatPatternsCartData = null;
let selectedPattern = {shareId: null, label: null, gctype: null, organismId: null, selectedWeights: []}; // This is used by the script that includes this file

// Add event listener to dropdown trigger
document.querySelector("#dropdown-pattern-lists > button.dropdown-trigger").addEventListener("click", (event) => {
    const item = event.currentTarget;
    item.closest(".dropdown").classList.toggle('is-active');
});

// Add event listeners to the pattern list category selectors
const categorySelectors = document.querySelectorAll('#dropdown-content-pattern-list-category .ul-li');
categorySelectors.forEach((element) => {
    element.addEventListener('click', (event) => {
        const category = event.target.dataset.category;
        setActivePatternCartCategory(category);

        categorySelectors.forEach((element) => {
            element.classList.remove('is-selected');
            element.classList.add('is-clickable');
        });

        event.target.classList.add('is-selected');
        event.target.classList.remove('is-clickable');
    });
});

document.getElementById('dropdown-pattern-list-proceed').addEventListener('click', (event) => {
    updatePatternListSelectorLabel();

    // close the dropdown
    document.getElementById('dropdown-pattern-lists').classList.remove('is-active');
});

// Add a click listener to the dropdown-pattern-list-cancel button
document.getElementById('dropdown-pattern-list-cancel').addEventListener('click', (event) => {
    // clear pattern lists and pattern list areas
    document.getElementById('dropdown-content-pattern-lists').innerHTML = '';
    document.getElementById('dropdown-content-weights').innerHTML = '';

    const categorySelectors = document.querySelectorAll('#dropdown-content-pattern-list-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    // clear the patterns-manually-entered input element
    document.getElementById('dropdown-pattern-list-search-input').value = '';
    document.getElementById('dropdown-pattern-list-selector-label').innerHTML = 'Quick search using pattern lists';

    // and finally the related pattern lists and patterns
    // Doing it this way to not break the Proxy object on projection.js
    selectedPattern.shareId = null;
    selectedPattern.label = null;
    selectedPattern.gctype = null;
    selectedPattern.organismId = null;
    selectedPattern.selectedWeights = [];
});

// Monitor key strokes after user types more than 2 characters in the dropdown-pattern-list-search-input box
document.getElementById('dropdown-pattern-list-search-input').addEventListener('keyup', (event) => {
    const search_term = event.target.value;

    if (search_term.length === 0) {
        // clear the gene list
        document.getElementById('dropdown-content-pattern-lists').innerHTML = '';
        document.getElementById('dropdown-content-weights').innerHTML = '';
        return;
    } else if (search_term.length <= 2) {
        return;
    }

    const categorySelectors = document.querySelectorAll('#dropdown-content-pattern-list-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    document.getElementById('dropdown-content-pattern-lists').innerHTML = '';
    const patternListItemTemplate = document.getElementById('tmpl-pattern-list-item');

    for (const cartType in patternsCartData) {
        for (const cart of patternsCartData[cartType]) {
            if (cart.label.toLowerCase().includes(search_term.toLowerCase())) {
                const row = patternListItemTemplate.content.cloneNode(true);
                createPatternListItem(row, cart);
            }
        }
    }
});

document.getElementById("dropdown-pattern-list-clear-weights").addEventListener("click", (event) => {
    // uncheck all the existing rows
    const rows = document.getElementsByClassName('dropdown-weight-item');

    if (rows[0].classList.contains('is-disabled')) return;  // If there is only one weight, we can't clear it

    for (const row of rows ) {
        row.classList.remove("is-selected");
        row.querySelector(".mdi").classList.replace("mdi-check", "mdi-plus");
    };

    selectedPattern.selectedWeights = [];
});

document.getElementById("dropdown-pattern-list-top5-weights").addEventListener("click", (event) => {
    // Get the labels of the first 5 weights and select them
    const rows = document.getElementsByClassName('dropdown-weight-item');
    const labels = Array.from(rows).slice(0, 5).map((row) => row.dataset.label);
    selectedPattern.selectedWeights = [];

    // Some labels may already be selected, so we need to clear them first
    for (const row of rows ) {
        row.classList.remove("is-selected");
        row.querySelector(".mdi").classList.replace("mdi-check", "mdi-plus");
    };

    selectPatternWeights(labels);

});
document.getElementById("dropdown-pattern-list-top10-weights").addEventListener("click", (event) => {
    // Get the labels of the first 10 weights and select them
    const rows = document.getElementsByClassName('dropdown-weight-item');
    const labels = Array.from(rows).slice(0, 10).map((row) => row.dataset.label);
    selectedPattern.selectedWeights = [];

    for (const row of rows ) {
        row.classList.remove("is-selected");
        row.querySelector(".mdi").classList.replace("mdi-check", "mdi-plus");
    };

    selectPatternWeights(labels);
});

document.getElementById("dropdown-pattern-list-top20-weights").addEventListener("click", (event) => {
    // Get the labels of the first 20 weights and select them
    const rows = document.getElementsByClassName('dropdown-weight-item');
    const labels = Array.from(rows).slice(0, 20).map((row) => row.dataset.label);
    selectedPattern.selectedWeights = [];

    for (const row of rows ) {
        row.classList.remove("is-selected");
        row.querySelector(".mdi").classList.replace("mdi-check", "mdi-plus");
    };


    selectPatternWeights(labels);
});
document.getElementById("dropdown-pattern-list-all-weights").addEventListener("click", (event) => {
    const rows = document.getElementsByClassName('dropdown-weight-item');
    const labels = Array.from(rows).map((row) => row.dataset.label);
    selectedPattern.selectedWeights = [];

    for (const row of rows ) {
        row.classList.remove("is-selected");
        row.querySelector(".mdi").classList.replace("mdi-check", "mdi-plus");
    };

    selectPatternWeights(labels);
});


const createPatternListItem = (item, cart) => {
    const gctype = cart.gctype;
    const num_genes = cart.num_genes;
    const text = `${cart.label} (${num_genes} genes)`;

    item.querySelector('.pattern-list-item-label').textContent = text;
    item.querySelector('.ul-li').dataset.shareId = cart.share_id;
    item.querySelector('.ul-li').dataset.label = cart.label;
    item.querySelector('.ul-li').dataset.gctype = gctype;
    item.querySelector('.ul-li').dataset.organismId = cart.organism_id;

    if (selectedPattern.shareId == cart.share_id) {
        item.querySelector('.ul-li').classList.add('is-selected');
        item.querySelector('.ul-li').classList.remove('is-clickable');
    } else {
        item.querySelector('.ul-li').classList.remove('is-selected');
        item.querySelector('.ul-li').classList.add('is-clickable');
    }

    // create tag
    item.querySelector('.tag').textContent = gctype;
    // give tags different colors based on the gctype
    if (gctype === "unweighted-list") {
        item.querySelector('.tag').classList.add('is-info');
    } else if (gctype === "weighted-list") {
        item.querySelector('.tag').classList.add('is-success');
    } else if (gctype === "labeled-list") {
        item.querySelector('.tag').classList.add('is-warning');
    }

    // if unweighted list, hide the arrow-right icon
    if (gctype === "unweighted-list") {
        item.querySelector('.dropdown-pattern-list-icon').classList.add('is-hidden');
    }

    document.getElementById('dropdown-content-pattern-lists').appendChild(item);

    // Get item after it's been added to the DOM
    const thisItem = document.querySelector(`.dropdown-pattern-list-item[data-share-id="${cart.share_id}"]`);

    // Event listener to select the pattern list and get the weights for the pattern.  Populate the weights dropdown with the weights.
    thisItem.addEventListener("click", async (event) => {
        // uncheck all the existing rows
        const rows = document.getElementsByClassName('dropdown-pattern-list-item');
        for (const row of rows) {
            row.classList.remove('is-selected');
            row.classList.add('is-clickable');
        };

        selectedPattern.shareId = event.currentTarget.dataset.shareId;
        selectedPattern.label = event.currentTarget.dataset.label;
        selectedPattern.gctype = event.currentTarget.dataset.gctype;
        selectedPattern.selectedWeights = [];

        event.currentTarget.classList.add('is-selected');
        event.currentTarget.classList.remove('is-clickable');

        // These buttons have no bearing on unweighted lists
        document.getElementById("weighted-shortcut-label").classList.add('is-hidden');
        document.getElementById("dropdown-pattern-list-clear-weights").classList.add('is-hidden');
        document.getElementById("dropdown-pattern-list-top5-weights").classList.add('is-hidden');
        document.getElementById("dropdown-pattern-list-top10-weights").classList.add('is-hidden');
        document.getElementById("dropdown-pattern-list-top20-weights").classList.add('is-hidden');
        document.getElementById("dropdown-pattern-list-all-weights").classList.add('is-hidden');
        if (gctype === "weighted-list") {
            document.getElementById("weighted-shortcut-label").classList.remove('is-hidden');
            document.getElementById("dropdown-pattern-list-clear-weights").classList.remove('is-hidden');
            document.getElementById("dropdown-pattern-list-top5-weights").classList.remove('is-hidden');
            document.getElementById("dropdown-pattern-list-top10-weights").classList.remove('is-hidden');
            document.getElementById("dropdown-pattern-list-top20-weights").classList.remove('is-hidden');
            document.getElementById("dropdown-pattern-list-all-weights").classList.remove('is-hidden');
        }


        await populatePatternWeights();
    });
}

/**
 * Fetches patterns data asynchronously and executes a callback function.
 * @param {Function} callback - The callback function to be executed after fetching patterns data.
 * @returns {Promise<void>} - A promise that resolves when the patterns data is fetched successfully.
 */
const fetchPatternsData = async (callback) => {
    try {
        patternsCartData = await apiCallsMixin.fetchGeneCarts({includeMembers: false});

        flatPatternsCartData = [...patternsCartData.domain_carts, ...patternsCartData.group_carts, ...patternsCartData.public_carts, ...patternsCartData.user_carts, ...patternsCartData.shared_carts, ...patternsCartData.recent_carts]
        document.getElementById('dropdown-pattern-lists').classList.remove('is-loading');
        document.getElementById('dropdown-pattern-lists').classList.remove('is-disabled');

        if (callback) {
            await callback();
        }

    } catch (error) {
        console.error(error);
    }
}

/**
 * Sets the active pattern cart category and updates the pattern list accordingly.
 * @param {string} category - The category of the pattern cart.
 */
const setActivePatternCartCategory = (category) => {
    // clear the pattern list
    document.getElementById('dropdown-pattern-list-search-input').value = '';

    const patternListItemTemplate = document.getElementById('tmpl-pattern-list-item');
    let data = null;

    document.getElementById('dropdown-content-pattern-lists').innerHTML = '';

    let recentChosen = false;

    switch (category) {
        case 'favorites':
            data = patternsCartData.domain_carts;
            break;
        case 'recent':
            recentChosen = true;
            break;
        case 'saved':
            data = patternsCartData.user_carts;
            break;
        case 'shared':
            data = patternsCartData.shared_carts;
            break;
    }

    // Prevent from breaking for now
    if (recentChosen) {
        return;
    }

    if (!data) {
        // User may be logged out and selected "saved"
        return;
    }

    for (const entry of data) {
        const row = patternListItemTemplate.content.cloneNode(true);
        createPatternListItem(row, entry);
    }
}

/**
 * Populates the weights dropdown with data fetched from the API.
 * @returns {Promise<void>} A promise that resolves once the weights dropdown is populated.
 */
const populatePatternWeights = async () => {
    // data is a list of weight and top/buttom genes (if weighted-list) and if binary weights
    const data = await apiCallsMixin.fetchPatternElementList(selectedPattern.shareId, selectedPattern.gctype)

    // All weights are selected by default
    selectedPattern.selectedWeights = data;

    // Use the weight info to populate the weights dropdown (tmpl-weight-item)
    document.getElementById('dropdown-content-weights').innerHTML = '';
    const weightListItemTemplate = document.getElementById('tmpl-weight-item');
    for (const weight of data) {
        const row = weightListItemTemplate.content.cloneNode(true);
        row.querySelector('.weight-item-label').textContent = weight.label;

        row.querySelector('.ul-li').dataset.label = weight.label;
        // These only show up for weighted lsits
        if (weight.top_up) {
            row.querySelector('.ul-li').dataset.top_up = weight.top_up;
            row.querySelector('.ul-li').dataset.top_down = weight.top_down;
        }

        // If there is just one weight, we are obviously going to select it
        if (data.length === 1) {
            row.querySelector('.ul-li').classList.add("is-disabled");
            row.querySelector('.icon').classList.add("is-hidden");
        }

        document.getElementById('dropdown-content-weights').appendChild(row);

        const thisRow = document.querySelector(`.dropdown-weight-item[data-label="${weight.label}"]`);

        // Event listener to select the weight and update the selectedPattern.selectedWeights
        // Multiple weights can be selected

        thisRow.addEventListener("click", (event) => {
            // create object from event.currentTarget.dataset
            const obj = {};
            for (const key in event.currentTarget.dataset) {
                obj[key] = event.currentTarget.dataset[key];
            }

            if (event.currentTarget.classList.contains('is-selected')) {
                event.currentTarget.classList.remove('is-selected');
                selectedPattern.selectedWeights = selectedPattern.selectedWeights.filter((weight) => weight.label !== obj.label);
                // change mdi-plus to mdi-minus
                event.currentTarget.querySelector('.mdi').classList.remove('mdi-check');
                event.currentTarget.querySelector('.mdi').classList.add('mdi-plus');
                return;
            }
            event.currentTarget.classList.add('is-selected');

            // This is done so that the Proxy object can detect the change
            const currentWeights = selectedPattern.selectedWeights;
            currentWeights.push(obj);
            selectedPattern.selectedWeights = currentWeights;

            // change mdi-minus to mdi-plus
            event.currentTarget.querySelector('.mdi').classList.remove('mdi-plus');
            event.currentTarget.querySelector('.mdi').classList.add('mdi-check');
        });
    }

}

/**
 * Updates the selected pattern list with the given share ID and updates the pattern list selector label.
 * @param {string} shareId - The share ID to set for the selected pattern.
 */
const selectPatternList = (shareId) => {
    const foundPattern = document.querySelector(`.dropdown-pattern-list-item[data-share-id="${shareId}"]`)
    if (!foundPattern) {
        console.error(`Pattern with share ID ${shareId} not found`);
        return;
    }
    selectedPattern.shareId = shareId;
    selectedPattern.label = document.querySelector(`.dropdown-pattern-list-item[data-share-id="${shareId}"]`).dataset.label;
    selectedPattern.gctype = document.querySelector(`.dropdown-pattern-list-item[data-share-id="${shareId}"]`).dataset.gctype;

    // click the pattern list to select it
    document.querySelector(`.dropdown-pattern-list-item[data-share-id="${shareId}"]`).click();

    updatePatternListSelectorLabel();
}

const selectPatternWeights = (labels) => {

    // uncheck all the existing rows
    const rows = document.getElementsByClassName('.dropdown-weight-item');
    for (const row of rows) {
        row.classList.remove('is-selected');
    };

    // select our labels
    for (const label of labels) {
        document.querySelector(`.dropdown-weight-item[data-label="${label}"]`).click();
    }

}

/**
 * Updates the pattern list select drop down label.
 */
const updatePatternListSelectorLabel = () => {
    const shareId = selectedPattern.shareId;
    document.querySelector('#dropdown-pattern-list-selector-label').innerHTML = shareId ? selectedPattern.label : 'Quick search using pattern sources';
}