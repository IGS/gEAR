"use strict";

// NOTE: This component depends on common.js and on Bulma CSS being imported in the parent HTML file

// Terminology:
// - patterns - The pattern source data
// - weights - The individual patterns

let patternsCartData = null;
let patternsCartLabelIndex = {};    // key is the share_id, value is the label

let selectedPattern = {shareId: null, label: null};

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

// Add a click listener to the dropdown-pattern-list-cancel button
document.querySelector('#dropdown-pattern-list-cancel').addEventListener('click', (event) => {
    // clear pattern lists and pattern list areas
    document.querySelector('#dropdown-content-pattern-lists').innerHTML = '';

    // reset the label container
    document.querySelector('#pattern-select-dropdown-dynamic-selections').innerHTML = '';

    const categorySelectors = document.querySelectorAll('#dropdown-content-pattern-list-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    // clear the patterns-manually-entered input element
    document.querySelector('#dropdown-pattern-list-search-input').value = '';
    document.querySelector('#dropdown-pattern-list-selector-label').innerHTML = 'Quick search using pattern Lists';

    // and finally the related pattern lists and patterns
    selectedPattern = {shareId: null, label: null};
});

// Monitor key strokes after user types more than 2 characters in the dropdown-pattern-list-search-input box
document.querySelector('#dropdown-pattern-list-search-input').addEventListener('keyup', (event) => {
    const search_term = event.target.value;

    if (search_term.length <= 2) {return}

    const categorySelectors = document.querySelectorAll('#dropdown-content-pattern-list-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    document.querySelector('#dropdown-content-pattern-lists').innerHTML = '';
    const pattern_list_item_template = document.querySelector('#tmpl-pattern-list-item');

    for (const cart_type in patternsCartData) {
        for (const cart of patternsCartData[cart_type]) {
            if (cart.label.toLowerCase().includes(search_term.toLowerCase())) {
                const row = pattern_list_item_template.content.cloneNode(true);
                row.querySelector('.pattern-list-item-label').textContent = cart.label;
                row.querySelector('.ul-li').dataset.shareId = cart.share_id;
                row.querySelector('.ul-li').dataset.patterns = patternsCartData[cart.share_id].join(',');

                if (selectedPattern.shareId == cart.share_id) {
                    row.querySelector('.ul-li').classList.add('is-selected');
                } else {
                    row.querySelector('.ul-li').classList.remove('is-selected');
                }

                document.querySelector('#dropdown-content-pattern-lists').appendChild(row);

                // Event listeners
                row.querySelector('.dropdown-pattern-list-item').click((event) => {
                    setActivePatternCart(row);
                });

            }
        }
    }
});

/**
 * Fetches pattern data from the server and performs necessary operations.
 * @param {Function} callback - Optional callback function to be executed after data is fetched.
 * @returns {Promise<void>} - A promise that resolves when the data is fetched and processed.
 */
const fetchPatternsData = async (callback) => {
    try {
        patternsCartData = await apiCallsMixin.fetchGeneCarts();
        document.querySelector('#dropdown-pattern-lists').classList.remove('is-loading');
        document.querySelector('#dropdown-pattern-lists').classList.remove('is-disabled');

        // Build the pattern cart label index for ease of use
        for (const cart_type in patternsCartData) {
            for (const cart of patternsCartData[cart_type]) {
                patternsCartLabelIndex[cart.share_id] = cart.label;
            }
        }

        if (callback) {
            callback();
        }

    } catch (error) {
        console.error(error);
    }
}

const setActivePatternCart = (cartRow) => {

    // populate the pattern list from this cart
    const patterns = cartRow.dataset.patterns.split(',');
    const patternItemTemplate = document.querySelector('#tmpl-pattern-item');

    for (const pattern of patterns.sort()) {
        const pattern_row = patternItemTemplate.content.cloneNode(true);
        pattern_row.querySelector('.pattern-item-label').textContent = pattern;
    }

    for (const patternDiv of document.querySelectorAll('.dropdown-pattern-item')) {
        const pattern = patternDiv.querySelector('.pattern-item-label').textContent;

        if (selectedPattern.label === pattern) {
            patternDiv.classList.add('is-selected');
            patternDiv.classList.remove('is-clickable');
        } else {
            patternDiv.classList.remove('is-selected');
            patternDiv.classList.add('is-clickable');

        }
    }
}

const setActivePatternCartCategory = (category) => {
    // clear the pattern list
    document.getElementById('dropdown-pattern-list-search-input').value = '';

    const patternListItemTemplate = document.querySelector('#tmpl-pattern-list-item');
    let data = null;

    document.querySelector('#dropdown-content-pattern-lists').innerHTML = '';

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

    for (const entry of data) {
        const gctype = entry.gctype;
        const num_genes = entry.num_genes;
        const text = `${entry.label} (${num_genes} genes)`;

        const row = patternListItemTemplate.content.cloneNode(true);
        row.querySelector('.pattern-list-item-label').textContent = text;
        row.querySelector('.ul-li').dataset.shareId = entry.share_id;

        if (selectedPattern.shareId === entry.share_id) {
            row.querySelector('.ul-li').classList.add('is-selected');
            row.querySelector('.ul-li').classList.remove('is-clickable');
        } else {
            row.querySelector('.ul-li').classList.remove('is-selected');
            row.querySelector('.ul-li').classList.add('is-clickable');
        }

        // create tag
        row.querySelector('.tag').textContent = gctype;
        // give tags different colors based on the gctype
        if (gctype === "unweighted-list") {
            row.querySelector('.tag').classList.add('is-info');
        } else if (gctype === "weighted-list") {
            row.querySelector('.tag').classList.add('is-success');
        } else if (gctype === "labeled-list") {
            row.querySelector('.tag').classList.add('is-warning');
        }

        document.querySelector('#dropdown-content-pattern-lists').appendChild(row);
    }
}

const selectPatternList = (shareId) => {
    selectedPattern.shareId = shareId;
    updatePatternListSelectorLabel();
}

const updatePatternListSelectorLabel = () => {
    // Updates the pattern list select drop down label

    const shareId = selectedPattern.shareId;
    if (shareId) {
        document.querySelector('#dropdown-pattern-list-selector-label').innerHTML = selectedPattern.label;
    } else {
        document.querySelector('#dropdown-pattern-list-selector-label').innerHTML = 'Quick search using pattern lists';
    }
}