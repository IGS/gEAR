"use strict";

import { apiCallsMixin, getCurrentUser } from "../../js/common.v2.js?v=a14c1ad";

export const datasetCollectionState = {
    "data": null,
    "labelIndex": {},
    "selectedShareId": null,
    "selectedLabel": null
}

// This many characters will be included and then three dots will be appended
const DatasetCollectionSelectorLabelMaxLength = 35;

// Add event listener to dropdown trigger
document.querySelector("#dropdown-dc > button.dropdown-trigger").addEventListener("click", (event) => {
    const item = event.currentTarget;
    // Close all other dropdowns
    const dropdowns = document.querySelectorAll('.dropdown.is-active');
    dropdowns.forEach((dropdown) => {
        if (dropdown !== item.closest('.dropdown')) {
            dropdown.classList.remove('is-active');
        }
    });
    item.closest(".dropdown").classList.toggle('is-active');
});

// Add event listeners to the dataset collection category selectors
const categorySelectors = document.querySelectorAll('#dropdown-content-dc-category .ul-li');
categorySelectors.forEach((element) => {
    element.addEventListener('click', (event) => {
        const category = event.target.dataset.category;
        setActiveDCCategory(category);

        categorySelectors.forEach((element) => {
            element.classList.remove('is-selected');
            element.classList.add('is-clickable');
        });

        event.target.classList.add('is-selected');
        event.target.classList.remove('is-clickable');
    });
});

// Add a click listener to the dancel button
document.querySelector('#dropdown-dc-cancel').addEventListener('click', (event) => {
    document.querySelector('#dropdown-dc-search-input').value = '';
    document.querySelector('#dropdown-content-dc').innerHTML = '';
    document.querySelector('#dropdown-dc').classList.remove('is-active');
});

// Monitor key strokes after user types more than 2 characters in the search box
document.querySelector('#dropdown-dc-search-input').addEventListener('keyup', (event) => {
    const search_term = event.target.value;

    if (search_term.length === 0) {
        document.querySelector('#dropdown-content-dc').innerHTML = '';
    }

    if (search_term.length <= 2) {return}

    const categorySelectors = document.querySelectorAll('#dropdown-content-dc-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    document.querySelector('#dropdown-content-dc').innerHTML = '';

    const dc_item_template = document.querySelector('#tmpl-dc');

    // build a label index for the dataset collections
    for (const category in datasetCollectionState.data) {
        // This data structure has mixed types - we only care about the arrayed categories
        if (!Array.isArray(datasetCollectionState.data[category])) {continue}

        for (const entry of datasetCollectionState.data[category]) {
            if (entry.label.toLowerCase().includes(search_term.toLowerCase())) {
                const row = dc_item_template.content.cloneNode(true);
                createDatasetCollectionListItem(row, entry);
            }
        }
    }
});

/**
 * Fetches dataset collections.
 *
 * @param {string} [shareId=null] - The share ID of the dataset collection.
 * @returns {Promise<void>} - A promise that resolves when the dataset collections are fetched.
 * @throws {Error} - If an error occurs during the fetch.
 */
export const fetchDatasetCollections = async (shareId=null) => {
    const layoutShareId = shareId || datasetCollectionState.selectedShareId || null;

    try {
        datasetCollectionState.data = await apiCallsMixin.fetchDatasetCollections({includeMembers: false, layoutShareId});

        document.querySelector('#dropdown-dc').classList.remove('is-loading');
        document.querySelector('#dropdown-dc').classList.remove('is-disabled');

        // build a label index for the dataset collections
        for (const category in datasetCollectionState.data) {
            // This data structure has mixed types - we only care about the arrayed categories
            if (!Array.isArray(datasetCollectionState.data[category])) {continue}

            for (const entry of datasetCollectionState.data[category]) {
                datasetCollectionState.labelIndex[entry.share_id] = entry.label;
            }
        }

        if (datasetCollectionState.data.selected) {
            selectDatasetCollection(datasetCollectionState.data.selected);
        }

    } catch (error) {
        console.error(error);
    }
}

/**
 * Selects a dataset collection based on the provided share ID.
 *
 * @param {string} share_id - The share ID of the dataset collection to be selected.
 */
export const selectDatasetCollection = (share_id) => {
    // reads the DC share_id passed and handles any UI and data updates to make
    //   it preselected

    const defaultLabel = "Choose a Dataset Collection";
    datasetCollectionState.selectedShareId = share_id || null;
    datasetCollectionState.selectedLabel = datasetCollectionState.labelIndex[datasetCollectionState.selectedShareId] || defaultLabel;

    updateDatasetCollectionSelectorLabel();
}

/**
 * Sets the active dataset collection category and updates the dropdown content accordingly.
 *
 * @param {string} category - The category to set as active. Possible values are 'domain', 'user', 'recent', 'group', and 'shared'.
 */
const setActiveDCCategory = (category) => {
    // clear the dataset collection search input and content
    document.querySelector('#dropdown-content-dc').innerHTML = '';
    document.querySelector('#dropdown-dc-search-input').value = '';

    const dc_item_template = document.querySelector('#tmpl-dc');
    let data = null;

    document.querySelector('#dropdown-content-dc').innerHTML = '';

    let recentChosen = false;

    switch (category) {
        case 'domain':
            data = datasetCollectionState.data.domain_layouts;
            break;
        case 'user':
            data = datasetCollectionState.data.user_layouts;
            break;
        case 'recent':
            //data = datasetCollectionState.data.recent_layouts;
            recentChosen = true;
            break;
        case 'group':
            data = datasetCollectionState.data.group_layouts;
            break;
        case 'shared':
            data = datasetCollectionState.data.shared_layouts;
            break;
    }

    if (recentChosen) {
        // prevent from breaking
        return;
    }

    if (!data) {
        // User may be logged out and selected "saved"
        return;
    }

    // sort the data by label before iterating
    data.sort((a, b) => {
        if (a.label.toUpperCase() < b.label.toUpperCase()) {return -1}
        if (a.label.toUpperCase() > b.label.toUpperCase()) {return 1}
        return 0;
    });

    for (const entry of data) {
        const row = dc_item_template.content.cloneNode(true);
        createDatasetCollectionListItem(row, entry);
    }
}

/**
 * Updates the label of the dataset collection selector based on the datasetCollectionState.selectedLabel.
 * If the datasetCollectionState.selectedLabel exceeds the length limit, it will be truncated and displayed with ellipsis.
 */
const updateDatasetCollectionSelectorLabel = () => {
    if (datasetCollectionState.selectedLabel.length > DatasetCollectionSelectorLabelMaxLength) {
        const truncated_label = `${datasetCollectionState.selectedLabel.substring(0, DatasetCollectionSelectorLabelMaxLength)}...`;
        document.querySelector('#dropdown-dc-selector-label').innerHTML = truncated_label;
    } else {
        document.querySelector('#dropdown-dc-selector-label').innerHTML = datasetCollectionState.selectedLabel;
    }
}

/**
 * Creates a dataset collection list item and updates its content based on the provided entry.
 * @param {HTMLElement} row - The row element representing the dataset collection list item.
 * @param {Object} entry - The entry object containing the data for the dataset collection.
 */
const createDatasetCollectionListItem = (row, entry) => {
    row.querySelector('.dc-item-label').textContent = entry.label;
    row.querySelector('.ul-li').dataset.shareId = entry.share_id;

    const tag_element = row.querySelector('.ul-li .dc-item-tag');

    if (entry.folder_label) {
        tag_element.textContent = entry.folder_label;
    } else {
        // we don't need the tag at all if there's no content for it
        tag_element.remove();
    }

    document.querySelector('#dropdown-content-dc').appendChild(row);

    // Create event listener to select the dataset collection
    const thisItem = document.querySelector(`.dropdown-dc-item[data-share-id="${entry.share_id}"]`);
    thisItem.addEventListener('click', (event) => {

        // uncheck all the existing rows
        const rows = document.getElementsByClassName('dropdown-dc-item');
        for (const row of rows) {
            row.classList.remove('is-selected');
        };

        const row_div = event.target.closest('div');
        row_div.classList.toggle('is-selected');
        datasetCollectionState.selectedShareId = row_div.dataset.shareId;
        datasetCollectionState.selectedLabel = datasetCollectionState.labelIndex[datasetCollectionState.selectedShareId];
        getCurrentUser().saveLayoutShareId(datasetCollectionState.selectedShareId);
        console.debug("Selected DC: ", datasetCollectionState.selectedShareId, datasetCollectionState.selectedLabel);

        updateDatasetCollectionSelectorLabel();

        document.querySelector('#dropdown-dc').classList.remove('is-active');
        document.querySelector('#dropdown-dc button').classList.remove('is-danger');
    });


}
