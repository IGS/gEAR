"use strict";

let dataset_collection_data = null;
let dataset_collection_label_index = {};

let selected_dc_share_id = null;
let selected_dc_label = null;

// This many characters will be included and then three dots will be appended
const DATASET_COLLECTION_SELECTOR_PROFILE_LABEL_LENGTH_LIMIT = 35;

document.addEventListener('DOMContentLoaded', () => {
    // Add event listener to dropdown trigger
    document.querySelector("#dropdown-dc > button.dropdown-trigger").addEventListener("click", (event) => {
        const item = event.currentTarget;
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
        for (const category in dataset_collection_data) {
            // This data structure has mixed types - we only care about the arrayed categories
            if (!Array.isArray(dataset_collection_data[category])) {continue}

            for (const entry of dataset_collection_data[category]) {
                if (entry.label.toLowerCase().includes(search_term.toLowerCase())) {
                    const row = dc_item_template.content.cloneNode(true);
                    createDatasetCollectionListItem(row, entry);
                }
            }
        }
    });
});

/**
 * Fetches dataset collections from the API.
 *
 * @param {boolean} includeMembers - Whether to include members in the dataset collections.
 * @param {Function} callback - Optional callback function to be called after fetching dataset collections.
 * @returns {Promise<void>} - A promise that resolves when the dataset collections are fetched.
 */
const fetchDatasetCollections = async (includeMembers=false, callback) => {
    const layoutShareId = selected_dc_share_id || null;

    try {
        dataset_collection_data = await apiCallsMixin.fetchDatasetCollections({includeMembers, layoutShareId});
        document.querySelector('#dropdown-dc').classList.remove('is-loading');
        document.querySelector('#dropdown-dc').classList.remove('is-disabled');

        // build a label index for the dataset collections
        for (const category in dataset_collection_data) {
            // This data structure has mixed types - we only care about the arrayed categories
            if (!Array.isArray(dataset_collection_data[category])) {continue}

            for (const entry of dataset_collection_data[category]) {
                dataset_collection_label_index[entry.share_id] = entry.label;
            }
        }

        if (dataset_collection_data.selected) {
            selectDatasetCollection(dataset_collection_data.selected);
        }

        if (callback) {
            callback();
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
const selectDatasetCollection = (share_id) => {
    // reads the DC share_id passed and handles any UI and data updates to make
    //   it preselected

    const defaultLabel = "Choose a Dataset Collection";
    selected_dc_share_id = share_id || null;
    selected_dc_label = dataset_collection_label_index[selected_dc_share_id] || defaultLabel;

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
            data = dataset_collection_data.domain_layouts;
            break;
        case 'user':
            data = dataset_collection_data.user_layouts;
            break;
        case 'recent':
            //data = dataset_collection_data.recent_layouts;
            recentChosen = true;
            break;
        case 'group':
            data = dataset_collection_data.group_layouts;
            break;
        case 'shared':
            data = dataset_collection_data.shared_layouts;
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
        if (a.label < b.label) {return -1}
        if (a.label > b.label) {return 1}
        return 0;
    });

    for (const entry of data) {
        const row = dc_item_template.content.cloneNode(true);
        createDatasetCollectionListItem(row, entry);
    }
}

/**
 * Updates the label of the dataset collection selector based on the selected_dc_label.
 * If the selected_dc_label exceeds the length limit, it will be truncated and displayed with ellipsis.
 */
const updateDatasetCollectionSelectorLabel = () => {
    if (selected_dc_label.length > DATASET_COLLECTION_SELECTOR_PROFILE_LABEL_LENGTH_LIMIT) {
        const truncated_label = `${selected_dc_label.substring(0, DATASET_COLLECTION_SELECTOR_PROFILE_LABEL_LENGTH_LIMIT)}...`;
        document.querySelector('#dropdown-dc-selector-label').innerHTML = truncated_label;
    } else {
        document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;
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
        selected_dc_share_id = row_div.dataset.shareId;
        selected_dc_label = dataset_collection_label_index[selected_dc_share_id];

        updateDatasetCollectionSelectorLabel();

        document.querySelector('#dropdown-dc').classList.remove('is-active');
        document.querySelector('#dropdown-dc button').classList.remove('is-danger');
    });


}
