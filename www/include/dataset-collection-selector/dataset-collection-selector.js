let dataset_collection_data = null;
let dataset_collection_label_index = {};

let selected_dc_share_id = null;
let selected_dc_label = null;

// This many characters will be included and then three dots will be appended
const DATASET_COLLECTION_SELECTOR_PROFILE_LABEL_LENGTH_LIMIT = 35;

document.addEventListener('DOMContentLoaded', () => {
    // Add event listeners to the gene list category selectors
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

    // Add event listeners to the DC selectors even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (event.target.classList.contains('dropdown-dc-item') ||
            event.target.classList.contains('dc-item-label') || 
            event.target.classList.contains('dc-item-tag')) {
            
            // uncheck all the existing rows
            const rows = document.querySelectorAll('.dropdown-dc-item');
            rows.forEach((row) => {
                row.classList.remove('is-selected');
            });

            const row_div = event.target.closest('div');
            row_div.classList.toggle('is-selected');
            selected_dc_share_id = row_div.dataset.shareId;
            selected_dc_label = dataset_collection_label_index[selected_dc_share_id];

            if (selected_dc_label.length > DATASET_COLLECTION_SELECTOR_PROFILE_LABEL_LENGTH_LIMIT) {
                let truncated_label = selected_dc_label.substring(0, DATASET_COLLECTION_SELECTOR_PROFILE_LABEL_LENGTH_LIMIT) + '...';
                document.querySelector('#dropdown-dc-selector-label').innerHTML = truncated_label;
            } else {
                document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;
            }

            document.querySelector('#dropdown-dc').classList.remove('is-active');
        }
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
                    row.querySelector('.dc-item-label').textContent = entry.label;
                    row.querySelector('.ul-li').dataset.shareId = entry.share_id;

                    let tag_element = row.querySelector('.ul-li .dc-item-tag');

                    if (entry.folder_label) {
                        tag_element.textContent = entry.folder_label;
                    } else {
                        // we don't need the tag at all if there's no content for it
                        tag_element.remove();
                    }

                    document.querySelector('#dropdown-content-dc').appendChild(row);
                }
            }
        }
    });
});

const fetchDatasetCollections = async (callback) => {
    try {
        dataset_collection_data = await apiCallsMixin.fetchDatasetCollections();
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

        callback();

    } catch (error) {
        console.error(error);
    }
}

const setActiveDCCategory = (category) => {
    // clear the gene list
    document.querySelector('#dropdown-content-dc').innerHTML = '';
    document.querySelector('#dropdown-dc-search-input').value = '';

    const dc_item_template = document.querySelector('#tmpl-dc');
    let data = null;
    
    document.querySelector('#dropdown-content-dc').innerHTML = '';

    switch (category) {
        case 'domain':
            data = dataset_collection_data.domain_layouts;
            break;
        case 'user':
            data = dataset_collection_data.user_layouts;
            break;
        case 'recent':
            //data = dataset_collection_data.recent_layouts;
            break;
        case 'group':
            data = dataset_collection_data.group_layouts;
            break;
        case 'shared':
            data = dataset_collection_data.shared_layouts;
            break;
    }

    // sort the data by label before iterating
    data.sort((a, b) => {
        if (a.label < b.label) {return -1}
        if (a.label > b.label) {return 1}
        return 0;
    });

    for (const entry of data) {
        const row = dc_item_template.content.cloneNode(true);
        row.querySelector('.dc-item-label').textContent = entry.label;
        row.querySelector('.ul-li').dataset.shareId = entry.share_id;

        let tag_element = row.querySelector('.ul-li .dc-item-tag');

        if (entry.folder_label) {
            tag_element.textContent = entry.folder_label;
        } else {
            // we don't need the tag at all if there's no content for it
            tag_element.remove();
        }

        document.querySelector('#dropdown-content-dc').appendChild(row);
    }
}