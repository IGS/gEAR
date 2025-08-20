"use strict";

import { apiCallsMixin } from "../../js/common.v2.js";

export const geneCollectionState = {
    "data": null,
    "labelIndex": {},
    "shareIdsToGenes": {},
    "selectedGeneLists": new Set(),
    "selectedGenes": new Set(),
    "manuallyEnteredGenes": new Set()
};

// Add event listener to dropdown trigger
document.querySelector("#dropdown-gene-lists > button.dropdown-trigger").addEventListener("click", (event) => {
    const item = event.currentTarget;

    // Close all other dropdowns
    const dropdowns = document.querySelectorAll('.dropdown.is-active');
    dropdowns.forEach((dropdown) => {
        if (dropdown !== item.closest('.dropdown')) {
            dropdown.classList.remove('is-active');
        }
    });

    item.closest(".dropdown").classList.toggle('is-active');

    // in case it was in errored state
    document.querySelector('#dropdown-gene-lists button').classList.remove('is-danger');
});

// Add event listeners to the gene list category selectors
const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li');
categorySelectors.forEach((element) => {
    element.addEventListener('click', (event) => {
        const category = event.target.dataset.category;
        setActiveGeneCartCategory(category);

        categorySelectors.forEach((element) => {
            element.classList.remove('is-selected');
            element.classList.add('is-clickable');
        });

        event.target.classList.add('is-selected');
        event.target.classList.remove('is-clickable');
    });
});

// Add event listeners to the gene list selectors even if they don't exist yet
// TODO: These click events should match what is in the dataset-collection and pattern-collection selectors
document.addEventListener('click', (event) => {
    // Map class names to actions
    const actionMap = {
        'gene-list-item-label': 'view',
        'dropdown-gene-list-item-right-selector': 'view',
        'dropdown-gene-list-item-add': 'add',
        'dropdown-gene-list-item-remove': 'remove'
    };

    // gene-list-item-label & dropdown-gene-list-item-right-selector both should only show the genes
    // dropdown-gene-list-item-add should add the entire cart

    // Find the first matching class in the map
    const matchedClass = Object.keys(actionMap).find(cls => event.target.classList.contains(cls));
    if (matchedClass) {
        const rowDiv = event.target.closest('div');
        setActiveGeneCart(rowDiv, actionMap[matchedClass]);
    }
});

// Add event listeners to the gene selectors even if they don't exist yet
document.addEventListener('click', (event) => {
    if (event.target.classList.contains('gene-list-item-add') ||
        event.target.classList.contains('gene-list-item-remove')) {

        const rowDiv = event.target.closest('div');
        const geneSymbol = rowDiv.querySelector('.gene-item-label').textContent;

        if (geneCollectionState.selectedGenes.has(geneSymbol)) {
            geneCollectionState.selectedGenes.delete(geneSymbol);
            rowDiv.querySelector('i.toggler').classList.replace("gene-list-item-remove", "gene-list-item-add")
        } else {
            geneCollectionState.selectedGenes.add(geneSymbol);
            rowDiv.querySelector('i.toggler').classList.replace("gene-list-item-add", "gene-list-item-remove")
        }

        rowDiv.classList.toggle('is-selected');
        rowDiv.querySelector('i.toggler').classList.toggle('mdi-plus');
        rowDiv.querySelector('i.toggler').classList.toggle('mdi-check');
    }
});

// Add a click listener to the dropdown-gene-list-cancel button
document.querySelector('#dropdown-gene-list-cancel').addEventListener('click', (event) => {
    // clear gene lists and gene list areas
    document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
    document.querySelector('#dropdown-content-genes').innerHTML = '';

    // reset the label container
    document.querySelector('#gene-select-dropdown-dynamic-selections').innerHTML = '';

    const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    // clear the searched gene list input element
    document.querySelector('#dropdown-gene-list-search-input').value = '';
    document.querySelector('#dropdown-gene-list-selector-label').innerHTML = 'Quick search using Gene Lists';

    // and finally the related gene lists and genes
    geneCollectionState.selectedGeneLists.clear();
    geneCollectionState.selectedGenes.clear();

    // Add back any manually-entered genes
    if (geneCollectionState.manuallyEnteredGenes.size > 0) {
        geneCollectionState.selectedGenes = new Set([...geneCollectionState.manuallyEnteredGenes,]);
    }
});

document.querySelector('#dropdown-gene-list-proceed').addEventListener('click', (event) => {
    updateGeneListSelectorLabel();

    // close the dropdown
    document.querySelector('#dropdown-gene-lists').classList.remove('is-active');
});

// Monitor key strokes after user types more than 2 characters in the dropdown-gene-list-search-input box
document.querySelector('#dropdown-gene-list-search-input').addEventListener('keyup', (event) => {
    const searchTerm = event.target.value;

    if (searchTerm.length === 0) {
        // clear the gene list
        document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
        document.querySelector('#dropdown-content-genes').innerHTML = '';
        return;
    } else if (searchTerm.length <= 2) {
        return;
    }

    const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li');
    categorySelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
    document.querySelector('#dropdown-content-genes').innerHTML = '';
    const geneListItemTemplate = document.querySelector('#tmpl-gene-list-item');

    let listShareIdsFound = new Set();

    for (const listType in geneCollectionState.data) {
        for (const cart of geneCollectionState.data[listType]) {
            if (cart.label.toLowerCase().includes(searchTerm.toLowerCase())) {

                if (! listShareIdsFound.has(cart.share_id)) {
                    const row = geneListItemTemplate.content.cloneNode(true);
                    row.querySelector('.gene-list-item-label').textContent = cart.label;
                    row.querySelector('.ul-li').dataset.shareId = cart.share_id;
                    row.querySelector('.ul-li').dataset.genes = geneCollectionState.shareIdsToGenes[cart.share_id].join(',');

                    if (geneCollectionState.selectedGeneLists.has(cart.share_id)) {
                        row.querySelector('i.toggler').classList.remove('mdi-plus');
                        row.querySelector('i.toggler').classList.add('mdi-check');
                        row.querySelector('.ul-li').classList.add('is-selected');
                    } else {
                        row.querySelector('i.toggler').classList.remove('mdi-check');
                        row.querySelector('i.toggler').classList.add('mdi-plus');
                        row.querySelector('.ul-li').classList.remove('is-selected');
                    }

                    document.querySelector('#dropdown-content-gene-lists').appendChild(row);

                    listShareIdsFound.add(cart.share_id);
                }
            }
        }
    }
});

/**
 * Fetches gene cart data from the server.
 *
 * @param {string|null} shareId - The share ID of the gene list. Default is null.
 * @returns {Promise<void>} - A promise that resolves when the gene cart data is fetched successfully.
 */
export const fetchGeneCartData = async (shareId=null) => {
    try {
        geneCollectionState.data = await apiCallsMixin.fetchGeneCarts({gcShareId: shareId, cartType: 'unweighted-list', includeMembers: false});
        document.querySelector('#dropdown-gene-lists').classList.remove('is-loading');
        document.querySelector('#dropdown-gene-lists').classList.remove('is-disabled');

        // Build the gene cart label index for ease of use
        for (const cartType in geneCollectionState.data) {
            for (const cart of geneCollectionState.data[cartType]) {
                geneCollectionState.labelIndex[cart.share_id] = cart.label;
                geneCollectionState.shareIdsToGenes[cart.share_id] = cart.genes;
            }
        }

    } catch (error) {
        console.error(error);
    }
}

const setActiveGeneCart = async (cartRow, mode) => {
    // TODO: this needs a spinner while it loads

    // clear the current gene list
    document.querySelector('#dropdown-content-genes').innerHTML = '';

    // populate the gene list from this cart
    let gene_list_member_data = await apiCallsMixin.fetchGeneCartMembers(cartRow.dataset.shareId);
    let genes = [];

    for (const member of gene_list_member_data.gene_symbols) {
        genes.push(member.label);
    }

    const gene_item_template = document.querySelector('#tmpl-gene-item');

    for (const gene of genes.sort()) {
        const gene_row = gene_item_template.content.cloneNode(true);
        gene_row.querySelector('.gene-item-label').textContent = gene;
        document.querySelector('#dropdown-content-genes').appendChild(gene_row);
    }

    // if adding or removing, update the inventory
    if (mode === 'add') {
        geneCollectionState.selectedGeneLists.add(cartRow.dataset.shareId);
        geneCollectionState.selectedGenes = new Set([...geneCollectionState.selectedGenes, ...genes, ...geneCollectionState.manuallyEnteredGenes]);
    } else if (mode === 'remove') {
        geneCollectionState.selectedGeneLists.delete(cartRow.dataset.shareId);

        for (const gene of genes) {
            geneCollectionState.selectedGenes.delete(gene);
        }
    }

    // now handle the coloring, icons and selection box based on the mode
    if (mode === 'add') {
        cartRow.querySelector('i.toggler').classList.remove('mdi-plus');
        cartRow.querySelector('i.toggler').classList.add('mdi-check')
        cartRow.querySelector('i.toggler').classList.remove('dropdown-gene-list-item-add');
        cartRow.querySelector('i.toggler').classList.add('dropdown-gene-list-item-remove');
        cartRow.classList.add('is-selected');

    } else if (mode === 'remove') {
        cartRow.querySelector('i.toggler').classList.remove('mdi-check');
        cartRow.querySelector('i.toggler').classList.add('mdi-plus');
        cartRow.querySelector('i.toggler').classList.add('dropdown-gene-list-item-add');
        cartRow.querySelector('i.toggler').classList.remove('dropdown-gene-list-item-remove');
        cartRow.classList.remove('is-selected');

    } else if (mode === 'view') {
        // do nothing
    }

    for (const gene_div of document.querySelectorAll('.dropdown-gene-item')) {
        const gene_symbol = gene_div.querySelector('.gene-item-label').textContent;

        if (geneCollectionState.selectedGenes.has(gene_symbol)) {
            gene_div.classList.add('is-selected');
            gene_div.querySelector('i.toggler').classList.remove('mdi-plus');
            gene_div.querySelector('i.toggler').classList.add('mdi-check');
        } else {
            gene_div.classList.remove('is-selected');
            gene_div.querySelector('i.toggler').classList.remove('mdi-check');
            gene_div.querySelector('i.toggler').classList.add('mdi-plus');
        }
    }

    // update the panel of selected items
    updateGeneListSelectionPanel();
}

const setActiveGeneCartCategory = (category) => {
    // clear the gene list
    document.getElementById('dropdown-content-genes').innerHTML = '';
    document.getElementById('dropdown-gene-list-search-input').value = '';

    const gene_list_item_template = document.querySelector('#tmpl-gene-list-item');
    let data = null;

    document.querySelector('#dropdown-content-gene-lists').innerHTML = '';

    switch (category) {
        case 'favorites':
            data = geneCollectionState.data.domain_carts;
            break;
        case 'recent':
            data = geneCollectionState.data.recent_carts;
            break;
        case 'saved':
            data = geneCollectionState.data.user_carts;
            break;
        case 'shared':
            data = geneCollectionState.data.shared_carts;
            break;
    }

    if (!data) {
        // User may be logged out and selected "saved"
        return;
    }

    for (const entry of data) {
        const row = gene_list_item_template.content.cloneNode(true);
        row.querySelector('.gene-list-item-label').textContent = entry.label;
        row.querySelector('.ul-li').dataset.shareId = entry.share_id;
        // This is from when we were pre-loading the gene lists, which isn't practical
        //row.querySelector('.ul-li').dataset.genes = geneCollectionState.shareIdsToGenes[entry.share_id].join(',');


        if (geneCollectionState.selectedGeneLists.has(entry.share_id)) {
            row.querySelector('i.toggler').classList.remove('mdi-plus');
            row.querySelector('i.toggler').classList.add('mdi-check');
            row.querySelector('.ul-li').classList.add('is-selected');
        } else {
            row.querySelector('i.toggler').classList.remove('mdi-check');
            row.querySelector('i.toggler').classList.add('mdi-plus');
            row.querySelector('.ul-li').classList.remove('is-selected');
        }

        document.getElementById('dropdown-content-gene-lists').appendChild(row);
    }
}

const selectGeneLists = (shareIds) => {
    console.debug('selectGeneLists', shareIds);
    console.debug(geneCollectionState.shareIdsToGenes);
    // reads the gene list share_ids passed and handles any UI and data updates to make
    //   them preselected
    for (const shareId of shareIds) {
        geneCollectionState.selectedGeneLists.add(shareId);
        geneCollectionState.selectedGenes = new Set([...geneCollectionState.selectedGenes, ...geneCollectionState.shareIdsToGenes[shareId], ...geneCollectionState.manuallyEnteredGenes]);
    }

    updateGeneListSelectorLabel();
}

const updateGeneListSelectionPanel = () => {
    const selection_box = document.querySelector('#gene-select-dropdown-dynamic-selections');

    // first empty it out, then populate it
    selection_box.innerHTML = '';

    for (const cart_share_id of geneCollectionState.selectedGeneLists) {
        selection_box.innerHTML += `<span class="tag is-info is-light is-small m-1">${geneCollectionState.labelIndex[cart_share_id]}</span>`;
    }
}

const updateGeneListSelectorLabel = () => {
    // Updates the gene list select drop down label based on the current state of the geneCollectionState.selectedGeneLists array
    const selected_cart_count = geneCollectionState.selectedGeneLists.size;

    if (selected_cart_count === 1) {
        // It's the only one
        const only_cart_id = Array.from(geneCollectionState.selectedGeneLists)[0];
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = geneCollectionState.labelIndex[only_cart_id];
    } else if (selected_cart_count > 1) {
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = `${selected_cart_count} gene lists selected`;
    } else {
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = 'Quick search using Gene Lists';
    }
}