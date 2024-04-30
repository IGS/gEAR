"use strict";

let gene_cart_data = null;
let gene_cart_label_index = {};

// Build this where key is share_id and values are arrays of gene symbols
let gene_cart_genes = {};

let selected_gene_lists = new Set();
let selected_genes = new Set();

document.addEventListener('DOMContentLoaded', () => {

    // Add event listener to dropdown trigger
    document.querySelector("#dropdown-gene-lists > button.dropdown-trigger").addEventListener("click", (event) => {
        const item = event.currentTarget;
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
    document.addEventListener('click', (event) => {
        // gene-list-item-label & dropdown-gene-list-item-right-selector both should only show the genes
        // dropdown-gene-list-item-add should add the entire cart

        if (event.target.classList.contains('gene-list-item-label') ||
            event.target.classList.contains('dropdown-gene-list-item-right-selector')) {

            const row_div = event.target.closest('div');
            setActiveGeneCart(row_div, 'view');
        } else if (event.target.classList.contains('dropdown-gene-list-item-add')) {
            const row_div = event.target.closest('div');
            setActiveGeneCart(row_div, 'add');
        } else if (event.target.classList.contains('dropdown-gene-list-item-remove')) {
            const row_div = event.target.closest('div');
            setActiveGeneCart(row_div, 'remove');
        }
    });

    // Add event listeners to the gene selectors even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (event.target.classList.contains('gene-list-item-add') ||
            event.target.classList.contains('gene-list-item-remove')) {

            const row_div = event.target.closest('div');
            const gene_symbol = row_div.querySelector('.gene-item-label').textContent;

            if (selected_genes.has(gene_symbol)) {
                selected_genes.delete(gene_symbol);
                row_div.querySelector('i.toggler').classList.replace("gene-list-item-remove", "gene-list-item-add")
            } else {
                selected_genes.add(gene_symbol);
                row_div.querySelector('i.toggler').classList.replace("gene-list-item-add", "gene-list-item-remove")
            }

            row_div.classList.toggle('is-selected');
            row_div.querySelector('i.toggler').classList.toggle('mdi-plus');
            row_div.querySelector('i.toggler').classList.toggle('mdi-check');
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

        // clear the genes-manually-entered input element
        document.querySelector('#dropdown-gene-list-search-input').value = '';
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = 'Quick search using Gene Lists';

        // and finally the related gene lists and genes
        selected_gene_lists.clear();
        selected_genes.clear();
    });

    document.querySelector('#dropdown-gene-list-proceed').addEventListener('click', (event) => {
        updateGeneListSelectorLabel();

        // close the dropdown
        document.querySelector('#dropdown-gene-lists').classList.remove('is-active');
    });

    // Monitor key strokes after user types more than 2 characters in the dropdown-gene-list-search-input box
    document.querySelector('#dropdown-gene-list-search-input').addEventListener('keyup', (event) => {
        const search_term = event.target.value;

        if (search_term.length === 0) {
            // clear the gene list
            document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
            document.querySelector('#dropdown-content-genes').innerHTML = '';
            return;
        } else if (search_term.length <= 2) {
            return;
        }

        const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li');
        categorySelectors.forEach((element) => {
            element.classList.remove('is-selected');
            element.classList.add('is-clickable');
        });

        document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
        document.querySelector('#dropdown-content-genes').innerHTML = '';
        const gene_list_item_template = document.querySelector('#tmpl-gene-list-item');

        for (const cart_type in gene_cart_data) {
            for (const cart of gene_cart_data[cart_type]) {
                if (cart.label.toLowerCase().includes(search_term.toLowerCase())) {
                    const row = gene_list_item_template.content.cloneNode(true);
                    row.querySelector('.gene-list-item-label').textContent = cart.label;
                    row.querySelector('.ul-li').dataset.shareId = cart.share_id;
                    row.querySelector('.ul-li').dataset.genes = gene_cart_genes[cart.share_id].join(',');

                    if (selected_gene_lists.has(cart.share_id)) {
                        row.querySelector('i.toggler').classList.remove('mdi-plus');
                        row.querySelector('i.toggler').classList.add('mdi-check');
                        row.querySelector('.ul-li').classList.add('is-selected');
                    } else {
                        row.querySelector('i.toggler').classList.remove('mdi-check');
                        row.querySelector('i.toggler').classList.add('mdi-plus');
                        row.querySelector('.ul-li').classList.remove('is-selected');
                    }

                    document.querySelector('#dropdown-content-gene-lists').appendChild(row);
                }
            }
        }
    });
});

const fetchGeneCartData = async (callback) => {
    try {
        gene_cart_data = await apiCallsMixin.fetchGeneCarts('unweighted-list');
        document.querySelector('#dropdown-gene-lists').classList.remove('is-loading');
        document.querySelector('#dropdown-gene-lists').classList.remove('is-disabled');

        // Build the gene cart label index for ease of use
        for (const cart_type in gene_cart_data) {
            for (const cart of gene_cart_data[cart_type]) {
                gene_cart_label_index[cart.share_id] = cart.label;
                gene_cart_genes[cart.share_id] = cart.genes;

                // remove the genes list from the original data structure so we don't
                //  use the memory twice
                delete cart.genes;
            }
        }

        if (callback) {
            callback();
        }

    } catch (error) {
        console.error(error);
    }
}

const setActiveGeneCart = (cart_row, mode) => {
    // clear the current gene list
    document.querySelector('#dropdown-content-genes').innerHTML = '';

    // populate the gene list from this cart
    const genes = cart_row.dataset.genes.split(',');
    const gene_item_template = document.querySelector('#tmpl-gene-item');

    for (const gene of genes.sort()) {
        const gene_row = gene_item_template.content.cloneNode(true);
        gene_row.querySelector('.gene-item-label').textContent = gene;
        document.querySelector('#dropdown-content-genes').appendChild(gene_row);
    }

    // if adding or removing, update the inventory
    if (mode === 'add') {
        selected_gene_lists.add(cart_row.dataset.shareId);
        selected_genes = new Set([...selected_genes, ...genes]);
    } else if (mode === 'remove') {
        selected_gene_lists.delete(cart_row.dataset.shareId);

        for (const gene of genes) {
            selected_genes.delete(gene);
        }
    }

    // now handle the coloring, icons and selection box based on the mode
    if (mode === 'add') {
        cart_row.querySelector('i.toggler').classList.remove('mdi-plus');
        cart_row.querySelector('i.toggler').classList.add('mdi-check')
        cart_row.querySelector('i.toggler').classList.remove('dropdown-gene-list-item-add');
        cart_row.querySelector('i.toggler').classList.add('dropdown-gene-list-item-remove');
        cart_row.classList.add('is-selected');

    } else if (mode === 'remove') {
        cart_row.querySelector('i.toggler').classList.remove('mdi-check');
        cart_row.querySelector('i.toggler').classList.add('mdi-plus');
        cart_row.querySelector('i.toggler').classList.add('dropdown-gene-list-item-add');
        cart_row.querySelector('i.toggler').classList.remove('dropdown-gene-list-item-remove');
        cart_row.classList.remove('is-selected');

    } else if (mode === 'view') {
        // do nothing
    }

    for (const gene_div of document.querySelectorAll('.dropdown-gene-item')) {
        const gene_symbol = gene_div.querySelector('.gene-item-label').textContent;

        if (selected_genes.has(gene_symbol)) {
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
    document.querySelector('#dropdown-content-genes').innerHTML = '';
    document.querySelector('#dropdown-gene-list-search-input').value = '';

    const gene_list_item_template = document.querySelector('#tmpl-gene-list-item');
    let data = null;

    document.querySelector('#dropdown-content-gene-lists').innerHTML = '';

    switch (category) {
        case 'favorites':
            data = gene_cart_data.domain_carts;
            break;
        case 'recent':
            data = gene_cart_data.recent_carts;
            break;
        case 'saved':
            data = gene_cart_data.user_carts;
            break;
        case 'shared':
            data = gene_cart_data.shared_carts;
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
        row.querySelector('.ul-li').dataset.genes = gene_cart_genes[entry.share_id].join(',');


        if (selected_gene_lists.has(entry.share_id)) {
            row.querySelector('i.toggler').classList.remove('mdi-plus');
            row.querySelector('i.toggler').classList.add('mdi-check');
            row.querySelector('.ul-li').classList.add('is-selected');
        } else {
            row.querySelector('i.toggler').classList.remove('mdi-check');
            row.querySelector('i.toggler').classList.add('mdi-plus');
            row.querySelector('.ul-li').classList.remove('is-selected');
        }

        document.querySelector('#dropdown-content-gene-lists').appendChild(row);
    }
}

const selectGeneLists = (share_ids) => {
    // reads the gene list share_ids passed and handles any UI and data updates to make
    //   them preselected
    for (const share_id of share_ids) {
        selected_gene_lists.add(share_id);
        selected_genes = new Set([...selected_genes, ...gene_cart_genes[share_id]]);
    }

    updateGeneListSelectorLabel();
}

const updateGeneListSelectionPanel = () => {
    const selection_box = document.querySelector('#gene-select-dropdown-dynamic-selections');

    // first empty it out, then populate it
    selection_box.innerHTML = '';

    for (const cart_share_id of selected_gene_lists) {
        selection_box.innerHTML += `<span class="tag is-info is-light is-small m-1">${gene_cart_label_index[cart_share_id]}</span>`;
    }
}

const updateGeneListSelectorLabel = () => {
    // Updates the gene list select drop down label based on the current state of the selected_gene_lists array
    const selected_cart_count = selected_gene_lists.size;

    if (selected_cart_count === 1) {
        // It's the only one
        const only_cart_id = Array.from(selected_gene_lists)[0];
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = gene_cart_label_index[only_cart_id];
    } else if (selected_cart_count > 1) {
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = `${selected_cart_count} gene lists selected`;
    } else {
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = 'Quick search using Gene Lists';
    }
}