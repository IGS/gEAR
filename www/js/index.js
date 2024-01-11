let gene_cart_data = null;
let gene_cart_label_index = {};

// For carts, key is share_id, value is array of genes symbol strings
let selected_carts = {};
let selected_genes = [];

document.addEventListener('DOMContentLoaded', () => {
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
        // dropdown-gene-list-item-add should add the entire cartsmpl

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

            if (selected_genes.includes(gene_symbol)) {
                selected_genes = selected_genes.filter((gene) => gene !== gene_symbol);
            } else {
                selected_genes.push(gene_symbol);
            }

            row_div.classList.toggle('is-selected');
            row_div.querySelector('i.toggler').classList.toggle('mdi-plus');
            row_div.querySelector('i.toggler').classList.toggle('mdi-minus');
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
        document.querySelector('#genes-manually-entered').value = '';
        document.querySelector('#dropdown-gene-list-search-input').value = '';
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = 'Quick search using Gene Lists';

        // and finally the related gene lists and genes
        selected_carts = {};
        selected_genes = [];
    });

    document.querySelector('#dropdown-gene-list-proceed').addEventListener('click', (event) => {
        const selected_cart_count = Object.keys(selected_carts).length;

        if (selected_cart_count === 1) {
            document.querySelector('#dropdown-gene-list-selector-label').innerHTML = gene_cart_label_index[selected_carts[0]];
        } else if (selected_cart_count > 1) {
            document.querySelector('#dropdown-gene-list-selector-label').innerHTML = `${selected_cart_count} gene lists selected`;
        } else {
            // Nothing to do
        }

        // close the dropdown
        document.querySelector('#dropdown-gene-lists').classList.remove('is-active');
    });

    // Minor key strokes after user types more than 2 characters in the dropdown-gene-list-search-input box
    document.querySelector('#dropdown-gene-list-search-input').addEventListener('keyup', (event) => {
        const search_term = event.target.value;
         
        if (search_term.length <= 2) {return}

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
                    row.querySelector('.ul-li').dataset.genes = cart.genes.join(',');
                    row.querySelector('.ul-li').dataset.shareId = cart.share_id;

                    if (cart.share_id in selected_carts) {
                        row.querySelector('i.toggler').classList.remove('mdi-plus');
                        row.querySelector('i.toggler').classList.add('mdi-minus');
                        row.querySelector('.ul-li').classList.add('is-selected');
                    } else {
                        row.querySelector('i.toggler').classList.remove('mdi-minus');
                        row.querySelector('i.toggler').classList.add('mdi-plus');
                        row.querySelector('.ul-li').classList.remove('is-selected');
                    }

                    document.querySelector('#dropdown-content-gene-lists').appendChild(row);
                }
            }
        }
    });
});

const fetchGeneCartData = async () => {
    console.log("Fetching gene cart data");
    console.log(CURRENT_USER);
    try {
        gene_cart_data = await apiCallsMixin.fetchGeneCarts('unweighted-list');
        document.querySelector('#dropdown-gene-lists').classList.remove('is-loading');
        document.querySelector('#dropdown-gene-lists').classList.remove('is-disabled');

        // Build the gene cart label index for ease of use
        for (const cart_type in gene_cart_data) {
            for (const cart of gene_cart_data[cart_type]) {                
                gene_cart_label_index[cart.share_id] = cart.label;
            }
        }
        
    } catch (error) {
        console.error(error);
    }
}

const populateUserHistoryTable = async () => {
    const numEntries = 5;

    // Load the spinner template
    const spinnerTemplate = document.querySelector('#user-history-loading');
    document.querySelector('#user-history-table-tbody').innerHTML = '';
    document.querySelector('#user-history-table-tbody').appendChild(spinnerTemplate.content.cloneNode(true));

    try {
        const data = await apiCallsMixin.fetchUserHistoryEntries(numEntries);
        const template = document.querySelector('#user-history-row');
        document.querySelector('#user-history-table-tbody').innerHTML = '';

        if (data.length === 0) {
            const noHistoryTemplate = document.querySelector('#user-history-no-entries');
            document.querySelector('#user-history-table-tbody').appendChild(noHistoryTemplate.content.cloneNode(true));
        } else {
            for (const entry of data) {
                const row = template.content.cloneNode(true);

                let formatted_category = entry.entry_category.replaceAll('_', ' ');
                formatted_category = formatted_category.charAt(0).toUpperCase() + formatted_category.slice(1);

                row.querySelector('.category').textContent = formatted_category;
                row.querySelector('.action-label').textContent = entry.label;
                row.querySelector('.date').textContent = entry.entry_date;
                row.querySelector('.url').setAttribute('href', entry.url);

                document.querySelector('#user-history-table-tbody').appendChild(row);
            }
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
        selected_carts[cart_row.dataset.shareId] = genes;
        selected_genes = [...new Set([...selected_genes, ...genes])];
    } else if (mode === 'remove') {
        delete selected_carts[cart_row.dataset.shareId];
        selected_genes = selected_genes.filter((gene) => !genes.includes(gene));
    }

    // now handle the coloring, icons and selection box based on the mode
    if (mode === 'add') {
        cart_row.querySelector('i.toggler').classList.remove('mdi-plus');
        cart_row.querySelector('i.toggler').classList.add('mdi-minus')
        cart_row.querySelector('i.toggler').classList.remove('dropdown-gene-list-item-add');
        cart_row.querySelector('i.toggler').classList.add('dropdown-gene-list-item-remove');
        cart_row.classList.add('is-selected');

    } else if (mode === 'remove') {
        cart_row.querySelector('i.toggler').classList.remove('mdi-minus');
        cart_row.querySelector('i.toggler').classList.add('mdi-plus');
        cart_row.querySelector('i.toggler').classList.add('dropdown-gene-list-item-add');
        cart_row.querySelector('i.toggler').classList.remove('dropdown-gene-list-item-remove');
        cart_row.classList.remove('is-selected');
        
    } else if (mode === 'view') {
        // do nothing
    }
    
    for (const gene_div of document.querySelectorAll('.dropdown-gene-item')) {
        let gene_symbol = gene_div.querySelector('.gene-item-label').textContent;

        if (selected_genes.includes(gene_symbol)) {
            gene_div.classList.add('is-selected');
            gene_div.querySelector('i.toggler').classList.remove('mdi-plus');
            gene_div.querySelector('i.toggler').classList.add('mdi-minus');            
        } else {
            gene_div.classList.remove('is-selected');
            gene_div.querySelector('i.toggler').classList.remove('mdi-minus');
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
            break;
        case 'saved':
            data = gene_cart_data.user_carts;
            break;
        case 'shared':
            data = gene_cart_data.shared_carts;
            break;
    }

    for (const entry of data) {
        const row = gene_list_item_template.content.cloneNode(true);
        row.querySelector('.gene-list-item-label').textContent = entry.label;
        row.querySelector('.ul-li').dataset.genes = entry.genes.join(',');
        row.querySelector('.ul-li').dataset.shareId = entry.share_id;

        if (entry.share_id in selected_carts) {
            row.querySelector('i.toggler').classList.remove('mdi-plus');
            row.querySelector('i.toggler').classList.add('mdi-minus');
            row.querySelector('.ul-li').classList.add('is-selected');
        } else {
            row.querySelector('i.toggler').classList.remove('mdi-minus');
            row.querySelector('i.toggler').classList.add('mdi-plus');
            row.querySelector('.ul-li').classList.remove('is-selected');
        }

        document.querySelector('#dropdown-content-gene-lists').appendChild(row);
    }
}

const updateGeneListSelectionPanel = () => {
    selection_box = document.querySelector('#gene-select-dropdown-dynamic-selections');

    // first empty it out, then populate it
    selection_box.innerHTML = '';

    for (const cart_share_id in selected_carts) {
        selection_box.innerHTML += `<span class="tag is-info is-light is-small m-1">${gene_cart_label_index[cart_share_id]}</span>`;
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        populateUserHistoryTable();
    }

    fetchGeneCartData();
    fetchDatasetCollections();
}
