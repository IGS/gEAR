let gene_cart_data = null;

document.addEventListener('DOMContentLoaded', () => {
    fetchGeneCartData();

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

    // Add even listeners to the gene list selectors even if they don't exist yet
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
        }
    });
});

const fetchGeneCartData = async () => {
    try {
        gene_cart_data = await apiCallsMixin.fetchGeneCarts('unweighted-list');
        console.log(gene_cart_data);

        document.querySelector('#dropdown-gene-lists').classList.remove('is-loading');
        document.querySelector('#dropdown-gene-lists').classList.remove('is-disabled');
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
    const genes = cart_row.dataset.genes.split(',');
    const gene_item_template = document.querySelector('#tmpl-gene-item');
    
    document.querySelector('#dropdown-content-genes').innerHTML = '';

    if (mode === 'add') {
        cart_row.querySelector('i.toggler').classList.remove('mdi-plus');
        cart_row.querySelector('i.toggler').classList.add('mdi-minus');
    }

    for (const gene of genes.sort()) {
        const gene_row = gene_item_template.content.cloneNode(true);
        gene_row.querySelector('.gene-item-label').textContent = gene;
        gene_row_div = gene_row.querySelector('.gene-item-label').closest('div');

        if (mode === 'view') {
            gene_row_div.classList.remove('is-selected');
        } else {
            gene_row_div.classList.add('is-selected');

            gene_row_div.querySelector('i').classList.remove('mdi-plus');
            gene_row_div.querySelector('i').classList.add('mdi-minus');
        }

        document.querySelector('#dropdown-content-genes').appendChild(gene_row);

    }

    const geneListSelectors = document.querySelectorAll('#dropdown-content-gene-lists .ul-li');
    geneListSelectors.forEach((element) => {
        element.classList.remove('is-selected');
        element.classList.add('is-clickable');
    });

    cart_row.classList.add('is-selected');
    cart_row.classList.remove('is-clickable');
}

const setActiveGeneCartCategory = (category) => {
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
        document.querySelector('#dropdown-content-gene-lists').appendChild(row);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        populateUserHistoryTable();
    }
}
