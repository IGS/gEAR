let gene_cart_data = null;

document.addEventListener('DOMContentLoaded', () => {
    fetchGeneCartData();

    // select all #dropdown-content-gene-list-category ul-li elements
    const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li');
    
    categorySelectors.forEach((element) => {
        console.log("Adding event listener to: " + element);

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
});

const fetchGeneCartData = async () => {
    try {
        gene_cart_data = await apiCallsMixin.fetchGeneCarts();
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

const setActiveGeneCartCategory = (category) => {
    const template = document.querySelector('#tmpl-gene-list-item');
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
        const row = template.content.cloneNode(true);
        row.querySelector('.gene-list-item-label').textContent = entry.label;
        document.querySelector('#dropdown-content-gene-lists').appendChild(row);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        populateUserHistoryTable();
    }
}
