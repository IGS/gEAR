'use strict';

let manually_entered_genes = [];

document.addEventListener('DOMContentLoaded', () => {

    // add event listener for when the dropdown-gene-list-search-input input box is changed
    document.querySelector('#genes-manually-entered').addEventListener('change', (event) => {
        const search_term_string = event.target.value;

        if (search_term_string.length > 0) {
            // split the string into an array of genes by spaces or commas
            manually_entered_genes = search_term_string.split(/[ ,]+/);
            selected_genes = [...new Set([...selected_genes, ...manually_entered_genes])];
        }
    });

    document.querySelector('#submit-expression-search').addEventListener('click', (event) => {
        const status = validateExpressionSearchForm();

        if (status) {
            // build the URL for a GET request
            let url = '/expression.html?';

            // add the manually-entered genes
            // TODO: need to combine selected_genes here to accommodate the case where a gene cart
            //  chosen but the individual genes removed.
            if (manually_entered_genes.length > 0) {
                url += `gene_symbol=${manually_entered_genes.join(',')}`;
            }

            // are we doing exact matches?
            if (document.querySelector('#gene-search-exact-match').checked) {
                url += '&gene_symbol_exact_match=1';
            }

            // get the value of the single-multi radio box
            const single_multi = document.querySelector('input[name="single-multi"]:checked').value;

            if (single_multi === 'single') {
                url += '&is_multigene=0';
            } else {
                url += '&is_multigene=1';
            }

            // add the gene lists
            //  TODO: This will only be for labeling purposes, since individual genes could have been
            //    deselected within
            if (Object.keys(selected_carts).length > 0) {
                url += `&gene_lists=${selected_carts.join(',')}`;
            }

            // add the dataset collections
            url += `&layout_id=${selected_dc_share_id}`;

            // now go there
            window.location.href = url;
        }
    });

    // Bulma gives the styling for tabs, but not the functionality
    const tabs = document.querySelectorAll('div.tabs li a');
    tabs.forEach((element) => {
        element.addEventListener('click', (event) => {
            const tab_id = element.closest('li').dataset.tabId;

            // First, hide all the tab-content elements
            document.querySelectorAll('.tabs-content li').forEach((element) => {
                if (element.dataset.tabId === tab_id) {
                    element.classList.add('is-active');
                } else {
                    element.classList.remove('is-active');
                }
            });
        });
    });
});

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

const validateExpressionSearchForm = () => {
    // User must have either selected a gene list or entered genes manually. Either of these
    // will populate the selected_genes array
    if (selected_genes.length + manually_entered_genes.length === 0) {
        createToast('Please enter at least one gene to proceed');
        return false;
    }

    // Check if the user has selected any dataset collections
    if (selected_dc_share_id === null) {
        createToast('Please select at least one dataset to proceed');
        return false;
    }

    return true;
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        populateUserHistoryTable();
    }

    fetchGeneCartData();
    fetchDatasetCollections();
}
