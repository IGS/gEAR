'use strict';

let manually_entered_genes = [];

document.addEventListener('DOMContentLoaded', () => {

    // if URL params are present, the user probably pasted a v1 URL into the browser
    //  so we should redirect them to the new expression URL but keep the same parameters
    if (window.location.search) {
        const url = new URL('/expression.html', window.location.origin);
        url.search = window.location.search;
        window.location = url.toString();
    }

    // handle when the dropdown-gene-list-search-input input box is changed
    document.getElementById('genes-manually-entered').addEventListener('change', (event) => {
        const search_term_string = event.target.value;
        let previously_manual_genes = manually_entered_genes;

        if (search_term_string.length > 0) {
            manually_entered_genes = search_term_string.split(/[ ,]+/);
        } else {
            manually_entered_genes = [];
        }

        // if any genes have been removed since last time, we need to remove them from the selected_genes array
        manually_entered_genes.forEach((gene) => {
            previously_manual_genes = previously_manual_genes.filter((g) => g !== gene);
        });

        previously_manual_genes.forEach((gene) => {
            selected_genes.delete(gene);
        });

        selected_genes = new Set([...selected_genes, ...manually_entered_genes]);
    });

    document.querySelector('#submit-expression-search').addEventListener('click', (event) => {
        const status = validateExpressionSearchForm();

        if (!status) {
            return;
        }

        // build the URL for a GET request
        const url = new URL('/expression.html', window.location.origin);

        // add the manually-entered genes
        // TODO: need to combine selected_genes here to accommodate the case where a gene cart
        //  chosen but the individual genes removed.
        const manuallyEnteredGenes =  Array.from(new Set([...selected_genes, ...manually_entered_genes]));
        if (manuallyEnteredGenes.length > 0) {
            url.searchParams.append('gene_symbol', manuallyEnteredGenes.join(','));
        }

        // are we doing exact matches?
        if (document.querySelector('#gene-search-exact-match').checked) {
            url.searchParams.append('gene_symbol_exact_match', '1');
        }

        // get the value of the single-multi radio box
        const singleMulti = document.querySelector('input[name="single-multi"]:checked').value;
        url.searchParams.append('is_multigene', singleMulti === 'single' ? '0' : '1');

        // add the gene lists
        //  TODO: This will only be for labeling purposes, since individual genes could have been
        //    deselected within
        if (selected_gene_lists.size > 0) {
            const geneCartShareIds = Array.from(selected_gene_lists);
            url.searchParams.append('gene_lists', geneCartShareIds.join(','));
        }

        // add the dataset collections
        url.searchParams.append('layout_id', selected_dc_share_id);

        // now go there
        window.location.href = url.toString();
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
    if (selected_genes.size + manually_entered_genes.length === 0) {
        createToast('Please enter at least one gene to proceed');
        document.querySelector('#dropdown-gene-lists button').classList.add('is-danger');
        return false;
    } else {
        document.querySelector('#dropdown-gene-lists button').classList.remove('is-danger');
    }

    // Check if the user has selected any dataset collections
    if (selected_dc_share_id === null) {
        createToast('Please select at least one dataset to proceed');
        document.querySelector('#dropdown-dc button').classList.add('is-danger');
        return false;
    } else {
        document.querySelector('#dropdown-dc button').classList.remove('is-danger');
    }

    return true;
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        populateUserHistoryTable();
    }

    document.getElementById("submit-expression-search").classList.add("is-loading");
    await Promise.all([
        fetchGeneCartData(),
        fetchDatasetCollections(false)
    ]);
    document.getElementById("submit-expression-search").classList.remove("is-loading");

    // Trigger the default dataset collection to be selected in the
    if (CURRENT_USER.layout_share_id) {
        selectDatasetCollection(CURRENT_USER.layout_share_id);
    }

}
