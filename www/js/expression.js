url_params_passed = false;

document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.querySelector('#page-header-label').textContent = 'Gene Expression Search';

    // add event listener for when the dropdown-gene-list-search-input input box is changed
    document.querySelector('#genes-manually-entered').addEventListener('change', (event) => {
        const search_term_string = event.target.value;

        if (search_term_string.length > 0) {
            // split the string into an array of genes by spaces or commas
            manually_entered_genes = search_term_string.split(/[ ,]+/);
            selected_genes = [...new Set([...selected_genes, ...manually_entered_genes])];
        }
    });

    // Handle passed URL parameters
    if (getUrlParameter('gene_symbol_exact_match') === 'true') {
        document.querySelector('#gene-search-exact-match').checked = true;
    }

    // add event listener for when the submit-expression-search button is clicked
    document.querySelector('#submit-expression-search').addEventListener('click', (event) => {
        const status = validateExpressionSearchForm();

        if (! status) {
            console.log("Aborting search");
            return;
        }

        console.log("Fetching gene annotations");
        fetchGeneAnnotations();
    });

});

const fetchGeneAnnotations = async (callback) => {
    try {
        const data = await apiCallsMixin.fetchGeneAnnotations(
            selected_genes.join(','),
            document.querySelector('#gene-search-exact-match').checked
        );
        console.log(data);



        /*
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
        */
    } catch (error) {
        console.error(error);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Wait until all pending API calls have completed before checking if we need to search
    const [cart_result, dc_result] = await Promise.all([
        fetchGeneCartData(parseGeneCartURLParams),
        fetchDatasetCollections(parseDatasetCollectionURLParams),
    ]);

    // Now, if URL params were passed and we have both genes and a dataset collection, 
    //  run the search
    if (url_params_passed) {
        if (selected_dc_share_id && selected_genes.length > 0) {
            document.querySelector('#submit-expression-search').click();
        }
    }
}

const parseGeneCartURLParams = () => {
    // handle manually-entered gene symbols
    const gene_symbols = getUrlParameter('gene_symbol');
    if (gene_symbols) {
        document.querySelector('#genes-manually-entered').value = gene_symbols.replaceAll(',', ' ');
        selected_genes = gene_symbols.split(',');

        url_params_passed = true;
    }

    // handle passed gene lists
    let gene_lists = [];
    if (getUrlParameter('gene_lists')) {
        gene_lists = getUrlParameter('gene_lists').split(',');
        selectGeneLists(gene_lists);
        url_params_passed = true;
    }

    // are we doing exact matches?
    const exact_match = getUrlParameter('gene_symbol_exact_match');
    if (exact_match === 'true') {
        document.querySelector('#gene-search-exact-match').checked = true;
    } else {
        document.querySelector('#gene-search-exact-match').checked = false;
    }
}

const parseDatasetCollectionURLParams = () => {
    // handle passed dataset collection
    const layout_id = getUrlParameter('layout_id');

    if (layout_id) {
        selected_dc_share_id = layout_id;
        selected_dc_label = dataset_collection_label_index[layout_id];
        document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;
    }
}

const validateExpressionSearchForm = () => {
    // User must have either selected a gene list or entered genes manually. Either of these
    // will populate the selected_genes array
    if (selected_genes.length === 0) {
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