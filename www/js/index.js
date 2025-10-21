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

    // For this page we want the gene collection dropdown to be right-aligned
    const geneCollectionDropdown = document.querySelector('#dropdown-gene-lists');
    geneCollectionDropdown.classList.remove('is-left');
    geneCollectionDropdown.classList.add('is-right');

    // For this page we want the dataset collection dropdown to be left-aligned
    const datasetCollectionDropdown = document.querySelector('#dropdown-dc');
    datasetCollectionDropdown.classList.remove('is-right');
    datasetCollectionDropdown.classList.add('is-left');

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

    loadOrganismList();
    populateDatasetSpinner();
});


document.getElementById("onboarding-btn").addEventListener("click", (event) => {
    // Start the onboarding process
    // NOTE: Click out of the intro will exit the intro
    introJs().setOptions({
    steps: [{
        intro: "Welcome to the main page."
    }, {
        element: document.querySelector("aside.menu"),
        intro: "Here you can navigate to different analysis tools as well as find or manage your gene lists and datasets.",
    }, {
        intro: "To perform a basic search, a collection of datasets, as well as one or more genes are necessary."
    }, {
        intro: "Let's start with searching for a gene in a collection of datasets. <br/><br/> When searching for genes, you can either select a gene list or enter genes manually.",
    }, {
        element: document.getElementById("genes-manually-entered"),
        intro: "You can enter genes manually here. You can enter multiple genes separated by commas or spaces.",
    }, {
        element: document.getElementById("dropdown-gene-lists"),
        intro: "You can also optionally select a gene list from the dropdown. <br/><br/> These genes will be added to the manually entered genes with the genes in the list.",
    }, {
        element: document.querySelector('#gene-search-exact-match').parentElement,
        intro: "You can choose to search for exact matches of the gene symbols. <br/><br/> If this is not checked, the search will match any gene with this text.",
    }, {
        element: document.querySelector('input[name="single-multi"][value="multi"]').parentElement,
        intro: "Currently plots specific to single genes will be shown. <br/><br/> If 'Multi-gene Display' is selected, plots that use all searched genes (such as heatmaps) will be shown instead.",
    }, {
        element: document.getElementById("dropdown-dc"),
        intro: "Next, select a dataset collection to search in. <br/><br/> These collections contain datasets under a certain theme, such as from a specific region, organism, or from a research paper",
    }, {
        element: document.getElementById("submit-expression-search"),
        intro: "Once you have selected your genes and dataset collections, click here to search for expression data.",
    }, {
        element: document.getElementById("user-history-table"),
        intro: "Finally, you can view your user history here. This will show you your recent searches and actions."
    }


    ]
    }).start();
});

const populateDatasetSpinner = async () => {
    let spinnerDatasets = [];
    const data = await apiCallsMixin.fetchDatasets({'limit': 10});
    spinnerDatasets = data.datasets;
   
    const template = document.querySelector('#dataset-card-template');
    const datasetSpinnerContainer = document.querySelector('#highlighted-datasets');
    datasetSpinnerContainer.innerHTML = '';

    let cardIdx = 0;

    for (const dataset of spinnerDatasets) {
        const card = template.content.cloneNode(true);
        //card.querySelector('.dataset-link').setAttribute('href', `/dataset.html?dataset_id=${dataset.dataset_id}`);
        //card.querySelector('.dataset-link').textContent = dataset.title;
        datasetSpinnerContainer.appendChild(card);
        cardIdx += 1;

        if (cardIdx >= 3) {
            break;
        }
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

/**
 * Loads the list of organisms from the server and populates the organism choices on the 
 * datasets tab
 * @function
 * @returns {void}
 */
const loadOrganismList = async () => {
    try {
        const data = await apiCallsMixin.fetchOrganismList();
        const organismChoices = document.getElementById("organism-choices"); // <select> element
        for (const organism of data.organisms) {
            const option = document.createElement("option");
            option.value = organism.id;
            option.textContent = organism.label;
            organismChoices.appendChild(option);
        }
    } catch (error) {
        logErrorInConsole(error);
        createToast("Failed to load organism list");
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
        fetchDatasetCollections()
    ]);
    document.getElementById("submit-expression-search").classList.remove("is-loading");

    // Trigger the default dataset collection to be selected in the
    if (CURRENT_USER.layout_share_id) {
        selectDatasetCollection(CURRENT_USER.layout_share_id);
    }

}
