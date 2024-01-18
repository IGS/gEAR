'use strict';

let urlParamsPassed = false;
let annotationData = null;

let isMultigene = false;

document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.querySelector('#page-header-label').textContent = 'Gene Expression Search';

    // add event listener for when the dropdown-gene-list-search-input input box is changed
    document.querySelector('#genes-manually-entered').addEventListener('change', (event) => {
        const searchTermString = event.target.value;

        if (searchTermString.length > 0) {
            // split the string into an array of genes by spaces or commas
            const manuallyEnteredGenes = searchTermString.split(/[ ,]+/);
            selected_genes = [...new Set([...selected_genes, ...manuallyEnteredGenes])];
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

        fetchGeneAnnotations();
    });

    // Add event listeners to the gene result list items even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (!event.target.classList.contains('gene-result-list-item')) {
            return;
        }
        const geneSymbol = event.target.innerHTML;

        // remove is-selected from all the existing rows, then add it to this one
        const rows = document.querySelectorAll('.gene-result-list-item');
        rows.forEach((row) => {
            row.classList.remove('is-selected');
        });

        event.target.classList.add('is-selected');
        selectGeneResult(geneSymbol);
    });

});

const fetchGeneAnnotations = async (callback) => {
    try {
        const annotationData = await apiCallsMixin.fetchGeneAnnotations(
            selected_genes.join(','),
            document.querySelector('#gene-search-exact-match').checked
        );

        document.querySelector('#gene-result-count').innerHTML = Object.keys(annotationData).length;

        if (Object.keys(annotationData).length === 0) {
            const noHistoryTemplate = document.querySelector('#tmpl-gene-result-none-found');
            document.querySelector('#gene-result-list').appendChild(noHistoryTemplate.content.cloneNode(true));
        } else {
            const template = document.querySelector('#tmpl-gene-result-item');
            document.querySelector('#gene-result-list').innerHTML = '';

            for (const gene_symbol in annotationData) {
                const row = template.content.cloneNode(true);
                row.querySelector('li').innerHTML = gene_symbol;
                document.querySelector('#gene-result-list').appendChild(row);

                // due to a python issue, at some point in depth the data becomes a string. Parse it.
                for (const organism_id in annotationData[gene_symbol]['by_organism']) {
                    const annot = JSON.parse(annotationData[gene_symbol]['by_organism'][organism_id][0]);
                    annotationData[gene_symbol]['by_organism'][organism_id] = annot;
                }
            }
        }
    } catch (error) {
        console.error(error);
    }
}

const fetchOrganisms = async (callback) => {
    try {
        const orgs = await apiCallsMixin.fetchOrganismList();
        const template = document.querySelector('#tmpl-organism-option');

        for (const organism of orgs['organisms']) {
            const row = template.content.cloneNode(true);
            row.querySelector('option').innerHTML = organism.label;
            row.querySelector('option').value = organism.id;
            document.querySelector('#organism-selector').appendChild(row);
        }

    } catch (error) {
        console.error(error);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Wait until all pending API calls have completed before checking if we need to search
    const [cartResult, dcResult] = await Promise.all([
        fetchGeneCartData(parseGeneCartURLParams),
        fetchDatasetCollections(parseDatasetCollectionURLParams),
        fetchOrganisms()
    ]);

    // Now, if URL params were passed and we have both genes and a dataset collection,
    //  run the search
    if (urlParamsPassed) {
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

        urlParamsPassed = true;
    }

    // handle passed gene lists
    let geneLists = [];
    if (getUrlParameter('gene_lists')) {
        geneLists = getUrlParameter('gene_lists').split(',');
        selectGeneLists(geneLists); // declared in gene-collection-selector.js
        urlParamsPassed = true;
    }

    // are we doing exact matches?
    const exactMatch = getUrlParameter('gene_symbol_exact_match');
    document.querySelector('#gene-search-exact-match').checked = exactMatch === 'true';

    // single or multiple gene view (convert to boolean)?
    const isMultigeneParam = getUrlParameter('is_multigene');
    isMultigene = isMultigeneParam === '1';
    if (isMultigene) {
        document.querySelector('#single-multi-multi').checked = true;
    } else {
        document.querySelector('#single-multi-single').checked = true;
    }
}

const parseDatasetCollectionURLParams = async () => {
    // handle passed dataset collection
    const layoutShareId = getUrlParameter('layout_id');

    if (!layoutShareId) {
        return;
    }

    selected_dc_share_id = layoutShareId;
    selected_dc_label = dataset_collection_label_index[layoutShareId];
    document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;

    const tileGrid = new TileGrid(layoutShareId, "#result-panel-grid");
    try {
        tileGrid.layout = await tileGrid.getLayout();
        tileGrid.tilegrid = tileGrid.generateTileGrid();
        tileGrid.applyTileGrid(isMultigene);
    }   catch (error) {
        logErrorInConsole(error);
    }
}

const selectGeneResult = (geneSymbol) => {
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
