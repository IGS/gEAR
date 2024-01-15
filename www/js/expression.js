document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.querySelector('#page-header-label').textContent = 'Gene Expression Search';

    // Handle passed URL parameters
    if (getUrlParameter('gene_symbol_exact_match') === 'true') {
        document.querySelector('#gene-search-exact-match').checked = true;
    }

});

const handlePageSpecificLoginUIUpdates = async (event) => {
    fetchGeneCartData(parseGeneCartURLParams);
    fetchDatasetCollections(parseDatasetCollectionURLParams);
}

const parseGeneCartURLParams = () => {
    // handle manually-entered gene symbols
    const gene_symbols = getUrlParameter('gene_symbol');

    if (gene_symbols) {
        document.querySelector('#genes-manually-entered').value = gene_symbols.replaceAll(',', ' ');
    }

    // handle passed gene lists
    let gene_lists = [];
    if (getUrlParameter('gene_lists')) {
        gene_lists = getUrlParameter('gene_lists').split(',');
        selectGeneLists(gene_lists);
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