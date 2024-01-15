document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.querySelector('#page-header-label').textContent = 'Gene Expression Search';

    // Handle passed URL parameters
 

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
        toggleOnGeneLists(gene_lists);
    }
}

const parseDatasetCollectionURLParams = () => {
    // handle passed dataset collections
    //const layout_id = getUrlParameter('layout_id');
}