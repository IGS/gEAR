document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.querySelector('#page-header-label').textContent = 'Gene Expression Search';

    // Handle passed URL parameters


});

const handlePageSpecificLoginUIUpdates = async (event) => {
    fetchGeneCartData();
    fetchDatasetCollections();
}

