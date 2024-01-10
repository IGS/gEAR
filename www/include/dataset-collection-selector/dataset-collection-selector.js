
let dataset_collection_data = null;


document.addEventListener('DOMContentLoaded', () => {
    fetchDatasetCollections();

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

    // Add event listeners to the gene list selectors even if they don't exist yet
    document.addEventListener('click', (event) => {
        // gene-list-item-label & dropdown-gene-list-item-right-selector both should only show the genes
        // dropdown-gene-list-item-add should add the entire cartsmpl

        if (event.target.classList.contains('gene-list-item-label') || 
            event.target.classList.contains('dropdown-gene-list-item-right-selector')) {
            
            const row_div = event.target.closest('div');
            setActiveGeneCart(row_div, 'view');
        } else if (event.target.classList.contains('dropdown-gene-list-item-add')) {
            const row_div = event.target.closest('div');
            setActiveGeneCart(row_div, 'add');
        } else if (event.target.classList.contains('dropdown-gene-list-item-remove')) {
            const row_div = event.target.closest('div');
            setActiveGeneCart(row_div, 'remove');
        }
    });

    // Add event listeners to the gene selectors even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (event.target.classList.contains('gene-list-item-add') || 
            event.target.classList.contains('gene-list-item-remove')) {
            
            const row_div = event.target.closest('div');
            const gene_symbol = row_div.querySelector('.gene-item-label').textContent;

            if (selected_genes.includes(gene_symbol)) {
                selected_genes = selected_genes.filter((gene) => gene !== gene_symbol);
            } else {
                selected_genes.push(gene_symbol);
            }

            row_div.classList.toggle('is-selected');
            row_div.querySelector('i.toggler').classList.toggle('mdi-plus');
            row_div.querySelector('i.toggler').classList.toggle('mdi-minus');
        }
    });

    // Add a click listener to the dropdown-gene-list-cancel button
    document.querySelector('#dropdown-gene-list-cancel').addEventListener('click', (event) => {
        // clear gene lists and gene list areas
        document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
        document.querySelector('#dropdown-content-genes').innerHTML = '';

        // reset the label container
        document.querySelector('#gene-select-dropdown-dynamic-selections').innerHTML = '';

        const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li'); 
        categorySelectors.forEach((element) => {
            element.classList.remove('is-selected');
            element.classList.add('is-clickable');
        });

        // clear the genes-manually-entered input element
        document.querySelector('#genes-manually-entered').value = '';
        document.querySelector('#dropdown-gene-list-search-input').value = '';
        document.querySelector('#dropdown-gene-list-selector-label').innerHTML = 'Quick search using Gene Lists';

        // and finally the related gene lists and genes
        selected_carts = {};
        selected_genes = [];
    });

    document.querySelector('#dropdown-gene-list-proceed').addEventListener('click', (event) => {
        const selected_cart_count = Object.keys(selected_carts).length;

        if (selected_cart_count === 1) {
            document.querySelector('#dropdown-gene-list-selector-label').innerHTML = gene_cart_label_index[selected_carts[0]];
        } else if (selected_cart_count > 1) {
            document.querySelector('#dropdown-gene-list-selector-label').innerHTML = `${selected_cart_count} gene lists selected`;
        } else {
            // Nothing to do
        }

        // close the dropdown
        document.querySelector('#dropdown-gene-lists').classList.remove('is-active');
    });

    // Minor key strokes after user types more than 2 characters in the dropdown-gene-list-search-input box
    document.querySelector('#dropdown-gene-list-search-input').addEventListener('keyup', (event) => {
        const search_term = event.target.value;
         
        if (search_term.length <= 2) {return}

        const categorySelectors = document.querySelectorAll('#dropdown-content-gene-list-category .ul-li'); 
        categorySelectors.forEach((element) => {
            element.classList.remove('is-selected');
            element.classList.add('is-clickable');
        });

        document.querySelector('#dropdown-content-gene-lists').innerHTML = '';
        document.querySelector('#dropdown-content-genes').innerHTML = '';
        const gene_list_item_template = document.querySelector('#tmpl-gene-list-item');
        
        for (const cart_type in gene_cart_data) {
            for (const cart of gene_cart_data[cart_type]) {
                if (cart.label.toLowerCase().includes(search_term.toLowerCase())) {                    
                    const row = gene_list_item_template.content.cloneNode(true);
                    row.querySelector('.gene-list-item-label').textContent = cart.label;
                    row.querySelector('.ul-li').dataset.genes = cart.genes.join(',');
                    row.querySelector('.ul-li').dataset.shareId = cart.share_id;

                    if (cart.share_id in selected_carts) {
                        row.querySelector('i.toggler').classList.remove('mdi-plus');
                        row.querySelector('i.toggler').classList.add('mdi-minus');
                        row.querySelector('.ul-li').classList.add('is-selected');
                    } else {
                        row.querySelector('i.toggler').classList.remove('mdi-minus');
                        row.querySelector('i.toggler').classList.add('mdi-plus');
                        row.querySelector('.ul-li').classList.remove('is-selected');
                    }

                    document.querySelector('#dropdown-content-gene-lists').appendChild(row);
                }
            }
        }
    });
});

const fetchDatasetCollections = async () => {
    try {
        const data = await apiCallsMixin.fetchDatasetCollections();
        dataset_collection_data = data;
        document.querySelector('#dropdown-dc').classList.remove('is-loading');
        document.querySelector('#dropdown-dc').classList.remove('is-disabled');
    } catch (error) {
        console.error(error);
    }
}