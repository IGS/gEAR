'use strict';

let url_params_passed = false;
let currently_selected_gene_symbol = null;
let currently_selected_org_id = "";
let is_multigene = false;

document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.querySelector('#page-header-label').textContent = 'Gene Expression Search';

    // handle when the dropdown-gene-list-search-input input box is changed
    document.querySelector('#genes-manually-entered').addEventListener('change', (event) => {
        const search_term_string = event.target.value;

        if (search_term_string.length > 0) {
            // split the string into an array of genes by spaces or commas
            const manually_entered_genes = search_term_string.split(/[ ,]+/);
            selected_genes = [...new Set([...selected_genes, ...manually_entered_genes])];
        }
    });

    // Handle passed URL parameters
    if (getUrlParameter('gene_symbol_exact_match') === 'true') {
        document.querySelector('#gene-search-exact-match').checked = true;
    }

    // add event listener for when the submit-expression-search button is clicked
    document.querySelector('#submit-expression-search').addEventListener('click', (event) => {
        // TODO: If I clear the gene symbol search, then click, the gene symbol is not updated and passes validation below - SAdkins

        const status = validateExpressionSearchForm();

        if (! status) {
            console.log("Aborting search");
            return;
        }

        fetchGeneAnnotations();
    });

    // handle when the organism-selector select box is changed
    document.querySelector('#organism-selector').addEventListener('change', (event) => {
        currently_selected_org_id = document.querySelector('#organism-selector').value;

        if (currently_selected_org_id === "") {
            showOrganismSelectorToolip();
            return;
        } else {
            hideOrganismSelectorToolip();
            updateAnnotationDisplay();
        }

    });

    // Add event listeners to the gene result list items even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (!event.target.classList.contains('gene-result-list-item')) {
            return;
        }
        const gene_symbol = event.target.innerHTML;

        // remove is-selected from all the existing rows, then add it to this one
        const rows = document.querySelectorAll('.gene-result-list-item');
        rows.forEach((row) => {
            row.classList.remove('is-selected');
        });

        event.target.classList.add('is-selected');
        selectGeneResult(gene_symbol);
    });
});

const fetchGeneAnnotations = async (callback) => {
    try {
        const annotation_data = await apiCallsMixin.fetchGeneAnnotations(
            selected_genes.join(','),
            document.querySelector('#gene-search-exact-match').checked
        );

        document.querySelector('#gene-result-count').innerHTML = Object.keys(annotation_data).length;

        if (Object.keys(annotation_data).length === 0) {
            const no_history_template = document.querySelector('#tmpl-gene-result-none-found');
            document.querySelector('#gene-result-list').appendChild(no_history_template.content.cloneNode(true));
        } else {
            const template = document.querySelector('#tmpl-gene-result-item');
            document.querySelector('#gene-result-list').innerHTML = '';

            for (const gene_symbol in annotation_data) {
                const row = template.content.cloneNode(true);
                row.querySelector('li').innerHTML = gene_symbol;
                document.querySelector('#gene-result-list').appendChild(row);

                // due to a python issue, at some point in depth the data becomes a string. Parse it.
                for (const organism_id in annotation_data[gene_symbol]['by_organism']) {
                    const annot = JSON.parse(annotation_data[gene_symbol]['by_organism'][organism_id][0]);
                    annotation_data[gene_symbol]['by_organism'][organism_id] = annot;
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
        createToast("There was an error fetching the organism list (" + error + ")");
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Wait until all pending API calls have completed before checking if we need to search
    const [cart_result, dc_result, org_result] = await Promise.all([
        fetchGeneCartData(parseGeneCartURLParams),
        fetchDatasetCollections(parseDatasetCollectionURLParams),
        fetchOrganisms()
    ]);

    // Now, if URL params were passed and we have both genes and a dataset collection,
    //  run the search
    if (url_params_passed) {
        if (selected_dc_share_id && selected_genes.length > 0) {
            document.querySelector('#submit-expression-search').click();
        }
    }
}

const hideOrganismSelectorToolip = () => {
    document.querySelector('#organism-selector-control').classList.remove('has-tooltip-top', 'has-tooltip-arrow', 'has-tooltip-active');
    document.querySelector('#organism-selector-control').removeAttribute('data-tooltip');
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
        selectGeneLists(gene_lists); // declared in gene-collection-selector.js
        url_params_passed = true;
    }

    // are we doing exact matches?
    const exact_match = getUrlParameter('gene_symbol_exact_match');
    document.querySelector('#gene-search-exact-match').checked = exact_match === 'true';

    // single or multiple gene view (convert to boolean)?
    const is_multigene_param = getUrlParameter('is_multigene');
    is_multigene = is_multigene_param === '1';
    if (is_multigene) {
        document.querySelector('#single-multi-multi').checked = true;
    } else {
        document.querySelector('#single-multi-single').checked = true;
    }
}

const parseDatasetCollectionURLParams = async () => {
    // handle passed dataset collection
    const layout_share_id = getUrlParameter('layout_id');

    if (!layout_share_id) {
        return;
    }

    selected_dc_share_id = layout_share_id;
    selected_dc_label = dataset_collection_label_index[layout_share_id];
    document.querySelector('#dropdown-dc-selector-label').innerHTML = selected_dc_label;

    const tile_grid = new TileGrid(layout_share_id, "#result-panel-grid");
    try {
        tile_grid.layout = await tile_grid.getLayout();
        tile_grid.tilegrid = tile_grid.generateTileGrid();
        tile_grid.applyTileGrid(is_multigene);
    }   catch (error) {
        logErrorInConsole(error);
    }
}

const selectGeneResult = (gene_symbol) => {
    const selected_organism_id = document.querySelector('#organism-selector').value;
    currently_selected_gene_symbol = gene_symbol;

    // if no organism is selected, display a tooltip to choose one
    if (selected_organism_id === "") {
        showOrganismSelectorToolip();
        return;
    }

    updateAnnotationDisplay();
}

const showOrganismSelectorToolip = () => {
    document.querySelector('#organism-selector-control').setAttribute('data-tooltip', 'Select an organism to view annotation');
    document.querySelector('#organism-selector-control').classList.add('has-tooltip-top', 'has-tooltip-arrow', 'has-tooltip-active');
}

const updateAnnotationDisplay = () => {
    // these make some of the syntax below shorter
    const gs = currently_selected_gene_symbol;
    const oid = currently_selected_org_id;

    // clear the external resource links
    document.querySelector('#external-resource-links').innerHTML = '';

    // if the selected organism is not in the annotation data, show a message
    if (annotation_data[gs]['by_organism'].hasOwnProperty(oid)) {
        const annotation = annotation_data[gs]['by_organism'][oid];
        //console.log(annotation);

        document.querySelector('#currently-selected-gene-product').innerHTML = " - " + annotation['product'];
        document.querySelector('#currently-selected-gene-product').classList.remove('is-hidden');

        let good_dbxref_count = 0;

        for (const dbxref of annotation['dbxrefs']) {
            if (dbxref['url'] !== null) {
                const dbxref_template = document.querySelector('#tmpl-external-resource-link');
                const row = dbxref_template.content.cloneNode(true);
                row.querySelector('a').innerHTML = dbxref['source'];
                row.querySelector('a').href = dbxref['url'];
                document.querySelector('#external-resource-links').appendChild(row);
                good_dbxref_count++;
            }
        }

        if (good_dbxref_count === 0) {
            const dbxref_template = document.querySelector('#tmpl-external-resource-link-none-found');
            const row = dbxref_template.content.cloneNode(true);
            document.querySelector('#external-resource-links').appendChild(row);
        }
        return;


    }
    document.querySelector('#currently-selected-gene-product').innerHTML = " - (annotation not available for this organism)";
    document.querySelector('#currently-selected-gene-product').classList.remove('is-hidden');

    const dbxref_template = document.querySelector('#tmpl-external-resource-link-none-found');
    const dbxref_template_row = dbxref_template.content.cloneNode(true);
    document.querySelector('#external-resource-links').appendChild(dbxref_template_row);

}

const validateExpressionSearchForm = () => {
    // User must have either selected a gene list or entered genes manually. Either of these
    // will populate the selected_genes array
    if (!selected_genes.length) {
        createToast('Please enter at least one gene to proceed');
        return false;
    }

    // Check if the user has selected any dataset collections
    if (!selected_dc_share_id) {
        createToast('Please select at least one dataset to proceed');
        return false;
    }

    return true;
}
