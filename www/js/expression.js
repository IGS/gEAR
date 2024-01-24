'use strict';

let url_params_passed = false;
let currently_selected_gene_symbol = null;
let currently_selected_org_id = "";
let is_multigene = false;
let annotation_data = null;
let manually_entered_genes = [];
let tilegrid = null;
let svg_scoring_method = 'gene';

/*
TODOs:
- hide the annotation panel when multi-gene searches are displayed
- check if the user has a stored default profile and select that one if none were passed (index page too)
- If I clear the gene symbol search, then click, the gene symbol is not updated and passes validation below - SAdkins
*/

document.addEventListener('DOMContentLoaded', () => {
    // Set the page header title
    document.getElementById('page-header-label').textContent = 'Gene Expression Search';

    // Set current sidebar menu item to active
	for (const elt of document.querySelectorAll("#primary_nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}

	document.querySelector("a[tool='search_expression'").classList.add("is-active");


    // handle when the dropdown-gene-list-search-input input box is changed
    document.querySelector('#genes-manually-entered').addEventListener('change', (event) => {
        const search_term_string = event.target.value;

        if (search_term_string.length > 0) {
            // split the string into an array of genes by spaces or commas
            manually_entered_genes = search_term_string.split(/[ ,]+/);
            selected_genes = [...new Set([...selected_genes, ...manually_entered_genes])];
        }
    });

    document.querySelector('#functional-annotation-toggle').addEventListener('click', (event) => {
        const annotation_panel = document.querySelector('#extended-annotation-panel');
        const toggle_icon = document.querySelector('#functional-annotation-toggle i');

        if (annotation_panel.classList.contains('is-hidden')) {
            annotation_panel.classList.remove('is-hidden');
            toggle_icon.classList.remove('mdi-chevron-down');
            toggle_icon.classList.add('mdi-chevron-up');

        } else {
            annotation_panel.classList.add('is-hidden');
            toggle_icon.classList.remove('mdi-chevron-up');
            toggle_icon.classList.add('mdi-chevron-down');
        }
    });

    // add event listener for when the submit-expression-search button is clicked
    document.querySelector('#submit-expression-search').addEventListener('click', async (event) => {
        const status = validateExpressionSearchForm();

        if (! status) {
            console.log("Aborting search");
            return;
        }

        // update multigene/single gene
        is_multigene = document.querySelector('#single-multi-multi').checked;

        try {
            const [annotRes, tilegridRes] = await Promise.allSettled([fetchGeneAnnotations(), setupTileGrid(selected_dc_share_id)]);
            tilegrid = tilegridRes.value;

            // auto-select the first gene in the list
            const first_gene = document.querySelector('.gene-result-list-item');
            if (!is_multigene && first_gene) {
                first_gene.click();
            }

        } catch (error) {
            logErrorInConsole(error);
        }
    });

    // handle when the organism-selector select box is changed
    document.querySelector('#organism-selector').addEventListener('change', (event) => {
        currently_selected_org_id = parseInt(document.querySelector('#organism-selector').value);

        if (currently_selected_org_id === "") {
            showOrganismSelectorToolip();
            document.querySelector('#set-default-organism').classList.add('is-hidden');
            return;
        } else {
            hideOrganismSelectorToolip();
            updateAnnotationDisplay();

            // If the user is logged in and doesn't have a default org ID or it's different from their current,
            //  show the control
            if (CURRENT_USER.session_id) {
                if (CURRENT_USER.default_org_id === undefined || CURRENT_USER.default_org_id !== currently_selected_org_id) {
                    document.querySelector('#set-default-organism').classList.remove('is-hidden');
                } else {
                    document.querySelector('#set-default-organism').classList.add('is-hidden');
                }
            }
        }
    });

    document.querySelector('#set-default-organism').addEventListener('click', (event) => {
        // we don't want to set to null, and the UI should have prevented this, but check just in case
        if (currently_selected_org_id !== "") {
            CURRENT_USER.default_org_id = currently_selected_org_id;
            apiCallsMixin.saveUserDefaultOrgId(CURRENT_USER);
            document.querySelector('#set-default-organism').classList.add('is-hidden');
        }
    });

    // Add event listeners to the gene result list items even if they don't exist yet
    document.addEventListener('click', (event) => {
        if (!event.target.classList.contains('gene-result-list-item')) {
            return;
        }

        const gene_symbol = event.target.textContent;
        document.querySelector('#currently-selected-gene').innerHTML = gene_symbol;

        // remove is-selected from all the existing rows, then add it to this one
        const rows = document.querySelectorAll('.gene-result-list-item');
        rows.forEach((row) => {
            row.classList.remove('is-selected');
        });

        event.target.classList.add('is-selected');
        selectGeneResult(gene_symbol);
    });

    // Change the svg scoring method when select element is changed
    document.getElementById('svg-scoring-method').addEventListener('change', (event) => {
        if (is_multigene) return;   // multigene does not use this

        svg_scoring_method = event.target.value;
        // Get gene symbol from currently selected list item
        let list_item = document.querySelector('.gene-result-list-item.is-selected');
        if (!list_item) {
            list_item = document.querySelector('.gene-result-list-item');
        }

        const gene_symbol = list_item.textContent;
        selectGeneResult(gene_symbol);
    });
});

const fetchGeneAnnotations = async (callback) => {
    try {
        annotation_data = await apiCallsMixin.fetchGeneAnnotations(
            selected_genes.join(','),
            document.querySelector('#gene-search-exact-match').checked
        );

        //console.log(annotation_data);

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

            // if this matches the user's default organism, select it
            if (CURRENT_USER.default_org_id === organism.id) {
                row.querySelector('option').selected = true;
                currently_selected_org_id = organism.id;
            }

            document.querySelector('#organism-selector').appendChild(row);
        }

    } catch (error) {
        createToast("There was an error fetching the organism list (" + error + ")");
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Wait until all pending API calls have completed before checking if we need to search
    try {
        // SAdkins note - Promise.all fails fast,
        // but Promise.allSettled waits until all resolve/reject and lets you know which ones failed
        const [cart_result, dc_result, org_result] = await Promise.all([
            fetchGeneCartData(parseGeneCartURLParams),
            fetchDatasetCollections(parseDatasetCollectionURLParams),
            fetchOrganisms()
        ]);
    } catch (error) {
        logErrorInConsole(error);
    }

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
    document.querySelector('#gene-search-exact-match').checked = exact_match === '1';

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
}

const selectGeneResult = (gene_symbol) => {
    const selected_organism_id = document.querySelector('#organism-selector').value;
    currently_selected_gene_symbol = gene_symbol;

    // if no organism is selected, display a tooltip to choose one
    if (selected_organism_id === "") {
        showOrganismSelectorToolip();
    } else {
        updateAnnotationDisplay();
    }

    // Other things can be called next, such as plotting calls
    if (tilegrid) {
        tilegrid.renderDisplays(currently_selected_gene_symbol, is_multigene, svg_scoring_method);
    }
}

const setupTileGrid = async (layout_share_id) => {
    const tilegrid = new TileGrid(layout_share_id, "#result-panel-grid");
    try {
        tilegrid.layout = await tilegrid.getLayout();
        await tilegrid.addAllDisplays();

        tilegrid.generateTileGrid(is_multigene);
        tilegrid.applyTileGrid(is_multigene);
        await tilegrid.addDefaultDisplays();

        // NOTE - the tilegrid.renderDisplays() call below can check and use the first array element of the selected_genes array if single_gene
        // But we are using a string for clarity.
        if (is_multigene) {
            // Don't render yet if a gene is not selected
            if (selected_genes.length) {
                await tilegrid.renderDisplays(selected_genes, is_multigene);
            }
        } else {
            // Don't render yet if a gene is not selected
            if (currently_selected_gene_symbol) {
                await tilegrid.renderDisplays(currently_selected_gene_symbol, is_multigene, svg_scoring_method);
            }
        }
    } catch (error) {
        logErrorInConsole(error);
    } finally {
        return tilegrid;
    }
}

const showOrganismSelectorToolip = () => {
    document.querySelector('#organism-selector-control').setAttribute('data-tooltip', 'Select an organism to view annotation');
    document.querySelector('#organism-selector-control').classList.add('has-tooltip-top', 'has-tooltip-arrow', 'has-tooltip-active');
}

const updateAnnotationDisplay = () => {
    // these make some of the syntax below shorter
    const gs = currently_selected_gene_symbol;
    const oid = currently_selected_org_id;

    // clear the external resource links and GO terms
    document.querySelector('#external-resource-links').innerHTML = '';
    document.querySelector('#go-terms').innerHTML = '';
    document.querySelector('#go-term-count').innerHTML = '';

    // if the selected organism is not in the annotation data, show a message
    if (! annotation_data[gs]['by_organism'].hasOwnProperty(oid)) {
        document.querySelector('#currently-selected-gene-product').innerHTML = " - (annotation not available for this organism)";
        document.querySelector('#currently-selected-gene-product').classList.remove('is-hidden');

        const dbxref_template = document.querySelector('#tmpl-external-resource-link-none-found');
        const dbxref_template_row = dbxref_template.content.cloneNode(true);
        document.querySelector('#external-resource-links').appendChild(dbxref_template_row);
        return;
    }

    // if we got this far, we have annotation for this one. let's display it
    const annotation = annotation_data[gs]['by_organism'][oid];
    document.querySelector('#annotation-panel-gene-symbol').innerHTML = gs;

    // Gene product
    document.querySelector('#currently-selected-gene-product').innerHTML = " - " + annotation['product'];
    document.querySelector('#currently-selected-gene-product').classList.remove('is-hidden');
    document.querySelector('#annotation-panel-gene-product').innerHTML = annotation['product'];

    // aliases and Ensembl ID
    if (annotation['aliases'].length > 1) {
        document.querySelector('#annotation-panel-gene-aliases').innerHTML = annotation['aliases'].join(', ');
    } else {
        document.querySelector('#annotation-panel-gene-aliases').innerHTML = "None found";
    }

    document.querySelector('#annotation-panel-gene-ensembl-release').innerHTML = annotation['ensembl_release'];
    document.querySelector('#annotation-panel-gene-ensembl-id').innerHTML = annotation['ensembl_id'];

    const ensembl_url = "https://www.ensembl.org/Multi/Search/Results?q=" + annotation['ensembl_id'] + ";site=ensembl";
    document.querySelector('#annotation-panel-gene-ensembl-id').setAttribute('href', ensembl_url);

    // External database references
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

    // GO terms
    document.querySelector('#go-term-count').innerHTML = '(' + annotation['go_terms'].length + ')';
    if (annotation['go_terms'].length === 0) {
        const go_term_template = document.querySelector('#tmpl-go-term-none-found');
        const row = go_term_template.content.cloneNode(true);
        document.querySelector('#go-terms').appendChild(row);
    } else {
        for (const go_term of annotation['go_terms']) {
            const go_term_template = document.querySelector('#tmpl-go-term');
            const go_term_url = "https://amigo.geneontology.org/amigo/search/ontology?q=" + go_term['go_id'];

            const row = go_term_template.content.cloneNode(true);
            row.querySelector('.go-term-id').innerHTML = go_term['go_id'];
            row.querySelector('.go-term-id').href = go_term_url;
            row.querySelector('.go-term-label').innerHTML = go_term['name'];
            document.querySelector('#go-terms').appendChild(row);
        }
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
    if (!selected_dc_share_id) {
        createToast('Please select at least one dataset to proceed');
        return false;
    }

    return true;
}
