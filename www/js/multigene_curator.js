/*
  SAdkins note - The styling on this page differs from other gEAR JS pages
  I was trying to follow the recommended Javascript style in using camelCase
  in addition using the StandardJS style linter ("semistandrd" variant)
  with some refactoring done with the P42+ VSCode extension (to modernize code)
  However, code inherited from common.js is still in snake_case rather than camelCase
*/

'use strict';
/* global $, axios, Plotly, CURRENT_USER, session_id, check_for_login */

let numObs = 0; // dummy value to initialize with

let obsFilters = {};
let sortCategories = {"primary": null, "secondary": null}; // Control sorting order of chosen categories
let genesFilter = [];

let plotConfig = {};  // Plot config that is passed to API or stored in DB

let selectedGenes = null;  // Genes selected in a plot (e.g. volcano with lasso tool)

let datasetId = null;
let displayId = null;
let obsLevels = null;
let obsNotUsed = null;
let geneSymbols = null;

const datasetTree = new DatasetTree({treeDiv: '#dataset_tree'});
const geneCartTree = new GeneCartTree({treeDiv: '#gene_cart_tree'});

const plotTypes = ['dotplot', 'heatmap', 'mg_violin', 'quadrant', 'volcano'];

const dotplotOptsIds = ["#obs_primary_container", "#obs_secondary_container"];
const heatmapOptsIds = ["#heatmap_options_container", "#obs_primary_container", "#obs_secondary_container"];
const quadrantOptsIds = ["#quadrant_options_container", "#de_test_container"];
const violinOptsIds = ["#violin_options_container", "#obs_primary_container", "#obs_secondary_container"];
const volcanoOptsIds = ["#volcano_options_container", "#de_test_container"];

// color palettes
const continuousPalettes = [
    {
        label: "Multi-color scales",
        options: [
            { value: "YlOrRd", text: "Yellow-Orange-Red" },
            { value: "Viridis", text: "Viridis" },
        ],
    },
    {
        label: "Single-color scales",
        options: [
            { value: "Greys", text: "Greyscale" },
            { value: "Blues", text: "Bluescale" },
            { value: "Purp", text: "Purplescale" },
        ],
    },
    {
        label: "Diverging Colorscales",
        options: [
            { value: "RdBu", text: "Red-Blue" },
            { value: "PiYG", text: "Pink-Yellow-Green" },
        ],
    },
];
const discretePalettes = ["alphabet", "vivid", "light24", "dark24"];

window.onload= async () => {
    // Hide further configs until a dataset is chosen.
    // Changing the dataset will start triggering these to show
    $('#plot_type_container').hide();
    $('#gene_container').hide();

    // check if the user is already logged in
    check_for_login();

    // Load gene carts and datasets before the dropdown appears
    await reloadTrees ();

    // Initialize plot types
     $('#plot_type_select').select2({
        placeholder: 'Choose how to plot',
        width: '25%'
    });

    // If brought here by the "gene search results" page, curate on the dataset ID that referred us
    const linkedDatasetId = getUrlParameter("dataset_id");
    if (linkedDatasetId) {
        $("#dataset").val(linkedDatasetId);
        try {
            // Had difficulties triggering a "select_node.jstree" event, so just add the data info here
            const tree_leaf = datasetTree.treeData.find(e => e.dataset_id === linkedDatasetId);
            $("#dataset").text(tree_leaf.text);
            $("#dataset").data("organism-id", tree_leaf.organism_id);
            $("#dataset").data("dataset-id", tree_leaf.dataset_id);
            $("#dataset").trigger('change');
        } catch {
            console.error(`Dataset id ${linkedDatasetId} was not returned as a public/private/shared dataset`);
        }
    }

    // Create observer to watch if user changes (ie. successful login does not refresh page)
    // See: https://developer.mozilla.org/en-US/docs/Web/API/MutationObserver

    // But we need to wait for navigation_bar to load first (in common.js) so do some polling
    // See: https://stackoverflow.com/q/38881301

    // Select the node that will be observed for mutations
    const targetNode = document.getElementById('loggedin_controls');
    const saferNode = document.getElementById("navigation_bar");   // Empty div until loaded
    // Create an observer instance linked to the callback function
    const observer = new MutationObserver(function(mutationList, observer) {
        if (targetNode) {
            reloadTrees();
            this.disconnect();  // Don't need to reload once the trees are updated
        }
    });
    // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
    // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
    observer.observe(targetNode || saferNode , { attributes: true });
};

// Call API to return plot JSON data
async function getData (datasetId, payload) {
    try {
        return await axios.post(`/api/plot/${datasetId}/mg_dash`, {
            ...payload
        })
    } catch (e) {

        const message = "There was an error in making this plot. Please contact the gEAR team using the 'Contact' button at the top of the page and provide as much information as possible.";
        const success = -1;
        const data = {message, success};
        return {data};
    }
}

// Call API to return a list of the dataset's gene symbols
async function fetchGeneSymbols (payload) {
    const { datasetId, analysis } = payload;
    const base = `./api/h5ad/${datasetId}/genes`;
    const query = analysis ? `?analysis=${analysis.id}` : '';

    const { data } = await axios.get(`${base}${query}`);
    return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
}

// Call API to return all observations
async function fetchH5adInfo (payload) {
    const { datasetId, analysis } = payload;
    const base = `./api/h5ad/${datasetId}`;
    const query = analysis ? `?analysis=${analysis.id}` : '';
    const { data } = await axios.get(`${base}${query}`);
    return data;
}

async function loadDatasets () {
    $('#pre_dataset_spinner').show();
    await $.ajax({
        type: 'POST',
        url: './cgi/get_h5ad_dataset_list.cgi',
        data: {
            session_id
        },
        dataType: 'json',
        success(data) {
            let counter = 0;

            // Populate select box with dataset information owned by the user
            const userDatasets = [];
            if (data.user.datasets.length > 0) {
                // User has some profiles
                $.each(data.user.datasets, (_i, item) => {
                    if (item) {
                        userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                    }
                });
            }
            // Next, add datasets shared with the user
            const sharedDatasets = [];
            if (data.shared_with_user.datasets.length > 0) {
                $.each(data.shared_with_user.datasets, (_i, item) => {
                    if (item) {
                        sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                    }
                });
            }
            // Now, add public datasets
            const domainDatasets = [];
            if (data.public.datasets.length > 0) {
                $.each(data.public.datasets, (_i, item) => {
                    if (item) {
                        domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                    }
                });
            }

            datasetTree.userDatasets = userDatasets;
            datasetTree.sharedDatasets = sharedDatasets;
            datasetTree.domainDatasets = domainDatasets;
            datasetTree.generateTree();
        },
        error(_xhr, _status, msg) {
            console.error(`Failed to load dataset list because msg: ${msg}`);
        }
    });
    $('#pre_dataset_spinner').hide();
}

// Draw plotly chart to image
async function drawPreviewImage (display) {
    // check if config has been stringified
    const config = typeof display.plotly_config === 'string' ? JSON.parse(display.plotly_config) : display.plotly_config;

    const { data } = await getData(datasetId, config);

    // If there was an error in the plot, put alert up
    if ( data.success < 1 ) {
        $(`#modal-display-img-${display.id} + .js-plot-error`).show();
        $(`#modal-display-img-${display.id} + .js-plot-error`).html(data.message);
        $(`#modal-display-${display.id}-loading`).hide();
        return;
    }

    const { plot_json: plotlyJson, plot_config: plotlyConfig } = data;
    Plotly.toImage(
        { ...plotlyJson, plotlyConfig },
        { height: 500, width: 500 }
    ).then(url => {
        $(`#modal-display-img-${display.id}`).prop('src', url);
    }).then(
        () => { $(`#modal-display-${display.id}-loading`).hide(); }
    );
}

// Invert a log function
function invertLogFunction(value, base=10) {
    return base ** value;
}

// Draw plotly chart in HTML
function drawChart (data, datasetId) {
    const targetDiv = `dataset_${datasetId}_h5ad`;
    const parentDiv = `dataset_${datasetId}`;
    const { plot_json: plotlyJson, plot_config: plotlyConfig, message, success } = data;

    // Since default plots are now added after dataset selection, wipe the plot when a new one needs to be drawn
    $(`#${targetDiv}`).empty()

    // If there was an error in the plot, put alert up
    if ( success < 1 || !plotlyJson.layout) {
        $(`#${parentDiv} .js-plot-error`).show();
        $(`#${parentDiv} .js-plot-error`).html(message);
        return;
    }

    // Make some complex edits to the plotly layout
    const plotType = $('#plot_type_select').select2('data')[0].id;
    if (plotType === 'heatmap') {
        setHeatmapHeightBasedOnGenes(plotlyJson.layout, genesFilter);
    } else if (plotType=== "mg_violin" && $("#stacked_violin").is(":checked")){
        adjustStackedViolinHeight(plotlyJson.layout);
    }

    const configMods = {
        responsive: true
    };

    const config = {
        ...plotlyConfig,
        ...configMods
    };
    Plotly.newPlot(targetDiv, plotlyJson.data, plotlyJson.layout, config);

    // Update plot with custom plot config stuff stored in plot_display_config.js
    const curator_conf = post_plotly_config.curator;
    for (const idx in curator_conf) {
        const conf = curator_conf[idx];
        // Get config (data and/or layout info) for the plot type chosen, if it exists
        if (conf.plot_type == plotType) {
            const update_data = "data" in conf ? conf.data : {};
            const update_layout = "layout" in conf ? conf.layout : {};
            Plotly.update(targetDiv, update_data, update_layout)
        }
    }


    // Show any warnings from the API call
    if (message && success === 2) {
        $(`#${parentDiv} .js-plot-warning`).show();
        $(`#${parentDiv} .js-plot-warning`).html(`<ul>${message}</ul>`);
    }

    // If plot data is selected, create the right-column table and do other misc things
    $(`#dataset_${datasetId}_h5ad`).on("plotly_selected", (_e, data) => {

        if (!(['volcano', 'quadrant'].includes(plotConfig.plot_type))) {
            return;
        }

        // Note: the jQuery implementation has slightly different arguments than what is in the plotlyJS implementation
        // We want 'data', which returns the eventData PlotlyJS events normally return
        selectedGenes = [];

        data.points.forEach((pt) => {
            selectedGenes.push({
                gene_symbol: pt.data.text[pt.pointNumber],
                x: pt.data.x[pt.pointNumber].toFixed(1),
                y: plotConfig.plot_type === "volcano" ? invertLogFunction(-pt.data.y[pt.pointNumber]).toExponential(2) : pt.data.y[pt.pointNumber].toFixed(2),
            });
        });

        // Sort in alphabetical order
        selectedGenes.sort();

        // Adjust headers to the plot type
        if (plotConfig.plot_type === "quadrant") {
            $("#gene_x").html('X LFC <i class="fa fa-sort" aria-hidden="true"></i>');
            $("#gene_y").html('Y LFC <i class="fa fa-sort" aria-hidden="true"></i>');
        } else {
            // volcano
            $("#gene_x").html('LFC <i class="fa fa-sort" aria-hidden="true"></i>');
            $("#gene_y").html('Pval <i class="fa fa-sort" aria-hidden="true"></i>');
        }

        const template = $.templates("#selected_genes_tmpl");
        const htmlOutput = template.render(selectedGenes);
        $("#selected_genes_c").html(htmlOutput);

		// Highlight table rows that match searched genes
		if (genesFilter.length) {
			// Select the first column (gene_symbols) in each row
			$("#selected_genes_c tr td:first-child").each(function() {
				const tableGene = $(this).text();
				genesFilter.forEach((gene) => {
                    if (gene.toLowerCase() === tableGene.toLowerCase() ) {
                        $(this).parent().addClass("table-success");
                    }
				});
			})
		}

		$("#saved_gene_cart_info_c").hide();
        $("#tbl_selected_genes").show();
    });
}

// Submit API request and draw the HTML
async function draw (datasetId, payload) {
    const {data } = await getData(datasetId, payload);
    drawChart(data, datasetId);
}

// If user changes, update genecart/profile trees
async function reloadTrees(){
    // Update dataset and genecart trees in parallel
    // Works if they were not populated or previously populated
    await Promise.all([loadDatasets(), loadGeneCarts()]);
}

// Taken from https://www.w3schools.com/howto/howto_js_sort_table.asp
function sortTable(n) {
	let table;
	let rows;
	let switching;
	let i;
	let x;
	let y;
	let shouldSwitch;
	let dir;
	let switchcount = 0;
	table = document.getElementById("tbl_selected_genes");

	switching = true;
	// Set the sorting direction to ascending:
	dir = "asc";
	/* Make a loop that will continue until
		no switching has been done: */
	while (switching) {
		// Start by saying: no switching is done:
		switching = false;
		rows = table.rows;
		/* Loop through all table rows (except the
		first, which contains table headers): */
		for (i = 1; i < rows.length - 1; i++) {
		// Start by saying there should be no switching:
		shouldSwitch = false;
		/* Get the two elements you want to compare,
			one from current row and one from the next: */
		x = rows[i].getElementsByTagName("td")[n];
		y = rows[i + 1].getElementsByTagName("td")[n];
		/* Check if the two rows should switch place,
			based on the direction, asc or desc: */
		if (dir == "asc") {
			// First column is gene_symbol... rest are numbers
			if (n === 0 && x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
			// If so, mark as a switch and break the loop:
			shouldSwitch = true;
			break;
			}
			if (Number(x.innerHTML) > Number(y.innerHTML)) {
			shouldSwitch = true;
			break;
			}
		} else if (dir == "desc") {
			if (n === 0 && x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
			// If so, mark as a switch and break the loop:
			shouldSwitch = true;
			break;
			}
			if (Number(x.innerHTML) < Number(y.innerHTML)) {
			shouldSwitch = true;
			break;
			}
		}
		}
		if (shouldSwitch) {
		/* If a switch has been marked, make the switch
			and mark that a switch has been done: */
		rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
		switching = true;
		// Each time a switch is done, increase this count by 1:
		switchcount++;
		} else {
		/* If no switching has been done AND the direction is "asc",
			set the direction to "desc" and run the while loop again. */
		if (switchcount == 0 && dir == "asc") {
			dir = "desc";
			switching = true;
		}
		}
	}
}

// Render the gene-selection dropdown menu
function createGeneDropdown (genes) {
    const tmpl = $.templates('#gene_dropdown_tmpl');
    const html = tmpl.render({ genes });
    $('#gene_dropdown_container').html(html);
    $('#gene_dropdown').select2({
        placeholder: 'To search, click to select or start typing some gene names',
        allowClear: true,
        width: 'resolve'
    });
    $('#gene_spinner').hide();
}

// Render the observation primary field HTML
function createObsPrimaryField (obsLevels) {
    const tmpl = $.templates('#obs_primary_tmpl'); // Get compiled template using jQuery selector for the script block
    const html = tmpl.render(obsLevels); // Render template using data - as HTML string
    $('#obs_primary_container').html(html); // Insert HTML string into DOM
}

// Render the observation secondary field
function createObsSecondaryField (obsLevels) {
    const tmpl = $.templates('#obs_secondary_tmpl');
    const html = tmpl.render(obsLevels);
    $('#obs_secondary_container').html(html);
}

// Render the observation filter dropdowns
function createObsFilterDropdowns (obsLevels) {
    const tmpl = $.templates('#obs_filters_tmpl');
    const html = tmpl.render(obsLevels);
    $('#obs_filters_container').html(html);
    $('select.js-obs-levels').select2({
        placeholder: 'Start typing to include groups from this category. Click "All" to use all groups',
        allowClear: true,
        width: 'resolve'
    });
}

// Render clusterbars for heatmap
function createObsClusterBarField (obsLevels) {
    const tmpl = $.templates('#obs_clusterbar_tmpl');
    const html = tmpl.render(obsLevels);
    $('#obs_clusterbar_container').html(html);
}

// Render the sortable list for the chosen category
function createObsSortable (obsLevel, scope) {
    const escapedObsLevel = $.escapeSelector(obsLevel);
    const propData = $(`#${escapedObsLevel}_dropdown`).select2('data');
    const sortData = propData.map((elem) => elem.id);
    const tmpl = $.templates('#obs_sortable_tmpl');
    const html = tmpl.render({ sortData });
    $(`#${scope}_sortable`).html(html);
    $(`#${scope}_sortable`).sortable();

    // Store category name to retrieve later.
    sortCategories[scope] = obsLevel;

    $(`#${scope}_order_label`).show();
}

// Render dropdowns specific to the dot plot
function createDotplotDropdowns (obsLevels) {
    createObsPrimaryField (obsLevels);
    $('#none_primary_group').hide();  // Primary is required
    createObsSecondaryField(obsLevels);
}

// Render dropdowns specific to the heatmap plot
function createHeatmapDropdowns (obsLevels) {
    createObsPrimaryField(obsLevels);
    createObsSecondaryField(obsLevels);
    createObsClusterBarField(obsLevels);

    // Initialize differential expression test dropdown
    $('#distance_select').select2({
        placeholder: 'Choose distance metric',
        width: '25%'
    });

    $('#cluster_obs_warning').hide();
}

// Render dropdowns specific to the quadrant plot
function createQuadrantDropdowns (obsLevels) {
    // Filter only categories with at least 3 conditions
    const goodObsLevels = {}
    for (const category in obsLevels) {
        if (obsLevels[category].length >= 3) {
            goodObsLevels[category] = obsLevels[category]
        }
    }

    const tmpl = $.templates('#select_conditions_tmpl');
    const html = tmpl.render(goodObsLevels);
    $('#quadrant_compare1_condition').html(html);
    $('#quadrant_compare1_condition').select2({
        placeholder: 'Select the first query condition.',
        width: '25%'
    });
    $('#quadrant_compare2_condition').html(html);
    $('#quadrant_compare2_condition').select2({
        placeholder: 'Select the second query condition.',
        width: '25%'
    });
    $('#quadrant_ref_condition').html(html);
    $('#quadrant_ref_condition').select2({
        placeholder: 'Select the reference condition.',
        width: '25%'
    });

    // Initialize differential expression test dropdown
    $('#de_test_select').select2({
        placeholder: 'Choose DE testing algorithm',
        width: '25%'
    });
}

// Render dropdowns specific to the violin plot
function createViolinDropdowns (obsLevels) {
    createObsPrimaryField (obsLevels);
    $('#none_primary_group').hide();  // Primary is required
    createObsSecondaryField(obsLevels);
}

// Render dropdowns specific to the volcano plot
function createVolcanoDropdowns (obsLevels) {
    // Filter only categories wtih at least 2 conditions
    const goodObsLevels = {}
    for (const category in obsLevels) {
        if (obsLevels[category].length >= 2) {
            goodObsLevels[category] = obsLevels[category]
        }
    }

    const tmpl = $.templates('#select_conditions_tmpl');
    const html = tmpl.render(goodObsLevels);
    $('#volcano_query_condition').html(html);
    $('#volcano_query_condition').select2({
        placeholder: 'Select the query condition.',
        width: '25%'
    });

    for (const category in goodObsLevels) {
        // Add the "union" option if the category has 3 or more groups.
        // "Union" in 2-group categories would just be the same as choosing the 2nd category
        if (goodObsLevels[category].length > 2) {
            goodObsLevels[category].push("Union of the rest of the groups");
        }
    }
    const html2 = tmpl.render(goodObsLevels);
    $('#volcano_ref_condition').html(html2);
    $('#volcano_ref_condition').select2({
        placeholder: 'Select the reference condition.',
        width: '25%'
    });

    // Initialize differential expression test dropdown
    $('#de_test_select').select2({
        placeholder: 'Choose DE testing algorithm',
        width: '25%'
    });

}

function curateObservations (obsLevels) {
    const obsNotUsed = [];  // Observations that will not be added to obsFilters due to issues (1 group, meaningless category)

    // Delete useless filters
    for (const property in obsLevels) {
        if (property === 'color' || property.endsWith('colors')) {
            obsNotUsed.push(property);
            delete obsLevels[property];
        } else if (obsLevels[property].length === 1) {
            obsNotUsed.push(property);
            delete obsLevels[property];
        }
    }
    return [obsLevels, obsNotUsed];
}

// Generate a list of saved plot displays the user has access in viewing.
async function loadSavedDisplays (datasetId, defaultDisplayId=null) {
    const datasetData = await fetchDatasetInfo(datasetId);
    const { owner_id: ownerId } = datasetData;
    const userDisplays = await fetchUserDisplays(CURRENT_USER.id, datasetId);
    // Do not duplicate user displays in the owner display area as it can cause HTML element issues
    const ownerDisplays = CURRENT_USER.id === ownerId ? [] : await fetchOwnerDisplays(ownerId, datasetId);

    // Filter displays to those only with multigene plot types
    const mgUserDisplays = userDisplays.filter(d => plotTypes.includes(d.plot_type));
    const mgOwnerDisplays = ownerDisplays.filter(d => plotTypes.includes(d.plot_type));

    //
    mgUserDisplays.forEach(display => {
        display.is_default = false;
        if (defaultDisplayId && display.id === defaultDisplayId ) {
            display.is_default = true;
        }
        drawPreviewImage(display);
    });
    mgOwnerDisplays.forEach(display => {
        display.is_default = false;
        if (defaultDisplayId && display.id === defaultDisplayId ) {
            display.is_default = true;
        }
        drawPreviewImage(display);
    });

    const displaysTmpl = $.templates('#saved_display_modal_tmpl');
    const displaysHtml = displaysTmpl.render({
        dataset_id: datasetId,
        user_displays: mgUserDisplays,
        owner_displays: mgOwnerDisplays
    });

    $('#saved_display_modal').html(displaysHtml);
}

// Populate the HTML config options based on what was in the plot
function loadDisplayConfigHtml (plotConfig) {
    // NOTE: The calling function also clicks "#reset_opts", so the options are rendered already
    // Populate filter-by dropdowns
    obsFilters = plotConfig.obs_filters;
    for (const property in obsFilters) {
        const escapedProperty = $.escapeSelector(property);
        $(`#${escapedProperty}_dropdown`).val(obsFilters[property]);
        $(`#${escapedProperty}_dropdown`).trigger('change');
    }

    const escapedPrimary = $.escapeSelector(plotConfig.primary_col);
    const escapedSecondary = $.escapeSelector(plotConfig.secondary_col);

    $(`#${escapedPrimary}_primary`).prop('checked', true).click();
    $(`#${escapedSecondary}_secondary`).prop('checked', true).click();

    $('#plot_title').val(plotConfig.plot_title);

    // Populate plot type-specific dropdowns and checkbox options
    switch ($('#plot_type_select').val()) {
    case 'dotplot':
        break;
    case 'heatmap':
        for (const field in plotConfig.clusterbar_fields) {
            const escapedField = $.escapeSelector(field);
            $(`#${escapedField}_clusterbar`).prop('checked', true).click();
        }
        $('#matrixplot').prop('checked', plotConfig.matrixplot);
        $('#center_around_zero').prop('checked', plotConfig.center_around_zero);
        $('#cluster_obs').prop('checked', plotConfig.cluster_obs);
        $('#cluster_genes').prop('checked', plotConfig.cluster_obs);
        $('#flip_axes').prop('checked', plotConfig.flip_axes);
        $('#distance_select').val(plotConfig.distance_metric);
        $('#distance_select').trigger('change');
        break;
    case 'mg_violin':
        $('#stacked_violin').prop('checked', plotConfig.stacked_violin);
        $('#violin_add_points').prop('checked', plotConfig.violin_add_points);
        break;
    case 'quadrant':
        $('#include_zero_foldchange').prop('checked', plotConfig.include_zero_fc);
        $("#quadrant_foldchange_cutoff").val(plotConfig.fold_change_cutoff);
        $("#quadrant_fdr_cutoff").val(plotConfig.fdr_cutoff);
        $('#quadrant_compare1_condition').val(plotConfig.compare1_condition);
        $('#quadrant_compare1_condition').trigger('change');
        $('#quadrant_compare2_condition').val(plotConfig.compare2_condition);
        $('#quadrant_compare2_condition').trigger('change');
        $('#quadrant_ref_condition').val(plotConfig.ref_condition);
        $('#quadrant_ref_condition').trigger('change');
        $('#de_test_select').val(plotConfig.de_test_algo);
        $('#de_test_select').trigger('change');
        break;
    default:
        // volcano
        $('#volcano_pvalue_threshold').val(plotConfig.pvalue_threshold);
        $('#volcano_lower_logfc_threshold').val(plotConfig.lower_logfc_threshold);
        $('#volcano_upper_logfc_threshold').val(plotConfig.upper_logfc_threshold);
        $('#adj_pvals').prop('checked', plotConfig.adj_pvals);
        $('#annot_nonsig').prop('checked', plotConfig.annotate_nonsignificant)
        $('#volcano_query_condition').val(plotConfig.query_condition);
        $('#volcano_query_condition').trigger('change');
        $('#volcano_ref_condition').val(plotConfig.ref_condition);
        $('#volcano_ref_condition').trigger('change');
        $('#de_test_select').val(plotConfig.de_test_algo);
        $('#de_test_select').trigger('change');
    }
}

// Load all saved gene carts for the current user
function loadGeneCarts () {
    const d = new $.Deferred();

    if (!session_id) {
        // User is not logged in. Hide gene carts container
        $('#gene_cart_container').hide();
        d.resolve();
    } else {
        $.ajax({
            url: './cgi/get_user_gene_carts.cgi',
            type: 'post',
            data: { session_id },
            dataType: 'json',
            success(data, _textStatus, _jqXHR) { // source https://stackoverflow.com/a/20915207/2900840
                const carts = {};
                const cartTypes = ['domain', 'user', 'group', 'shared', 'public'];
                let cartsFound = false;

                for (const ctype of cartTypes) {
                    carts[ctype] = [];

                    if (data[`${ctype}_carts`].length > 0) {
                        cartsFound = true;

                        //User has some profiles
                        $.each(data[`${ctype}_carts`], (_i, item) => {
                            carts[ctype].push({value: item.id, text: item.label });
                        });
                    }
                }

                geneCartTree.domainGeneCarts = carts.domain;
                geneCartTree.userGeneCarts = carts.user;
                geneCartTree.groupGeneCarts = carts.group;
                geneCartTree.sharedGeneCarts = carts.shared;
                geneCartTree.publicGeneCarts = carts.public;
                geneCartTree.generateTree();

                if (! cartsFound ) {
                    $('#gene_cart_container').show();
                }

                d.resolve();
            },
            error(_jqXHR, _errorThrown) {
                // display_error_bar(jqXHR.status + ' ' + errorThrown.name);
                d.fail();
            }
        });
    }
    d.promise();
}

function saveGeneCart () {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: CURRENT_USER.session_id
        , label: $("#gene_cart_name").val()
        , gctype: "unweighted-list"
        , organism_id:  $("#dataset").data('organism-id')
        , is_public: 0
    });

    selectedGenes.forEach((sg) => {
        const gene = new Gene({
            id: sg.gene_id, // TODO: prop never defined... could make = gene_symbol
            gene_symbol: sg.gene_symbol,
        });
        gc.add_gene(gene);
    });

    gc.save(updateUIAfterGeneCartSaveSuccess, updateUIAfterGeneCartSaveFailure);
}

function updateUIAfterGeneCartSaveSuccess(gc) {
	$("#saved_gene_cart_info_c > h3").html(`Cart: ${gc.label}`);
    $('#saved_gene_cart_info_c > h3').addClass('text-success');
	$("#gene_cart_member_count").html(gc.genes.length);
	$("#saved_gene_cart_info_c").show();
}

function updateUIAfterGeneCartSaveFailure(gc) {
    $('#saved_gene_cart_info_c > h3').text('Issue with saving gene cart.');
    $('#saved_gene_cart_info_c > h3').addClass('text-danger');
    $('#saved_gene_cart_info_c').show();
}

function fetchDatasetInfo (datasetId) {
    return $.ajax({
        url: './cgi/get_dataset_info.cgi',
        type: 'POST',
        data: { dataset_id: datasetId },
        dataType: 'json'
    });
}

function fetchUserDisplays (userId, datasetId) {
    return $.ajax({
        url: './cgi/get_dataset_displays.cgi',
        type: 'POST',
        data: { user_id: userId, dataset_id: datasetId },
        dataType: 'json'
    });
}
function fetchOwnerDisplays (ownerId, datasetId) {
    return $.ajax({
        url: './cgi/get_dataset_displays.cgi',
        type: 'POST',
        data: { user_id: ownerId, dataset_id: datasetId },
        dataType: 'json'
    });
}

function getDefaultDisplay (datasetId) {
    return $.ajax({
        url: './cgi/get_default_display.cgi',
        type: 'POST',
        data: {
            user_id: CURRENT_USER.id,
            dataset_id: datasetId,
            is_multigene: 1
        },
        dataType: 'json'
    });
}

function downloadSelectedGenes() {
	// Builds a file in memory for the user to download.  Completely client-side.
	// plot_data contains three keys: x, y and symbols
	// build the file string from this

    const plotType = plotConfig.plot_type;

    // Adjust headers to the plot type
    let xLabel;
    let yLabel;

    if (plotType === "quadrant") {
        const query1 = plotConfig.compare1_condition.split(';-;')[1];
        const query2 = plotConfig.compare2_condition.split(';-;')[1];
        const ref = plotConfig.ref_condition.split(';-;')[1];
        xLabel = `${query1} vs ${ref} Log2FC`;
        yLabel = `${query2} vs ${ref} Log2FC`;
    } else {
        const query = plotConfig.query_condition.split(';-;')[1];
        let ref = plotConfig.ref_condition.split(';-;')[1];

        ref = ref === "Union of the rest of the groups" ? "rest" : ref;

        // volcano
        xLabel = `${query} vs ${ref} Log2FC`;
        yLabel = `${query} vs ${ref} p-value`;
    }
	let fileContents = `gene_symbol\t${xLabel}\t${yLabel}\n`;

    // Entering genes and info now.
    selectedGenes.forEach((gene) => {
        fileContents +=
            `${gene.gene_symbol}\t`
            + `${gene.x}\t`
            + `${gene.y}\n`
	});

	const element = document.createElement("a");
	element.setAttribute(
		"href",
		`data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(fileContents)}`
	);
	element.setAttribute("download", "selected_genes.tsv");
	element.style.display = "none";
	document.body.appendChild(element);
	element.click();
	document.body.removeChild(element);
}

$('#dataset').change(async function () {
    datasetId = $('#dataset').val();
    displayId = null;

    $('#post_dataset_spinner').show();

    // Obtain default display ID for this dataset
    const {default_display_id: defaultDisplayId} = await getDefaultDisplay(datasetId);

    // Populate saved displays modal
    loadSavedDisplays(datasetId, defaultDisplayId);

    $('#load_saved_plots').show();
    $('#plot_type_container').show();

    $('#gene_container').show();
    $('#gene_spinner').show();
    $('#gene_cart_clear').click();
    $('#categories_not_used').hide();
    $('#no_observations_error').hide();

    // Create promises to get genes and observations for this dataset
    const geneSymbolsPromise = fetchGeneSymbols({ datasetId, analysis: undefined })
    const h5adPromise =  fetchH5adInfo({ datasetId, analysis: undefined });

    // Execute both in parallel
    geneSymbols = await geneSymbolsPromise;
    const data = await h5adPromise;

    createGeneDropdown(geneSymbols);  // gene_spinner hidden here

    $('#post_dataset_spinner').hide();

    // Cannot cluster with just one gene (because function is only available
    // in dash.clustergram which requires 2 or more genes in plot)
    // Adding a gene will trigger a change to enable the property
    $("#cluster_obs").prop("disabled", true);
    $("#cluster_genes").prop("disabled", true);

    // Catch datasets with no observations (before filtering), and indicate this dataset cannot be used.
    if (! Object.keys(data.obs_levels).length) {
        $('#no_observations_error').show();
    }

    // Get categorical observations for this dataset
    [obsLevels, obsNotUsed] = curateObservations(data.obs_levels);
    numObs = data.num_obs;

    // Determine if at least one category has at least two groups (for volcanos) or at least three groups (for quadrants)
    // If condition is not met, disable the plot option
    let hasTwoGroupCategory = false;
    let hasThreeGroupCategory = false;

    for (const category in obsLevels) {
        if (obsLevels[category].length >= 3) {
            hasThreeGroupCategory = true;
        }
        if (obsLevels[category].length >= 2) {
            hasTwoGroupCategory = true;
        }
        if (hasTwoGroupCategory && hasThreeGroupCategory) {
            break;
        }
    }
    $('#volcano_opt').prop("disabled", false);
    if (!hasTwoGroupCategory) {
        $('#volcano_opt').prop("disabled", true);
    }
    $('#quadrant_opt').prop("disabled", false);
    if (!hasThreeGroupCategory) {
        $('#quadrant_opt').prop("disabled", true);
    }

    $('#create_plot').show();
    $("#create_plot").prop("disabled", true);
    $('#reset_opts').show();

    // If a plot type was already selected,
    // reset the options so configs are populated for the current dataset
    if ($('#plot_type_select').val() ) {
        $("#create_plot").prop("disabled", false);
        $('#reset_opts').click();
    }
});

$("#create_gene_cart").on("click", () => {
    $("#create_gene_cart_dialog").show("fade");
});

$("#cancel_save_gene_cart").on("click", () => {
    $("#create_gene_cart_dialog").hide("fade");
    $("#gene_cart_name").val("");
});

$("#gene_cart_name").on("input", function () {
    if ($(this).val() == "") {
        $("#save_gene_cart").prop("disabled", true);
    } else {
        $("#save_gene_cart").prop("disabled", false);
    }
});

$("#save_gene_cart").on("click", () => {
    $("#save_gene_cart").prop("disabled", true);

    if (CURRENT_USER) {
        saveGeneCart();
    } else {
        window.alert("You must be signed in to do that.");
    }
    $("#save_gene_cart").prop('disabled', false);
});

$("#download_plot").on("click", () => {
    Plotly.downloadImage(
        `dataset_${datasetId}_h5ad`, {
            width: $('#plot_download_width').val()
            , height: $('#plot_download_height').val()
            , scale: $('#plot_download_scale').val()
        }
    )
});

// Load user's gene carts
$('#gene_cart').change(function () {
    const geneCartId = $(this).val();
    const params = { session_id, gene_cart_id: geneCartId };
    const d = new $.Deferred(); // Causes editable to wait until results are returned
    // User is not logged in
    if (!session_id) {
        d.resolve();
    } else {
        // User is logged in

        $('#gene_spinner').show();
        // Get the gene cart members and populate the gene symbol search bar
        $.ajax({
            url: './cgi/get_gene_cart_members.cgi',
            type: 'post',
            data: params,
            success(data, _newValue, _oldValue) {
                if (data.success === 1) {
                    // Append gene symbols to search bar
                    const geneCartSymbols = [];

                    // format gene symbols into search string
                    $.each(data.gene_symbols, (_i, item) => {
                        geneCartSymbols.push(item.label);
                    });

                    const geneCartSymbolsLowerCase = geneCartSymbols.map(x => x.toLowerCase());
                    const geneSymbolsLowerCase = geneSymbols.map(x => x.toLowerCase());

                    // Get genes from gene cart that are present in dataset's genes.  Preserve casing of dataset's genes.
                    const intersection = geneSymbols.filter(x => geneCartSymbolsLowerCase.includes(x.toLowerCase()));

                    // Get genes from gene cart that are not in dataset's genes
                    const difference = geneCartSymbols.filter(x => !geneSymbolsLowerCase.includes(x.toLowerCase()));

                    $('#gene_dropdown').val(intersection);
                    $('#gene_dropdown').trigger('change');

                    $('#genes_not_found').hide();
                    if (difference.length > 0) {
                        const differenceString = difference.join(', ');
                        $('#genes_not_found').text(`The following gene cart genes were not found in this dataset: ${differenceString}`);
                        $('#genes_not_found').show();
                    }
                }
                d.resolve();
            }
        });
        $('#gene_spinner').hide();
    }
    return d.promise();
});

// Some options are specific to certain plot types
$('#plot_type_select').change(() => {
    $('#reset_opts').click();  // Reset all options
    $("#create_plot").prop("disabled", false);

    dotplotOptsIds.forEach(id => {
        $(id).hide();
    })
    heatmapOptsIds.forEach(id => {
        $(id).hide();
    });
    quadrantOptsIds.forEach(id => {
        $(id).hide();
    })
    violinOptsIds.forEach(id => {
        $(id).hide();
    });
    volcanoOptsIds.forEach(id => {
        $(id).hide();
    });

    switch ($('#plot_type_select').val()) {
    case 'dotplot':
        dotplotOptsIds.forEach(id => {
            $(id).show();
        })
        $("#gene_selection_help").text("Choose the genes to include in plot.");
        break;
    case 'heatmap':
        heatmapOptsIds.forEach(id => {
            $(id).show();
        });
        $("#gene_selection_help").text("Choose the genes to include in plot.");
        break;
    case 'mg_violin':
        violinOptsIds.forEach(id => {
            $(id).show();
        });
        $("#gene_selection_help").text("Choose the genes to include in plot.");
        break;
    case 'quadrant':
        quadrantOptsIds.forEach(id => {
            $(id).show();
        })
        $("#gene_selection_help").text("OPTIONAL: Gene selection is optional for this plot type. Selected genes are annotated in the plot.");
        break;
    default:
        // volcano
        volcanoOptsIds.forEach(id => {
            $(id).show();
        });
        $("#gene_selection_help").text("OPTIONAL: Gene selection is optional for this plot type. Selected genes are annotated in the plot.");
    }
});

$(document).on('change', '#cluster_obs', () => {
    $('#cluster_obs_warning').hide();
    if ($('#cluster_obs').is(':checked') && numObs >= 1000){
        $('#cluster_obs_warning').show();
    }
});

$(document).on('change', '#gene_dropdown', () => {
    const genesFilter = $('#gene_dropdown').select2('data').map((elem) => elem.id);
    // Show warning if too many genes are entered
    $("#too_many_genes_warning").hide();
    if (genesFilter.length > 10) {
        $("#too_many_genes_warning").text(`There are currently ${genesFilter.length} genes to be plotted. Be aware that with some plots, a high number of genes can make the plot congested or unreadable.`);
        $("#too_many_genes_warning").show();
    }

    // Cannot cluster columns with just one gene (because function is only available
    // in dash.clustergram which requires 2 or more genes in plot)
    if (genesFilter.length > 1) {
        $("#cluster_obs").prop("disabled", false);
        $("#cluster_genes").prop("disabled", false);
        return;
    }
    $("#cluster_obs").prop("disabled", true);
    $("#cluster_obs").prop("checked", false);
    $("#cluster_genes").prop("disabled", true);
    $("#cluster_genes").prop("checked", false);
});

// When a column is chosen, populate the sortable list
$(document).on('change', 'input[name="obs_primary"]', () => {
    const obsLevel = $('input[name="obs_primary"]:checked').val();
    if (obsLevel !== "none") {
        createObsSortable(obsLevel, "primary");
        $('#primary_order_label').show();
        $('#primary_order_help').show();
        return;
    }
    $('#primary_sortable').empty();
    $('#primary_order_label').hide();
    $('#primary_order_help').hide();
});

$(document).on('change', 'input[name="obs_secondary"]', () => {
    const obsLevel = $('input[name="obs_secondary"]:checked').val();
    if ( obsLevel !== "none") {
        createObsSortable(obsLevel, "secondary");
        $('#secondary_order_label').show();
        $('#secondary_order_help').show();
        return;
    }
    $('#secondary_sortable').empty();
    $('#secondary_order_label').hide();
    $('#secondary_order_help').hide();
});

// Determine if condition has no groups so the sort container will be disabled or not.
$(document).on('change', '.js-obs-levels', function () {
    const id = this.id; // This is not escaped
    const group = id.replace('_dropdown', '');
    const escapedGroup = $.escapeSelector(group);
    const propData = $(`#${escapedGroup}_dropdown`).select2('data');
    const props = propData.map((elem) => elem.id);
    $(`#${escapedGroup}_primary`).prop("disabled", false);
    $(`#${escapedGroup}_secondary`).prop("disabled", false);
    // Update sortables with current filters list
    if ($(`#${escapedGroup}_primary`).is(":checked")) {
        createObsSortable(group, "primary");

    }
    if ($(`#${escapedGroup}_secondary`).is(":checked")) {
        createObsSortable(group, "secondary");
    }
    // Disable sorting until it is known that filters have length
    if (!props.length) {
        $(`#${escapedGroup}_primary`).prop("disabled", true);
        $(`#${escapedGroup}_secondary`).prop("disabled", true);
    }
});

$(document).on('click', '#create_plot', async () => {

    // Reset plot errors and warnings for both plots
    $('.js-plot-error').empty().hide();
    $('.js-plot-warning').empty().hide();

    plotConfig = {};

    const plotType = $('#plot_type_select').select2('data')[0].id;
    plotConfig.plot_type = plotType;

    // Update filters based on selection
    obsFilters = {};
    for (const property in obsLevels) {
        const escapedProperty = $.escapeSelector(property);
        const propData = $(`#${escapedProperty}_dropdown`).select2('data');
        obsFilters[property] = propData.map((elem) => elem.id);

        // If no groups for an observation are selected, delete filter
        if (!obsFilters[property].length) {
            delete obsFilters[property];
        }
    }
    plotConfig.obs_filters = obsFilters;

    if (!plotType) {
        window.alert('Please select a plot type.');
        return;
    }

    plotConfig.gene_symbols = genesFilter = $('#gene_dropdown').select2('data').map((elem) => elem.id);

    const sortOrder = {};
    const categoriesUsed = [];
    // Grab the sorted order of the list and convert to array
    if (sortCategories.primary) {
        sortOrder[sortCategories.primary] = $('#primary_sortable').sortable("toArray", {attribute:"value"});
        categoriesUsed.push(sortCategories.primary);
    }
    // This should be rare, but just use the primary order if both are the same category
    if (sortCategories.secondary && !(categoriesUsed.includes(sortCategories.secondary))) {
        sortOrder[sortCategories.secondary] = $('#secondary_sortable').sortable("toArray", {attribute:"value"});
    }
    plotConfig.sort_order = sortOrder;

    if ($('input[name="obs_primary"]:checked').val() !== "none"){
        plotConfig.primary_col = $('input[name="obs_primary"]:checked').val();
    }
    if ($('input[name="obs_secondary"]:checked').val() !== "none"){
        plotConfig.secondary_col = $('input[name="obs_secondary"]:checked').val();
    }

    plotConfig.plot_title = $('#plot_title').val();

    // Add specific plotConfig options depending on plot type
    switch (plotType) {
    case 'dotplot':
        if ((plotConfig.gene_symbols).length < 1) {
            window.alert('At least one gene must be provided.');
            return;
        }
        if (!plotConfig.primary_col) {
            window.alert("Must select at least a primary category to aggregate groups by for dot plots.");
            return;
        }
        break;
    case 'heatmap':
        if ((plotConfig.gene_symbols).length < 2) {
            window.alert("Must select at least 2 genes to generate a heatmap");
            return;
        }
        plotConfig.clusterbar_fields = [];
        $('input[name="obs_clusterbar"]:checked').each( (idx, elem) => {
            plotConfig.clusterbar_fields.push($(elem).val());
        });
        plotConfig.matrixplot = $('#matrixplot').is(':checked');
        plotConfig.center_around_zero = $('#center_around_zero').is(':checked');
        plotConfig.cluster_obs = $('#cluster_obs').is(':checked');
        plotConfig.cluster_genes = $('#cluster_genes').is(':checked');
        plotConfig.flip_axes = $('#flip_axes').is(':checked');
        plotConfig.distance_metric = $('#distance_select').select2('data')[0].id;
        break;
    case 'mg_violin':
        if ((plotConfig.gene_symbols).length < 1) {
            window.alert('At least one gene must be provided.');
            return;
        }
        if (!plotConfig.primary_col) {
            window.alert("Must select at least a primary category to aggregate groups by for violin plots.");
            return;
        }
        plotConfig.stacked_violin = $('#stacked_violin').is(':checked');
        plotConfig.violin_add_points = $('#violin_add_points').is(':checked');
        break;
    case 'quadrant':
        plotConfig.include_zero_fc = $('#include_zero_foldchange').is(':checked');
        plotConfig.fold_change_cutoff = Number($("#quadrant_foldchange_cutoff").val());
        plotConfig.fdr_cutoff = Number($("#quadrant_fdr_cutoff").val());
        if (! $('#de_test_select').select2('data')[0].id) {
            window.alert('Must select a DE statistical test.');
            return;
        }
        plotConfig.de_test_algo = $('#de_test_select').select2('data')[0].id;
        plotConfig.compare1_condition = $('#quadrant_compare1_condition').select2('data')[0].id;
        plotConfig.compare2_condition = $('#quadrant_compare2_condition').select2('data')[0].id;
        plotConfig.ref_condition = $('#quadrant_ref_condition').select2('data')[0].id;
        if (!(plotConfig.compare1_condition && plotConfig.compare2_condition && plotConfig.ref_condition)) {
            window.alert('All comparision conditions must be selected to generate a quadrant plot.');
            return;
        }
        const condition1Key = plotConfig.compare1_condition.split(';-;')[0];
        const condition2Key = plotConfig.compare2_condition.split(';-;')[0];
        const refQuadrantKey = plotConfig.ref_condition.split(';-;')[0];

        const condition1Val = plotConfig.compare1_condition.split(';-;')[1];
        const condition2Val = plotConfig.compare2_condition.split(';-;')[1];
        const refQuadrantVal = plotConfig.ref_condition.split(';-;')[1];
        if ((condition1Key !== condition2Key) && (condition1Key !== refQuadrantKey)) {
            window.alert('Please choose 3 conditions from the same observation group.');
            return;
        }

        if ((condition1Val === refQuadrantVal) || (condition2Val === refQuadrantVal) || (condition1Val === condition2Val)) {
            window.alert('Please choose 3 different conditions.');
            return;
        }

        // If condition category was filtered, the selected groups must be present
        if (condition1Key in obsFilters
            && !(obsFilters[condition1Key].includes(condition1Val)
            && obsFilters[condition2Key].includes(condition2Val)
            && obsFilters[condition1Key].includes(refQuadrantVal))) {
            window.alert('One of the selected conditions is also chosen to be filtered out. Please adjust.');
            return;
        }

        break;
    default:
        // volcano
        plotConfig.pvalue_threshold = $('#volcano_pvalue_threshold').val();
        plotConfig.lower_logfc_threshold = $('#volcano_lower_logfc_threshold').val();
        plotConfig.upper_logfc_threshold = $('#volcano_upper_logfc_threshold').val();
        plotConfig.adj_pvals = $('#adj_pvals').is(':checked');
        plotConfig.annotate_nonsignificant = $('#annot_nonsig').is(':checked');
        if (! $('#de_test_select').select2('data')[0].id) {
            window.alert('Must select a DE statistical test.');
            return;
        }
        plotConfig.de_test_algo = $('#de_test_select').select2('data')[0].id;
        plotConfig.query_condition = $('#volcano_query_condition').select2('data')[0].id;
        plotConfig.ref_condition = $('#volcano_ref_condition').select2('data')[0].id;
        // Validation related to the conditions
        if (!(plotConfig.query_condition && plotConfig.ref_condition)) {
            window.alert('Both comparision conditions must be selected to generate a volcano plot.');
            return;
        }
        const queryKey = plotConfig.query_condition.split(';-;')[0];
        const refKey = plotConfig.ref_condition.split(';-;')[0];
        const queryVal = plotConfig.query_condition.split(';-;')[1];
        const refVal = plotConfig.ref_condition.split(';-;')[1];
        if (queryKey !== refKey) {
            window.alert('Please choose 2 conditions from the same observation group.');
            return;
        }

        if (queryVal === refVal) {
            window.alert('Please choose 2 different conditions.');
            return;
        }

        // If condition category was filtered, the selected groups must be present
        if (queryKey in obsFilters
            && !(obsFilters[queryKey].includes(queryVal)
            && obsFilters[queryKey].includes(refVal))
            && refVal !== "Union of the rest of the groups") {
            window.alert('One of the selected conditions is also chosen to be filtered out. Please adjust.');
            return;
        }
    }


    // Render dataset plot HTML
    const plotTemplate = $.templates('#dataset_plot_tmpl');
    const plotHtml = plotTemplate.render({ dataset_id: datasetId });
    $('#dataset_plot').html(plotHtml);


	$("#tbl_selected_genes").hide();
    if (["quadrant", "volcano"].includes(plotType)) {
        $('#dataset_plot').removeClass("col").addClass("col-9");
        $('#genes_list_bar').show();
    } else {
        $('#dataset_plot').addClass("col").removeClass("col-9");
        $('#genes_list_bar').hide();
    }

    // Draw the updated chart
    $('#plot_spinner').show();
    await draw(datasetId, plotConfig);
    $('#plot_spinner').hide();

    // Show plot options and disable selected genes button (since genes are not selected anymore)
    $('#post_plot_options').show();
});

// If "all" button is clicked, populate dropdown with all groups in this observation
$(document).on('click', '.js-all', function () {
    const id = this.id;
    const group = id.replace('_all', '');
    const escapedGroup = $.escapeSelector(group);

    $(`#${escapedGroup}_dropdown`).val(obsLevels[group]);
    $(`#${escapedGroup}_dropdown`).trigger('change'); // This actually triggers select2 to show the dropdown vals
});

// If "all" button is clicked, populate dropdown with all groups in this observation
$(document).on('click', '.js-clear', function () {
    const id = this.id;
    const group = id.replace('_clear', '');
    const escapedGroup = $.escapeSelector(group);

    $(`#${escapedGroup}_dropdown`).val('');
    $(`#${escapedGroup}_dropdown`).trigger('change'); // This actually triggers select2 to clear the dropdown vals
});

// Clear gene cart
$(document).on('click', '#gene_cart_clear', () => {
    $('#gene_dropdown').val('');  // Clear genes
    $('#gene_dropdown').trigger('change');
    $('#gene_cart').text('Choose gene cart'); // Reset gene cart text
    $('#gene_cart').val('');
    $('#genes_not_found').hide();
    $('#too_many_genes_warning').hide();
});

// Reset observation filters choices to be empty
$(document).on('click', '#reset_opts', async function () {
    $('#options_container').show();
    $('#options_spinner').show();

    // Reset sorting order
    sortCategories = {"primary": null, "secondary": null};
    $('#primary_order_label').hide();
    $('#primary_order_help').hide();
    $('#secondary_order_label').hide();
    $('#secondary_order_help').hide();

    // Update fields dependent on dataset observations
    createObsFilterDropdowns(obsLevels);
    $('#categories_not_used').hide();
    if (obsNotUsed && obsNotUsed.length) {
        const obsNotUsedString = obsNotUsed.join(", ")
        $('#categories_not_used').show();
        $('#categories_not_used').html(`The following observation categories were excluded either due to having one group or being considered to have no meaning: ${obsNotUsedString}`);
    }

    switch ($('#plot_type_select').val()) {
    case 'dotplot':
        createDotplotDropdowns(obsLevels);
        break;
    case 'heatmap':
        createHeatmapDropdowns(obsLevels);
        break;
    case 'mg_violin':
        createViolinDropdowns(obsLevels);
        break;
    case 'quadrant':
        createQuadrantDropdowns(obsLevels);
        break;
    default:
        // volcano
        createVolcanoDropdowns(obsLevels);
    }

    $('.js-all').click();  // Include all groups for every category (filter nothing)

    $('#options_spinner').hide();
});

// Save plot
$(document).on('click', '#save_display_btn', async function () {
    $('#saved_plot_confirmation').hide();
    $('#saved_plot_confirmation').removeClass('text-success');
    $('#saved_plot_confirmation').removeClass('text-danger');

    const plotType = $('#plot_type_select').select2('data')[0].id;

    const payload = {
        id: null, // Want to save as a new display
        dataset_id: datasetId,
        user_id: CURRENT_USER.id,
        label: $('#display_name').val(),
        plot_type: plotType,
        plotly_config: JSON.stringify({
            // depending on display type, this object will
            // have different properties
            ...plotConfig
        })
    };

    const res = await $.ajax({
        url: './cgi/save_dataset_display.cgi',
        type: 'POST',
        data: payload,
        dataType: 'json'
    });

    if (res?.success) {
        let msg = 'Plot successfully saved'

        if ($("#save_as_default_check").is(':checked') && res.display_id) {
            displayId = res.display_id;
            const res2 = await $.ajax({
                url: './cgi/save_default_display.cgi',
                type: 'POST',
                data: {
                    user_id: CURRENT_USER.id,
                    dataset_id: datasetId,
                    display_id: displayId,
                    is_multigene: 1
                },
                dataType: 'json'
            });

            if (res2?.success) {
                // Swap current default buttons
                $('.js-current-default')
                    .prop('disabled', false)
                    .addClass('js-save-default')
                    .addClass('btn-purple')
                    .removeClass('btn-secondary')
                    .removeClass('js-current-default')
                    .text("Make Default");
                $(`#${displayId}_default`)
                    .prop('disabled', true)
                    .removeClass('js-save-default')
                    .removeClass('btn-purple')
                    .addClass('btn-secondary')
                    .addClass('js-current-default')
                    .text("Default");
                msg += " and was set as the default display.";
            } else {
                msg += " but there was an issue saving as the default display.";
            }
        } else {
            msg += " but not set as default display."
        }

        $('#saved_plot_confirmation').text(msg);
        $('#saved_plot_confirmation').addClass('text-success');
    } else {
        $('#saved_plot_confirmation').text('There was an issue saving the plot');
        $('#saved_plot_confirmation').addClass('text-danger');
    }
    $('#saved_plot_confirmation').show();

    // Update saved displays modal so new plot is included
    loadSavedDisplays(datasetId);
});

// Load display information back into the curator page
$(document).on('click', '.js-load-display', async function () {
    const id = this.id;
    displayId = id.replace('_load', '');

    const display = await $.ajax({
        url: './cgi/get_dataset_display.cgi',
        type: 'POST',
        data: { display_id: displayId },
        dataType: 'json'
    });

    plotConfig = display.plotly_config;

    // Load plot type
    $('#plot_type_select').val(display.plot_type);
    $('#plot_type_select').trigger('change');

    // Load gene symbols
    geneSymbols = plotConfig.gene_symbols;
    $('#gene_dropdown').val(geneSymbols);
    $('#gene_dropdown').trigger('change');

    // Load config options
    loadDisplayConfigHtml(plotConfig);

    // Hide modal box
    $('#load_plots_modal').modal('hide');

    // Draw the updated chart
    $('#plot_spinner').show();
    const plotTemplate = $.templates('#dataset_plot_tmpl');
    const plotHtml = plotTemplate.render({ dataset_id: datasetId });
    $('#dataset_plot').html(plotHtml);

    if (["quadrant", "volcano"].includes(display.plot_type)) {
        $('#dataset_plot').removeClass("col").addClass("col-10");
        $('#genes_list_bar').show();
    } else {
        $('#dataset_plot').addClass("col").removeClass("col-10");
        $('#genes_list_bar').hide();
    }
    await draw(datasetId, plotConfig);
    $('#plot_spinner').hide();

    // Show plot options
    $('#post_plot_options').show();

});

// Delete user display
$(document).on('click', '.js-delete-display', async function () {
    $('#delete_display_confirmation').hide();
    $('#delete_display_confirmation').removeClass('alert-success');
    $('#delete_display_confirmation').removeClass('alert-danger');

    const id = this.id;
    const displayId = id.replace('_delete', '');

    const res = await $.ajax({
        url: './cgi/delete_dataset_display.cgi',
        type: 'POST',
        data: { id: displayId, user_id: CURRENT_USER.id },
        dataType: 'json'
    });

    $('#load_plots_modal').modal('hide');
    $('#delete_display_confirmation').show();
    if (res?.success) {
        $('#delete_display_confirmation_text').text('Display was successfully deleted.');
        $('#delete_display_confirmation').addClass('alert-success');
    } else {
        $('#delete_display_confirmation_text').text('There was an issue deleting the saved display.');
        $('#delete_display_confirmation').addClass('alert-danger');
    }

    // Update saved displays, now that display has been deleted
    loadSavedDisplays(datasetId);
});

// Save this particular display as the user's default display
$(document).on('click', '.js-save-default', async function () {
    $('#saved_default_confirmation').hide();
    $('#saved_default_confirmation').removeClass('alert-success');
    $('#saved_default_confirmation').removeClass('alert-danger');

    const id = this.id;
    const displayId = id.replace('_default', '');

    const res = await $.ajax({
        url: './cgi/save_default_display.cgi',
        type: 'POST',
        data: {
            user_id: CURRENT_USER.id,
            dataset_id: datasetId,
            display_id: displayId,
            is_multigene: 1
        },
        dataType: 'json'
    });

    // Swap current default buttons
    $('.js-current-default')
        .prop('disabled', false)
        .addClass('js-save-default')
        .addClass('btn-purple')
        .removeClass('btn-secondary')
        .removeClass('js-current-default')
        .text("Make Default")
    $(`#${displayId}_default`)
        .prop('disabled', true)
        .removeClass('js-save-default')
        .removeClass('btn-purple')
        .addClass('btn-secondary')
        .addClass('js-current-default')
        .text("Default")

    $('#saved_default_confirmation').show();
    if (res?.success) {
        $('#saved_default_confirmation_text').text('Display successfully saved as your new default.');
        $('#saved_default_confirmation').addClass('alert-success');
    } else {
        $('#saved_default_confirmation_text').text('There was an issue setting the default display.');
        $('#saved_default_confirmation').addClass('alert-danger');
    }
});
