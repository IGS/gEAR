"use strict";

let currentAnalysis = null;
let clickedMarkerGenes = new Set();
let enteredMarkerGenes = new Set();
let currentLabel = null;
let analysisLabels = new Set();

// TODO:  Check font sizes on all instruction blocks
// TODO:  Check if mitochrondrial QC actually returned anything
// TODO:  Complete work on limiting the gene count
// TODO:  Louvain options are escaping their box
// TODO:  Make sure all plotting buttons either disable or show something else while the compute runs


/**
 * Represents a dataset tree.
 *
 * @class
 * @constructor
 * @param {Object} options - The options for the dataset tree.
 * @param {HTMLElement} options.element - The element to render the dataset tree.
 * @param {HTMLElement} options.searchElement - The element for searching the dataset tree.
 * @param {Function} options.selectCallback - The callback function to be called when a dataset is selected.
 */
const datasetTree = new DatasetTree({
    element: document.getElementById("dataset-tree")
    , searchElement: document.getElementById("dataset-query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current-dataset-c").classList.remove("is-hidden");
        document.getElementById("current-dataset").textContent = e.node.title;
        document.getElementById("current-dataset-post").textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        datasetId = newDatasetId;

        // TODO: reset_workbench() here
    })
});

/**
 * Loads the dataset tree by fetching dataset information from the curator API.
 * Populates the userDatasets, sharedDatasets, and domainDatasets arrays with dataset information.
 * Generates the dataset tree using the generateTree method of the datasetTree object.
 * @throws {Error} If there is an error fetching the dataset information.
 */
const loadDatasetTree = async () => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const {data: datasetData} = await apiCallsMixin.fetchAllDatasets();

        let counter = 0;

        // Create data structure with  dataset information owned by the user
        if (datasetData.user.datasets.length > 0) {
            // User has some profiles
            for (const item of datasetData.user.datasets) {
                if (item) {
                    userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        // Next, add datasets shared with the user
        if (datasetData.shared_with_user.datasets.length > 0) {
            for (const item of datasetData.shared_with_user.datasets) {
                if (item) {
                    sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        // Now, add public datasets
        if (datasetData.public.datasets.length > 0) {
            for (const item of datasetData.public.datasets) {
                if (item) {
                    domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        datasetTree.userDatasets = userDatasets;
        datasetTree.sharedDatasets = sharedDatasets;
        datasetTree.domainDatasets = domainDatasets;
        datasetTree.generateTree();
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch datasets. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }

}

/**
 * Updates the manual marker gene entries based on the provided gene string.
 *
 * @param {string} geneString - The gene string to update the marker gene entries with.
 */
const updateManualMarkerGeneEntries = (geneString) => {
    enteredMarkerGenes = new Set();
    const geneSyms = geneString.split(',');

    for (let geneSym of geneSyms) {
        geneSym = geneSym.trim();
        if (geneSym) {
            enteredMarkerGenes.add(geneSym);
        }
    }

    const counterSet = new Set([...enteredMarkerGenes, ...clickedMarkerGenes]);
    document.querySelector(UI.markerGenesUniqueCountElt).textContent = counterSet.size;
    document.querySelector(UI.markerGenesEnteredCountElt).textContent = enteredMarkerGenes.size;
}

/**
 * Handles page-specific login UI updates (after the login event is triggered).
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are complete.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {
	document.getElementById("page-header-label").textContent = "Single Cell Workbench";

    for (const elt of document.querySelectorAll("#primary-nav .menu-list a.is-active")) {
        elt.classList.remove("is-active");
    }

    document.querySelector("a[tool='sc_workbench'").classList.add("is-active");


    const sessionId = CURRENT_USER.session_id;
    if (! sessionId ) {
        createToast("Not logged in so saving analyses is disabled.", "is-warning");
        // TODO: Other actions
    }

	try {
		await loadDatasetTree()
        // If brought here by the "gene search results" page, curate on the dataset ID that referred us
        const urlParams = new URLSearchParams(window.location.search);
        if (urlParams.has("dataset_id")) {
            const linkedDatasetId = urlParams.get("dataset_id");
            try {
                // find DatasetTree node and trigger "activate"
                const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
                foundNode.setActive(true);
                datasetId = linkedDatasetId;
            } catch (error) {
                createToast(`Dataset id ${linkedDatasetId} was not found as a public/private/shared dataset`);
                throw new Error(error);
            }
        }
	} catch (error) {
		logErrorInConsole(error);
	}

}

/* Event listenters for elements already loaded */

document.querySelector(UI.visualizeMarkerGenesBtnElt).addEventListener("click", async (event) => {
    await performMarkerGeneVisualization();
});

document.querySelector(UI.markerGenesManuallyEnteredElt).addEventListener("keyup", (event) => {
    const geneString = event.target.value;
    updateManualMarkerGeneEntries(geneString);
});

/* -------------------------------------------------------- */

window.onload=() => {
    $('[data-toggle="tooltip"]').tooltip()

    $('.tooltoggle').bootstrapToggle('disable');

    $('#asvg_flavor_tooltip').tooltip({
        placement: 'bottom',
        html: true,
        trigger: 'hover',
        // I cannot get this to work, even though I can in this pen:
        //  https://jsbin.com/gebiju/edit?html,js,output
        delay: { "show": 100, "hide": 2000 }
    });


    $( "#dataset_id" ).on('change', async () => {
        show_working("Loading dataset");

        $('.js-next-step').hide();

        if (currentAnalysis != null) {
            reset_workbench();
        }

        currentAnalysis = new Analysis({'dataset_id': $("#dataset_id").val(),
                                         'type': 'primary',
                                         'dataset_is_raw': true});

        $(".initial_instructions").hide();
        // Technically these could load asynchronously, but logically the progress logs make more sense sequentially
        await get_dataset_info($("#dataset_id").val());
        load_preliminary_figures($("#dataset_id").val());

    });

    // This series traps the <enter> key within each form and fires a click event on
    //  the corresponding button instead.
    catch_enter_and_instead_click('primary_filter_options_f', 'btn_apply_primary_filter');
    catch_enter_and_instead_click('qc_mito_options_f', 'btn_do_analysis_qc_by_mito');
    catch_enter_and_instead_click('asvg_options_f', 'btn_do_analysis_select_variable_genes');
    catch_enter_and_instead_click('pca_options_f', 'btn_pca_run');
    catch_enter_and_instead_click('tsne_options_f', 'btn_tsne_run');
    catch_enter_and_instead_click('louvain_options_f', 'btn_louvain_run');
    catch_enter_and_instead_click('marker_genes_options_f', 'btn_marker_genes_run');

    $("#btn_apply_primary_filter").on('click', function() {
        apply_primary_filter();
    });

    $("#btn_asvg_save").on('click', function() {
        run_analysis_select_variable_genes(1);
        $("#btn_do_analysis_select_variable_genes").hide();
        $(this).text('Saved');
        $(this).prop("disabled", true);
    });

    $("#btn_compare_genes_run").on('click', function() {
        check_dependencies_and_run(run_analysis_compare_genes);
    });

    $("#btn_download_marker_genes").on('click', function() {
        download_marker_genes_table();
    });

    $(".btn_delete_analysis").on('click', function() {
        currentAnalysis.delete();
    });

    $("#btn_do_analysis_qc_by_mito").on('click', function() {
        run_analysis_qc_by_mito(0);
    });

    $("#btn_do_analysis_select_variable_genes").on('click', function() {
        run_analysis_select_variable_genes(0);
    });

    $("#btn_louvain_run").on('click', function() {
        check_dependencies_and_run(run_analysis_louvain);
    });

    $("#btn_make_public_copy").on('click', function() {
        currentAnalysis.make_public_copy();
    });

    $("#btn_marker_genes_run").on('click', function() {
        check_dependencies_and_run(run_analysis_marker_genes);
    });

    $("#btn_pca_run").on('click', function() {
        check_dependencies_and_run(run_analysis_pca);
    });

    $("#btn_pca_top_genes").on('click', function() {
       check_dependencies_and_run(run_analysis_pca_top_genes);
    });

    $("#btn_louvain_rerun_with_groups").on('click', function() {
        $('#marker_genes_group_labels td.group_user_label input').removeClass('duplicate');
        if ($("#louvain_merge_clusters").is(':checked') ) {
            check_dependencies_and_run(run_analysis_louvain);
            return;
        }

        const duplicate_count = currentAnalysis.marker_genes.count_and_highlight_duplicates();
        if (duplicate_count == 0) {
            check_dependencies_and_run(run_analysis_louvain);
        }
    });

    $("#btn_qbm_save").on('click', function() {
        run_analysis_qc_by_mito(1);
        $("#btn_do_analysis_qc_by_mito").hide();
    });

    $("#btn_save_analysis").on('click', function() {
        $(this).prop("disabled", true);
        $(this).text('Saving ...');
        currentAnalysis.save_to_user_area();
    });

    $("#btn_tsne_run").on('click', function() {
        check_dependencies_and_run(run_analysis_tsne);
    });

    $('#cg_download_table_f').on('click', () => {
        const qry_id = $('#query_cluster').val();
        const ref_id = $('#reference_cluster').val();

        download_table_as_excel('compare_genes_table_f',
                                `cluster_comparison_${qry_id}_vs_` +
                                ref_id + '.xls');
    });

    $('#cg_download_table_r').on('click', () => {
        const qry_id = $('#query_cluster').val();
        const ref_id = $('#reference_cluster').val();

        download_table_as_excel('compare_genes_table_r',
                                `cluster_comparison_${ref_id}_vs_` +
                                qry_id + '.xls');
    });

    $("#cg_show_table_f").on('click', () => {
        if ($("#compare_genes_table_f").is(":visible")) {
            $("#compare_genes_table_f").hide(500);
            $("#cg_show_table_f").html(' Show table');
            $("#cg_show_table_f").removeClass('fa-eye-slash');
            $("#cg_show_table_f").addClass('fa-eye');
            return;
        }
        $("#compare_genes_table_f").show(500);
        $("#cg_show_table_f").html(' Hide table');
        $("#cg_show_table_f").addClass('fa-eye-slash');
        $("#cg_show_table_f").removeClass('fa-eye');
    });

    $("#cg_show_table_r").on('click', () => {
        if ($("#compare_genes_table_r").is(":visible")) {
            $("#compare_genes_table_r").hide(500);
            $("#cg_show_table_r").html(' Show table');
            $("#cg_show_table_r").removeClass('fa-eye-slash');
            $("#cg_show_table_r").addClass('fa-eye');
            return;
        }
        $("#compare_genes_table_r").show(500);
        $("#cg_show_table_r").html(' Hide table');
        $("#cg_show_table_r").addClass('fa-eye-slash');
        $("#cg_show_table_r").removeClass('fa-eye');
    });

    $('.tooltoggle').change( function() {
        const analysis_block_id = `#analysis_${$(this).data('analysis-name')}`;

        if ( $(this).prop('checked') ) {
            $(analysis_block_id).show();
        } else {
            $(analysis_block_id).hide();
        }
    });

    $('select#analysis_id').on('change', function() {
        // The first analysis ID is the blank 'New' one
        if ($(this).find(':selected').data('analysis-id') == "0") {
            reset_workbench();

            $("#primary_analysis_notification").hide();
            $("#analysis_action_c").hide();
            $("#analysis_status_info").text("");
            $("#analysis_status_info_c").hide();
            $("#btn_make_public_copy").hide();
            $("#btn_delete_saved_analysis").hide();
            $("#btn_delete_unsaved_analysis").hide();

            currentAnalysis = new Analysis({'dataset_id': $("#dataset_id").val(),
            'type': 'primary',
            'dataset_is_raw': true});
            // Need to reload prelim step so qc_by_mito toggle will work
            $("#dataset_info").show();
            load_preliminary_figures($("#dataset_id").val());

            $('#analysis_id').selectpicker('refresh');
            $("#stored_analyses_c").show(10);

            return;
        }
        show_working("Loading stored analysis");

        $('#new_analysis_label_c').hide();
        reset_workbench();

        const analysis_type = $(this).find(':selected').data('analysis-type');

        load_stored_analysis($(this).find(':selected').data('analysis-id'),
            analysis_type,
            $(this).find(':selected').data('dataset-id')
        );

        if (analysis_type == 'user_unsaved') {
            $("#primary_analysis_notification").hide();
            $("#analysis_action_c").show();
            $("#analysis_status_info_c").hide();
            $("#btn_make_public_copy").hide();
            $("#btn_delete_saved_analysis").hide();
            $("#btn_delete_unsaved_analysis").show();
            return;
        }
        if (analysis_type == 'user_saved') {
            $("#primary_analysis_notification").hide();
            $("#analysis_action_c").hide();
            $("#analysis_status_info").text("This analysis is stored in your profile.");
            $("#analysis_status_info_c").show();
            $("#btn_make_public_copy").show();
            $("#btn_delete_saved_analysis").show();
            $("#btn_delete_unsaved_analysis").hide();
            return;
        }
        if (analysis_type == 'public') {
            $("#primary_analysis_notification").hide();
            $("#analysis_action_c").hide();
            $("#analysis_status_info").text("Changes made to this public analysis will spawn a local copy within your profile.");
            $("#analysis_status_info_c").show();
            $("#btn_make_public_copy").hide();
            $("#btn_delete_saved_analysis").hide();
            $("#btn_delete_unsaved_analysis").hide();
            return;
        }
        if (analysis_type == 'primary') {
            $("#primary_analysis_notification").show();
            $("#analysis_action_c").hide();
            $("#analysis_status_info_c").hide();
            $("#btn_make_public_copy").hide();
            $("#btn_delete_saved_analysis").hide();
            $("#btn_delete_unsaved_analysis").hide();
        }
    });

    $('select#compare_genes_method').on('change', function() {
        /*
           p-value correction method. Used only for ‘t-test’, ‘t-test_overestim_var’,
           and ‘wilcoxon’ methods.

           The only other option is 'logreg'
        */
        if (this.value === 'logreg') {
            $("#compare_genes_corr_method_c").hide();
        } else {
            $("#compare_genes_corr_method_c").show();
        }
    });

    $("#marker_genes_table").mouseover(function(e) {
        $(this).css("cursor", "pointer");
    });

    $("#marker_genes_table").click(function(e) {
        let clicked_cell = $(e.target).closest("td");
        const goi = clicked_cell.text().trim(); // goi = genes of interest

        if ($(e.target).hasClass("js-col-idx")) {
            clicked_cell = $(e.target).closest("th");
        }

        // If row index was clicked, operate on whole row. Otherwise, just on individual cells.
        if (clicked_cell.hasClass('js-row-idx')) {
            const row_cells = clicked_cell.siblings();
            // note - jQuery map, not array.prototype map
            const gois = row_cells.map((i, el) => el.innerText.trim()).get();
            if (clicked_cell.hasClass('highlighted')) {
                // unhighlight all cells in row
                row_cells.removeClass("highlighted");
                clicked_cell.removeClass("highlighted");
                gois.forEach(el => {
                    clickedMarkerGenes.delete(el);
                    currentAnalysis.remove_gene_of_interest(el);
                });
            } else {
                row_cells.addClass("highlighted");
                clicked_cell.addClass("highlighted");
                gois.forEach(el => {
                    clickedMarkerGenes.add(el);
                    currentAnalysis.add_gene_of_interest(el);
                });
            }
        } else if (clicked_cell.hasClass("js-col-idx")){
            // get nth element in each row +1 (for the row idx)
            const clicked_header_idx = clicked_cell.index()+1;
            const col_cells = $('table tr td:nth-child(' + clicked_header_idx + ')');
            // note - jQuery map, not array.prototype map
            const gois = col_cells.map((i, el) => el.innerText.trim()).get();
            if (clicked_cell.hasClass('highlighted')) {
                // unhighlight all cells in row
                col_cells.removeClass("highlighted");
                clicked_cell.removeClass("highlighted");
                gois.forEach(el => {
                    clickedMarkerGenes.delete(el);
                    currentAnalysis.remove_gene_of_interest(el);
                });
            } else {
                col_cells.addClass("highlighted");
                clicked_cell.addClass("highlighted");
                gois.forEach(el => {
                    clickedMarkerGenes.add(el);
                    currentAnalysis.add_gene_of_interest(el);
                });
            }
        } else if (clicked_cell.hasClass('highlighted')) {
            clicked_cell.removeClass("highlighted");
            clickedMarkerGenes.delete(goi);
            currentAnalysis.remove_gene_of_interest(goi);
        } else {
            clicked_cell.addClass("highlighted");
            clickedMarkerGenes.add(goi);
            currentAnalysis.add_gene_of_interest(goi);
        }

        // Occasionally an empty string finds its way in here, which can throw off counts.
        clickedMarkerGenes.delete('');


        $('#marker_genes_selected_count').text(clickedMarkerGenes.size);
        const counter_set = new Set([...enteredMarkerGenes, ...clickedMarkerGenes]);
        $('#marker_genes_unique_count').text(counter_set.size);
    });

    $("#new_analysis_label_cancel").click(function(e) {
        // Set the label back to what it currently is, then hide
        $("#new_analysis_label").val(currentAnalysis.label);
        $("#new_analysis_label_c").hide(500);
    });

    $("#new_analysis_label_save").click(function(e) {
        currentAnalysis.label = $("#new_analysis_label").val();
        currentAnalysis.save();
        $("#new_analysis_label_c").hide(500);
    });

    $("#btn_sbs_tsne_run").click(async function(e) {
        /*
          When the user clicks to search a gene for tSNE display,
          generate a tSNE image based on the default layout.
        */
        e.preventDefault();
        $("#sbs_tsne_plot_c").empty();
        $("#sbs_tsne_gene_not_found").hide();
        $('#sbs_tsne_plot_c').show();

        show_working("Generating tSNE display");
        const ds = currentAnalysis.dataset;
        const dsp = new DatasetPanel({...ds});

        const data = await dsp.get_embedded_tsne_display(dsp.id);
        dsp.config = data.plotly_config;

        const img = document.createElement('img');
        img.className = 'img-fluid';

        get_tsne_image_data($("#sbs_tsne_gene_symbol").val(), dsp.config).then(
            data => {
                if (typeof data === 'object' || typeof data === "undefined") {
                    // If this is true, there was an error
                    $("#sbs_tsne_gene_not_found").show();
                } else {
                    // place image
                    img.src = `data:image/png;base64,${data}`
                    document.getElementById('sbs_tsne_plot_c').appendChild(img);
                }
            }
        );

        $("#btn_sbs_tsne_run").prop('disabled', false);

        done_working("tSNE visualized", false);
    });

    $(".show_analysis_renamer").click(function(e) {
        $("#new_analysis_label").val(currentAnalysis.label);
        $("#new_analysis_label_c").show(500);
    });

    // Handling all buttons manually
    $("button").click(function(e) {
        e.preventDefault();
    });


    $('#marker_genes_manually_entered').focus(function() {
        reset_manual_marker_gene_entries();
    });

    $('#marker_genes_manually_entered').change(function() {
        process_manual_marker_gene_entries();
    });


    $('#new_analysis_label').focus(function() {
        currentLabel = $(this).val();
    });
    $('#new_analysis_label').keyup(function() {
        if ( analysisLabels.has($(this).val()) ) {
            // it's also OK if the current value was what it was when the user started to edit
            if ( $(this).val() != currentLabel ) {
                // turn red, disable save, and show dup message
                $(this).addClass('duplicate');
                $('#new_analysis_label_save').prop("disabled", true);
                $('#duplicate_label_warning').show(500);
            }
            return;

        }
        // clear duplication message, remove red, and enable save
        $('#duplicate_label_warning').hide(500);
        $(this).removeClass('duplicate');
        $('#new_analysis_label_save').prop("disabled", false);
    });

	$("#save_weighted_gene_cart").on("click", () => {
		$("#save_weighted_gene_cart").prop("disabled", true);
		if (CURRENT_USER) {
			save_weighted_gene_cart();
		} else {
			alert("You must be signed in to do that.");
		}
	});

    $("#save_marker_gene_cart").on("click", () => {
		$("#save_marker_gene_cart").prop("disabled", true);
		if (CURRENT_USER) {
			save_marker_gene_cart();
		} else {
			alert("You must be signed in to do that.");
		}
	});

    // Create observer to watch if user changes (ie. successful login does not refresh page)
    // See: https://developer.mozilla.org/en-US/docs/Web/API/MutationObserver

    // But we need to wait for navigation_bar to load first (in common.js) so do some polling
    // See: https://stackoverflow.com/q/38881301

    // Select the node that will be observed for mutations
    const target_node = document.getElementById('loggedin_controls');
    const safer_node = document.getElementById("navigation_bar");   // Empty div until loaded
    // Create an observer instance linked to the callback function
    const observer = new MutationObserver(function(mutationList, observer) {
        if (target_node) {
            populate_dataset_selection();
            this.disconnect();  // Don't need to reload once the trees are updated
        }
    });
    // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
    // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
    observer.observe(target_node || safer_node , { attributes: true });
}

function save_marker_gene_cart() {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: CURRENT_USER.session_id,
        label: $("#marker_gene_cart_name").val(),
        gctype: 'unweighted-list',
        organism_id: $("#dataset_id").data('organism-id'),
        is_public: 0
    });

    currentAnalysis.genes_of_interest.forEach((gene_id) => {
        const gene = new Gene({
            //id: gene_id,    // TODO: figure out how to get ensembl ID for this
            gene_symbol: gene_id,
        });
        gc.add_gene(gene);
    });

    gc.save(update_ui_after_marker_gene_cart_save_success, update_ui_after_marker_gene_cart_save_failure);

}

function update_ui_after_marker_gene_cart_save_success(gc) {
	$("#saved_marker_gene_cart_info_c > p").html(`Cart "${gc.label}" successfully saved.`);
	$("#saved_marker_gene_cart_info_c > p").removeClass("text-danger").addClass("text-success");
	$("#saved_marker_gene_cart_info_c").show();
    done_working("Saved marker gene cart", false);
}

function update_ui_after_marker_gene_cart_save_failure(gc, message) {
	$("#saved_marker_gene_cart_info_c > p").html("There was an issue saving the marker gene cart.");
	$("#saved_marker_gene_cart_info_c > p").removeClass("text-success").addClass("text-danger");
	$("#saved_marker_gene_cart_info_c").show();
    report_error(`Error saving gene cart: ${gc.label}`);
    report_error(message);
    $('#save_marker_gene_cart').attr("disabled", false);
}

function save_weighted_gene_cart() {
    // Return PC data from Anndata object
    $.ajax({
        type: "POST",
        url: "./cgi/get_PCs_from_anndata.cgi",
        data: {
            'dataset_id': currentAnalysis.dataset_id,
            'analysis_id': currentAnalysis.id,
            'analysis_type': currentAnalysis.type,
            'session_id': currentAnalysis.user_session_id,
        },
        datatype: "json",
    }).done((data) => {
        if (! data.success) {
            report_error(`Error getting PCs for saving: ${data.msg}`);
            $('#save_weighted_gene_cart').attr("disabled", false);
            return;
        }

        const weight_labels = data.pc_data.columns;

        // must have access to USER_SESSION_ID
        const gc = new WeightedGeneCart({
            session_id: CURRENT_USER.session_id,
            label: $("#weighted_gene_cart_name").val(),
            gctype: 'weighted-list',
            organism_id: $("#dataset_id").data('organism-id'),
            is_public: 0
        }, weight_labels
        );

        data.pc_data.index.forEach((gene_id, i) => {
            const weights = data.pc_data.data[i];
            const gene = new WeightedGene({
                id: gene_id,
                gene_symbol: data.gene_symbols[i]
            }, weights
            );
            gc.add_gene(gene);
        });

        gc.save(update_ui_after_weighted_gene_cart_save_success, update_ui_after_weighted_gene_cart_save_failure);
    }).fail((xhr, status, msg) => {
        report_error(`Error saving gene cart: ${msg}`);
        $('#save_weighted_gene_cart').attr("disabled", false);
    });

}

function update_ui_after_weighted_gene_cart_save_success(gc) {
	$("#saved_weighted_gene_cart_info_c > p").html(`Cart "${gc.label}" successfully saved.`);
	$("#saved_weighted_gene_cart_info_c > p").removeClass("text-danger").addClass("text-success");
	$("#saved_weighted_gene_cart_info_c").show();
    done_working("Saved weighted gene cart", false);
}

function update_ui_after_weighted_gene_cart_save_failure(gc, message) {
	$("#saved_weighted_gene_cart_info_c > p").html("There was an issue saving the weighted gene cart.");
	$("#saved_weighted_gene_cart_info_c > p").removeClass("text-success").addClass("text-danger");
	$("#saved_weighted_gene_cart_info_c").show();
    report_error(`Error saving gene cart: ${gc.label}`);
    report_error(message);
    $('#save_weighted_gene_cart').attr("disabled", false);
}

function get_tsne_image_data(gene_symbol, config) {
    config.colorblind_mode = CURRENT_USER.colorblind_mode;
    // then craziness: https://stackoverflow.com/a/48980526
    return axios.post(`/api/plot/${currentAnalysis.dataset_id}/tsne`, {
        gene_symbol,
        analysis: currentAnalysis,
        colorize_legend_by: config.colorize_legend_by,
        plot_type: 'tsne_static',
        plot_by_group: config.plot_by_group,
        max_columns: config.max_columns,
        horizontal_legend: config.horizontal_legend,
        x_axis: config.x_axis,
        y_axis: config.y_axis,
        analysis_owner_id: currentAnalysis.user_id,
        colors: config.colors,
        colorblind_mode: config.colorblind_mode,
        // helps stop caching issues
        timestamp: new Date().getTime()
    }).then(response => {
        return response.data.image
    }); // end axios
}

function apply_primary_filter() {
    show_working("Applying dataset filters");
    $("#dataset_info .empty_on_change").empty();

    // If a user redoes the primary filters, it's assumed they want to do this on the primary
    //  datasource again.  This allows them to lessen the stringency of a filter.
    const original_analysis_type = currentAnalysis.type
    currentAnalysis.type = 'primary'

    if ($('#filter_cells_lt_n_genes_selected').is(":checked")) {
        currentAnalysis.primary_filter.filter_cells_lt_n_genes_selected = true
        currentAnalysis.primary_filter.filter_cells_lt_n_genes = $("#filter_cells_lt_n_genes").val();
    } else {
        currentAnalysis.primary_filter.filter_cells_lt_n_genes_selected = false
    }


    if ($('#filter_cells_gt_n_genes_selected').is(":checked")) {
        currentAnalysis.primary_filter.filter_cells_gt_n_genes_selected = true
        currentAnalysis.primary_filter.filter_cells_gt_n_genes = $("#filter_cells_gt_n_genes").val();
    } else {
        currentAnalysis.primary_filter.filter_cells_gt_n_genes_selected = false
    }

    if ($('#filter_genes_lt_n_cells_selected').is(":checked")) {
        currentAnalysis.primary_filter.filter_genes_lt_n_cells_selected = true
        currentAnalysis.primary_filter.filter_genes_lt_n_cells = $("#filter_genes_lt_n_cells").val();
    } else {
        currentAnalysis.primary_filter.filter_genes_lt_n_cells_selected = false
    }

    if ($('#filter_genes_gt_n_cells_selected').is(":checked")) {
        currentAnalysis.primary_filter.filter_genes_gt_n_cells_selected = true
        currentAnalysis.primary_filter.filter_genes_gt_n_cells = $("#filter_genes_gt_n_cells").val();
    } else {
        currentAnalysis.primary_filter.filter_genes_gt_n_cells_selected = false
    }

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_apply_primary_filter.cgi",
        data: {'analysis_id': currentAnalysis.id,
               'analysis_type': currentAnalysis.type,
               'dataset_id': currentAnalysis.dataset_id,
               'session_id': currentAnalysis.user_session_id,
               'using_primary_datasource': currentAnalysis.using_primary_datasource,
               'filter_cells_lt_n_genes': currentAnalysis.primary_filter.filter_cells_lt_n_genes,
               'filter_cells_gt_n_genes': currentAnalysis.primary_filter.filter_cells_gt_n_genes,
               'filter_genes_lt_n_cells': currentAnalysis.primary_filter.filter_genes_lt_n_cells,
               'filter_genes_gt_n_cells': currentAnalysis.primary_filter.filter_genes_gt_n_cells
              },
        dataType: "json",
        success: function(data) {
            currentAnalysis.primary_filter.filtered_gene_count = data['n_genes'];
            currentAnalysis.primary_filter.filtered_cell_count = data['n_obs'];

            if (data['success'] == 1) {
                $("#selected_dataset_shape_filtered").html(data['n_genes'] + " genes x " + data['n_obs'] + " obs");
                $("#selected_dataset_shape_filtered_c").show(500);

                if (original_analysis_type == 'primary') {
                    currentAnalysis.type = 'user_unsaved'
                } else {
                    currentAnalysis.type = original_analysis_type
                }

                currentAnalysis.primary_filter.calculated = true;
                currentAnalysis.primary_filter.update_ui(currentAnalysis);
                done_working("Data filters applied");
                $('#primary_filter_options_c .js-next-step').show();  // Show that next toggle can be clicked
            } else {
                if (data['n_genes'] == 0) {
                    report_error("Filter reduced genes to 0. Try less stringent cutoffs");
                } else if (data['n_obs'] == 0) {
                    report_error("Filter reduced genes to 0. Try less stringent cutoffs");
                } else {
                    report_error("There was an error filtering this dataset");
                }
            }

            $('#btn_apply_primary_filter').attr("disabled", false);
        },
        error: function(xhr, status, msg) {
            report_error("Error filtering dataset: " + msg);
            $('#btn_apply_primary_filter').attr("disabled", false);
        }
    });
}

function catch_enter_and_instead_click(form_id, button_id) {
    /*
      Each analysis tool has a series of form elements and at least one submit button.  This
      function allows a quick way to override the default form submission when the user
      hits <enter> and instead simulates the click of a chosen button.
      */
    $('#' + form_id).on('submit', function(e) {
        e.preventDefault();
    });

    $('#' + form_id).keypress(function (e) {
        if (e.which == 13) {
            $("#" + button_id).click();
        }
    });
}

function check_dependencies_and_run(callback, opts) {
    if (currentAnalysis.type == 'public') {
        currentAnalysis.copy_to_user_unsaved(callback, opts);
    } else {
        callback(opts);
    }
}

function done_working(msg, do_save) {
    $("#action_log li.working").remove();
    if (msg) {
        $("<li>" + msg + "</li>").prependTo("#action_log");
    }

    if (do_save == null || do_save == true) {
        if (currentAnalysis.type != 'primary') {
            if ( $('select#analysis_id').find(':selected').data('analysis-id') == "0") {
                currentAnalysis.label = $('#new_analysis_label').val();
            }
            currentAnalysis.save();

            if (currentAnalysis.type == 'user_unsaved') {
                $("#analysis_action_c").show();
                $("#analysis_status_info_c").hide();
            } else if (currentAnalysis.type == 'user_saved') {
                $("#analysis_action_c").hide();
                $("#analysis_status_info_c").show();
            }
        }
    }
}

function download_marker_genes_table() {
    /*
      This builds a file in-memory for the user to download which is a the marker genes table
      in tab-delimited form.
     */
    var file_contents = '';

    // Do the header row
    var row = [];
    $("table#marker_genes_table thead th").each(function(){
        if ( row.length == 0) {
            row.push('');
        } else {
            row.push($(this).text());
        }
    });
    file_contents += row.join("\t") + "\n";

    // Now all the other rows
    $("table#marker_genes_table tbody tr").each(function(){
        row = [];

        $(this).find('td').each(function(){
            row.push($(this).text());
        });

        file_contents += row.join("\t") + "\n";
    });

    const element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(file_contents));
    element.setAttribute('download', 'marker_genes.xls');
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

//TODO: move this into a generic utils module
async function get_dataset_info(dataset_id) {
    $("#stored_analyses_c").hide();

    await $.ajax({
        type: "POST",
        url: "./cgi/get_dataset_info.cgi",
        data: {'dataset_id': dataset_id, 'include_shape': 1},
        dataType: "json"
    }).done((data) => {
        const ds = new Dataset(data);

        $('#new_analysis_label').val(currentAnalysis.label);

        // TODO: make Analysis.dataset the replacement for Analysis.dataset_id
        currentAnalysis.dataset = ds;

        update_selected_dataset(ds);
        $("#dataset_info").show();
        $("#analysis_list_c").show();
        analysisLabels = currentAnalysis.get_saved_analyses_list(ds.id, 0, 'sc_workbench');
        done_working();
    }).fail((xhr, status, msg) => {
        report_error("Failed to access dataset");
        report_error(`Failed ID was: ${dataset_id} because msg: ${msg}`);
    });
}

function load_preliminary_figures(dataset_id) {
    $("#stored_analyses_c").hide();

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_preview_primary_filter.cgi",
        data: {'dataset_id': currentAnalysis.dataset_id, 'analysis_id': currentAnalysis.id,
               'analysis_type': currentAnalysis.type, 'session_id': currentAnalysis.user_session_id
              },
        dataType: "json"
    }).done((data) => {
        $("#primary_initial_plot_loading_c").hide();

        if (data['success'] == 1) {
            $('#primary_initial_violin_c').html(`<a target="_blank" href="./datasets/${dataset_id}.prelim_violin.png"><img src="./datasets/` + dataset_id + '.prelim_violin.png" class="img-fluid img-zoomed" /></a>');
            $('#primary_initial_scatter_c').html(`<a target="_blank" href="./datasets/${dataset_id}.prelim_n_genes.png"><img src="./datasets/` + dataset_id + '.prelim_n_genes.png" class="img-fluid img-zoomed" /></a>');
            done_working("Prelim plots displayed");
        } else {
            $('#primary_initial_violin_c').html('Preliminary plots not yet generated. Continue your analysis.');
            done_working("Prelim figures missing.");
        }

        $("#primary_initial_plot_c").show(500);
    }).fail((xhr, status, msg) => {
        report_error("Failed to access dataset");
        report_error(`Failed ID was: ${dataset_id} because msg: ${msg}`);
    });
}

function load_stored_analysis(analysis_id, analysis_type, dataset_id) {
    $.ajax({
        type: "POST",
        url: "./cgi/get_stored_analysis.cgi",
        data: {'analysis_id': analysis_id,
               'analysis_type': analysis_type,
               'session_id': CURRENT_USER.session_id,
               'dataset_id': dataset_id
              },
        dataType: "json",
        success: function(data) {
            const ana = Analysis.load_from_json(data);
            ana.dataset = currentAnalysis.dataset;
            currentAnalysis = ana;

            if (data['tsne']['tsne_calculated']) {
                $("#analysis_sbs_tsne").show();
            } else {
                $("#analysis_sbs_tsne").hide();
            }

            done_working("Analysis loaded", false);
        },
        error: function(xhr, status, msg) {
            report_error("Failed to load stored analysis: " + msg);
        }
    });
}

$(document).on("build_jstrees", () => populate_dataset_selection());


function process_manual_marker_gene_entries() {
    currentAnalysis.genes_of_interest = new Set([...enteredMarkerGenes, ...clickedMarkerGenes]);
    // Only allow saving of gene cart if genes are selected
    $("#save_marker_gene_cart").prop("disabled", true);
    if (currentAnalysis.genes_of_interest.size) {
		$("#save_marker_gene_cart").prop("disabled", false);

    }
}

function reset_manual_marker_gene_entries() {
    clickedMarkerGenes = new Set();

    // remember which GOI are from the table
    $.each(
        $('#marker_genes_table td.highlighted'),
        function() {
            clickedMarkerGenes.add($(this).text());
        }
    );
}

function reset_workbench() {
    // Performs all the steps needed to reset the workbench if the input dataset changes

    // handle any data structures
    currentAnalysis.reset();

    // Reset the analysis tool fields and some display labels
    $("form.reset_on_change").trigger("reset");
    $("p.reset_on_change").hide();
    $("span.empty_on_change").empty();
    $("div.empty_on_change").empty();
    $("tbody.empty_on_change").empty();
    $("input#new_analysis_label").val('');
    $('#top_genes strong').empty();

    // Hide any non-analysis-flow steps
    $("#analysis_sbs_tsne").hide();
    $("#group_labels_c").hide();

    // Toggle the tool buttons to hide them in the UI
    $('.tooltoggle').bootstrapToggle('off');

    // Disable those analysis steps which have previous requirements
    $('.tooltoggle').bootstrapToggle('disable');
}


