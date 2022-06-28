current_analysis = null;
var clicked_marker_genes = new Set();
var entered_marker_genes = new Set();
var current_label = null;
var analysis_labels = new Set();

// TODO:  Check font sizes on all instruction blocks
// TODO:  Check if mitochrondrial QC actually returned anything
// TODO:  Complete work on limiting the gene count
// TODO:  Louvain options are escaping their box
// TODO:  Make sure all plotting buttons either disable or show something else while the compute runs

const dataset_tree = new DatasetTree({treeDiv: '#dataset_tree'});

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


    $( "#dataset_id" ).on('change', () => {
        show_working("Loading dataset");

        $('.js-next-step').hide();

        if (current_analysis != null) {
            reset_workbench();
        }

        current_analysis = new Analysis({'dataset_id': $("#dataset_id").val(),
                                         'type': 'primary',
                                         'dataset_is_raw': true});

        $(".initial_instructions").hide();
        get_dataset_info($("#dataset_id").val());
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
        current_analysis.delete();
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
        current_analysis.make_public_copy();
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
        var duplicate_count = current_analysis.marker_genes.count_and_highlight_duplicates();
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
        current_analysis.save_to_user_area();
    });

    $("#btn_tsne_run").on('click', function() {
        check_dependencies_and_run(run_analysis_tsne);
    });

    $('#cg_download_table_f').on('click', function() {
        var qry_id = $('#query_cluster').val();
        var ref_id = $('#reference_cluster').val();

        download_table_as_excel('compare_genes_table_f',
                                'cluster_comparison_' + qry_id + '_vs_' +
                                ref_id + '.xls');
    });

    $('#cg_download_table_r').on('click', function() {
        var qry_id = $('#query_cluster').val();
        var ref_id = $('#reference_cluster').val();

        download_table_as_excel('compare_genes_table_r',
                                'cluster_comparison_' + ref_id + '_vs_' +
                                qry_id + '.xls');
    });

    $("#cg_show_table_f").on('click', function() {
        if ($("#compare_genes_table_f").is(":visible")) {
            $("#compare_genes_table_f").hide(500);
            $("#cg_show_table_f").html(' Show table');
            $("#cg_show_table_f").removeClass('fa-eye-slash');
            $("#cg_show_table_f").addClass('fa-eye');
        } else{
            $("#compare_genes_table_f").show(500);
            $("#cg_show_table_f").html(' Hide table');
            $("#cg_show_table_f").addClass('fa-eye-slash');
            $("#cg_show_table_f").removeClass('fa-eye');
        }
    });

    $("#cg_show_table_r").on('click', function() {
        if ($("#compare_genes_table_r").is(":visible")) {
            $("#compare_genes_table_r").hide(500);
            $("#cg_show_table_r").html(' Show table');
            $("#cg_show_table_r").removeClass('fa-eye-slash');
            $("#cg_show_table_r").addClass('fa-eye');
        } else{
            $("#compare_genes_table_r").show(500);
            $("#cg_show_table_r").html(' Hide table');
            $("#cg_show_table_r").addClass('fa-eye-slash');
            $("#cg_show_table_r").removeClass('fa-eye');
        }
    });

    $('.tooltoggle').change( function() {
        var analysis_block_id = '#analysis_' + $(this).data('analysis-name');

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

            $("#analysis_action_c").hide();
            $("#analysis_status_info").text("");
            $("#analysis_status_info_c").hide();
            $("#btn_make_public_copy").hide();
            $("#btn_delete_saved_analysis").hide();
            $("#btn_delete_unsaved_analysis").hide();
        } else {
            show_working("Loading stored analysis");

            $('#new_analysis_label_c').hide();
            reset_workbench();

            var analysis_type = $(this).find(':selected').data('analysis-type');

            load_stored_analysis($(this).find(':selected').data('analysis-id'),
                                 analysis_type,
                                 $(this).find(':selected').data('dataset-id'));

            if (analysis_type == 'user_unsaved') {
                $("#primary_analysis_notification").hide();
                $("#analysis_action_c").show();
                $("#analysis_status_info_c").hide();
                $("#btn_make_public_copy").hide();
                $("#btn_delete_saved_analysis").hide();
                $("#btn_delete_unsaved_analysis").show();
            } else if (analysis_type == 'user_saved') {
                $("#primary_analysis_notification").hide();
                $("#analysis_action_c").hide();
                $("#analysis_status_info").text("This analysis is stored in your profile.");
                $("#analysis_status_info_c").show();
                $("#btn_make_public_copy").show();
                $("#btn_delete_saved_analysis").show();
                $("#btn_delete_unsaved_analysis").hide();
            } else if (analysis_type == 'public') {
                $("#primary_analysis_notification").hide();
                $("#analysis_action_c").hide();
                $("#analysis_status_info").text("Changes made to this public analysis will spawn a local copy within your profile.");
                $("#analysis_status_info_c").show();
                $("#btn_make_public_copy").hide();
                $("#btn_delete_saved_analysis").hide();
                $("#btn_delete_unsaved_analysis").hide();
            } else if (analysis_type == 'primary') {
                $("#primary_analysis_notification").show();
                $("#analysis_action_c").hide();
                $("#analysis_status_info_c").hide();
                $("#btn_make_public_copy").hide();
                $("#btn_delete_saved_analysis").hide();
                $("#btn_delete_unsaved_analysis").hide();
            }
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
        const clicked_cell = $(e.target).closest("td");
        const goi = clicked_cell.text().trim();

        // If row index was clicked, operate on whole row. Otherwise, just on individual cells.
        if (clicked_cell.hasClass('js-row-idx')) {
            const row_cells = clicked_cell.siblings();
            // note - jQuery map, not array.prototype map
            const gois = row_cells.map((i, el) => el.innerText.trim()).get();
            if (clicked_cell.hasClass('highlighted')) {
                row_cells.removeClass("highlighted");
                clicked_cell.removeClass("highlighted");
                gois.forEach(el => {
                    clicked_marker_genes.delete(el);
                    current_analysis.remove_gene_of_interest(el);
                });
            } else {
                row_cells.addClass("highlighted");
                clicked_cell.addClass("highlighted");
                gois.forEach(el => {
                    clicked_marker_genes.add(el);
                    current_analysis.add_gene_of_interest(el);
                });
            }
        } else if (clicked_cell.hasClass('highlighted')) {
            clicked_cell.removeClass("highlighted");
            clicked_marker_genes.delete(goi);
            current_analysis.remove_gene_of_interest(goi);
        } else {
            clicked_cell.addClass("highlighted");
            clicked_marker_genes.add(goi);
            current_analysis.add_gene_of_interest(goi);
        }

        // Occasionally an empty string finds its way in here, which can throw off counts.
        clicked_marker_genes.delete('');


        $('#marker_genes_selected_count').text(clicked_marker_genes.size);
        const counter_set = new Set([...entered_marker_genes, ...clicked_marker_genes]);
        $('#marker_genes_unique_count').text(counter_set.size);
    });

    $("#new_analysis_label_cancel").click(function(e) {
        // Set the label back to what it currently is, then hide
        $("#new_analysis_label").val(current_analysis.label);
        $("#new_analysis_label_c").hide(500);
    });

    $("#new_analysis_label_save").click(function(e) {
        current_analysis.label = $("#new_analysis_label").val();
        current_analysis.save();
        $("#new_analysis_label_c").hide(500);
    });

    $("#btn_visualize_marker_genes").click(function(e) {
        $("#marker_genes_dotplot_c").empty();
        $("#marker_genes_violin_c").empty();
        show_working("Generating marker gene visualization");

        $.ajax({
            type: "POST",
            url: "./cgi/h5ad_generate_marker_gene_visualization.cgi",
            data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
                   'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
                   // This is necessary because stringify won't support a Set directly
                   'marker_genes': JSON.stringify([...current_analysis.genes_of_interest])
                  },
            dataType: "json",
            success: function(data) {
                if (data['success'] == 1) {
                    var params = {
                        'analysis_id': current_analysis.id,
                        'analysis_name': 'dotplot_goi',
                        'analysis_type': current_analysis.type,
                        'dataset_id': current_analysis.dataset_id,
                        'session_id': current_analysis.user_session_id,
                        // this saves the user from getting a cached image each time
                        datetime: (new Date()).getTime()
                    }

                    current_analysis.place_analysis_image(
                        {'params': params, 'title': 'Marker genes (dotplot)', 'target': '#marker_genes_dotplot_c'});

                    params['analysis_name'] = 'stacked_violin_goi'
                    current_analysis.place_analysis_image(
                        {'params': params, 'title': 'Marker genes (stacked violin)', 'target': '#marker_genes_violin_c'});

                    $("#marker_genes_results_c").show(500);

                    done_working("Marker genes visualized", false);
                } else {
                    done_working("Marker gene visualization failed.", false);
                }
            },
            error: function(xhr, status, msg) {
                report_error("Error visualizing marker genes");
            }
        });
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
        var ds = current_analysis.dataset;
        var dsp = new DatasetPanel({...ds});

        const data = await dsp.get_embedded_tsne_display(dsp.id);
        dsp.config = data.plotly_config;

        var img = document.createElement('img');
        img.className = 'img-fluid';

        get_tsne_image_data($("#sbs_tsne_gene_symbol").val(), dsp.config).then(
            data => {
                if (typeof data === 'object' || typeof data === "undefined") {
                    // If this is true, there was an error
                    $("#sbs_tsne_gene_not_found").show();
                } else {
                    // place image
                    img.src = "data:image/png;base64," + data
                    document.getElementById('sbs_tsne_plot_c').appendChild(img);
                }
            }
        );

        $("#btn_sbs_tsne_run").prop('disabled', false);

        done_working("tSNE visualized", false);
    });

    $(".show_analysis_renamer").click(function(e) {
        $("#new_analysis_label").val(current_analysis.label);
        $("#new_analysis_label_c").show(500);
    });

    // Handling all buttons manually
    $("button").click(function(e) {
        e.preventDefault();
    });


    $('#marker_genes_manually_entered').focus(function() {
        reset_manual_marker_gene_entries();
    });
    $('#marker_genes_manually_entered').keyup(function() {
        update_manual_marker_gene_entries($(this).val());
    });
    $('#marker_genes_manually_entered').change(function() {
        process_manual_marker_gene_entries($(this).val());
    });


    $('#new_analysis_label').focus(function() {
        current_label = $(this).val();
    });
    $('#new_analysis_label').keyup(function() {
        if ( analysis_labels.has($(this).val()) ) {
            // it's also OK if the current value was what it was when the user started to edit
            if ( $(this).val() != current_label ) {
                // turn red, disable save, and show dup message
                $(this).addClass('duplicate');
                $('#new_analysis_label_save').prop("disabled", true);
                $('#duplicate_label_warning').show(500);
            }

        } else {
            // clear duplication message, remove red, and enable save
            $('#duplicate_label_warning').hide(500);
            $(this).removeClass('duplicate');
            $('#new_analysis_label_save').prop("disabled", false);
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

function get_tsne_image_data(gene_symbol, config) {
    // then craziness: https://stackoverflow.com/a/48980526
    return axios.post(`/api/plot/${current_analysis.dataset_id}/tsne`, {
        gene_symbol,
        analysis: current_analysis,
        colorize_legend_by: config.colorize_legend_by,
        plot_type: 'tsne_static',
        plot_by_group: config.plot_by_group,
        max_columns: config.max_columns,
        horizontal_legend: config.horizontal_legend,
        x_axis: config.x_axis,
        y_axis: config.y_axis,
        analysis_owner_id: current_analysis.user_id,
        colors: config.colors,
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
    original_analysis_type = current_analysis.type
    current_analysis.type = 'primary'

    if ($('#filter_cells_lt_n_genes_selected').is(":checked")) {
        current_analysis.primary_filter.filter_cells_lt_n_genes_selected = true
        current_analysis.primary_filter.filter_cells_lt_n_genes = $("#filter_cells_lt_n_genes").val();
    } else {
        current_analysis.primary_filter.filter_cells_lt_n_genes_selected = false
    }


    if ($('#filter_cells_gt_n_genes_selected').is(":checked")) {
        current_analysis.primary_filter.filter_cells_gt_n_genes_selected = true
        current_analysis.primary_filter.filter_cells_gt_n_genes = $("#filter_cells_gt_n_genes").val();
    } else {
        current_analysis.primary_filter.filter_cells_gt_n_genes_selected = false
    }

    if ($('#filter_genes_lt_n_cells_selected').is(":checked")) {
        current_analysis.primary_filter.filter_genes_lt_n_cells_selected = true
        current_analysis.primary_filter.filter_genes_lt_n_cells = $("#filter_genes_lt_n_cells").val();
    } else {
        current_analysis.primary_filter.filter_genes_lt_n_cells_selected = false
    }

    if ($('#filter_genes_gt_n_cells_selected').is(":checked")) {
        current_analysis.primary_filter.filter_genes_gt_n_cells_selected = true
        current_analysis.primary_filter.filter_genes_gt_n_cells = $("#filter_genes_gt_n_cells").val();
    } else {
        current_analysis.primary_filter.filter_genes_gt_n_cells_selected = false
    }

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_apply_primary_filter.cgi",
        data: {'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type,
               'dataset_id': current_analysis.dataset_id,
               'session_id': current_analysis.user_session_id,
               'using_primary_datasource': current_analysis.using_primary_datasource,
               'filter_cells_lt_n_genes': current_analysis.primary_filter.filter_cells_lt_n_genes,
               'filter_cells_gt_n_genes': current_analysis.primary_filter.filter_cells_gt_n_genes,
               'filter_genes_lt_n_cells': current_analysis.primary_filter.filter_genes_lt_n_cells,
               'filter_genes_gt_n_cells': current_analysis.primary_filter.filter_genes_gt_n_cells
              },
        dataType: "json",
        success: function(data) {
            current_analysis.primary_filter.filtered_gene_count = data['n_genes'];
            current_analysis.primary_filter.filtered_cell_count = data['n_obs'];

            if (data['success'] == 1) {
                $("#selected_dataset_shape_filtered").html(data['n_genes'] + " genes x " + data['n_obs'] + " obs");
                $("#selected_dataset_shape_filtered_c").show(500);

                if (original_analysis_type == 'primary') {
                    current_analysis.type = 'user_unsaved'
                } else {
                    current_analysis.type = original_analysis_type
                }

                current_analysis.primary_filter.calculated = true;
                current_analysis.primary_filter.update_ui(current_analysis);
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
    if (current_analysis.type == 'public') {
        current_analysis.copy_to_user_unsaved(callback, opts);
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
        if (current_analysis.type != 'primary') {
            if ( $('select#analysis_id').find(':selected').data('analysis-id') == "0") {
                current_analysis.label = $('#new_analysis_label').val();
            }
            current_analysis.save();

            if (current_analysis.type == 'user_unsaved') {
                $("#analysis_action_c").show();
                $("#analysis_status_info_c").hide();
            } else if (current_analysis.type == 'user_saved') {
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

    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(file_contents));
    element.setAttribute('download', 'marker_genes.xls');
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

//TODO: move this into a generic utils module
function get_dataset_info(dataset_id) {
    $("#stored_analyses_c").hide();

    $.ajax({
        type: "POST",
        url: "./cgi/get_dataset_info.cgi",
        data: {'dataset_id': dataset_id, 'include_shape': 1},
        dataType: "json",
        success: function(data) {
            var ds = new Dataset(data);

            $('#new_analysis_label').val(current_analysis.label);

            // TODO: make Analysis.dataset the replacement for Analysis.dataset_id
            current_analysis.dataset = ds;

            update_selected_dataset(ds);
            $("#dataset_info").show();
            $("#analysis_list_c").show();
            analysis_labels = current_analysis.get_saved_analyses_list(ds.id, 0);
            done_working();
        },
        error: function(xhr, status, msg) {
            report_error("Failed to access dataset");
            report_error("Failed ID was: " + dataset_id + " because msg: " + msg);
        }
    });
}

function load_preliminary_figures(dataset_id) {
    $("#stored_analyses_c").hide();

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_preview_primary_filter.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id
              },
        dataType: "json",
        success: function(data) {
            $("#primary_initial_plot_loading_c").hide();

            if (data['success'] == 1) {
                $('#primary_initial_violin_c').html('<a target="_blank" href="./datasets/' + dataset_id + '.prelim_violin.png"><img src="./datasets/' + dataset_id + '.prelim_violin.png" class="img-fluid img-zoomed" /></a>');
                $('#primary_initial_scatter_c').html('<a target="_blank" href="./datasets/' + dataset_id + '.prelim_n_genes.png"><img src="./datasets/' + dataset_id + '.prelim_n_genes.png" class="img-fluid img-zoomed" /></a>');
                done_working("Prelim plots displayed");
            } else {
                $('#primary_initial_violin_c').html('Preliminary plots not yet generated. Continue your analysis.');
                done_working("Prelim figures missing.");
            }

            $("#primary_initial_plot_c").show(500);
        },
        error: function(xhr, status, msg) {
            report_error("Failed to access dataset");
            report_error("Failed ID was: " + dataset_id + " because msg: " + msg);
        }
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
            ana = Analysis.load_from_json(data);
            ana.dataset = current_analysis.dataset;
            current_analysis = ana;

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

async function populate_dataset_selection() {
    $('#pre_dataset_spinner').show();
    await $.ajax({
        type: "POST",
        url: "./cgi/get_h5ad_dataset_list.cgi",
        data: {
            'session_id': CURRENT_USER.session_id,
            'for_page': 'analyze_dataset',
            'include_dataset_id': getUrlParameter('dataset_id')
        },
        dataType: "json",
        success(data) {
            let counter = 0
            // Populate select box with dataset information owned by the user
            const user_datasets = [];
            if (data.user.datasets.length > 0) {
              // User has some profiles
              $.each(data.user.datasets, (_i, item) => {
                if (item) {
                    user_datasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
              });
            }
            // Next, add datasets shared with the user
            const shared_datasets = [];
            if (data.shared_with_user.datasets.length > 0) {
              // User has some profiles
              $.each(data.shared_with_user.datasets, (_i, item) => {
                if (item) {
                    shared_datasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id  });
                }
              });
            }
            // Now, add public datasets
            const domain_datasets = [];
            if (data.public.datasets.length > 0) {
              // User has some profiles
              $.each(data.public.datasets, (_i, item) => {
                  if (item) {
                    domain_datasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id  });
                  }
              });
            }

            dataset_tree.userDatasets = user_datasets;
            dataset_tree.sharedDatasets = shared_datasets;
            dataset_tree.domainDatasets = domain_datasets;
            dataset_tree.generateTree();

            // was there a requested dataset ID already?
            const dataset_id = getUrlParameter('dataset_id');
            if (dataset_id !== undefined) {
                $("#dataset_id").val(dataset_id);
                try {
                    $('#dataset_id').text(dataset_tree.treeData.find(e => e.dataset_id === dataset_id).text);
                    $("#dataset_id").trigger("change");
                } catch {
                    console.error(`Dataset id ${dataset_id} was not returned as a public/private/shared dataset`);
                }
            }

        },
        error(xhr, status, msg) {
            report_error(`Failed to load dataset list because msg: ${msg}`);
        }
    });
    $('#pre_dataset_spinner').hide();
}

function process_manual_marker_gene_entries(gene_str) {
    current_analysis.genes_of_interest = new Set([...entered_marker_genes, ...clicked_marker_genes]);
}

function report_error(msg) {
    $("#action_log li.working").remove();
    if (msg) {
        $("<li class='failure'>" + msg + "</li>").prependTo("#action_log");
    }
}

function reset_manual_marker_gene_entries() {
    clicked_marker_genes = new Set();

    // remember which GOI are from the table
    $.each(
        $('#marker_genes_table td.highlighted'),
        function() {
            clicked_marker_genes.add($(this).text());
        }
    );
}

function reset_workbench() {
    // Performs all the steps needed to reset the workbench if the input dataset changes

    // handle any data structures
    current_analysis.reset();

    // Reset the analysis tool fields and some display labels
    $("form.reset_on_change").trigger("reset");
    $("p.reset_on_change").hide();
    $("span.empty_on_change").empty();
    $("div.empty_on_change").empty();
    $("tbody.empty_on_change").empty();
    $("input#new_analysis_label").val('');

    // Hide any non-analysis-flow steps
    $("#analysis_sbs_tsne").hide();

    // Toggle the tool buttons to hide them in the UI
    $('.tooltoggle').bootstrapToggle('off');

    // Disable those analysis steps which have previous requirements
    $('.tooltoggle').bootstrapToggle('disable');
}

function run_analysis_compare_genes() {
    show_working("Computing comparison");
    $("#compare_genes_results_c .empty_on_change").empty();

    var compute_gene_comparison = 1;

    if (current_analysis.gene_comparison.calculated == 0) {
        compute_gene_comparison = 1;
    }

    $.ajax({
        type: "POST",
        // TODO: Drop the DPI here. These images are huge.
        url: "./cgi/h5ad_compare_genes.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'n_genes': $("#compare_genes_n_genes").val(),
               'compute_gene_comparison': compute_gene_comparison,
               'group_labels': JSON.stringify(current_analysis.group_labels),
               'query_cluster': $("#query_cluster").val(),
               'reference_cluster': $("#reference_cluster").val(),
               'method': $("#compare_genes_method").val(),
               'corr_method': $("#compare_genes_corr_method").val()
              },
        dataType: "json",
        success: function(data) {
            if (data['success'] == 1) {
                $("#analysis_compare_genes .tool_instructions").hide();
                $('#compare_genes_results_c').show(500);
                current_analysis.gene_comparison.calculated = true;
                current_analysis.gene_comparison.n_genes = $("#compare_genes_n_genes").val();
                current_analysis.gene_comparison.query_cluster = $("#query_cluster").val();
                current_analysis.gene_comparison.reference_cluster = $("#reference_cluster").val();
                current_analysis.gene_comparison.method = $("#compare_genes_method").val();
                current_analysis.gene_comparison.corr_method = $("#compare_genes_corr_method").val();
                current_analysis.gene_comparison.update_ui_with_results(current_analysis, data);
                done_working("Comparison computed");
            } else {
                done_working("Comparison failed.");
            }
        },
        error: function(xhr, status, msg) {
            report_error("Error running comparison");
        }
    });
}

function run_analysis_louvain() {
    show_working("Computing Louvain clusters");
    $("#analysis_louvain div.image_result_c").empty();

    var compute_louvain = true;

    current_analysis.group_labels = [];
    $("input[name='group_labels[]']").each(function() {
        current_analysis.group_labels.push($(this).val());
    });
    // update gene comparison options to include new labels
    current_analysis.gene_comparison.populate_group_selectors(current_analysis.group_labels);

    // TODO: check parameters to be sure we don't need to recluster
    if (current_analysis.louvain.calculated == true) {
        if (current_analysis.louvain.n_neighbors == $("#louvain_n_neighbors").val() &&
            current_analysis.louvain.resolution == $("#louvain_resolution").val()) {
            compute_louvain = false;
        }
    }

    var plot_tsne = 0;
    var plot_umap = 0;

    if ($("#dimensionality_reduction_method_tsne").is(":checked")) {
        plot_tsne = 1;
    }

    if ($("#dimensionality_reduction_method_umap").is(":checked")) {
        plot_umap = 1;
    }

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_generate_louvain.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'resolution': $("#louvain_resolution").val(),
               'compute_louvain': compute_louvain,
               'plot_tsne': plot_tsne,
               'plot_umap': plot_umap,
               'group_labels': JSON.stringify(current_analysis.group_labels)
              },
        dataType: "json",
        success: function(data) {
            if (data['success'] == 1) {
                current_analysis.louvain.calculated = true;
                current_analysis.louvain.n_neighbors = $("#louvain_n_neighbors").val();
                current_analysis.louvain.resolution = $("#louvain_resolution").val();
                current_analysis.louvain.plot_tsne = plot_tsne;
                current_analysis.louvain.plot_umap = plot_umap;
                current_analysis.louvain.update_ui(current_analysis);

                $('#btn_louvain_run').attr("disabled", false);
                done_working("Louvain clusters computed");
                $('#louvain_run_c .js-next-step').show();  // Show that next toggle can be clicked
            } else {
                $('#btn_louvain_run').attr("disabled", false);
                done_working("Louvain cluster compute failed.");
            }
        },
        error: function(xhr, status, msg) {
            report_error("Error generating clusters");
            $('#btn_louvain_run').attr("disabled", false);
        }
    });
}

function run_analysis_marker_genes() {
    show_working("Computing marker genes");
    $('#btn_marker_genes_run').attr("disabled", true);
    $("#marker_genes_plot_c").empty();
    $("#marker_genes_table_c .empty_on_change").empty();

    $('#marker_genes_manually_entered').val('');
    $('#marker_genes_selected_count').text(0);
    $('#marker_genes_entered_count').text(0);
    $('#marker_genes_unique_count').text(0);
    current_analysis.genes_of_interest = new Set();
    clicked_marker_genes = new Set();
    entered_marker_genes = new Set();

    let compute_marker_genes = true;
    // If marker genes were already calculated, just reuse the results... unless they change the number they want.
    if (current_analysis.marker_genes.calculated == true
        && $("#marker_genes_n_genes").val() == current_analysis.marker_genes.n_genes) {
        compute_marker_genes = false;
    }

    $.ajax({
        type: "POST",
        // TODO: Drop the DPI here. These images are huge.
        url: "./cgi/h5ad_find_marker_genes.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'n_genes': $("#marker_genes_n_genes").val(),
               'compute_marker_genes': compute_marker_genes
              },
        dataType: "json",
        success: function(data) {
            if (data['success'] == 1) {
                $("#analysis_marker_genes .tool_instructions").hide(500);
                current_analysis.marker_genes.calculated = true;
                current_analysis.marker_genes.n_genes = $("#marker_genes_n_genes").val();
                current_analysis.marker_genes.update_ui_with_results(current_analysis, data);
                current_analysis.group_labels = data['group_labels'].map(x => x.group_label);

                $('#btn_marker_genes_run').attr("disabled", false);
                done_working("Marker genes computed");
                $('#marker_genes_options_c .js-next-step').show();  // Show that next toggle can be clicked
            } else {
                $('#btn_marker_genes_run').attr("disabled", false);
                done_working("Marker genes compute failed.");
            }
        },
        error: function(xhr, status, msg) {
            report_error("Error reporting marker genes");
            $('#btn_marker_genes_run').attr("disabled", false);
        }
    });
}

async function run_analysis_pca_top_genes() {
    show_working('Computing top genes for principal components');
    $("#pca_top_genes_c").empty();

    const response = await $.ajax({
        type: "POST",
        url: "./api/analysis/plotTopGenesPCA",
        data: {
            'dataset_id': current_analysis.dataset_id,
            'analysis_id': current_analysis.id,
            'analysis_type': current_analysis.type,
            'session_id': current_analysis.user_session_id,
            'pcs': $("#top_pca_genes").val(),
        },
        dataType: "json"
    }).catch(() => report_error("Error finding top genes"));

    if (response.success == 1) {
        const params = {
            'analysis_id': current_analysis.id,
            'analysis_name': 'pca_loadings',
            'analysis_type': current_analysis.type,
            'dataset_id': current_analysis.dataset_id,
            'session_id': current_analysis.user_session_id,
            datetime: (new Date()).getTime()
        };
        $("#pca_top_genes_c").show(500);
        current_analysis.place_analysis_image({
            'params': params,
            'title': 'Top PCA Genes',
            'target': '#pca_top_genes_c'
        });
        $("#pca_top_genes_c").show(500);
        done_working("Top PCA genes visualized", false);
    } else if (response.success == 0 ) {
        report_error(response.message);
    } else {
        report_error("Error finding top genes");
    }

    $('#btn_pca_top_genes').attr("disabled", false);
}

function run_analysis_pca() {
    show_working("Computing principal components and variance plot");
    $("#pca_missing_gene_c").hide(500);
    $("#analysis_pca div.image_result_c").empty();

    var compute_pca = true;

    if (current_analysis.pca.calculated == true) {
        compute_pca = false;
    }

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_generate_pca.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'genes_to_color': $("#pca_genes_to_color").val(),
               'compute_pca': compute_pca
              },
        dataType: "json",
        success: function(data) {
            if (data['success'] == 1) {
                current_analysis.pca.calculated = true;
                current_analysis.pca.genes_to_color = $("#pca_genes_to_color").val();
                current_analysis.pca.update_ui_with_results(current_analysis, data);
                done_working("PCA and variance computed");
                $("#pca_options_g").show();
                $('#pca_options_c .js-next-step').show();  // Show that next toggle can be clicked
            } else {
                $("#pca_missing_gene").text(data['missing_gene']);
                $("#pca_missing_gene_c").show(500);
                done_working("PCA: requested gene missing");
            }

            $('#btn_pca_run').attr("disabled", false);
        },
        error: function(xhr, status, msg) {
            report_error("Error running PCA");
            $('#btn_pca_run').attr("disabled", false);
        }
    });
}

function run_analysis_qc_by_mito(save_dataset) {
    if (save_dataset == 0) {
        $("#btn_qbm_save").prop('disabled', true);
        $("#btn_qbm_save").show();
        show_working("Analyzing mitochondrial genes");
    } else {
        show_working("Applying mitochondrial filter");
        $("#btn_qbm_save").text('Saving');
        $("#btn_qbm_save").prop('disabled', true);
    }

    // clear the images
    $("#qbm_violin_c").empty();
    $("#qbm_scatter_percent_mito_c").empty();
    $("#qbm_scatter_n_genes_c").empty();

    current_analysis.qc_by_mito.gene_prefix = $("#qbm_gene_prefix").val();
    current_analysis.qc_by_mito.filter_mito_perc = $('#qbm_filter_mito_perc').val();
    current_analysis.qc_by_mito.filter_mito_count = $('#qbm_filter_mito_count').val();

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_qc_by_mito.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'genes_prefix': current_analysis.qc_by_mito.gene_prefix,
               'filter_mito_perc': current_analysis.qc_by_mito.filter_mito_perc,
               'filter_mito_count': current_analysis.qc_by_mito.filter_mito_count,
               'save_dataset': save_dataset
              },
        dataType: "json",
        success: function(data) {
            if (data['success'] == 1) {
                current_analysis.qc_by_mito.n_genes = data['n_genes'];
                current_analysis.qc_by_mito.n_obs = data['n_obs'];

                if (save_dataset == 1) {
                    $("#btn_qbm_save").prop('disabled', true);
                    current_analysis.qc_by_mito.calculated = true;
                } else {
                    $("#btn_qbm_save").prop('disabled', false);
                    current_analysis.qc_by_mito.calculated = false;
                }

                current_analysis.qc_by_mito.update_ui_with_results(
                    current_analysis, data, current_analysis.qc_by_mito.calculated);
                done_working("Mitochondrial plot displayed");
                $('#qc_mito_options_c .js-next-step').show();  // Show that next toggle can be clicked
            } else {
                done_working("QC by mito failed.");
            }

            $('#btn_do_analysis_qc_by_mito').attr("disabled", false);
        },
        error: function(xhr, status, msg) {
            report_error("Error doing QC analysis");
            $('#btn_do_analysis_qc_by_mito').attr("disabled", false);
        }
    });
}

function run_analysis_tsne() {
    show_working("Computing tSNE/UMAP and generating plot");
    $("#tsne_missing_gene_c").hide(500);
    $("#tsne_plot_c").empty();
    $("#umap_plot_c").empty();

    // Anytime we run this there are three things which might need to be computed depending on what
    //  has happened.  neighborhood, tSNE and UMAP need to be done if it's the first time or if
    //  the settings have changed.
    // Also, both tSNE and UMAP can be replotted (with different gene coloring) and not recomputed

    var compute_neighbors = 1;
    var compute_tsne = 0;
    var compute_umap = 0;
    var plot_tsne = 0;
    var plot_umap = 0;

    if ($("#dimensionality_reduction_method_tsne").is(":checked")) {
        compute_tsne = 1;
        plot_tsne = 1;
    }

    if ($("#dimensionality_reduction_method_umap").is(":checked")) {
        compute_umap = 1;
        plot_umap = 1;
    }

    // We don't have to recompute if none of the plotting parameters have changed.  This allows
    //  us to just do something like recolor.
    if (current_analysis.tsne.neighbors_calculated == 1) {
        if (current_analysis.tsne.n_pcs == $("#tsne_n_pcs").val() &&
            current_analysis.tsne.n_neighbors == $("#dredux_n_neighbors").val() &&
            current_analysis.tsne.random_state == $("#tsne_random_state").val()) {
            compute_neighbors = 0;
        }
    }

    // don't recompute tSNE if we already have and nothing has changed
    if (current_analysis.tsne.tsne_calculated == 1) {
        if (current_analysis.tsne.n_pcs == $("#tsne_n_pcs").val() &&
            current_analysis.tsne.random_state == $("#tsne_random_state").val()) {
            compute_tsne = 0;
        }
    }

    // Don't recompute UMAP if we already have
    if (current_analysis.tsne.umap_calculated == 1) {
        compute_umap = 0;
    }

    // If plotting of either is requested and neighbors are being recomputed, we need to
    //  also recompute tSNE/UMAP
    if (compute_neighbors == 1) {
        if (plot_umap == 1) {
            compute_umap = 1;
        }

        if (plot_tsne == 1) {
            compute_tsne = 1;
        }
    }

    var use_scaled = false;
    if ($("#tsne_use_scaled").is(":checked")) {
        use_scaled = true;
    }

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_generate_tsne.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'genes_to_color': $("#tsne_genes_to_color").val(),
               'n_pcs': $("#tsne_n_pcs").val(),
               'n_neighbors': $("#dredux_n_neighbors").val(),
               'random_state': $("#tsne_random_state").val(),
               'use_scaled': use_scaled,
               'compute_neighbors': compute_neighbors,
               'compute_tsne': compute_tsne,
               'compute_umap': compute_umap,
               'plot_tsne': plot_tsne,
               'plot_umap': plot_umap
              },
        dataType: "json",
        success: function(data) {
            if (data['success'] == 1) {
                $("#analysis_tsne .tool_instructions").hide(500);
                current_analysis.tsne.neighbors_calculated = true;
                current_analysis.tsne.n_pcs = $("#tsne_n_pcs").val();
                current_analysis.tsne.n_neighbors = $("#dredux_n_neighbors").val();
                current_analysis.tsne.random_state = $("#tsne_random_state").val();
                current_analysis.tsne.genes_to_color = $("#tsne_genes_to_color").val();
                current_analysis.tsne.plot_tsne = plot_tsne;
                current_analysis.tsne.plot_umap = plot_umap;
                current_analysis.tsne.update_ui(current_analysis);

                if (compute_tsne == 1) {
                    current_analysis.tsne.tsne_calculated = 1;
                }

                if (compute_umap == 1) {
                    current_analysis.tsne.umap_calculated = 1;
                }

                $('#btn_tsne_run').attr("disabled", false);
                done_working("tSNE/UMAP computed and displayed");
                $('#tsne_options_c .js-next-step').show();  // Show that next toggle can be clicked
            } else {
                $("#tsne_missing_gene").text(data['missing_gene']);
                $("#tsne_missing_gene_c").show(500);
                $('#btn_tsne_run').attr("disabled", false);
                done_working("tSNE generation failed.");
            }
        },
        error: function(xhr, status, msg) {
            report_error("Error generating tSNE");
            $('#btn_tsne_run').attr("disabled", false);
        }
    });
}

function run_analysis_select_variable_genes(save_dataset) {
    if (save_dataset == 0) {
        show_working("Analyzing variable genes");
    } else {
        show_working("Saving variable genes");
        current_analysis.select_variable_genes.calculated = true;
    }

    $("#asvg_plot_c").empty();
    $("#asvg_plot_norm_c").empty();
    $("#btn_asvg_save").hide();

    current_analysis.select_variable_genes.norm_counts_per_cell = $("#asvg_norm_counts_per_cell").val();
    current_analysis.select_variable_genes.flavor = $("#asvg_flavor").val();
    current_analysis.select_variable_genes.n_top_genes = $("#asvg_n_top_genes").val();
    current_analysis.select_variable_genes.min_mean = $("#asvg_min_mean").val();
    current_analysis.select_variable_genes.max_mean = $("#asvg_max_mean").val();
    current_analysis.select_variable_genes.min_dispersion = $("#asvg_min_dispersion").val();

    if ($('#asvg_regress_out').is(":checked")) {
        current_analysis.select_variable_genes.regress_out = true;
    } else {
        current_analysis.select_variable_genes.regress_out = false;
    }

    if ($('#asvg_scale_unit_variance').is(":checked")) {
        current_analysis.select_variable_genes.scale_unit_variance = true;
    } else {
        current_analysis.select_variable_genes.scale_unit_variance = false;
    }

    $.ajax({
        type: "POST",
        url: "./cgi/h5ad_identify_variable_genes.cgi",
        data: {'dataset_id': current_analysis.dataset_id, 'analysis_id': current_analysis.id,
               'analysis_type': current_analysis.type, 'session_id': current_analysis.user_session_id,
               'norm_counts_per_cell': current_analysis.select_variable_genes.norm_counts_per_cell,
               'flavor': current_analysis.select_variable_genes.flavor,
               'n_top_genes': current_analysis.select_variable_genes.n_top_genes,
               'min_mean': current_analysis.select_variable_genes.min_mean,
               'max_mean': current_analysis.select_variable_genes.max_mean,
               'min_dispersion': current_analysis.select_variable_genes.min_dispersion,
               'regress_out': current_analysis.select_variable_genes.regress_out,
               'scale_unit_variance': current_analysis.select_variable_genes.scale_unit_variance,
               'save_dataset': save_dataset
              },
        dataType: "json",
        success: function(data) {
            if (save_dataset == 1) {
                current_analysis.select_variable_genes.calculated = true;
            } else {
                current_analysis.select_variable_genes.calculated = false;
            }

            current_analysis.select_variable_genes.update_ui_with_results(
                current_analysis, data, current_analysis.select_variable_genes.calculated);
            done_working("Variable genes image created");

            $("#asvg_result_count").text("(" + data['n_genes'] + ")");
            $(".asvg_save_options").show();
            $("#btn_asvg_save").show();
            $('#btn_do_analysis_select_variable_genes').attr("disabled", false);
            $('#top_genes').html("Suggested highly-variable genes:<br><strong>" + data['top_genes'] + "<strong>");
            $('#top_genes').show();
            $('#asvg_options_c .js-next-step').show();  // Show that next toggle can be clicked
        },
        error: function(xhr, status, msg) {
            report_error("Error identifying variable genes");
            $('#btn_do_analysis_select_variable_genes').attr("disabled", false);
        }
    });
}

function show_working(msg) {
    $("<li class='working'><img src='./img/loading_search.gif' alt='Working' /><span>" +
      msg + "</span></li>").prependTo("#action_log");
}

function update_manual_marker_gene_entries(gene_str) {
    entered_marker_genes = new Set();
    var gene_syms = gene_str.split(',');

    for (gene_sym of gene_syms) {
        gene_sym = gene_sym.trim();
        if (gene_sym) {
            entered_marker_genes.add(gene_sym);
        }
    }

    var counter_set = new Set([...entered_marker_genes, ...clicked_marker_genes]);
    $('#marker_genes_unique_count').text(counter_set.size);
    $('#marker_genes_entered_count').text(entered_marker_genes.size);
}

function update_selected_dataset(ds) {
    $("#selected_dataset_title").html(ds['title']);
    $("#selected_dataset_shape_initial").html(ds.shape());
    $("<li>Changed dataset to: " + ds.title + "</li>").prependTo("#action_log");
}
