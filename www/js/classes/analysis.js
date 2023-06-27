"use strict";

/*
  Classes representing overall analysis (pipeline) elements and their related/child
  classes.

  TODO: Separate this into analysis.js and ui-panel-analysis.js
*/

// requires common.js

class Analysis {
    constructor ({id, dataset_id, dataset_is_raw, label, type, vetting, user_session_id, genes_of_interest, group_labels} = {}) {
        this.id = id == null ? uuid() : id;
        this.user_session_id = user_session_id == null ? CURRENT_USER.session_id : user_session_id;

        this.dataset_id = dataset_id;
        // This next one isn't fully implemented yet, but should be.
        // Currently, the only place it gets populated is analyze_dataset.js:get_dataset_info()
        // this.dataset = dataset

        // should be either 'primary', 'public', 'user_saved', 'user_unsaved'
        this.type = type;

        // should be either null, 'owner', 'gear' or 'community'
        this.vetting = vetting;

        // free-form label to display for this analysis
        var new_label = 'Unlabeled ' + common_datetime();
        this.label = label == null ? new_label : label;

        // is this raw or log-transformed?
        this.dataset_is_raw = dataset_is_raw == null ? true : dataset_is_raw;

        // By keeping track of these, we can prevent re-running them for things
        this.louvain = new AnalysisStepLouvain();
        this.marker_genes = new AnalysisStepMarkerGenes();
        this.pca = new AnalysisStepPCA();
        this.tsne = new AnalysisSteptSNE();
        this.primary_filter = new AnalysisStepPrimaryFilter();
        this.qc_by_mito = new AnalysisStepQCByMito();
        this.select_variable_genes = new AnalysisStepSelectVariableGenes();
        this.gene_comparison = new AnalysisStepCompareGenes();

        this.group_labels = group_labels == null ? [] : group_labels;
        this.genes_of_interest = genes_of_interest == null ? new Set() : genes_of_interest;

        if (this.genes_of_interest instanceof Array) {
            this.genes_of_interest = new Set(this.genes_of_interest);

        // checks for an empty object like {}
        } else if (Object.entries(this.genes_of_interest).length === 0 && this.genes_of_interest.constructor === Object) {
            this.genes_of_interest = new Set();
        }
    }

    add_gene_of_interest(gene_symbol) {
        gene_symbol = gene_symbol.trim();
        this.genes_of_interest.add(gene_symbol);

        if (this.genes_of_interest.length == 1) {
            $("#btn_visualize_marker_genes").show(300);
        }
    }

    copy_to_user_unsaved(callback, opts) {
        /* This function takes an analysis and makes a copy of it to the current user's
           unsaved area.  The most common use case is to copy the public analysis of another
           user so changes can be made.
        */
        var this_analysis = this;
        var new_analysis_id = uuid();

        $.ajax({
            type: "POST",
            url: "./cgi/copy_dataset_analysis.cgi",
            data: {session_id: CURRENT_USER.session_id,
                   dataset_id: this.dataset_id,
                   source_analysis_id: this.id,
                   dest_analysis_id: new_analysis_id,
                   dest_analysis_type: 'user_unsaved',
                   source_analysis_type:this.type},
            dataType: "json",
            success: function(data) {
                this_analysis.type = 'user_unsaved';
                this_analysis.id = new_analysis_id;
                this_analysis.user_session_id = CURRENT_USER.session_id;

                $("#analysis_action_c").show();
                $("#analysis_status_info_c").hide();

                this_analysis.get_saved_analyses_list(this_analysis.dataset_id, new_analysis_id);

                if (callback) {
                    callback(opts);
                }
            },
            error: function(xhr, status, msg) {
                report_error("Error creating sandbox copy of the analysis: " + msg);
            }
        });
    }

    delete() {
        var this_analysis = this;

        $.ajax({
            type: "POST",
            url: "./cgi/delete_dataset_analysis.cgi",
            data: {session_id: CURRENT_USER.session_id,
                   dataset_id: this.dataset_id,
                   analysis_id: this.id,
                   analysis_type: this.type},
            dataType: "json",
            success: function(data) {
                if (data['success'] == 1) {
                    // Trigger the selection of a 'New' analysis
                    $('#analysis_id option[data-analysis-id="0"]').prop("selected", true).change();
                    this_analysis.get_saved_analyses_list(this_analysis.dataset_id, 0);

                } else {
                    console.log("Deletion failed: " + data['error']);
                }
            },
            error: function(xhr, status, msg) {
                report_error("Error deleting analysis: " + msg);
            }
        });
    }

    get_saved_analyses_list(dataset_id, selected_analysis_id) {
        $.ajax({
            type: "POST",
            url: "./cgi/get_stored_analysis_list.cgi",
            data: {'dataset_id': dataset_id, 'session_id': CURRENT_USER.session_id},
            dataType: "json",
            success: function(data) {
                var empty_analysis_list_html = $('#analyses_list_empty_tmpl').html();
                var al = new Set();

                if (data['primary'].length > 0) {
                    var primary_list_tmpl = $.templates("#analyses_list_tmpl");
                    var primary_list_html = primary_list_tmpl.render(data['primary']);
                    $("#analyses_primary").html(primary_list_html);

                    for (var anum in data['primary']) {
                        al.add(data['primary'][anum]['label']);
                    }
                } else {
                    $("#analyses_primary").html(empty_analysis_list_html);
                }

                if (data['user_unsaved'].length > 0) {
                    var unsaved_list_tmpl = $.templates("#analyses_list_tmpl");
                    var unsaved_list_html = unsaved_list_tmpl.render(data['user_unsaved']);
                    $("#analyses_unsaved").html(unsaved_list_html);

                    for (var anum in data['user_unsaved']) {
                        al.add(data['user_unsaved'][anum]['label']);
                    }
                } else {
                    $("#analyses_unsaved").html(empty_analysis_list_html);
                }

                if (data['user_saved'].length > 0) {
                    var saved_list_tmpl = $.templates("#analyses_list_tmpl");
                    var saved_list_html = saved_list_tmpl.render(data['user_saved']);
                    $("#analyses_saved").html(saved_list_html);

                    for (var anum in data['user_saved']) {
                        al.add(data['user_saved'][anum]['label']);
                    }
                } else {
                    $("#analyses_saved").html(empty_analysis_list_html);
                }

                if (data['public'].length > 0) {
                    var public_list_tmpl = $.templates("#analyses_list_tmpl");
                    var public_list_html = public_list_tmpl.render(data['public']);
                    $("#analyses_public").html(public_list_html);

                    for (var anum in data['public']) {
                        al.add(data['public'][anum]['label']);
                    }
                } else {
                    $("#analyses_public").html(empty_analysis_list_html);
                }

                // preselect any analysis ID
                $('#analysis_id option[data-analysis-id="' + selected_analysis_id + '"]').attr("selected", "selected");

                $('#analysis_id').selectpicker('refresh');
                $("#stored_analyses_c").show(10);

                // global in analyze_dataset.js
                analysis_labels = al;
            },
            error: function(xhr, status, msg) {
                report_error("Failed to access dataset");
                report_error("Failed ID was: " + dataset_id + " because msg: " + msg);
            }
        });
    }

    static load_from_json(data) {
        var ana = new Analysis({
            'id': data['id'],
            'dataset_id': data['dataset_id'],
            'dataset_is_raw': data['dataset_is_raw'],
            'label': data['label'],
            'type': data['type'],
            'user_session_id': data['user_session_id'],
            'group_labels': data['group_labels'],
            'genes_of_interest': data['genes_of_interest']
        });

        if (ana.type == 'primary') {
            // If showing a primary display we only want to show marker genes and gene comparison
            //  tools
            $("#toggle_marker_genes").bootstrapToggle('enable');
            $("#toggle_marker_genes").bootstrapToggle('on');
            $("#dataset_info").hide();

        } else {
            ana.primary_filter = AnalysisStepPrimaryFilter.load_from_json(ana, data['primary_filter']);
            ana.primary_filter.update_ui(ana);

            ana.qc_by_mito = AnalysisStepQCByMito.load_from_json(ana, data['qc_by_mito']);
            ana.select_variable_genes = AnalysisStepSelectVariableGenes.load_from_json(ana, data['select_variable_genes']);
            ana.pca = AnalysisStepPCA.load_from_json(ana, data['pca']);

            ana.tsne = AnalysisSteptSNE.load_from_json(ana, data['tsne']);
            ana.tsne.update_ui(ana);

            ana.louvain = AnalysisStepLouvain.load_from_json(ana, data['louvain']);
            ana.louvain.update_ui(ana);

            ana.marker_genes = AnalysisStepMarkerGenes.load_from_json(ana, data['marker_genes']);
            ana.gene_comparison = AnalysisStepCompareGenes.load_from_json(ana, data['gene_comparison']);
        }

        return ana;
    }

    make_public_copy() {
        var this_analysis = this;

        $.ajax({
            type: "POST",
            url: "./cgi/copy_dataset_analysis.cgi",
            data: {session_id: CURRENT_USER.session_id,
                   dataset_id: this.dataset_id,
                   source_analysis_id: this.id,
                   dest_analysis_id: this.id,
                   source_analysis_type:this.type,
                   dest_analysis_type:'public'
                  },
            dataType: "json",
            success: function(data) {
                this_analysis.type = 'user_saved';
                $("#analysis_action_c").hide();
                $("#analysis_status_info").text("Changes made to this public analysis will spawn a local copy within your profile.");
                $("#analysis_status_info_c").show();
                $("#btn_make_public_copy").hide();
                $("#btn_delete_saved_analysis").hide();
                $("#btn_delete_unsaved_analysis").hide();

                this_analysis.get_saved_analyses_list(this_analysis.dataset_id, this_analysis.id);
            },
            error: function(xhr, status, msg) {
                report_error("Error saving analysis state: " + msg);
            }
        });
    }

    place_analysis_image({params, title, target = []} = {}) {
        /*
          Grabs an image from the server, usually generated temporarily by a
          module like scanpy, and places it into the target location as a binary stream.
          This prevents us from having to keep a lot of temporary images on the server
          and also ensures the user never gets a cached image.
         */
        var param_str = jQuery.param(params);
        var img_src = './cgi/get_analysis_image.cgi?' + param_str;
        $("<a target='_blank' href='" + img_src + "'><img src='" + img_src + "' class='img-fluid img-zoomed' alt='" + title + "' /></a>").appendTo(target);
    }

    remove_gene_of_interest(gene_symbol) {
        this.genes_of_interest.delete(gene_symbol);

        if (this.genes_of_interest.length == 0) {
            $("#btn_visualize_marker_genes").hide(300);
        }
    }

    reset() {
        /*
          Starts over as a completely new Analysis instance, with new ID, and resets
          existing component
         */
        this.louvain.reset();
        this.marker_genes.reset();
        this.pca.reset();
        this.tsne.reset();
        this.primary_filter.reset();
        this.qc_by_mito.reset();
        this.select_variable_genes.reset();
        this.gene_comparison.reset();
        this.dataset_is_raw = true;
        this.id = uuid();
        this.label = null;
    }

    save() {
        /*
          Saves the current analysis parameters to disk.
         */
        var state = JSON.stringify(this);
        var this_analysis = this;

        $.ajax({
            type: "POST",
            url: "./cgi/save_dataset_analysis.cgi",
            data: {session_id: CURRENT_USER.session_id,
                   dataset_id: this.dataset_id,
                   analysis_id: this.id,
                   analysis_type: this.type,
                   analysis_vetting: this.vetting,
                   label: this.label,
                   state: state},
            dataType: "json",
            success: function(data) {
                this_analysis.get_saved_analyses_list(this_analysis.dataset_id, this_analysis.id);
            },
            error: function(xhr, status, msg) {
                report_error("Error saving analysis state: " + msg);
            }
        });
    }

    save_to_user_area() {
        var this_analysis = this;

        $.ajax({
            type: "POST",
            url: "./cgi/copy_dataset_analysis.cgi",
            data: {session_id: CURRENT_USER.session_id,
                   dataset_id: this.dataset_id,
                   source_analysis_id: this.id,
                   dest_analysis_id: this.id,
                   source_analysis_type:this.type,
                   dest_analysis_type:'user_saved'
                  },
            dataType: "json",
            success: function(data) {
                this_analysis.type = 'user_saved';
                $("#btn_save_analysis").text("Saved");
                $("#analysis_action_c").hide();
                $("#analysis_status_info").text("This analysis is stored in your profile.");
                $("#analysis_status_info_c").show(500);
                $("#btn_delete_saved_analysis").show();
                $("#btn_make_public_copy").show();
                $("#new_analysis_label_c").hide();

                this_analysis.get_saved_analyses_list(this_analysis.dataset_id, this_analysis.id);
            },
            error: function(xhr, status, msg) {
                report_error("Error saving analysis state: " + msg);
            }
        });
    }
}

class AnalysisStepCompareGenes {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepCompareGenes();
        if (data === undefined) {
            data = {
                "calculated": false,
                "n_genes": 0,
                "query_cluster": null,
                "reference_cluster": null,
                "method": 't-test_overestim_var',
                "corr_method": 'benjamini-hochberg',
                "table_json_f": null,
                "table_json_r": null
            };
        }

        step.calculated = data['calculated'];
        step.n_genes = data['n_genes'];
        step.query_cluster = data['query_cluster'];
        step.reference_cluster = data['reference_cluster'];
        step.method = 't-test_overestim_var';
        step.corr_method = 'benjamini-hochberg';

        if (data.hasOwnProperty('method')) {
            step.method = data['method'];
        }

        if (data.hasOwnProperty('corr_method')) {
            step.corr_method = data['corr_method'];
        }

        if (data.hasOwnProperty('table_json')) {
            step.table_json = data['table_json'];
        }

        if (step.calculated == true) {
            step.update_ui_with_results(ana, data, true);
        }

        return step;
    }

    populate_group_selectors(group_labels) {
        var selectors = [];
        var group_num = 0;

        for (var label of group_labels) {
            //selectors.push({label:label, value:group_num });
            selectors.push({label:label, value:label });
            group_num += 1;
        }

        var list_tmpl = $.templates("#cluster_list_tmpl");
        var list_html = list_tmpl.render(selectors);
        $("#query_cluster_options").html(list_html);
        $("#reference_cluster_options").html(list_html);

        $("#query_cluster").val(this.query_cluster);
        $("#reference_cluster").val(this.reference_cluster);
    }

    populate_comparison_table(ana, table_id, table_json) {
        // add the table rows
        var body_tmpl = $.templates("#compare_genes_table_body_tmpl");
        var body_html = body_tmpl.render(table_json);
        $("#" + table_id + " tbody").html(body_html);
    }

    reset() {
        this.calculated = false;
        this.reset_ui();
        this.n_genes = 0;
        this.query_cluster = null;
        this.reference_cluster = null;
        this.method = 't-test_overestim_var';
        this.corr_method = 'benjamini-hochberg';
        this.table_json_f = null;
        this.table_json_r = null;
    }

    reset_ui() {
        $("#marker_genes_n_genes").val('');
        $("#query_cluster").val('');
        $("#reference_cluster").val('all-reference-clusters');
        $("#compare_genes_method").val('t-test_overestim_var');
        $("#compare_genes_corr_method").val('benjamini-hochberg');
        $('#analysis_compare_genes .tool_instructions').show();
        $('#compare_genes_results_c').hide();
    }

    update_ui_with_results(ana, data, results_saved) {
        $('#analysis_compare_genes .tool_instructions').hide();
        $('#compare_genes_results_c').show(500);

        if (this.calculated == true) {
            $('#toggle_compare_genes').bootstrapToggle('on');

            $("#compare_genes_n_genes").val(this.n_genes);
            $("#compare_genes_method").val(this.method);
            $("#compare_genes_corr_method").val(this.corr_method);
            // the query and reference cluster values have to be loaded after the option lists

            var params = {
                'analysis_id': ana.id,
                'analysis_name': 'rank_genes_groups_' + data['cluster_label'] + '_comp_ranked',
                'analysis_type': ana.type,
                'dataset_id': ana.dataset_id,
                'session_id': ana.user_session_id,
                // this saves the user from getting a cached image each time
                datetime: (new Date()).getTime()
            }

            ana.place_analysis_image(
                {'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_ranked_c'});

            params['analysis_name'] = 'rank_genes_groups_' + data['cluster_label'] + '_' + this.query_cluster + '_comp_violin'
            ana.place_analysis_image(
                {'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_violin_c'});

            if (data.hasOwnProperty('table_json_f')) {
                this.table_json_f = data['table_json_f'];
                data['table_json_f'] = JSON.parse(data['table_json_f']);
                this.populate_comparison_table(ana, 'compare_genes_table_f', data['table_json_f']['data']);
            }

            if (data.hasOwnProperty('table_json_r')) {
                this.table_json_r = data['table_json_r'];
                data['table_json_r'] = JSON.parse(data['table_json_r']);
                this.populate_comparison_table(ana, 'compare_genes_table_r', data['table_json_r']['data']);
            }

            if (this.reference_cluster != 'all-reference-clusters') {
                params['analysis_name'] ='rank_genes_groups_' + data['cluster_label'] + '_comp_ranked_rev'
                ana.place_analysis_image(
                    {'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_ranked_rev_c'});

                params['analysis_name'] = 'rank_genes_groups_' + data['cluster_label'] + '_' + this.reference_cluster + '_comp_violin_rev'
                ana.place_analysis_image(
                    {'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_violin_rev_c'});
            }
        }
    }
}

class AnalysisStepLouvain {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepLouvain();
        if (data === undefined) return step;

        step.calculated = data['calculated'];
        step.n_neighbors = data['n_neighbors'];
        step.resolution = data['resolution'];
        step.plot_tsne = data['plot_tsne'];
        step.plot_umap = data['plot_umap'];

        return step;
    }

    reset() {
        this.calculated = false;
        this.n_neighbors = null;
        this.resolution = null;
        this.plot_umap = 0;
        this.plot_tsne = 0;
        this.reset_ui();
    }

    reset_ui() {
        // this should match the value in the form element
        $("#louvain_resolution").val(1.3);
    }

    update_ui(ana) {
        $('#analysis_louvain .tool_instructions').hide(500);

        if (this.calculated == true) {
            $('#toggle_louvain').bootstrapToggle('on');
            $('#louvain_n_neighbors').val(this.n_neighbors);
            $('#louvain_resolution').val(this.resolution);

            var params = {
                'analysis_id': ana.id,
                'analysis_name': 'tsne_louvain',
                'analysis_type': ana.type,
                'dataset_id': ana.dataset_id,
                'session_id': ana.user_session_id,
                // this saves the user from getting a cached image each time
                datetime: (new Date()).getTime()
            }

            if (this.plot_tsne == 1) {
                ana.place_analysis_image(
                    {'params': params, 'title': 'Louvain groups', 'target': '#louvain_tsne_plot_c'});
            }

            if (this.plot_umap == 1) {
                params['analysis_name'] = 'umap_louvain'
                ana.place_analysis_image(
                    {'params': params, 'title': 'Louvain groups', 'target': '#louvain_umap_plot_c'});
            }

            // clustering enables marker gene identification
            $("#toggle_marker_genes").bootstrapToggle('enable');
        }
    }
}

class AnalysisStepMarkerGenes {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepMarkerGenes();
        if (data === undefined) step;

        step.calculated = data['calculated']
        step.n_genes = data['n_genes']
        step.group_names = data['group_names']

        if (step.calculated == true) {
            $('#toggle_marker_genes').bootstrapToggle('on');
            $('#marker_genes_n_genes').val(step.n_genes);
            step.update_ui_with_results(ana, data);
        }

        return step;
    }

    count_and_highlight_duplicates() {
        var all_values = [];
        var dup_values = [];

        // first remove any duplicate-labeled ones
        $('#marker_genes_group_labels td.group_user_label input').removeClass('duplicate');

        $('#marker_genes_group_labels td.group_user_label input').each(function() {
            var cluster_label = $(this).val();

            // this means it WAS found
            if ( $.inArray( cluster_label, all_values ) > -1 ) {
                dup_values.push(cluster_label);
                $(this).addClass('duplicate');
            } else {
                all_values.push(cluster_label);
            }
        });

        return dup_values.length;
    }

    populate_marker_genes_labels(ana, data) {
        this.group_labels = [];

        // If the user has saved labels before, put them in the table here.  Else leave it
        //  as the top gene.
        let i = 0;
        if (ana.group_labels.length > 0) {
            for (i=0; i < ana.group_labels.length; i++) {
                data['group_labels'][i]['new_group_label'] = ana.group_labels[i];
                this.group_labels.push(ana.group_labels[i]);
            }
        } else {
            for (i=0; i < data['group_labels'].length; i++) {
                data['group_labels'][i]['new_group_label'] = data['group_labels'][i]['genes'];
                // For the overall labels, do the group number rather than gene since that's what's
                //  displayed by scanpy in the images
                this.group_labels.push(i);
            }
        }

        // show the abbreviated table in the louvain analysis block
        const marker_genes_group_labels_tmpl = $.templates("#marker_genes_group_labels_tmpl");
        const marker_genes_group_labels_html = marker_genes_group_labels_tmpl.render(data['group_labels']);
        $("#marker_genes_group_labels tbody").html(marker_genes_group_labels_html);
        $("#group_labels_c").show();
    }

    populate_marker_genes_table(ana, params) {
        var mg_analysis = this;

        $.ajax({
            type: "POST",
            // TODO: Drop the DPI here. These images are huge.
            url: "./cgi/h5ad_find_marker_genes.cgi",
            data: {'dataset_id': params['dataset_id'], 'analysis_id': params['analysis_id'],
                   'analysis_type': params['analysis_type'], 'session_id': params['session_id'],
                   'n_genes': this.n_genes,
                   'compute_marker_genes': false
                  },
            dataType: "json",
            success: function(data) {
                if (data['success'] == 1) {
                    // add the table header
                    var marker_genes_header_tmpl = $.templates("#marker_genes_table_head_tmpl");
                    var marker_genes_header_html = marker_genes_header_tmpl.render(data['table']['columns']);
                    $("#marker_genes_table thead tr").html("<th>&nbsp;</th>" + marker_genes_header_html);

                    // add the table rows
                    var marker_genes_body_tmpl = $.templates("#marker_genes_table_body_tmpl");
                    var marker_genes_body_html = marker_genes_body_tmpl.render(data['table']['rows']);
                    $("#marker_genes_table tbody").html(marker_genes_body_html);

                    $('#btn_download_marker_genes').show();
                    mg_analysis.populate_marker_genes_labels(ana, data);
                    var group_labels = data['group_labels'].map(x => x.group_label);
                    ana.gene_comparison.populate_group_selectors(group_labels);
                }
            }
        });
    }

    reset() {
        this.calculated = false;
        this.n_genes = false;
        this.reset_ui();
    }

    reset_ui() {
        $("#marker_genes_n_genes").val(5);
    }

    update_ui_with_results(ana, data) {
        var params = {
            'analysis_id': ana.id,
            'analysis_name': 'rank_genes_groups_' + data['cluster_label'],
            'analysis_type': ana.type,
            'dataset_id': ana.dataset_id,
            'session_id': ana.user_session_id,
            // this saves the user from getting a cached image each time
            datetime: (new Date()).getTime()
        }

        $("#marker_genes_table_header").show();
        ana.place_analysis_image(
            {'params': params, 'title': 'Marker genes', 'target': '#marker_genes_plot_c'});

        if (data['table']) {
            // TODO, reduce this into one call from populate_marker_genes_table()
            // add the table header
            var mg_analysis = this;
            var marker_genes_header_tmpl = $.templates("#marker_genes_table_head_tmpl");
            var marker_genes_header_html = marker_genes_header_tmpl.render(data['table']['columns']);
            $("#marker_genes_table thead tr").html("<th>&nbsp;</th>" + marker_genes_header_html);

            // add the table rows
            var marker_genes_body_tmpl = $.templates("#marker_genes_table_body_tmpl");
            var marker_genes_body_html = marker_genes_body_tmpl.render(data['table']['rows']);
            $("#marker_genes_table tbody").html(marker_genes_body_html);

            $('#btn_download_marker_genes').show();
            mg_analysis.populate_marker_genes_labels(ana, data);
            var group_labels = data['group_labels'].map(x => x.group_label);
            ana.gene_comparison.populate_group_selectors(group_labels);
        } else {
            this.populate_marker_genes_table(ana, params);
        }

        $('#marker_genes_visualization_c').show(500);

        // marker gene calculation enables cluster comparison
        $("#toggle_compare_genes").bootstrapToggle('enable');
    }
}

class AnalysisStepPCA {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepPCA();
        if (data === undefined) step;

        step.calculated = data['calculated'];
        step.genes_to_color = data['genes_to_color'];

        if (step.calculated == true) {
            $('#toggle_pca').bootstrapToggle('on');
            $('#pca_genes_to_color').val(step.genes_to_color);
            $('#pca_options_g').show();
            step.update_ui_with_results(ana, data);
        }

        return step;
    }

    reset() {
        this.calculated = false;
        this.genes_to_color = false;
        this.reset_ui();
    }

    reset_ui() {
        $("#pca_genes_to_color").val('');
        $("#top_pca_genes").val('');
    }

    update_ui_with_results(ana, data) {
        $('#analysis_pca .tool_instructions').hide(500);

        var params = {
             analysis_id: ana.id,
             analysis_name: 'pca',
             'analysis_type': ana.type,
             'dataset_id': ana.dataset_id,
             'session_id': ana.user_session_id,
             // this saves the user from getting a cached image each time
             datetime: (new Date()).getTime()
         }

        ana.place_analysis_image(
            {'params': params, 'title': 'PCA scatter', 'target': '#pca_scatter_c'});

        params['analysis_name'] = 'pca_variance_ratio';
        ana.place_analysis_image(
            {'params': params, 'title': 'PCA variance', 'target': '#pca_variance_c'});

        // PCA calculation enables tSNE/UMAP
        $("#toggle_tsne").bootstrapToggle('enable');
    }
}

class AnalysisStepPrimaryFilter {
    constructor() {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepPrimaryFilter();
        if (data === undefined) return step;
        step.calculated = data['calculated'];

        step.filter_cells_gt_n_genes = data['filter_cells_gt_n_genes'];
        step.filter_cells_gt_n_genes_selected = data['filter_cells_gt_n_genes_selected'];
        step.filter_cells_lt_n_genes = data['filter_cells_lt_n_genes'];
        step.filter_cells_lt_n_genes_selected = data['filter_cells_lt_n_genes_selected'];
        step.filter_genes_gt_n_cells = data['filter_genes_gt_n_cells'];
        step.filter_genes_gt_n_cells_selected = data['filter_genes_gt_n_cells_selected'];
        step.filter_genes_lt_n_cells = data['filter_genes_lt_n_cells'];
        step.filter_genes_lt_n_cells_selected = data['filter_genes_lt_n_cells_selected'];
        step.filtered_gene_count = data['filtered_gene_count'];
        step.filtered_cell_count = data['filtered_cell_count'];

        return step;
    }

    reset() {
        this.calculated = false;
        this.filter_cells_gt_n_genes = null
        this.filter_cells_gt_n_genes_selected = false
        this.filter_cells_lt_n_genes = null
        this.filter_cells_lt_n_genes_selected = false
        this.filter_genes_gt_n_cells = null
        this.filter_genes_gt_n_cells_selected = false
        this.filter_genes_lt_n_cells = null
        this.filter_genes_lt_n_cells_selected = false
        this.filtered_gene_count = null
        this.filtered_cell_count = null
        this.reset_ui();
    }

    reset_ui() {
        $("#filter_cells_lt_n_genes").val(300);
        $("#filter_cells_gt_n_genes").val('');
        $("#filter_genes_lt_n_cells").val(3);
        $("#filter_genes_gt_n_cells").val('');

        $('#filter_cells_lt_n_genes_selected').prop('checked', false);
        $('#filter_cells_gt_n_genes_selected').prop('checked', false);
        $('#filter_genes_lt_n_cells_selected').prop('checked', false);
        $('#filter_genes_gt_n_cells_selected').prop('checked', false);

        $(".primary_initial_plot").show();
        $("#primary_initial_plot_c").hide();

    }

    update_ui(ana) {
        if (this.calculated == true) {
            if (this.filter_cells_lt_n_genes_selected == true) {
                $('#filter_cells_lt_n_genes_selected').prop('checked', true);
            } else {
                $('#filter_cells_lt_n_genes_selected').prop('checked', false);
            }

            if (this.filter_cells_gt_n_genes_selected == true) {
                $('#filter_cells_gt_n_genes_selected').prop('checked', true);
            } else {
                $('#filter_cells_gt_n_genes_selected').prop('checked', false);
            }

            if (this.filter_genes_lt_n_cells_selected == true) {
                $('#filter_genes_lt_n_cells_selected').prop('checked', true);
            } else {
                $('#filter_genes_lt_n_cells_selected').prop('checked', false);
            }

            if (this.filter_genes_gt_n_cells_selected == true) {
                $('#filter_genes_gt_n_cells_selected').prop('checked', true);
            } else {
                $('#filter_genes_gt_n_cells_selected').prop('checked', false);
            }

            $('#filter_cells_lt_n_genes').val(this.filter_cells_lt_n_genes);
            $('#filter_cells_gt_n_genes').val(this.filter_cells_gt_n_genes);
            $('#filter_genes_lt_n_cells').val(this.filter_genes_lt_n_cells);
            $('#filter_genes_gt_n_cells').val(this.filter_genes_gt_n_cells);
            $("#selected_dataset_shape_filtered").html(this.filtered_gene_count + " genes x " + this.filtered_cell_count + " obs");
            $("#selected_dataset_shape_filtered_c").show(500);

            var params = {
                'analysis_id': ana.id,
                'analysis_name': 'highest_expr_genes',
                'analysis_type': ana.type,
                'dataset_id': ana.dataset_id,
                'session_id': ana.user_session_id,
                // this saves the user from getting a cached image each time
                datetime: (new Date()).getTime()
            }

            ana.place_analysis_image(
                {'params': params, 'title': 'Highest expressed genes', 'target': '#primary_top_genes_c'});

            $(".primary_initial_plot").hide();
            $("#primary_top_genes_plot_c").show();

            $("#toggle_qc_by_mito").bootstrapToggle('enable');
            $("#toggle_select_variable_genes").bootstrapToggle('enable');
        }
    }
}

class AnalysisStepQCByMito {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepQCByMito();
        if (data === undefined) return step;

        step.calculated = data['calculated'];
        step.gene_prefix = data['gene_prefix'];
        step.filter_mito_perc = data['filter_mito_perc'];
        step.filter_mito_count = data['filter_mito_count'];
        step.n_genes = data['n_genes'];
        step.n_obs = data['n_obs'];

        if (step.calculated == true) {
            $('#toggle_qc_by_mito').bootstrapToggle('on');
            $('#qbm_gene_prefix').val(step.gene_prefix);
            $('#qbm_filter_mito_perc').val(step.filter_mito_perc);
            $('#qbm_filter_mito_count').val(step.filter_mito_count);

            step.update_ui_with_results(ana, data, true);
        }

        return step;
    }

    reset() {
        this.calculated = false;
        this.gene_prefix = null;
        this.filter_mito_perc = null;
        this.filter_mito_count = null;
        this.n_genes = null;
        this.n_obs = null;
        this.reset_ui();
    }

    reset_ui() {
        $("#qbm_gene_prefix").val('mt-');
        $("#btn_qbm_save").hide();
        $("#btn_qbm_save").prop('disabled', true);
        $("#btn_qbm_save").text('Save these genes');

        $("#btn_do_analysis_qc_by_mito").show();
    }

    update_ui_with_results(ana, data, results_saved) {
        if (results_saved == true) {
            $("#btn_do_analysis_qc_by_mito").hide();
            $("#btn_qbm_save").prop('disabled', true);
            $("#btn_qbm_save").text('Saved');
            $("#btn_qbm_save").show();
            $("#qbm_post_shape").show(500);
        }

        $("#analysis_qc_by_mito .tool_instructions").hide(500);
        $("#qbm_gene_count").text(this.n_genes);
        $("#qbm_obs_count").text(this.n_obs);

        var params = {
            'analysis_id': ana.id,
            'analysis_name': 'violin_qc_by_mito',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset_id,
            'session_id': ana.user_session_id,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.place_analysis_image(
            {'params': params, 'title': 'QC by mito - Violin',
             'target': '#qbm_violin_c'});

        params['analysis_name'] = 'scatter_percent_mito'
        ana.place_analysis_image(
            {'params': params, 'title': 'QC by mito - Scatter percent mito',
             'target': '#qbm_scatter_percent_mito_c'});

        params['analysis_name'] = 'scatter_n_genes'
        ana.place_analysis_image(
            {'params': params, 'title': 'QC by mito - Scatter N genes',
             'target': '#qbm_scatter_n_genes_c'});
    }
}

class AnalysisStepSelectVariableGenes {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisStepSelectVariableGenes();
        if (data === undefined ) return step;

        step.calculated = data['calculated'];
        step.norm_counts_per_cell = data['norm_counts_per_cell'];
        step.flavor = data['flavor'];
        step.n_top_genes = data['n_top_genes']
        step.min_mean = data['min_mean'];
        step.max_mean = data['max_mean'];
        step.min_dispersion = data['min_dispersion'];
        step.regress_out = data['regress_out'];
        step.scale_unit_variance = data['scale_unit_variance'];

        if (step.calculated == true) {
            step.update_ui_with_results(ana, data, true);
        }

        return step;
    }

    reset() {
        this.calculated = false;
        this.norm_counts_per_cell = null;
        this.flavor = 'seurat';
        this.n_top_genes = null;
        this.min_mean = null;
        this.max_mean = null;
        this.min_dispersion = null;
        this.regress_out = true;
        this.scale_unit_variance = true;
        this.reset_ui();
    }

    reset_ui() {
        $('#asvg_norm_counts_per_cell').val('1e4');
        $('#asvg_flavor').val('seurat');
        $('#asvg_n_top_genes').val('');
        $('#asvg_min_mean').val(0.0125);
        $('#asvg_max_mean').val(3);
        $('#asvg_min_dispersion').val(0.5);
        $('#asvg_regress_out').prop('checked', true);
        $('#asvg_scale_unit_variance').prop('checked', true);

        $("#btn_asvg_save").hide();
        $("#btn_asvg_save").text('Save these genes');
        $("#btn_asvg_save").prop("disabled", false);
        $("#btn_do_analysis_select_variable_genes").show();
        $('.asvg_save_options').hide();
    }

    update_ui_with_results(ana, data, results_saved) {
        if (this.calculated == true) {
            $('#toggle_select_variable_genes').bootstrapToggle('on');
            $('#asvg_norm_counts_per_cell').val(this.norm_counts_per_cell);
            $('#asvg_flavor').val(this.flavor);
            $('#asvg_n_top_genes').val(this.n_top_genes);
            $('#asvg_min_mean').val(this.min_mean);
            $('#asvg_max_mean').val(this.max_mean);
            $('#asvg_min_dispersion').val(this.min_dispersion);

            if (this.regress_out == true) {
                $('#asvg_regress_out').prop('checked', true);
            } else {
                $('#asvg_regress_out').prop('checked', false);
            }

            if (this.scale_unit_variance == true) {
                $('#asvg_scale_unit_variance').prop('checked', true);
            } else {
                $('#asvg_scale_unit_variance').prop('checked', false);
            }

            $("#btn_asvg_save").show();
            $("#btn_asvg_save").text('Saved');
            $("#btn_asvg_save").prop("disabled", true);
            $("#btn_do_analysis_select_variable_genes").hide();

            // Variable gene calculation enables PCA
            $("#toggle_pca").bootstrapToggle('enable');

            $('.asvg_save_options').show();
        }

        var params = {
            'analysis_id': ana.id,
            'analysis_name': 'filter_genes_dispersion',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset_id,
            'session_id': ana.user_session_id,

            // this saves the user from getting a cached image each time
            datetime: (new Date()).getTime()
        }

        ana.place_analysis_image(
            {'params': params, 'title': 'Variable genes', 'target': '#asvg_plot_c'});

        $("#analysis_select_variable_genes .tool_instructions").hide(500);
    }
}

class AnalysisSteptSNE {
    constructor () {
        this.reset();
    }

    static load_from_json(ana, data) {
        var step = new AnalysisSteptSNE();
        if (data === undefined) return step;

        step.calculated = data['calculated'];
        step.neighbors_calculated = data['neighbors_calculated'];
        step.tsne_calculated = data['tsne_calculated'];
        step.umap_calculated = data['umap_calculated'];
        step.genes_to_color = data['genes_to_color'];
        step.n_pcs = data['n_pcs'];
        step.n_neighbors = data['n_neighbors'];
        step.random_state = data['random_state'];
        step.plot_tsne = data['plot_tsne'];
        step.plot_umap = data['plot_umap'];

        return step;
    }

    reset() {
        this.neighbors_calculated = false;
        this.tsne_calculated = false;
        this.umap_calculated = false;
        this.genes_to_color = false;
        this.n_pcs = null;
        this.n_neighbors = null;
        this.random_state = null;
        this.plot_umap = 0;
        this.plot_tsne = 0;
        this.reset_ui();
    }

    reset_ui() {
        $('#tsne_genes_to_color').val('');
        $('#dredux_n_neighbors').val('');
        $('#tsne_n_pcs').val('');
        $('#tsne_random_state').val(2);
        $('#dimensionality_reduction_method_tsne').prop('checked', false);
        $('#dimensionality_reduction_method_umap').prop('checked', true);
    }

    update_ui(ana) {
        $('#analysis_tsne .tool_instructions').hide(500);

        //if (this.neighbors_calculated == true) {
        if (this.neighbors_calculated || this.tsne_calculated || this.umap_calculated) {
            $('#toggle_tsne').bootstrapToggle('on');
            $('#tsne_genes_to_color').val(this.genes_to_color);
            $('#tsne_n_pcs').val(this.n_pcs);
            $('#dredux_n_neighbors').val(this.n_neighbors);
            $('#tsne_random_state').val(this.random_state);

            if (this.plot_tsne == 1) {
                $('#dimensionality_reduction_method_tsne').prop('checked', true);
            } else {
                $('#dimensionality_reduction_method_tsne').prop('checked', false);
            }

            if (this.plot_umap == 1) {
                $('#dimensionality_reduction_method_umap').prop('checked', true);
            } else {
                $('#dimensionality_reduction_method_umap').prop('checked', false);
            }

            var params = {
                analysis_id: ana.id,
                analysis_name: 'tsne',
                'analysis_type': ana.type,
                'dataset_id': ana.dataset_id,
                'session_id': ana.user_session_id,
                // this saves the user from getting a cached image each time
                datetime: (new Date()).getTime()
            }

            if (this.plot_tsne == 1) {
                ana.place_analysis_image(
                    {'params': params, 'title': 'tSNE plot', 'target': '#tsne_plot_c'});
            }

            if (this.plot_umap == 1) {
                params['analysis_name'] = 'umap'
                ana.place_analysis_image(
                    {'params': params, 'title': 'UMAP plot', 'target': '#umap_plot_c'});
            }

            // dimensionality reduction enables clustering
            $("#toggle_louvain").bootstrapToggle('enable');
        } else {
            console.log("tSNE not calculated");
        }

    }
}
