// TODO
//  - Make existing plots disappear right when the user hits 'Plot' to redraw

plot_data = null;
selected_data = null;

window.onload=function() {
    // check if the user is already logged in
    check_for_login();
    session_id = Cookies.get('gear_session_id');

    $( ".btn-apply-filter" ).on('click', function() {
        $('.initial_instructions').hide();
        $('#myChart').html('');
        $('#error_loading_c').hide();
        $('#plot_loading').show();
        load_comparison_graph();
    });

    /***** gene cart stuff *****/
    $('#create_gene_cart').on('click', function () {
        $('#create_gene_cart_dialog').show('fade');
    });
    $('#cancel_save_gene_cart').on('click', function () {
        $('#create_gene_cart_dialog').hide('fade');
        $('#gene_cart_name').val('');
    });
    $("#gene_cart_name").on('input', function() {
        if ($(this).val() == '') {
            $("#save_gene_cart").prop("disabled", true);
        } else {
            $("#save_gene_cart").prop("disabled", false);
        }
    });
    $("#save_gene_cart").on('click', function () {
        $("#save_gene_cart").prop("disabled", true);

        if (CURRENT_USER) {
            save_gene_cart();
        } else {
            alert("You must be signed in to do that.");
        }

    });
    /***** end gene cart stuff *****/

    $( "#dataset_id" ).on('change', function() {
        populate_condition_selection_control();
    });

    $("#statistical_test").on('change', function() {
        if ($("#statistical_test").val()) {
            $("#test_pval_cutoff").prop("disabled", false);
        } else {
            $("#test_pval_cutoff").prop("disabled", true);
        }
    });

    // initially disable the condition selectors
    $('#dataset1_conditions').attr('disabled', 'disabled');
    $('#dataset2_conditions').attr('disabled', 'disabled');
}

function download_selected_genes() {
    // Builds a file in memory for the user to download.  Completely client-side.
    // plot_data contains three keys: x, y and symbols
    // build the file string from this
    if ( $('#log_base').val() == 'raw' ) {
        file_contents = "gene_symbol\t" + $('#dataset1_conditions').val() + "\t" +
            $('#dataset2_conditions').val() + "\n";
    } else {
        file_contents = "gene_symbol\t" + $('#dataset1_conditions').val() + " (log" + $('#log_base').val() + ")\t" +
            $('#dataset2_conditions').val() + " (log" + $('#log_base').val() +")\n";
    }

    selected_data.points.forEach(function(pt) {
        // Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
        file_contents += plot_data['symbols'][pt.pointNumber] + "\t" +
                         pt.x.toFixed(1) + "\t" + pt.y.toFixed(1) + "\n";
    });

    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(file_contents));
    element.setAttribute('download', 'selected_genes.xls');
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

function load_comparison_graph() {
    $.ajax({
        url : './cgi/get_dataset_comparison.cgi',
        type: "POST",
        data: { 'dataset1_id': $('#dataset_id').val(),
                'dataset1_condition': $('#dataset1_conditions').val(),
                'dataset2_condition': $('#dataset2_conditions').val(),
                'fold_change_cutoff': $('#fold_change_cutoff').val(),
                'std_dev_num_cutoff': $('#std_dev_num_cutoff').val(),
                'log_transformation': $('#log_base').val(),
                'statistical_test'  : $('#statistical_test').val(),
                'test_pval_cutoff  ': $('#test_pval_cutoff').val()
              },
        dataType: "json",
        success: function(data, textStatus, jqXHR) {
            if ( data['success'] == 1 ) {
                $('#fold_change_std_dev').html(data['fold_change_std_dev']);
                plot_data = data
                plot_data_to_graph(data)
            } else {
                console.log("CGI unsuccessful");
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            // Handle graphing failures
            $('#plot_loading').hide();
            $('#ticket_datasetx_id').text($('#dataset_id').val());
            $('#ticket_datasetx_condition').text($('#dataset1_conditions').val());
            $('#ticket_datasety_id').text($('#dataset_id').val());
            $('#ticket_datasety_condition').text($('#dataset2_conditions').val());
            $('#error_loading_c').show();
        }
    });
}

function populate_condition_selection_control() {
    dataset_id = $('#dataset_id').val();
    $('#dataset1_conditions').attr('disabled', 'disabled');
    $('#dataset2_conditions').attr('disabled', 'disabled');
    $('#dataset1_conditions').html("<option>Loading ... </option>");
    $('#dataset2_conditions').html("<option>Loading ... </option>");

    $.ajax({
        url : './cgi/get_condition_list.cgi',
        type: "POST",
        data : { 'dataset_id': dataset_id },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if (data['success'] === 0) {
                display_error_bar(data['error']);

            } else if (data['success'] === -1) {
                $('#dataset1_conditions > option')
                    .html('This dataset is not ready to compare.')
                    .attr('selected', 'selected');
            } else {
                var selectorTmpl = $.templates('#dataset_condition_options');
                var selectorHtml = selectorTmpl.render(data['conditions']);
                $('#dataset1_conditions').html(selectorHtml);
                $('#dataset2_conditions').html(selectorHtml);

                $('#dataset1_conditions').removeAttr('disabled');
                $('#dataset2_conditions').removeAttr('disabled');

                if (data['has_replicates'] == 1) {
                    $("#statistical_test_label").html("");
                    $("#statistical_test").attr("disabled", false);
                } else {
                    $("#statistical_test_label").html("Not applicable since this dataset has no replicates");
                    $("#statistical_test").attr("disabled", true);
                }
            }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
	          console.log('textStatus= ', textStatus);
	          console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}

function populate_dataset_selection_controls() {
    $.ajax({
        type: "POST",
        url: "./cgi/get_h5ad_dataset_list.cgi",
        data: {
            'session_id': CURRENT_USER.session_id,
            'for_page': 'compare_dataset',
            'include_dataset_id': getUrlParameter('dataset_id')
        },
        dataType: "json",
        success: function(data) {
            if (data['user']['datasets'].length > 0) {
                var user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var user_dataset_list_html = user_dataset_list_tmpl.render(data['user']['datasets']);
                $("#dataset_ids_user").html(user_dataset_list_html);
            } else {
                $("#dataset_id .user_initial").html("Not logged in");
            }

            if (data['shared_with_user']['datasets'].length > 0) {
                var shared_with_user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var shared_with_user_dataset_list_html = shared_with_user_dataset_list_tmpl.render(data['shared_with_user']['datasets']);
                $("#dataset_ids_shared_with_user").html(shared_with_user_dataset_list_html);
            }

            if (data['public']['datasets'].length > 0) {
                var public_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var public_dataset_list_html = public_dataset_list_tmpl.render(data['public']['datasets']);
                $("#dataset_ids_public").html(public_dataset_list_html);
            }

            // was there a requested dataset ID already?
            var dataset_id = getUrlParameter('dataset_id');
            if (dataset_id !== undefined) {
                $('#dataset_id').val(dataset_id);
                $( "#dataset_id" ).trigger( "change" );
            }
        },
        error: function(xhr, status, msg) {
            report_error("Failed to load dataset list because msg: " + msg);
        }
    });
}

function plot_data_to_graph(data) {
    $('#tbl_selected_genes').hide();
    $('#selection_methods_c').show();

    var point_labels = [];
    var perform_ranking = false;

    if ( $("#statistical_test").val()) {
        perform_ranking = true;
    }

    var plotdata = null;

    if (perform_ranking) {
        var pval_cutoff = parseFloat($("#test_pval_cutoff").val());
        var passing = {'x': [], 'y': [], 'labels': [], 'id':[], 'pvals':[]};
        var failing = {'x': [], 'y': [], 'labels': [], 'id':[], 'pvals':[]};

        for (i = 0; i < data['x'].length; ++i) {
            var this_pval = parseFloat(data['pvals_adj'][i]);

            if (this_pval <= pval_cutoff) {
                // good scoring match
                passing['x'].push(data['x'][i]);
                passing['y'].push(data['y'][i]);
                passing['labels'].push("Gene symbol: " + data['symbols'][i] + "   P-value: " + this_pval.toPrecision(6));
                passing['id'].push(data['symbols'][i]);
                passing['pvals'].push(data['pvals_adj'][i]);

            } else {
                // this one didn't pass the p-value cutoff
                failing['x'].push(data['x'][i]);
                failing['y'].push(data['y'][i]);
                failing['labels'].push("Gene symbol: " + data['symbols'][i] + "   P-value: " + this_pval.toPrecision(6));
                failing['id'].push(data['symbols'][i]);
                failing['pvals'].push(data['pvals_adj'][i]);

            }
        }

        if ($("input[name='stat_action']:checked").val() == 'colorize') {
            plotdata = [
                {
                    id: passing['id'],
                    pvals: passing['pvals'],
                    x: passing['x'],
                    y: passing['y'],
                    mode: 'markers',
                    name: "Passed cutoff",
                    type: 'scatter',
                    text: passing['labels'],
                    marker: {
                        color: '#FF0000',
                        size: 4
                    }
                },
                {
                    id: failing['id'],
                    pvals: failing['pvals'],
                    x: failing['x'],
                    y: failing['y'],
                    mode: 'markers',
                    name: "Did not pass cutoff",
                    type: 'scatter',
                    text: failing['labels'],
                    marker: {
                        color: '#A1A1A1',
                        size: 4
                    }
                }
            ];
        } else {
            plotdata = [
                {
                    id: passing['id'],
                    pvals: passing['pvals'],
                    x: passing['x'],
                    y: passing['y'],
                    mode: 'markers',
                    type: 'scatter',
                    text: passing['labels'],
                    marker: {
                        color: '#2F103E',
                        size: 4
                    }
                }
            ];
        }
    } else {
        for (i = 0; i < data['symbols'].length; ++i) {
            point_labels.push("Gene symbol: " + data['symbols'][i]);
        }

        plotdata = [
            {
                id: data['symbols'],
                pvals: data['pvals_adj'],
                x: data['x'],
                y: data['y'],
                mode: 'markers',
                type: 'scatter',
                text: point_labels,
                marker: {
                    color: '#2F103E',
                    size: 4
                }
            }
        ];
    }

    var layout = {
        title: '',
        xaxis: {title: $('#dataset_id option:selected').text() + ' - ' + $('#dataset1_conditions option:selected').text(),
                type: ''},
        yaxis: {title: $('#dataset_id option:selected').text() + ' - ' + $('#dataset2_conditions option:selected').text(),
                type: ''},
        margin: {t: 20},
        hovermode: 'closest',
        dragmode: 'select'
    }

    $('#plot_loading').hide();
    var graphDiv = document.getElementById('myChart');
    Plotly.newPlot(graphDiv, plotdata, layout, {showLink: false});
    $('#selected_label').hide();
    $('#controls_label').show();

    graphDiv.on('plotly_selected', function(eventData) {
        selected_data = eventData;
        selected_gene_data = [];
        eventData.points.forEach(function(pt) {
            // Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
            // Each trace has its own "pointNumber" ids so gene symbols and pvalues needed to be passed in for each plotdata trace
            selected_gene_data.push({'gene_symbol': pt.data.id[pt.pointNumber],
                                     'pval_adj': pt.data.pvals[pt.pointNumber],
                                     'x': pt.x.toFixed(1), 'y': pt.y.toFixed(1)});
        });
        var template = $.templates('#selected_genes_tmpl');
        var htmlOutput = template.render(selected_gene_data);
        $('#selected_genes_c').html(htmlOutput);

        // toggle visibilities
        $('.selection_instructions').hide();
        $('#saved_gene_cart_info_c').hide();
        $('#tbl_selected_genes').show();

        $('#controls_label').hide();
        $('#selected_label').show();

        if ( $('#log_base').val() == 'raw' ) {
            $('#tbl_selected_genes_transformation_row').hide();
        } else {
            $('#table_transformation_label').text('Log' + $('#log_base').val());
            $('#tbl_selected_genes_transformation_row').show();
        }
    });

    window.onresize = function() {
        Plotly.Plots.resize(graphDiv);
    }
};

function save_gene_cart() {
    // must have access to USER_SESSION_ID
    var gc = new GeneCart({'session_id': CURRENT_USER.session_id, 'label': $('#gene_cart_name').val()});

    selected_data.points.forEach(function(pt) {
        var gene = new Gene({'id': plot_data['gene_ids'][pt.pointNumber], 'gene_symbol': plot_data['symbols'][pt.pointNumber]});
        gc.add_gene(gene);
    });

    gc.save(update_ui_after_gene_cart_save);
}

function update_ui_after_gene_cart_save(gc) {
    $('#create_gene_cart_dialog').hide('fade');
    $('#saved_gene_cart_info_c > h3').html('Cart: ' + gc.label);
    $('#gene_cart_member_count').html(gc.genes.length);
    $('#saved_gene_cart_info_c').show();
}
