"use strict";

var obs_filters = {};
var genes_filters = [];
var supplementary_genes_filters = [];

var dataset_id = null;
var groupbyed_obs = null;
var obs_levels = null;
var gene_symbols = null;    //TODO: get encoded string from POST and decode to array to pre-populate


// Async to ensure data is fetched before proceeding
(async () => {

    // check if the user is already logged in
    await check_for_login();
    session_id = Cookies.get("gear_session_id");

    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    })

    $('#dataset_select').select2({
        placeholder: 'To search, click to select or start typing a dataset name',
    });
    await populate_datasets();

    // Initialize plot types
    $('#plot_type_select').select2({
        placeholder: 'Choose how to plot',
    });
    $('#plot_type_select').hide();

    $("#dataset_select").change(async function() {
        dataset_id = $("#dataset_select").select2('data')[0].id;

        $('#plot_type_select').show();

        // Get genes for this dataset
        gene_symbols = await fetch_gene_symbols({dataset_id, undefined});
        create_gene_dropdown(gene_symbols);

        // Get categorical observations for this dataset
        const data = await fetch_h5ad_info({dataset_id, undefined});
        obs_levels = curate_observations(data.obs_levels);
        create_obs_groupby_field(obs_levels);
        create_obs_dropdowns(obs_levels);

        // Ensure genes and observation columns dropdown toooltip shows
        $(function () {
            $('[data-toggle="tooltip"]').tooltip()
        })

    })

    $(document).on('click', "#update_plot", function() {

        // Remove supplementary plot and reset its genes filter
        if (supplementary_genes_filters.length) {
            supplementary_genes_filters = []
            $('#supplementary_plot').remove()
        }

        // Render dataset plot HTML
        const plot_template = $.templates("#dataset_plot_tmpl");
        const plot_html = plot_template.render({dataset_id: dataset_id});
        $('#dataset_plot').html(plot_html);

        var plot_type = $("#plot_type_select").select2('data')[0].id;

        var groupby_filter = $('input[name="obs_groupby"]:checked').val();

        // Update filters based on selection
        obs_filters = {};
        for (var property in obs_levels) {
            var prop_data = $(`#${property}_dropdown`).select2('data');
            obs_filters[property] = prop_data.map(function(elem) {
                return elem.id;
            });

            // If no groups for an observation are selected, delete filter
            if (!obs_filters[property].length) {
                delete obs_filters[property];
            }
        }

        if (!Object.keys(obs_filters).length) {
            alert("At least one observation must have categories filtered.");
            return;
        }

        genes_filters = $("#gene_dropdown").select2('data').map(function(elem) {
            return elem.id;
        })

        if (!genes_filters.length) {
            alert("At least one gene must be provided.");
            return;
        }

	    const cluster_cols = $('#cluster_cols').is(':checked');
	    const payload = { cluster_cols, groupby_filter };

        // Draw the updated chart
        draw(dataset_id, plot_type, genes_filters, obs_filters, payload);

    });

    // Cannot groupby if observations are going to be clustered
    $('#cluster_cols').change(function(){
        if(this.checked) {
            $('.obs_groupby').prop('checked', false);
            $('#obs_groupby_container').hide();
        } else {
            $('#obs_groupby_container').show();
        }
    });

    // Some options are specific to certain plot types
    $('#plot_type_select').change(function(){
        switch ($('#plot_type_select').val()) {
            case "heatmap":
                $('#obs_checkbox_container').show();
                $('#obs_groupby_container').show();
                break;
            case "violin":
                $('#cluster_cols').prop('checked', false);
                $('#obs_checkbox_container').hide();
                $('#obs_groupby_container').show();
                break;
            default:
                $('#cluster_cols').prop('checked', false);
                $('#obs_checkbox_container').hide();
                $('.obs_groupby').prop('checked', false);
                $('#obs_groupby_container').hide();
        }

    });

    // If "all" button is clicked, populate dropdown with all groups in this observation
    $(document).on('click', ".all", function() {
        var id = this.id;
        var group = id.replace("_all", "");

        $(`#${group}_dropdown`).val(obs_levels[group]);
        $(`#${group}_dropdown`).trigger('change');  // This actually triggers select2 to show the dropdown vals
    });

    // If gene is clicked in plot display supplementary violin plot
    $(document).on('click', "g.y5tick text a", function() {
        var gene = $(this).text();

        // Add or remove gene depending on if it is already in array
        const index = supplementary_genes_filters.indexOf(gene);
        if (index === -1) {
            supplementary_genes_filters.push(gene);
            $(this).parent().css("fill", "crimson");
        } else {
            supplementary_genes_filters.splice(index, 1);
            $(this).parent().css("fill", "rgb(42, 63, 95)"); // original default fill color
        }

        // Render supplementary plot HTML
        if (supplementary_genes_filters.length) {
            // Render supplementary plot HTML
            const plot_template = $.templates("#supplementary_plot_tmpl");
            const plot_html = plot_template.render({dataset_id: dataset_id});
            $('#supplementary_plot').html(plot_html);
            // Draw the supplementary chart
            var groupby_filter = $('input[name="obs_groupby"]:checked').val();
            var payload = {groupby_filter}
            draw(dataset_id, "violin", supplementary_genes_filters, obs_filters, payload, true);
        }
    });

    $(document).on('click', "#reset_obs", async function() {
        // Get categorical observations for this dataset
        const data = await fetch_h5ad_info({dataset_id, undefined});
        obs_levels = curate_observations(data.obs_levels);
        create_obs_groupby_field(obs_levels);
        create_obs_dropdowns(obs_levels);

        // Update list of observations to groupby by
        create_obs_groupby_field(obs_levels);
        create_obs_dropdowns(obs_levels);
    });

})();


async function get_data(dataset_id, plot_type, genes_filters, obs_filters, payload) {
    return await axios.post(`/api/plot/${dataset_id}/mg_dash`, {
        ...payload,
        plot_type: plot_type,
        gene_symbols: genes_filters,
        obs_filters: obs_filters,
    });
}

async function fetch_gene_symbols(payload) {
    const { dataset_id, analysis } = payload;
    const base = `./api/h5ad/${dataset_id}/genes`;
    const query = analysis ? `?analysis=${analysis.id}` : "";

    const { data } = await axios.get(`${base}${query}`);
    return data.gene_symbols;
}

async function fetch_h5ad_info(payload) {
    const { dataset_id, analysis } = payload;
    const base = `./api/h5ad/${dataset_id}`;
    const query = analysis ? `?analysis=${analysis.id}` : "";
    const { data } = await axios.get(`${base}${query}`);
    return data;
}

async function populate_datasets() {
    $.ajax({
        type: "POST",
        url: "./cgi/get_h5ad_dataset_list.cgi",
        data: {
            session_id: undefined,  // TODO: define with user login info
        },
        dataType: "json",
        success: function (data) {
            // Populate select box with dataset information owned by the user
            if (data["user"]["datasets"].length > 0) {
                var user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var user_dataset_list_html = user_dataset_list_tmpl.render(
                    data["user"]["datasets"]
                );
                $("#dataset_ids_user").html(user_dataset_list_html);
            } else {
                $("#dataset_id .user_initial").html("Not logged in");
            }

            // Next, add datasets shared with the user
            if (data["shared_with_user"]["datasets"].length > 0) {
                var shared_with_user_dataset_list_tmpl = $.templates(
                    "#dataset_list_tmpl"
                );
                var shared_with_user_dataset_list_html = shared_with_user_dataset_list_tmpl.render(
                    data["shared_with_user"]["datasets"]
                );
                $("#dataset_ids_shared_with_user").html(
                    shared_with_user_dataset_list_html
                );
            }

            // Now, add public datasets
            if (data["public"]["datasets"].length > 0) {
                var public_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var public_dataset_list_html = public_dataset_list_tmpl.render(
                    data["public"]["datasets"]
                );
                $("#dataset_ids_public").html(public_dataset_list_html);
            }

        },
        error: function (xhr, status, msg) {
            console.error("Failed to load dataset list because msg: " + msg);
        },
        });
}

// Draw plotly chart in HTML
function draw_chart(data, dataset_id, supplementary=false) {
    const target_div = supplementary ? `dataset_${dataset_id}_secondary` : `dataset_${dataset_id}_h5ad`;
    const {
        plot_json,
        plot_config
    } = data;

    var layout_mods = {
        height: target_div.clientHeight,
        width: target_div.clientWidth,
    };

    // Overwrite plot layout and config values with custom ones from display
    var layout = {
        ...plot_json.layout,
        ...layout_mods,
    };

    var config_mods = {
        responsive: false,
    };

    const config = {
        ...plot_config,
        ...config_mods,
    };
    Plotly.newPlot(target_div, plot_json.data, layout, config);
}

// Submit API request and draw the HTML
async function draw(dataset_id, plot_type, genes_filters, obs_filters, payload, supplementary=false) {
    const {
        data
    } = await get_data(dataset_id, plot_type, genes_filters, obs_filters, payload);
    draw_chart(data, dataset_id, supplementary);
}

function create_gene_dropdown(genes) {
    var tmpl = $.templates("#gene_dropdown_tmpl");
    var data = {genes: genes}
    var html = tmpl.render(data);
    $("#gene_dropdown_container").html(html);
    $('select.genes').select2({
        placeholder: 'To search, click to select or start typing some gene names',
        allowClear: true
    });
}

function create_obs_groupby_field(obs_levels) {
    var tmpl = $.templates("#obs_groupby_tmpl");  // Get compiled template using jQuery selector for the script block
    var html = tmpl.render(obs_levels);                 // Render template using data - as HTML string
    $("#obs_groupby_container").html(html);       // Insert HTML string into DOM
}

function create_obs_dropdowns(obs_levels) {
    var tmpl = $.templates("#obs_dropdowns_tmpl");
    var html = tmpl.render(obs_levels);
    $("#obs_dropdowns_container").html(html);
    $('select.obs_levels').select2({
        placeholder: 'Start typing to filter categories. Click "All" to use all categories',
        allowClear: true
    });
}

function curate_observations(obs_levels) {
    // Delete useless filters
    for (var property in obs_levels) {
        if (property === "color" || property.endsWith("colors")) {
            delete obs_levels[property];
        } else if (obs_levels[property].length === 1) {
            delete obs_levels[property];
        }
    }
    return obs_levels
}
