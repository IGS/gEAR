//const dataset_id = "b0420910-a0fa-e920-152d-420b6275d3af";

var obs_filters = {};
var genes_filters = [];

var dataset_id = null;
var sorted_obs = null;
var obs_levels = null;
var gene_symbols = null;    //TODO: get encoded string from POST and decode to array to pre-populate

// Async to ensure data is fetched before proceeding
(async () => {

    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    })

    $('#dataset_select').select2({
        placeholder: 'To search, click to select or start typing a dataset name',
    });
    datasets = await populate_datasets()

    // Initialize plot types
    $('#plot_type_select').select2({
        placeholder: 'Choose how to plot',
    });
    $('#plot_type_select').hide();

    $("#dataset_select").change(async function() {
        dataset_id = $("#dataset_select").select2('data')[0].id;

        dataset_id = "b0420910-a0fa-e920-152d-420b6275d3af";

        $('#plot_type_select').show();

        // Get genes for this dataset
        gene_symbols = await fetch_gene_symbols({dataset_id, undefined});
        make_gene_dropdown(gene_symbols);
        $('select.genes').select2({
            placeholder: 'To search, click to select or start typing some gene names',
            allowClear: true
        });

        // Get categorical observations for this dataset
        const data = await fetch_h5ad_info({dataset_id, undefined});
        obs_levels = data.obs_levels;
        // Delete useless filters
        for (var property in obs_levels) {
            if (property === "color" || property.endsWith("colors")) {
                delete obs_levels[property];
            } else if (obs_levels[property].length === 1) {
                delete obs_levels[property];
            }
        }
        make_obs_sort_dropdown(obs_levels);
        make_obs_dropdowns(obs_levels);
        $('select.obs_levels').select2({
            placeholder: 'Start typing to filter categories or leave empty for no filter on this category. Click "x" to remove filter.',
            allowClear: true
        });

        // Ensure genes and observation columns dropdown toooltip shows
        $(function () {
            $('[data-toggle="tooltip"]').tooltip()
        })

        // Render dataset plot HTML
        const plot_template = $.templates("#dataset_plot_tmpl");
        const plot_html = plot_template.render({dataset_id: dataset_id});
        $('#dataset_plot').html(plot_html);

    })

    // If close button is clicked for a obs group, remove the property from list of observations, and hide the div
    $(".close").click(function() {
        var id = this.id;
        var group = id.replace("_close", "");

        delete obs_levels[group];
        $(`${group}_div`).hide();

        // Update list of observations to sort by
        make_obs_sort_dropdown(obs_levels);
    });

    $("#update_plot").click(function() {

        var plot_type = $("#plot_type_select").select2('data')[0].id;

        var sort_filter = $('input[name="obs_sort"]:checked').val();

        // Update filters based on selection
        obs_filters = {};
        for (var property in obs_levels) {
            prop_data = $(`#${property}_dropdown`).select2('data');
            obs_filters[property] = prop_data.map(function(elem) {
                return elem.id;
            });

            // If no groups for an observation are selected, use all of them
            if (!obs_filters[property].length) {
                obs_filters[property] = obs_levels[property];
            }
        }

        genes_filters = $("#gene_dropdown").select2('data').map(function(elem) {
            return elem.id;
        })

	const cluster_cols = $('#cluster_cols').is(':checked');
	const payload = { cluster_cols, sort_filter };

        // Draw the updated chart
        draw(dataset_id, plot_type, genes_filters, obs_filters, payload);

    });

    // Cannot sort if observations are going to be clustered
    $('#cluster_cols').change(function(){
        if(this.checked) {
            $('.obs_sort').prop('checked', false);
            $('#obs_sort_dropdown_container').hide();
        } else {
            $('#obs_sort_dropdown_container').show();
        }
    });

    // Some options are heatmap-specific
    $('#plot_type_select').change(function(){
        if($('#plot_type_select').val() == "heatmap") {
            $('#obs_checkbox_container').show();
            $('#obs_sort_dropdown_container').show();
        } else {
            $('#cluster_cols').prop('checked', false);
            $('#obs_checkbox_container').hide();
            $('.obs_sort').prop('checked', false);
            $('#obs_sort_dropdown_container').hide();
        }
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

function draw_chart(data, dataset_id) {

    const target_div = `dataset_${dataset_id}_h5ad`;
    const {
        plot_json,
        plot_config
    } = data;

    $(`#dataset_${dataset_id}`).append(template());

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

function template() {
    const template = `
    <div
      id='dataset_${dataset_id}_h5ad'
      class="h5ad-container"
      style="position: relative;">
    </div>
  `;
    return template;
}

async function draw(dataset_id, plot_type, genes_filters, obs_filters, payload) {
    const {
        data
    } = await get_data(dataset_id, plot_type, genes_filters, obs_filters, payload);
    draw_chart(data, dataset_id);
}

function make_gene_dropdown(genes) {
    var tmpl = $.templates("#gene_dropdown_tmpl"); // Get compiled template using jQuery selector for the script block
    var data = {genes: genes}
    var html = tmpl.render(data);             // Render template using data - as HTML string
    $("#gene_dropdown_container").html(html);           // Insert HTML string into DOM
}

function make_obs_sort_dropdown(obs_levels) {
    var tmpl = $.templates("#obs_sort_dropdown_tmpl"); // Get compiled template using jQuery selector for the script block
    var html = tmpl.render(obs_levels);             // Render template using data - as HTML string
    $("#obs_sort_dropdown_container").html(html);           // Insert HTML string into DOM
}

function make_obs_dropdowns(obs_levels) {
    var tmpl = $.templates("#obs_dropdowns_tmpl"); // Get compiled template using jQuery selector for the script block
    var html = tmpl.render(obs_levels);             // Render template using data - as HTML string
    $("#obs_dropdowns_container").html(html);           // Insert HTML string into DOM
}
