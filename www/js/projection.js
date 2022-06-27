let SELECTED_PATTERN_LABEL = null;

window.onload=() => {
    $("#projection_source").on('change', () => {
        $('#projectr_intro').hide();
        source = $('#projection_source').val();

        if (source == 'dataset') {
            $('#set_of_patterns_c').hide();
            $('#dataset_id_c').show();
            $('#target_dataset_id_c').show();
        } else if (source == 'repository') {
            $('#dataset_id_c').hide();
            $('#target_dataset_id_c').hide();
            $('#set_of_patterns_c').show();
        } else {
            $('#dataset_id_c').hide();
            $('#target_dataset_id_c').hide();
            $('#set_of_patterns_c').hide();
        }
    });

    $("#set_of_patterns").on('change', () => {
        $.ajax({
            type: "POST",
            url: "./cgi/get_pattern_element_list.cgi",
            data: {
                'file_name': $('#set_of_patterns').val(),
            },
            dataType: "json",
            success: (data) => {
                const pattern_elements_tmpl = $.templates("#pattern_elements_tmpl");
                const pattern_elements_html = pattern_elements_tmpl.render(data);
                $("#projection_pattern_elements").html(pattern_elements_html);
                $("#projection_pattern_elements_c").show();
            },
            error(xhr, status, msg) {
                report_error(`Failed to load dataset list because msg: ${msg}`);
            }
        });
    });

    $("#target_dataset_id").on('change', () => {
        $('#submitter_c').show();
    });

    $(document).on('click', 'li.projection_pattern_element', function(e) {
        // set only this one as active
        $("li.projection_pattern_element").removeClass("active");
        $(this).addClass("active");

        SELECTED_PATTERN_LABEL = $(this).data('label');

        $("#target_dataset_id_c").show();
        $("#projection_c").prop("disabled", false);
        $("#projection_c").show();
    });

    $(document).on('click', '#btn_project', async function(e) {
        e.preventDefault();
        $(this).attr("disabled", true);
        const datasetId = $("#target_dataset_id").val();
        const config = {
            "scope":source
            , "input_value": $('#set_of_patterns').val()
            , "pattern_value": SELECTED_PATTERN_LABEL
        }

       /* const {default_display_id: defaultDisplayId} = await getDefaultDisplay(datasetId);
        const display = await $.ajax({
            url: './cgi/get_dataset_display.cgi',
            type: 'POST',
            data: { display_id: defaultDisplayId },
            dataType: 'json'
        });

        plotConfig = display.plotly_config;

        config.plot_config = plotConfig;
        */

        const { data } = await runProjectR(datasetId, config);
        $(this).attr("disabled", false);
        drawChart(data, datasetId);

    });

};

// Call API to return plot JSON data
async function runProjectR (datasetId, payload) {
    try {
        return await axios.post(`/api/projectr/${datasetId}`, {
            ...payload
        })
    } catch (e) {

        const message = "There was an error in making this projection. Please contact the gEAR team using the 'Contact' button at the top of the page and provide as much information as possible.";
        const success = -1;
        const data = {message, success};
        return {data};
    }
}

// Draw plotly chart in HTML
function drawChart (data, datasetId) {
    const targetDiv = "projection_view";
    const { plot_json: plotlyJson, plot_config: plotlyConfig, message, success } = data;

    // Since default plots are now added after dataset selection, wipe the plot when a new one needs to be drawn
    $(`#${targetDiv}`).empty()

    const configMods = {
        responsive: true
    };

    const config = {
        ...plotlyConfig,
        ...configMods
    };
    Plotly.newPlot(targetDiv, plotlyJson.data, plotlyJson.layout, config);
}

function populate_dataset_selection() {
    // NOTE: Called in common.js
    $.ajax({
        type: "POST",
        url: "./cgi/get_h5ad_dataset_list.cgi",
        data: {
            'session_id': CURRENT_USER.session_id,
            'for_page': 'projection',
            'include_dataset_id': getUrlParameter('dataset_id')
        },
        dataType: "json",
        success: (data) => {
            if (data['user']['datasets'].length > 0) {
                const user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                const user_dataset_list_html = user_dataset_list_tmpl.render(data['user']['datasets']);
                $("#dataset_ids_user").html(user_dataset_list_html);
                $("#target_dataset_ids_user").html(user_dataset_list_html);
            }

            if (data['shared_with_user']['datasets'].length > 0) {
                const shared_with_user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                const shared_with_user_dataset_list_html = shared_with_user_dataset_list_tmpl.render(data['shared_with_user']['datasets']);
                $("#dataset_ids_shared_with_user").html(shared_with_user_dataset_list_html);
                $("#target_dataset_ids_shared_with_user").html(shared_with_user_dataset_list_html);
            }

            if (data['public']['datasets'].length > 0) {
                const public_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                const public_dataset_list_html = public_dataset_list_tmpl.render(data['public']['datasets']);
                $("#dataset_ids_public").html(public_dataset_list_html);
                $("#target_dataset_ids_public").html(public_dataset_list_html);
            }

            // was there a requested dataset ID already?
            const dataset_id = getUrlParameter('dataset_id');

            if (dataset_id !== undefined) {
                $('#dataset_id').val(dataset_id);
                $( "#dataset_id" ).trigger( "change" );
            }
        },
        error(xhr, status, msg) {
            report_error(`Failed to load dataset list because msg: ${msg}`);
        }
    });
}

function populate_pattern_selection() {
    // NOTE: Called in common.js
    $.ajax({
        type: "POST",
        url: "./cgi/get_projection_pattern_list.cgi",
        data: {
            'session_id': CURRENT_USER.session_id,
        },
        dataType: "json",
        success: (data) => {
            const pattern_list_tmpl = $.templates("#dataset_list_tmpl");    // recycling the template... same output
            const pattern_list_html = pattern_list_tmpl.render(data);
            $("#set_of_patterns").html(pattern_list_html);
        },
        error(xhr, status, msg) {
            report_error(`Failed to load dataset list because msg: ${msg}`);
        }
    });
}

function getDefaultDisplay (datasetId) {
    return $.ajax({
        url: './cgi/get_default_display.cgi',
        type: 'POST',
        data: {
            user_id: CURRENT_USER.id,
            dataset_id: datasetId,
            is_multigene: 0 //TODO: Change
        },
        dataType: 'json'
    });
}

function report_error(msg) {
    if (msg) {
        $(`<li class='failure'>${msg}</li>`).prependTo("#action_log");
    }
}
