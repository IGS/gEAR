window.onload=function() {
    // check if the user is already logged in
    check_for_login();
    session_id = Cookies.get('gear_session_id');

    $("#projection_source").on('change', function() {
        $('#projectr_intro_text').hide();
        source = $('#projection_source').val();

        if (source == 'dataset') {
            $('#repository_pattern_c').hide();
            $('#dataset_id_c').show();
            $('#target_dataset_id_c').show();
        } else if (source == 'repository') {
            $('#dataset_id_c').hide();
            $('#target_dataset_id_c').hide();
            $('#repository_pattern_c').show();
        } else {
            $('#dataset_id_c').hide();
            $('#target_dataset_id_c').hide();
            $('#repository_pattern_c').hide();
        }
    });
};

function populate_dataset_selection() {
    $.ajax({
        type: "POST",
        url: "./cgi/get_h5ad_dataset_list.cgi",
        data: {
            'session_id': CURRENT_USER.session_id,
            'for_page': 'projection',
            'include_dataset_id': getUrlParameter('dataset_id')
        },
        dataType: "json",
        success: function(data) {
            if (data['user']['datasets'].length > 0) {
                var user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var user_dataset_list_html = user_dataset_list_tmpl.render(data['user']['datasets']);
                $("#dataset_ids_user").html(user_dataset_list_html);
                $("#target_dataset_ids_user").html(user_dataset_list_html);
            }

            if (data['shared_with_user']['datasets'].length > 0) {
                var shared_with_user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var shared_with_user_dataset_list_html = shared_with_user_dataset_list_tmpl.render(data['shared_with_user']['datasets']);
                $("#dataset_ids_shared_with_user").html(shared_with_user_dataset_list_html);
                $("#target_dataset_ids_shared_with_user").html(shared_with_user_dataset_list_html);
            }

            if (data['public']['datasets'].length > 0) {
                var public_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
                var public_dataset_list_html = public_dataset_list_tmpl.render(data['public']['datasets']);
                $("#dataset_ids_public").html(public_dataset_list_html);
                $("#target_dataset_ids_public").html(public_dataset_list_html);
            }

            // was there a requested dataset ID already?
            var dataset_id = getUrlParameter('dataset_id');
            var share_id = getUrlParameter('share_id');

            if (dataset_id !== undefined) {
                $('#dataset_id').val(dataset_id);
                $( "#dataset_id" ).trigger( "change" );
            }

            // Update the select box
            $('select#dataset_id').selectpicker('refresh');
            $('select#target_dataset_id').selectpicker('refresh');
        },
        error: function(xhr, status, msg) {
            report_error("Failed to load dataset list because msg: " + msg);
        }
    });
}
