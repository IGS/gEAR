var animation_time = 200;

window.onload=function() {
    // check if the user is already logged in
    check_for_login();

    $("#initial_instructions_bar").on('click', function() {
        if ($("#initial_instructions_body").is(":visible")) {
            $("#initial_instructions_body").hide(animation_time);
        } else {
            $("#initial_instructions_body").show(animation_time);
        }
    });

    $("#initial_instructions_closer i").on('click', function() {
        $("#initial_instructions_c").hide();
    });
};  // end window onloads

function display_error_bar(msg) {
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' +
      '<strong>Fail. </strong> Sorry, something went wrong.  Please contact us with this message if you need help.' +
      '</p>' +
      '<p style="text-align: center;">(<em>Error: ' + msg + '</em>)</p>' +
      '</div>').show();
}


