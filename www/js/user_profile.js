/*
 This script relies on the source having also included the
 common.js within this project
*/

window.onload=function() {
    sleep(500).then(() => {
        // If CURRENT_USER is defined at this point, add information as placeholder test
        if (CURRENT_USER) {
            $('#email').attr('placeholder', CURRENT_USER.email);
            $('#institution').attr('placeholder', CURRENT_USER.institution);
            $('#wantUpdates').prop('checked', false);
            if (CURRENT_USER.updates_wanted == 1) {
                $('#wantUpdates').prop('checked', true);
            }
        }
    })

};

$(document).on('keyup', '#newPassword, #repeatPassword', function () {
    // Originally was adding "match/nomatch" classes to #passwordMatch
    // but for some reason, class does not want to be removed once added
    if ($('#newPassword').val() == $('#repeatPassword').val()) {
        $('#passwordMatch').html('Matching').css('color', 'green');
        $("#passwordInvalidTooltip").remove();
        $('#newPassword').removeClass('is-invalid');
    } else {
        $('#passwordMatch').html('Not Matching').css('color', 'red');
    }
});

// Toggle if password is visible or not
$('#showPassword').on('click', function() {
    var x = document.getElementById("newPassword");
    if (x.type === "password") {
    x.type = "text";
    } else {
    x.type = "password";
    }
});

// Save changed user settings
$(document).on('click', "#submit_preferences", function(e) {

    // Prevent page refresh
    e.preventDefault();

    // ensure password matches
    if ( $('#newPassword').val() !== $('#repeatPassword').val()) {
        $('#newPassword').addClass('is-invalid');
        $('#newPassword').after('<div class="invalid-tooltip" id="passwordInvalidTooltip">Passwords must match.</div>')
        //console.log("passwords did not match");
        return;
    }

    // disable the submit button and put the UID into the form
    $("#submit_preferences").text("Processing");
    $("#submit_preferences").attr("disabled", true);
    $('#help_id').val(CURRENT_USER.help_id);

    var formData = $("#settings_form").serializeArray();
    formData.push({
        'name':'scope',
        'value':'settings_change'
    },
    {
        'name':'user_id',
        'value':CURRENT_USER.id
    });

    $.ajax({
        url: './cgi/save_user_account_changes.cgi',
        type: "POST",
        data : formData,
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            if (data['success'] == 1) {
                var msg = "Settings have been saved!";
                alert_user_success(msg);
            } else {
                var msg = "An error occurred while updating preferences.  Please try again.";
                alert_user_error(msg);
	          }
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('jqXHR: ', jqXHR);
            console.log('textStatus: ', textStatus);
            console.log('errorThrown: ', errorThrown);

            var msg = "Unable to update preferences. Please try again.";
            alert_user_error(msg);
        }
    }); //end ajax

    //reset submit button area no matter the outcome
    $("#submit_preferences").text("Update Profile");
    $("#submit_preferences").attr("disabled", false);
});

function alert_user_error(message){
    $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message"><strong>Oops! </strong> ' + message + '</p></div>').show();
};

function alert_user_success(message){
    $('.alert-container').html('<div class="alert alert-success alert-dismissible" role="alert">' +
      '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      '<p class="alert-message">' + message + '</p></div>').show();
};
