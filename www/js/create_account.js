window.onload=function() {

    $('#btn_account_creation_cancel').click(function(e) {
        e.preventDefault();

        // redirect to home page
        window.location.href = "./index.html";
    });

    $(document).on('click', '#btn_account_creation_submit', function(e) {
        e.preventDefault();
        var formData = $("#account_creation").serializeArray();

        // Validate form's completion. Exit if it contains errors and alert user
        if (validate_account_creation_form(formData) == false){
            return false;
        }

        $.ajax({
            url : './cgi/create_account.cgi',
            type: "POST",
            data : formData,
            dataType:"json",
            success: function(data, textStatus, jqXHR) {
                // -1 means the account email existed already
                if (data['session_id'] == -1) {
                    $("#email_already_exists").show();

                    // ? - Any other value is the session ID
                } else {
                    CURRENT_USER.email = $('#inputEmail').val();
                    CURRENT_USER.session_id = data['session_id'];
                    $('span.user_logged_in').text(CURRENT_USER.user_name);

                    // $('#login_controls').hide();
                    $('#login_controls').attr("style", "display: none !important");

                    $('#loggedin_controls').show();
                    // https://github.com/js-cookie/js-cookie
                    Cookies.set('gear_session_id', CURRENT_USER.session_id, { expires: 7 });

                    // now redirect to the home page
                    window.location.href = "./index.html";
                }
            },
            error: function (jqXHR, textStatus, errorThrown) {
      			$('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
      				                       '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
      				                       '<p class="alert-message"><strong>Oops! </strong>Something went wrong.</p>' +
      				                       '<p>Please e-mail jorvis@gmail.com for help creating your account.</p></div>').show();
            }
        });
    });
};

function validate_account_creation_form(formData) {
    $(".account-creation-error-message").hide();
    var is_form_valid = true;

    var user_name = "";
    var email = "";
    var password1 = "";
    var password2 = "";
    $.each(formData, function(i, item){
        if (item.name == 'inputName') {
            user_name = item.value;
        }
        if (item.name == 'inputEmail') {
            email = item.value;
        }
        if (item.name == 'inputPassword') {
            password1 = item.value;
        }
        if (item.name == 'retypePassword') {
            password2 = item.value;
        }
    });

    //User name: Check length
    if ( user_name.length < 2 ) {
        $("#name_invalid").show();
        is_form_valid = false;
    }

    //Email: Check format
    var email_regex = /^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$/;
    if ( !email_regex.test(email) ) {
      $("#email_invalid").show();
      is_form_valid = false;
    }

    //Password: 1) Check length 2) Does the retype match
    if ( password1.length > 2  ) {
        if ( password1 != password2 ) {
          $("#password_mismatch").show();
          is_form_valid = false;
        }
    } else {
        $("#password_invalid").show();
        is_form_valid = false;
    }

    return is_form_valid;
}
