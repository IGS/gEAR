window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Create an account';

    document.addEventListener('click', function(e) {

        if (e.target.id === 'btn-account-creation-submit') {
            e.preventDefault();
            var formData = new FormData(document.getElementById('account_creation'));

            // Validate form's completion. Exit if it contains errors and alert user
            if (validate_account_creation_form(formData) == false){
                return false;
            }

            console.log("Would have submitted the account creation form");
            return false;

            var xhr = new XMLHttpRequest();
            xhr.open('POST', './cgi/create_account.cgi', true);
            xhr.onload = function() {
                if (xhr.status === 200) {
                    var data = JSON.parse(xhr.responseText);
                    // -1 means the account email existed already
                    if (data['session_id'] == -1) {
                        document.getElementById('email_already_exists').style.display = 'block';
                    } else {
                        CURRENT_USER.email = document.getElementById('inputEmail').value;
                        CURRENT_USER.session_id = data['session_id'];
                        document.querySelector('span.user_logged_in').textContent = CURRENT_USER.user_name;

                        document.getElementById('login_controls').style.display = 'none';
                        document.getElementById('loggedin_controls').style.display = 'block';

                        // https://github.com/js-cookie/js-cookie
                        Cookies.set('gear_session_id', CURRENT_USER.session_id, { expires: 7 });

                        // now redirect to the home page
                        window.location.href = './index.html';
                    }
                } else {
                    document.querySelector('.alert-container').innerHTML = '<div class="alert alert-danger alert-dismissible" role="alert">' +
                                                                           '<button type="button" class="close close-alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                                                                           '<p class="alert-message"><strong>Oops! </strong>Something went wrong.</p>' +
                                                                           '<p>Please e-mail jorvis@gmail.com for help creating your account.</p></div>';
                    document.querySelector('.alert-container').style.display = 'block';
                }
            };
            xhr.send(formData);
        }
    });

    document.getElementById('first-last').addEventListener('blur', function() {
        validateFirstLast();
    });

    document.getElementById('email').addEventListener('blur', function() {
        validateEmail();
    });

    document.getElementById('password1').addEventListener('keyup', function() {
        validatePassword('typing');
    });

    document.getElementById('password2').addEventListener('blur', function() {
        validatePassword2();
    });
};

const handlePageSpecificLoginUIUpdates = async (event) => {
    // Nothing to do here at the moment
}

function validateEmail() {
    const email = document.getElementById('email').value;
    const email_regex = /^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$/;
    if (! email_regex.test(email) ) {
        // add the is-danger class to the input
        document.getElementById('email').classList.add('is-danger');
        document.getElementById('email-error-message').classList.remove('is-hidden');
        document.getElementById('email-alert-icon').classList.remove('is-hidden');
        return false;
    } else {
        document.getElementById('email').classList.remove('is-danger');
        document.getElementById('email-error-message').classList.add('is-hidden');
        document.getElementById('email-alert-icon').classList.add('is-hidden');
        return true;
    }

    // TODO: Also check if this e-mail is already registered
}

function validateFirstLast() {
    const first_last = document.getElementById('first-last').value;

    if ( first_last.length > 2 && first_last.includes(' ')) {
        document.getElementById('first-last').classList.remove('is-danger');
        document.getElementById('first-last-error-message').classList.add('is-hidden');
        document.getElementById('first-last-alert-icon').classList.add('is-hidden');
        return true;
    } else {
        document.getElementById('first-last').classList.add('is-danger');
        document.getElementById('first-last-error-message').classList.remove('is-hidden');
        document.getElementById('first-last-alert-icon').classList.remove('is-hidden');
        return false;
    }
}

function validatePassword(mode) {
    const password1 = document.getElementById('password1').value;
    const password2 = document.getElementById('password2').value;

    if (mode == 'typing') {
        // Check the list of password requirements
        if (password1.length < 8) {
            validatePasswordToggleRequirement('character-limit', 'fail');
        } else {
            validatePasswordToggleRequirement('character-limit', 'pass');
        }

        if (password1.match(/[A-Z]/) == null) {
            validatePasswordToggleRequirement('upper-char', 'fail');
        } else {
            validatePasswordToggleRequirement('upper-char', 'pass');
        }

        if (password1.match(/[a-z]/) == null) {
            validatePasswordToggleRequirement('lower-char', 'fail');
        } else {
            validatePasswordToggleRequirement('lower-char', 'pass');
        }

        if (password1.match(/[0-9]/) == null) {
            validatePasswordToggleRequirement('number', 'fail');
        } else {
            validatePasswordToggleRequirement('number', 'pass');
        }

        if (password1.match(/[^A-Za-z0-9]/) == null) {
            validatePasswordToggleRequirement('special-char', 'fail');
        } else {
            validatePasswordToggleRequirement('special-char', 'pass');
        }

    } else {
        if (password1 != password2) {
            document.getElementById('password1').classList.add('is-danger');
            document.getElementById('password2').classList.add('is-danger');
            document.getElementById('password1-error-message').innerHTML = 'Passwords do not match';
            document.getElementById('password2-error-message').innerHTML = 'Passwords do not match';
            document.getElementById('password1-error-message').classList.remove('is-hidden');
            document.getElementById('password2-error-message').classList.remove('is-hidden');
            document.getElementById('password1-alert-icon').classList.remove('is-hidden');
            document.getElementById('password2-alert-icon').classList.remove('is-hidden');
            return false;
        } else {
            document.getElementById('password1').classList.remove('is-danger');
            document.getElementById('password2').classList.remove('is-danger');
            document.getElementById('password1-error-message').classList.add('is-hidden');
            document.getElementById('password2-error-message').classList.add('is-hidden');
            document.getElementById('password1-alert-icon').classList.add('is-hidden');
            document.getElementById('password2-alert-icon').classList.add('is-hidden');
            return true;
        }
    }

    return true;
}

function validatePassword2() {
    // Here we only care about if the two passwords match, and we don't want to be 
    //  annoying about it while the user is still typing
    if (document.getElementById('password1').value != document.getElementById('password2').value) {
        document.getElementById('password1').classList.add('is-danger');
        document.getElementById('password2').classList.add('is-danger');
        document.getElementById('password1-error-message').innerHTML = 'Passwords do not match';
        document.getElementById('password2-error-message').innerHTML = 'Passwords do not match';
        document.getElementById('password1-error-message').classList.remove('is-hidden');
        document.getElementById('password2-error-message').classList.remove('is-hidden');
        document.getElementById('password1-alert-icon').classList.remove('is-hidden');
        document.getElementById('password2-alert-icon').classList.remove('is-hidden');
        return false;
    } else {
        document.getElementById('password1').classList.remove('is-danger');
        document.getElementById('password2').classList.remove('is-danger');
        document.getElementById('password1-error-message').classList.add('is-hidden');
        document.getElementById('password2-error-message').classList.add('is-hidden');
        document.getElementById('password1-alert-icon').classList.add('is-hidden');
        document.getElementById('password2-alert-icon').classList.add('is-hidden');
    }
}

function validatePasswordToggleRequirement(requirement, state) {
    // state is either 'pass' or 'fail'   
    let selectorString = '#pc-' + requirement + ' i';

    if (state == 'pass') {
        document.querySelector(selectorString).classList.remove('mdi-emoticon-sad-outline');
        document.querySelector(selectorString).classList.add('mdi-check-bold');
    } else {
        document.querySelector(selectorString).classList.remove('mdi-check-bold');
        document.querySelector(selectorString).classList.add('mdi-emoticon-sad-outline');
    }
}

function validate_account_creation_form(formData) {
    if (validateFirstLast() == false) {
        return false;
    }

    if (validateEmail() == false) {
        return false;
    }

    if (validatePassword('submit') == false) {
        return false;
    }

    // if we made it this far, things are good
    return true;
}
