'use strict';

import { convertToFormData, initCommonUI, registerPageSpecificLoginUIUpdates } from "./common.v2.js?v=056890d";
// Pre-initialize some stuff
initCommonUI();

let help_id = null;

window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Lost password recovery';

    document.addEventListener('click', async function(e) {

        if (e.target.id === 'btn-password-update-submit') {
            e.preventDefault();

            // disable the button so the user doesn't click it again
            document.getElementById('btn-password-update-submit').classList.add('is-loading');

            // Validate form's completion. Exit if it contains errors and alert user
            const form_data_valid = await validatePasswordUpdateForm();
            if (form_data_valid == false){
                document.getElementById('btn-password-update-submit').classList.remove('is-loading');
                return false;
            } else {
                // Submit the form
                const {data} = await axios.post('./cgi/update_password.cgi', convertToFormData({
                    'password': document.getElementById('password1').value,
                    'help_id': help_id
                }));

                document.getElementById('btn-password-update-submit').classList.remove('is-loading');

                if (data["success"]) {
                    document.getElementById('reset-form').classList.add('is-hidden');
                    document.getElementById('reset-success').classList.remove('is-hidden');
                } else {
                    document.getElementById('password-update-status-message').innerHTML = 'There was a problem updating your password. ' + data['error'] + '. Please try again later or contact us.';
                    document.getElementById('password-update-status-message').classList.remove('is-hidden');
                }

            }
        }
    });

    document.getElementById('initial-submit').addEventListener('click', function() {
        document.getElementById('initial-submit').classList.add('is-loading');
        document.getElementById('initial-submit').disabled = true;

        // Did they enter a valid email address?
        if (validateEmail(document.getElementById('email').value)) {
            sendVerificationEmail();
        }
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
    help_id = getUrlParameter('help_id');

    if (help_id) {
        document.getElementById('reset-form').classList.remove('is-hidden');
    } else {
        document.getElementById('initial-form').classList.remove('is-hidden');
    }
}
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

async function sendVerificationEmail(verification_uuid) {
    console.debug("Sending verification email");

    const current_url = window.location.href;
    const current_page = current_url.lastIndexOf("/");
    const destination_page = `${current_url.substring(0, current_page)}/forgot_password.html`;

    const {data} = await axios.post('./cgi/send_email.cgi', convertToFormData({
        'email': document.getElementById('email').value,
        'scope': 'forgot_password',
        'destination_page': destination_page
    }));

    document.getElementById('initial-submit').classList.remove('is-loading');

    if (data["success"]) {
        document.getElementById('email-success-message').classList.remove('is-hidden');
    } else {
        document.getElementById('email-error-message').innerHTML = 'There was a problem sending the e-mail. Please try again later.';
        document.getElementById('email').classList.add('is-danger');
        document.getElementById('email-error-message').classList.remove('is-hidden');
        document.getElementById('email-alert-icon').classList.remove('is-hidden');
        document.getElementById('initial-submit').disabled = false;
    }

    return Boolean(data["success"]);
}

const validateEmail = (emailAddress) =>  {
    // TODO: Make this part of common.js along with the create_account.js methods
    const email = emailAddress.trim();
    const email_regex = /^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$/;
    if (email_regex.test(email) ) {
        document.getElementById('email').classList.remove('is-danger');
        document.getElementById('email-error-message').classList.add('is-hidden');
        document.getElementById('email-alert-icon').classList.add('is-hidden');
        return true;
    } else {
        document.getElementById('email-error-message').innerHTML = 'Please enter a valid e-mail address';
        document.getElementById('email').classList.add('is-danger');
        document.getElementById('email-error-message').classList.remove('is-hidden');
        document.getElementById('email-alert-icon').classList.remove('is-hidden');
        return false;
    }
}

function validatePassword(mode) {
    /*
    This always checks password complexity requirements, but if mode is
    'submit' it will also check if the two passwords match.
    */

    const password1 = document.getElementById('password1').value;
    const password2 = document.getElementById('password2').value;

    let complexity_met = true;

    // Check the list of password requirements
    if (password1.length < 8) {
        validatePasswordToggleRequirement('character-limit', 'fail');
        complexity_met = false;
    } else {
        validatePasswordToggleRequirement('character-limit', 'pass');
    }

    if (password1.match(/[A-Z]/) == null) {
        validatePasswordToggleRequirement('upper-char', 'fail');
        complexity_met = false;
    } else {
        validatePasswordToggleRequirement('upper-char', 'pass');
    }

    if (password1.match(/[a-z]/) == null) {
        validatePasswordToggleRequirement('lower-char', 'fail');
        complexity_met = false;
    } else {
        validatePasswordToggleRequirement('lower-char', 'pass');
    }

    if (password1.match(/[0-9]/) == null) {
        validatePasswordToggleRequirement('number', 'fail');
        complexity_met = false;
    } else {
        validatePasswordToggleRequirement('number', 'pass');
    }

    if (password1.match(/[^A-Za-z0-9]/) == null) {
        validatePasswordToggleRequirement('special-char', 'fail');
        complexity_met = false;
    } else {
        validatePasswordToggleRequirement('special-char', 'pass');
    }

    if (complexity_met == false) {
        document.getElementById('password1').classList.add('is-danger');
        document.getElementById('password1-error-message').innerHTML = 'Password does not meet complexity requirements';
        document.getElementById('password1-error-message').classList.remove('is-hidden');
        document.getElementById('password1-alert-icon').classList.remove('is-hidden');
        return false;
    }

    document.getElementById('password1').classList.remove('is-danger');
    document.getElementById('password1-error-message').classList.add('is-hidden');
    document.getElementById('password1-alert-icon').classList.add('is-hidden');

    if (mode != 'submit') {
        return;
    }

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
    }
    document.getElementById('password1').classList.remove('is-danger');
    document.getElementById('password2').classList.remove('is-danger');
    document.getElementById('password1-error-message').classList.add('is-hidden');
    document.getElementById('password2-error-message').classList.add('is-hidden');
    document.getElementById('password1-alert-icon').classList.add('is-hidden');
    document.getElementById('password2-alert-icon').classList.add('is-hidden');
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
    const selectorString = `#pc-${requirement} i`;

    if (state == 'pass') {
        document.querySelector(selectorString).classList.remove('mdi-emoticon-sad-outline');
        document.querySelector(selectorString).classList.add('mdi-check-bold');
    } else {
        document.querySelector(selectorString).classList.remove('mdi-check-bold');
        document.querySelector(selectorString).classList.add('mdi-emoticon-sad-outline');
    }
}

async function validatePasswordUpdateForm() {
    if (!validatePassword('submit')) {
        return false;
    }

    // if we made it this far, things are good
    return true;
}