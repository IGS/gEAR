'use strict';

import { convertToFormData, initCommonUI } from "./common.v2.js?v=f0a5adc";

let verification_uuid = null;

window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Create an account';

    document.addEventListener('click', async function(e) {

        if (e.target.id === 'btn-account-creation-submit') {
            e.preventDefault();

            // disable the button so the user doesn't click it again
            document.getElementById('btn-account-creation-submit').classList.add('is-loading');

            // Validate form's completion. Exit if it contains errors and alert user
            const form_data_valid = await validateAccountCreationForm();
            if (form_data_valid == false){
                document.getElementById('btn-account-creation-submit').classList.remove('is-loading');
                return false;
            }

            // generate a UUID for the user
            verification_uuid = crypto.randomUUID();
            const email_sent = await sendVerificationEmail(verification_uuid);

            if (email_sent == false) {
                // TODO: Handle UI display here.
                alert("There was an error sending the verification email. Please try again later.");
                document.getElementById('btn-account-creation-submit').classList.remove('is-loading');
                return false;
            }

            // hide the account info and show the verification info
            document.getElementById('account-info-c').classList.add('is-hidden');
            document.getElementById('email-verification-c').classList.remove('is-hidden');

            return false;

        } else if (e.target.id === 'btn-email-verification-submit') {
            e.preventDefault();

            // This returns false no matter what. await/sync (or developer) issue, so UI code shifted within
            const account_created = await createAccount(verification_uuid);
            return false;
        }

    });

    document.getElementById('first-last').addEventListener('blur', function() {
        validateFirstLast();
    });

    document.getElementById('email').addEventListener('blur', async function() {
        await validateEmail();
    });

    document.getElementById('password1').addEventListener('keyup', function() {
        validatePassword('typing');
    });

    document.getElementById('password2').addEventListener('blur', function() {
        validatePassword2();
    });
};

async function createAccount(verification_uuid) {
    // disable the button so it's not clicked again
    document.getElementById('btn-email-verification-submit').classList.add('is-loading');

    // get the value of the colorblind mode checkbox, if it's checked
    const colorblind_mode = document.getElementById('colorblind-mode').checked ? 'yes' : 0;

    // get the value of the email_updates checkbox, if it's checked
    const email_updates = document.getElementById('email-updates').checked ? 'yes' : 0;

    const {data} = await axios.post('./cgi/create_account.cgi', convertToFormData({
        'first-last': document.getElementById('first-last').value,
        'institution': document.getElementById('institution').value,
        'email': document.getElementById('email').value,
        'password': document.getElementById('password1').value,
        'verification_code_long': verification_uuid,
        'verification_code_short': document.getElementById('verification-code').value,
        'colorblind_mode': colorblind_mode,
        'email_updates': email_updates
    }));

    if (data['success']) {
        document.getElementById('account-info-c').classList.add('is-hidden');
        document.getElementById('email-verification-c').classList.add('is-hidden');
        document.getElementById('account-creation-success-c').classList.remove('is-hidden');

        // populate the login form, then log the user in
        document.getElementById('user-email').value = document.getElementById('email').value;
        document.getElementById('user-password').value = document.getElementById('password1').value;
        doLogin(false);
    } else {
        document.getElementById('account-info-c').classList.add('is-hidden');
        document.getElementById('email-verification-c').classList.add('is-hidden');
        document.getElementById('account-creation-failure-c').classList.remove('is-hidden');
        console.error('error: ' + data['error'])
        document.getElementById('btn-email-verification-submit').classList.remove('is-loading');
        return false;
    }

    return Boolean(data["success"]);
}

async function sendVerificationEmail(verification_uuid) {
    const {data} = await axios.post('./cgi/send_email.cgi', convertToFormData({
        'email': document.getElementById('email').value,
        'scope': 'user_verification',
        'verification_code_long': verification_uuid,
    }));

    return Boolean(data["success"]);
}

async function validateEmail() {
    const email = document.getElementById('email').value;
    const email_regex = /^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$/;
    if (email_regex.test(email) ) {
        document.getElementById('email').classList.remove('is-danger');
        document.getElementById('email-error-message').classList.add('is-hidden');
        document.getElementById('email-alert-icon').classList.add('is-hidden');
    } else {
        document.getElementById('email-error-message').innerHTML = 'Please enter a valid e-mail address';
        document.getElementById('email').classList.add('is-danger');
        document.getElementById('email-error-message').classList.remove('is-hidden');
        document.getElementById('email-alert-icon').classList.remove('is-hidden');
        return false;
    }

    // Also check if this e-mail is already registered
    const {data} = await axios.post('./cgi/check_existing_email.cgi', convertToFormData({
        'email': email,
    }));

    if (data['email_exists'] == 1) {
        document.getElementById('email-error-message').innerHTML = 'This e-mail address is already registered';
        document.getElementById('email').classList.add('is-danger');
        document.getElementById('email-error-message').classList.remove('is-hidden');
        return false;
    }

    return true;
}

function validateFirstLast() {
    const first_last = document.getElementById('first-last').value;

    if ( first_last.length > 2 && first_last.includes(' ')) {
        document.getElementById('first-last').classList.remove('is-danger');
        document.getElementById('first-last-error-message').classList.add('is-hidden');
        document.getElementById('first-last-alert-icon').classList.add('is-hidden');
        return true;
    }
    document.getElementById('first-last').classList.add('is-danger');
    document.getElementById('first-last-error-message').classList.remove('is-hidden');
    document.getElementById('first-last-alert-icon').classList.remove('is-hidden');
    return false;
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

async function validateAccountCreationForm() {
    // function returns bool
    if (!validateFirstLast()) {
        return false;
    }

    // function returns bool
    if (!(await validateEmail())) {
        return false;
    }

    if (!validatePassword('submit')) {
        return false;
    }

    // if we made it this far, things are good
    return true;
}

// Pre-initialize some stuff
await initCommonUI();
