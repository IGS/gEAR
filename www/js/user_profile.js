'use strict';

import { convertToFormData, createToast, getCurrentUser, initCommonUI, logErrorInConsole, registerPageSpecificLoginUIUpdates } from "./common.v2.js?v=670b2ed";

let passwordUpdateWanted = false;

// When password and repeated password are not the same, add a tooltip
for (const classElt of document.getElementsByClassName("js-password")) {
    classElt.addEventListener("keyup", () => {
        const newPasswordElt = document.getElementById("new-password");
        const repeatPasswordElt = document.getElementById("repeat-password");

        if (newPasswordElt.value !== repeatPasswordElt.value) {
            repeatPasswordElt.classList.add("is-danger");
            repeatPasswordElt.classList.remove("is-success");
            document.getElementById("password-match").classList.add("is-hidden");
            document.getElementById("password-no-match").classList.remove("is-hidden");
            // disable submit button
            document.getElementById("submit-preferences").disabled = true;
            return;
        }
        repeatPasswordElt.classList.add("is-success");
        repeatPasswordElt.classList.remove("is-danger");
        document.getElementById("password-match").classList.remove("is-hidden");
        document.getElementById("password-no-match").classList.add("is-hidden");
        // enable submit button
        document.getElementById("submit-preferences").disabled = false;
    });
}

// Toggle if password is visible or not
document.getElementById("show-password").addEventListener("click", () => {
    const newPasswordElt = document.getElementById("new-password");
    newPasswordElt.type = newPasswordElt.type === "password" ? "text" : "password";
});

document.getElementById("update-password-toggle").addEventListener("click", () => {
    document.getElementById("update-password-toggle").classList.add("is-hidden");
    document.getElementById("password-update-section").classList.remove("is-hidden");
    document.getElementById("new-password").focus();
    passwordUpdateWanted = true;
});

// submit the form
document.getElementById("submit-preferences").addEventListener("click", async (event) => {
    // Prevent page refresh
    event.preventDefault();

    event.target.classList.add("is-loading");

    const formData ={
        scope: "settings_change",
        session_id: getCurrentUser()?.session_id,
        help_id: getCurrentUser()?.help_id,
        email: document.getElementById("email").value,
        institution: document.getElementById("institution").value,
        colorblind_mode: document.getElementById("colorblind-mode").checked,
        want_updates: document.getElementById("want-updates").checked,
        new_password: document.getElementById("new-password").value,
    };

    try {
        const {data} = await axios.post('./cgi/save_user_account_changes.cgi', convertToFormData(formData));
        if (data?.success === 1) {
            var msg = "Settings have been saved!";
            createToast(msg, "is-success");
        } else {
            var msg = "An error occurred while updating preferences.  Please try again.";
            throw new Error(msg);
        }
    } catch (error) {
        const msg = "Unable to update preferences. Please try again.";
        logErrorInConsole(msg);
        createToast(msg);
    } finally {
        event.target.classList.remove("is-loading");
    }

});

/* --- Entry point --- */
/**
 * Handles page-specific login UI updates.
 *
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are complete.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// User settings has no "active" state for the sidebar
	document.getElementById("page-header-label").textContent = "User Profile";
	for (const elt of document.querySelectorAll("#primary-nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}

    const sessionId = getCurrentUser()?.session_id;

	if (! sessionId ) {
        document.getElementById("not-logged-in-msg").classList.remove("is-hidden");
        return;
    }

    // Get the user's settings from the server
    document.getElementById("email").value = getCurrentUser()?.email;
    document.getElementById("institution").value = getCurrentUser()?.institution;
    document.getElementById("colorblind-mode").checked = getCurrentUser()?.colorblind_mode;
    document.getElementById("want-updates").checked = getCurrentUser()?.updates_wanted;

    // Forcibly empty the password fields so the browser's don't autofill.
    //  Then wait and do it again.
    //  If only they would respect the autocomplete="off" attribute...
    document.getElementById("new-password").value = "";
    document.getElementById("repeat-password").value = "";

    setTimeout(() => {
        document.getElementById("new-password").value = "";
        document.getElementById("repeat-password").value = "";
    }, 500);
};
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

// Pre-initialize some stuff
await initCommonUI();