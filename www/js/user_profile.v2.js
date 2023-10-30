'use strict';

/* Creates a Toast-style message in the upper-corner of the screen. */
const createToast = (msg, levelClass="is-danger") => {
    const template = `
    <div class="notification js-toast ${levelClass} animate__animated animate__fadeInUp animate__faster">
        <button class="delete"></button>
        ${msg}
    </div>
    `
    const html = generateElements(template);

    const numToasts = document.querySelectorAll(".js-toast.notification").length;

    if (document.querySelector(".js-toast.notification")) {
        // If .js-toast notifications are present, append under final notification
        // This is to prevent overlapping toast notifications
        document.querySelector(".js-toast.notification:last-of-type").insertAdjacentElement("afterend", html);
        // Position new toast under previous toast with CSS
        html.style.setProperty("top", `${(numToasts * 70) + 30}px`);
    } else {
        // Otherwise prepend to top of main content
        document.getElementById("main_c").prepend(html);
    }

    // This should get the newly added notification since it is now the first
    document.querySelector(".js-toast.notification .delete").addEventListener("click", (event) => {
        const notification = event.target.closest(".js-toast.notification");
        notification.remove(notification);
    });

    // For a success message, remove it after 3 seconds
    if (levelClass === "is-success") {
        const notification = document.querySelector(".js-toast.notification:last-of-type");
        notification.classList.remove("animate__fadeInUp");
        notification.classList.remove("animate__faster");
        notification.classList.add("animate__fadeOutDown");
        notification.classList.add("animate__slower");
    }
}

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// User settings has no "active" state for the sidebar
	document.querySelector("#header_bar .navbar-item").textContent = "User Profile";
	for (const elt of document.querySelectorAll("#primary_nav > aside > ul > li > a")) {
		elt.classList.remove("is-active");
	}

    const sessionId = CURRENT_USER.session_id;

	if (! sessionId ) {
        document.getElementById("not_logged_in_msg").classList.remove("is-hidden");
        return;
    }

    // Get the user's settings from the server
    document.getElementById("email").value = CURRENT_USER.email;
    document.getElementById("institution").value = CURRENT_USER.institution;
    document.getElementById("colorblind_mode").checked = CURRENT_USER.colorblind_mode;
    document.getElementById("want_updates").checked = CURRENT_USER.updates_wanted;

};

// When password and repeated password are not the same, add a tooltip
for (const classElt of document.getElementsByClassName("js-password")) {
    classElt.addEventListener("keyup", () => {
        const newPasswordElt = document.getElementById("new_password");
        const repeatPasswordElt = document.getElementById("repeat_password");

        if (newPasswordElt.value !== repeatPasswordElt.value) {
            repeatPasswordElt.classList.add("is-danger");
            repeatPasswordElt.classList.remove("is-success");
            document.getElementById("password_match").classList.add("is-hidden");
            document.getElementById("password_no_match").classList.remove("is-hidden");
            // disable submit button
            document.getElementById("submit_preferences").disabled = true;
            return;
        }
        repeatPasswordElt.classList.add("is-success");
        repeatPasswordElt.classList.remove("is-danger");
        document.getElementById("password_match").classList.remove("is-hidden");
        document.getElementById("password_no_match").classList.add("is-hidden");
        // enable submit button
        document.getElementById("submit_preferences").disabled = false;
    });
}

// Toggle if password is visible or not
document.getElementById("show_password").addEventListener("click", () => {
    const newPasswordElt = document.getElementById("new_password");
    newPasswordElt.type = newPasswordElt.type === "password" ? "text" : "password";
});

// submit the form
document.getElementById("submit_preferences").addEventListener("click", async (event) => {
    // Prevent page refresh
    event.preventDefault();

    event.target.classList.add("is-loading");

    const formData ={
        scope: "settings_change",
        session_id: CURRENT_USER.session_id,
        help_id: CURRENT_USER.help_id,
        email: document.getElementById("email").value,
        institution: document.getElementById("institution").value,
        colorblind_mode: document.getElementById("colorblind_mode").checked,
        want_updates: document.getElementById("want_updates").checked,
        new_password: document.getElementById("new_password").value,
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