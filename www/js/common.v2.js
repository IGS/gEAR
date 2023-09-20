document.addEventListener('DOMContentLoaded', () => {
    // Add listeners for any elements which have generic close controls
    // NOTE: ".delete" has Bulma CSS associated, which overwrites mdi icon classes
    (document.querySelectorAll('.notification .delete') || []).forEach(($delete) => {
        // const $notification = $delete.parentNode;
        // SAdkins - modifying to find ".notification" in case .delete is not immediate child. Otherwise partial removal could happen.
        const $notification = $delete.closest(".notification");

        $delete.addEventListener('click', () => {
            $notification.parentNode.removeChild($notification);
        });
    });
});

const checkForLogin = async () => {
    // success:  returns session_id
    // failure:  returns null
    const session_id = Cookies.get('gear_session_id');
    CURRENT_USER = new User();

    if (! session_id || session_id === "undefined") {
        // no cookie found, so user is not logged in
        handleLoginUIUpdates();
    } else {
        const payload = { session_id }
        // we found a cookie, so see if it's a valid session
        try {
            const {data} = await axios.post("/cgi/get_session_info.cgi", convertToFormData(payload));

            if (data.success) {
                CURRENT_USER = new User({session_id, ...data});
                CURRENT_USER.setDefaultProfile();
                document.querySelector('#navbarBasicExample > div.navbar-end > div > div.navbar-item.has-dropdown.is-hoverable > a').textContent = CURRENT_USER.user_name;
                handleLoginUIUpdates();
            } else {
                // session_id is invalid, so remove cookie
                Cookies.remove('gear_session_id');
                throw new Error(`Invalid session_id: ${session_id}`);
            }
        } catch (error) {
            console.error(error);
        }
    }

    // Now that session_id has been obtained, we can trigger events that depend on it
    trigger(document, "build_jstrees");

    return session_id;
}


/* Generate a DocumentFragment based on an HTML template. Returns htmlCollection that can be appended to a parent HTML */
const generateElements = (html) => {
    const template = document.createElement('template');
    template.innerHTML = html.trim();
    return template.content.children[0];
}

const handleLoginUIUpdates = () => {
    // pass
}

// Equivalent to jQuery "trigger" (https://youmightnotneedjquery.com/#trigger_native)
const trigger = (el, eventType) => {
    if (typeof eventType === 'string' && typeof el[eventType] === 'function') {
        el[eventType]();
    } else {
        const event =
        typeof eventType === 'string'
            ? new Event(eventType, {bubbles: true})
            : eventType;
        el.dispatchEvent(event);
    }
}

// Convert POST payload to FormData so that POSTing to CGI scripts that use cgi.FieldStorage will work
const convertToFormData = (object) => {
    // Source -> https://stackoverflow.com/a/66611630
    // NOTE: When using FormData do not set headers to application/json
    const formData = new FormData();
    for (const key of Object.keys(object)) {
        if (typeof object[key] !== 'object') formData.append(key, object[key])
        else formData.append(key, JSON.stringify(object[key]))
    };
    return formData;
}

// If "step" sections are clicked, expand that section and collaps the others
const jsSteps = document.getElementsByClassName("js-step");
const resetSteps = (event) => {
    // If clicked area has no parentNode (i.e. clicked target is destroyed and recreated), just leave be
    if (! event.target.parentNode) return;

    // If clicked on existing section, leave it
    const currentStep = event.target.closest(".js-step");
    if (currentStep.classList.contains("step-active")) return;

    // Make every step uniform
    for (const jsStep of jsSteps) {
        // Reset active step
        jsStep.classList.remove("step-active");

        // Reset the bg color
        jsStep.classList.remove("has-background-light");
        jsStep.classList.add("has-background-white-bis");

        // Reset the rightmost arrow
        const rightIcon = jsStep.querySelector("span.is-pulled-right i");
        rightIcon.classList.remove("mdi-chevron-up");
        rightIcon.classList.add("mdi-chevron-down");

        // Collapse all "non-header" parts of step (in CSS file)
    }

    // Modify this step to look active
    currentStep.classList.remove("has-background-white-bis");
    currentStep.classList.add("has-background-light");
    currentStep.querySelector("span.is-pulled-right i").classList.remove("mdi-chevron-down");
    currentStep.querySelector("span.is-pulled-right i").classList.add("mdi-chevron-up");
    // Adding "step-active" expands all "non-header" parts of step (in CSS file)
    currentStep.classList.add("step-active")

}

for (const jsStep of jsSteps) {
    // Add "capture=true" to trigger parent before children events
    // This helps with manually triggering click events to open up the next section.
    jsStep.addEventListener("click", resetSteps, {capture: true});
}