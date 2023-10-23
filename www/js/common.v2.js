let CURRENT_USER;

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

    // modal code from https://bulma.io/documentation/components/modal/

    // Add a click event on various child elements to close the parent modal
    (document.querySelectorAll('.modal-background, .modal-close, .modal-card-head .delete, .modal-card-foot .button') || []).forEach(($close) => {
        const $target = $close.closest('.modal');

        $close.addEventListener('click', () => {
            closeModal($target);
        });
    });

    // Add a keyboard event to close all modals
    document.addEventListener('keydown', (event) => {
        if (event.code === 'Escape') {
            closeAllModals();
        }
    });

    // Add click event for collapsable left navbar text
    const navbar_toggler = document.querySelector('#navbar-toggler');
    navbar_toggler.addEventListener('click', (event) => {

        // toggle is-collapsed on all span elements of class icon-text-part
        (document.querySelectorAll('.icon-text-part') || []).forEach(($menuLabel) => {
            $menuLabel.classList.toggle('is-collapsed');
        });

        // hide the citation element
        document.querySelector("#citation_c").classList.add("is-hidden")
    });



    checkForLogin();
});

// Functions to open and close a modal
function openModal($el) {
    $el.classList.add('is-active');
}

function closeModal($el) {
    $el.classList.remove('is-active');
}

function closeAllModals() {
    (document.querySelectorAll('.modal') || []).forEach(($modal) => {
        closeModal($modal);
    });
}

/*************************************************************************************
 Code related to the login process, which is available in the header across all pages.
*************************************************************************************/

document.getElementById('submit-login').addEventListener('click', function(event) {
    event.preventDefault(); // Prevent the default form submission

    // reset any UI elements
    document.getElementById('user-email-help').classList.add('is-hidden');
    document.getElementById('user-password-help').classList.add('is-hidden');

    // try to log in
    doLogin();
});

document.getElementById('submit-logout').addEventListener('click', function(event) {
    // Clear session information and redirect to home page
    Cookies.remove('gear_session_id');
    CURRENT_USER = undefined;
    window.location.replace('./index.html');
});

const checkForLogin = async () => {
    const session_id = Cookies.get('gear_session_id');
    CURRENT_USER = new User();

    if (! session_id || session_id === "undefined") {
        // no cookie found, so user is not logged in
        handleLoginUIUpdates();
    } else {
        const payload = { session_id }
        // we found a cookie, so see if it's a valid session
        try {
            const {data} = await axios.post("/cgi/get_session_info.v2.cgi", convertToFormData(payload));

            if (data.success) {
                CURRENT_USER = new User({session_id, ...data});
                //CURRENT_USER.setDefaultProfile();
                document.getElementById('current-user-name').textContent = data.user_name;
                handleLoginUIUpdates();

            } else {
                // session_id is invalid, so remove cookie
                Cookies.remove('gear_session_id');
                throw new Error(`Invalid session_id`);
            }
        } catch (error) {
            console.error(error);
        }
    }
}

const doLogin = async () => {
    const formdata = new FormData(document.getElementById("login-form"));
    const payload = new URLSearchParams(formdata);
    const {data} = await axios.post("/cgi/login.v2.cgi", payload);

    if (data.session_id == 0) {
        // user name wasn't found at all
        document.getElementById('user-email-help').classList.remove('is-hidden');

    } else if (data.session_id == -1) {
        // User was found, password was incorrect
        document.getElementById('user-password-help').classList.remove('is-hidden');

    } else if (data.session_id.toString().length > 10) {
        // Looks like we got a session ID
        document.getElementById('current-user-name').textContent = data.user_name;
        hideNotLoggedInElements();
        showLoggedInElements();

        // create the user object
        CURRENT_USER = new User({...data});

        // set the cookie
        Cookies.set('gear_session_id', CURRENT_USER.session_id, { expires: 7 });

    } else {
        // Something went muy wrong

    }
}

const hideLoggedInElements = () => {
    document.querySelectorAll('.logged-in').forEach(element => element.style.display = 'none');
}

const hideNotLoggedInElements = () => {
    document.querySelectorAll('.not-logged-in').forEach(element => element.style.display = 'none');
}

const showLoggedInElements = () => {
    document.querySelectorAll('.logged-in').forEach(element => element.style.display = '');
}

const showNotLoggedInElements = () => {
    document.querySelectorAll('.not-logged-in').forEach(element => element.style.display = '');
}

const handleLoginUIUpdates = () => {
    // So that all elements don't initially show while login is checked, we
    //  show/hide elements first then parent container
    if (CURRENT_USER.session_id) {
        hideNotLoggedInElements();
        showLoggedInElements();
    } else {
        hideLoggedInElements();
        showNotLoggedInElements();
    }
    document.querySelector("#navbar-login-controls").classList.remove("is-hidden");

    trigger(document, handlePageSpecificLoginUIUpdates);
}

/*************************************************************************************
   End of login-related code
*************************************************************************************/

/* Generate a DocumentFragment based on an HTML template. Returns htmlCollection that can be appended to a parent HTML */
const generateElements = (html) => {
    const template = document.createElement('template');
    template.innerHTML = html.trim();
    return template.content.children[0];
}

// Equivalent to jQuery "trigger" (https://youmightnotneedjquery.com/#trigger_native)
const trigger = (el, eventType) => {

    if (typeof eventType === 'string' && typeof el[eventType] === 'function') {
        el[eventType]();
    } else if (typeof eventType === 'function') {
        eventType();
    } else {
        const event =
        typeof eventType === 'string'
            ? new Event(eventType, {bubbles: true})
            : eventType;
        el.dispatchEvent(event);
    }
}

const logErrorInConsole = (error) => {
    if (error.response) {
        // The request was made and the server responded with a status code
        // that falls out of the range of 2xx
        console.error(error.response.data);
        console.error(error.response.status);
        console.error(error.response.headers);
    } else if (error.request) {
        // The request was made but no response was received
        // `error.request` is an instance of XMLHttpRequest in the browser and an instance of
        // http.ClientRequest in node.js
        console.error(error.request);
    } else {
        // Something happened in setting up the request that triggered an Error
        console.error('Error', error.message);
    }
    console.error(error.config);
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

const collapseJsStep = (stepElt) => {
        // Reset active step
        stepElt.classList.remove("step-active");

        // Reset the bg color
        stepElt.classList.remove("has-background-light");
        stepElt.classList.add("has-background-white-bis");

        // Reset the rightmost arrow
        const rightIcon = stepElt.querySelector("span.is-pulled-right i");
        rightIcon.classList.remove("mdi-chevron-up");
        rightIcon.classList.add("mdi-chevron-down");
}

// If "step" sections are clicked, expand that section and collapse the others
const jsSteps = document.getElementsByClassName("js-step");
const resetSteps = (event) => {
    // If clicked area has no parentNode (i.e. clicked target is destroyed and recreated), just leave be
    if (! event.target.parentNode) return;

    // ? Instead of adding and removing ".step-active", should we toggle ".is-hidden" on the collapsable parts?

    // If clicked on existing section, collapse it
    const currentStep = event.target.closest(".js-step");
    if (currentStep.classList.contains("step-active")) {
        // If collapsable contents were clicked, do nothing
        if (event.target.closest(".js-step-collapsable")) {
            return;
        }

        // Otherwise, collapse the section
        collapseJsStep(currentStep);
        return
    };

    // Make every step uniform
    for (const jsStep of jsSteps) {
        collapseJsStep(jsStep);
        // Collapse all "non-header" parts of step (see common.css)
    }

    // NOTE: If manually clicked and then click event fires the element will collapse, which is not user-friendly.
    // So ensure all triggered click event happens after all API stuff is done.

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
