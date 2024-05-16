let CURRENT_USER;
let SIDEBAR_COLLAPSED = false;
let SITE_PREFS = null;

// Handle unhandled promise rejections (general catch-all for anything that wasn't caught)
// https://developer.mozilla.org/en-US/docs/Web/API/Window/unhandledrejection_event
window.addEventListener("unhandledrejection", function (event) {
    createToast("Something went wrong. Please contact the gEAR team and provide steps to reproduce.");
});

document.addEventListener('DOMContentLoaded', () => {
    // load the site preferences JSON file, then call any functions which need it
    getDomainPreferences().then((result) => {
        SITE_PREFS = result;
        loadPlugins();

        // Replace "gEAR" in the page title with the "domain_short_display_label" from the site preferences
        const pageTitle = document.title;
        document.title = pageTitle.replace("gEAR", SITE_PREFS.domain_short_display_label);
    });


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

    // If the user has agreed to the site's beta status, don't show the modal
    const betaSiteModal = document.getElementById('beta-site-modal');
    const betaCookie = Cookies.get('gear_beta_agreed');
    if (betaCookie != "true") {
        betaSiteModal.classList.add('is-active');
    }

    /**
     * Temporary code to handle the warning modal while in beta mode
     */
    document.getElementById('beta-modal-agree').addEventListener('click', () => {
        Cookies.set('gear_beta_agreed', 'true', { expires: 7 });
        betaSiteModal.classList.remove('is-active');
    });

/**
 * Controls for the left navbar visibility
 */
    const navbarElementsToAnimate = document.querySelectorAll('.icon-text-part');

    /**
     * Toggles a class on the selected elements.
     *
     * @param {string} selector - The CSS selector for the elements to toggle the class on.
     * @param {string} className - The class name to toggle.
     * @param {boolean} add - Whether to add or remove the class. If true, the class will be added; if false, the class will be removed.
     */
    const toggleClass = (selector, className, add) => {
        document.querySelectorAll(selector).forEach((element) => {
            if (add) {
                element.classList.add(className);
            } else {
                element.classList.remove(className);
            }
        });
    }

    /**
     * Handles tooltips for elements with the class 'icon-image-part'.
     * @param {boolean} add - Indicates whether to add or remove tooltips.
     */
    const handleTooltips = (add) => {
        document.querySelectorAll('span.icon-image-part').forEach((element) => {
            if (add) {
                const text = element.parentNode.querySelector('span.icon-text-part').textContent;
                element.setAttribute('data-tooltip', text);
                element.classList.add('has-tooltip-right', 'has-tooltip-arrow');
            } else {
                element.removeAttribute('data-tooltip');
                element.classList.remove('has-tooltip-right', 'has-tooltip-arrow');
            }
        });
    }

    /**
     * Hides the navbar elements with animation.
     */
    const hideNavbarElementsWithAnimation = () => {
        toggleClass('span.menu-label-text', 'is-hidden', true);
        toggleClass('span.icon-text-part', 'is-hidden', true);
        toggleClass('#navbar-logo-normal', 'is-hidden', true);
        toggleClass('#navbar-logo-small', 'is-hidden', false);
        toggleClass('#citation-c', 'is-hidden', true);
        toggleClass("#navbar-toggler i", "mdi-arrow-collapse-left", false);
        toggleClass("#navbar-toggler i", "mdi-arrow-collapse-right", true);

        navbarElementsToAnimate.forEach((element) => {
            element.classList.add('hidden-sidenavbar');
            element.classList.remove('shown-sidenavbar');
        });

        handleTooltips(true);
    }

    /**
     * Shows the navbar elements with animation.
     */
    const showNavbarElementsWithAnimation = () => {
        toggleClass('span.menu-label-text', 'is-hidden', false);
        toggleClass('span.icon-text-part', 'is-hidden', false);
        toggleClass('#navbar-logo-normal', 'is-hidden', false);
        toggleClass('#navbar-logo-small', 'is-hidden', true);
        toggleClass('#citation-c', 'is-hidden', false);
        toggleClass("#navbar-toggler i", "mdi-arrow-collapse-left", true);
        toggleClass("#navbar-toggler i", "mdi-arrow-collapse-right", false);

        navbarElementsToAnimate.forEach((element) => {
            element.classList.add('shown-sidenavbar');
            element.classList.remove('hidden-sidenavbar');
        });

        handleTooltips(false);
    }

    const navbarToggler = document.querySelector('#navbar-toggler');

    navbarToggler.addEventListener('click', (event) => {
        if (SIDEBAR_COLLAPSED == false) {
            hideNavbarElementsWithAnimation();
            SIDEBAR_COLLAPSED = true;
            Cookies.set('gear_sidebar_collapsed', true, { expires: 7 });
            return;
        }
        showNavbarElementsWithAnimation();
        SIDEBAR_COLLAPSED = false;
        Cookies.set('gear_sidebar_collapsed', false, { expires: 7 });
    });

    // now, if the page was initially loaded check and see if this has already been toggled via a cookie
    const sidebar_cookie = Cookies.get('gear_sidebar_collapsed');
    if (sidebar_cookie == "true") {
        hideNavbarElementsWithAnimation();
        SIDEBAR_COLLAPSED = true;
    } else {
        showNavbarElementsWithAnimation();
        SIDEBAR_COLLAPSED = false;
    }

    document.querySelector('#epiviz-panel-designer-link').addEventListener('click', (event) => {
        createToast("This feature is not yet available.", "is-warning");
    });

/**
 * / End controls for the left navbar visibility
 */

    /*************************************************************************************
     Code related to the login process, which is available in the header across all pages.
    *************************************************************************************/

    document.getElementById('submit-login').addEventListener('click', (event) => {
        event.preventDefault(); // Prevent the default form submission

        // reset any UI elements
        document.getElementById('user-email-help').classList.add('is-hidden');
        document.getElementById('user-password-help').classList.add('is-hidden');

        // try to log in
        doLogin();
    });

    document.getElementById('submit-logout').addEventListener('click', (event) => {
        // Clear session information and redirect to home page
        Cookies.remove('gear_session_id');
        CURRENT_USER = undefined;
        window.location.replace('./index.html');
    });


    checkForLogin();
});

const getDomainPreferences = async () => {
    const response = await fetch('/site_domain_prefs.json');
    return response.json();
}

/**
 * Retrieves the value of a URL parameter.
 * @param {string} sParam - The name of the parameter to retrieve.
 * @returns {string|null} - The value of the parameter, or null if it doesn't exist.
 */
const getUrlParameter = (sParam) => {
    const urlParams = new URLSearchParams(window.location.search);
    if (urlParams.has(sParam)) {
        return urlParams.get(sParam);
    }
    return null;
}

/**
 * Opens a modal by adding the 'is-active' class to the specified element.
 * @param {HTMLElement} $el - The element to open as a modal.
 */
const openModal = ($el) => {
    $el.classList.add('is-active');
}

/**
 * Closes the modal by removing the 'is-active' class from the specified element.
 * @param {HTMLElement} $el - The element representing the modal.
 */
const closeModal = ($el) => {
    $el.classList.remove('is-active');
}

/**
 * Closes all modals on the page.
 */
const closeAllModals = () => {
    (document.querySelectorAll('.modal') || []).forEach(($modal) => {
        closeModal($modal);
    });
}

/**
 * Checks if the user is logged in and performs necessary UI updates.
 * @returns {Promise<void>} A promise that resolves once the login check is complete.
 */
const checkForLogin = async () => {
    const session_id = Cookies.get('gear_session_id');
    apiCallsMixin.sessionId = session_id;
    CURRENT_USER = new User();

    if (! session_id || session_id === "undefined") {
        // no cookie found, so user is not logged in
        handleLoginUIUpdates();
    } else {
        // we found a cookie, so see if it's a valid session
        try {
            const data = await apiCallsMixin.getSessionInfo();

            if (data.success) {
                CURRENT_USER = new User({session_id, ...data});
                CURRENT_USER.setDefaultProfile();
                document.getElementById('current-user-name').textContent = data.user_name;
                handleLoginUIUpdates();

                apiCallsMixin.colorblindMode = CURRENT_USER.colorblind_mode;

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

/**
 * Performs the login process.
 * @async
 * @function doLogin
 * @returns {Promise<void>}
 */
const doLogin = async (do_reload) => {
    const formdata = new FormData(document.getElementById("login-form"));
    const data = await apiCallsMixin.login(formdata);

    if (do_reload === undefined) {
        do_reload = true;
    }

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
        apiCallsMixin.sessionId = CURRENT_USER.session_id;
        apiCallsMixin.colorblindMode = CURRENT_USER.colorblind_mode;

        // refresh the page
        if (do_reload) {
            window.location.reload();
        }

    } else {
        // Something went wrong
        console.error("Invalid session_id returned from login.cgi");

    }
}

/**
 * Disables and hides an element.
 * @param {HTMLElement} element - The element to disable and hide.
 * @param {boolean} hideParent - Indicates whether to hide the parent element instead of the element itself.
 */
const disableAndHideElement = (element, hideParent=false) => {
    element.disabled = true;
    if (hideParent) {
        // This is useful for where the parent is a ".control" and
        // removing the element itself still leaves a gap.
        element.parentNode.classList.add('is-hidden');
    } else {
        element.classList.add('is-hidden');
    }
}

/**
 * Enables and shows the specified element.
 * @param {HTMLElement} element - The element to enable and show.
 * @param {boolean} showParent - Indicates whether to show the parent element instead of the element itself.
 */
const enableAndShowElement = (element, showParent=false) => {
    element.disabled = false;
    if (showParent) {
        element.parentNode.classList.remove('is-hidden');
    } else {
        element.classList.remove('is-hidden');
    }
}

/**
 * Hides elements with the class 'logged-in'.
 */
const hideLoggedInElements = () => {
    document.querySelectorAll('.logged-in').forEach(element => element.classList.add('is-hidden'));
}

/**
 * Hides elements with the class 'not-logged-in'.
 */
const hideNotLoggedInElements = () => {
    document.querySelectorAll('.not-logged-in').forEach(element => element.classList.add('is-hidden'));
}

/**
 * Shows the logged-in elements by setting their display property to an empty string.
 */
const showLoggedInElements = () => {
    document.querySelectorAll('.logged-in').forEach(element => element.classList.remove('is-hidden'));
}

/**
 * Shows the elements that are only visible when the user is not logged in.
 */
const showNotLoggedInElements = () => {
    document.querySelectorAll('.not-logged-in').forEach(element => element.classList.remove('is-hidden'));
}

/**
 * Handles the UI updates for the login functionality.
 */
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

/**
 * Triggers an event on the given element.
 * @param {HTMLElement} el - The element on which to trigger the event.
 * @param {string|Function|Event} eventType - The type of event to trigger, or a function to execute as an event handler, or a custom event object.
 */
const trigger = (el, eventType) => {
    // Equivalent to jQuery "trigger" (https://youmightnotneedjquery.com/#trigger_native)
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

/**
 * Logs the error details to the console.
 * @param {Error} error - The error object.
 */
const logErrorInConsole = (error) => {
    if (error.response) {
        // The request was made and the server responded with a status code
        // that falls out of the range of 2xx
        console.error('Response Error:', {
            data: error.response.data,
            status: error.response.status,
            headers: error.response.headers,
        });
    } else if (error.request) {
        // The request was made but no response was received
        // `error.request` is an instance of XMLHttpRequest in the browser and an instance of
        // http.ClientRequest in node.js
        console.error('Request Error:', error.request);
    } else {
        // Something happened in setting up the request that triggered an Error
    }

    if (error.config) {
        console.error('Config:', error.config);
    }

    console.error('Stack Trace:', error);
}

/**
 * Copies the specified text to the clipboard.
 * @param {string} text - The text to be copied.
 * @returns {Promise<boolean>} - A promise that resolves to true if the text was successfully copied, false otherwise.
 */
const copyToClipboard = async (text) => {
    try {
        await navigator.clipboard.writeText(text);
        return true;
    } catch (ex) {
        console.warn("Copy to clipboard failed.", ex);
        return false;
    }
}

// Convert POST payload to FormData so that POSTing to CGI scripts that use cgi.FieldStorage will work
/**
 * Converts an object into FormData.
 * @param {Object} object - The object to be converted.
 * @returns {FormData} - The converted FormData object.
 */
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

/**
 * Creates a toast notification with the given message and level class.
 * @param {string} msg - The message to display in the toast notification.
 * @param {string} [levelClass="is-danger"] - The level class for the toast notification. Defaults to "is-danger".
 */
const createToast = (msg, levelClass="is-danger") => {
    const toast = document.createElement("div");
    toast.classList.add("notification", "js-toast", levelClass, "animate__animated", "animate__fadeInUp", "animate__faster");
    const toastButton = document.createElement("button");
    toastButton.classList.add("delete");
    toastButton.addEventListener("click", (event) => {
        const notification = event.target.closest(".js-toast.notification");
        notification.remove(notification);
    });
    toast.appendChild(toastButton);
    toast.appendChild(document.createTextNode(msg));

    const numToasts = document.querySelectorAll(".js-toast.notification").length;

    if (document.querySelector(".js-toast.notification")) {
        // If .js-toast notifications are present, append under final notification
        // This is to prevent overlapping toast notifications
        const notifications = document.querySelectorAll(".js-toast.notification")
        const notification = notifications[notifications.length - 1];
        notification.insertAdjacentElement("afterend", toast);
        // Position new toast under previous toast with CSS
        toast.style.setProperty("top", `${(numToasts * 70) + 30}px`);
    } else {
        // Otherwise prepend to top of main content
        document.getElementById("main-c").prepend(toast);
    }

    // For a success message, remove it after 3 seconds
    if (levelClass === "is-success") {
        const notifications = document.querySelectorAll(".js-toast.notification")
        const notification = notifications[notifications.length - 1];
        notification.classList.remove("animate__fadeInUp");
        notification.classList.remove("animate__faster");
        notification.classList.add("animate__fadeOutDown");
        notification.classList.add("animate__slower");
    }

    // remove the toast
    toast.addEventListener("animationend", (event) => {
        if (event.animationName === "fadeOutDown") {
            event.target.remove();
        }
    });
}

const loadPlugin = (pluginName, pageName, pluginImportBasename) => {
    const head = document.getElementsByTagName('head')[0];
    const body = document.getElementsByTagName('body')[0];

    // create a hidden element at the end of the document and put it there
    const pluginImportHtmlUrl = `./plugins/${pluginName}/${pageName}`;
    const pluginHtmlElement = document.createElement('div');
    pluginHtmlElement.id = `${pluginName}_html_c`;  // TODO: rename in kebab-case, which requires fixing in plugin HTML and JS files

    fetch(pluginImportHtmlUrl)
        .then(response => {
            return response.text();
        })
        .then(data => {
            body.append(pluginHtmlElement);
            document.getElementById(pluginHtmlElement.id).innerHTML = data;
        });

    const pluginImportCssUrl = `./plugins/${pluginName}/${pluginImportBasename}.css`;
    const style = document.createElement('link');
    style.href = pluginImportCssUrl;
    style.type = 'text/css';
    style.rel = 'stylesheet';
    head.append(style);

    const pluginImportJsUrl = `./plugins/${pluginName}/${pluginImportBasename}.js`;
    const script = document.createElement('script');
    script.src = pluginImportJsUrl;
    script.type = 'text/javascript';
    head.append(script);
}

const loadPlugins = () => {
    const pageName = location.pathname.replace(/^\//, '') || 'index.html';
    const pluginImportBasename = pageName.replace(/\.[^/.]+$/, "");

    for (const [pluginName, pluginPageNames] of Object.entries(SITE_PREFS['enabled_plugins'])) {
        if (pluginPageNames.includes(pageName)) {
            loadPlugin(pluginName, pageName, pluginImportBasename);
        }
    }
}

/**
 * Collapses a JS step element.
 *
 * @param {HTMLElement} stepElt - The step element to collapse.
 */
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
/**
 * Resets the steps and collapses the sections based on the event target.
 * If the clicked area has no parentNode, the function will exit.
 * If the clicked target is an existing section, it will collapse the section.
 * If the clicked target is a collapsable content, it will do nothing.
 * Otherwise, it will collapse all steps and make the clicked step active.
 *
 * @param {Event} event - The event object triggered by the click.
 */
const resetSteps = (event) => {
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

// Generates an RFC4122 version 4 compliant UUID
const uuid = () => {
    return 'xxxxxxxx-xxxx-xxxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
      const r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
      return v.toString(16);
    });
}

for (const jsStep of jsSteps) {
    // Add "capture=true" to trigger parent before children events
    // This helps with manually triggering click events to open up the next section.
    jsStep.addEventListener("click", resetSteps, {capture: true});
}

/* API Calls Mixin */

/* NOTES
Any axios methods that impolement these calls, must provide their own success/error handling
? Some "classes" JS files like "genecart" are already abstracted away... should we keep their methods or move them here?
*/

/**
 * Mixin containing various API calls for interacting with datasets, displays, and analyses.
 * @mixin
 */
const apiCallsMixin = {
    sessionId: null,
    colorblindMode: null,

    // TODO: change all server side scripts with "layout id" to "layout share id" or "collection share id"
    // TODO: change all "gene cart" function names to "gene list"


    /**
     * Adds a display to the collection.
     *
     * @param {string} layoutShareId - The layout share ID.
     * @param {string} displayId - The display ID.
     * @returns {Promise<any>} - A promise that resolves to the response data.
     */
    async addDisplayToCollection(layoutShareId, displayId) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId, display_id: displayId};
        const {data} = await axios.post("cgi/add_display_to_layout.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Creates a new dataset collection with the specified name.
     * @param {string} collectionName - The name of the collection to add.
     * @returns {Promise<any>} - A promise that resolves with the response data.
     */
    async createDatasetCollection(collectionName) {
        const payload = {session_id: this.sessionId, layout_name: collectionName};
        const {data} = await axios.post("cgi/add_layout.cgi", convertToFormData(payload));
        return data;
    },

    /**
     * Checks for projection using the specified dataset ID, pattern source, algorithm, and other options.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} patternSource - The source of the pattern.
     * @param {string} algorithm - The algorithm to use for projection.
     * @returns {Promise<any>} - A promise that resolves to the projection data.
     */
    async checkForProjection(datasetId, patternSource, algorithm) {
        const payload = { session_id: this.sessionId, genecart_id: patternSource, algorithm };
        const {data} = await axios.post(`api/projectr/${datasetId}/output_file`, payload);
        return data;
    },
    /**
     * Deletes a dataset with the specified dataset ID.
     *
     * @param {string} datasetId - The ID of the dataset to be deleted.
     * @returns {Promise<any>} - A promise that resolves with the response data from the server.
     */
    async deleteDataset(datasetId) {
        const payload = {session_id: this.sessionId, dataset_id: datasetId};
        const {data} = await axios.post("cgi/remove_dataset.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Deletes a display from a collection.
     *
     * @param {string} layoutShareId - The layout share ID.
     * @param {string} displayId - The display ID.
     * @returns {Promise<any>} - A promise that resolves with the response data.
     */
    async deleteDisplayFromCollection(layoutShareId, displayId) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId, display_id: displayId};
        const {data} = await axios.post("cgi/remove_display_from_layout.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Deletes a display.
     * @param {string} displayId - The ID of the display to be deleted.
     * @returns {Promise<null>} - A promise that resolves to null.
     */
    async deleteDisplay(displayId) {
        const payload = {session_id: this.sessionId, id: displayId};
        await axios.post("/cgi/delete_dataset_display.cgi", convertToFormData(payload));
        return null;
    },
    /**
     * Deletes a dataset collection with the specified layout ID.
     *
     * @param {string} layoutShareId - The share ID of the dataset collection.
     * @returns {Promise<any>} - A promise that resolves with the response data from the server.
     */
    async deleteDatasetCollection(layoutShareId) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId};
        const {data} = await axios.post("cgi/remove_layout.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Deletes a gene list.
     *
     * @param {string} geneListId - The ID of the gene list to be deleted.
     * @returns {Promise<any>} - A promise that resolves to the response data from the server.
     */
    async deleteGeneList(geneListId) {
        const payload = {session_id: this.sessionId, gene_list_id: geneListId};
        const {data} = await axios.post("/cgi/remove_gene_cart.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches aggregations for a given dataset and analysis.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis.
     * @param {object} filters - The filters to apply to the aggregations.
     * @returns {Promise<object>} - The aggregated data.
     */
    async fetchAggregations(datasetId, analysisId, filters) {
        const payload = {session_id: this.sessionId, dataset_id: datasetId, analysis_id: analysisId, filters};
        const {data} = await axios.post(`/api/h5ad/${datasetId}/aggregations`, payload);
        return data;
    },
    /**
     * Fetches analyses for a given dataset.
     * @param {string} datasetId - The ID of the dataset.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchAnalyses(datasetId) {
        const {data} = await axios.get(`./api/h5ad/${datasetId}/analyses`)
        return data;
    },
    /**
     * Fetches the available plot types for a given dataset and analysis.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis.
     * @param {boolean} [isMultigene=false] - Indicates whether the plot types are for a multigene analysis.
     * @returns {Promise<any>} - A promise that resolves to the available plot types data.
     */
    async fetchAvailablePlotTypes(datasetId, analysisId, isMultigene=false){
        const flavor = isMultigene ? "mg_availableDisplayTypes" : "availableDisplayTypes";
        const payload = {session_id: this.sessionId, dataset_id: datasetId, analysis_id: analysisId};
        const {data} = await axios.post(`/api/h5ad/${datasetId}/${flavor}`, payload);
        return data;
    },
    /**
     * Fetches dashboard data from the server.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysis - The type of analysis.
     * @param {string} plotType - The type of plot.
     * @param {object} plotConfig - The configuration for the plot.
     * @param {object} [otherOpts={}] - Additional options for the request.
     * @returns {Promise<object>} - The fetched data.
     */
    async fetchDashData(datasetId, analysis, plotType, plotConfig, otherOpts={}) {
        // NOTE: gene_symbol should already be already passed to plotConfig
        const payload = { ...plotConfig, plot_type: plotType, analysis, colorblind_mode: this.colorblindMode };
        const {data} = await axios.post(`/api/plot/${datasetId}/mg_dash`, payload, otherOpts);
        return data;
    },
    /**
     * Fetches all datasets asynchronously.
     * @returns {Promise<any>} The fetched data.
     */
    async fetchAllDatasets() {
        const payload = {session_id: this.sessionId};
        const {data} = await axios.post("cgi/get_h5ad_dataset_list.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches dataset comparison data.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {object} filters - The filters to apply to the dataset.
     * @param {string} compareKey - The key to compare the dataset.
     * @param {string} conditionX - The X condition for comparison.
     * @param {string} conditionY - The Y condition for comparison.
     * @param {number} foldChangeCutoff - The fold change cutoff value.
     * @param {number} stDevNumCutoff - The standard deviation number cutoff value.
     * @param {number} logBase - The base for logarithmic transformation.
     * @param {string} statisticalTestAction - The statistical test action to perform.
     * @returns {Promise<any>} The dataset comparison data.
     */
    async fetchDatasetComparison(datasetId, filters, compareKey, conditionX, conditionY, foldChangeCutoff, stDevNumCutoff, logBase, statisticalTestAction) {
        const payload = {
            dataset_id: datasetId,
            obs_filters: filters,
            compare_key: compareKey,
            condition_x: conditionX,
            condition_y: conditionY,
            fold_change_cutoff: foldChangeCutoff,
            std_dev_num_cutoff: stDevNumCutoff,
            log_transformation: logBase,
            statistical_test: statisticalTestAction
        };
		const {data} = await axios.post("cgi/get_dataset_comparison.cgi", convertToFormData(payload));
		return data;
    },
    /**
     * Fetches dataset information from the server.
     * @param {string} datasetId - The ID of the dataset to fetch information for.
     * @returns {Promise<Object>} - A promise that resolves to the dataset information.
     */
    async fetchDatasetInfo(datasetId) {
        const params = {dataset_id: datasetId, include_shape: 1};
        const {data} = await axios.post("cgi/get_dataset_info.cgi", convertToFormData(params));
        return data;
    },
    /**
     * Fetches the dataset collection members for a given layout ID.
     *
     * @param {string} layoutShareId - The share ID of the layout.
     * @returns {Promise<any>} - A promise that resolves to the dataset collection members.
     */
    async fetchDatasetCollectionMembers(layoutShareId) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId};
        const {data} = await axios.post("cgi/get_users_layout_members.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches dataset collections.
     *
     * @param {string|null} layoutShareId - The layout share ID.
     * @param {string|null} noDomain - If 1, the domain layout will not be included in the results.
     * @returns {Promise<any>} The response data.
     */
    async fetchDatasetCollections(layoutShareId=null, noDomain=0) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId, no_domain: noDomain};
        const {data} = await axios.post("cgi/get_user_layouts.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches the display image for a dataset.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} displayId - The ID of the display.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchDatasetDisplayImage(datasetId, displayId) {
        // POST due to payload variables being sensitive
        const payload = {dataset_id: datasetId, display_id: displayId};
        const {data} = await axios.post("/cgi/get_dataset_display_image.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches dataset displays.
     * @param {string} datasetId - The ID of the dataset.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchDatasetDisplays(datasetId) {
        const payload = {session_id: this.sessionId, dataset_id: datasetId};
        const {data} = await axios.post("/cgi/get_dataset_displays.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches dataset list information for the user, such as metadata and default display settings.
     * This can be used for a single dataset id or a layout id.
     * @param {Object} requestConfig - The request configuration.
     * @returns {Promise<Object>} The dataset list information.
     */
    async fetchDatasetListInfo(requestConfig) {
        const payload = {session_id: this.sessionId, ...requestConfig};
        const {data} = await axios.post("/cgi/get_dataset_list.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches datasets based on the provided search criteria.
     *
     * @param {Object} searchCriteria - The search criteria for filtering datasets.
     * @returns {Promise} A promise that resolves to the fetched datasets.
     */
    async fetchDatasets(searchCriteria) {
        const payload = {session_id: this.sessionId, ...searchCriteria};
        const {data} = await axios.post("/cgi/search_datasets.cgi", convertToFormData(payload));
        return data;
    },

    /**
     * Fetches the default display for a dataset.
     * @param {string} datasetId - The ID of the dataset.
     * @param {boolean} [isMultigene=false] - Indicates if the dataset is multigene.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchDefaultDisplay(datasetId, isMultigene=false) {
        const payload = {session_id: this.sessionId, dataset_id: datasetId, is_multigene: isMultigene};
        const {data} = await axios.post("/cgi/get_default_display.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches a single display by ID.
     * @param {string} displayId - The ID of the display.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchDisplay(displayId) {
        const payload = {session_id: this.sessionId, display_id: displayId};
        const {data} = await axios.post("/cgi/get_dataset_display.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Fetches the Epiviz display data for a given dataset, gene symbol, and genome.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} geneSymbol - The gene symbol.
     * @param {string} genome - The genome.
     * @param {Object} [otherOpts={}] - Additional options for the axios request.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchEpivizDisplay(datasetId, geneSymbol, genome, otherOpts={}) {

        const urlParams = new URLSearchParams();
        urlParams.append('gene', geneSymbol);
        urlParams.append('genome', genome);

        const {data} = await axios.get(`/api/plot/${datasetId}/epiviz?${urlParams.toString()}`, otherOpts);
        return data;
    },
    /**
     * Fetches annotations for the passed gene symbols
     *
     * @param {string[]} geneSymbols - The gene symbols to search for. Comma-separated.
     * @param {boolean} exactMatch - Indicates whether to perform an exact match search.
     * @param {string} layoutShareId - The layout share ID.
     * @param {boolean} isMulti - Indicates whether multiple gene symbols are being searched.
     * @returns {Promise<any>} - The fetched gene annotations.
     */
    async fetchGeneAnnotations(geneSymbols, exactMatch, layoutShareId, isMulti) {
        const payload = { session_id: this.sessionId, search_gene_symbol: geneSymbols, exact_match: exactMatch, is_multi: isMulti, layout_share_id: layoutShareId };
        const {data} = await axios.post(`/cgi/search_genes.cgi`, convertToFormData(payload));
        return data;
    },

    /**
     * Fetches gene lists based on the provided search criteria.
     *
     * @param {Object} searchCriteria - The search criteria for fetching gene lists.
     * @returns {Promise} - A promise that resolves to the fetched gene lists.
     */
    async fetchGeneLists(searchCriteria) {
        const payload = {session_id: this.sessionId, ...searchCriteria};
        const {data} = await axios.post(`/cgi/search_gene_carts.cgi`, convertToFormData(payload));
        return data;
    },
    /**
     * Fetches the members of a gene cart.
     * @param {string} geneCartId - The ID of the gene cart.
     * @returns {Promise<any>} - A promise that resolves to the data of the gene cart members.
     */
    async fetchGeneCartMembers(geneCartId) {
        const payload = { session_id: this.sessionId, gene_cart_id: geneCartId };
        const {data} = await axios.post(`/cgi/get_gene_cart_members.cgi`, convertToFormData(payload));
        return data;
    },
    /**
     * Fetches gene carts based on the specified cart type.
     * @param {string} cartType - The type of gene cart to fetch.
     * @returns {Promise<any>} - A promise that resolves to the fetched gene carts data.
     */
    async fetchGeneCarts(cartType=null) {
        const payload = {session_id: this.sessionId};
        if (cartType) {
            payload.cart_type = cartType;
        }

        const {data} = await axios.post(`/cgi/get_user_gene_carts.cgi`, convertToFormData(payload));
        return data;
    },
    /**
     * Fetches gene symbols for a given dataset and analysis.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysisId - The ID of the analysis (optional).
     * @returns {Promise<Array<string>>} - A promise that resolves to an array of gene symbols.
     */
    async fetchGeneSymbols(datasetId, analysisId) {
        let url = `./api/h5ad/${datasetId}/genes`;
        if (analysisId) url += `?analysis_id=${analysisId}`;
        const {data} = await axios.get(url);
        return data;
    },
    /**
     * Fetches H5ad info from the server.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} [analysisId] - The ID of the analysis (optional).
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchH5adInfo(datasetId, analysisId) {
        let url = `/api/h5ad/${datasetId}`
        if (analysisId) url += `?analysis_id=${analysisId}`;
        const {data} = await axios.get(url);
        return data;
    },
    /**
     * Fetches info on all the organisms in the database
     * @returns {Promise<any>} - A promise that resolves to the list of organisms
     */
    async fetchOrganismList() {
        const {data} = await axios.post(`/cgi/get_organism_list.cgi`);
        return data;
    },
    /**
     * Fetches orthologs for a given dataset and gene symbols.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string[]} geneSymbols - An array of gene symbols.
     * @returns {Promise<any>} - A promise that resolves to the fetched ortholog data.
     */
    async fetchOrthologs(datasetId, geneSymbols, geneOrganismId=null) {
        const payload = { session_id: this.sessionId, gene_symbols:geneSymbols, gene_organism_id: geneOrganismId };
        const {data} = await axios.post(`/api/h5ad/${datasetId}/orthologs`, payload);
        return data;
    },
    /**
     * Fetches the pattern element list from the server.
     * @param {string} sourceID - The ID of the source.
     * @param {string} scope - The scope of the pattern element list.
     * @returns {Promise<any>} - A promise that resolves to the pattern element list data.
     */
    async fetchPatternElementList(sourceID, scope) {
        const payload = { source_id: sourceID, scope };
        const {data} = await axios.post(`/cgi/get_pattern_element_list.cgi`, convertToFormData(payload));
        return data;
    },
    /**
     * Fetches pattern weighted genes from the server.
     * @param {string} sourceId - The ID of the source.
     * @param {string} patternId - The ID of the pattern.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchPatternWeightedGenes(sourceId, patternId) {
        const payload = {source_id: sourceId, pattern_id: patternId};
        const {data} = await axios.post(`/cgi/get_pattern_weighted_genes.cgi`, convertToFormData(payload));
        return data;
    },
    /**
     * Fetches Plotly data for a given dataset, analysis, plot type, and plot configuration.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysis - The analysis type.
     * @param {string} plotType - The type of plot.
     * @param {object} plotConfig - The configuration object for the plot.
     * @param {object} [otherOpts={}] - Additional options for the request.
     * @returns {Promise<object>} - The fetched Plotly data.
     */
    async fetchPlotlyData(datasetId, analysis, plotType, plotConfig, otherOpts={}) {
        // NOTE: gene_symbol should already be already passed to plotConfig
        const payload = { ...plotConfig, plot_type: plotType, analysis, colorblind_mode: this.colorblindMode };
        const {data} = await axios.post(`/api/plot/${datasetId}`, payload, otherOpts);
        return data;
    },

    /**
     * Fetches a projection from the server.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} projectionId - The ID of the projection.
     * @param {string} patternSource - The source of the pattern.
     * @param {string} algorithm - The algorithm used for projection.
     * @param {string} gctype - The type of gene cart.
     * @param {Object} [otherOpts={}] - Additional options for the request.
     * @returns {Promise} - A promise that resolves with the fetched data.
     */
    async fetchProjection(datasetId, projectionId, patternSource, algorithm, gctype, otherOpts={}) {
        const payload = { session_id: this.sessionId, projection_id: projectionId, genecart_id: patternSource, algorithm, scope: gctype};
        const {data} = await axios.post(`api/projectr/${datasetId}`, payload, otherOpts);
        return data;
    },

    /**
     * Fetches SVG data for a given dataset, gene symbol, and projection ID.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} geneSymbol - The gene symbol.
     * @param {string} projectionId - The ID of the projection.
     * @param {Object} [otherOpts={}] - Additional options for the GET request.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchSvgData(datasetId, geneSymbol, projectionId=null, otherOpts={}) {

        // create URL Params for the GET request
        const urlParams = new URLSearchParams();
        urlParams.append('gene', geneSymbol);
        if (projectionId) {
            urlParams.append('projection_id', projectionId);
        }

        // fetch the data
        const {data} = await axios.get(`/api/plot/${datasetId}/svg?${urlParams.toString()}`, otherOpts);

        return data;
    },

    /**
     * Fetches the TSNE image for a given dataset, analysis, plot type, plot configuration, and other options.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} analysis - The analysis type.
     * @param {string} plotType - The type of plot.
     * @param {object} plotConfig - The plot configuration.
     * @param {object} [otherOpts={}] - Additional options for the request.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchTsneImage(datasetId, analysis, plotType, plotConfig, otherOpts={}) {
        // NOTE: gene_symbol should already be already passed to plotConfig
        const payload = { ...plotConfig, plot_type: plotType, analysis, colorblind_mode: this.colorblindMode };
        const {data} = await axios.post(`/api/plot/${datasetId}/tsne`, payload, otherOpts);
        return data;
    },
    /**
     * Fetches user history entries.
     * @param {number} numEntries - The number of entries to fetch.
     * @returns {Promise<any>} - A promise that resolves to the fetched data.
     */
    async fetchUserHistoryEntries(numEntries) {
        const payload = { session_id: this.sessionId, num_entries: numEntries };
        const {data} = await axios.post("/cgi/get_user_history_entries.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Retrieves session information.
     * @returns {Promise<Object>} The session information.
     */
    async getSessionInfo() {
        const payload = {session_id: this.sessionId};
        const {data} = await axios.post("/cgi/get_session_info.v2.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Logs in the user.
     * @param {FormData} formData - The form data containing the login credentials.
     * @returns {Promise<any>} - A promise that resolves to the login response data.
     */
    async login(formData) {
        const payload = new URLSearchParams(formData);
        const {data} = await axios.post("/cgi/login.v2.cgi", payload);
        return data;
    },

    /**
     * Renames a dataset collection.
     *
     * @param {string} layoutShareId - The ID of the layout share.
     * @param {string} collectionName - The new name for the collection.
     * @returns {Promise<any>} - A promise that resolves with the response data.
     */
    async renameDatasetCollection(layoutShareId, collectionName) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId, layout_name: collectionName};
        const {data} = await axios.post("cgi/rename_layout.cgi", convertToFormData(payload));
        return data;
    },

    /**
     * Saves the dataset display with the specified parameters.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} displayId - The ID of the display.
     * @param {string} label - The label for the display.
     * @param {string} plotType - The type of plot for the display.
     * @param {object} plotConfig - The configuration object for the plot.
     * @returns {Promise<any>} - A promise that resolves to the saved data.
     */
    async saveDatasetDisplay(datasetId, displayId, label, plotType, plotConfig) {
        const payload = {
            id: displayId,
            dataset_id: datasetId,
            session_id: this.sessionId,
            label,
            plot_type: plotType,
            plotly_config: JSON.stringify({
                ...plotConfig,  // depending on display type, this object will have different properties
            }),
        };
        if (!displayId) delete payload.id;  // Prevent passing in "null" as a string.

        const {data} = await axios.post("/cgi/save_dataset_display.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Saves the changes made to the dataset information.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} visibility - The visibility of the dataset. 1 for public, 0 for private.
     * @param {string} isDownloadable - Indicates if the dataset is downloadable. 1 for downloadable, 0 for not downloadable.
     * @param {string} title - The title of the dataset.
     * @param {string} pubmedId - The PubMed ID of the dataset.
     * @param {string} geoId - The GEO ID of the dataset.
     * @param {string} lDesc - The long description of the dataset.
     * @returns {Promise<any>} - A promise that resolves to the response data.
     */
    async saveDatasetInfoChanges(datasetId, visibility, isDownloadable, title, pubmedId, geoId, lDesc) {
        const payload = {session_id: this.sessionId, dataset_id: datasetId, visibility, is_downloadable: isDownloadable, title, pubmed_id: pubmedId, geo_id: geoId, ldesc: lDesc};
        const {data} = await axios.post("/cgi/save_datasetinfo_changes.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Saves the default display for a dataset.
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} displayId - The ID of the display.
     * @param {boolean} [isMultigene=false] - Indicates if the display is for multiple genes.
     * @returns {Promise<any>} - A promise that resolves to the saved data.
     */
    async saveDefaultDisplay(datasetId, displayId, isMultigene=false) {
        const payload = {display_id: displayId, session_id: this.sessionId, dataset_id: datasetId, is_multigene: isMultigene};
        const {data} = await axios.post("/cgi/save_default_display.cgi", convertToFormData(payload));
        return data;
    },

    /**
     * Saves the dataset collection arrangement.
     *
     * @param {string} layoutShareId - The layout share ID.
     * @param {Object} layoutArrangement - The layout arrangement object.
     * @returns {Promise<any>} - A promise that resolves with the response data.
     */
    async saveDatasetCollectionArrangement(layoutShareId, layoutArrangement) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId, layout_arrangement: JSON.stringify(layoutArrangement)};
        const {data} = await axios.post("/cgi/save_layout_arrangement.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Saves the changes made to gene list information.
     *
     * @param {string} gcId - The ID of the gene list.
     * @param {string} visibility - The visibility of the gene list. 1 for public, 0 for private.
     * @param {string} title - The title of the gene list.
     * @param {string} organismId - The ID of the organism.
     * @param {string} ldesc - The description of the gene list.
     * @returns {Promise<any>} - A promise that resolves to the response data.
     */
    async saveGeneListInfoChanges(gcId, visibility, title, organismId, ldesc) {
        const payload = {session_id: this.sessionId, gc_id: gcId, visibility, title, organism_id: organismId, ldesc};
        const {data} = await axios.post("/cgi/save_genecart_changes.cgi", convertToFormData(payload));
        return data;
    },
    /**
         * Saves the default organism for all annotation views.
         * @param {Object} user - The JS object for the currently-logged in user
         * @returns {Promise<any>} - A promise that resolves to the saved data.
         */
    async saveUserDefaultOrgId(user) {
        const payload = {session_id: this.sessionId, default_org_id: user.default_org_id};
        const {data} = await axios.post("/cgi/save_user_default_organism.cgi", convertToFormData(payload));
        return data;
    },
    /**
     * Sets the user's primary dataset collection.
     *
     * @param {string} layoutShareId - The share ID of the layout.
     * @returns {Promise<any>} - A promise that resolves with the response data from the server.
     */
    async setUserPrimaryDatasetCollection(layoutShareId) {
        const payload = {session_id: this.sessionId, layout_share_id: layoutShareId};
        const {data} = await axios.post("/cgi/set_primary_layout.cgi", convertToFormData(payload));
        return data;
    }

}
