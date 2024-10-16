'use strict';

const status2Class = {
    pending : "is-warning"
    , loading : "is-info"
    , completed : "is-success"
    , canceled : "is-dark"  // should not be seen ideally
    , failed : "is-danger"
}

const finishElements = {
    importFinish: document.getElementById("import-finish")
    , partialSuccess: document.getElementById("partial-success")
    , completeSuccess: document.getElementById("complete-success")
}

const POLL_TIMEOUT = 10_000; // 10 seconds

// Error messages
const FETCH_ERROR_MESSAGE = "Could not fetch submission.";
const FETCH_DATASET_ERROR_MESSAGE = "Could not fetch submission dataset.";
const SESSION_ERROR_MESSAGE = "Must be logged in to import a new submission. Please login and refresh page.";
const DB_ERROR_MESSAGE = "Something went wrong with saving this submission to the database. Please contact gEAR support.";


const generateUUIDFromString = (string) => {
    const hash = crypto.createHash('sha256').update(string).digest('hex');
    const uuid = `${hash.substring(0, 8)}-${hash.substring(8, 12)}-4${hash.substring(13, 16)}-8${hash.substring(17, 20)}-${hash.substring(20, 32)}`;
    return uuid;
  }

/**
 * Checks if the user is logged in by verifying the session ID.
 *
 * @returns {boolean} True if the user is logged in, false otherwise.
 */
const isUserLoggedIn = () => {
    return Boolean(CURRENT_USER.session_id);
};

/**
 * Updates the message cell in a table with a specified message and styling based on the level.
 *
 * @param {string} datasetId - The ID of the dataset, used to locate the message cell in the table.
 * @param {string} msg - The message to be displayed in the cell.
 * @param {string} [level] - The optional level of the message, used to apply additional styling (e.g., 'info', 'warning', 'error').
 */
const printTblMessage = (datasetId, msg, level) => {
    // ? Should I make this a hover element since messages may be long?
    const cell = document.querySelector(`#dataset-${datasetId} .js-messages-cell`);
    cell.classList = ["js-messages-cell"];  // Reset classes
    cell.textContent = msg;
    if (level) {
        cell.classList.add(`has-text-${level}-dark`)
    }
};

/* Create a dataset permalink URL from the stored share ID */
/**
 * Generates a permalink URL for a dataset based on the provided share ID.
 *
 * @param {string} shareId - The unique identifier for the dataset to be shared.
 * @returns {string} - The generated permalink URL for the dataset.
 */
const createDatasetPermalinkUrl = (shareId) => {
    const currentUrl = window.location.href;
    const currentPage = currentUrl.lastIndexOf("nemoarchive_import");
    return `${currentUrl.substring(0, currentPage)}p?s=${shareId}`;
};

/**
 * Extracts the name of the first entity with the entityType "file_set" from the given array.
 *
 * @param {Array} jd - The array of entities to search through.
 * @returns {string} - The name of the first entity with the entityType "file_set".
 */
const grabSubmissionId = (jd) => {
    return jd.find(entity => entity.entityType === "file_set").name;
};

/**
 * Filters and returns entities of type "file" from the provided array.
 *
 * @param {Array} jd - The array of entities to filter.
 * @returns {Array} - An array of entities where the entityType is "file".
 */
const getFileEntities = (jd) => {
    return jd.filter(entity => entity.entityType === "file");
};

/**
 * Generates a DOM element from an HTML string.
 *
 * @param {string} html - The HTML string to convert into a DOM element.
 * @returns {Element} The first child element of the generated content.
 */
const generateElements = (html) => {
    const template = document.createElement('template');
    template.innerHTML = html.trim();
    return template.content.children[0];
};

/**
 * Fetches the dataset display information based on the provided display ID.
 *
 * @param {string} displayId - The ID of the display to fetch information for.
 * @returns {Promise<Object>} A promise that resolves to the dataset display information.
 */
const getDatasetDisplay = async (displayId) => {
    const params = new URLSearchParams({"display_id": displayId});

    const data = await axios.get(`./cgi/get_dataset_display.cgi?${params}`);
    return data;
};

/**
 * Fetches the default display ID for a given dataset.
 *
 * @param {string} datasetId - The ID of the dataset for which to fetch the default display.
 * @returns {Promise<string>} - A promise that resolves to the default display ID.
 * @throws {Error} - Throws an error if the default display ID could not be fetched.
 */
const getDefaultDisplay = async (datasetId) => {
    try {
        // POST due to payload variables being sensitive
        const {data} =  await apiCallsMixin.fetchDefaultDisplay(datasetId)
        const {default_display_id: defaultDisplayId} = data;
        return defaultDisplayId;
    } catch (error) {
        logErrorInConsole(error);
        const msg = `Could not fetch default display for dataset ${datasetId}.`
        throw new Error(msg);
    }
};

/**
 * Fetches a submission by its ID.
 *
 * @param {string} submissionId - The ID of the submission to fetch.
 * @returns {Promise<Object>} The submission data.
 * @throws {Error} If the submission could not be fetched.
 */
const getSubmission = async (submissionId) => {

    try {
        const response = await axios.get(`/api/submissions/${submissionId}`);
        const { data } = response;

        if (!data.success) {
            throw new Error(data.message);
        }

        return data;
    } catch (error) {
        logErrorInConsole(error);
        throw new Error(FETCH_ERROR_MESSAGE);
    }
};

/**
 * Fetches the submission dataset from the given URL.
 *
 * @param {string} href - The URL to fetch the submission dataset from.
 * @returns {Promise<Object>} The data from the fetched submission dataset.
 * @throws {Error} If the fetch operation fails or the response indicates failure.
 */
const getSubmissionDataset = async (href) => {
    // ? Should I use submission_id and dataset_id as args instead?
    try {
        const response = await axios.get(href);
        const { data } = response;

        if (!data.success) {
            throw new Error(data.message);
        }

        return data;
    } catch (error) {
        logErrorInConsole(error);
        throw new Error(FETCH_DATASET_ERROR_MESSAGE);
    }
};

/**
 * Initializes new submission-related entries in the database.
 *
 * @param {Array} fileEntities - The file entities to process.
 * @param {string} submissionId - The ID of the submission.
 * @returns {Promise<void>}
 */
const initializeNewSubmission = async (fileEntities, submissionId) => {
    const { datasets } = await processSubmission(fileEntities, submissionId);
    const payload = { action: "reset_steps" };

    for (const datasetId in datasets) {
        try {
            await axios.put(`${datasets[datasetId].href}/status`, payload, {
                headers: { "Content-Type": "application/json" }
            });

            // Set up the dataset row in the submission table
            // ? Consider pagination
            initializeDatasetRow(datasets[datasetId]);

            // Log info about any dataset that could not be initialized
            if (datasets[datasetId].error_message) {
                printTblMessage(datasetId, datasets[datasetId].error_message, "danger")
            }

        } catch (error) {
            logErrorInConsole(error);
            // Log info about any dataset that could not be initialized
            const msg = `Could not initialize new submission for dataset ${datasetId}.`;
            throw new Error(msg);
        }
    }
};

/**
 * Initializes the submission element with the given submission ID.
 *
 * @param {string} submissionId - The ID of the submission.
 */
const initializeSubmissionElement = (submissionId) => {
    const submissionElt = document.getElementById("submission-title");
    submissionElt.dataset.submissionId = submissionId;
    submissionElt.addEventListener("click", handleSubmissionLink);
};

/**
 * Initializes a dataset row in the HTML table with the provided dataset information.
 *
 * @param {Object} dataset - The dataset object containing necessary information.
 * @param {string} dataset.datasetId - The unique identifier for the dataset.
 * @param {string} dataset.shareId - The share identifier for the dataset.
 * @param {string} dataset.identifier - The identifier string for the dataset, which includes the namespace.
 */
const initializeDatasetRow = (dataset) => {
    const {datasetId, shareId, identifier} = dataset;
    const namespace = identifier.split("nemo:")[1];
    const identifierUrl = `https://assets.nemoarchive.org/${namespace}`;

    // Create row template
    const template = document.getElementById("dataset-row-tmpl").content.cloneNode(true);
    const row = template.querySelector("tr");
    row.id = `dataset-${datasetId}`;
    row.dataset.shareId = shareId;
    row.dataset.datasetId = datasetId;

    // first cell
    const identifierCell = row.querySelector(".js-identifier-cell");
    const identifierLink = identifierCell.querySelector("a");
    identifierLink.href = identifierUrl;
    identifierLink.textContent = identifier;
    const permalink = row.querySelector(".js-identifier-permalink");

    // second cell
    //const progress = row.querySelector(".js-progress-cell");

    // third cell
    //const currentStep = row.querySelector(".js-current-step-cell");

    // fourth cell
    //const messages = row.querySelector(".js-messages-cell");

    // fifth cell
    //const obsLevels = row.querySelector(".js-obslevels-cell");

    // Append row to table
    const parent = document.querySelector("#submission-datasets tbody");
    const htmlCollection = generateElements(template);
    parent.append(htmlCollection);
};

/**
 * Checks if the import process is complete based on the given status.
 *
 * @param {Object} status - The status object containing the state of the import process.
 * @param {string} status.make_umap - The status of the make_umap step, which can be "pending", "loading", or other states.
 * @returns {boolean} - Returns true if the import process is complete, otherwise false.
 */
const isImportComplete = (status) => {
    /* Returns True if last step finished */
    // i.e. a non-runnable state. Incomplete steps are reset to "pending" when an import is attempted
    if (status.make_umap == "pending" || status.make_umap == "loading") return false;
    return true;
};

/**
 * Launches the new submission import process.
 *
 * @param {string} submissionId - The ID of the submission.
 * @param {Array} fileEntities - The file entities to import.
 * @param {Array} sampleEntities - The sample entities to import.
 */
const launchSubmissionImport = (submissionId, fileEntities) => {
    const params = { file_metadata: fileEntities, action: "import" };
    axios.post(`/api/submissions/${submissionId}`, params, {
        headers: { "Content-Type": "application/json" },
    });
};

/**
 * Populates the non-supported element with a value from the URL parameter 'notsupported'.
 * If the 'notsupported' parameter is greater than 0, it displays an additional paragraph element.
 *
 * @function
 */
const populateNonSupportedElement = () => {
    const notSupported = getUrlParameter('notsupported') || 0;
    const notSupportedEl = document.getElementById("not-supported");
    notSupportedEl.textContent = notSupported;
    if (parseInt(notSupported) > 0) {
        const notSupportedP = document.getElementById("not-supported-p");
        notSupportedP.style.display = "block";
    }
};

/* Populate select dropdown of observation categorical metadata */
const populateObsDropdown = async (datasetId) => {
    const parent = document.querySelector(`#dataset-${datasetId} .js-obslevels-cell`);
    try {
        const {data} = await axios.get(`/api/h5ad/${datasetId}`);
        const obsLevels = Object.keys(data.obs_levels);
        if (!obsLevels.length) {
            parent.textContent = "Not applicable.";
            return
        }
        parent.classList.add("is-loading", "select");
        const optionElts = obsLevels.map(cat => `<option>${cat}</option>`);
        const template = `<select>${optionElts.join("")}</select>`;
        const htmlCollection = generateElements(template);
        parent.append(htmlCollection);
        parent.classList.remove("is-loading");
    } catch (error) {
        logErrorInConsole(error);
        const msg = `Could not fetch categorical observations for dataset ${datasetId}.`;
        createToast(msg, "is-warning");
    }
};

/**
 * Populates the element with the ID "submission-id" with the provided submission ID.
 *
 * @param {string} submissionId - The submission ID to be displayed in the element.
 */
const populateSubmissionId = (submissionId) => {
    const submissionEl = document.getElementById("submission-id");
    submissionEl.textContent = submissionId;
};

/**
 * Processes a submission by either retrieving an existing submission or creating a new one.
 *
 * @param {Array} fileEntities - An array of file entities to be processed.
 * @param {string} submissionId - The ID of the submission to be processed.
 * @returns {Promise<Object>} - A promise that resolves to the submission data.
 *
 * @throws {Error} - Throws an error if the submission retrieval or creation fails.
 */
const processSubmission = async (fileEntities, submissionId) => {
    const submissionElt = document.getElementById("submission-title");

    try {
        // Going to attempt to get existing submission (and resume)
        const {data} = await axios.get(`/api/submissions/${submissionId}`)
        if (!data.success) throw new Error(data.message);

        // Append extra dataset information
        for (const datasetId in data.datasets) {
            data.datasets[datasetId] = await getSubmissionDataset(data.datasets[datasetId].href);
        }

        // Prepopulate collection name if previously added
        if (data.collection_name) {
            document.getElementById("collection-name").value = data.collection_name;
        }
        // If this is a different user viewing the submission, disable naming the layout
        if (! data.is_submitter) {
            document.getElementById("collection-name").disabled = "disabled";
        }

        submissionElt.dataset.layout_share_id = data.layout_share_id; // Store layout_share_id for future retrieval

        return data

    } catch {
        const isRestricted = 0;    // hardcoded for now

        // If existing submission does not exist, create new submission
        const postParams = {"submission_id": submissionId, "is_restricted": isRestricted };
        const {data} = await axios.post("/api/submissions", postParams, {
            headers: {"Content-Type": "application/json"},
        });
        if (!data.success) throw new Error(data.message);

        submissionElt.dataset.layout_share_id = data.layout_share_id; // Store layout_share_id for future retrieval


        // NOTE: at this point we have no routes for our submission
        // Go through the projected dataset routes, and create the ones that do not already exist
        // If the dataset exists (for another submission), then associate with this submission

        await Promise.allSettled([...fileEntities].map( async (entity) => {
            // Create a UUID based on the identifier to use as the dataset_id
            const identifier = entity.attributes.identifier;
            if (!identifier) {
                // Log error and continue to next entity
                console.warn(`No identifier found for entity ${entity.name}`);
                // increment the not supported counter
                const notSupportedElt = document.getElementById("not-supported");
                notSupportedElt.textContent = parseInt(notSupportedElt.textContent) + 1;
                return;
            }

            const datasetId = generateUUIDFromString(identifier);

            const sdParams = {"dataset_id":datasetId, "identifier":identifier, "is_restricted":isRestricted}

            const href = `/api/submissions/${submissionId}/datasets/${datasetId}`

            // GET route first
            try {
                data.datasets[datasetId] = await getSubmissionDataset(href);
            } catch {
                // Submission dataset did not already exist... create it
                try {
                    const {data: sdData} = await axios.post(`/api/submissions/${submissionId}/datasets`, sdParams, {
                        headers: {"Content-Type": "application/json"}
                    });
                    if (sdData.success) {
                        data.datasets[datasetId] = sdData
                    } else {
                        data.datasets[datasetId].error_message = sdData.message;
                    }

                } catch (error) {
                    data.datasets[datasetId].error_message = error.message;
                }
            }

            // Save dataset as a new SubmissionMember
            await axios.put(`${postJsonRes.datasets[datasetId].href}/members`, {
                headers: {"Content-Type": "application/json"}
            })
            // TODO: Create a rollback script in case of failure of submission or datasets;
        }));
        return data;
    }
};

const saveNewDisplay = async (datasetId, plotConfig, plotType, label) => {
    const data = await apiCallsMixin.saveDatasetDisplay(datasetId, null, label, plotType, plotConfig)
    return data.display_id;
};

/**
 * Saves the default display for a given dataset.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} displayId - The ID of the display to set as default.
 * @returns {Promise<string>} - A promise that resolves to the ID of the default display.
 */
const saveDefaultDisplay = async (datasetId, displayId) => {
    const data = await apiCallsMixin.saveDefaultDisplay(datasetId, displayId);
    return data.default_display_id;
};

/**
 * Updates the layout name for a given submission.
 *
 * @param {string} submissionId - The ID of the submission to update.
 * @param {string} collectionName - The new layout name for the submission.
 * @returns {Promise<Object>} The response data from the server.
 */
const updateLayoutName = async (submissionId, collectionName) => {
    const payload = {submission_id: submissionId, layout_name: collectionName}
    const {data} = await axios.post("./cgi/nemoarchive_update_submission_layout_name.cgi", convertToFormData(payload));
    return data;
};

/**
 * Displays a permalink for a dataset in the specified table row.
 *
 * @param {string} datasetId - The unique identifier of the dataset.
 */
const showDatasetPermalink = (datasetId) => {
    const tableRow = document.getElementById(`dataset-${datasetId}`);
    const permalinkRow = tableRow.querySelector(".js-identifier-permalink");
    const shareId = tableRow.dataset.shareId;
    const permalinkUrl =  createDatasetPermalinkUrl(shareId);
    const template = `<a href=${permalinkUrl} target="_blank">Visualize</a>`;
    permalinkRow.innerHTML = template;
};

/**
 * Determines the current step based on the provided statuses.
 *
 * @param {Object} statuses - An object representing the statuses of various steps.
 * @param {string} statuses.make_umap - Status of the "make_umap" step.
 * @param {string} statuses.convert_to_h5ad - Status of the "convert_to_h5ad" step.
 * @param {string} statuses.convert_metadata - Status of the "convert_metadata" step.
 * @param {string} statuses.pulled_to_vm - Status of the "pulled_to_vm" step.
 * @returns {string} - The current step, which can be "loading", "failed", "make_umap", "convert_to_h5ad", "convert_metadata", or "pulled_to_vm".
 */
const getCurrentStep = (statuses) => {
    // NOTE: This is under the assumption only one status is loading or failed.
    // Is anything currently loading or failed.
    for (const status in statuses) {
        if (statuses[status] == "loading") return status;
        if (statuses[status] == "failed") return status;
    }

    // Find last complete step
    if (statuses.make_umap == "completed") return "make_umap"
    // Just report what was the most recent step
    if (statuses.convert_to_h5ad == "completed") return "convert_to_h5ad"
    if (statuses.convert_metadata == "completed") return "convert_metadata"
    return "pulled_to_vm"

};

/**
 * Updates the status of a dataset's progress table and progress bar.
 *
 * @param {string} datasetId - The unique identifier for the dataset.
 * @param {Object} statuses - An object containing the statuses of each step.
 * @throws {Error} If the current step is not found in the step order.
 */
const updateTblStatus = (datasetId, statuses) => {
    const currentStep = getCurrentStep(statuses);
    const stepOrder = ["pulled_to_vm", "convert_metadata", "convert_to_h5ad", "make_umap"];
    const stepIndex = stepOrder.indexOf(currentStep);
    if (stepIndex == -1) {
        throw new Error(`Could not find index of ${currentStep}`);
    }
    const currentStepStatus = statuses[currentStep];
    const statusCapitalized = currentStepStatus.charAt(0).toUpperCase() + currentStepStatus.slice(1)
    const progressModifier = currentStepStatus == "completed" ? 1 : 0;

    const datasetRow = document.getElementById(`dataset-${datasetId}`);

    // Update progress bar
    const progressBarElt = datasetRow.querySelector(".js-progress-cell")
    const progressMaxValue = (stepIndex + progressModifier) * 25;

    if (currentStepStatus == "loading") {
        // Dynamically update bar to look like it's "loading"
        progressBarElt.value = 0;
        setInterval(() => {
            progressBarElt.value = progressBarElt.value + 1;
            if (progressBarElt.value == progressMaxValue) {
                setTimeout(() => {
                    progressBarElt.value = 0;
                }, 500);
            }
        },25);
    } else {
        progressBarElt.value = progressMaxValue;
    }

    // ? classList is supposed to be read-only so unsure why this works
    progressBarElt.classList = [`progress ${status2Class[statuses[currentStep]]}`];
    progressBarElt.textContent = statusCapitalized;

    const step2Name = {
        "pulled_to_vm": "Pull Files"
        , "convert_metadata": "Validate Metadata"
        , "convert_to_h5ad": "Convert to H5AD"
        , "make_umap": "Make UMAP"
    }

    // Update current step information
    const currentStepElt = datasetRow.querySelector(".js-current-step-cell");
    currentStepElt.textContent = progressBarElt.value == 100 ? "Complete" : step2Name[currentStep];
};

const updateDisplayInLayout = async (submissionId, oldDisplayId, newDisplayId) => {
    const {data} = await axios.post("./cgi/nemoarchive_update_display_in_layout.cgi", convertToFormData({
        submission_id: submissionId
        , old_display_id: oldDisplayId
        , new_display_id: newDisplayId
    }));
    return data;
}

/**
 * Handles the submission link by generating a shareable URL with the submission ID
 * and copying it to the clipboard. Notifies the user whether the URL was successfully
 * copied or not.
 *
 * @function handleSubmissionLink
 */
const handleSubmissionLink = () => {
    const submissionElt = document.getElementById("submission-title");
    const submissionId = submissionElt.dataset.submissionId;
    const currentUrl = window.location.href;
    const currentPage = currentUrl.lastIndexOf("?");
    if (currentPage === -1) currentPage = currentUrl.length;    // If nemoarchive import url with no params, use whole URL
    const share_url = `${currentUrl.substring(0, currentPage)}?submission_id=${submissionId}`;

    if (copyToClipboard(share_url)) {
        createToast("URL copied to clipboard", "is-info");
    } else {
        createToast(`Failed to copy to clipboard. URL: ${share_url}`, "is-warning" );
    }
};

/**
 * Redirects the user to a gene search page in a new tab.
 *
 * @param {string} layoutShareId - The ID of the layout share.
 * @param {string} [gene] - The gene to search for (optional).
 */
const redirectToGeneSearch = (layoutShareId, gene) => {
    const geneAddition = (gene) ? `&g=${gene}` : "";
    window.open(`${window.origin}/p?l=${layoutShareId}${geneAddition}`, "_blank");
};

/* Poll submission status */
const pollSubmission = async (submissionId) => {
    // Source -> https://javascript.info/long-polling
    const { is_finished: isFinished, datasets } = await getSubmission(submissionId);

    const datasetStatus = {}

    // Update the statuses
    await Promise.allSettled([...Object.keys(datasets)].map( async (datasetId) => {
        const datasetInfo = await getSubmissionDataset(datasets[datasetId].href);
        datasetStatus[datasetId] = datasetInfo.status;

        // Set up the dataset row in the submission table
        updateTblStatus(datasetId, datasetInfo.status)
        printTblMessage(datasetId, datasetInfo.message, "danger")

        if (datasetInfo.status.convert_to_h5ad == "completed") {
            // Populated categorical observations if possible
            populateObsDropdown(datasetId);

            // show permalink (though no curations will be made until make_umap finishes)
            showDatasetPermalink(datasetId);
        }

        const importComplete = isImportComplete(datasetInfo.status);
        if (!importComplete) {
            return;
        }



        if (datasetInfo.status.make_umap == "completed") {
            // After the first dataset is finished, we can View Datasets if desired
            finishElements.importFinish.classList.remove("is-hidden");
        }
        // But if final step did not complete, only a partial success
        finishElements.partialSuccess.classList.remove("is-hidden");
    }));

    const isAllComplete = Object.keys(datasetStatus).every( (datasetId) => isImportComplete(datasetStatus[datasetId]) );

    // If submission is complete (nothing in runnable state), then exit polling.
    if (isAllComplete && finishElements.partialSuccess.classList.contains("is-hidden")) {
        finishElements.completeSuccess.classList.remove("is-hidden");
    }

    //if (isFinished) return;   // This (from submission object) seems to execute early if things go bad, which can lead to confusion
    if (isAllComplete) return;

    // Pull again after a brief timeout
    setTimeout(() => {pollSubmission(submissionId)}, POLL_TIMEOUT);
};

/**
 * Handles the submission creation process.
 *
 * @param {string} jsonUrl - The URL of the JSON data to import.
 */
const createSubmission = async (jsonUrl) => {
    // ! For debugging, the JSON blobs in the GCP bucket have a 30 minute token limit for accessing.
    // ! Currently we are also using local test.json as the nemoarchive API specs are evolving.

    if (!isUserLoggedIn()) {
        createToast(SESSION_ERROR_MESSAGE);
        return;
    }

    try {
        const { data } = await axios.get(jsonUrl);
        //const {data} = await axios.get("nemoarchive_import/test.json");

        const submissionId = grabSubmissionId(data);

        populateSubmissionId(submissionId);
        initializeSubmissionElement(submissionId);
        populateNonSupportedElement();

        const fileEntities = getFileEntities(data);

        // Initialize db information about this submission
        // Check db if datasets were loaded in previous submissions
        // (initialize_new_submission.cgi)
        await initializeNewSubmission(fileEntities, submissionId);

        const emailDiv = document.getElementById("email-on-success-div");
        emailDiv.classList.remove("is-hidden");

        launchSubmissionImport(submissionId, fileEntities);

        // Poll datasets for updates. Function calls itself until importing is finished.
        // If import has an error, it should at least poll once to show final statuses
        await pollSubmission(submissionId);

        emailDiv.disabled = true;
    } catch (error) {
        logErrorInConsole(error);
        createToast(DB_ERROR_MESSAGE);
    }
};

/**
 * Asynchronously views a submission by its parameter, updating the UI with submission details.
 * This function does not resume the submission but checks its status.
 *
 * @param {string} submissionParam - The parameter identifying the submission to view.
 * @returns {Promise<void>} - A promise that resolves when the submission has been fully processed.
 *
 * @throws {Error} - Throws an error if there is an issue fetching the submission datasets.
 */
const viewSubmission = async (submissionParam) => {
    // This is meant purely to check a submission status, but will not resume it (view-only mode)
    const submissionElt = document.getElementById("submission-title");
    submissionElt.dataset.submissionId = submissionParam;
    populateSubmissionId(submissionParam);
    try {
        const { layout_share_id: layoutShareId, collection_name: collectionName, is_submitter: isSubmitter, datasets } = await getSubmission(submissionParam);

        // Prepopulate collection name if previously added
        if (collectionName) {
            document.getElementById("collection-name").value = collectionName;
        }
        // If this is a different user viewing the submission, disable naming the layout
        if (! isSubmitter) {
            document.getElementById("collection-name").disabled = "disabled";
        }

        // Set up the rows
        for (const datasetId in datasets) {
            try {
                const {data} = await axios.get(datasets[datasetId].href);
                // Set up the dataset row in the submission table
                initializeDatasetRow(data);
            } catch (error) {
                logErrorInConsole(error);
                const msg = `Could not fetch submission dataset ${datasetId}.`;
                throw new Error(msg);
            }

        }

        submissionElt.dataset.layout_share_id = layoutShareId; // Store layout_share_id for future retrieval
        submissionElt.addEventListener("click", handleSubmissionLink);

        // Poll datasets for updates. Function calls itself until importing is finished.
        return pollSubmission(submissionParam);
    } catch (error) {
        createToast(`Something went wrong with viewing submission ${submissionParam} - ${error}. Please contact gEAR support.`);
        console.error(error);
    }
};

/* Create UMAPs, save as a layout, and navigate to gene expression results */
const handleViewDatasets = async () => {
    const submissionElt = document.getElementById("submission-title");
    const submissionId = submissionElt.dataset.submissionId;

    const gene = document.getElementById("global-gene-selection").value;

    const plotType = "umap_static";
    const label = "nemoanalytics import default plot";

    // Save new layout for submisison if it does not exist
    const collectionName = document.getElementById("collection-name").value;
    if (collectionName) {
        await updateLayoutName(submissionId, collectionName)
    }

    const tableRows = document.querySelectorAll("#submission-datasets tbody tr");
    await Promise.allSettled([...tableRows].map( async (row) => {
        const datasetId = row.dataset.datasetId;

        if (!(isImportComplete(datasetId))) {
            console.info(`Dataset ${datasetId} did not finish. Cannot save displays or to layout.`)
            return;
        }

        // If select dropdown is not present, then there is no category to plot by.
        const obsLevel = row.querySelector(".js-obslevels-cell");
        const category = (obsLevel.textContent == "Categories not found") ? null : obsLevel.querySelector("select").value;

        // Determine if we need to make or update default displays for this dataset.
        // Either because a display did not previously exist, or a default was created but the user selected a category.
        const defaultDisplayId = await getDefaultDisplay(datasetId)
        if (!defaultDisplayId) {
            return;
        }
        const display = await getDatasetDisplay(defaultDisplayId);
        if (display.label === label && category) {
            const displayConfig = display.plotly_config;
            if (category) displayConfig.colorize_legend_by = category;
            const newLabel = `${label} with category`;
            if (gene) displayConfig.gene_symbol = gene;
            const displayId = await saveNewDisplay(datasetId, displayConfig, plotType, newLabel);
            await saveDefaultDisplay(datasetId, displayId);

            // Update the layout with the new default display for this dataset.
            await updateDisplayInLayout(submissionId, defaultDisplayId, displayId);
        }
        return;
    }));

    // Redirect to layout share ID
    redirectToGeneSearch(submissionElt.dataset.layout_share_id, gene);
};

/**
 * Handles the email button click event.
 * Sends a request to subscribe to email updates for the submission.
 * Notifies the user if the subscription is successful or alerts if it fails.
 *
 * @async
 * @function handleEmailButton
 * @returns {Promise<void>} A promise that resolves when the email subscription request is complete.
 */
const handleEmailButton = async () => {
    const submissionElt = document.getElementById("submission-title");
    const submissionId = submissionElt.dataset.submissionId;

    // Send email when importing has finished. Some RabbitMQ stuff will happen
    const {data} = await axios.put(`/api/submissions/${submissionId}/email`, {
        headers: {"Content-Type": "application/json"},
        });


    if (data.success) {
        createToast("A email will be sent when importing has finished. It is safe to close the tab or browser.", "is-info", true);
        return;
    }
    createToast("Could not subscribe to email updates. Please contact gEAR support.");
    return;

};

/* --- Entry point --- */
/**
 * Handles page-specific login UI updates.
 * @param {Event} event - The event object.
 * @returns {Promise<void>} - A promise that resolves when the UI updates are completed.
 */
const handlePageSpecificLoginUIUpdates = async (event) => {
    // Essentially a "view-only" mode.
    const submissionParam = getUrlParameter("submission_id")
    if (submissionParam) {
        console.info("viewing exising submission");
        return viewSubmission(submissionParam);
    }

    // URI and it's component is encoded so need to decode all that
    // https://thisthat.dev/encode-uri-vs-encode-uri-component/

    // Cannot use decoding from getUrlParameter as the entire url needs to be decoded first
    // Otherwise, the components within the url are dropped
    const encodedJsonUrl = getUrlParameter('url', false);
    if (encodedJsonUrl) {
        const jsonUrl = decodeURIComponent(decodeURI(encodedJsonUrl));
        console.info("staged URL: " + jsonUrl);

        // Create a brand new submission from the JSON pulled in
        return createSubmission(jsonUrl);

        // ? Figure out what to do if counts are per sample rather than per file
    }

    // Neither path above was passed.  Alert user
    createToast("You accessed this page without providing an import URL or an existing submission ID.")

    document.getElementById("view-datasets").addEventListener("click", handleViewDatasets);
    document.getElementById("email-button").addEventListener("click", handleEmailButton);

/* enables 'delete' button - https://bulma.io/documentation/elements/notification/#javascript-example */
    (document.querySelectorAll('.notification .delete') || []).forEach(($delete) => {
        const $notification = $delete.parentNode;
        $delete.addEventListener('click', () => {
            $notification.classList.add("is-hidden");
            const $notificationBody = $notification.querySelector("span");
            $notificationBody.textContent = "";
        });
    });
};

