'use strict';

// STATUS ELEMENTS
const tdPending = '<td class="has-text-white has-background-warning-dark">Pending</td>';
const tdLoading = '<td class="has-text-white has-background-info-dark">Loading</td>';
const tdCompleted = '<td class="has-text-white has-background-success-dark">Completed</td>';
const tdCanceled = '<td class="has-text-white has-background-dark">Canceled</td>';
const tdFailed = '<td class="has-text-white has-background-danger-dark">Failed</td>';

const status2Element = {
    pending : tdPending
    , loading : tdLoading
    , completed : tdCompleted
    , canceled : tdCanceled
    , failed : tdFailed
}

const finishElements = {
    importFinish: document.getElementById("import_finish")
    , partialSuccess: document.getElementById("partial_success")
    , completeSuccess: document.getElementById("complete_success")
}

const POLL_TIMEOUT = 10_000; // 10 seconds

const alert = (msg, fade=false) => {
    // TODO: Add some fade effect
    const alertElt = document.getElementById("main_alert");
    const alertBody = alertElt.querySelector("span");
    alertElt.classList.remove("is-hidden");
    alertBody.textContent = msg;
}

const notify = (msg, fade=false) => {
    const notifyElt = document.getElementById("main_notify");
    const notifyBody = notifyElt.querySelector("span");
    notifyElt.classList.remove("is-hidden");
    notifyBody.textContent = msg;
    if (fade) {
        setTimeout(() => {
            notifyElt.classList.add('is-hidden');
            }, 2000);

    }
}

const convertToFormData = (object) => {
    // Source -> https://stackoverflow.com/a/66611630
    // NOTE: When using FormData do not set headers to application/json
    const formData = new FormData();
    Object.keys(object).forEach(key => {
        if (typeof object[key] !== 'object') formData.append(key, object[key])
        else formData.append(key, JSON.stringify(object[key]))
    })
    return formData;
}

const printTblMessage = (datasetId, msg, level) => {
    // ? Should I make this a hover element since messages may be long?
    const cell = document.getElementById(`${datasetId}-messages`);
    cell.classList = [];
    cell.textContent = msg;
    if (level) {
        cell.classList.add(`has-text-${level}-dark`)
    }
}

/* Add this dataset to this layout */
const addDatasetToLayout = async (datasetId, layoutId) => {
    const response = await fetch("./cgi/add_dataset_to_layout.cgi", {
        method: "POST",
        body: convertToFormData({session_id, dataset_id: datasetId, layout_id: layoutId})
        });
    return await response.json();
}

/* Create a dataset permalink URL from the stored share ID */
const createDatasetPermalinkUrl = (shareId) => {
    const currentUrl = window.location.href;
    const currentPage = currentUrl.lastIndexOf("nemoarchive_import");
    return `${currentUrl.substring(0, currentPage)}p?s=${shareId}`;
}

/* Get submission ID from JSON */
const grabSubmissionId = (jd) => {
    return jd.find(entity => entity.entityType === "file_set").name;
};

/* Get array of all "file" entityType objects */
const getFileEntities = (jd) => {
    return jd.filter(entity => entity.entityType === "file");
}

/* Get array of all "sample" entityType objects */
const getSampleEntities = (jd) => {
    return jd.filter(entity => entity.entityType === "sample");
}

/* Get sample attributes when passed in a sample name. Accepts a list of sample entites */
const getSampleAttributesByName = (jd, name) => {
    return jd.find(entity => entity.name === name).attributes;
}

/* Generate a DocumentFragment based on an HTML template. Returns htmlCollection */
const generateElements = (html) => {
    const template = document.createElement('template');
    template.innerHTML = html.trim();
    return template.content.children[0];
}

/* Get specific dataset display information */
const getDatasetDisplay = async (displayId) => {
    const params = new URLSearchParams({"display_id": displayId});

    const response = await fetch(`./cgi/get_dataset_display.cgi?${params}`);
    return await response.json();
}

/* Get a saved default display for the given dataset and user */
const getDefaultDisplay = async (datasetId) => {
    const params = new URLSearchParams({"user_id":CURRENT_USER.id, "dataset_id": datasetId});

    const response = await fetch(`./cgi/get_default_display.cgi?${params}`);
    const jsonRes = await response.json();
    const {default_display_id: defaultDisplayId} = jsonRes;
    return defaultDisplayId;
}

/* Retrieve existing submission */
const getSubmission = async (submissionId) => {
    // Going to attempt to get existing submission
    const response = await fetch(`/api/submissions/${submissionId}`);
    // NOTE: This is actually OK for situations where we want to check the submission exists before creating a new one
    if (!response?.ok) throw new Error(response.statusText);
    const jsonRes = await response.json();
    if (!jsonRes.success) throw new Error(jsonRes.message);

    // Still want to show table even with a saving failure to indicate something went awry with at least one dataset
    return jsonRes;
}

/* Retrieve existing submission_dataset */
const getSubmissionDataset = async (href) => {
    // ? Should I use submission_id and dataset_id as args instead?
    const response = await fetch(href);
    if (!response?.ok) {
        throw new Error(response.statusText);
    }
    const jsonRes = await response.json();
    if (!jsonRes.success) {
        throw new Error(jsonRes.message);
    }
    return jsonRes;
}

/* Create new submission-related entries in database */
const initializeNewSubmission = async (fileEntities, submissionId) => {
    const { datasets } = await processSubmission(fileEntities, submissionId);

    for (const datasetId in datasets) {
        // Reset any non-complete, non-loading steps to pending
        const putResponse = await fetch(`${datasets[datasetId].href}/status`, {
            method: "PUT"
            , headers: {"Content-Type": "application/json"}
            , body: JSON.stringify({"action":"reset_steps"})
        })
        if (!putResponse?.ok) {
            throw new Error(putResponse.statusText);
        }
        // Set up the dataset row in the submission table
        initializeDatasetRow(datasets[datasetId]);

        // Log info about any dataset that could not be initialized
        if (datasets[datasetId].failed_datasets) {
            printTblMessage(datasetId, datasets[datasetId].failed_datasets, "danger")
        }
    }
}

/* Set up the dataset row in the submission table */
const initializeDatasetRow = (dataset) => {
    const {dataset_id:datasetId, share_id:shareId, identifier} = dataset;
    const namespace = identifier.split("nemo:")[1];
    const identifierUrl = `https://assets.nemoarchive.org/${namespace}`;

    const template = `
    <tr id="dataset-${datasetId}" data-share_id="${shareId}" data-dataset_id="${datasetId}">
        <td><a class="has-text-link" href="${identifierUrl}" target="_blank">${identifier}</a></td>
        <td id="${datasetId}-pulled-to-vm" class="has-text-white has-background-warning-dark">Pending</td>
        <td id="${datasetId}-validate-metadata" class="has-text-white has-background-warning-dark">Pending</td>
        <td id="${datasetId}-convert-to-h5ad" class="has-text-white has-background-warning-dark">Pending</td>
        <td id="${datasetId}-make-tsne" class="has-text-white has-background-warning-dark">Pending</td>
        <td id="${datasetId}-messages"></td>
        <td id="${datasetId}-obslevels">Not ready</td>
        <td id="${datasetId}-permalink">Not ready</td>
    </tr>
    `

    // Append row to table
    const parent = document.querySelector("#submission_datasets tbody");
    const htmlCollection = generateElements(template);
    parent.append(htmlCollection);
}

const isImportComplete = (status) => {
    /* Returns True if step finished */
    // i.e. a non-runnable state. Incomplete steps are reset to "pending" when an import is attempted
    if (status.make_tsne == "pending" || status.make_tsne == "loading") return false;
    return true;
}

/* Populate and show number of non-supported files, if any exist */
const populateNonSupportedElement = () => {
    const notSupported = getUrlParameter('notsupported') || 0;
    const notSupportedEl = document.querySelector("#not_supported");
    notSupportedEl.textContent = notSupported;
    if (parseInt(notSupported) > 0) {
        const notSupportedP = document.querySelector("#not_supported_p");
        notSupportedP.style.display = "block";
    }
}

/* Populate select dropdown of observation categorical metadata */
const populateObsDropdown = async (datasetId) => {
    const parent = document.getElementById(`${datasetId}-obslevels`);
    const response = await fetch(`/api/h5ad/${datasetId}`);
    const jsonData = await response.json();
    const obsLevels = Object.keys(jsonData.obs_levels);
    if (!obsLevels.length) {
        const parent = document.getElementById(`${datasetId}-obslevels`);
        parent.textContent = "No categories found.";
        return
    }
    parent.classList.add("is-loading", "select");
    const optionElts = obsLevels.map(cat => `<option>${cat}</option>`);
    const template = `<select>${optionElts.join("")}</select>`;
    const htmlCollection = generateElements(template);
    parent.append(htmlCollection);
    parent.classList.remove("is-loading");

}

/* Populate submission ID in DOM */
const populateSubmissionId = (submissionId) => {
    const submissionEl = document.getElementById("submission_id");
    submissionEl.textContent = submissionId;
}

/* Retrieve existing submission or add new one. Returns datasets from submission */
const processSubmission = async (fileEntities, submissionId) => {
    const submissionElt = document.getElementById("submission_title");

    try {
        // Going to attempt to get existing submission (and resume)
        const getResponse = await fetch(`/api/submissions/${submissionId}`)
        if (!getResponse?.ok) {
            throw new Error(getResponse.statusText);
        }
        const getJsonRes = await getResponse.json();
        if (!getJsonRes.success) {
            throw new Error(getJsonRes.message);
        }

        // Append extra dataset information
        for (const datasetId in getJsonRes.datasets) {
            getJsonRes.datasets[datasetId] = await getSubmissionDataset(getJsonRes.datasets[datasetId].href);
        }

        // Prepopulate collection name if previously added
        if (getJsonRes.collection_name) {
            document.getElementById("collection_name").value = getJsonRes.collection_name;
        }
        // If this is a different user viewing the submission, disable naming the layout
        if (! getJsonRes.is_submitter) {
            document.getElementById("collection_name").disabled = "disabled";
        }

        submissionElt.dataset.layout_share_id = getJsonRes.layout_share_id; // Store layout_share_id for future retrieval

        return getJsonRes

    } catch {
        const isRestricted = 0;    // hardcoded for now

        // If existing submission does not exist, create new submission
        const postParams = {"submission_id": submissionId, "is_restricted": isRestricted };
        const postResponse = await fetch("/api/submissions", {
            method: "POST",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify(postParams)
            });
        if (!postResponse?.ok) {
            throw new Error(postResponse.statusText);
        }

        const postJsonRes = await postResponse.json();
        if (!postJsonRes.success) {
            throw new Error(postJsonRes.message);
        }

        submissionElt.dataset.layout_share_id = getJsonRes.layout_share_id; // Store layout_share_id for future retrieval


        // NOTE: at this point we have no routes for our submission
        // Go through the projected dataset routes, and create the ones that do not already exist
        // If the dataset exists (for another submission), then associate with this submission

        await Promise.allSettled([...fileEntities].map( async (entity) => {
            console.log(entity);
            const datasetId = entity.attributes.id;
            const identifier = entity.attributes.identifier;
            const sdParams = {"dataset_id":datasetId, "identifier":identifier, "is_restricted":isRestricted}

            const href = `/api/submissions/${submissionId}/datasets/${datasetId}`

            // GET route first
            try {
                postJsonRes.datasets[datasetId] = await getSubmissionDataset(href);
            } catch {
                // Submission dataset did not already exist... create it
                const sdPostResponse = await fetch(`/api/submissions/${submissionId}/datasets`, {
                    method:"POST"
                    , headers: {"Content-Type": "application/json"}
                    , body: JSON.stringify(sdParams)
                });
                if (! sdPostResponse?.ok) {
                    postJsonRes.datasets[datasetId].failed_datasets = sdPostResponse.statusText;
                }
                const sdPostJsonRes = await sdPostResponse.json();
                if (sdPostJsonRes.success) {
                    postJsonRes.datasets[datasetId] = sdPostJsonRes
                } else {
                    postJsonRes.datasets[datasetId].failed_datasets = sdPostResponse.message;
                }
            }

            // Save dataset as a new SubmissionMember
            await fetch(`${postJsonRes.datasets[datasetId].href}/members`, {
                method:"PUT"
                , headers: {"Content-Type": "application/json"}
            })
            // TODO: Create a rollback script in case of failure of submission or datasets;
        }));
        return postJsonRes;
    }
}

const saveNewDisplay = async (datasetId, plotConfig, plotType, label) => {
    const response = await fetch("./cgi/save_dataset_display.cgi", {
        method: "POST",
        body: convertToFormData({user_id: CURRENT_USER.id
            , dataset_id: datasetId
            , plotly_config: plotConfig
            , plot_type: plotType
            , label
            })
        });
    const jsonRes = await response.json();
    return jsonRes.display_id;
}

const saveDefaultDisplay = async (datasetId, displayId) => {
    const response = await fetch("./cgi/save_default_display.cgi", {
        method: "POST",
        body: convertToFormData({user_id: CURRENT_USER.id, dataset_id: datasetId, display_id: displayId})
        });
    const jsonRes = await response.json();
    return jsonRes.default_display_id;
}

const updateLayoutName = async (submissionId, collectionName) => {
    const response = await fetch("./cgi/nemoarchive_update_submission_layout_name.cgi", {
        method: "POST",
        body: convertToFormData({submission_id: submissionId, layout_name: collectionName})
        });
    return await response.json();
}

const showDatasetPermalink = (datasetId) => {
    const tableRow = document.getElementById(`dataset-${datasetId}`);
    const permalinkRow = document.getElementById(`${datasetId}-permalink`);
    const shareId = tableRow.dataset.share_id;
    const permalinkUrl =  createDatasetPermalinkUrl(shareId);
    const template = `<a href=${permalinkUrl} target="_blank">View Dataset</a>`;
    permalinkRow.innerHTML = template;
}

/* Update the statuses of the table */
const updateTblStatus = (datasetId, statuses) => {
    // Initially I had the status HTML elements in <div> to be added as innerHTML.
    // However Bulma adds padding to the <td> so the bg-colors look weird.
    // Resolving this by creating a new element, adding the ID, and replacing the old with the new.

    const pulledToVm = document.getElementById(`${datasetId}-pulled-to-vm`);
    const newPulledToVm = generateElements(status2Element[statuses.pulled_to_vm]);
    newPulledToVm.setAttribute("id", `${datasetId}-pulled-to-vm`);
    pulledToVm.replaceWith(newPulledToVm);

    const validateMetadata = document.getElementById(`${datasetId}-validate-metadata`);
    const newValidateMetadata = generateElements(status2Element[statuses.convert_metadata]);
    newValidateMetadata.setAttribute("id", `${datasetId}-validate-metadata`);
    validateMetadata.replaceWith(newValidateMetadata);

    const convertToH5ad = document.getElementById(`${datasetId}-convert-to-h5ad`);
    const newConvertToH5ad = generateElements(status2Element[statuses.convert_to_h5ad]);
    newConvertToH5ad.setAttribute("id", `${datasetId}-convert-to-h5ad`);
    convertToH5ad.replaceWith(newConvertToH5ad);

    const makeTSNE = document.getElementById(`${datasetId}-make-tsne`);
    const newMakeTSNE = generateElements(status2Element[statuses.make_tsne]);
    newMakeTSNE.setAttribute("id", `${datasetId}-make-tsne`);
    makeTSNE.replaceWith(newMakeTSNE);

}

const handleSubmissionLink = () => {
    const submissionElt = document.getElementById("submission_title");
    const submissionId = submissionElt.dataset.submission_id;
    const currentUrl = window.location.href;
    const currentPage = currentUrl.lastIndexOf("?");
    if (currentPage === -1) currentPage = currentUrl.length;    // If nemoarchive import url with no params, use whole URL
    const share_url = `${currentUrl.substring(0, currentPage)}?submission_id=${submissionId}`;

    if (copyToClipboard(share_url)) {
        notify("URL copied to clipboard", true);
    } else {
        notify(`Failed to copy to clipboard. URL: ${share_url}`, true );
    }
}

const redirect_to_gene_search = (layoutShareId, gene) => {
    const gene_addition = (gene) ? `&g=${gene}` : "";
    window.open(`${window.origin}/p?l=${layoutShareId}${gene_addition}`, "_blank");
}

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
        }

        const importSuccess = isImportComplete(datasetInfo.status);
        if (!importSuccess) {
            return;
        }
        // show permalink
        showDatasetPermalink(datasetId);

        if (datasetInfo.status.make_tsne == "completed") {
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
}

/* Create a brand new submission from JSON contents */
const createSubmission = async (jsonUrl) => {
    // ! For debugging, the JSON blobs in the GCP bucket have a 30 minute token limit for accessing.
    // ! Currently we are also using local test.json as the nemoarchive API specs are evolving.

    // Super hacky but sometimes this loads before common.check_for_login completes.
    setTimeout(() => {CURRENT_USER}, 500);
    if (! CURRENT_USER?.id) {
        alert("Must be logged in to import a new submission. Please login and refresh page.");
        throw("Must be logged in to import a new submission");
    }

    // https://github.github.io/fetch
    //const urlResponse = await fetch(jsonUrl);
    const urlResponse = await fetch("nemoarchive_import/test4.json");
    const jsonData = await urlResponse.json();

    const submissionId = await grabSubmissionId(jsonData);
    populateSubmissionId(submissionId);
    const submissionElt = document.getElementById("submission_title");
    submissionElt.dataset.submission_id = submissionId;
    submissionElt.addEventListener("click", () => handleSubmissionLink());
    populateNonSupportedElement();

    const fileEntities = getFileEntities(jsonData);
    const sampleEntities = getSampleEntities(jsonData);

    // Initialize db information about this submission
    // Check db if datasets were loaded in previous submissions
    // (initialize_new_submission.cgi)
    try {
        await initializeNewSubmission(fileEntities, submissionId);
    } catch (error) {
        alert("Something went wrong with saving this submission to database. Please contact gEAR support.");
        console.error(error);
        return;
    }

    const emailDiv = document.getElementById("email_on_success_div")
    emailDiv.classList.remove("is-hidden");

    // Launch off the new submission import. Since we are polling for status, we do not need to block here.
    const params = { "file_metadata": fileEntities, "sample_metadata":sampleEntities, "action":"import"};
    const response = fetch(`/api/submissions/${submissionId}`, {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify(params)
        });

    // Poll datasets for updates. Function calls itself until importing is finished.
    // If import has an error, it should at least poll once to show final statuses
    await pollSubmission(submissionId);

    // Submission is finished and API has responded... email updates won't be triggered at this point
    emailDiv.classList.add("is-hidden");

    return;
}

/* View a submission in "view-only" mode */
const viewSubmission = async (submissionParam) => {
    // This is meant purely to check a submission status, but will not resume it
    const submissionElt = document.getElementById("submission_title");
    submissionElt.dataset.submission_id = submissionParam;
    populateSubmissionId(submissionParam);
    const { layout_share_id: layoutShareId, collection_name: collectionName, is_submitter: isSubmitter, datasets } = await getSubmission(submissionParam);

    // Prepopulate collection name if previously added
    if (collectionName) {
        document.getElementById("collection_name").value = collectionName;
    }
    // If this is a different user viewing the submission, disable naming the layout
    if (! isSubmitter) {
        document.getElementById("collection_name").disabled = "disabled";
    }

    // Set up the rows
    for (const datasetId in datasets) {
        const response = await fetch(datasets[datasetId].href);
        if (!response?.ok) {
            throw new Error(response.statusText);
        }
        const datasetInfo = await response.json();

        // Set up the dataset row in the submission table
        initializeDatasetRow(datasetInfo)
    }

    submissionElt.dataset.layout_share_id = layoutShareId; // Store layout_share_id for future retrieval
    submissionElt.addEventListener("click", handleSubmissionLink);

    // Poll datasets for updates. Function calls itself until importing is finished.
    return pollSubmission(submissionParam);
}

/* Create UMAPs, save as a layout, and navigate to gene expression results */
const handleViewDatasets = async () => {
    const submissionElt = document.getElementById("submission_title");
    const submissionId = submissionElt.dataset.submission_id;

    const gene = document.getElementById("global_gene_selection").value;

    const plotType = "tsne_static";
    const label = "nemoanalytics import default plot";

    // Save new layout for submisison if it does not exist
    const collectionName = document.getElementById("collection_name").value;
    if (collectionName) {
        await updateLayoutName(submissionId, collectionName)
    }

    const tableRows = document.querySelectorAll("#submission_datasets tbody tr");
    await Promise.allSettled([...tableRows].map( async (row) => {
        const datasetId = row.dataset.dataset_id;

        if (!(isImportComplete(datasetId))) {
            console.info(`Dataset ${datasetId} did not finish. Cannot save displays or to layout.`)
            return;
        }

        // If select dropdown is not present, then there is no category to plot by.
        const obsLevel = document.getElementById(`${datasetId}-obslevels`);
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
        }
        return;
    }));

    // Redirect to layout share ID
    redirect_to_gene_search(submissionElt.dataset.layout_share_id, gene);
}

/* Send email when importing has finished. Some RabbitMQ stuff will happen */
const handleEmailButton = async () => {
    const submissionElt = document.getElementById("submission_title");
    const submissionId = submissionElt.dataset.submission_id;

    const response = await fetch(`/api/submissions/${submissionId}/email`, {
        method: "PUT",
        headers: {"Content-Type": "application/json"},
        });

    const jsonRes = await response.json();

    if (jsonRes.success) {
        notify("A email will be sent when importing has finished. It is safe to close the tab or browser.");
        return;
    }
    alert("Could not subscribe to email updates. Please contact gEAR support.");
    return

}

window.onload = () => {

    // Essentially a "view-only" mode.
    const submissionParam = getUrlParameter("submission_id")
    if (submissionParam) {
        return viewSubmission(submissionParam);
    }

    // URI and it's component is encoded so need to decode all that
    // https://thisthat.dev/encode-uri-vs-encode-uri-component/

    // Cannot use decoding from getUrlParameter as the entire url needs to be decoded first
    // Otherwise, the components within the url are dropped
    const encodedJsonUrl = getUrlParameter('url', false);
    const jsonUrl = decodeURIComponent(decodeURI(encodedJsonUrl));
    console.log("staged URL: " + jsonUrl);

    // Create a brand new submission from the JSON pulled in
    createSubmission(jsonUrl);

    // ? Figure out what to do if counts are per sample rather than per file
};

document.getElementById("view_datasets").addEventListener("click", handleViewDatasets);
document.getElementById("email_button").addEventListener("click", handleEmailButton);

/* enables 'delete' button - https://bulma.io/documentation/elements/notification/#javascript-example */
document.addEventListener('DOMContentLoaded', () => {
    (document.querySelectorAll('.notification .delete') || []).forEach(($delete) => {
        const $notification = $delete.parentNode;
        $delete.addEventListener('click', () => {
            $notification.classList.add("is-hidden");
            const $notificationBody = $notification.querySelector("span");
            $notificationBody.textContent = "";
        });
    });
});

