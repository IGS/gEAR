'use strict';

// STATUS ELEMENTS
const tdPending = '<td class="has-text-white has-background-warning-dark">Pending</td>';
const tdLoading = '<td class="has-text-white has-background-info-dark">Loading</td>';
const tdCompleted = '<td class="has-text-white has-background-success-dark">Completed</td>';
const tdCanceled = '<td class="has-text-white has-background-dark">Canceled</td>';
const tdFailed = '<td class="has-text-white has-background-danger-dark">Failed</td>';

const status2Element = {
    "pending" : tdPending
    , "loading" : tdLoading
    , "completed" : tdCompleted
    , "canceled" : tdCanceled
    , "failed" : tdFailed
}

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
        notifyElt.classList.replace('show', 'hide');
        notifyElt.classList.replace('hide', 'is-hidden');
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

/* Given a table header "step" name, return the cell index position */
const getIndexOfTableStep = (stepName) => {
    const tableHead = document.querySelectorAll("#submission_datasets thead tr th");
    const foundElt = [...tableHead].filter(elt => elt.textContent === stepName);
    if (foundElt.length) {
        return foundElt[0].cellIndex;
    }
    console.error(`Could not find step ${stepName} in table headers`);
    return -1;
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

const getFileMetadata = async (attributes) => {
    const {id : datasetId, identifier} = attributes
    if (! identifier) {
        throw new Error(`No identifier found for dataset ${datasetId}`)
    }
    // Call NeMO Archive assets API using identifier
    //const response = await fetch(`https://nemoarchive.org/asset/derived/${identifier}`);
    const response = await fetch(`/api/mock_identifier/${identifier}`);
    if (!response?.ok) {
        throw new Error(response.statusText);
    }
    const jsonRes = await response.json();
    if (!jsonRes.success) {
        printTblMessage(datasetId, jsonRes.message, "danger");
        throw new Error(jsonRes.message);
    }
    return jsonRes.metadata
}

/* Retrieve existing submission */
const getSubmission = async(submissionId) => {
    // Going to attempt to get existing submission
    const response = await fetch(`/api/submissions/${submissionId}`)
    if (!response?.ok) {
        throw new Error(response.statusText);
    }
    const jsonRes = await response.json();
    if (!jsonRes.success) {
        throw new Error(jsonRes.message);
    }

    for (const datasetId in jsonRes.datasets) {
        const response = await fetch(datasets[datasetId].href);
        const datasetInfo = await response.json();

        // Set up the dataset row in the submission table
        intitializeDatasetRow(datasetInfo)
    }

    // Still want to show table even with a saving failure to indicate something went awry with at least one dataset
    return jsonRes;
}

/* Create new submission-related entries in database */
const initializeNewSubmission = async (fileEntities, submissionId) => {
    const { datasets } = await processSubmission(fileEntities, submissionId);
    for (const datasetId in datasets) {
        // Link submission and dataset together if they have not been
        await fetch(datasets[datasetId].href, {
            method:"PUT"
            , body:JSON.stringify({"action":"save_submission_member"})
        })

        const response = await fetch(datasets[datasetId].href);
        const datasetInfo = await response.json();

        // Reset any non-complete, non-loading steps to pending
        const putResponse = await fetch(datasets[datasetId].href, {
            method: "PUT"
            , body: json.stringify({"action":"reset_steps"})
        })
        const {newStatus} = await response.json();
        datasetInfo.status = newStatus

        // Set up the dataset row in the submission table
        updateDatasetRow(datasetInfo)
    }
}

/* Set up the dataset row in the submission table */
const updateDatasetRow = (dataset) => {
    const {dataset_id:datasetId, identifier, share_id:shareId} = dataset;
    const namespace = identifier.split("nemo:")[1];
    const identifierUrl = `https://assets.nemoarchive.org/${namespace}`;

    const {
        pulled_to_vm: pulledToVm
        , convert_metadata: validateMetadata
        , convert_to_h5ad: convertToH5ad
        } = datasets[datasetId].status;

    const template = `
    <tr id="dataset-${datasetId}" data-share_id="${shareId}" data-dataset_id="${datasetId}">
        <td><a class="has-text-link" href="${identifierUrl}" target="_blank">${identifier}</a></td>
        ${status2Element[pulledToVm]}
        ${status2Element[validateMetadata]}
        ${status2Element[convertToH5ad]}
        <td id="${datasetId}-messages"></td>
        <td id="${datasetId}-obslevels">Not ready</td>
        <td id="${datasetId}-permalink">Not ready</td>
    </tr>
    `
    const parent = document.querySelector("#submission_datasets tbody");
    const htmlCollection = generateElements(template);
    parent.append(htmlCollection);
}

const isImportComplete = (datasetId) => {
    /* Returns True if step finished */
    const columnIndex = getIndexOfTableStep("Convert to H5AD");
    if (!shouldStepRun(datasetId, columnIndex)) return true;
    return false
    // ? Consolidate with first steps of "ConvertH5ad"
}

/* Convert metadata into a JSON file */
const makeDefaultDisplay = async (datasetId, category, gene) => {
    const params = {"dataset_id": datasetId, "category":category, "gene":gene, "session_id":session_id, "action":'make_display'}
    const response = await fetch(`/api/submission_dataset/${datasetId}`, {
            method: "POST",
            body: json.stringify(params)
            });
    const jsonRes = await response.json();
    return jsonRes.plot_config;
}

/* Populate and show number of non-supported files, if any exist */
const populateNonSupportedElement = () => {
    const notSupported = getUrlParameter('notsupported') || 0;
    const notSupportedEl = document.querySelector("#not_supported");
    notSupportedEl.textContent = notSupported;
    if (parseInt(notSupported) === 0) {
        const notSupportedP = document.querySelector("#not_supported_p");
        notSupportedP.style.display = "none";
    }
}

/* Populate select dropdown of observation categorical metadata */
const populateObsDropdown = async (datasetId) => {
    const parent = document.getElementById(`${datasetId}-obslevels`);
    parent.classList.add("is-loading", "select");
    const response = await fetch(`/api/h5ad/${datasetId}`);
    const jsonData = await response.json();
    const obsLevels = Object.keys(jsonData.obs_levels);

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
    // Going to attempt to get existing submission (and resume)
    const getResponse = await fetch(`/api/submissions/${submissionId}?${getParams}`)
    if (!getResponse?.ok) {
        throw new Error(getResponse.statusText);
    }
    const getJsonRes = await getResponse.json();
    if (getJsonRes.success) {
        return getJsonRes;
    }

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
        printTblMessage(datasetId, jsonRes.message, "danger");
        throw new Error(postJsonRes.message);
    }

    // NOTE: at this point we have no routes for our submission
    // Go through the projected dataset routes, and create the ones that do not already exist

    for (const entity in fileEntities) {
        const datasetId = entity.attributes.id;
        const identifier = entity.attributes.identifier;
        const sdParams = {"dataset_id":datasetId, "identifier":identifier, "is_restricted":isRestricted}

        // GET route first
        const sdGetResponse = await fetch(`/api/submissions/${submissionId}/datasets/${datasetId}`);
        if (!sdGetResponse?.ok) {
            // Submission dataset did not already exist... create it
            const sdPostResponse = await fetch(`/api/submissions/${submissionId}/datasets`, {
                method:"POST"
                , body: JSON.stringify(sdParams)
            });
            const sdPostJsonRes = await sdPostResponse.json();
            if (sdPostJsonRes.success) {
                postJsonRes.datasets[datasetId].href = sdPostJsonRes.href
            }
        } else {
            const sdGetJsonRes = await sdGetResponse.json();
            if (sdGetJsonRes.success) {
                postJsonRes.datasets[datasetId].href = sdGetJsonRes.href
            }
        }
    }

    return postJsonRes;
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

const saveNewLayout = async (submissionId, collectionName) => {
    const response = await fetch("./cgi/nemoarchive_add_submission_layout.cgi", {
        method: "POST",
        body: convertToFormData({session_id, submission_id: submissionId, layout_name: collectionName})
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

const setTblLoadingStatus = (datasetId, columnIndex) => {
    const tableRow = document.getElementById(`dataset-${datasetId}`);
    tableRow.children[columnIndex].outerHTML = tdLoading;
}

/* Return true if this step should run (is pending), false otherwise */
const shouldStepRun = (datasetId, columnIndex) => {
    const tableRow = document.getElementById(`dataset-${datasetId}`);
    return (tableRow.children[columnIndex].outerHTML === tdPending);
}

/* Update the statuses of the table based on the current job's failure state */
const updateTblStatus = (datasetId, columnIndex, isSuccess) => {
    //TODO: Eventually get this info from the database
    const updatedStatus = isSuccess ? tdCompleted : tdFailed;
    const tableRow = document.getElementById(`dataset-${datasetId}`);
    // Change status. In case of failure, "cancel" all subsequent tasks
    tableRow.children[columnIndex].outerHTML = updatedStatus;
    if (!isSuccess) {
        for (i = columnIndex+1; i < tableRow.children.length; i++) {
            tableRow.children[i].outerHTML = tdCanceled;
        }
    }
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

/* Create a brand new submission from JSON contents */
const createSubmission = async (jsonUrl, importFinish, partialSuccess) => {
    // https://github.github.io/fetch
    //const urlResponse = await fetch(jsonUrl);
    const urlResponse = await fetch("nemoarchive_import/test.json");
    const jsonData = await urlResponse.json();

    const submissionId = await grabSubmissionId(jsonData);
    populateSubmissionId(submissionId);
    const submissionElt = document.getElementById("submission_title");
    submissionElt.dataset.submission_id = submissionId;
    submissionElt.addEventListener("click", () => handleSubmissionLink());
    populateNonSupportedElement();

    const postParams = { "metadata": jsonData, "action":"import" };
    const postResponse = await fetch(`/api/submissions/${submissionId}`, {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify(postParams)
        });
    if (!postResponse?.ok) {
        throw new Error(postResponse.statusText);
    }

    const postJsonRes = await postResponse.json();
    if (!postJsonRes.success) {
        printTblMessage(datasetId, jsonRes.message, "danger");
        throw new Error(postJsonRes.message);
    }

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
    }

    await Promise.allSettled(fileEntities.map(async (entity) => {
        const { attributes: fileAttributes } = entity;
        // Merge in sample attributes
        const sampleAttributes = getSampleAttributesByName(sampleEntities, fileAttributes.sample_id);
        const allAttributes = { dataset: fileAttributes, sample: sampleAttributes };
        allAttributes.name = entity.name;

        // Create JSON from metadata
        const fileMetadata = await getFileMetadata(allAttributes.dataset);

        // Merge all API metadata and metadata pulled from the portal
        // TODO: Pull sample just for specific dataset
        const allMetadata = {
            dataset: { ...allAttributes.dataset, ...fileMetadata.dataset },
            sample: { ...allAttributes.sample, ...fileMetadata.sample }
        };

        const datasetId = allMetadata.dataset.id;

        // If import has already finished, let's just focus on the display
        if (isImportComplete(datasetId) && await getDefaultDisplay(datasetId)) {
            return true;
        }

        const params = { "metadata": allMetadata};
        const response = await fetch(`/api/submission_dataset/${datasetId}`, {
            method: "POST",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify(params)
            });
        if (!response?.ok) {
            throw new Error(postResponse.statusText);
        }
        const jsonRes = await response.json();
        if (!jsonRes.success) {
            printTblMessage(datasetId, jsonRes.message, "danger");
            partialSuccess.classList.remove("is-hidden");
            return false
        }

        // Save the newly created display as a default.
        const plotType = "tsne_static";
        const label = "nemoanalytics import default plot";
        const displayId = await saveNewDisplay(datasetId, jsonRes.plot_config, plotType, label);
        await saveDefaultDisplay(datasetId, displayId);

        // enable select-obs dropdown
        if (["MEX"].includes(allMetadata.dataset.filetype)) {
            const parent = document.getElementById(`${datasetId}-obslevels`);
            parent.textContent = "Categories not found";
        } else {
            populateObsDropdown(datasetId);
        }
        // show permalink
        showDatasetPermalink(datasetId);

        // After the first dataset is finished, we can View Datasets if desired
        importFinish.classList.remove("is-hidden");
        return true
    }));
}

/* View a submission in "view-only" mode */
const viewSubmission = async (submissionParam, importFinish) => {
    // This is meant purely to check a submission status, but will not resume it
    const submissionElt = document.getElementById("submission_title");
    submissionElt.dataset.submission_id = submissionParam;
    populateSubmissionId(submissionParam);
    const { layout_share_id: layoutShareId, datasets } = await getSubmission(submissionParam);
    submissionElt.dataset.layout_share_id = layoutShareId; // Store layout_share_id for future retrieval
    submissionElt.addEventListener("click", handleSubmissionLink);

    for (const datasetId in datasets) {
        const importSuccess = isImportComplete(datasetId);
        if (importSuccess) {
            //TODO: Hide obsLevels column
            // show permalink
            showDatasetPermalink(datasetId);
        }
    }

    importFinish.classList.remove("is-hidden");
    return;
}

/* Create UMAPs, save as a layout, and navigate to gene expression results */
const handleViewDatasets = async () => {
    const submissionElt = document.getElementById("submission_title");
    const submissionId = submissionElt.dataset.submission_id;

    const gene = document.getElementById("global_gene_selection").value;

    // If submission already has an associated layout, we do not want to change it.
    // Just redirect to the gene search
    if (submissionElt.dataset.layout_share_id) {
        redirect_to_gene_search(submissionElt.dataset.layout_share_id, gene);
        return;
    }

    const plotType = "tsne_static";
    const label = "nemoanalytics import default plot";

    // Save new layout for submisison if it does not exist
    const collectionName = document.getElementById("collection_name").value;
    const tableRows = document.querySelectorAll("#submission_datasets tbody tr");
    try {
        const {layout_id: layoutId, layout_share_id: layoutShareId} = await saveNewLayout(submissionId, collectionName)
        await Promise.allSettled([...tableRows].map( async (row) => {
            const datasetId = row.dataset.dataset_id;

            if (!(isImportComplete(datasetId))) {
                console.info(`Dataset ${datasetId} did not finish. Cannot save displays or to layout.`)
                return;
            }

            // Add completed datasets to layout (easy since they are in submission)
            await addDatasetToLayout(datasetId, layoutId)

            // If select dropdown is not present, then there is no category to plot by.
            const obsLevel = document.getElementById(`${datasetId}-obslevels`);
            const category = (obsLevel.textContent == "Categories not found") ? null : obsLevel.querySelector("select").value;

            // Determine if we need to make or update default displays for this dataset.
            // Either because a display did not previously exist, or a default was created but the user selected a category.
            const defaultDisplayId = await getDefaultDisplay(datasetId)
            if (defaultDisplayId) {
                const display = await getDatasetDisplay(defaultDisplayId);
                if (display.label === label && category) {
                    const displayConfig = display.plotly_config;
                    displayConfig.colorize_legend_by = category;
                    const newLabel = label + " with category";
                    if (gene) displayConfig.gene_symbol = gene;
                    const displayId = await saveNewDisplay(datasetId, displayConfig, plotType, newLabel);
                    await saveDefaultDisplay(datasetId, displayId);
                }
                return;
            }

            // If user chooses a category, make a new tSNE_static display
            if (!category) {
                return;
            }
            // No default display was found. Need to make one.
            // ! This step can take a while. Probably should add a loading icon and notify user.
            const plotConfig = await makeDefaultDisplay(datasetId, category, gene);
            const displayId = await saveNewDisplay(datasetId, plotConfig, plotType, label);
            await saveDefaultDisplay(datasetId, displayId);
        }));

        // Redirect to layout share ID
        redirect_to_gene_search(layoutShareId, gene);

    } catch (e) {
        console.warn(`Submission ${submissionId} already has layout ID. Not saving new datasets.`)
        // ? Should we return the layout ID and redirect anyways?
    }
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
    return

}

window.onload = () => {
    // ! For debugging, the JSON blobs in the GCP bucket have a 30 minute token limit for accessing.
    // ! Currently we are also using local test.json as the nemoarchive API specs are evolving.

    if (! CURRENT_USER?.id) {
        throw("Must be logged in to use this tool")
    }

    const importFinish = document.getElementById("import_finish");
    const partialSuccess = document.getElementById("partial_success");
    const completeSuccess = document.getElementById("complete_success");

    // Essentially a "view-only" mode.
    const submissionParam = getUrlParameter("submission_id")
    if (submissionParam) {
        return viewSubmission(submissionParam, importFinish);
    }


    // URI and it's component is encoded so need to decode all that
    // https://thisthat.dev/encode-uri-vs-encode-uri-component/

    // Cannot use decoding from getUrlParameter as the entire url needs to be decoded first
    // Otherwise, the components within the url are dropped
    const encodedJsonUrl = getUrlParameter('url', false);
    const jsonUrl = decodeURIComponent(decodeURI(encodedJsonUrl));
    console.log("staged URL: " + jsonUrl);

    // Create a brand new submission from the JSON pulled in
    createSubmission(jsonUrl, importFinish, partialSuccess);

    // If everything went well, change to a complete success
    if (partialSuccess.classList.contains("is-hidden")) {
        completeSuccess.classList.remove("is-hidden");
    }

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

