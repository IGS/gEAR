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

/* Create a dataset permalink URL from the stored share ID */
const createPermalinkUrl = (shareId) => {
    const currentUrl = window.location.href;
    const currentPage = currentUrl.lastIndexOf("nemoarchive_import");
    return currentUrl.substring(0, currentPage) + 'p?s=' + shareId;
}

/* Get submission ID from JSON */
const grabSubmissionId = (jd) => {
    return jd.find(entity => entity.entityType === "file_set").name;
};

/* Get array of all "file" entityType objects */
const getFileEntities = (jd) => {
    return jd.filter(entity => entity.entityType === "file")
}

/* Get array of all "sample" entityType objects */
const getSampleEntities = (jd) => {
    return jd.filter(entity => entity.entityType === "sample")
}

/* Get sample attributes when passed in a sample name. Accepts a list of sample entites */
const getSampleAttributesByName = (jd, name) => {
    return jd.find(entity => entity.name === name).attributes
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
    const foundElt = [...tableHead].filter( elt => elt.textContent === stepName);
    if (foundElt.length) {
        return foundElt[0].cellIndex;
    }
    console.error(`Could not find step ${stepName} in table headers`);
    return -1;
}

/* Convert metadata into a JSON file */
const convertMetadata = (attributes) => {
    const columnIndex = getIndexOfTableStep("Convert Metadata");
    if (stepIsComplete(columnIndex)) return;
    setTblLoadingStatus(attributes.id, columnIndex);
    const isSuccess = async () => {
        try {
            const params = {"dataset_id": attributes.id, "filetype":attributes.filetype, session_id};
            const response = await fetch("/cgi/convert_nemoarchive_metadata_to_json.cgi?", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify(params)
            });
            if (!response.ok) {
                throw new Error(response.statusText);
            }
            const jsonRes = await response.json();
            if (!jsonRes.success) {
                throw new Error(jsonRes.message);
            }
            return true;
        } catch (error) {
            return false;
        }
    }
    // Update status of file in table
    updateTblStatus(datasetId, columnIndex, isSuccess);
}

/* Convert file set to H5AD and write to final destination */
const convertToH5ad = (attributes) => {
    const columnIndex = getIndexOfTableStep("Convert to H5AD");
    if (stepIsComplete(columnIndex)) return;
    setTblLoadingStatus(attributes.id, columnIndex);
    const isSuccess = async () => {
        try {
            const params = new URLSearchParams({"dataset_id": attributes.id, "filetype":attributes.filetype});
            const response = await fetch("/cgi/convert_and_write_h5ad.cgi?" + params);
            if (!response.ok) {
                throw new Error(response.statusText);
            }
            const jsonRes = await response.json();
            if (!jsonRes.success) {
                throw new Error(jsonRes.message);
            }
            return true;
        } catch (error) {
            return false;
        }
    }
    // Update status of file in table
    updateTblStatus(datasetId, columnIndex, isSuccess);
}

const getFileMetadata = async (attributes) => {
    const {dataset_id : datasetId, identifier} = attributes
    if (! identifier) {
        throw new Error(`No identifier found for dataset ${datasetId}`)
    }
    // Call NeMO Archive assets API using identifier
    //const response = await fetch(`https://nemoarchive.org/asset/derived/${identifier}`);
    const response = await fetch('/api/mock_identifier');
    if (!response.ok) {
        throw new Error(response.statusText);
    }
    const jsonRes = await response.json();
    if (!jsonRes.success) {
        throw new Error(jsonRes.message);
    }
}

/* Create new submission-related entries in database */
const initializeNewSubmission = async (fileEntities, submissionId) => {
    const datasets = await processSubmission(fileEntities, submissionId)
    // Still want to show table even with a saving failure to indicate something went awry with at least one dataset
    initializeSubmissionTable(fileEntities, datasets)

}

/* Initialize submission table to show files */
const initializeSubmissionTable = (fileEntities, datasets) => {
    for (const entity of fileEntities) {
        const { attributes, name } = entity;
        const datasetId = attributes.id;

        const {
            pulled_to_vm: pulledToVm
            , convert_metadata: convertMetadata
            , convert_to_h5ad: convertToH5ad
            } = datasets[datasetId].dataset_status

        console.log(datasets[datasetId].dataset_status)

        const template = `
        <tr id="submission-${datasetId}">
            <td class="dataset-id">${datasetId}</td>
            <td>${name}</td>
            ${status2Element[pulledToVm]}
            ${status2Element[convertMetadata]}
            ${status2Element[convertToH5ad]}
            <td id="${datasetId}_obslevels">Not ready</td>
            <td id="${datasetId}_permalink">Not ready</td>
        </tr>
        `
        const parent = document.querySelector("#submission_datasets tbody");
        const htmlCollection = generateElements(template);
        parent.append(htmlCollection);
    }
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
    const parent = document.querySelector(`#${datasetId}_obslevels`);
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
const populateSubmissionId = (jsonData) => {
    const submissionId = grabSubmissionId(jsonData);
    const submissionEl = document.querySelector("#submission_id");
    submissionEl.textContent = submissionId;
    return submissionId;
}

/* Retrieve existing submission or add new one. Returns datasets from submission */
const processSubmission = async (fileEntities, submissionId) => {
    // Going to attempt to get existing submission
    const getResponse = await fetch(`/api/submission/${submissionId}`)
    if (!getResponse.ok) {
        throw new Error(getResponse.statusText);
    }
    const getJsonRes = await getResponse.json();
    if (getJsonRes.success) {
        return getJsonRes.datasets;
    }

    // If existing submission does not exist, create new submission
    const params = { "file_entities": fileEntities, "submission_id": submissionId };
    const postresponse = await fetch("/api/submission", {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify(params)
        });
    if (!postresponse.ok) {
        throw new Error(postresponse.statusText);
    }

    const postJsonRes = await postresponse.json();
    if (!postJsonRes.success) {
        throw new Error(postJsonRes.message);
    }
    return postJsonRes.datasets;
}

/* Pull all files to VM */
const pullFileSetToVm = async (attributes) => {
    const columnIndex = getIndexOfTableStep("Pull to VM");
    if (stepIsComplete(columnIndex)) return;
    setTblLoadingStatus(attributes.id, columnIndex);
    // Pull component files from GCP to VM
    const isSuccess = await pullComponentFilesToVm(attributes);

    // Update status of file in table
    updateTblStatus(attributes.id, columnIndex, isSuccess);
}

/* For a given file set, pull all component files to VM.
   Fail-fast if anything is not successful */
const pullComponentFilesToVm = async (attributes) => {
    try {
        const { component_fields, id } = attributes
        await Promise.all(component_fields.map(async (component) => {
            const path = attributes[component];
            const params = new URLSearchParams({ "bucket_path": path, "dataset_id": id });
            const response = await fetch("/cgi/pull_nemoarchive_gcp_files_to_vm.cgi?" + params);
            if (!response.ok) {
                throw new Error(response.statusText);
            }
            const jsonRes = await response.json();
            if (!jsonRes.success) {
                throw new Error(jsonRes.message);
            }
        }));
        return true;
    } catch (error) {
        return false;
    }
}

const showPermalink = (datasetId) => {
    const permalinkRow = document.querySelector(`${datasetId}_permalink`)
    // get permalink
    const shareId = $(this).attr('value');
    const permalinkUrl =  createPermalinkUrl(shareId);
    const template = `<a href=${permalinkUrl} target="_blank">View Dataset</a>`;
    const htmlCollection = generateElements(template);
    parent.append(htmlCollection);
}

const setTblLoadingStatus = (datasetId, columnIndex) => {
    const tableRow = document.querySelector(`#submission-${datasetId}`);
    // Change status. In case of failure, "cancel" all subsequent tasks
    tableRow.children[columnIndex].outerHTML = tdLoading;
}

/* Update the statuses of the table based on the current job's failure state */
const updateTblStatus = (datasetId, columnIndex, isSuccess) => {
    //TODO: Eventually get this info from the database
    const updatedStatus = isSuccess ? tdCompleted : tdFailed;
    const tableRow = document.querySelector(`#submission-${datasetId}`);
    // Change status. In case of failure, "cancel" all subsequent tasks
    tableRow.children[columnIndex].outerHTML = updatedStatus;
    if (!isSuccess) {
        for (i = columnIndex+1; i < tableRow.children.length; i++) {
            tableRow.children[i].outerHTML = tdCanceled;
        }
    }
}

/* Create UMAPs, save as a layout, and navigate to gene expression results */
const handleViewDatasets = async () => {
    const tableRows = document.querySelectorAll("#submission_datasets tbody tr");
    await Promise.allSettled([...tableRows].map( async (row) => {
        // Create UMAPs
        const datasetId = row.querySelector(".dataset-id").textContent;
        const obsLevel = row.querySelector(`#${datasetId}_obslevels select`).value();
    }));
}

/* Send email when importing has finished. Some RabbitMQ stuff will happen */
const handleEmailButton = () => {

}

window.onload = async () => {
    // ! For debugging, the JSON blobs in the GCP bucket have a 30 minute token limit for accessing.

    if (!CURRENT_USER.id) {
        throw("Must be logged in to use this tool")
    }

    // URI and it's component is encoded so need to decode all that
    // https://thisthat.dev/encode-uri-vs-encode-uri-component/

    // Cannot use decoding from getUrlParameter as the entire url needs to be decoded first
    // Otherwise, the components within the url are dropped
    const encodedJsonUrl = getUrlParameter('url', false);
    const jsonUrl = decodeURIComponent(decodeURI(encodedJsonUrl));
    console.log("staged URL: " + jsonUrl);

    // https://github.github.io/fetch
    //const urlResponse = await fetch(jsonUrl);
    const urlResponse = await fetch("nemoarchive_import/test.json")
    const jsonData = await urlResponse.json();

    const submissionId = await populateSubmissionId(jsonData);
    populateNonSupportedElement();

    const fileEntities = getFileEntities(jsonData);
    const sampleEntities = getSampleEntities(jsonData);


    // Initialize db information about this submission
    // Check db if datasets were loaded in previous submissions
    // (initialize_new_submission.cgi)
    try {
        await initializeNewSubmission(fileEntities, submissionId);
    } catch (error) {
        const mainAlert = document.querySelector("#main_alert");
        mainAlert.textContent = "Something went wrong with saving this submission to database. Please contact gEAR support.";
        mainAlert.classList.remove("is-hidden");
        console.error(error);
    }

    throw("stopping here");

    await Promise.allSettled(fileEntities.map( async (entity) => {
        const { attributes } = entity;
        // Merge in sample attributes
        const sampleAttributes = getSampleAttributesByName(sampleEntities, attributes.sample_id);
        const allAttributes = {...attributes, ...sampleAttributes};
        allAttributes.name = entity.name;

        // Pull files into VM
        await pullFileSetToVm(allAttributes);
        // Create JSON from metadata
        await getFileMetadata(allAttributes);
        // convertMetadata(allAttributes);
        // Convert to H5AD while validating metadata and move to final destination
        // const datasetInfo = convertToH5ad(allAttributes);

        // enable select-obs dropdown
        if (! ["MEX"].includes(allAttributes.filetype)) {
            populateObsDropdown(allAttributes.id)
        }
        // show permalink
        showPermalink(allAttributes.id)
    }));

    // Finalize submission.load_status depending on everything succeeded or not
    // if at least one is successfully imported, show the layout button

    // ? Figure out what to do if counts are per sample rather than per file



    // ! If dataset importer is not dataset owner (submitter in nemoarchive), link to owner
    // ! If owner not created, create owner account where flag is set that entry not created by owner
};

document.getElementById("view_datasets").addEventListener("click", () => handleViewDatasets);
document.getElementById("email_button").addEventListener("click", () => handleEmailButton);
