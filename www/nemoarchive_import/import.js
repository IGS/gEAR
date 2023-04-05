// STATUS ELEMENTS
const tdPending = '<td class="has-text-white has-background-warning-dark">Pending</td>';
const tdLoading = '<td class="has-text-white has-background-info-dark">Loading</td>';
const tdCompleted = '<td class="has-text-white has-background-success-dark">Complete</td>';
const tdCanceled = '<td class="has-text-white has-background-dark">Canceled</td>';
const tdFailed = '<td class="has-text-white has-background-danger-dark">Failed</td>';

const grabSubmissionId = (jd) => {
    return jd.filter( entity => entity.entityType === "file_set")[0].name;
};

const getFileEntities = (jd) => {
    return jd.filter(entity => entity.entityType == "file")
}

const generateElements = (html) => {
    /* Generate a DocumentFragment based on an HTML template. Returns htmlCollection */
    const template = document.createElement('template');
    template.innerHTML = html.trim();
    return template.content.children[0];
}

const populateFileSubmissionTable = (fileEntities) => {
    // Populate table to show files
    for (entity of fileEntities) {
        const { attributes, name } = entity;
        const datasetId = attributes.id;
        // TODO: handle situations where dataset ID is null (https://developer.mozilla.org/en-US/docs/Web/API/Crypto/randomUUID)
        const parent = document.querySelector("#submission_datasets tbody");
        const htmlCollection = generateElements(`
        <tr id="submission-${datasetId}">
            <td>${datasetId}</td>
            <td>${name}</td>
            ${tdPending}
            ${tdPending}
        </tr>
        `);
        parent.append(htmlCollection);
    }
}

const populateNonSupportedElement = () => {
    const notSupported = getUrlParameter('notsupported') || 0;
    const notSupportedEl = document.querySelector("#not_supported");
    notSupportedEl.textContent = notSupported;
    if (!(parseInt(notSupported))) {
        const notSupportedP = document.querySelector("#not_supported_p");
        notSupportedP.style.display = "none";
    }
}

const populateSubmissionId = (jsonData) => {
    const submissionId = grabSubmissionId(jsonData);
    const submissionEl = document.querySelector("#submission_id");
    submissionEl.textContent = submissionId;
}

const pullFilesToVm = async (fileEntities) => {
    /* Pull all files to VM */
    await Promise.allSettled(fileEntities.map( async (entity) => {
        const { attributes } = entity;
        const datasetId = attributes.id;
        // Pull component files to VM to upload & convert (Joshua's revision)
        const isSuccess = await pullComponentFilesToVm(attributes, attributes.component_fields, datasetId);
        // Update status of file in table
        const columnIndex = 2;
        updateTblStatus(datasetId, columnIndex, isSuccess);
    }));
}

const pullComponentFilesToVm = (attributes, compoonentFields, datasetId) => {
    /* For a given file set, pull all component files to VM. Fail-fast if anything is not successful */
    // TODO: Would love to figure out how to use async/await.catch for this
    return Promise.all(compoonentFields.map( async (component) => {
        const path = attributes[component];
        const params = new URLSearchParams({"bucket_path":path, "dataset_id":datasetId});
            const response = await fetch("/cgi/pull_nemoarchive_gcp_files_to_vm.cgi?" + params)
            if (!response.ok) {
                throw new Error(response.statusText);
            }
            const jsonRes = await response.json();
            if (!jsonRes.success) {
                throw new Error(jsonRes.message);
            }
    })).then((values) => {
        // All jobs passed.
        return true;
    }).catch((error) => {
        // At least one failed
        return false;
    });

}

const updateTblStatus = (datasetId, columnIndex, isSuccess) => {
    /* Update the statuses of the table based on the current job's failure state */
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

window.onload = async () => {
    // ! For debugging, the JSON blobs in the GCP bucket have a 30 minute token limit for accessing.

    // URI and it's component is encoded so need to decode all that
    // https://thisthat.dev/encode-uri-vs-encode-uri-component/

    // Cannot use decoding from getUrlParameter as the entire url needs to be decoded first
    // Otherwise, the components within the url are dropped
    const encodedJsonUrl = getUrlParameter('url', decode=false);
    const jsonUrl = decodeURIComponent(decodeURI(encodedJsonUrl));
    console.log("staged URL: " + jsonUrl);

    // https://github.github.io/fetch
    const urlResponse = await fetch(jsonUrl);
    const jsonData = await urlResponse.json();

    populateSubmissionId(jsonData);
    populateNonSupportedElement();

    // Check db if datasets were loaded in previous submissions


    // Pull files into VM
    const fileEntities = getFileEntities(jsonData);
    populateFileSubmissionTable(fileEntities);
    pullFilesToVm(fileEntities);

    // ? Figure out what to do if counts are per sample rather than per file

    // Convert to H5AD and move to final destination

    // ! If dataset importer is not dataset owner (submitter in nemoarchive), link to owner
    // ! If owner not created, create owner account where flag is set that entry not created by owner

};
