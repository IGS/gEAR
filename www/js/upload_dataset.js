'use strict';

import { apiCallsMixin, convertToFormData, createToast, getCurrentUser, guid, initCommonUI, openModal } from "./common.v2.js?v=670b2ed";

/* --- constants and variables --- */

let datasetUid = null;
let shareUid = null;
let datasetFormat = null;   // set when user chooses a dataset type
let spatialFormat = null;   // set when user chooses a spatial platform (if applicable)
let performPrimaryAnalysis = true

let processingStatus = null;
const processingStatusCheckInterval = 10; // seconds

const requiredMetadataFields = ['metadata-title', 'metadata-summary', 'metadata-dataset-type',
    'metadata-contact-name', 'metadata-annotation-source', 'metadata-annotation-version',
    'metadata-contact-name', 'metadata-contact-email',
    'metadata-taxon-id', 'metadata-organism'
];

const optionalMetadataFields = ['metadata-contact-institute', 'metadata-platform-id',
    'metadata-instrument', 'metadata-library-selection', 'metadata-library-source',
    'metadata-geo-id', 'metadata-library-strategy', 'metadata-pubmed-id'
];

/* --- Functions and Classes --- */
/**
 * Sends a POST request to add a primary analysis to the current dataset upload.
 * Displays a success toast if the operation is successful, or a warning toast if not.
 *
 * @async
 * @function addPrimaryAnalysisToDataset
 * @returns {Promise<void>} Resolves when the operation is complete and the toast is shown.
 */
const addPrimaryAnalysisToDataset = async () => {
    const {data} = await axios.post('./cgi/add_primary_analysis_to_dataset_upload.cgi', convertToFormData({
        share_uid: shareUid,
        dataset_format: datasetFormat,
        session_id: getCurrentUser()?.session_id,
    }));

    document.getElementById('finalize-migrating-primary-analysis-li').classList.remove("is-hidden");
    if (data.success) {
        createToast('Primary analysis added successfully','is-success');

        // If dataset was not single-cell or spatial, then we cannot have a primary analysis
        if (!data.perform_primary_analysis) {
            performPrimaryAnalysis = false;
            document.getElementById('finalize-migrating-primary-analysis-li').classList.add("is-hidden");
        }
    } else {
        // This is non-fatal, so just show a warning toast
        createToast('Error adding primary analysis to uploaded dataset');
        processingStatus = "error";
    }
}

/**
 * Checks the current processing status of the dataset by making an asynchronous request
 * to the server. Updates the UI elements to reflect the current status, message, and progress.
 * If processing is complete, enables the dataset processing submit button.
 *
 * @async
 * @function checkDatasetProcessingStatus
 * @returns {Promise<void>} Resolves when the status check and UI updates are complete.
 */
const checkDatasetProcessingStatus = async () => {
    const {data} = await axios.post('./cgi/check_dataset_processing_status.cgi', convertToFormData({
        share_uid: shareUid,
        session_id: getCurrentUser()?.session_id
    }));

    processingStatus = data.status;
    document.getElementById('step-process-dataset-status').textContent = processingStatus.charAt(0).toUpperCase() + processingStatus.slice(1);
    document.getElementById('step-process-dataset-status-message').textContent = data.message;
    document.getElementById('dataset-processing-progress').value = data.progress;

    // TODO: Handle the different statuses here
    if (processingStatus === 'complete') {
        await addPrimaryAnalysisToDataset();

        // If still complete after the primary analysis step, enable the next step button
        if (processingStatus === "complete") {
            document.getElementById('dataset-processing-submit').disabled = false;
        } else {
            document.getElementById('step-process-dataset-status').textContent = processingStatus.charAt(0).toUpperCase() + processingStatus.slice(1);
            document.getElementById('step-process-dataset-status-message').textContent = "Error adding primary analysis to dataset";
        }
    }
}

/**
 * Deletes an upload in progress for a given share and dataset.
 *
 * Sends a POST request to the server to delete the specified upload in progress.
 * On success, reloads the list of uploads in progress. On failure, displays a warning toast.
 *
 * @async
 * @function
 * @param {string} shareUid - The unique identifier for the share.
 * @param {string} datasetId - The unique identifier for the dataset.
 * @returns {Promise<void>} Resolves when the operation is complete.
 */
const deleteUploadInProgress = async (shareUid, datasetId) => {
    const {data} = await axios.post('./cgi/delete_upload_in_progress.cgi', convertToFormData({
        share_uid: shareUid,
        dataset_id: datasetId,
        session_id: getCurrentUser()?.session_id
    }));

    if (data.success) {
        loadUploadsInProgress();
    } else {
        createToast('Error deleting upload in progress', data.message, 'is-warning');
    }
}

/**
 * Finalizes the dataset upload process by sending form data to the server,
 * updating UI elements to reflect the status of metadata loading, h5ad migration,
 * and user data migration, and enabling the next step if successful.
 *
 * @async
 * @function finalizeUpload
 * @returns {Promise<void>} Resolves when the upload finalization process is complete.
 */
const finalizeUpload = async () => {

    const payload = {
        dataset_uid: datasetUid,
        share_uid: shareUid,
        session_id: getCurrentUser()?.session_id,
        dataset_format: datasetFormat,
        perform_analysis_migration: performPrimaryAnalysis ? 1 : 0,
        dataset_visibility: document.querySelector('input[name=dataset-visibility]:checked').value
    }

    const data = await apiCallsMixin.finalizeExpressionUpload(payload);

    if (data['metadata_loaded']) {
        document.getElementById('finalize-storing-metadata').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-storing-metadata').classList.add('mdi-checkbox-marked');
    } else {
        document.getElementById('finalize-storing-metadata').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-storing-metadata').classList.add('mdi-skull-scan');
    }

    if (data['h5ad_migrated']) {
        document.getElementById('finalize-migrating-h5ad').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-migrating-h5ad').classList.add('mdi-checkbox-marked');
    } else {
        document.getElementById('finalize-migrating-h5ad').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-migrating-h5ad').classList.add('mdi-skull-scan');
    }

    if (data['userdata_migrated']) {
        document.getElementById('finalize-migrating-userdata').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-migrating-userdata').classList.add('mdi-checkbox-marked');
    } else {
        document.getElementById('finalize-migrating-userdata').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-migrating-userdata').classList.add('mdi-skull-scan');
    }

    if (data['primary_analysis_migrated']) {
        document.getElementById('finalize-migrating-primary-analysis').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-migrating-primary-analysis').classList.add('mdi-checkbox-marked');
    } else {
        document.getElementById('finalize-migrating-primary-analysis').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-migrating-primary-analysis').classList.add('mdi-skull-scan');
    }

    if (data.success) {
        // Not actually doing anything further with access rights after the initial insert,
        // so can set it as successful here
        document.getElementById('finalize-setting-access').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-setting-access').classList.add('mdi-checkbox-marked');
        document.getElementById('dataset-finalize-next-step').disabled = false;
        return;
    }
    msg = data.message || 'Error finalizing dataset upload';

    console.error(`ERROR: ${msg}`);
    document.getElementById('dataset-finalize-status-message').innerText = msg;
    document.getElementById('dataset-finalize-status-message-c').classList.remove('is-hidden');
}

/**
 * Asynchronously populates the metadata form fields with values parsed from an uploaded metadata file.
 * Utilizes the `apiCallsMixin.parseMetadataFile` method to extract metadata from the uploaded file,
 * and fills the corresponding form fields if parsing is successful. Handles both success and error states,
 * updating the UI accordingly (form fields, status messages, and button states).
 *
 * @async
 * @function populateMetadataFormFromFile
 * @returns {Promise<void>} Resolves when the form has been populated or an error message has been displayed.
 */
const populateMetadataFormFromFile = async () => {
    const formData = new FormData(document.getElementById('metadata-upload-form'));
    const data = await apiCallsMixin.parseMetadataFile(formData);
    const button = document.getElementById('metadata-upload-submit');

    // Fill out the form with the data
    if (data.success) {
        document.getElementsByName('metadata-title')[0].value = data.metadata.title.value;
        document.getElementsByName('metadata-summary')[0].value = data.metadata.summary.value;
        document.getElementsByName('metadata-annotation-source')[0].value = data.metadata.annotation_source.value;
        document.getElementsByName('metadata-annotation-version')[0].value = data.metadata.annotation_release_number.value;
        document.getElementsByName('metadata-geo-id')[0].value = data.metadata.geo_accession.value;

        document.getElementsByName('metadata-contact-name')[0].value = data.metadata.contact_name.value;
        document.getElementsByName('metadata-contact-email')[0].value = data.metadata.contact_email.value;
        document.getElementsByName('metadata-contact-institute')[0].value = data.metadata.contact_institute.value;
        document.getElementsByName('metadata-taxon-id')[0].value = data.metadata.sample_taxid.value;
        document.getElementsByName('metadata-organism')[0].value = data.metadata.sample_organism.value;
        document.getElementsByName('metadata-platform-id')[0].value = data.metadata.platform_id.value;
        document.getElementsByName('metadata-instrument')[0].value = data.metadata.instrument_model.value;
        document.getElementsByName('metadata-library-selection')[0].value = data.metadata.library_selection.value;
        document.getElementsByName('metadata-library-source')[0].value = data.metadata.library_source.value;
        document.getElementsByName('metadata-library-strategy')[0].value = data.metadata.library_strategy.value;
        document.getElementsByName('metadata-pubmed-id')[0].value = data.metadata.pubmed_id.value;

        // Handle the metadata-dataset-type select box
        const datasetTypeSelect = document.getElementsByName('metadata-dataset-type')[0];
        for (let i = 0; i < datasetTypeSelect.options.length; i++) {
            if (datasetTypeSelect.options[i].value === data.metadata.dataset_type.value) {
                datasetTypeSelect.selectedIndex = i;
                break;
            }
        }

        document.getElementById('metadata-upload-status-message').textContent = "Form populated with uploaded metadata";
        button.disabled = false;

    } else {
        // Display error message, and clear out the file selected
        document.getElementById('metadata-upload-status-message').textContent = data.message;
        document.getElementById('metadata-file-name').textContent = 'No file selected';
    }

    // change spinner back to submit button
    document.getElementById('metadata-upload-status').classList.remove('is-hidden');
    button.classList.remove('is-loading');
}

/**
 * Asynchronously fetches GEO metadata based on the user-provided GEO ID and populates
 * corresponding form fields with the retrieved data. If no data is found, displays a status message.
 * Also manages the loading state of the GEO lookup button.
 *
 * @async
 * @function getGeoData
 * @returns {Promise<void>} Resolves when the GEO data has been fetched and the form updated.
 */
const getGeoData = async () => {
    const geoId = document.getElementsByName('metadata-geo-id')[0].value;
    const geoData = await apiCallsMixin.fetchGeoData(geoId);

    if (Object.keys(geoData).length === 0) {
        document.getElementById('metadata-geo-lookup-status').classList.remove('is-hidden');
    } else {
        document.getElementById('metadata-geo-lookup-status').classList.add('is-hidden');
        document.getElementsByName('metadata-contact-name')[0].value = geoData.contact_name;
        document.getElementsByName('metadata-contact-email')[0].value = geoData.contact_email;
        document.getElementsByName('metadata-contact-institute')[0].value = geoData.contact_institute;
        document.getElementsByName('metadata-taxon-id')[0].value = geoData.taxid_ch1;
        document.getElementsByName('metadata-organism')[0].value = geoData.organism_ch1;
        document.getElementsByName('metadata-platform-id')[0].value = geoData.platform_id;
        document.getElementsByName('metadata-instrument')[0].value = geoData.instrument_model;
        document.getElementsByName('metadata-library-selection')[0].value = geoData.library_selection;
        document.getElementsByName('metadata-library-source')[0].value = geoData.library_source;
        document.getElementsByName('metadata-library-strategy')[0].value = geoData.library_strategy;
        document.getElementsByName('metadata-pubmed-id')[0].value = geoData.pubmed_id;
    }

    const button = document.getElementById('metadata-geo-lookup');
    button.disabled = false;
    button.classList.remove('is-loading');
}

/**
 * Navigates to a specific step in the dataset upload process, updating the UI to reflect the current step.
 * Handles step marker icons and classes for visual feedback, manages step content visibility,
 * and triggers dataset processing status checks when appropriate.
 *
 * @param {string} step - The label of the step to navigate to. Must be one of:
 *   'enter-metadata', 'upload-dataset', 'process-dataset', 'finalize-dataset', 'curate-dataset'.
 */
const stepTo = (step) => {
    // TODO: switch to using the stepper-fxns.js functions (and unify the two stepper implementations)

    const stepLabels = ['enter-metadata', 'upload-dataset', 'process-dataset',
        'finalize-dataset', 'curate-dataset'
    ];
    let stepReached = false;

    // Walk forward in the steps and handle each
    for (const label of stepLabels) {
        const stepLi = document.getElementById(`step-${label}`);
        const stepMarker = stepLi.firstElementChild;
        const stepIcon = stepMarker.firstElementChild;

        // If the step is the one we want, add the check icon
        if (label === step) {
            stepMarker.classList.add('is-light');
            stepLi.classList.add('is-active');
            stepIcon.firstElementChild.classList.remove('mdi-check-bold');
            stepIcon.firstElementChild.classList.add('mdi-wrench');
            stepReached = true;

        // If not the current step, handle if it's before or after the current
        } else if (stepReached) {
            // These are the steps markers after the current one
            stepLi.classList.remove('is-active');
            stepMarker.classList.remove('is-light');
            stepIcon.firstElementChild.classList.remove('mdi-wrench');
            stepIcon.firstElementChild.classList.remove('mdi-check-bold');
        } else {
            // These are the step markers before the current one
            stepLi.classList.remove('is-active');
            stepMarker.classList.remove('is-light');
            stepIcon.firstElementChild.classList.remove('mdi-wrench');
            stepIcon.firstElementChild.classList.add('mdi-check-bold');
        }
    }

    // if the step is process-dataset, we need to check on the status
    if (step === 'process-dataset') {
        // Check the status immediately, then set an interval to keep doing it.
        checkDatasetProcessingStatus();

        setInterval(() => {
            if (processingStatus !== 'complete' && processingStatus !== 'error') {
                checkDatasetProcessingStatus();
            }
        }, processingStatusCheckInterval * 1000);
    }

    // Inactivate all step contents, then display the one we want
    document.querySelectorAll('.step-c').forEach((item) => {
        item.classList.add('is-hidden');
    });

    document.getElementById('step-' + step + '-c').classList.remove('is-hidden');

    // Scroll to the top of the page
    window.scrollTo(0, 0);
}

/**
 * Loads the list of uploads currently in progress for the current user and updates the UI accordingly.
 *
 * - Fetches uploads in progress from the server using the user's session ID.
 * - Populates the submissions-in-progress table with the retrieved uploads.
 * - Adds event listeners to handle resuming or deleting uploads.
 * - Shows or hides relevant UI sections based on whether uploads are present.
 * - Displays a warning toast if an error occurs during the fetch.
 *
 * @async
 * @function loadUploadsInProgress
 * @returns {Promise<void>} Resolves when the uploads in progress have been loaded and the UI updated.
 */
const loadUploadsInProgress = async () => {
    const {data} = await axios.post('./cgi/get_uploads_in_progress.cgi', convertToFormData({
        session_id: getCurrentUser()?.session_id
    }));

    if (data.success) {
        const tableBody = document.getElementById('submissions-in-progress-table-tbody');

        if (data.uploads.length > 0) {
            const template = document.getElementById('submission-history-row');
            tableBody.innerHTML = '';

            for (const upload of data.uploads) {
                const clone = template.content.cloneNode(true);
                clone.querySelector('tr').dataset.shareId = upload.share_id;
                clone.querySelector('tr').dataset.datasetId = upload.dataset_id;
                clone.querySelector('tr').dataset.loadStep = upload.load_step;
                clone.querySelector('tr').dataset.performPrimaryAnalysis = upload.perform_primary_analysis; // true/false
                clone.querySelector('tr').dataset.datasetIsSpatial = upload.dataset_is_spatial; // true/false

                clone.querySelector('.submission-share-id').textContent = upload.share_id;
                clone.querySelector('.submission-status').textContent = upload.status;
                clone.querySelector('.submission-title').textContent = upload.title;
                clone.querySelector('.submission-dataset-type').textContent = upload.dataset_type;
                tableBody.appendChild(clone);
            };

            // Add click listeners for submissions-in-progress-table-tbody rows we just added
            // First, the resume button
            for (const row of document.querySelectorAll('.submission-history-row .submission-resume')) {
                row.addEventListener('click', (event) => {
                    const row = event.target.closest('tr');

                    shareUid = row.dataset.shareId;
                    const step = row.dataset.loadStep;

                    if (row.dataset.datasetId) {
                        datasetUid = row.dataset.datasetId;
                    }

                    performPrimaryAnalysis = row.dataset.performPrimaryAnalysis;
                    if (performPrimaryAnalysis) {
                        document.getElementById('finalize-migrating-primary-analysis-li').classList.add("is-hidden");
                    }

                    if (row.dataset.datasetIsSpatial) {
                        datasetFormat = 'spatial';
                    }

                    // Do we want to dynamically load the next step or page refresh for it?
                    //  If dynamic we have to reset all the forms.
                    stepTo(step);

                    document.getElementById('submissions-in-progress').classList.add('is-hidden');
                    document.getElementById('submission-c').classList.remove('is-hidden');
                });
            };

            // Now the delete button
            for (const row of document.querySelectorAll('.submission-history-row .submission-delete')) {
                row.addEventListener('click', (event) => {
                    // reset row to be the parent tr element
                    const row = event.target.closest('tr');

                    shareUid = row.dataset.shareId;
                    const dataset_id = row.dataset.datasetId;
                    deleteUploadInProgress(shareUid, dataset_id);
                });
            };

            document.getElementById('submissions-in-progress').classList.remove('is-hidden');
        } else {
            // remove the last row of the table and hide the submissions in progress section
            tableBody.innerHTML = '';
            document.getElementById('submissions-in-progress').classList.add('is-hidden');

            document.getElementById('submission-c').classList.remove('is-hidden');
        }

    } else {
        createToast(`Error loading uploads in progress: ${data.message}`, 'is-warning');
    }
}

/**
 * Converts a field name by removing the 'metadata-' prefix, replacing hyphens with spaces,
 * and capitalizing the first letter.
 *
 * @param {string} field - The field name to prettify.
 * @returns {string} The prettified field name.
 */
const prettifyFieldName = (field) => {
    field = field.replace('metadata-', '');
    field = field.replaceAll('-', ' ');
    return field.charAt(0).toUpperCase() + field.slice(1);
}

/**
 * Asynchronously collects metadata from form fields, sends it to the server for storage,
 * and updates the UI based on the server response.
 *
 * Sends a POST request with dataset metadata to the backend CGI script.
 * On success, advances the UI to the next step and scrolls to the top of the page.
 * On failure, displays an error toast notification.
 *
 * @async
 * @function storeMetadata
 * @returns {Promise<void>} Resolves when the operation is complete.
 */
const storeMetadata = async () => {
    const {data} = await axios.post('./cgi/store_expression_metadata.cgi', convertToFormData({
        dataset_uid: datasetUid,
        share_uid: shareUid,
        session_id: getCurrentUser()?.session_id,
        title: document.getElementsByName('metadata-title')[0].value,
        summary: document.getElementsByName('metadata-summary')[0].value,
        dataset_type: document.getElementsByName('metadata-dataset-type')[0].value,
        annotation_source: document.getElementsByName('metadata-annotation-source')[0].value,
        annotation_version: document.getElementsByName('metadata-annotation-version')[0].value,
        geo_id: document.getElementsByName('metadata-geo-id')[0].value,
        contact_name: document.getElementsByName('metadata-contact-name')[0].value,
        contact_email: document.getElementsByName('metadata-contact-email')[0].value,
        contact_institute: document.getElementsByName('metadata-contact-institute')[0].value,
        taxon_id: document.getElementsByName('metadata-taxon-id')[0].value,
        organism: document.getElementsByName('metadata-organism')[0].value,
        platform_id: document.getElementsByName('metadata-platform-id')[0].value,
        instrument: document.getElementsByName('metadata-instrument')[0].value,
        library_selection: document.getElementsByName('metadata-library-selection')[0].value,
        library_source: document.getElementsByName('metadata-library-source')[0].value,
        library_strategy: document.getElementsByName('metadata-library-strategy')[0].value,
        pubmed_id: document.getElementsByName('metadata-pubmed-id')[0].value
    }));

    if (data.success) {
        // UI for next step:
        /*
        // For the current step:
        <span class="steps-marker">
            <span class="icon">
            <i class="mdi mdi-check-bold"></i>
            </span>
        </span>
        // For the next step:
        <span class="steps-marker is-light">
            <span class="icon">
            <i class="mdi mdi-wrench"></i>
            </span>
        </span>
        */

        stepTo('upload-dataset');

        // Scroll to the top of the page
        window.scrollTo(0, 0);

    } else {
        // Show error with a toast
        createToast('Error saving metadata', data.message, 'is-warning');
    }
}

/**
 * Handles uploading a dataset file to the server using XMLHttpRequest.
 * Collects form data including dataset identifiers, user session, format, and file input,
 * then sends it to the backend endpoint. Updates the UI to reflect upload progress,
 * handles server response, and initiates dataset processing upon successful upload.
 *
 * Side Effects:
 * - Updates progress bar and status messages in the DOM.
 * - Disables/enables submit button and toggles loading/status classes.
 * - Calls `processDataset()` and advances UI steps on success.
 *
 * Dependencies:
 * - Assumes existence of global variables: `datasetUid`, `shareUid`, `datasetFormat`.
 * - Requires DOM elements with IDs: 'dataset-file-input', 'dataset-upload-progress',
 *   'dataset-upload-status-message', 'dataset-upload-submit', 'dataset-upload-status'.
 * - Requires functions: `getCurrentUser()`, `processDataset()`, `stepTo()`.
 * @async
 * @function uploadDataset
 * @returns {Promise<void>} Resolves when the upload process is complete.
 */
const uploadDataset = () => {
    const formData = new FormData();
    formData.append('share_uid', shareUid);
    formData.append('session_id', getCurrentUser()?.session_id);
    formData.append('dataset_format', datasetFormat);
    if (spatialFormat) {
        formData.append('spatial_format', spatialFormat);
    }
    formData.append('dataset_file', document.getElementById('dataset-file-input').files[0]);

    const xhr = new XMLHttpRequest();
    xhr.open('POST', './cgi/store_expression_dataset.cgi', true);

    xhr.upload.onprogress = (event) => {
        if (event.lengthComputable) {
            const percentComplete = (event.loaded / event.total) * 100;
            document.getElementById('dataset-upload-progress').value = percentComplete;
        }
    };

    xhr.onload = () => {
        const response = JSON.parse(xhr.responseText);

        document.getElementById('dataset-upload-status-message').textContent = '';
        document.getElementById('dataset-upload-submit').classList.remove('is-loading');
        document.getElementById('dataset-upload-status').classList.remove('is-hidden');

        if (response.success) {
            document.getElementById('dataset-upload-status-message').textContent = 'Dataset uploaded successfully. Processing beginning momentarily ...';
            document.getElementById('dataset-upload-status').classList.remove('is-hidden');

            processDataset();

            // Wait a few seconds, then move to the next step. The process script
            // (called above) will run for a long time and be monitored separately
            setTimeout(() => {
                stepTo('process-dataset');
            }, 3000);
            return;

        }
        document.getElementById('dataset-upload-status-message').textContent = response.message;
        document.getElementById('dataset-upload-submit').disabled = false;
    };

    xhr.send(formData);
}

/**
 * Processes the uploaded dataset by sending form data to the server.
 * Appends share UID, dataset format, and current user session ID to the form data,
 * then sends a POST request to the server endpoint for processing.
 * Handles the server response, but relies on external status checking for further actions.
 *
 * @async
 * @function
 * @returns {Promise<void>} Resolves when the request is complete.
 */
const processDataset = async () => {
    const formData = new FormData();
    formData.append('share_uid', shareUid);
    formData.append('dataset_format', datasetFormat);
    if (spatialFormat) {
        formData.append('spatial_format', spatialFormat);
    }
    formData.append('session_id', getCurrentUser()?.session_id);

    const xhr = new XMLHttpRequest();
    xhr.open('POST', './cgi/process_uploaded_expression_dataset.cgi', true);

    xhr.onload = () => {
        const response = JSON.parse(xhr.responseText);

        if (response.success) {
            // Nothing really to do here since status checking happening elsewhere
        }
    }

    xhr.send(formData);
}

/**
 * Validates the metadata form by checking required fields for values and enforcing SQL character length limits.
 * Highlights fields with errors by adding the 'is-danger' class.
 *
 * @returns {Object} An object mapping field names to error messages for fields that failed validation.
 */
const validateMetadataForm = () => {
    const erroredFields = {};

    // Each element with a name in required_metadata_fields must have a value
    for (const field of requiredMetadataFields) {
        const element = document.getElementsByName(field)[0];

        if (element.value) {
            element.classList.remove('is-danger');
        } else {
            element.classList.add('is-danger');
            erroredFields[field] = 'Requires a value';
        }
    }

    // Check SQL length limitations
    const fieldCharacterLimits = {
        'metadata-title': 255,
        'metadata-summary': 65535,
        'metadata-annotation-source': 20,
        //'metadata-annotation-version': int,
        'metadata-contact-name': 100,
        'metadata-contact-email': 100,
        'metadata-contact-institute': 255,
        'metadata-platform-id': 255,
        'metadata-instrument': 255,
        'metadata-library-selection': 255,
        'metadata-library-source': 255,
        'metadata-geo-id': 50,
        'metadata-library-strategy': 65535,
        'metadata-pubmed-id': 20
    };

    for (const field in fieldCharacterLimits) {
        const element = document.getElementsByName(field)[0];
        if (element.value.length > fieldCharacterLimits[field]) {
            element.classList.add('is-danger');
            erroredFields[field] = `Exceeds ${fieldCharacterLimits[field]} character limit`;
        } else {
            element.classList.remove('is-danger');
        }
    }

    return erroredFields;
}

/**
 * Initializes the upload dataset page by:
 * - Checking if the user is logged in and displaying the appropriate UI elements.
 * - Loading uploads in progress for logged-in users.
 * - Setting the page header label.
 * - Generating unique identifiers for the dataset and sharing.
 *
 * Assumes the existence of global functions `getCurrentUser()`, `loadUploadsInProgress()`, and `guid()`,
 * as well as global variables `datasetUid` and `shareUid`.
 */
const initPage = () => {
    if (getCurrentUser()?.session_id) {
        document.getElementById('logged-in-c').classList.remove('is-hidden');
        loadUploadsInProgress();
    } else {
        document.getElementById('not-logged-in-c').classList.remove('is-hidden');
    }

    // Set the page title
    document.getElementById('page-header-label').textContent = 'Upload an expression dataset';

    // Generate the UID that will be used for this submission
    datasetUid = guid('long');
    shareUid = guid('short');
}

/* --- Initialization logic --- */

await initCommonUI();
await initPage();

/* --- Event listeners --- */

// Add click listeners for all buttons of class 'format-selector'
const formatSelectorElts = document.getElementsByClassName('format-selector');
for (const btn of formatSelectorElts) {
    btn.addEventListener('click', (event) => {
        // Reset each as selectable
        for (const element of formatSelectorElts) {
            if (element.disabled) {
                continue;
            }
            // set the classList on this button to only be 'mdi' and 'mdi-cancel'
            const icon = element.querySelector('span.icon i');
            icon.classList.remove(...icon.classList);
            icon.classList.add('mdi', 'mdi-checkbox-blank-outline');

            element.querySelector('span.format-status').textContent = 'Choose';
        };

        // Now set things for the one actually clicked
        btn.querySelector('span.icon i').classList.remove('mdi', 'mdi-checkbox-blank-outline');
        btn.querySelector('span.icon i').classList.add('mdi', 'mdi-checkbox-outline');
        btn.querySelector('span.format-status').textContent = 'Selected';
        datasetFormat = btn.dataset.format;

        // If the format is spatial, change text to say "Zarr store"
        const migrateH5adSpan = document.getElementById("finalize-migrating-h5ad-text");
        migrateH5adSpan.textContent = 'Migrating H5AD file';
        if (datasetFormat === 'spatial') {
            migrateH5adSpan.textContent = 'Migrating Zarr store';
        }

    });
};

document.getElementById('new-submission-toggle').addEventListener('click', (event) => {
    event.preventDefault();

    document.getElementById('submissions-in-progress').classList.add('is-hidden');
    document.getElementById('submission-c').classList.remove('is-hidden');
});

document.getElementById('dataset-processing-submit').addEventListener('click', (event) => {
    event.preventDefault();

    stepTo('finalize-dataset');
});

document.getElementById('dataset-curate-submit').addEventListener('click', (event) => {
    event.preventDefault();

    const url = `/dataset_curator.html?dataset_id=${datasetUid}`;
    window.location.href = url;
});

document.getElementById('metadata-form-submit').addEventListener('click', (event) => {
    event.preventDefault();
    const erroredFields = validateMetadataForm();

    if (erroredFields && Object.keys(erroredFields).length === 0) {
        // Form looks good!
        document.getElementById('errored-field-list-c').classList.add('is-hidden');
        storeMetadata();
        return;

    }
    const erroredFieldsUl = document.getElementById('errored-field-list');
    erroredFieldsUl.innerHTML = '';

    // iterate over the errored fields and display them
    for (const field in erroredFields) {
        const fieldLabel = prettifyFieldName(field);
        const fieldMsg = erroredFields[field];

        const li = document.createElement('li');
        li.textContent = `${fieldLabel}: ${fieldMsg}`;
        erroredFieldsUl.appendChild(li);
    };

    document.getElementById('errored-field-list-c').classList.remove('is-hidden');
});

document.getElementById('dataset-file-input').addEventListener('change', (event) => {
    // Was a file selected?
    if (event.currentTarget.files.length > 0) {
        document.getElementById('dataset-upload-submit').disabled = false;
        const file = event.currentTarget.files[0];
        document.getElementById('dataset-file-name').textContent = file.name;
        return;
    }
    document.getElementById('dataset-upload-submit').disabled = true;
    document.getElementById('dataset-file-name').textContent = 'No file selected';
});

document.getElementById('dataset-finalize-submit').addEventListener('click', (event) => {
    event.preventDefault();
    event.currentTargetdisabled = true;
    document.getElementById('finalize-dataset-status-c').classList.remove('is-hidden');

    finalizeUpload();
});

document.getElementById('dataset-finalize-next-step').addEventListener('click', (event) => {
    event.preventDefault();

    // Move to the next step
    stepTo('curate-dataset');
});

document.getElementById('metadata-file-input').addEventListener('change', (event) => {
    const metadataUploadSubmit = document.getElementById('metadata-upload-submit');
    // Was a file selected?
    if (event.currentTarget.files.length > 0) {
        metadataUploadSubmit.disabled = false;
        const file = event.currentTarget.files[0];
        document.getElementById('metadata-file-name').textContent = file.name;
        document.getElementsByName('metadata-dataset-id')[0].value = datasetUid;
        return;
    }
    metadataUploadSubmit.disabled = true;
    document.getElementById('metadata-file-name').textContent = 'No file selected';
});

document.getElementById('dataset-upload-submit').addEventListener('click', (event) => {
    event.preventDefault();
    // make sure they chose a format
    if (!datasetFormat) {
        document.getElementById('dataset-upload-status-message').textContent = 'Please choose a format above first.';
        document.getElementById('dataset-upload-status').classList.remove('is-hidden');
        return;
    }

    if (datasetFormat === "spatial") {
        const spatialFormatSelect = document.getElementById('select-spatial-platform');
        if (spatialFormatSelect.value === '') {
            document.getElementById('dataset-upload-status-message').textContent = 'Please choose a spatial platform above first.';
            document.getElementById('dataset-upload-status').classList.remove('is-hidden');
            return;
        }
        spatialFormat = spatialFormatSelect.value;
    } else {
        spatialFormat = null;   // safeguard
    }

    // change submit button to spinner
    const button = event.currentTarget;
    button.disabled = true;
    button.classList.add('is-loading');
    document.getElementById('dataset-upload-status').classList.add('is-hidden');
    uploadDataset();
});

document.getElementById('metadata-upload-submit').addEventListener('click', (event) => {
    // change submit button to spinner
    event.preventDefault();
    const button = event.currentTarget;
    button.disabled = true;
    button.classList.add('is-loading');
    document.getElementById('metadata-upload-status').classList.add('is-hidden');
    populateMetadataFormFromFile();
});

document.getElementById('metadata-geo-lookup').addEventListener('click', (event) => {
    event.preventDefault();
    const button = event.currentTarget
    button.disabled = true;
    button.classList.add('is-loading');
    const getData = getGeoData();
});

document.getElementById('select-spatial-platform').addEventListener('change', (e) => {
    const platform = e.target.value;
    const reqsSpan = document.getElementById('spatial-requirements');

    if (platform === '') {
        reqsSpan.classList.add('is-hidden');
        document.getElementById("btn-spatial-format-selector").disabled = true;
    } else {
        reqsSpan.classList.remove('is-hidden');
        document.getElementById("btn-spatial-format-selector").disabled = false;
    }
});

// Show modal of spatial requirements if clicked
document.getElementById("spatial-requirements").addEventListener("click", (e) => {
    e.preventDefault();
    const modalId = e.target.dataset.target;
    const modalElt = document.getElementById(modalId);
    openModal(modalElt);

    const platform = document.getElementById("select-spatial-platform").value;
    // Get the text content for the selected option
    const platformText = document.querySelector(`#select-spatial-platform option[value="${platform}"]`).textContent;
    // populate the modal
    document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
        document.getElementById(`${platform}-spatial-reqs`).innerHTML;

    document.querySelector(`#${modalId} .modal-card-title`).textContent = platformText;
});



/*  From Shaun, used to toggle element's stickiness */
// if scrolling makes #summar-s go above top of screen, make it sticky
/*
window.addEventListener("scroll", (event) => {
    const summary = document.querySelector("<<whatever you want here>>");
    summary.classList.remove("stick-to-top");
    if (summary.getBoundingClientRect().top < 0) {
        summary.classList.add("stick-to-top");
    }
});
*/