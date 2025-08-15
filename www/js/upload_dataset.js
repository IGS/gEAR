'use strict';

import { apiCallsMixin, createToast, getCurrentUser, registerPageSpecificLoginUIUpdates } from './common.v2.js';

let dataset_uid = null;
let share_uid = null;
let dataset_format = null;

let processing_status = null;
const processing_status_check_interval = 10; // seconds

let required_metadata_fields = ['metadata-title', 'metadata-summary', 'metadata-dataset-type',
    'metadata-contact-name', 'metadata-annotation-source', 'metadata-annotation-version',
    'metadata-contact-name', 'metadata-contact-email',
    'metadata-taxon-id', 'metadata-organism'
];

let optional_metadata_fields = ['metadata-contact-institute', 'metadata-platform-id',
    'metadata-instrument', 'metadata-library-selection', 'metadata-library-source',
    'metadata-geo-id', 'metadata-library-strategy', 'metadata-pubmed-id'
];

window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Upload an expression dataset';

    // Generate the UID that will be used for this submission
    dataset_uid = guid('long');
    share_uid = guid('short');

    // Add click listeners for all buttons of class 'format-selector'
    document.querySelectorAll('.format-selector').forEach((btn) => {
        btn.addEventListener('click', (event) => {

            // Reset each as selectable
            document.querySelectorAll('.format-selector').forEach((element) => {
                if (! element.disabled) {
                    // set the classList on this button to only be 'mdi' and 'mdi-cancel'
                    let icon = element.querySelector('span.icon i');
                    icon.classList.remove(...icon.classList);
                    icon.classList.add('mdi', 'mdi-checkbox-blank-outline');

                    element.querySelector('span.format-status').textContent = 'Choose';
                }
            });

            // Now set things for the one actually clicked
            btn.querySelector('span.icon i').classList.remove('mdi', 'mdi-checkbox-blank-outline');
            btn.querySelector('span.icon i').classList.add('mdi', 'mdi-checkbox-outline');
            btn.querySelector('span.format-status').textContent = 'Selected';
            dataset_format = btn.dataset.format;
        });
    });

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

        let url = `/dataset_curator.html?dataset_id=${dataset_uid}`;
        window.location.href = url;
    });

    document.getElementById('metadata-form-submit').addEventListener('click', (event) => {
        event.preventDefault();
        let errored_fields = validateMetadataForm();

        if (errored_fields && Object.keys(errored_fields).length === 0) {
            // Form looks good!
            document.getElementById('errored-field-list-c').classList.add('is-hidden');
            storeMetadata();

        } else {
            const errored_fields_ul = document.getElementById('errored-field-list');
            errored_fields_ul.innerHTML = '';

            // iterate over the errored fields and display them
            for (const field in errored_fields) {
                let field_label = prettifyFieldName(field);
                let field_msg = errored_fields[field];

                const li = document.createElement('li');
                li.textContent = `${field_label}: ${field_msg}`;
                errored_fields_ul.appendChild(li);
            };

            document.getElementById('errored-field-list-c').classList.remove('is-hidden');
        }
    });

    document.getElementById('dataset-file-input').addEventListener('change', (event) => {
        // Was a file selected?
        if (event.target.files.length > 0) {
            document.getElementById('dataset-upload-submit').disabled = false;
            const file = event.target.files[0];
            document.getElementById('dataset-file-name').textContent = file.name;
        } else {
            document.getElementById('dataset-upload-submit').disabled = true;
            document.getElementById('dataset-file-name').textContent = 'No file selected';
        }
    });

    document.getElementById('dataset-finalize-submit').addEventListener('click', (event) => {
        event.preventDefault();
        document.getElementById('dataset-finalize-submit').disabled = true;
        document.getElementById('finalize-dataset-status-c').classList.remove('is-hidden');

        finalizeUpload();
    });

    document.getElementById('dataset-finalize-next-step').addEventListener('click', (event) => {
        event.preventDefault();

        // Move to the next step
        stepTo('curate-dataset');
    });

    document.getElementById('metadata-file-input').addEventListener('change', (event) => {
        // Was a file selected?
        if (event.target.files.length > 0) {
            document.getElementById('metadata-upload-submit').disabled = false;
            const file = event.target.files[0];
            document.getElementById('metadata-file-name').textContent = file.name;
            document.getElementsByName('metadata-dataset-id')[0].value = dataset_uid;
        } else {
            document.getElementById('metadata-upload-submit').disabled = true;
            document.getElementById('metadata-file-name').textContent = 'No file selected';
        }
    });

    document.getElementById('dataset-upload-submit').addEventListener('click', (event) => {
        event.preventDefault();

        // make sure they chose a format
        if (!dataset_format) {
            document.getElementById('dataset-upload-status-message').textContent = 'Please choose a format above first.';
            document.getElementById('dataset-upload-status').classList.remove('is-hidden');
            return;
        }

        // change submit button to spinner
        let button = document.getElementById('dataset-upload-submit');
        button.disabled = true;
        button.classList.add('is-loading');
        document.getElementById('dataset-upload-status').classList.add('is-hidden');
        uploadDataset();
    });

    document.getElementById('metadata-upload-submit').addEventListener('click', (event) => {
        // change submit button to spinner
        event.preventDefault();

        let button = document.getElementById('metadata-upload-submit');
        button.disabled = true;
        button.classList.add('is-loading');
        document.getElementById('metadata-upload-status').classList.add('is-hidden');
        populateMetadataFormFromFile();
    });

    document.getElementById('metadata-geo-lookup').addEventListener('click', (event) => {
        event.preventDefault();
        let button = document.getElementById('metadata-geo-lookup');
        button.disabled = true;
        button.classList.add('is-loading');
        let geo_data = getGeoData();
    });
};

const checkDatasetProcessingStatus = async () => {
    const {data} = await axios.post('./cgi/check_dataset_processing_status.cgi', convertToFormData({
        share_uid: share_uid,
        session_id: getCurrentUser().session_id
    }));

    processing_status = data.status;
    document.getElementById('step-process-dataset-status').textContent = processing_status.charAt(0).toUpperCase() + processing_status.slice(1);
    document.getElementById('step-process-dataset-status-message').textContent = data.message;
    document.getElementById('dataset-processing-progress').value = data.progress;

    // TODO: Handle the different statuses here
    if (processing_status === 'complete') {
        document.getElementById('dataset-processing-submit').disabled = false;
    }
}

const deleteUploadInProgress = async (share_uid, dataset_id) => {
    const {data} = await axios.post('./cgi/delete_upload_in_progress.cgi', convertToFormData({
        share_uid: share_uid,
        dataset_id: dataset_id,
        session_id: getCurrentUser().session_id
    }));

    if (data.success) {
        loadUploadsInProgress();
    } else {
        createToast('Error deleting upload in progress', data.message, 'is-warning');
    }
}

const finalizeUpload = async () => {
    let formData = new FormData();
    formData.append('share_uid', share_uid);
    formData.append('session_id', getCurrentUser().session_id);
    formData.append('dataset_uid', dataset_uid);
    formData.append('dataset_format', dataset_format);

    const dataset_visibility = document.querySelector('input[name=dataset-visibility]:checked').value;
    formData.append('dataset_visibility', dataset_visibility);

    const data = await apiCallsMixin.finalizeExpressionUpload(formData);

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

    if (data.success) {
        // Not actually doing anything further with access rights after the initial insert,
        // so can set it as successful here
        document.getElementById('finalize-setting-access').classList.remove('mdi-checkbox-blank-outline');
        document.getElementById('finalize-setting-access').classList.add('mdi-checkbox-marked');
        document.getElementById('dataset-finalize-next-step').disabled = false;
    } else {
        console.log("ERROR");
        console.log(data);
        document.getElementById('dataset-finalize-status-message').innerText = data.message;
        document.getElementById('dataset-finalize-status-message-c').classList.remove('is-hidden');
    }
}

const populateMetadataFormFromFile = async () => {
    const formData = new FormData(document.getElementById('metadata-upload-form'));
    const data = await apiCallsMixin.parseMetadataFile(formData);
    let button = document.getElementById('metadata-upload-submit');

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
        let dataset_type_select = document.getElementsByName('metadata-dataset-type')[0];
        for (let i = 0; i < dataset_type_select.options.length; i++) {
            if (dataset_type_select.options[i].value === data.metadata.dataset_type.value) {
                dataset_type_select.selectedIndex = i;
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

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (getCurrentUser().session_id) {
        document.getElementById('logged-in-c').classList.remove('is-hidden');
        loadUploadsInProgress();
    } else {
        document.getElementById('not-logged-in-c').classList.remove('is-hidden');
    }
}
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

const getGeoData = async () => {
    const geo_id = document.getElementsByName('metadata-geo-id')[0].value;
    const geo_data = await apiCallsMixin.fetchGeoData(geo_id);

    if (Object.keys(geo_data).length === 0) {
        document.getElementById('metadata-geo-lookup-status').classList.remove('is-hidden');
    } else {
        document.getElementById('metadata-geo-lookup-status').classList.add('is-hidden');
        document.getElementsByName('metadata-contact-name')[0].value = geo_data.contact_name;
        document.getElementsByName('metadata-contact-email')[0].value = geo_data.contact_email;
        document.getElementsByName('metadata-contact-institute')[0].value = geo_data.contact_institute;
        document.getElementsByName('metadata-taxon-id')[0].value = geo_data.taxid_ch1;
        document.getElementsByName('metadata-organism')[0].value = geo_data.organism_ch1;
        document.getElementsByName('metadata-platform-id')[0].value = geo_data.platform_id;
        document.getElementsByName('metadata-instrument')[0].value = geo_data.instrument_model;
        document.getElementsByName('metadata-library-selection')[0].value = geo_data.library_selection;
        document.getElementsByName('metadata-library-source')[0].value = geo_data.library_source;
        document.getElementsByName('metadata-library-strategy')[0].value = geo_data.library_strategy;
        document.getElementsByName('metadata-pubmed-id')[0].value = geo_data.pubmed_id;
    }

    let button = document.getElementById('metadata-geo-lookup');
    button.disabled = false;
    button.classList.remove('is-loading');
}

const stepTo = (step) => {
    const step_labels = ['enter-metadata', 'upload-dataset', 'process-dataset',
        'finalize-dataset', 'curate-dataset'
    ];
    let step_reached = false;

    // Walk forward in the steps and handle each
    for (const label of step_labels) {
        let step_li = document.getElementById('step-' + label);
        let step_marker = step_li.firstElementChild;
        let step_icon = step_marker.firstElementChild;

        // If the step is the one we want, add the check icon
        if (label === step) {
            step_marker.classList.add('is-light');
            step_li.classList.add('is-active');
            step_icon.firstElementChild.classList.remove('mdi-check-bold');
            step_icon.firstElementChild.classList.add('mdi-wrench');
            step_reached = true;

        // If not the current step, handle if it's before or after the current
        } else {
            if (step_reached) {
                // These are the steps markers after the current one
                step_li.classList.remove('is-active');
                step_marker.classList.remove('is-light');
                step_icon.firstElementChild.classList.remove('mdi-wrench');
                step_icon.firstElementChild.classList.remove('mdi-check-bold');
            } else {
                // These are the step markers before the current one
                step_li.classList.remove('is-active');
                step_marker.classList.remove('is-light');
                step_icon.firstElementChild.classList.remove('mdi-wrench');
                step_icon.firstElementChild.classList.add('mdi-check-bold');
            }
        }
    }

    // if the step is process-dataset, we need to check on the status
    if (step === 'process-dataset') {
        // Check the status immediately, then set an interval to keep doing it.
        checkDatasetProcessingStatus();

        setInterval(() => {
            if (processing_status !== 'complete' && processing_status !== 'error') {
                checkDatasetProcessingStatus();
            }
        }, processing_status_check_interval * 1000);
    }

    // Inactivate all step contents, then display the one we want
    document.querySelectorAll('.step-c').forEach((item) => {
        item.classList.add('is-hidden');
    });

    document.getElementById('step-' + step + '-c').classList.remove('is-hidden');

    // Scroll to the top of the page
    window.scrollTo(0, 0);
}

const loadUploadsInProgress = async () => {
    const {data} = await axios.post('./cgi/get_uploads_in_progress.cgi', convertToFormData({
        session_id: getCurrentUser().session_id
    }));

    if (data.success) {
        if (data.uploads.length > 0) {
            const template = document.getElementById('submission-history-row');
            document.querySelector('#submissions-in-progress-table-tbody').innerHTML = '';

            data.uploads.forEach((upload) => {
                let clone = template.content.cloneNode(true);
                clone.querySelector('tr').dataset.shareId = upload.share_id;
                clone.querySelector('tr').dataset.datasetId = upload.dataset_id;
                clone.querySelector('tr').dataset.loadStep = upload.load_step;

                clone.querySelector('.submission-share-id').textContent = upload.share_id;
                clone.querySelector('.submission-status').textContent = upload.status;
                clone.querySelector('.submission-title').textContent = upload.title;
                clone.querySelector('.submission-dataset-type').textContent = upload.dataset_type;
                document.querySelector('#submissions-in-progress-table-tbody').appendChild(clone);
            });

            // Add click listeners for submissions-in-progress-table-tbody rows we just added
            // First, the resume button
            document.querySelectorAll('.submission-history-row .submission-resume').forEach((row) => {
                row.addEventListener('click', (event) => {
                    let row = event.target.closest('tr');

                    share_uid = row.dataset.shareId;
                    const step = row.dataset.loadStep;

                    if (row.dataset.datasetId) {
                        dataset_uid = row.dataset.datasetId;
                    }

                    // Do we want to dynamically load the next step or page refresh for it?
                    //  If dynamic we have to reset all the forms.
                    stepTo(step);

                    document.getElementById('submissions-in-progress').classList.add('is-hidden');
                    document.getElementById('submission-c').classList.remove('is-hidden');
                });
            });

            // Now the delete button
            document.querySelectorAll('.submission-history-row .submission-delete').forEach((row) => {
                row.addEventListener('click', (event) => {
                    // reset row to be the parent tr element
                    let row = event.target.closest('tr');

                    share_uid = row.dataset.shareId;
                    const dataset_id = row.dataset.datasetId;
                    deleteUploadInProgress(share_uid, dataset_id);
                });
            });

            document.getElementById('submissions-in-progress').classList.remove('is-hidden');
        } else {
            // remove the last row of the table and hide the submissions in progress section
            document.querySelector('#submissions-in-progress-table-tbody').innerHTML = '';
            document.getElementById('submissions-in-progress').classList.add('is-hidden');

            document.getElementById('submission-c').classList.remove('is-hidden');
        }

    } else {
        createToast('Error loading uploads in progress: ' + data.message, 'is-warning');
    }
}

const prettifyFieldName = (field) => {
    field = field.replace('metadata-', '');
    field = field.replaceAll('-', ' ');
    field = field.charAt(0).toUpperCase() + field.slice(1);
    return field;
}

const storeMetadata = async () => {
    const {data} = await axios.post('./cgi/store_expression_metadata.cgi', convertToFormData({
        dataset_uid: dataset_uid,
        share_uid: share_uid,
        session_id: getCurrentUser().session_id,
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

const uploadDataset = () => {
    const formData = new FormData();
    formData.append('dataset_uid', dataset_uid);
    formData.append('share_uid', share_uid);
    formData.append('session_id', getCurrentUser().session_id);
    formData.append('dataset_format', dataset_format);
    formData.append('dataset_file', document.getElementById('dataset-file-input').files[0]);

    const xhr = new XMLHttpRequest();
    xhr.open('POST', './cgi/store_expression_dataset.cgi', true);

    xhr.upload.onprogress = function(event) {
        if (event.lengthComputable) {
            const percentComplete = (event.loaded / event.total) * 100;
            document.getElementById('dataset-upload-progress').value = percentComplete;
        }
    };

    xhr.onload = function() {
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

        } else {
            document.getElementById('dataset-upload-status-message').textContent = response.message;
            document.getElementById('dataset-upload-submit').disabled = false;
        }
    };

    xhr.send(formData);
}

const processDataset = async () => {
    const formData = new FormData();
    formData.append('share_uid', share_uid);
    formData.append('dataset_format', dataset_format);
    formData.append('session_id', getCurrentUser().session_id);

    const xhr = new XMLHttpRequest();
    xhr.open('POST', './cgi/process_uploaded_expression_dataset.cgi', true);

    xhr.onload = function() {
        const response = JSON.parse(xhr.responseText);

        if (response.success) {
            // Nothing really to do here since status checking happening elsewhere
        }
    }

    xhr.send(formData);
}

const validateMetadataForm = () => {
    let errored_fields = {};

    // Each element with a name in required_metadata_fields must have a value
    for (const field of required_metadata_fields) {
        const element = document.getElementsByName(field)[0];

        if (!element.value) {
            element.classList.add('is-danger');
            errored_fields[field] = 'Requires a value';
        } else {
            element.classList.remove('is-danger');
        }
    }

    // Check SQL length limitations
    const field_character_limits = {
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

    for (const field in field_character_limits) {
        const element = document.getElementsByName(field)[0];
        if (element.value.length > field_character_limits[field]) {
            element.classList.add('is-danger');
            errored_fields[field] = `Exceeds ${field_character_limits[field]} character limit`;
        } else {
            element.classList.remove('is-danger');
        }
    }

    return errored_fields;
}

/**
 * Generates a GUID string
 * @returns {String} The generated GUID
 * @example af8a8416-6e18-a307-bd9c-f2c947bbb3aa
 * @author Slavik Meltser (slavik@meltser.info)
 * @link http://slavik.meltser.info/?p=142
 */
function guid(uid_length) {
    function _p8(s) {
        var p = (Math.random().toString(16)+"000000000").substr(2,8);
        return s ? "-" + p.substr(0,4) + "-" + p.substr(4,4) : p ;
    }
    if (uid_length == 'long') {
        return _p8() + _p8(true) + _p8(true) + _p8();
    }
    if (uid_length == 'short') {
        return _p8();
    }
}


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