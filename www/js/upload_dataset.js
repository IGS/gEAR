/*
 
*/

'use strict';

let dataset_uid = null;
let share_uid = null;

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

    document.getElementById('metadata-form-submit').addEventListener('click', (event) => {
        event.preventDefault();
        let errored_fields = validateMetadataForm();

        if (errored_fields.length === 0) {
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

    document.getElementById('metadata-file-input').addEventListener('change', (event) => {
        const file = event.target.files[0];
        document.getElementById('metadata-file-name').textContent = file.name;
        document.getElementsByName('metadata-session-id')[0].value = CURRENT_USER.session_id;
        document.getElementsByName('metadata-dataset-id')[0].value = dataset_uid;
    });

    document.getElementById('metadata-upload-submit').addEventListener('click', (event) => {
        document.getElementById('metadata-upload-form').submit();
    });

    document.getElementById('metadata-geo-lookup').addEventListener('click', (event) => {
        event.preventDefault();
        let button = document.getElementById('metadata-geo-lookup');
        button.disabled = true;
        button.classList.add('is-loading');
        let geo_data = getGeoData();
    });
};

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        document.getElementById('logged-in-c').classList.remove('is-hidden');
    } else {
        document.getElementById('not-logged-in-c').classList.remove('is-hidden');
    }
}

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
        session_id: CURRENT_USER.session_id,
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

        let metadata_step_li = document.getElementById('step-enter-metadata');
        let metadata_step_marker = metadata_step_li.firstElementChild;
        let metadata_step_icon = metadata_step_marker.firstElementChild;

        let upload_step_li = document.getElementById('step-upload-dataset');
        let upload_step_marker = upload_step_li.firstElementChild;
        let upload_step_icon = upload_step_marker.firstElementChild;

        metadata_step_li.classList.remove('is-active');
        upload_step_li.classList.add('is-active');

        metadata_step_marker.classList.remove('is-light');
        metadata_step_icon.firstElementChild.classList.remove('mdi-wrench');
        metadata_step_icon.firstElementChild.classList.add('mdi-check-bold');

        upload_step_marker.classList.add('is-light');
        upload_step_icon.firstElementChild.classList.add('mdi-wrench');

        document.getElementById('step-upload-metadata-c').classList.add('is-hidden');
        document.getElementById('step-upload-dataset-c').classList.remove('is-hidden');

        // Scroll to the top of the page
        window.scrollTo(0, 0);

    } else {
        alert('Failed to store metadata');
    }
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