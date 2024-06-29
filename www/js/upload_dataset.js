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
        let fields_missing = validateMetadataForm();

        if (fields_missing.length === 0) {
            document.getElementById('missing-field-list-c').classList.add('is-hidden');
        } else {
            const missing_fields_ul = document.getElementById('missing-field-list');
            missing_fields_ul.innerHTML = '';
            fields_missing.forEach((field) => {
                // we want to transform each field string to something more readable
                field = field.replace('metadata-', '');
                field = field.replaceAll('-', ' ');
                field = field.charAt(0).toUpperCase() + field.slice(1);

                const li = document.createElement('li');
                li.textContent = field;
                missing_fields_ul.appendChild(li);
            });

            document.getElementById('missing-field-list-c').classList.remove('is-hidden');
        }
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
    console.log(geo_data);

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

    let button = document.getElementById('metadata-geo-lookup');
    button.disabled = false;
    button.classList.remove('is-loading');
}

const validateMetadataForm = () => {
    let fields_missing = [];

    // Each element with a name in required_metadata_fields must have a value
    for (const field of required_metadata_fields) {
        const element = document.getElementsByName(field)[0];

        if (!element.value) {
            element.classList.add('is-danger');
            fields_missing.push(field);
        } else {
            element.classList.remove('is-danger');
        }
    }

    return fields_missing;
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