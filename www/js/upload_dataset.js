/*
 
*/

'use strict';

let dataset_uid = null;
let share_uid = null;

window.onload=function() {
    // Set the page title
    document.getElementById('page-header-label').textContent = 'Upload an expression dataset';

    // Generate the UID that will be used for this submission
    dataset_uid = guid('long');
    share_uid = guid('short');

    

};

const handlePageSpecificLoginUIUpdates = async (event) => {
    if (CURRENT_USER.session_id) {
        document.getElementById('logged-in-c').classList.remove('is-hidden');
    } else {
        document.getElementById('not-logged-in-c').classList.remove('is-hidden');
    }
}

/**
 * Generates a GUID string.
 * @returns {String} The generated GUID.
 * @example af8a8416-6e18-a307-bd9c-f2c947bbb3aa
 * @author Slavik Meltser (slavik@meltser.info).
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