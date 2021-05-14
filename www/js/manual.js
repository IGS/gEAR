window.onload=function() {
    // check if the user is already logged in
    check_for_login();
    session_id = Cookies.get('gear_session_id');

    // This matches the choices in #list-tab
    //  'basic' corresponds to '#list-basic-choice', 'curation' to '#list-curation-choice', etc.
    var doc_link = getUrlParameter('doc');
    if (doc_link) {
        $("#list-" + doc_link + "-choice").trigger('click');
    }
}
