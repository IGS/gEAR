window.onload=function() {
    // check if the user is already logged in
    //check_for_login();
    //session_id = Cookies.get('gear_session_id');

    // Load any events for which the user is already registered
    $.ajax({
        url : './cgi/get_event_registration_list.cgi',
        type: "POST",
        data : { 'session_id': session_id, },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            // Loop through the events on the page and update status of each
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });
}
