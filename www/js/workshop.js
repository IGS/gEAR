window.onload=function() {
    // check if the user is already logged in
    //check_for_login();
    //session_id = Cookies.get('gear_session_id');
    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    })

    // Load any events for which the user is already registered
    $.ajax({
        url : './cgi/get_event_registration_list.cgi',
        type: "POST",
        data : { 'session_id': session_id,
                 'min_event_id': 1,
                 'max_event_id': 8
               },
        dataType:"json",
        success: function(data, textStatus, jqXHR) {
            // Loop through the events on the page and update status of each
            console.log(data);
        },
        error: function (jqXHR, textStatus, errorThrown) {
            console.log('textStatus= ', textStatus);
            console.log('errorThrown= ', errorThrown);
            display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        }
    });

    // Button states         - Color/class
    //  Register (normal)    - green/success
    //  Register (waitlist)  - blue/primary
    //  Full                 - grey/secondary
    //  Unregister           - red/danger
}
