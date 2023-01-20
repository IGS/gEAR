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
            console.log("Data:");
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
    $(document).on('click', "ul.session_actions button", (e) => {
        e.preventDefault();

        let button = $(e.target);
        let button_type = $(button).parent().attr("class");
        let event_id = $(button).data('event-id');

        if (button_type == 'register') {
            $.ajax({
                url : './cgi/set_user_event_registration.cgi',
                type: "POST",
                data : { 'session_id': session_id,
                         'event_id': event_id,
                         'registration_status': 1
                       },
                dataType:"json",
                success: function(data, textStatus, jqXHR) {
                    if (data.success == 0) {
                        display_error_bar(data.msg);
                    } else {
                        // show only the unregister button now
                        $("ul#session_" + event_id + " li").hide();
                        $("ul#session_" + event_id + " li.unregister").show();
                    }
                },
                error: function (jqXHR, textStatus, errorThrown) {
                    console.log('textStatus= ', textStatus);
                    console.log('errorThrown= ', errorThrown);
                    display_error_bar(jqXHR.status + ' ' + errorThrown.name);
                }
            });

        } else if (button_type == 'register_wait') {

        } else if (button_type == 'unregister' || button_type == 'unregister_wait') {
            $.ajax({
                url : './cgi/set_user_event_registration.cgi',
                type: "POST",
                data : { 'session_id': session_id,
                         'event_id': event_id,
                         'registration_status': 0
                       },
                dataType:"json",
                success: function(data, textStatus, jqXHR) {
                    if (data.success == 0) {
                        display_error_bar(data.msg);
                    } else {
                        // show only the unregister button now
                        $("ul#session_" + event_id + " li").hide();

                        // show register or register_waitlist depending on seats available
                        $("ul#session_" + event_id + " li.register").show();
                    }
                },
                error: function (jqXHR, textStatus, errorThrown) {
                    console.log('textStatus= ', textStatus);
                    console.log('errorThrown= ', errorThrown);
                    display_error_bar(jqXHR.status + ' ' + errorThrown.name);
                }
            });
        } 

    });

}
