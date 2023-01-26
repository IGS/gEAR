window.onload=function() {
    // check if the user is already logged in
    //check_for_login();
    //session_id = Cookies.get('gear_session_id');
    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    })

    if (! session_id) {
        $("p.account_needed_message").show();
    }

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
            for (let event_id in data) {
                let event = data[event_id];

                $("#session_" + event_id + "_attendees").html(event['attendees']);
                $("#session_" + event_id + "_seats_left").html(event['max_attendees'] - event['attendees']);
                $("#session_" + event_id + "_max_attendees").text(event['max_attendees']);
                $("#session_" + event_id + "_waitlist_current").html(event['attendees'] - event['max_attendees']);
                $("#session_" + event_id + "_waitlist_size").html(event['waitlist_size']);

                if (event['user_attending'] == 1) {
                    $("ul#session_" + event_id + " li.unregister").show();
                        
                } else if (event['user_waitlisted'] == 1) {
                    $("ul#session_" + event_id + " li.unregister_wait").show();

                } else {
                    if (event['attendees'] < event['max_attendees']) {
                        // user can register for a seat
                        $("ul#session_" + event_id + " li.register").show();
                        
                    } else if (event['attendees'] < (event['max_attendees'] + event['waitlist_size'])) {
                        // user can be added to waitlist
                        $("ul#session_" + event_id + " li.register_wait").show();

                    } else if (event['attendees'] >= (event['max_attendees'] + event['waitlist_size'])) {
                        // event is full
                        $("ul#session_" + event_id + " li.event_full").show();
                        
                    } else {
                        console.log("Unhandled event: " + event_id);
                    }
                }
            }
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
                        $("ul#session_" + event_id + " li.unregister_wait").show();
                    }
                },
                error: function (jqXHR, textStatus, errorThrown) {
                    console.log('textStatus= ', textStatus);
                    console.log('errorThrown= ', errorThrown);
                    display_error_bar(jqXHR.status + ' ' + errorThrown.name);
                }
            });
            
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
                        let event = data['event'];
                        if (data['current_count'] < event['max_attendees']) {
                            $("#session_" + event_id + "_seats_left").html(event['max_attendees'] - data['current_count']);
                            $("ul#session_" + event_id + " li.register").show();
                        } else if (data['current_count'] < (event['max_attendees'] + event['waitlist_size'])) {
                            $("#session_" + event_id + "_waitlist_current").html(data['current_count'] - event['max_attendees']);
                            $("ul#session_" + event_id + " li.register_wait").show();
                        } else {
                            $("ul#session_" + event_id + " li.event_full").show();
                        }
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
