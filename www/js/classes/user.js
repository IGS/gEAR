"use strict";

class User {
    constructor ({id, name, email, institution, updates_wanted, is_admin,
                  help_id, date_created, session_id} = {}) {
        this.id = id;
        this.user_name = name;
        this.email = email;
        this.institution = institution;
        this.updates_wanted = updates_wanted;
        // Note: Store is_admin as boolean in database rather than 0/1?
        this.is_admin = (is_admin === 1) ? true : false;
        this.help_id = help_id;
        this.date_created = date_created;
        this.session_id = session_id;

        // derived
        this.profile = undefined;
    }

    set_default_profile() {
        session_id = Cookies.get('gear_session_id');

        if (!session_id) {
            //User not logged in.
            if (! this.profile) {
                //User did not change selected profile
                this.profile = "Hearing (site default)"
            }

            //Set profile to 'Hearing (default)'
            // TODO: Get this UI call out of here.
            $('#selected_profile').text(this.profile);
        } else {
            //User is logged in.
            if (! this.profile) {
                //Get user's default profile
                this.profile = Cookies.get('gear_default_domain');
            }
            if (! this.profile) {
                this.profile = "Hearing (site default)"
            }

            // Selected profile is empty if user has a domain selected as primary.
            // Populate it with that domain
            if ( $('#selected_profile').text() == 'Empty') {
                $('#selected_profile').text(this.profile);
                $('#search_param_profile').text(this.profile);
            }
        }
    }
}
