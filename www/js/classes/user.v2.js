"use strict";

class User {
    constructor ({id, name, email, institution, colorblind_mode, updates_wanted, is_admin,
                  help_id, date_created, session_id} = {}) {
        this.id = id;
        this.user_name = name;
        this.email = email;
        this.institution = institution;
        this.colorblind_mode = colorblind_mode === 1;
        this.updates_wanted = updates_wanted === 1;
        // Note: Store is_admin as boolean in database rather than 0/1?
        this.is_admin = is_admin === 1;
        this.help_id = help_id;
        this.date_created = date_created;
        this.session_id = session_id;

        // derived
        this.profile = undefined;
    }

    set_default_profile() {
        const profileElt = document.getElementById("selected_profile");
        const searchParamElt = document.getElementById("search_param_profile");

        if (!session_id) {
            //User not logged in.
            if (! this.profile) {
                //User did not change selected profile
                this.profile = "Hearing (site default)"
            }

            //Set profile to 'Hearing (default)'
            // TODO: Get this UI call out of here.
            profileElt.textContent = this.profile;
            return;
        }
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
        if (profileElt.textContent == 'Empty') {
            profileElt.textContent = this.profile;
            searchParamElt.textContent = this.profile;
        }
    }
}
