"use strict";

class User {
    constructor ({id, user_name, email, institution, colorblind_mode, updates_wanted, is_admin,
                  help_id, date_created, session_id, default_profile_share_id, default_org_id,
                  layout_share_id} = {}) {

        // SAdkins - I would love to make these properties camelCase as per JS convention, but that would be an effort to fix in all code that uses this class.
        this.user_name = user_name;
        this.email = email;
        this.institution = institution;
        this.colorblind_mode = colorblind_mode === 1;
        this.updates_wanted = updates_wanted === 1;
        // ?: Store is_admin as boolean in database rather than 0/1?
        this.is_admin = is_admin === 1;
        this.help_id = help_id;
        this.date_created = date_created;
        this.session_id = session_id;
        this.default_org_id = default_org_id;
        this.layout_share_id = layout_share_id;

        // derived
        this.default_profile_share_id = default_profile_share_id
    }

    setDefaultProfile() {

         // TODO: Set all this via a config file as it can vary by site and DB
        // This corresponds to "Hearing (site default)"
        this.default_profile_share_id = "f64f9c22";

        if (!this.session_id) {
            return;
        }
        //User is logged in.
        if (Cookies.get('gear_default_domain')) {
            //Get user's default profile
            this.default_profile_share_id = Cookies.get('gear_default_domain');
        }

    }
}
