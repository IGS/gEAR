"use strict";

class User {
    constructor ({id, user_name, email, institution, colorblind_mode, updates_wanted, is_admin,
                  help_id, date_created, session_id, default_org_id, layout_share_id} = {}) {

        this.user_name = user_name;
        this.email = email;
        this.institution = institution;
        this.colorblind_mode = colorblind_mode === 1;
        this.updates_wanted = updates_wanted === 1;
        this.is_admin = is_admin === 1;
        this.help_id = help_id;
        this.date_created = date_created;
        this.session_id = session_id;
        this.default_org_id = default_org_id;
        this.layout_share_id = layout_share_id;

        // TODO: Check any calls like this and have them save to the DB
        //  Cookies.set('gear_default_domain');
    }
}
