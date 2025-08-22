"use strict";

import { apiCallsMixin } from '../common.v2.js?v=3dd4d96';

export class User {
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
    }

    // Write a function to save the layout_share_id to the database via a CGI script
    saveLayoutShareId(shareId) {
        this.layout_share_id = shareId;
        apiCallsMixin.setUserPrimaryDatasetCollection(shareId);
    }
}
