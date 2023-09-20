"use strict";

/*
  Requirements
  - JQuery
  - JSRender
  - gear/common.js (and that check_for_login() has been called.)
  - gear/classes/dataset.js
*/

// Does not yet require common.js
class DatasetCollectionPanel {
    constructor({ datasets, layout_id, layout_label } = {}) {
        this.datasets = datasets;
        this.layout_id = layout_id;
        this.layout_label = layout_label;
        this.search_performed = false;

        if (!this.datasets) {
            this.datasets = new Array();
        }

    }

    load_frames({
        dataset_id = null,
        multigene = false,
        projection = false,
    } = {}) {
        /*
              Queries the database to get the list of datasets in the user's current
              view.  Initializes the dataset frame panels with placeholders for each
              dataset.
            */

        this.reset();

        // we have to do this because 'this' gets scoped out within the AJAX call
        const dsc_panel = this;
        $.ajax({
            url: "./cgi/get_dataset_list.cgi",
            type: "POST",
            async: false, // Adding so datasets are updated before the set_layout() AJAX call happens
            data: {
                session_id,
                permalink_share_id: dataset_id,
                exclude_pending: 1,
                default_domain: this.layout_label,
                layout_id: dsc_panel.layout_id,
            },
            dataType: "json"
        }).done((data) => {
            $.each(data["datasets"], (_i, ds) => {
                // Choose single-gene or multigene grid-width
                const grid_width = multigene ? ds.mg_grid_width : ds.grid_width;

                this.reset_abort_controller();
                const dsp = new DatasetPanel(ds, grid_width, multigene, projection, dsc_panel.controller);

                if (dsp.load_status == "completed") {
                    // reformat the date
                    dsp.date_added = new Date(dsp.date_added);

                    // Insert line-breaks into dataset descriptions if it doesn't look
                    // like HTML already
                    if (dsp.ldesc && !/<\/?[a-z][\s\S]*>/i.test(dsp.ldesc)) {
                        dsp.ldesc = dsp.ldesc.replace(/(\r\n|\n\r|\r|\n)/g, "<br />");
                    }
                }
                dsc_panel.datasets.push(dsp);
            });

            $("span#domain_choice_info_count").text(dsc_panel.datasets.length);

            if (dataset_id) {
                const permalinkViewTmpl = $.templates("#tmpl_permalink_info");
                const permalinkViewHtml = permalinkViewTmpl.render(data["datasets"]);
                $("#permalink_info").html(permalinkViewHtml);
                dsc_panel.datasets.forEach((ds) => (ds.zoomed = true));
            }
            const listViewTmpl = $.templates("#tmpl_datasetbox");
            const listViewHtml = listViewTmpl.render(dsc_panel.datasets);
            $("#dataset_grid").html(listViewHtml);

            // Anytime this happens, the annotation pane needs to be (re-)loaded
            try {
                if (annotation_panel) {
                    annotation_panel = new FunctionalAnnotationPanel();
                }
            } catch (e) {
                console.warn(e);
            }

        }).fail((jqXHR, textStatus, errorThrown) => {
            display_error_bar(`${jqXHR.status} ${errorThrown.name}`);
        });
    }

    reset() {
        this.datasets = [];
        $("#dataset_grid").empty();
        this.reset_abort_controller();
    }

    reset_abort_controller() {
        if (this.controller && !this.performing_projection) {
            this.controller.abort(); // Cancel any previous axios requests (such as drawing plots for a previous dataset)
        }
        this.controller = new AbortController(); // Create new controller for new set of frames

        for (const dataset of this.datasets) {
            dataset.controller = this.controller;
        }
    }

    set_layout(
        layout_id,
        layout_share_id,
        layout_label,
        do_load_frames,
        multigene = false,
        projection = false
    ) {
        /*
              Updates this object, the user's stored cookie, and the database and UI labels
            */
        Cookies.set("gear_default_domain", layout_label);

        this.layout_id = layout_id;
        this.layout_share_id = layout_share_id;
        this.layout_label = layout_label;

        // According to this, we shouldn't be manually setting the label:
        //  https://github.com/vitalets/x-editable/issues/332#issuecomment-22379178
        $("#selected_profile").text(layout_label);

        $("#selected_profile").val(layout_id);
        $("#search_param_profile").text(layout_label);

        if (do_load_frames) {
            this.load_frames({ multigene, projection });
        }

        // If a user is logged in, we need to save to the db also
        if (CURRENT_USER.session_id) {
            CURRENT_USER.profile = layout_label;

            //Update user's primary layout profile
            $.ajax({
                url: "./cgi/set_primary_layout.cgi",
                type: "post",
                data: { session_id: CURRENT_USER.session_id, layout_id: layout_id }
            }).done((data) => {
                if (!(data.success == 1)) {
                    $(".alert-container")
                        .html(
                            '<div class="alert alert-danger alert-dismissible" role="alert">' +
                            '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                            '<p class="alert-message"><strong>Oops! </strong> ' +
                            data.error +
                            "</p></div>"
                        )
                        .show();
                }
            }).fail((jqXHR, _textStatus, errorThrown) => {
                display_error_bar(`${jqXHR.status} ${errorThrown.name}`);
            });
        }

        //Was a search already performed?
        if (this.search_performed && !projection) {
            // User has already searched, automatically update datasets and gene searches
            update_datasetframes_generesults();
        }
    }

    // Single-gene displays
    update_by_search_result(entry) {
        this.reset_abort_controller();

        for (const dataset of this.datasets) {
            if (
                typeof entry !== "undefined" &&
                dataset.organism_id in entry.by_organism
            ) {
                // If working with actual genes, ensure dataset and entry's organisms match for annotation purposes.
                const gene = JSON.parse(entry.by_organism[dataset.organism_id][0]);
                const gene_symbol = gene.gene_symbol;
                dataset.draw({ gene_symbol });
            } else {
                if (dataset.active_display) dataset.active_display.clear_display();
                dataset.show_no_match();
            }
        }
    }

    // Multigene displays
    // Only executes in "gene" mode
    update_by_all_results(entries) {
        this.reset_abort_controller();

        for (const dataset of this.datasets) {
            if (typeof entries !== "undefined") {
                // TODO: Do something with "by_organism" like single-gene "update_by_search_result"
                // 'entries' is array of gene_symbols
                dataset.draw_mg({ gene_symbols: entries });
            } else {
                if (dataset.active_display) dataset.active_display.clear_display();
                dataset.show_no_match();
            }
        }
    }
}
