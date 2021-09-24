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
    constructor({datasets, layout_id, layout_label} = {}) {
        this.datasets = datasets;
        this.layout_id = layout_id;
        this.layout_label = layout_label;

        if (! this.datasets) {
            this.datasets = new Array();
        }
    }

    load_frames({share_id=null, multigene=false} = {}) {
        /*
          Queries the database to get the list of datasets in the user's current
          view.  Initializes the dataset frame panels with placeholders for each
          dataset.
        */
        this.reset();

        // we have to do this because 'this' gets scoped out within the AJAX call
        var dsc_panel = this;

        $.ajax({
            url : './cgi/get_dataset_list.cgi',
            type: "POST",
            data : { 'session_id': session_id, 'permalink_share_id': share_id,
                     'exclude_pending': 1, 'default_domain': this.layout_label,
                     'layout_id': dsc_panel.layout_id },
            dataType: "json",
            success: function(data, textStatus, jqXHR) {
                //console.log(data);
                $.each( data['datasets'], function(i, ds) {
                    // Choose single-gene or multigene grid-width
                    let grid_width = (multigene) ? ds["mg_grid_width"] : ds["grid_width"];
                    var dsp = new DatasetPanel( ds, grid_width, multigene );
                    console.table(dsp);
                    
                    if (dsp.load_status == 'completed') {
                        // reformat the date
                        dsp.date_added = new Date(dsp.date_added);

                        // Insert line-breaks into dataset descriptions if it doesn't look
                        // like HTML already
                        if (dsp.ldesc && ! /<\/?[a-z][\s\S]*>/i.test(dsp.ldesc)) {
                            dsp.ldesc = dsp.ldesc.replace(/(\r\n|\n\r|\r|\n)/g, '<br />');
                        }
                    }

                    dsc_panel.datasets.push(dsp);
                });

                $('span#domain_choice_info_count').text(dsc_panel.datasets.length);

                if (share_id) {
                    const permalinkViewTmpl = $.templates("#tmpl_permalink_info");
                    const permalinkViewHtml = permalinkViewTmpl.render(data['datasets']);
                    $("#permalink_info").html(permalinkViewHtml);
                    const listViewTmpl = $.templates("#tmpl_datasetbox");
                    dsc_panel.datasets.forEach(ds => ds.zoomed = true);
                    const listViewHtml = listViewTmpl.render(dsc_panel.datasets);
                    $('#dataset_grid').html(listViewHtml);
                } else {
                    const listViewTmpl = $.templates("#tmpl_datasetbox");
                    const listViewHtml = listViewTmpl.render(dsc_panel.datasets);
                    $('#dataset_grid').html(listViewHtml);
                }
            },
            error: function (jqXHR, textStatus, errorThrown) {
                display_error_bar(jqXHR.status + ' ' + errorThrown.name);
            }
        });
    }

    reset() {
        this.datasets = [];
        $('#dataset_grid').empty();
    }

    set_layout(layout_id, layout_label, do_load_frames, multigene=false) {
        /*
          Updates this object, the user's stored cookie, and the database and UI labels
        */
        var d = new $.Deferred();

        Cookies.set('gear_default_domain', layout_label);

        this.layout_id = layout_id;
        this.layout_label = layout_label;

        $('#selected_profile').text(layout_label);
        $('#selected_profile').val(layout_id);
        $('#search_param_profile').text(layout_label);

        if (do_load_frames) {
            this.load_frames({multigene});
        }

        // If a user is logged in, we need to save to the db also
        if (CURRENT_USER.session_id) {
            CURRENT_USER.profile = layout_label;

            //Update user's primary layout profile
            $.ajax({
                url: './cgi/set_primary_layout.cgi',
                type: 'post',
                data: {'session_id': CURRENT_USER.session_id, 'layout_id': layout_id},
                success: function(data, newValue, oldValue) {
                    if (data['success'] == 1) {
                        //Was a search already performed?
                        if ($('#search_gene_symbol').val()) {
                            // User has already searched, automatically update datasets and gene searches
                            update_datasetframes_generesults();
                        }
                    } else {
                        $('.alert-container').html('<div class="alert alert-danger alert-dismissible" role="alert">' +
                                                   '<button type="button" class="close close-alert" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                                                   '<p class="alert-message"><strong>Oops! </strong> ' + data["error"] + '</p></div>').show();
                    }
                    d.resolve();
                }
            });
        } else {
            d.resolve();
        }

        return d.promise();
    }

    // Single-gene displays
    update_by_search_result(entry) {
        for (var dataset of this.datasets) {
            if (typeof entry !== 'undefined' &&
                dataset.organism_id in entry['by_organism']) {
                var gene = JSON.parse(entry['by_organism'][dataset.organism_id][0]);
                var gene_symbol = gene.gene_symbol;
                dataset.draw({'gene_symbol':gene_symbol});
            } else {
                if (dataset.display) dataset.display.clear_display();
                dataset.show_no_match();
            }
        }
    }

    // Multigene displays
    update_by_all_results(entries) {
        for (var dataset of this.datasets) {
            if (typeof entries !== 'undefined' ) {
                // TODO: Do something with "by_organism" like single-gene "update_by_search_result"
                // 'entries' is array of gene_symbols
                dataset.draw_mg({'gene_symbols':entries});
            } else {
                if (dataset.display) dataset.display.clear_display();
                dataset.show_no_match();
            }
        }
    }
}
