"use strict";

/*
  This class treats a panel representing a gene's functional annotation as a class,
  controlling which annotation data (gene product name, gene symbol, etc.) is displayed.

  Requirements
  - JQuery
  - JSRender
*/

// Does not yet require common.js

class FunctionalAnnotationPanel {
    constructor({annotation} = {}) {
        this.annotation = annotation;
        this.max_organism_id = 0;
        this.organism_ids_present = new Set();
        this.active_organism_id = null;

        for (const dsp of dataset_collection_panel.datasets) {
            this.organism_ids_present.add(dsp.organism_id);

            if (dsp.organism_id > this.max_organism_id) {
                this.max_organism_id = dsp.organism_id;
            }
        }
        
        $(document).on('click', 'button.icon-organism', function(e) {
            annotation_panel.set_active_organism($(this).data('organism-id'));
        });
    }

    add_custom_external_link(ad, label) {
        /*
          Given the JSON-parsed annotation attributes structure this will add
          additional links which need special handling and are not stored
          explicitly in our database.

          This currently isn't working because the DVD doesn't support SSL.  Leaving 
          this here to finish once they do.

          Currently being handled instead within geardb.Gene.load_dbxref_links()
        */
        const url = 'http://deafnessvariationdatabase.org/api?type=gene&method=exists&format=bool&terms=' + ad.gene_symbol;

        // Just trying two different methods.
        if (label == 'DVD2') {
            fetch(url, {
                method: 'GET',
                mode: 'no-cors',
                headers: {
                    'Content-Type': 'application/json'
                },
            }).then(
                response => {
                    console.log(response.json());                    
                }
            ).catch(
                error => {
                    console.log(error);
                }
            );
        } else if (label == 'DVD') {
            // callback=? addition is to try to force JSONP
            $.getJSON(url, function(data) {
                console.log(data);
            });
        }
    }

    autoselect_organism() {
        /*
          Given a full annotation structure this displays the first one.  Here, 'first' is arbitrary
          so we'll start in this priority, skipping any for which annotation isn't found or for 
          which there are no datasets of that organism in this profile:

          mouse > human > zebrafish > chicken

          This really just happens numerically from the organism ID starting at 1 and increasing.

          Returns the ID of the selected organism.
         */
        this.organism_ids_present.forEach(org_id => {
            if (typeof this.annotation !== 'undefined' && 
                typeof this.annotation['by_organism'][org_id] !== 'undefined' &&
                this.organism_ids_present.has(org_id)) {

                if (! this.active_organism_id) {
                    this.set_active_organism(org_id);
                }
            }
        });
        
        return null;
    }

    remove_annotation_attributes() {
        /* Clears all the annotation display fields, usually called before updating them
           so no attributes are carried over accidentally from the previous gene selection
        */
        $('#annot-product').text('');
        $('#gene_details_header_title').text('');
        $('#annot-gene_sym').text('');
        $('#annot-ensembl_id a').text('');
        $('#annot-ensembl_id a').attr('href', '#');
        $('#annot-ensembl_release').text('');
        $('#selected_gene_symbol').text('');
        $('#go_term_count').text('');
        $('#annot-go_terms').html('');
        $('#annot-aliases').html('');
        $('#links_out').html('');
    }

    reset() {
        /*
          Resets/clears the annotation panel
         */
        //this.set_active_organism_button();
        this.remove_annotation_attributes();
    }

    set_active_organism(organism_id) {
        console.log("Setting active organism (" + organism_id + ") with these data:");
        console.log(this.annotation);
        /*
          Updates the annotation panel based on selected organism.  This does 
          both the button array and annotation fields
        */
        this.remove_annotation_attributes();
        this.active_organism_id = organism_id;
        this.set_active_organism_button(organism_id);
        var ad = JSON.parse(this.annotation['by_organism'][organism_id][0]);
        this.update_annotation_attributes(ad);
    }

    set_active_organism_button(organism_id) {
        /*
          Controls the panel of organism buttons and which is the selected one.
        */
        $('.icon-organism').prop('disabled', true);

        const icon_buttons = document.querySelectorAll('.icon-organism');
        for (const btn_el of icon_buttons) {
            btn_el.setAttribute('title', btn_el.dataset.label + " - Not found in this collection of datasets");        

            let icon_image = btn_el.querySelector('img');
            if (icon_image) {
                icon_image.setAttribute("src", "./img/icons/org-" + icon_image.parentElement.dataset.organismId + "-outline-64.svg");
            }
        }

        for (var i = 1; i <= this.max_organism_id; i++) {
            // reset any custom labels
            let org_btn = document.querySelector("button#annot_organism_" + i);
            
            if (typeof this.annotation !== 'undefined' && this.organism_ids_present.has(i)) {
                document.querySelector('#annot_organism_' + i).disabled = false;
                org_btn.setAttribute('title', org_btn.dataset.label);

                if (organism_id == i) {
                    $('#annot_organism_' + i + ' img').prop('src', './img/icons/org-' + i + '-dark-64.svg');
                } else {
                    $('#annot_organism_' + i + ' img').prop('src', './img/icons/org-' + i + '-light-64.svg');
                }
            }
        }
    }

    update_annotation_attributes(ad) {
        /*
          Populates the annotation fields in the panel based on the passed organism/gene annotation
        */
        $('#annot-product').text(ad.product);
        $('#gene_details_header_title').text(ad.product);
        $('#annot-gene_sym').text(ad.gene_symbol);
        $('#annot-ensembl_id a').text(ad.ensembl_id);
        $('#annot-ensembl_id a').attr('href', 'https://www.ensembl.org/id/' + ad.ensembl_id);
        $('#annot-ensembl_release').text(ad.ensembl_release);
        $('#selected_gene_symbol').text(ad.gene_symbol);
        $('#searched_gene').show();

        var go_term_count = ad.go_terms.length;
        $('#go_term_count').text(go_term_count);

        var go_tmpl = $.templates("#tmpl_go_terms");
        var go_html = go_tmpl.render(ad.go_terms);
        $('#annot-go_terms').html(go_html);

        var links_tmpl = $.templates("#tmpl_external_links");
        var links_html = links_tmpl.render(ad.dbxrefs);
        $('#links_out').html(links_html);

        var aliases_tmpl = $.templates("#tmpl_aliases");
        var aliases_html = aliases_tmpl.render(ad.aliases);
        $('#annot-aliases').html(aliases_html);
    }
};
