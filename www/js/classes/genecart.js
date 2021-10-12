"use strict";

class GeneCart {
    constructor ({id, session_id, label, organism_id, share_id, is_public,
                  genes = [], gctype, ldesc} = {}) {
        this.id = id;
        this.session_id = session_id;
        this.label = label;
        this.organism_id = organism_id;
        this.share_id = share_id;
        this.is_public = is_public;
        this.genes = genes;
        this.gctype = gctype;
        this.ldesc = ldesc;
    }

    add_cart_to_db(callback, gc) {
        /*
          This method is to save a cart after it has been built in the 
          standard way, setting attributes on an instantiated object.
        */
        $.ajax({
            type: "POST",
            url: "./cgi/save_new_genecart.cgi",
            data: JSON.stringify(this),
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(data) {
                if (callback) {
                    this.id = data['id']
                    callback(this);
                }
            },
            error: function(msg) {
                console.log("error: " + msg);
            }
        });
    }

    add_cart_to_db_from_form(callback, form_data) {
        /*
          This method is to save a cart by submitting form data.  Once
          completed the object properties are filled in and the callback
          is executed.
        */
        $.ajax({
            type: "POST",
            method: "POST",
            url: "./cgi/save_new_genecart.cgi",
            data: form_data,
            contentType: false,
           // contentType: 'multipart/form-data',
            processData: false,
            cache: false,
            success: function(data) {
                if (callback) {
                    this.id = data['id']
                    callback(this);
                }
            },
            error: function(msg) {
                console.log("error: " + msg);
            }
        });
    }
    
    add_gene(gene) {
        // Pass a Gene object
        this.genes.push(gene)
    }

    save(callback) {
        /*
        If the 'id' is empty, it's assumed to be new, so an INSERT is
        performed.  Otherwise, if ID is populated this does an
        UPDATE of the existing cart.

        If a callback function is passed it is called after a successful
        save and the gene cart object is passed to it.
        */
        if (this.id) {
            this.update_cart_in_db(callback);
        } else {
            this.add_cart_to_db(callback, this);
        }
    }
}

