"use strict";

class GeneCart {
    constructor ({id, session_id, label, genes = []} = {}) {
        this.id = id;
        this.session_id = session_id;
        this.label = label;
        this.genes = genes;
    }

    add_cart_to_db(callback, gc) {
        $.ajax({
            type: "POST",
            url: "./cgi/save_new_genecart.cgi",
            data: JSON.stringify(this),
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            success: function(data) {
                if (callback) {
                    gc.id = data['id']
                    callback(gc);
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

