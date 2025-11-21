"use strict";

import { convertToFormData } from "../common.v2.js?v=670b2ed";

export class GeneCart {
    constructor ({id, session_id, label, organism_id, share_id, is_public, is_domain,
                  genes = [], gctype, ldesc} = {}) {
        this.id = id;
        this.session_id = session_id;
        this.label = label;
        this.organism_id = organism_id;
        this.share_id = share_id;
        this.is_public = is_public;
        this.is_domain = is_domain;
        this.genes = genes;
        this.gctype = gctype;
        this.ldesc = ldesc;
    }

    async addCartToDb(callback, errCallback=null) {
        /*
          This method is to save a cart after it has been built in the
          standard way, setting attributes on an instantiated object.
        */

        try {
            const response = await axios.post("./cgi/save_new_genecart_json.cgi", this);
            const {data} = response;
            this.id = data.id
            if (callback) callback(this);
        } catch (error) {
            let msg;
            if (error.response?.data) {
                msg = typeof error.response.data === 'string'
                    ? error.response.data
                    : (error.response.data.message || JSON.stringify(error.response.data));
            } else {
                msg = error?.message || "Something went wrong saving genecart to database.";
            }
            console.error(msg);
            if (errCallback) {
                errCallback(this, msg);
            }
        }
    }

    async addCartToDbFromForm(payload, callback, errCallback=null) {
        /*
          This method is to save a cart by submitting form data.  Once
          completed the object properties are filled in and the callback
          is executed.
        */

        try {

            const formData = convertToFormData(payload);

            // Properly add the new_cart_file to the form data
            if (payload?.new_cart_file) {
                formData.delete("new_cart_file");   // Remove the old one if it exists

                formData.append("new_cart_file", payload.new_cart_file);
            }

            const response = await axios.post("./cgi/save_new_genecart_form.cgi", formData);
            const {data} = response;

            this.id = data.id
            if (callback) callback(this);
        } catch (error) {
            let msg;
            if (error.response?.data) {
                msg = typeof error.response.data === 'string'
                    ? error.response.data
                    : (error.response.data.message || JSON.stringify(error.response.data));
            } else {
                msg = error?.message || "Something went wrong saving genecart to database.";
            }
            console.error(msg);
            if (errCallback) {
                errCallback(this, msg);
            }
        }
    }

    addGene(gene) {
        // Pass a Gene object
        this.genes.push(gene)
    }

    async save(callback=()=>{}, errCallback=()=>{}) {
        /*
        If the 'id' is empty, it's assumed to be new, so an INSERT is
        performed.  Otherwise, if ID is populated this does an
        UPDATE of the existing cart.

        If a callback function is passed it is called after a successful
        save and the gene cart object is passed to it.
        */
        if (this.id) {
            this.updateCartInDb(callback, errCallback);
        } else {
            await this.addCartToDb(callback, errCallback);
        }
    }

    updateCartInDb(callback, errCallback=null) {
        alert("Not implemented yet");
        return false;
    }
}

export class WeightedGeneCart extends GeneCart {
    constructor ({...args} = {}, weightLabels) {
        super(args);
        this.weight_labels = weightLabels || [];
    }
}

class LabeledGeneCart extends GeneCart {
    constructor ({...args} = {}, labels) {
        super(args);
        this.labels = labels || [];
        this.genecarts = [];
    }
}

