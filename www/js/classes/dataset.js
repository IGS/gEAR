"use strict";

class Dataset {
    constructor ({dataset_id, id=dataset_id, is_public, date_added, dtype, geo_id, ldesc, load_status, math_format, organism_id,
                  pubmed_id, schematic_image, share_id, tags=[], title, user_id, gene_count, obs_count,
                  has_h5ad, platform_id, instrument_model, library_selection, library_source, library_strategy,
                  contact_email, contact_institute, contact_name, annotation_source, plot_default,
                  annotation_release, is_permalink, links} = {}) {
        this.id = id;
        this.is_permalink = is_permalink;
        this.date_added = date_added;
        this.dtype = dtype;
        this.geo_id = geo_id;
        this.is_public = is_public;
        this.ldesc = ldesc;
        this.load_status = load_status;
        this.has_h5ad = has_h5ad;
        this.math_format = math_format;
        this.organism_id = organism_id;
        this.user_id = user_id;
        this.pubmed_id = pubmed_id;
        this.schematic_image = schematic_image;
        this.share_id = share_id;
        this.tags = tags;
        this.title = title;
        this.platform_id = platform_id;
        this.instrument_model = instrument_model;
        this.library_selection = library_selection;
        this.library_source = library_source;
        this.library_strategy = library_strategy;
        this.contact_email = contact_email;
        this.contact_institute = contact_institute;
        this.contact_name = contact_name;
        this.annotation_source = annotation_source;
        this.plot_default = plot_default;
        this.annotation_release = annotation_release;
        this.links = links;

        // derived
        this.gene_count = gene_count;
        this.obs_count = obs_count;
    }

    shape() {
        return this.gene_count + " genes x " + this.obs_count + " obs";
    }

}
