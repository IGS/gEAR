"use strict";
/*
  Requirements
  - JQuery
  - JSRender
  - gear/common.js (and that check_for_login() has been called.)
  - gear/classes/dataset.js
*/

// If multiple Bootstrap modals are created, multiple "modal-backdrop"
// divs are created, which can break site functionality.
// This only keeps a single instance of it.
// Source: https://stackoverflow.com/a/44588254
$(document).on("shown.bs.modal", ".modal", () => {
    if ($(".modal-backdrop").length > -1) {
        $(".modal-backdrop").not(":first").remove();
    }
});

class DatasetPanel extends Dataset {
    constructor(
        { grid_position, ...args },
        grid_width,
        multigene = false,
        projection = false,
        colorblind_mode = false,
        controller = null,
    ) {
        super(args);
        this.grid_position = grid_position;
        this.grid_width = grid_width ? grid_width : args.grid_width;
        this.multigene = Number(multigene); // "true/false" works but keeping consistent with the other bool properties
        this.projection = Number(projection);
        this.config = null;
        this.displays = null;
        this.default_display_id = null;
        this.active_display_id = null;
        this.zoomed = false;
        const single_or_multi = this.multigene ? "multi" : "single";
        const genes_or_projection = this.projection ? "projection" : "genes";
        this.primary_key = `${this.id}_${this.grid_position}_${genes_or_projection}_${single_or_multi}`;
        this.projection_id = null;
        this.h5ad_info = null;
        //this.links = args.links;
        //this.linksfoo = "foo";
        this.controller = controller;   // AbortController
    }

    // Call API to return observation information on this dataset
    fetch_h5ad_info (dataset_id, analysis) {
        const other_opts = {}
        if (this.controller) {
            other_opts.signal = this.controller.signal;
        }

        const base = `./api/h5ad/${dataset_id}`;
        const query = analysis ? `?analysis=${analysis.id}` : '';
        return axios.get(`${base}${query}`, other_opts);
    }

    get_dataset_displays(user_id, dataset_id) {
        return $.ajax({
            url: "./cgi/get_dataset_displays.cgi",
            type: "POST",
            data: { user_id, dataset_id },
            dataType: "json",
        });
    }

    get_dataset_display(display_id) {
        return $.ajax({
            url: "./cgi/get_dataset_display.cgi",
            type: "POST",
            data: { display_id },
            dataType: "json",
        });
    }

    get_default_display(user_id, dataset_id, is_multigene = 0) {
        return $.ajax({
            url: "./cgi/get_default_display.cgi",
            type: "POST",
            data: { user_id, dataset_id, is_multigene },
            dataType: "json",
        });
    }

    get_embedded_tsne_display(dataset_id) {
        return $.ajax({
            url: "./cgi/get_embedded_tsne_display.cgi",
            type: "POST",
            data: { dataset_id },
            dataType: "json",
        });
    }

    async draw_chart(gene_symbol, display_id) {
        const data = await this.get_dataset_display(display_id);

        let zoom = false;
        if (this.display?.zoomed) {
            zoom = true;
        }

        data.primary_key = this.primary_key;
        data.controller = this.controller;
        data.projection_id = this.projection_id;
        data.colorblind_mode = CURRENT_USER.colorblind_mode;

        let display;
        if (
            ["bar", "scatter", "violin", "line", "contour", "tsne_dynamic", "tsne/umap_dynamic"].includes(data.plot_type)
        ) {
            display = new PlotlyDisplay(data);
        } else if (data.plot_type === "svg") {
            display = new SVGDisplay(data, this.grid_width);
        } else if (
            ["tsne", "tsne_static", "umap_static", "pca_static"].includes(data.plot_type)
        ) {
            display = new TsneDisplay(data, gene_symbol);
        } else if (data.plot_type === "epiviz") {
            display = new EpiVizDisplay(data, gene_symbol);
        }
        this.display = display;

        // We first draw with the display then zoom in, so whenever
        // gene search is updated, the other displays behind the zoomed
        // is updated too.
        display.draw(gene_symbol);
        if (zoom) display.zoom_in();
    }

    async draw_mg_chart(gene_symbols, display_id) {
        const data = await this.get_dataset_display(display_id);

        let zoom = false;
        if (this.display?.zoomed) {
            zoom = true;
        }

        data.primary_key = this.primary_key;
        data.controller = this.controller;
        data.projection_id = this.projection_id;
        data.colorblind_mode = CURRENT_USER.colorblind_mode;

        let display;
        if (
            ["dotplot", "heatmap", "quadrant", "mg_violin", "volcano"].includes(data.plot_type)
        ) {
            display = new MultigeneDisplay(data, gene_symbols);
        }
        this.display = display;

        // We first draw with the display then zoom in, so whenever
        // gene search is updated, the other displays behind the zoomed
        // is updated too.
        if (this.has_h5ad) {
            display.draw_mg(gene_symbols);
        } else {
            this.show_error(
                "This dataset type does not currently support curated multigene displays."
            );
        }
        if (zoom) display.zoom_in();
    }

    // Draw single-gene plots
    async draw({ gene_symbol } = {}) {
        if (this.display) this.display.clear_display();

        this.show_loading();

        // cache gene_symbol so we can use it to redraw with different display
        this.gene_symbol = gene_symbol;

        if (this.display?.gene_symbol) {
            // Ensure that this display is a single-gene display
            this.draw_chart(gene_symbol, this.display.id);
        } else {
            // first time searching gene and displays have not been loaded
            const { default_display_id } = await this.get_default_display(
                CURRENT_USER.id,
                this.id
            );

            if (default_display_id) {
                this.default_display_id = default_display_id;
                this.draw_chart(gene_symbol, default_display_id);

                // cache all owner/user displays for this panel;
                // If user is the owner, do not duplicate their displays as it can cause the HTML ID to duplicate
                const owner_displays =
                    CURRENT_USER.id === this.user_id
                        ? []
                        : await this.get_dataset_displays(this.user_id, this.id);
                this.owner_displays = owner_displays;

                const user_displays = await this.get_dataset_displays(
                    CURRENT_USER.id,
                    this.id
                );
                this.user_displays = user_displays;

                this.register_events();
            } else {
                // No default display, this really shouldn't happen because
                // owners should always have atleast done this after upload
                this.show_error(
                    "No default display. Create one in the dataset curator."
                );
            }
        }
    }

    // Draw multigene plots
    /**
     * Initialize dash display.
     * @param {Array} gene_symbols - Array of gene symbols
     */
    async draw_mg({ gene_symbols } = {}) {
        if (this.display) this.display.clear_display();

        this.show_loading();

        // cache gene_symbol so we can use it to redraw with different display
        this.gene_symbols = gene_symbols;

        // For datasets with no h5ad (like Epiviz), we cannot call the multigene API to create a display
        // So return an explanation for this.
        if (!this.has_h5ad) {
            this.show_error(
                "This dataset type does not currently support curated multigene displays."
            );
            return;
        }

        // All multigene plot types require a categorical observation to plot from. So there would be no curations for these datasets anyways.
        if (! this.h5ad_info) {
            try {
                const res = await this.fetch_h5ad_info(this.id, undefined )
                this.h5ad_info = res.data;
            } catch (err) {
                if (err.name == "CanceledError") {
                    console.info("Canceled fetching h5ad info for previous request");
                    return;
                }
                this.show_error("Could not retrieve observation info for this dataset.");
                return;
            }
        }
        if (!Object.keys(this.h5ad_info.obs_levels).length) {
            this.show_error(
                "This dataset does not have any categorical observations to plot from."
            );
            return;
        }

        if (this.display?.gene_symbols) {
            // Ensure this display is a multigene display
            this.draw_mg_chart(gene_symbols, this.display.id);
        } else {
            // first time searching gene and displays have not been loaded
            const { default_display_id } = await this.get_default_display(
                CURRENT_USER.id,
                this.id,
                1
            );

            if (default_display_id) {
                this.default_display_id = default_display_id;
                this.draw_mg_chart(gene_symbols, default_display_id);

                // cache all owner/user displays for this panel;
                // If user is the owner, do not duplicate their displays as it can cause the HTML ID to duplicate
                const owner_displays =
                    CURRENT_USER.id === this.user_id
                        ? []
                        : await this.get_dataset_displays(this.user_id, this.id);
                this.owner_displays = owner_displays;

                const user_displays = await this.get_dataset_displays(
                    CURRENT_USER.id,
                    this.id
                );
                this.user_displays = user_displays;

                this.register_events();
            } else {
                // No default display, this really shouldn't happen because
                // owners should always have atleast done this after upload
                this.show_error(
                    "No default display. Create one in the multigene curator."
                );
            }
        }
    }

    async redraw(display_id) {
        // This check is here in case there was no
        // default display rendered, and a user tries
        // to toggle to a different display. There would
        // be no display to clear.
        if (this.display) this.display.clear_display();

        this.show_loading();
        if (this.multigene) {
            this.draw_mg_chart(this.gene_symbols, display_id);
        } else {
            this.draw_chart(this.gene_symbol, display_id);
        }
    }

    // Generate single-gene preview images (Plotly) or plots
    draw_preview_images(display) {
        // check if config has been stringified
        const config =
            typeof display.plotly_config == "string"
                ? JSON.parse(display.plotly_config)
                : display.plotly_config;
        const { gene_symbol } = config;

        display.controller = this.controller;
        display.colorblind_mode = CURRENT_USER.colorblind_mode;

        if (gene_symbol) {
            if (
                ["violin", "bar", "line", "scatter", "contour", "tsne_dynamic", "tsne/umap_dynamic"].includes(display.plot_type)
            ) {
                const d = new PlotlyDisplay(display);
                d.get_data(gene_symbol).then(({ data }) => {
                    const { plot_json } = data;
                    Plotly.toImage(
                        { ...plot_json, ...{static_plot:true} },
                        { height: 500, width: 500 }
                    ).then((url) => {
                        $(`#modal-display-img-${display.id}`).attr("src", url);
                    });
                });
            } else if (display.plot_type === "svg") {
                const target = `modal-display-${display.id}`;
                const d = new SVGDisplay(display, null, null, target);
                d.get_data(gene_symbol).then(({ data }) => {
                    d.draw_chart(data);
                });
            } else {
                // tsne
                const target = `modal-display-${display.id}`;
                const d = new TsneDisplay(display, gene_symbol, null, target);
                d.draw(gene_symbol);
            }
            $(`#modal-display-${display.id}-loading`).hide();
        } else {
            // Hide the container box for the multigene plot
            $(`#modal-display-${display.id}`).hide();
        }
    }

    // Generate multigene preview plots
    draw_preview_images_mg(display) {
        // check if config has been stringified
        const config =
            typeof display.plotly_config === "string"
                ? JSON.parse(display.plotly_config)
                : display.plotly_config;
        const { gene_symbols } = config;

        display.controller = this.controller;
        display.colorblind_mode = CURRENT_USER.colorblind_mode;

        // Draw preview image
        if (gene_symbols) {
            const d = new MultigeneDisplay(display, gene_symbols);
            d.get_data(gene_symbols).then(({ data }) => {
                const { plot_json } = data;
                Plotly.toImage(
                    { ...plot_json, ...{static_plot:true} },
                    { height: 500, width: 500 }
                ).then((url) => {
                    $(`#modal-display-img-${display.id}`).attr("src", url);
                });
            });
            $(`#modal-display-${display.id}-loading`).hide();
        } else {
            // Hide the container box for the single-gene plot
            $(`#modal-display-${display.id}`).hide();
        }
    }

    register_events() {
        const primary_key = this.primary_key;

        // redraw plot when user changes display
        $(`#dataset_${primary_key}`).on("changePlot", (e, display_id) =>
            this.redraw(display_id)
        );

        // zoom event
        $(`#dataset_${primary_key}_zoom_in_control`).click((event) => {
            $(`#dataset_grid`).fadeOut(() => {
                this.display.zoom_in();
            });
        });

        // info event
        $(`#dataset_${primary_key}_info_launcher`).click(() => {
            // SAdkins - found bug where if single-gene/multigene is toggled,
            // two modals will be overlayed.  Dispose of the earlier one.
            //$(`#dataset_${primary_key}_info`).modal('dispose');

            const panel = dataset_collection_panel.datasets.find(
                (d) => d.id == this.id
            );

            // get default display id for this display?

            const { id, title, ldesc, schematic_image } = panel;
            const infobox_tmpl = $.templates("#tmpl_infobox");
            const infobox_html = infobox_tmpl.render({
                dataset_id: id,
                primary_key: this.primary_key,
                title,
                ldesc,
                schematic_image,
            });
            $("#modals_c").html(infobox_html);
            $(`#dataset_${primary_key}_info`).modal("show");
        });

        // Displays modal events
        $(`#dataset_${primary_key}_displays`).click(() => {
            // SAdkins - found bug where if single-gene/multigene is toggled,
            // two modals will be overlayed.  Dispose of the earlier one.
            //$(`#dataset_${primary_key}_displays_modal`).modal('dispose');

            const panel = dataset_collection_panel.datasets.find(
                (d) => d.id == this.id
            );
            const { id, title, user_displays, owner_displays } = panel;

            const displays_tmpl = $.templates("#tmpl_displays_box");
            const displays_html = displays_tmpl.render({
                dataset_id: id,
                primary_key: this.primary_key,
                title,
                user_displays,
                owner_displays,
                active_display_id: this.active_display_id || this.default_display_id,
            });

            $("#modals_c").html(displays_html);

            // change plot when user clicks on a different display
            const display_panel = this;
            $("#modals_c div.modal-display").click(function (event) {
                $("#modals_c img.active-display").removeClass("active-display");
                const display_id = this.id.split("-").pop();
                display_panel.active_display_id = display_id;
                $(`#dataset_${primary_key}`).trigger("changePlot", display_id);
                $(`#modals_c img#modal-display-img-${display_id}`).addClass(
                    "active-display"
                );
            });

            user_displays.forEach((display) => {
                if (multigene) {
                    this.draw_preview_images_mg(display);
                } else {
                    this.draw_preview_images(display);
                }
            });
            owner_displays.forEach((display) => {
                if (multigene) {
                    this.draw_preview_images_mg(display);
                } else {
                    this.draw_preview_images(display);
                }
            });
            $(`#dataset_${primary_key}_displays_modal`).modal("show");
        });
    }

    show_error(msg) {
        $(`#${this.primary_key}_dataset_status_c h2`).text(msg);
        $(`#dataset_${this.primary_key} .dataset-status-container`).show();
        $(`#dataset_${this.primary_key} .plot-container`).hide();
    }

    show_loading() {
        $(`#${this.primary_key}_dataset_status_c h2`).text("Loading...");
        $(`#dataset_${this.primary_key} .dataset-status-container`).show();

        if (this.dtype !== "svg-expression") {
            if (
                this.dtype === "image-static-standard" ||
                this.dtype === "image-static"
            ) {
                $(`#dataset_${this.primary_key} .image-static-container`).hide();
            } else {
                $(`#dataset_${this.primary_key} .plot-container`).hide();
            } // TODO: Hide other container plot types...
        }
    }

    show_no_match() {
        $(`#${this.primary_key}_dataset_status_c h2`).text("Gene not found");
        $(`#dataset_${this.primary_key} .dataset-status-container`).show();
        $(`#dataset_${this.primary_key} .plot-container`).hide();
    }

    show_plot() {
        $(`#dataset_${this.primary_key} .dataset-status-container`).hide();
        $(`#dataset_${this.primary_key} .plot-container`).show();
    }

    show_no_default() {
        const lastForwardSlash = window.location.href.lastIndexOf("/");
        const baseURL = window.location.href.slice(0, lastForwardSlash);
        const link = `${baseURL}/dataset_curator.html#/dataset/${this.id}/displays`;
        $(`#${this.primary_key}_dataset_status_c h2`).html(
            `<a href="${link}"> No default display for this dataset. <br> Click to curate how you\'d like it displayed</a>`
        );
        $(`#dataset_${this.primary_key} .dataset-status-container`).show();
        $(`#dataset_${this.primary_key} .plot-container`).hide();
    }

    // Call API to return plot JSON data
    async run_projectR(projection_source, is_pca, scope) {
        const dataset_id = this.id;
        const payload = {
            genecart_id: projection_source,
            is_pca,
        };
        const other_opts = {}
        if (this.controller) {
            other_opts.signal = this.controller.signal;
        }

        try {
            const response = await axios.post(`api/projectr/${dataset_id}/output_file`, payload, other_opts);
            // If file was not found, put some loading text in the plot
            if (! response.data.projection_id) {
                this.show_loading({
                    info:"Plot generation may take a few minutes as projections need to be generated beforehand."
                });
            }
            payload.projection_id = response.data.projection_id ? response.data.projection_id : null;
        } catch (e) {
            this.show_error(e.message);
            throw(e.message);
        }

        payload.scope = scope;

        let message = "There was an error projecting patterns onto this dataset.";
        try {
            const { data } = await axios.post(`/api/projectr/${dataset_id}`, payload, other_opts);
            if (data.message) {
                message = data.message;
            }
            if (data.success < 1) {
                throw message; // will be caught below
            }
            this.projection_id = data.projection_id;
            this.show_info(message)
        } catch (e) {
            if (e.name == "CanceledError") {
                console.info("Canceled previous projectR request.");
                return;
            }
            this.show_error(message);
            throw message
        }
    }
}
