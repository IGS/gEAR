/** Base class representing a display */
class Display {
    /**
     * Initialize display.
     * @param {Object} - Display data.
     */
    constructor({
        id,
        dataset_id,
        user_id,
        label,
        plot_type,
        primary_key
    }) {
        this.id = id;
        this.dataset_id = dataset_id;
        this.user_id = user_id;
        this.label = label;
        this.plot_type = plot_type;
        this.primary_key = primary_key; // Combination of dataset_id, grid_position, and single/multigene indicator
        this.data = null;
        this.first_draw = true; // Keep track of if this is the original draw so that effects are not doubly applied.
        this.zoomed = false;
    }
    zoom_in() {
        // data is already fetched, so all we need to do
        // is redraw the data in a bigger container.
        this.zoomed = true;
        const zoomedViewTmpl = $.templates('#tmpl_dataset_zoomed');
        const dataset_panel = dataset_collection_panel.datasets.find(
            dataset => dataset.id === this.dataset_id
        );

        const panel_data = {
            ...dataset_panel
        };
        const zoomedViewHtml = zoomedViewTmpl.render(panel_data);
        $('#dataset_zoomed').html(zoomedViewHtml);

        $('#dataset_zoom_out_btn').click(() => {
            $('#dataset_zoomed').fadeOut(() => {
                this.zoomed = false;
                $('#dataset_grid').fadeIn();

                // now we want to redraw display panels
                // so they are stay contained inside their container
                dataset_collection_panel.datasets.forEach(dataset => {
                    const display = dataset.display;
                    // check for display data, it could be null
                    // if gene was not found.
                    if (display?.data) {
                        // TODO: Handle here if it's a 0-white color range
                        display.draw_chart(display.data);
                        // Exit status 2 is status to show plot but append warning message
                        if (display.data.message && display.data.success === 2)
                            display.show_warning(display.data.message);
                    }
                });
            });
        });

        $('#dataset_zoomed').fadeIn();
        this.draw_zoomed();
    }
    /**
     *  Draw the visualization.
     * @param {string} gene_symbol - Gene Symbol to visualize.
     */
    async draw(gene_symbol) {
        this.gene_symbol = gene_symbol;
        const {
            data
        } = await this.get_data(gene_symbol);
        if (data.success === -1) {
            this.show_error(data.message);
            return;
        }
        this.data = data;
        if (this.zoomed) {
            this.draw_zoomed();
        } else {
            let attempts_left = 3;

            while (attempts_left) {
                var draw_success = this.draw_chart(data);
                //console.log(this.dataset_id + " - Drawing attempts left: " + attempts_left + " success: " + draw_success);
                attempts_left -= 1;

                // if it didn't work, wait one second and try again
                if (draw_success) {
                    attempts_left = 0;
                } else if (attempts_left) {
                    // else the scope is lost in the anon function below
                    var that = this;
                    setTimeout(
                        () => {
                            draw_success = that.draw_chart(data);
                            attempts_left -= 1;
                        }, 1000);
                }
            }
        }
        // Exit status 2 is status to show plot but append warning message
        if (data.message && data.success === 2)
        this.show_warning(this.data.message);
    }

    /**
     *  Draw the multigene visualization.
     * @param {string} gene_symbols - Gene Symbols to visualize.
     */
     async draw_mg(gene_symbols) {
        this.gene_symbols = gene_symbols;
        const {
            data
        } = await this.get_data(gene_symbols);
        if (data.success === -1) {
            this.show_error(data.message);
            return;
        }
        this.data = data;
        if (this.zoomed) {
            this.draw_zoomed();
        } else {
            let attempts_left = 3;

            while (attempts_left) {
                var draw_success = this.draw_chart(data);
                //console.log(this.dataset_id + " - Drawing attempts left: " + attempts_left + " success: " + draw_success);
                attempts_left -= 1;

                // if it didn't work, wait one second and try again
                if (draw_success) {
                    attempts_left = 0;
                } else if (attempts_left) {
                    // else the scope is lost in the anon function below
                    var that = this;
                    setTimeout(
                        () => {
                            draw_success = that.draw_chart(data);
                            attempts_left -= 1;
                        }, 1000);
                }
            }
        }
        // Exit status 2 is status to show plot but append warning message
        if (data.message && data.success === 2)
        this.show_warning(this.data.message);
    }
    /**
     * Hides the display container
     */
    hide_loading() {
        $(`#dataset_${this.primary_key} .dataset-status-container`).hide();
    }

    show_loading() {
        $(`#${this.primary_key}_dataset_status_c h2`).text('Loading...');
        $(`#dataset_${this.primary_key} .dataset-status-container`).show();

        if (this.dtype !== 'svg-expression') {
            if (
                this.dtype === 'image-static-standard' ||
                this.dtype === 'image-static'
            ) {
                $(`#dataset_${this.primary_key} .image-static-container`).hide();
            } else {
                $(`#dataset_${this.primary_key} .plot-container`).hide();
            } // TODO: Hide other container plot types...
        }
    }

    /**
     *  Display the error from the server.
     * @param {string} msg - Message returned back from server
     */
    show_error(msg) {
        $(`#${this.primary_key}_dataset_status_c h2`).text(msg);
        $(`#dataset_${this.primary_key} .dataset-status-container`).show();
        $(`#dataset_${this.primary_key} .plot-container`).hide();
    }

    /**
     * Show warning overlay above plot
     */
    show_warning(msg) {}

    /**
     * Show the plot container.
     */
    show() {
        $(`#dataset_${this.primary_key}_h5ad.plot-container`).show();
    }

    /**
     * Add label that shows full text of axis labels upon hover of axis
     */
    create_hover_info_area(target_div) {
        $(`#${target_div}`).append(`<div class="hoverarea" ></div>`);
    }
}

/**
 * Class representing an EpiViz display
 * @extends Display
 */
class EpiVizDisplay extends Display {
    /**
     * Initialize tSNE
     * This subclass of display takes an extra argument, gene symbol.
     * This is because we draw this one differently and don't query
     * the server for data, but query the server for the image.
     * @param {Object} Data - Display data
     * @param {string} gene_symbol - Gene symbol to visualize
     * @param {String} projection_csv - Basename of CSV file containing projection data
     */
    constructor({
        plotly_config,
        ...args
    }, gene_symbol, projection_csv, target) {
        super(args);
        const config = plotly_config;
        this.gene_symbol = gene_symbol;
        this.econfig = config;
        this.extendRangeRatio = 10;
        this.projection_csv = projection_csv;   // Maybe use eventually but I don't think we can project on Epiviz displays

        const genes_track = plotly_config.tracks["EPIVIZ-GENES-TRACK"];
        if (genes_track.length > 0) {
            const gttrack = genes_track[0];
            this.genome = gttrack.measurements ? gttrack.measurements[0].id : gttrack.id[0].id;
        }

        this.target = target ? target : `dataset_${this.primary_key}`;

        this.epiviztemplate = this.epiviz_template(config);
    }
    /**
     * Not sure if this one is applicable for this display since
     * EpiViz pulls data for its own tracks based on the config.
     *
     * @param {string} gene_symbol - Gene symbol to visualize.
     */
    get_data(gene_symbol) {
        const base = `/api/plot/${this.dataset_id}/epiviz`;
        const query = `?gene=${gene_symbol}&genome=${this.genome}`;
        return axios.get(`${base}${query}`);
        // return null;
    }

    draw_zoomed() {
        this.clear_display();
        $(`#dataset_zoomed div#dataset_${this.primary_key}`).append(this.template());
    }

    /**
     * Draw EpiViz chart.
     */
    draw_chart() {
        // this.clear_display();
        // this.data.start = parseInt(this.data.start);
        // this.data.end = parseInt(this.data.end);
        // epiviz container doesn't exist
        if ($(`#${this.target}_epiviznav`).length == 0) {
            const target_div = `#${this.target}`;
            $(target_div).append(this.template());
        } else {
            // epiviz container already exists, so only update gneomic position in the browser
            const epiviznav = document.querySelector(`#${this.target}_epiviznav`);
            epiviznav.setAttribute("chr", this.data.chr);
            const nstart = this.data.start - Math.round((this.data.end - this.data.start) * this.extendRangeRatio);
            const nend = this.data.end + Math.round((this.data.end - this.data.start) * this.extendRangeRatio)
            epiviznav.setAttribute("start", nstart);
            epiviznav.setAttribute("end", nend);
            epiviznav.range = epiviznav.getGenomicRange(this.data.chr, nstart, nend);
        }

        this.hide_loading();

        this.show();

        return true;
    }

    epiviz_template(config) {
        let epiviztemplate = "";
        for (const track in config.tracks) {
            const track_config = config.tracks[track];
            track_config.forEach((tc) => {
                var temp_track = `<${track} slot='charts' `;
                temp_track += Object.keys(tc).includes("id") ? ` dim-s='${JSON.stringify(tc.id)}' ` : ` measurements='${JSON.stringify(tc.measurements)}' `;

                if (tc.colors != null) {
                    temp_track += ` chart-colors='${JSON.stringify(tc.colors)}' `;
                }

                if (tc.settings != null) {
                    temp_track += ` chart-settings='${JSON.stringify(tc.settings)}' `;
                }

                temp_track += ` style='min-height:200px;'></${track}> `;

                epiviztemplate += temp_track;
            })
        }

        return epiviztemplate
    }

    /**
     * HTML template for EpiViz
     */
    template() {

        // the chr, start and end should come from query - map gene to genomic position.

        return `
        <div id='${this.target}_epiviz' class='epiviz-container'>
            <script src="https://cdn.jsdelivr.net/gh/epiviz/epiviz-chart/cdn/renderingQueues/renderingQueue.js"></script>
            <script src="https://cdn.jsdelivr.net/gh/epiviz/epiviz-chart/cdn/webcomponentsjs/webcomponents-lite.js"></script>

            <link rel="import" href="https://cdn.jsdelivr.net/gh/epiviz/epiviz-chart/cdn/epiviz-components-gear.html">

            <epiviz-data-source provider-type="epiviz.data.WebServerDataProvider"
            id='${this.target}_epivizds'
            provider-id="fileapi"
            provider-url="${this.econfig.dataserver}">
            </epiviz-data-source>
            <epiviz-navigation
            hide-chr-input
            hide-search
            hide-add-chart
            show-viewer
            id='${this.target}_epiviznav'
            chr='${this.data.chr}'
            start=${this.data.start - Math.round((this.data.end - this.data.start) * this.extendRangeRatio)}
            end=${this.data.end + Math.round((this.data.end - this.data.start) * this.extendRangeRatio)}
            viewer=${`/epiviz.html?dataset_id=${this.id}&chr=${this.data.chr}&start=${this.data.start}&end=${this.data.end}`}
            >
            ${this.epiviztemplate}
            </epiviz-navigation>
        </div>
        `;
    }

    clear_display() {
        // $(`#${this.target}_epiviz epiviz-navigation`).remove();
        if (this.zoomed) $(`#dataset_${this.primary_key}_epiviz_zoomed`).remove();
    }
}

/**
 * Class representing a display drawn with Plotly
 * @extends Display
 */
class PlotlyDisplay extends Display {
    /**
     * Initialize plotly display.
     * @param {Object} Data - Data used to draw violin, bar, or line.
     * @param {String} projection_csv - Basename of CSV file containing projection data
     */
    constructor({
        plotly_config,
        ...args
    }, projection_csv) {
        super(args);
        const {
            x_axis,
            y_axis,
            z_axis,
            x_min,
            y_min,
            x_max,
            y_max,
            x_title,
            y_title,
            vlines,
            hide_x_labels,
            hide_y_labels,
            hide_legend,
            point_label,
            color_name,
            facet_row,
            facet_col,
            marker_size,
            jitter,
            colors,
            color_palette,
            reverse_palette,
            order,
            analysis,
        } = plotly_config;
        this.x_axis = x_axis;
        this.y_axis = y_axis;
        this.z_axis = z_axis;
        this.x_min = x_min;
        this.y_min = y_min;
        this.x_max = x_max;
        this.y_max = y_max;
        this.x_title = x_title;
        this.y_title = y_title;
        this.vlines = vlines;
        this.hide_x_labels = hide_x_labels;
        this.hide_y_labels = hide_y_labels;
        this.hide_legend = hide_legend;
        this.point_label = point_label;
        this.color_name = color_name;
        this.facet_row = facet_row;
        this.facet_col = facet_col;
        this.marker_size = marker_size;
        this.jitter = jitter;
        this.colors = colors;
        this.color_palette = color_palette;
        this.reverse_palette = reverse_palette;
        this.order = order;
        this.analysis = analysis;
        this.projection_csv = projection_csv;
    }
    clear_display() {
        $(`#dataset_${this.primary_key}_h5ad`).remove();
        if (this.zoomed) $(`#dataset_${this.primary_key}_h5ad_zoomed`).remove();
    }
    clear_zoomed() {
        $(`#dataset_${this.primary_key}_h5ad_zoomed`).remove();
    }
    /**
     * Query the server for data to draw plotly chart.
     * @param {string} gene_symbol - Gene symbol to visualize.
     */
    get_data(gene_symbol) {
        return axios.post(`/api/plot/${this.dataset_id}`, {
            plot_type: this.plot_type,
            analysis_owner_id: this.user_id,
            gene_symbol,
            x_axis: this.x_axis,
            y_axis: this.y_axis,
            z_axis: this.z_axis,
            x_min: this.x_min,
            y_min: this.y_min,
            x_max: this.x_max,
            y_max: this.y_max,
            x_title: this.x_title,
            y_title: this.y_title,
            vlines: this.vlines,
            point_label: this.point_label,
            hide_x_labels: this.hide_x_labels,
            hide_y_labels: this.hide_y_labels,
            hide_legend: this.hide_legend,
            color_name: this.color_name,
            facet_row: this.facet_row,
            facet_col: this.facet_col,
            marker_size: this.marker_size,
            jitter: this.jitter,
            colors: this.colors,
            color_palette: this.color_palette,
            reverse_palette: this.reverse_palette,
            order: this.order,
            analysis: this.analysis,
            projection_csv: this.projection_csv
        });
    }
    /**
     * Draw chart.
     * @param {*} data - Plotly data.
     */
    draw_chart(data) {
        this.clear_display();

        const target_div = `dataset_${this.primary_key}_h5ad`;
        const {
            plot_json,
        } = data;

        this.hide_loading();
        $(`#dataset_${this.primary_key}`).append(this.template());
        this.show();

        // Get config
        const index_conf = post_plotly_config.index;
        const plot_config = get_plotly_updates(index_conf, this.plot_type, "config");

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const update_layout = get_plotly_updates(index_conf, this.plot_type, "layout");

        Plotly.newPlot(target_div, plot_json.data, plot_json.layout, plot_config);
        Plotly.relayout(target_div, update_layout)

        //truncateAxisLabels();
        //this.create_hover_info_area(target_div);

        return true;
    }

    draw_zoomed() {
        this.clear_display();
        const target_div = `dataset_${this.primary_key}_h5ad_zoomed`;
        // const target_div = document
        //   .querySelector(`#dataset_zoomed div#dataset_${this.dataset_id}`);

        $(`#dataset_zoomed div#dataset_${this.primary_key}`).append(
            this.zoomed_template()
        );

        this.hide_loading();
        const {
            plot_json,
        } = this.data;

        const index_conf = post_plotly_config.index;
        const plot_config = get_plotly_updates(index_conf, this.plot_type, "config");

        Plotly.newPlot(target_div, plot_json.data, plot_json.layout, plot_config);

        // Do not need to add plot_display_config stuff since this.data has the updates already
    }

    show() {
        $(`#dataset_${this.primary_key} .plot-container`).show();
    }
    template() {
        return `
        <div
          id='dataset_${this.primary_key}_h5ad'
          class="h5ad-container"
          style="position: relative;">
        </div>
        `;
    }
    zoomed_template() {
        return `
        <div
          style='height:70vh;'
          id='dataset_${this.primary_key}_h5ad_zoomed'
          class="h5ad-container">
        </div>
        `;
    }

    /**
     * Display warning above the plot
     */
    show_warning(msg) {

        const hover_msg = " Hover to see warning.";

        const dataset_selector = $( `#dataset_${this.primary_key}_h5ad` );
        const warning_template = `
        <div class='dataset-warning bg-warning' id='dataset_${this.primary_key}_h5ad_warning'>
            <i class='fa fa-exclamation-triangle'></i>
            <span id="dataset_${this.primary_key}_h5ad_msg">${hover_msg}</span>
        </div>`;

        // Add warning above the chart
        // NOTE: must add to DOM before making selector variables
        dataset_selector.prepend(warning_template);

        const warning_selector = $( `#dataset_${this.primary_key}_h5ad_warning` );
        const msg_selector = $( `#dataset_${this.primary_key}_h5ad_msg` );

        // Add some CSS to warning to keep at top of container and not push display down
        warning_selector.css('position', 'absolute').css('z-index', '2');

        // Add hover events (via jQuery)
        warning_selector.mouseover(() => {
            msg_selector.html(msg);
        });
        warning_selector.mouseout(() => {
            msg_selector.text(hover_msg);
       });
    }
}

/**
 * Class representing a multigene display drawn with Dash
 * @extends Display
 */
 class MultigeneDisplay extends Display {
    /**
     * Initialize dash display.
     * @param {Object} Data - Data used to draw any multigene plot
     * @param {Array} gene_symbols - Array of gene symbols
     * @param {String} projection_csv - Basename of CSV file containing projection data
     */
    constructor({
        plotly_config,
        ...args
    }, gene_symbols, projection_csv) {
        super(args);
        const {
            primary_col,
            secondary_col,
            sort_order,
            obs_filters,
            clusterbar_fields,
            matrixplot,
            cluster_obs,
            cluster_genes,
            flip_axes,
            distance_metric,
            adj_pvals,
            annotate_nonsignificant,
            include_zero_fc,
            fold_change_cutoff,
            fdr_cutoff,
            de_test_algo,
            pvalue_threshold,
            lower_logfc_threshold,
            upper_logfc_threshold,
            query_condition,
            compare1_condition,
            compare2_condition,
            ref_condition,
            stacked_violin,
            violin_add_points,
            plot_title,
            legend_title,
            analysis,   // Analysis
        } = plotly_config;
        this.gene_symbols = gene_symbols;
        this.primary_col = primary_col;
        this.secondary_col = secondary_col;
        this.sort_order = sort_order;
        this.obs_filters = obs_filters;
        this.clusterbar_fields = clusterbar_fields;
        this.matrixplot = matrixplot;
        this.cluster_obs = cluster_obs;
        this.cluster_genes = cluster_genes;
        this.flip_axes = flip_axes;
        this.distance_metric = distance_metric,
        this.adj_pvals = adj_pvals;
        this.annot_nonsig = annotate_nonsignificant;
        this.include_zero_fc = include_zero_fc;
        this.fold_change_cutoff = fold_change_cutoff;
        this.fdr_cutoff = fdr_cutoff;
        this.de_test_algo = de_test_algo;
        this.pvalue_threshold = pvalue_threshold;
        this.lower_logfc_threshold = lower_logfc_threshold;
        this.upper_logfc_threshold = upper_logfc_threshold;
        this.query_condition = query_condition;
        this.compare1_condition = compare1_condition;
        this.compare2_condition = compare2_condition;
        this.ref_condition = ref_condition;
        this.stacked_violin = stacked_violin;
        this.violin_add_points = violin_add_points;
        this.plot_title = plot_title;
        this.legend_title = legend_title;
        this.analysis = analysis;
        this.projection_csv = projection_csv;
    }
    clear_display() {
        $(`#dataset_${this.primary_key}_mg`).remove();
        if (this.zoomed) $(`#dataset_${this.primary_key}_mg_zoomed`).remove();
    }
    clear_zoomed() {
        $(`#dataset_${this.primary_key}_mg_zoomed`).remove();
    }

    /**
     * Query the server for data to draw dash chart.
     * @param {string} gene_symbols - Gene symbols to visualize.
     */
    async get_data(gene_symbols) {
        return axios.post(`/api/plot/${this.dataset_id}/mg_dash`, {
            gene_symbols,
            analysis: this.analysis,
            analysis_owner_id: this.user_id,
            primary_col: this.primary_col,
            secondary_col: this.secondary_col,
            sort_order: this.sort_order,
            plot_type: this.plot_type,
            obs_filters: this.obs_filters,
            clusterbar_fields: this.clusterbar_fields,
            matrixplot: this.matrixplot,
            cluster_obs: this.cluster_obs,
            cluster_genes: this.cluster_genes,
            flip_axes: this.flip_axes,
            distance_metric: this.distance_metric,
            adj_pvals: this.adj_pvals,
            annotate_nonsignificant: this.annot_nonsig,
            include_zero_fc: this.include_zero_fc,
            fold_change_cutoff: this.fold_change_cutoff,
            fdr_cutoff: this.fdr_cutoff,
            de_test_algo: this.de_test_algo,
            pvalue_threshold: this.pvalue_threshold,
            lower_logfc_threshold: this.lower_logfc_threshold,
            upper_logfc_threshold: this.upper_logfc_threshold,
            query_condition: this.query_condition,
            compare1_condition: this.compare1_condition,
            compare2_condition: this.compare2_condition,
            ref_condition: this.ref_condition,
            stacked_violin: this.stacked_violin,
            violin_add_points: this.violin_add_points,
            plot_title: this.plot_title,
            legend_title: this.legend_title,
            projection_csv: this.projection_csv,
        });
    }
    /**
     * Draw chart.
     * @param {*} data - Plotly data.
     */
    draw_chart(data) {
        this.clear_display();

        const target_div = `dataset_${this.primary_key}_mg`;
        const {
            plot_json,
        } = data;

        this.hide_loading();
        $(`#dataset_${this.primary_key}`).append(this.template());
        this.show();

        // Get config
        const index_conf = post_plotly_config.index;
        const plot_config = get_plotly_updates(index_conf, this.plot_type, "config");

        if (this.plot_type == "heatmap" && this.first_draw) {
            // These functions modify `this.data` by reference, so we need to ensure we only call them once.
            // Subsequent draw_chart calls will have the correct data, post-modification
            this.first_draw = false;
            adjustExpressionColorbar(plot_json.data);
            adjustClusterColorbars(plot_json.data);
        }

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const update_layout = get_plotly_updates(index_conf, this.plot_type, "layout");

        Plotly.newPlot(target_div, plot_json.data, plot_json.layout, plot_config);
        Plotly.relayout(target_div, update_layout)

        //truncateAxisLabels();
        //this.create_hover_info_area(target_div);

        return true;
    }

    draw_zoomed() {
        this.clear_display();
        const target_div = `dataset_${this.primary_key}_mg_zoomed`;
        // const target_div = document
        //   .querySelector(`#dataset_zoomed div#dataset_${this.dataset_id}`);

        $(`#dataset_zoomed div#dataset_${this.primary_key}`).append(
            this.zoomed_template()
        );

        this.hide_loading();
        const {
            plot_json,
        } = this.data;

        // Get config
        const index_conf = post_plotly_config.index;
        const plot_config = get_plotly_updates(index_conf, this.plot_type, "config");

        Plotly.newPlot(target_div, plot_json.data, plot_json.layout, plot_config);

        // Do not need to add plot_display_config stuff since this.data has the updates already

    }

    show() {
        $(`#dataset_${this.primary_key} .plot-container`).show();
    }
    template() {
        return `
        <div
          id='dataset_${this.primary_key}_mg'
          class="h5ad-container"
          style="position: relative;">
        </div>
        `;
    }
    zoomed_template() {
        return `
        <div
          style='max-width:96%; height:70vh;'
          id='dataset_${this.primary_key}_mg_zoomed'
          class="h5ad-container">
        </div>
        `;
    }

    /**
     * Display warning above the plot
     */
    show_warning(msg) {

        const hover_msg = " Hover to see warning.";

        const dataset_selector = $( `#dataset_${this.primary_key}_mg` );
        const warning_template = `
        <div class='dataset-warning bg-warning' id='dataset_${this.primary_key}_mg_warning'>
            <i class='fa fa-exclamation-triangle'></i>
            <span id="dataset_${this.primary_key}_mg_msg">${hover_msg}</span>
        </div>`;

        // Add warning above the chart
        // NOTE: must add to DOM before making selector variables
        dataset_selector.prepend(warning_template);

        const warning_selector = $( `#dataset_${this.primary_key}_mg_warning` );
        const msg_selector = $( `#dataset_${this.primary_key}_mg_msg` );

        // Add some CSS to warning to keep at top of container and not push display down
        warning_selector.css('position', 'absolute').css('z-index', '2');

        // Add hover events (via jQuery)
        warning_selector.mouseover(() => {
            msg_selector.html(msg);
        });
        warning_selector.mouseout(() => {
            msg_selector.text(hover_msg);
       });
    }

}

/**
 * Class representing an SVG display.
 * @extends Display
 */
class SVGDisplay extends Display {
    /**
     * Initialize SVG
     * @param {Object} data - SVG display data
     * @param {number} grid_width - UI Panel width
     * @param {String} projection_csv - Basename of CSV file containing projection data
     */
    constructor({
        plotly_config,
        ...args
    }, grid_width, projection_csv, target) {
        super(args);
        this.grid_width = grid_width;
        this.projection_csv = projection_csv;

        this.target = target ? target : `dataset_${this.primary_key}`;

        const config = plotly_config;

        let high_color;
        let mid_color;
        let low_color;
        // some cases where colors was not saved
        if (config.colors && Object.entries(config.colors).length !== 0) {
            low_color = config.colors.low_color;
            mid_color = config.colors.mid_color;
            high_color = config.colors.high_color;
        } else {
            // default
            low_color = '#e7d1d5';
            mid_color = null;
            high_color = '#401362';
        }

        this.low_color = low_color;
        this.mid_color = mid_color;
        this.high_color = high_color;

        this.clear_display();
        $(`#${this.target}`).append(this.template());
        this.fetch_svg_paths();
    }
    clear_display() {
        $(`#${this.target}_svg_cc`).remove();
        if (this.zoomed) $(`#${this.target}_svg_cc_zoomed`).remove();
    }

    clear_zoomed() {
        if (this.zoomed) $(`#${this.target}_svg_cc_zoomed`).remove();
    }
    /**
     * Get SVG data for coloring svg.
     * @param {string} gene_symbol - Gene symbol to visualize.
     */
    get_data(gene_symbol) {
        let url = `/api/plot/${this.dataset_id}/svg?gene=${gene_symbol}`;
        if (this.projection_csv) {
            url += `&projection_csv=${this.projection_csv}`;
        }
        return axios.get(url);
    }
    /**
     * Fetch the svg files from the server and cache them.
     */
    fetch_svg_paths() {
        const id = `#${this.target}_svg_c`;
        const svg = $(id);
        const snap = Snap(id);
        const imgSrc = svg.data('path');
        Snap.load(imgSrc, svg => {
            //console.log("Snap loaded:");
            //console.log(svg);
            const paths = svg.selectAll('path, circle, rect, ellipse');
            svgs[this.dataset_id] = paths;
            snap.append(svg);
        });
    }
    /**
     * Fetch the svg file from the server for zoomed.
     */
    fetch_svg_paths_and_draw_zoomed() {
        const id = `#${this.target}_svg_c_zoomed`;
        const svg = $(id);
        const snap = Snap(id);
        const imgSrc = svg.data('path');
        Snap.load(imgSrc, svg => {
            const paths = svg.selectAll('path, circle, rect, ellipse');
            this.zoomed_paths = paths;
            snap.append(svg);
            this.draw_chart(this.data, true);
        });
    }
    /**
     * Draw chart by coloring the svg paths. There are 3
     * ways to color, by gene, dataset, or tissue.
     * @param {Object} data - SVG data.
     */
    draw_chart(data, zoomed = false) {
        const score = data.scores[SCORING_METHOD];
        // for those fields which have no reading, a specific value is sometimes put in instead
        // These are colored a neutral color
        const NA_FIELD_PLACEHOLDER = -0.012345679104328156;
        const NA_FIELD_COLOR = '#808080';

        const paths = zoomed ? this.zoomed_paths : svgs[this.dataset_id];

        const {
            data: expression
        } = data;

        if (SCORING_METHOD === 'gene' || SCORING_METHOD === 'dataset') {
            const {
                min,
                max
            } = score;
            var color_range = null;

            // are we doing a three- or two-color gradient?
            if (this.mid_color) {
                if (min >= 0) {
                    // All values greater than 0, do right side of three-color
                    color_range = d3
                        .scaleLinear()
                        .domain([min, max])
                        .range([this.mid_color, this.high_color]);
                } else if (max <= 0) {
                    // All values under 0, do left side of three-color
                    color_range = d3
                        .scaleLinear()
                        .domain([min, max])
                        .range([this.low_color, this.mid_color]);
                } else {
                    // We have a good value range, do the three-color
                    color_range = d3
                        .scaleLinear()
                        .domain([min, 0, max])
                        .range([this.low_color, this.mid_color, this.high_color]);
                }
            } else {
                color_range = d3
                    .scaleLinear()
                    .domain([min, max])
                    .range([this.low_color, this.high_color]);
            }

            const color = color_range;
            const tissues = Object.keys(data.data);

            // Sometimes path isn't defined yet - latency/async issue??
            if (! paths) {
                return false;
            }

            paths.forEach(path => {
                const tissue_classes = path.node.className.baseVal.split(' ');

                tissue_classes.forEach(tissue => {
                    if (tissues.includes(tissue)) {
                        if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                            path.attr('fill', NA_FIELD_COLOR);
                        } else {
                            path.attr('fill', color(expression[tissue]));
                        }

                        if (!this.target.includes('modal')) {
                            const tooltip = d3
                                .select('#tip')
                                .attr('class', 'tooltip')
                                .style('opacity', 0);
                            path.mouseover(() => {
                                // TODO:
                                // Changing this toggle on index doesn't work
                                const math = $(`#math_menu_${this.id}`).attr('data-math');
                                let score;
                                // Apply math transformation to expression score
                                if (math == 'log2') {
                                    score = df.format('.2f')(Math.log2(expression[tissue]));
                                } else if (math == 'log10') {
                                    score = d3.format('.2f')(Math.log10(expression[tissue]));
                                } else {
                                    //math == 'raw'
                                    score = d3.format('.2f')(expression[tissue]);
                                }
                                const tmpl_tooltip = $.templates('#tmpl_tooltip');
                                const tooltip_html = tmpl_tooltip.render({
                                    tissue,
                                    score,
                                });
                                tooltip.html(tooltip_html);
                                $('#tip').show();
                                tooltip
                                    .transition()
                                    .duration(1)
                                    .style('opacity', 1);
                            });
                            path.mouseout(() => {
                                tooltip
                                    .transition()
                                    .duration(1)
                                    .style('opacity', 0);
                                $('#tip').hide();
                            });
                        }
                    }
                });
                // Draw the tooltip
            });
        } else {
            // tissues scoring
            const tissues = Object.keys(score);

            var color_range = {};

            if (this.mid_color) {
                tissues.forEach(tissue => {
                    let {
                        min,
                        max
                    } = score[tissue];

                    if (min >= 0) {
                        color_range[tissue] = d3
                            .scaleLinear()
                            .domain([min, max])
                            .range([this.mid_color, this.high_color]);
                    } else if (max <= 0) {
                        color_range[tissue] = d3
                            .scaleLinear()
                            .domain([min, max])
                            .range([this.low_color, this.mid_color]);
                    } else {
                        color_range[tissue] = d3
                            .scaleLinear()
                            .domain([min, 0, max])
                            .range([this.low_color, this.mid_color, this.high_color]);
                    }
                });
            } else {
                tissues.forEach(tissue => {
                    let {
                        min,
                        max
                    } = score[tissue];

                    color_range[tissue] = d3
                        .scaleLinear()
                        .domain([min, max])
                        .range([this.low_color, this.high_color]);
                });
            }

            const color = color_range;

            paths.forEach(path => {
                const tissue_classes = path.node.className.baseVal.split(' ');

                tissue_classes.forEach(tissue => {
                    if (tissue && color[tissue]) {
                        const color_scale = color[tissue];
                        path.attr('fill', color_scale(expression[tissue]));

                        if (!this.target.includes('modal')) {
                            const tooltip = d3
                                .select('#tip')
                                .attr('class', 'tooltip')
                                .style('opacity', 0);
                            path.mouseover(() => {
                                // TODO:
                                // Changing this toggle on index doesn't work
                                const math = $(`#math_menu_${this.id}`).attr('data-math');
                                let score;
                                // Apply math transformation to expression score
                                if (math == 'log2') {
                                    score = df.format('.2f')(Math.log2(expression[tissue]));
                                } else if (math == 'log10') {
                                    score = d3.format('.2f')(Math.log10(expression[tissue]));
                                } else {
                                    //math == 'raw'
                                    score = d3.format('.2f')(expression[tissue]);
                                }
                                const tmpl_tooltip = $.templates('#tmpl_tooltip');
                                const tooltip_html = tmpl_tooltip.render({
                                    tissue,
                                    score,
                                });
                                tooltip.html(tooltip_html);
                                $('#tip').show();
                                tooltip
                                    .transition()
                                    .duration(1)
                                    .style('opacity', 1);
                            });
                            path.mouseout(() => {
                                tooltip
                                    .transition()
                                    .duration(1)
                                    .style('opacity', 0);
                                $('#tip').hide();
                            });
                        }
                    }
                });
            });
        }
        this.hide_loading();
        this.show();
        if (!this.target.includes('modal')) this.draw_legend(data, zoomed);

        return true;
    }
    /**
     * Draw linear gradient as legend
     * @param {object} data - SVG data
     */
    draw_legend(data, zoomed = false) {
        let target;
        target = zoomed ? `dataset_${this.primary_key}_svg_cc_zoomed` : `${this.target}_svg_c`;
        const node = document.getElementById(target);
        // Create our legend svg
        const legend = d3
            .select(node)
            .append('svg')
            .style('position', 'absolute')
            .attr('width', '100%')
            .attr('class', 'svg-gradient-container');
        const defs = legend.append('defs');
        // Define our gradient shape
        const linear_gradient = defs
              .append('linearGradient')
              .attr('id', `${this.target}-linear-gradient${zoomed ? '_zoomed' : ''}`)
              .attr('x1', '0%')
              .attr('y1', '0%')
              .attr('x2', '100%')
              .attr('y2', '0%');

        // TODO: Issues to resolve here.  The 'atf4' gene in this datasets:
        //  The Adult Cochlea Response to PTS-Inducing Noise - Summary View
        // Has a data range of around -0.34 up but here the min is returning as
        // 0. Is this being converted to an into somewhere?
        const score = data.scores[SCORING_METHOD];
        const {
            min,
            max
        } = score;

        // Create the gradient points for either three- or two-color gradients
        if (this.mid_color) {
            // Even if a midpoint is called for, it doesn't make sense if the values are
            //  all less than or all greater than 0
            if (min >= 0) {
                linear_gradient
                    .append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', this.mid_color);
                linear_gradient
                    .append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', this.high_color);
            } else if (max <= 0) {
                linear_gradient
                    .append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', this.low_color);
                linear_gradient
                    .append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', this.mid_color);
            } else {
                // This means we've got a good distribution of min under 0 and max above
                //  it, so we can do a proper three-color range
                // midpoint offset calculation, so the mid color is at 0
                //var mid_offset = (1 - (min / max - min))*100;
                const mid_offset = (Math.abs(min)/(max + Math.abs(min)))*100;

                linear_gradient
                    .append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', this.low_color);
                linear_gradient
                    .append('stop')
                    .attr('offset', `${mid_offset}%`)
                    .attr('stop-color', this.mid_color);
                linear_gradient
                    .append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', this.high_color);
            }
        } else {
            linear_gradient
                .append('stop')
                .attr('offset', '0%')
                .attr('stop-color', this.low_color);
            linear_gradient
                .append('stop')
                .attr('offset', '100%')
                .attr('stop-color', this.high_color);
        }

        const {
            width
        } = node.getBoundingClientRect();
        // Draw they rectangle using the linear gradient
        legend
            .append('rect')
            .attr('width', width / 2)
            .attr('y', 10)
            .attr('x', width / 4)
            .attr('height', 10)
            .style(
                'fill',
                `url(#${this.target}-linear-gradient${zoomed ? '_zoomed' : ''}`
            );

        const xScale = d3
              .scaleLinear()
              .domain([min, max])
              .range([0, width / 2]);

        const xAxis = d3
            .axisBottom()
            .ticks(3)
            .scale(xScale);
        legend
            .append('g')
            .attr('class', 'axis')
            .attr('transform', `translate(${width / 4}, 20)`)
            .call(xAxis);
    }
    show() {
        $(`#${this.target}_svg_cc`).show();
        if (this.zoomed) $(`#${this.target}_svg_cc_zoomed`).show();
    }
    template() {
        const svg_class =
            this.grid_width == 8 ?
            'grid-width-8' :
            this.grid_width == 12 ?
            'grid-width-12' :
            '';
        const template = `
      <div
        id="${this.target}_svg_cc"
        style="position: relative; display:none; ${
          this.target.includes('modal') ? 'height:100%;' : ''
        }"
        class="h5ad-svg-container ${svg_class}">
        <div
          id="${this.target}_svg_c"
          class="svg-content" data-dataset-id="${this.dataset_id}"
          data-path="datasets_uploaded/${this.dataset_id}.svg">
        </div>
      </div>
    `;
        return template;
    }

    draw_zoomed() {
        this.clear_zoomed();
        $(`#dataset_zoomed div#${this.target}`).append(this.zoomed_template());
        this.fetch_svg_paths_and_draw_zoomed();
    }

    zoomed_template() {
        return `
        <div
        id="dataset_${this.primary_key}_svg_cc_zoomed"
        style='display:none'
        class="h5ad-svg-container-zoomed">
            <div
                id="dataset_${this.primary_key}_svg_c_zoomed"
                class="svg-content" data-dataset-id="${this.dataset_id}"
                data-path="datasets_uploaded/${this.dataset_id}.svg">
            </div>
        </div>
        `;
    }

    /**
     * Display warning above the plot
     */
     show_warning(msg) {

        const hover_msg = " Hover to see warning.";

        const dataset_selector = $( `#dataset_${this.primary_key}_svg_cc` );
        const warning_template = `
        <div class='dataset-warning bg-warning' id='dataset_${this.primary_key}_svg_warning'>
            <i class='fa fa-exclamation-triangle'></i>
            <span id="dataset_${this.primary_key}_svg_msg">${hover_msg}</span>
        </div>`;

        // Add warning below the chart
        // NOTE: must add to DOM before making selector variables
        dataset_selector.append(warning_template);

        const warning_selector = $( `#dataset_${this.primary_key}_svg_warning` );
        const msg_selector = $( `#dataset_${this.primary_key}_svg_msg` );

        // Add some CSS to warning to keep at top of container and not push display down
        warning_selector.css('position', 'absolute').css('bottom', '0').css('z-index', '2');

        // Add hover events (via jQuery)
        warning_selector.mouseover(() => {
            msg_selector.html(msg);
        });
        warning_selector.mouseout(() => {
            msg_selector.text(hover_msg);
       });
    }
}

/**
 * Class representing a TSNE display
 * @extends Display
 */
class TsneDisplay extends Display {
  /**
   * Initialize tSNE
   * This subclass of display takes an extra argument, gene symbol.
   * This is because we draw this one differently and don't query
   * the server for data, but query the server for the image.
   * @param {Object} Data - Display data
   * @param {string} gene_symbol - Gene symbol to visualize
   * @param {String} projection_csv - Basename of CSV file containing projection data
   */
  constructor({ plotly_config, ...args }, gene_symbol, projection_csv, target) {
    super(args);
    const config = plotly_config;
    this.gene_symbol = gene_symbol;
    this.analysis_id = config.analysis ? config.analysis.id : null;
    this.colors = JSON.stringify(config.colors);
    this.order = JSON.stringify(config.order);
    this.colorize_legend_by = config.colorize_legend_by;
    this.skip_gene_plot = config.skip_gene_plot;
    this.plot_by_group = config.plot_by_group;
    this.max_columns = config.max_columns;
    this.x_axis = config.x_axis;
    this.y_axis = config.y_axis;

    this.projection_csv = projection_csv;

    this.target = target ? target : `dataset_${this.primary_key}`;
  }
  /**
   * Get data for the tsne. This return the tsne image, and
   * is used mainly to check if the request is successful
   * before appending img to the DOM.
   * @param {string} gene_symbol - Gene symbol to visualize.
   */
   get_data(gene_symbol) {
      return axios.get(`/api/plot/${this.dataset_id}/tsne`, {
          params: {
              gene: gene_symbol,
              analysis: this.analysis_id,
              plot_type: this.plot_type,
              colorize_by: this.colorize_legend_by,
              skip_gene_plot: this.skip_gene_plot,
              plot_by_group: this.plot_by_group,
              max_columns: this.max_columns,
              x_axis: this.x_axis,
              y_axis: this.y_axis,
              analysis_owner_id: this.user_id,
              colors: this.colors,
              order: this.order,
              // helps stop caching issues
              timestamp: new Date().getTime(),
              projection_csv: this.projection_csv
          }
      });
  }

  draw_zoomed() {
    this.draw_zoomed_chart();
  }
  /**
   * Draw tSNE chart.
   */
  draw_chart() {
    this.clear_display();
    const target_div = `#${this.target}`;
    const target_div_img = `#${this.target}_tsne img`;
    $(target_div).append(this.template());
    $(target_div_img).attr('src', `data:image/png;base64,${this.data.image}`);
    this.hide_loading();
    this.show();
    return true;
  }

  // Essentially recycled the "draw_chart" function
  draw_zoomed_chart() {
    this.clear_display();
    const target_div = `#dataset_zoomed div#${this.target}`;
    const target_div_img = `${target_div}_tsne_zoomed img`;
    $(target_div).append(this.zoomed_template());
    $(target_div_img).attr('src', `data:image/png;base64,${this.data.image}`);
    this.hide_loading();
    this.show();
  }

  /**
   * HTML template representing the tSNE image
   */
    template() {
        return `
        <div id='${this.target}_tsne' class='img-static-container' style="position: relative;">
            <img style='max-width:96%; max-height:40em;'></img>
        </div>
        `;
    }

    zoomed_template() {
        return `
        <div id='${this.target}_tsne_zoomed' class='img-static-container'>
            <img style='max-width:96%; max-height:70em;'></img>
        </div>
        `;
    }

    clear_display() {
        $(`#${this.target}_tsne`).remove();
        if (this.zoomed) $(`#dataset_${this.primary_key}_tsne_zoomed`).remove();
    }

    /**
     * Display warning above the plot
     */
     show_warning(msg) {

        const hover_msg = " Hover to see warning."

        const dataset_selector = $( `#dataset_${this.primary_key}_tsne` );
        const warning_template = `
        <div class='dataset-warning bg-warning' id='dataset_${this.primary_key}_tsne_warning'>
            <i class='fa fa-exclamation-triangle'></i>
            <span id="dataset_${this.primary_key}_tsne_msg">${hover_msg}</span>
        </div>`;

        // Add warning above the chart
        // NOTE: must add to DOM before making selector variables
        dataset_selector.prepend(warning_template);

        const warning_selector = $( `#dataset_${this.primary_key}_tsne_warning` );
        const msg_selector = $( `#dataset_${this.primary_key}_tsne_msg` );

        // Add some CSS to warning to keep at top of container and not push display down
        warning_selector.css('position', 'absolute').css('z-index', '2');

        // Add hover events (via jQuery)
        warning_selector.mouseover(() => {
            msg_selector.html(msg);
        });
        warning_selector.mouseout(() => {
            msg_selector.text(hover_msg);
       });
    }
}

// *** General functions ***

// Get updates and additions to plot from the plot_display_config JS object
function get_plotly_updates(conf_area, plot_type, category) {
    let updates = {};
    for (const idx in conf_area) {
        const conf = conf_area[idx];
        // Get config (data and/or layout info) for the plot type chosen, if it exists
        if (conf.plot_type == "all" || conf.plot_type == plot_type) {
            const update = category in conf ? conf[category] : {};
            updates = {...updates, ...update};    // Merge updates
        }
    }
    return updates;
}
