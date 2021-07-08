'use strict';

/*
  Requirements
  - JQuery
  - JSRender
  - gear/common.js (and that check_for_login() has been called.)
  - gear/classes/dataset.js
*/

class DatasetPanel extends Dataset {
  constructor({ grid_position, grid_width, ...args }) {
    super(args);
    this.grid_position = grid_position;
    this.grid_width = grid_width;
    this.config = null;
    this.displays = null;
    this.default_display_id = null;
    this.active_display_id = null;
    this.zoomed = false;
  }

  get_dataset_displays(user_id, dataset_id) {
    return $.ajax({
      url: './cgi/get_dataset_displays.cgi',
      type: 'POST',
      data: { user_id, dataset_id },
      dataType: 'json',
    });
  }

  get_dataset_display(display_id) {
    return $.ajax({
      url: './cgi/get_dataset_display.cgi',
      type: 'POST',
      data: { display_id },
      dataType: 'json',
    });
  }

  get_default_display(user_id, dataset_id) {
    return $.ajax({
      url: './cgi/get_default_display.cgi',
      type: 'POST',
      data: { user_id, dataset_id },
      dataType: 'json',
    });
  }

  get_embedded_tsne_display(dataset_id) {
    return $.ajax({
      url: './cgi/get_embedded_tsne_display.cgi',
      type: 'POST',
      data: { dataset_id },
      dataType: 'json',
    });
  }

  async draw_chart(gene_symbol, display_id) {
    const data = await this.get_dataset_display(display_id)

    let zoom = false;
    if (this.display && this.display.zoomed) {
      zoom = true;
    }

    let display;
    if (
      data.plot_type === 'bar' ||
      data.plot_type === 'scatter' ||
      data.plot_type === 'violin' ||
      data.plot_type === 'line' ||
      data.plot_type === 'contour' ||
      data.plot_type === 'tsne_dynamic' ||  // legacy
      data.plot_type === 'tsne/umap_dynamic'
    ) {
      display = new PlotlyDisplay(data);
    } else if (data.plot_type === 'svg') {
      display = new SVGDisplay(data, this.grid_width);
    } else if (data.plot_type === 'tsne' ||
      data.plot_type === 'tsne_static' ||
      data.plot_type === 'umap_static' ||
      data.plot_type === 'pca_static') {
      display = new TsneDisplay(data, gene_symbol);
    } else if (data.plot_type === 'epiviz') {
      display = new EpiVizDisplay(data, gene_symbol);
    }
    this.display = display;

    // We first draw with the display then zoom in, so whenever
    // gene search is updated, the other displays behind the zoomed
    // is updated too.
    display.draw(gene_symbol);
    if (zoom) display.zoom_in();
  }

  // async draw_mg_chart(gene_symbols, display_id)
  async draw_mg_chart(gene_symbols) {
    const data = {dataset_id: this.id, user_id: this.user_id, plot_type: "heatmap", plotly_config:{}};
    //const data = await this.get_dataset_display(display_id)

    let zoom = false;
    if (this.display && this.display.zoomed) {
      zoom = true;
    }
    let display;
    if (data.plot_type === 'heatmap' ||
      data.plot_type === 'mg_violin' ||
      data.plot_type === 'volcano'
    ) {
      display = new DashMGDisplay(data, gene_symbols)
    }
    this.display = display;

    // We first draw with the display then zoom in, so whenever
    // gene search is updated, the other displays behind the zoomed
    // is updated too.
    display.draw(gene_symbols);
    if (zoom) display.zoom_in();
  }

  // Draw single-gene plots
  async draw({ gene_symbol } = {}) {
    if (this.display) this.display.clear_display();

    this.show_loading();

    // cache gene_symbol so we can use it to redraw with different display
    this.gene_symbol = gene_symbol;

    if (this.display) {
      this.draw_chart(gene_symbol, this.display.id);

    } else {
      // first time searching gene and displays have not been loaded
      const { default_display_id } = await this.get_default_display(CURRENT_USER.id, this.id);

      if (default_display_id) {
        this.default_display_id = default_display_id;
        this.draw_chart(gene_symbol, default_display_id);

        // cache all owner/user displays for this panel;
        const owner_displays = await this.get_dataset_displays(this.user_id, this.id);
        this.owner_displays = owner_displays;

        const user_displays = await this.get_dataset_displays(CURRENT_USER.id, this.id);
        this.user_displays = user_displays;

        this.register_events();
      } else {
        // No default display, this really shouldn't happen because
        // owners should always have atleast done this after upload
        this.show_error(
          'No default display. Create one in the dataset curator.'
        );
      }
    }
  }

  // Draw multigene plots
  // NOTE: Currently hardcoded to draw a default heatmap plot
  async draw_mg({ gene_symbols } = {}) {
    if (this.display) this.display.clear_display();

    this.show_loading();

    // cache gene_symbol so we can use it to redraw with different display
    this.gene_symbols = gene_symbols;

    // NOTE: Forcing multigene charts to be drawn for now
    this.draw_mg_chart(gene_symbols);

    /* NOTE: No default displays saved for multigene plots yet.
    if (this.display) {
      this.draw_mg_chart(gene_symbols, this.display.id);
    } else {
      // first time searching gene and displays have not been loaded
      const { default_display_id } = await this.get_default_display(CURRENT_USER.id, this.id);

      if (default_display_id) {
        this.default_display_id = default_display_id;
        this.draw_mg_chart(gene_symbols, default_display_id);

        // cache all owner/user displays for this panel;
        const owner_displays = await this.get_dataset_displays(this.user_id, this.id);
        this.owner_displays = owner_displays;

        const user_displays = await this.get_dataset_displays(CURRENT_USER.id, this.id);
        this.user_displays = user_displays;

        this.register_events();
      } else {
        // No default display, this really shouldn't happen because
        // owners should always have atleast done this after upload
        this.show_error(
          'No default display. Create one in the dataset curator.'
        );
      }
    }
    */
  }

  async redraw(display_id) {
    // This check is here in case there was no
    // default display rendered, and a user tries
    // to toggle to a different display. There would
    // be no display to clear.
    if (this.display) this.display.clear_display();

    this.show_loading();
    this.draw_chart(this.gene_symbol, display_id);
  }

  register_events() {
    // redraw plot when user changes display
    $(`#dataset_${this.id}`).on('changePlot', (e, display_id) =>
      this.redraw(display_id)
    );

    // zoom event
    $(`#dataset_${this.id}_zoom_in_control`).click(event => {
      $(`#dataset_grid`).fadeOut(() => {
        this.display.zoom_in();
      });
    });

    // info event
    $(`#dataset_${this.id}_info_launcher`).click(() => {
      const panel = dataset_collection_panel.datasets.find(
        d => d.id == this.id
      );

      // get default display id for this display?

      const { id, title, ldesc, schematic_image } = panel;
      const infobox_tmpl = $.templates('#tmpl_infobox');
      const infobox_html = infobox_tmpl.render({
        dataset_id: id,
        title,
        ldesc,
        schematic_image,
      });
      $('#modals_c').html(infobox_html);
      $(`#dataset_${this.id}_info`).modal({ show: true });
    });

    // Displays modal events
    $(`#dataset_${this.id}_displays`).click(() => {
      const panel = dataset_collection_panel.datasets.find(
        d => d.id == this.id
      );
      const { id, title, user_displays, owner_displays } = panel;

      const displays_tmpl = $.templates('#tmpl_displays_box');
      const displays_html = displays_tmpl.render({
        dataset_id: id,
        title,
        user_displays,
        owner_displays,
        active_display_id: this.active_display_id || this.default_display_id,
      });

      $('#modals_c').html(displays_html);

      // change plot when user clicks on a different display
      const dataset_id = this.id;
      let display_panel = this;
      $('#modals_c div.modal-display').click(function(event) {
        $('#modals_c img.active-display').removeClass('active-display');
        const display_id = this.id.split('-').pop();
        display_panel.active_display_id = display_id;
        $(`#dataset_${dataset_id}`).trigger('changePlot', display_id);
        $(`#modals_c img#modal-display-img-${display_id}`).addClass(
          'active-display'
        );
      });

      user_displays.forEach(display => {
        // check if config has been stringified
        let gene_symbol;
        if (typeof display.plotly_config == 'string') {
          const config = JSON.parse(display.plotly_config);
          gene_symbol = config.gene_symbol;
        } else {
          const config = display.plotly_config;
          gene_symbol = config.gene_symbol;
        }

        if (
          display.plot_type === 'violin' ||
          display.plot_type === 'bar' ||
          display.plot_type === 'line' ||
          display.plot_type === 'scatter' ||
          display.plot_type === 'contour' ||
          display.plot_type === 'tsne_dynamic' ||  // legacy
          display.plot_type === 'tsne/umap_dynamic'
        ) {
          const d = new PlotlyDisplay(display);
          d.get_data(gene_symbol).then(({ data }) => {
            const { plot_json, plot_config } = data;
            Plotly.toImage(
              { ...plot_json, plot_config },
              { height: 500, width: 500 }
            ).then(url => {
              $(`#modal-display-img-${display.id}`).attr('src', url);
            });
          });
        } else if (display.plot_type === 'svg') {
          const target = `modal-display-${display.id}`;
          const d = new SVGDisplay(display, null, target);
          d.get_data(gene_symbol).then(({ data }) => {
            d.draw_chart(data);
          });
        } else {
          // tsne
          const target = `modal-display-${display.id}`;
          const d = new TsneDisplay(display, gene_symbol, target);
          d.draw(gene_symbol);
        }
      });
      owner_displays.forEach(display => {
        let config;
        if (typeof display.plotly_config == 'string') {
            config = JSON.parse(display.plotly_config);
        } else {
            config = display.plotly_config;
        }
        const { gene_symbol } = config;

        if (
          display.plot_type === 'violin' ||
          display.plot_type === 'bar' ||
          display.plot_type === 'line' ||
          display.plot_type === 'scatter' ||
          display.plot_type === 'contour' ||
          display.plot_type === 'tsne_dynamic' ||  // legacy
          display.plot_type === 'tsne/umap_dynamic'
        ) {
            const d = new PlotlyDisplay(display);
            d.get_data(gene_symbol).then(
                ({ data }) => {
                    const { plot_json, plot_config } = data;
                    Plotly.toImage(
                        { ...plot_json, plot_config },
                        { height: 500, width: 500 }
                    ).then(
                        url => { $(`#modal-display-img-${display.id}`).attr('src', url) }
                    ).then(
                        () => { $(`#modal-display-${display.id}-loading`).hide() }
                    );
                }
            );
        } else if (display.plot_type === 'svg') {
            const target = `modal-display-${display.id}`;
            const d = new SVGDisplay(display, null, target);
            d.get_data(gene_symbol).then(
                ({ data }) => {
                    d.draw_chart(data).then(
                        () => { $(`#modal-display-${display.id}-loading`).hide() }
                    )
                }
            );
        } else {
            // tsne
            const target = `modal-display-${display.id}`;
            const d = new TsneDisplay(display, gene_symbol, target);

            d.draw(gene_symbol);
            $(`#modal-display-${display.id}-loading`).hide();
        }
        // TODO: Add multigene plots as well, or break into new function
      });
      $(`#dataset_${this.id}_displays_modal`).modal({ show: true });
    });
  }

  show_error(msg) {
    $('#' + this.id + '_dataset_status_c h2').text(msg);
    $('#dataset_' + this.id + ' .dataset-status-container').show();
    $('#dataset_' + this.id + ' .plot-container').hide();
  }

  show_loading() {
    $('#' + this.id + '_dataset_status_c h2').text('Loading...');
    $('#dataset_' + this.id + ' .dataset-status-container').show();

    if (this.dtype === 'svg-expression') {
      // this.hide_Svg();
    } else if (
      this.dtype === 'image-static-standard' ||
      this.dtype === 'image-static'
    ) {
      $('#dataset_' + this.id + ' .image-static-container').hide();
    } else {
      $('#dataset_' + this.id + ' .plot-container').hide();
    } // TODO: Hide other container plot types...
  }

  show_no_match() {
    $('#' + this.id + '_dataset_status_c h2').text('Gene not found');
    $('#dataset_' + this.id + ' .dataset-status-container').show();
    $('#dataset_' + this.id + ' .plot-container').hide();
  }

  show_plot() {
    $('#dataset_' + this.id + ' .dataset-status-container').hide();
    $('#dataset_' + this.id + ' .plot-container').show();
  }

  show_no_default() {
    const lastForwardSlash = window.location.href.lastIndexOf('/');
    const baseURL = window.location.href.slice(0, lastForwardSlash);
    const link = `${baseURL}/dataset_curator.html#/dataset/${this.id}/displays`;
    $('#' + this.id + '_dataset_status_c h2').html(
      `<a href="${link}"> No default display for this dataset. <br> Click to curate how you\'d like it displayed</a>`
    );
    $('#dataset_' + this.id + ' .dataset-status-container').show();
    $('#dataset_' + this.id + ' .plot-container').hide();
  }
}
