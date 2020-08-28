/** Base class representing a display */
class Display {
  /**
   * Initialize display.
   * @param {Object} - Display data.
   */
  constructor({ id, dataset_id, user_id, label, plot_type }) {
    this.id = id;
    this.dataset_id = dataset_id;
    this.user_id = user_id;
    this.label = label;
    this.plot_type = plot_type;
    this.data = null;
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

    const panel_data = { ...dataset_panel };
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
          if (display && display.data) display.draw_chart(display.data);
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
    const { data } = await this.get_data(gene_symbol);
    if (data.success === -1) {
      this.show_error(data.message);
    } else {
      this.data = data;
      if (this.zoomed) {
        this.draw_zoomed();
      } else {
        this.draw_chart(data);
      }
    }
  }
  /**
   * Hides the display container
   */
  hide_loading() {
    $('#dataset_' + this.dataset_id + ' .dataset-status-container').hide();
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

  /**
   *  Display the error from the server.
   * @param {string} msg - Message returned back from server
   */
  show_error(msg) {
    $('#' + this.dataset_id + '_dataset_status_c h2').text(msg);
    $('#dataset_' + this.dataset_id + ' .dataset-status-container').show();
    $('#dataset_' + this.dataset_id + ' .plot-container').hide();
  }

  /**
   * Show the plot container.
   */
  show() {
    $(`#dataset_${this.dataset_id}_h5ad.plot-container`).show();
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
   */
  constructor({ plotly_config, ...args }, gene_symbol, target) {
    super(args);
    const config = plotly_config;
    this.gene_symbol = gene_symbol;
    this.econfig = config;
    this.extendRangeRatio = 10;

    var genes_track = plotly_config["tracks"]["EPIVIZ-GENES-TRACK"];
    if (genes_track.length > 0) {
      var gttrack = genes_track[0];
      if (gttrack.measurements) {
        this.genome = gttrack.measurements[0]["id"];
      }
      else {
        this.genome = gttrack.id[0]["id"];
      }
    }

    if (target) {
      this.target = target;
    } else {
      this.target = `dataset_${this.dataset_id}`;
    }

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
    $(`#dataset_zoomed div#dataset_${this.dataset_id}`).append(this.template());
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
      this.hide_loading();
      this.show();
    }
    else {
      // epiviz container already exists, so only update gneomic position in the browser
      const epiviznav = document.querySelector(`#${this.target}_epiviznav`);
      epiviznav.setAttribute("chr", this.data.chr);
      const nstart = this.data.start - Math.round((this.data.end - this.data.start) * this.extendRangeRatio);
      const nend = this.data.end + Math.round((this.data.end - this.data.start) * this.extendRangeRatio)
      epiviznav.setAttribute("start", nstart);
      epiviznav.setAttribute("end", nend);
      epiviznav.range = epiviznav.getGenomicRange(this.data.chr, nstart, nend);
      this.hide_loading();
      this.show();
    }
  }

  epiviz_template(config) {

    var epiviztemplate = "";
    for (const track in config.tracks) {
      const track_config = config.tracks[track];
      track_config.forEach(function(tc) {
        var temp_track = "<" + track + " slot='charts' ";
        if (Object.keys(tc).includes("id")){
          temp_track += " dim-s='" + JSON.stringify(tc.id) + "' ";
        } else {
          temp_track += " measurements='" + JSON.stringify(tc.measurements) + "' ";
        }

        if (tc.colors != null) {
          temp_track += " chart-colors='" + JSON.stringify(tc.colors) + "' ";
        }

        if (tc.settings != null) {
          temp_track += " chart-settings='" + JSON.stringify(tc.settings) + "' ";
        }

        temp_track += " style='min-height:200px;'></" + track + "> ";

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
          viewer=${"/epiviz.html?dataset_id=" + this.id + "&chr=" + this.data.chr + "&start=" + this.data.start + "&end=" + this.data.end}
        >
          ${this.epiviztemplate}
        </epiviz-navigation>
      </div>
    `;
  }

  clear_display() {
    // $(`#${this.target}_epiviz epiviz-navigation`).remove();
    if (this.zoomed) $(`#dataset_${this.dataset_id}_epiviz_zoomed`).remove();
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
   */
  constructor({ plotly_config, ...args }) {
    super(args);
    const {
      x_axis,
      y_axis,
      hide_x_labels,
      hide_y_labels,
      point_label,
      color_name,
      facet_row,
      facet_col,
      marker_size,
      jitter,
      colors,
      order,
      analysis,
      showlegend
    } = plotly_config;
    this.x_axis = x_axis;
    this.y_axis = y_axis;
    this.hide_x_labels = hide_x_labels;
    this.hide_y_labels = hide_y_labels;
    this.point_label = point_label;
    this.color_name = color_name;
    this.facet_row = facet_row;
    this.facet_col = facet_col;
    this.marker_size = marker_size;
    this.jitter = jitter;
    this.colors = colors;
    this.order = order;
    this.analysis = analysis;
    this.showlegend = showlegend;
  }
  clear_display() {
    $(`#dataset_${this.dataset_id}_h5ad`).remove();
    if (this.zoomed) $(`#dataset_${this.dataset_id}_h5ad_zoomed`).remove();
  }
  clear_zoomed() {
    $(`#dataset_${this.dataset_id}_h5ad_zoomed`).remove();
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
      point_label: this.point_label,
      hide_x_labels: this.hide_x_labels,
      hide_y_labels: this.hide_y_labels,
      color_name: this.color_name,
      facet_row: this.facet_row,
      facet_col: this.facet_col,
      marker_size: this.marker_size,
      jitter: this.jitter,
      colors: this.colors,
      order: this.order,
      analysis: this.analysis,
      showlegend: this.showlegend
    });
  }
  /**
   * Draw chart.
   * @param {*} data - Plotly data.
   */
  draw_chart(data) {
    this.clear_display();

    const target_div = `dataset_${this.dataset_id}_h5ad`;
    const { plot_json, plot_config } = data;

    this.hide_loading();
    $(`#dataset_${this.dataset_id}`).append(this.template());
    this.show();

    var layout_mods = {
      autosize: true,
      showlegend: this.showlegend,
      margin: { l: 20, t: 20 },
      height: 350,
    };

    // Overwrite plot layout and config values with custom ones from display
    var layout = {
      ...plot_json.layout,
      ...layout_mods,
    };
    // These are updated outside 'layout_mods' as unpacking the nested object
          // overwrites the x-axis and y-axis from 'plot_json_layout'
    layout.xaxis.automargin = true;
    layout.yaxis.automargin = true;

    var config_mods = {
        responsive: false,
    };

    const config = {
      ...plot_config,
      ...config_mods,
    };

    Plotly.newPlot(target_div, plot_json.data, layout, config);
  }

  draw_zoomed() {
    this.clear_display();
    const target_div = `dataset_${this.dataset_id}_h5ad_zoomed`;
    // const target_div = document
    //   .querySelector(`#dataset_zoomed div#dataset_${this.dataset_id}`);

    $(`#dataset_zoomed div#dataset_${this.dataset_id}`).append(
      this.zoomed_template()
    );

    this.hide_loading();
    const { plot_json, plot_config } = this.data;
    Plotly.newPlot(target_div, { ...plot_json, plot_config });
  }

  show() {
    $('#dataset_' + this.dataset_id + ' .plot-container').show();
  }
  template() {
    const template = `
        <div
          id='dataset_${this.dataset_id}_h5ad'
          class="h5ad-container">
        </div>
      `;
    return template;
  }
  zoomed_template() {
    const template = `
        <div
          style='height:70vh'
          id='dataset_${this.dataset_id}_h5ad_zoomed'
          class="h5ad-container">
        </div>
      `;
    return template;
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
   */
  constructor({ plotly_config, ...args }, grid_width, target) {
    super(args);
    this.grid_width = grid_width;

    if (target) {
      this.target = target;
    } else {
      this.target = `dataset_${this.dataset_id}`;
    }

    let config = plotly_config;

    let high_color, low_color;
    // some cases where colors was not saved
    if (config.colors && Object.entries(config.colors).length !== 0) {
      low_color = config.colors.low_color;
      high_color = config.colors.high_color;
    } else {
      // default
      low_color = '#e7d1d5';
      high_color = '#401362';
    }

    this.low_color = low_color;
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
    return axios.get(`/api/plot/${this.dataset_id}/svg?gene=${gene_symbol}`);
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
      let paths = svg.selectAll('path, circle');
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
      let paths = svg.selectAll('path, circle');
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

    let paths;
    if (zoomed) {
      paths = this.zoomed_paths;
    } else {
      paths = svgs[this.dataset_id];
    }

    const { data: expression } = data;

    if (SCORING_METHOD === 'gene' || SCORING_METHOD === 'dataset') {
      const { min, max } = score;
      const color = d3
        .scaleLinear()
        .domain([min, max])
        .range([this.low_color, this.high_color]);

      const tissues = Object.keys(data.data);
      paths.forEach(path => {
        const tissue = path.node.className.baseVal;
        if (tissues.includes(tissue)) {
          path.attr('fill', color(expression[tissue]));

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
        // Draw the tooltip
      });
    } else {
      // tissues scoring
      const tissues = Object.keys(score);
      const color = {};
      tissues.forEach(tissue => {
        let { min, max } = score[tissue];
        color[tissue] = d3
          .scaleLinear()
          .domain([min, max])
          .range([this.low_color, this.high_color]);
      });
      paths.forEach(path => {
        const tissue = path.node.className.baseVal;
        if (tissue && color[tissue]) {
          let color_scale = color[tissue];
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
    }
    this.hide_loading();
    this.show();
    if (!this.target.includes('modal')) this.draw_legend(data, zoomed);
  }
  /**
   * Draw linear gradient as legend
   * @param {object} data - SVG data
   */
  draw_legend(data, zoomed = false) {
    let target;
    if (zoomed) {
      target = `dataset_${this.dataset_id}_svg_cc_zoomed`;
    } else {
      target = `${this.target}_svg_cc`;
    }
    const node = document.getElementById(target);
    // Create our legend svg
    const legend = d3
      .select(node)
      .append('svg')
      .style('position', 'absolute')
      .attr('width', '100%');
    const defs = legend.append('defs');
    // Define our gradient shape
    const linear_gradient = defs
      .append('linearGradient')
      .attr('id', `${this.target}-linear-gradient${zoomed ? '_zoomed' : ''}`)
      .attr('x1', '0%')
      .attr('y1', '0%')
      .attr('x2', '100%')
      .attr('y2', '0%');
    // Create the two gradient points for
    // both the low and high colors
    linear_gradient
      .append('stop')
      .attr('offset', '0%')
      .attr('stop-color', this.low_color);
    linear_gradient
      .append('stop')
      .attr('offset', '100%')
      .attr('stop-color', this.high_color);
    const { width } = node.getBoundingClientRect();
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
    // Create the ticks
    const score = data.scores[SCORING_METHOD];
    const { min, max } = score;
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
    let svg_class =
      this.grid_width == 8
        ? 'grid-width-8'
        : this.grid_width == 12
        ? 'grid-width-12'
        : '';
    const template = `
      <div
        id="${this.target}_svg_cc"
        style="display:none; ${
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
    const template = `
    <div
      id="dataset_${this.dataset_id}_svg_cc_zoomed"
      style='display:none'
      class="h5ad-svg-container-zoomed">
      <div
        id="dataset_${this.dataset_id}_svg_c_zoomed"
        class="svg-content" data-dataset-id="${this.dataset_id}"
        data-path="datasets_uploaded/${this.dataset_id}.svg">
      </div>
    </div>
  `;
    return template;
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
   */
  constructor({ plotly_config, ...args }, gene_symbol, target) {
    super(args);
    const config = plotly_config;
    this.gene_symbol = gene_symbol;
    this.analysis_id = config.analysis_id;
    this.colors = JSON.stringify(config.colors);
    this.colorize_legend_by = config.colorize_legend_by;
    this.x_axis = config.x_axis;
    this.y_axis = config.y_axis;

    if (target) {
      this.target = target;
    } else {
      this.target = `dataset_${this.dataset_id}`;
    }
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
              colorize_by: this.colorize_legend_by,
              x_axis: this.x_axis,
              y_axis: this.y_axis,
              analysis_owner_id: this.user_id,
              colors: this.colors,
              // helps stop caching issues
              timestamp: new Date().getTime()
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
    $(target_div_img).attr('src', "data:image/png;base64," + this.data);
    this.hide_loading();
    this.show();
  }

  // Essentially recycled the "draw_chart" function
  draw_zoomed_chart() {
    this.clear_display();
    const target_div = `#dataset_zoomed div#${this.target}`;
    const target_div_img = `${target_div}_tsne_zoomed img`;
    $(target_div).append(this.zoomed_template());
    $(target_div_img).attr('src', "data:image/png;base64," + this.data);
    this.hide_loading();
    this.show();
  }

  /**
   * HTML template representing the tSNE image
   */
  template() {
   return `
      <div id='${this.target}_tsne' class='img-static-container'>
        <img style='max-width:96%; max-height:20em;'></img>
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
    if (this.zoomed) $(`#dataset_${this.dataset_id}_tsne_zoomed`).remove();
  }
}
