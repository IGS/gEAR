'use strict';

/*
 This script relies on the source having also included the
 common.js within this project (for login purposes)
*/

/* global $ Vue */

Vue.use(Vuex);
Vue.use(VueRouter);
Vue.use(BootstrapVue);
Vue.use(BootstrapVueIcons);

// Install VeeValidate components globally
Vue.component("ValidationObserver", VeeValidate.ValidationObserver);
Vue.component("ValidationProvider", VeeValidate.ValidationProvider);

window.onload=() => {
  const datasetTitle = Vue.component("dataset-title", {
    template: `
      <b-row v-if='title' class='justify-content-md-center id="dataset-title"'>
        <div class="mt-5 col-12">
        <!--
          <h2 class='font-weight-light'>You are curating <span class='font-weight-bold'>{{ title }}</span></h2>
          <hr />
          -->
        </div>
      </b-row>
    `,
    computed: {
      ...Vuex.mapState(["dataset_id", "title"]),
    },
  });

  const addDisplayBtn = Vue.component("add-display-btn", {
    template: `
      <b-card
        @mouseover="isHovering = true"
        @mouseout="isHovering = false"
        @click="routeToAddDisplay"
        :class='classed'
      >
       <div class='mt-4 pt-5 card-body text-center id="add-display-btn"'>
         <div class='display-4'>
           <i class="fa fa-plus"></i>
         </div>
         <p class='card-text'>Add Display</p>
       </div>
     </b-card>
    `,
    data() {
      return {
        isHovering: false,
        displayCard: true,
      };
    },
    computed: {
      ...Vuex.mapState(["dataset_id"]),
      classed() {
        return {
          hovering: this.isHovering,
          "display-card": true,
          elevation: true,
          "border-0": true,
        };
      },
    },
    methods: {
      routeToAddDisplay() {
        this.$router.push(`displays/new`);
      },
    },
  });

  const tsneChart = Vue.component("tsne-chart", {
    template: `
     <div>
        <div v-if='is_loading' class='elevation border-0 mb-5 sticky-chart'>
          <img class='ml-5 m-5' style="width:50px;height:50px;" src='../img/loading_search.gif'></img>
        </div>
        <template v-if='plot_params_ready || preconfigured'>
          <img :id="preview_id"></img>
          <b-alert :show="success > 1 && !preconfigured" variant="warning" dismissible>
          {{ message }}
          </b-alert>

          <b-alert :show="success < 0" variant="danger" dismissible>
          {{ message }}
          </b-alert>
        </template>
      </div>
    `,
    props: {
      preconfigured: {
        default: false,
      },
      display_data: {},
    },
    data() {
      return {
        // local variables to consolidate Vuex version of these, and user/owner versions
        display_tsne_is_loading: true,
        display_image_data: null,
        // stolen from https://gist.github.com/gordonbrander/2230317
        preview_id: `tsne_preview_${Math.random().toString(36).substr(2, 9)}`,
      };
    },
    computed: {
      ...Vuex.mapState([
        "dataset_id",
        "analysis",
        "config",
        "image_data",
        "success",
        "message",
        "tsne_is_loading",
      ]),
      is_loading() {
        return (this.preconfigured && this.display_tsne_is_loading == true) ||
        this.tsne_is_loading;
      },
      plot_params_ready() {
        return this.config.x_axis && this.config.y_axis;
      },
    },
    created() {
      if (this.preconfigured) {
        this.display_tsne_is_loading = this.tsne_is_loading;
        this.display_image_data = this.image_data;
        this.draw_image();
      }
    },
    watch: {
      //TODO: Redo all this
      display_image_data() {
        if (this.display_image_data) {
          // applies to tSNE plots within user/owner display previews
          $(`#${this.preview_id}`).addClass("img-fluid");
        }
        $(`#${this.preview_id}`).css({ "max-height": "205px" });
        $(`#${this.preview_id}`).attr(
          "src",
          `data:image/png;base64,${this.display_image_data}`
        );
      },
      // This is the Vuex store and is for tsne plots within the DatasetDisplay component
      image_data() {
        if (!this.preconfigured && this.image_data) {
          $(`#${this.preview_id}`).addClass("img-fluid");
          $(`#${this.preview_id}`).attr(
            "src",
            `data:image/png;base64,${this.image_data}`
          );
        }
      },
    },
    methods: {
      async draw_image() {
        const { dataset_id } = this;
        const { plot_type } = this.display_data;
        const config = this.display_data.plotly_config;
        const { analysis } = config;
        const analysis_owner_id = this.display_data.user_id;

        // This has to be separate from "fetch_tsne_image" because in user/owner displays, different image data may be returned
        const payload = { ...config, plot_type, analysis, analysis_owner_id };
        const { data } = await axios.post(`/api/plot/${dataset_id}/tsne`, payload);

        this.display_tsne_is_loading = false;
        this.display_image_data = data.image;
      },
    },
  });

  const svgChart = Vue.component("svg-chart", {
    template: `
      <div>
        <div v-if='!display_data' class='elevation border-0 mb-5'>
          <b-card-body no-body
            class='elevation border-0 mb-5'>
            <div ref='chart'></div>
          </b-card-body>
          <b-alert :show="success > 1" variant="warning" dismissible>
          {{ message }}
          </b-alert>

          <b-alert :show="success < 0" variant="danger" dismissible>
          {{ message }}
          </b-alert>
        </div>
        <div v-else ref='chart' style='height:230px'></div>
      </div>
    `,
    props: {
      chart_data: {
        default: null,
      },
      low_color: {
        default: "",
      },
      mid_color: {
        default: "",
      },
      high_color: {
        default: "",
      },
      display_data: {
        default: null,
      },
    },
    data() {
      return {
        loading: false,
        paths: [],
        svg: {},
        scoring_method: "gene", // need to make this a toggle
        success: 0,
        message: "",
      };
    },
    computed: {
      ...Vuex.mapState(["dataset_id"]),
    },
    watch: {
      async svg(svg) {
        this.paths = svg.selectAll("path, circle");
        svg.select("svg").attr({
          width: "100%",
          height: this.display_data ? "200px" : "",
        });
        const snap = Snap(this.$refs.chart);
        snap.append(svg);

        if (this.display_data) {
          const { plotly_config } = this.display_data;
          const payload = { ...plotly_config };
          const { gene_symbol } = payload;
          const { data } = await axios.get(
            `/api/plot/${this.dataset_id}/svg?gene=${gene_symbol}`
          );
          this.color_svg(data);
        } else {
          this.color_svg();
        }
      },
      low_color() {
        this.color_svg();
      },
      mid_color() {
        this.color_svg();
      },
      high_color() {
        this.color_svg();
      },
      chart_data() {
        this.color_svg();
        this.success = this.chart_data.success;
        this.message = this.chart_data.message;
      },
    },
    async created() {
      const svg_path = `datasets_uploaded/${this.dataset_id}.svg`;
      Snap.load(svg_path, (svg) => {
        this.svg = svg;
      });

      if (this.display_data) {
        const { plotly_config } = this.display_data;
        const { gene_symbol } = { ...plotly_config };
        const { data } = await axios.get(
          `/api/plot/${this.dataset_id}/svg?gene=${gene_symbol}`
        );
        this.color_svg(data);
      }
    },
    methods: {
      color_svg(data) {
        let chart_data;
        let low_color;
        let mid_color;
        let high_color;
        if (data) {
          const { plotly_config } = this.display_data;
          const { colors } = { ...plotly_config };
          low_color = colors.low_color;
          mid_color = colors.mid_color;
          high_color = colors.high_color;
          chart_data = data;
        } else {
          chart_data = this.chart_data;
          low_color = this.low_color;
          mid_color = this.mid_color;
          high_color = this.high_color;
        }

        const score = chart_data.scores[this.scoring_method];
        const paths = this.paths;
        const { data: expression } = chart_data;
        if (
          this.scoring_method === "gene" ||
          this.scoring_method === "dataset"
        ) {
          const { min, max } = score;
          let color = null;
          // are we doing a three- or two-color gradient?
          if (mid_color) {
              if (min >= 0) {
                  // All values greater than 0, do right side of three-color
                  color = d3
                      .scaleLinear()
                      .domain([min, max])
                      .range([mid_color, high_color]);
              } else if (max <= 0) {
                  // All values under 0, do left side of three-color
                  color = d3
                      .scaleLinear()
                      .domain([min, max])
                      .range([low_color, mid_color]);
              } else {
                  // We have a good value range, do the three-color
                  color = d3
                      .scaleLinear()
                      .domain([min, 0, max])
                      .range([low_color, mid_color, high_color]);
              }
          } else {
              color = d3
                  .scaleLinear()
                  .domain([min, max])
                  .range([low_color, high_color]);
          }

          const tissues = Object.keys(chart_data.data);

          paths.forEach((path) => {
            const tissue = path.node.className.baseVal;
            if (tissues.includes(tissue)) {
              path.attr("fill", color(expression[tissue]));
            }
          });
        }
      },
    },
  });

  const plotlyChart = Vue.component("plotly-chart", {
    template: `
      <div>
        <div v-if='!display_data' class='elevation border-0 mb-5 sticky-chart'>
          <div v-if='!is_there_data'>
            <img class='ml-5 m-5' style="width:50px;height:50px;" src='../img/loading_search.gif'></img>
          </div>
          <div v-else ref='chart'>
            <img class='img-fluid' v-if='img' :src='imgData'></img>
          </div>
          <b-alert :show="success > 1" variant="warning" dismissible>
          {{ message }}
          </b-alert>

          <b-alert :show="success < 0" variant="danger" dismissible>
          {{ message }}
          </b-alert>
        </div>
        <div v-else ref='chart'>
          <div v-if='!imgData' class='col align-middle text-center mt-5 pt-4'>
            <p class=''>Loading</p>
          </div>
          <div v-else>
            <img class='img-fluid' v-if='img' :src='imgData'></img>
          </div>
        </div>
      </div>
    `,
    props: {
      data: {
        default: null,
      },
      display_data: {
        default: null,
      },
      img: Boolean,
    },
    data() {
      return {
        loading: false,
        imgData: "",
        success: 0,
        message: ""
      };
    },
    async created() {
      if (this.display_data) {
        this.loading = true;
        const { plotly_config, plot_type } = this.display_data;
        const payload = { ...plotly_config, plot_type };
        const { data } = await axios.post(
          `/api/plot/${this.dataset_id}`,
          payload
        );
        this.loading = false;
        this.draw_chart(data);
      }
    },
    mounted() {
      if (this.is_there_data) {
        this.draw_chart();
      }
    },
    computed: {
      ...Vuex.mapState([
        "dataset_id"
      , "plot_type"
    ]),
      is_there_data() {
        if (this.data === null) return false;
        return (
          Object.entries(this.data).length !== 0 &&
          this.data.constructor === Object
        );
      },
    },
    watch: {
      data() {
        // because data isn't included in template it
        // is not reactive. So here we explicitly watch
        // for changes and redraw
        this.draw_chart();
      },
    },
    methods: {
      draw_chart(data) {
        // Reset some params in case of failures
        this.imgData = '';
        this.loading = true;

        // Reset the div to plot a new image
        if (!this.img) {
          (this.$refs.chart).innerHTML = "";
        }

        const { plot_json } = data ? data : this.data;
        if (data) {
          this.success = data.success;
          this.message = data.message;
        } else {
          this.success = this.data.success;
          this.message = this.data.message;
        }
        if (this.success >= 1 ) {
          if (this.img) {
            Plotly.toImage({ ...plot_json, ...{static_plot:true} }).then((url) => {
              this.imgData = url;
            });
          } else {
            const curator_conf = post_plotly_config.curator;
            const plot_config = this.get_plotly_updates(curator_conf, this.plot_type, "config");
            Plotly.newPlot(this.$refs.chart, plot_json.data, plot_json.layout, plot_config);
            // Update plot with custom plot config stuff stored in plot_display_config.js
            const update_layout = this.get_plotly_updates(curator_conf, this.plot_type, "layout")
            Plotly.relayout(this.$refs.chart, update_layout)
          }
        }
        this.loading = false;
      },
      get_plotly_updates(conf_area, plot_type, category) {
        // Get updates and additions to plot from the plot_display_config JS object
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
    },
  });

  const userDisplays = Vue.component("user-displays", {
    template: `
      <b-row class='mt-4 justify-content-md-center'>
        <div class='col-12'>
          <h5 class=''>Your Displays</h5>
          <div class='row'>
            <div class=''>
            <add-display-btn class='m-2'></add-display-btn>
            </div>
            <div
              v-for="display_data in user_displays"
              :key='display_data.id'
              class='display-card m-2'
            >
              <transition name="fade" mode="out-in">
                <b-card
                  class='display-card elevation border-0'
                  footer-bg-variant="white"
                  header-bg-variant="white"
                  body-class='border-0 p-0 col'
                  header-class='border-0 pb-0'
                  footer-class='border-0 pt-0'
                >
                  <div slot='header'>
                  <p class='float-right'>
                    <b-badge
                      pill
                      class='purple'
                      variant="dark">{{ display_data.plotly_config.gene_symbol }}
                    </b-badge>
                    <b-badge
                      pill
                      class='purple'
                      v-if='display_data.is_default' variant="dark">Default
                    </b-badge>
                  </p>
                  </div>
                  <template v-if="display_data.plot_type === 'svg'">
                    <svg-chart :display_data='display_data'></svg-chart>
                  </template>
                  <template v-else-if="is_type_tsne(display_data)">
                    <tsne-chart
                      :display_data='display_data'
                      :preconfigured=true
                    ></tsne-chart>
                  </template>
                  <template v-else>
                    <plotly-chart :display_data='display_data' :img='true'></plotly-chart>
                  </template>
                  <div slot='footer' class='text-right'>
                    <b-button
                      style='color:red'
                      size='sm'
                      variant='link'
                      @click="delete_display(display_data.id)"
                    >Delete
                    </b-button>
                    <b-button
                      size='sm'
                      variant='link'
                      style='color:#562a6f'
                      @click='save_as_default(display_data.id)'
                    >Make Default</b-button>
                    <b-button
                      style='color:#562a6f'
                      variant='link'
                      size='sm'
                      @click="edit_display(display_data.id)"
                    >Edit
                    </b-button>
                  </div>
                </b-card>
            </transition>

            </div>
          </div><!-- end .row -->
        </div> <!-- end .col12 -->
      </b-row>
    `,
    components: {
      plotlyChart,
      svgChart,
      tsneChart,
      addDisplayBtn,
    },
    data() {
      return {
        // deleteModal: false,
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState(["user", "dataset_id", "default_display_id", "config"]),
      ...Vuex.mapGetters(["user_displays"]),
      get_config_analysis_id() {
        return typeof this.config.analysis == "undefined" ||
        this.config.analysis == null ? null : this.config.analysis.id;
      },
      user_id() {
        return this.user.id;
      },
    },
    methods: {
      ...Vuex.mapActions([
        "fetch_user_displays",
        "fetch_default_display",
        "remove_display",
        "update_default_display_id",
      ]),
      is_type_tsne(display_data) {
        return (
          display_data.plot_type === "tsne_static" ||
          display_data.plot_type === "umap_static" ||
          display_data.plot_type === "pca_static" ||
          display_data.plot_type === "tsne"
        );
      },
      edit_display(display_id) {
        // this.$router.replace(`/dataset/${this.dataset_id}/displays/${display_id}/edit`)
        this.$router.push(`displays/${display_id}/edit`);
      },
      get_default_display() {
        const user_id = this.user_id;
        const dataset_id = this.dataset_id;

        return $.ajax({
          url: "./cgi/get_default_display.cgi",
          type: "POST",
          data: { user_id, dataset_id },
          dataType: "json",
        });
      },
      async save_as_default(display_id) {
        const payload = {
          display_id,
          user_id: this.user.id,
          dataset_id: this.dataset_id,
        };
        await $.ajax({
          url: "./cgi/save_default_display.cgi",
          type: "POST",
          data: payload,
          dataType: "json",
        });

        this.update_default_display_id({ display_id });
      },
      async delete_display(display_id) {
        const payload = {
          id: display_id,
          user_id: this.user_id,
        };

        const res = await $.ajax({
          url: "./cgi/delete_dataset_display.cgi",
          type: "POST",
          data: payload,
          dataType: "json",
        });

        if (res.success) {
          // update our displays so it is removed from the
          // cards that are currently rendered
          this.remove_display({ display_id });
        }
      },
    },
    created() {
      if (this.dataset_id) {
        const user_id = this.user_id;
        const dataset_id = this.dataset_id;

        this.fetch_default_display({ user_id, dataset_id });
        this.fetch_user_displays({ user_id, dataset_id });
      }
    },
    beforeMount() {
      if (this.dataset_id) {
        const user_id = this.user_id;
        const dataset_id = this.dataset_id;

        this.fetch_default_display({ user_id, dataset_id });
        this.fetch_user_displays({ user_id, dataset_id });
      }
    },
  });
  const ownerDisplays = Vue.component("owner-displays", {
    template: `
      <b-row class='mt-4 justify-content-md-center'>
        <div class='col-12'>
          <h5 class=''>Owner's Displays</h5>
          <div class='row'>
            <div
              v-for="display_data in owner_displays"
              :key='display_data.id'
              class='display-card m-2'
            >
              <b-card
                class='display-card elevation border-0'
                footer-bg-variant="white"
                header-bg-variant="white"
                body-class='border-0 p-0 col'
                header-class='border-0 pb-0'
                footer-class='border-0 pt-0'
              >
                <div slot='header'>
                <p class='float-right'>
                  <b-badge
                    pill
                    class='purple'
                    variant="dark">{{ display_data.plotly_config.gene_symbol }}
                  </b-badge>
                  <b-badge
                    pill
                    class='purple'
                    v-if='display_data.is_default' variant="dark">Default
                  </b-badge>
                </p>
                </div>
                <template v-if="display_data.plot_type === 'svg'">
                  <svg-chart :display_data='display_data'></svg-chart>
                </template>
                <template v-else-if="is_type_tsne(display_data)">
                  <tsne-chart
                    :display_data='display_data'
                    :preconfigured=true
                  ></tsne-chart>
                </template>
                <template v-else>
                  <plotly-chart :display_data='display_data' :img='true'></plotly-chart>
                </template>
                <div slot='footer' class='text-right'>
                  <b-button
                    size='sm'
                    variant='link'
                    style='color:#562a6f'
                    @click='save_as_default(display_data.id)'
                  >Make Default
                  </b-button>
                </div>
              </b-card>
            </div>
          </div>
        </div>
      </b-row>
    `,
    components: {
      plotlyChart,
      svgChart,
      tsneChart,
    },
    data() {
      return {
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState([
        "user",
        "owner_id",
        "dataset_id",
        "default_display_id",
        "config",
      ]),
      ...Vuex.mapGetters(["is_user_owner", "owner_displays"]),
    },
    methods: {
      ...Vuex.mapActions([
        "fetch_owner_displays",
        "fetch_default_display",
        "update_default_display_id",
      ]),
      async save_as_default(display_id) {
        const payload = {
          display_id,
          user_id: this.user.id,
          dataset_id: this.dataset_id,
        };

        await $.ajax({
          url: "./cgi/save_default_display.cgi",
          type: "POST",
          data: payload,
          dataType: "json",
        });

        this.update_default_display_id({ display_id });
      },
      is_type_tsne(display_data) {
        return (
          display_data.plot_type === "tsne_static" ||
          display_data.plot_type === "umap_static" ||
          display_data.plot_type === "pca_static" ||
          display_data.plot_type === "tsne"
        );
      },
    },
    created() {
      if (this.dataset_id) {
        const owner_id = this.owner_id;
        const dataset_id = this.dataset_id;

        this.fetch_owner_displays({ owner_id, dataset_id });
      }
    },
  });

  const datasetCurator = Vue.component("dataset-curator", {
    template: `
      <b-container fluid>
        <dataset-title></dataset-title>
        <b-container fluid>
          <transition name="fade" mode="out-in">
            <router-view></router-view>
          </transition>
        </b-container>
      </b-container>
    `,
    props: ["dataset_id", "user"],
    components: {
      datasetTitle,
    },
    watch: {
      user() {
        location.reload();
      },
      dataset_id() {
        // not elegant, but watch for ids to change
        // and force a reload...
        location.reload();
      },
    },
    created() {
      this.fetch_dataset_info(this.dataset_id);
    },
    methods: {
      ...Vuex.mapActions(["fetch_dataset_info"]),
    },
  });

  const datasetDisplays = Vue.component("dataset-displays", {
    template: `
      <div id='dataset-displays'>
        <user-displays></user-displays>
        <owner-displays v-if='owner_id && !is_user_owner'></owner-displays>
      </div>
    `,
    computed: {
      ...Vuex.mapState(["owner_id"]),
      ...Vuex.mapGetters(["is_user_owner"]),
    },
    components: {
      userDisplays,
      ownerDisplays,
    },
  });

  const chooseDisplayType = Vue.component("choose-display-type", {
    template: `
      <div id="choose-display-type">
        <b-row>
        <b-col>
          <b-icon class="float-right" v-b-tooltip.hover="'Select plot type.  Some choices may not be available depending on observations in dataset.'" icon="question-circle-fill"></b-icon>
          <h3>Display Type</h3>
            <b-form-group>
              <b-form-select
                :disabled='loading || available_plot_types.success === -1 || display_options.length === 0'
                :value="plot_type"
                @input="update_plot_type($event)"
                :options="display_options"
              >
                <option v-if='loading' slot="first" :value="null">Loading available displays...</option>
                <option v-else-if='available_plot_types.success === -1' slot="first">{{ available_plot_types.message }}</option>
                <option v-else-if='display_options.length === 0' :value='null'>No available display types</option>
                <option v-else slot="first" :value="null">Choose display type</option>
              </b-form-select>
            </b-form-group>
          </b-col>
          </b-row>
          <hr v-if='plot_type'>
      </div>
    `,
    props: ['analysis_id'],
    computed: {
      ...Vuex.mapState([
        "user",
        "plot_type",
        "dataset_id",
        "available_plot_types",
        "dataset_type",
      ]),
      display_options() {
        const display_options = Object.keys(this.available_plot_types)
          .filter((type) => this.available_plot_types[type])
          .map((type, i) => type);

         return display_options;
      },
      loading() {
        return (
          Object.entries(this.available_plot_types).length === 0 &&
          this.available_plot_types.constructor === Object
        );
      },
    },
    created() {
      const user_id = this.user.id;
      const session_id = this.user.session_id;
      const dataset_id = this.dataset_id;
      const analysis_id = this.analysis_id;

      this.fetch_available_plot_types({ user_id, session_id, dataset_id, analysis_id });
    },
    methods: {
      ...Vuex.mapActions(["fetch_available_plot_types", "set_plot_type"]),
      update_plot_type(plot_type) {
        this.set_plot_type(plot_type);
        this.$emit("input", plot_type);
      },
    },
    watch: {
      analysis_id(analysis_id) {
          // When analysis ID changes, update plot types dropdown list
          const user_id = this.user.id;
          const session_id = this.user.session_id;
          const dataset_id = this.dataset_id;

          this.fetch_available_plot_types({ user_id, session_id, dataset_id, analysis_id });
      },
    }
  });

  const verticalLine = Vue.component("vertical-line", {
    template: `
      <b-input-group size="sm" prepend="x =">
        <b-form-input type="number" v-model="vl_pos" @keyup="$emit('update:vl_pos', vl_pos);">
        </b-form-input>
        <b-input-group-append>
          <b-form-select size="sm" :options="options" v-model="vl_style" @change="$emit('update:vl_style', vl_style);">
          </b-form-select>
        </b-input-group-append>
      </b-input-group>
    `,
    props: { vl: Object },
    data() {
      return {
        vl_pos: null,
        vl_style: "solid",
        options: [
          { value: "solid", text: "Solid" },
          { value: "dash", text: "Dashed" },
          { value: "dot", text: "Dotted" },
          { value: "dashdot", text: "Dash/Dot" },
        ],
      };
    },
  });

  const PlotlyArguments = Vue.component("plotly-arguments", {
    template: `
      <div id="plotly-arguments">
        <b-row>
          <b-col>
          <h3>Display Parameters</h3>
          <b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">X</label>
                <!-- NOTE: b-icon has a 'title' property which conflicts with the 'title' property v-b-tooltip uses -->
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for X axis.'" icon="question-circle-fill"></b-icon>
                <b-form-select :options='columns' v-model='x_axis' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
                <b-form-checkbox v-model='hide_x_labels'>
                  Hide X Tickmarks
                </b-form-checkbox>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Y</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for Y axis.'" icon="question-circle-fill"></b-icon>
                <b-form-select :options='columns' v-model='y_axis' size='sm'>
                  <template slot="first">
                  <option value="raw_value">expression</option>
                  </template>
                </b-form-select>
                <b-form-checkbox v-model='hide_y_labels'>
                  Hide Y Tickmarks
                </b-form-checkbox>
              </b-form-group>
              <b-form-group label-align-sm="right" v-if="['contour'].includes(plot_type)">
                <label class="mb-0">Z</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for Z axis.'" icon="question-circle-fill"></b-icon>
                <b-form-select :options='columns' v-model='z_axis' size='sm'>
                  <template slot="first">
                  <option value="raw_value">expression</option>
                  </template>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Color</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for coloring'" icon="question-circle-fill"></b-icon>
                <b-form-select :options='columns' v-model='color_name' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                <option value="raw_value">expression</option>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right" v-if="['scatter', 'tsne_dynamic', 'tsne/umap_dynamic'].includes(plot_type)">
                <label class="mb-0">Size</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for marker size'" icon="question-circle-fill"></b-icon>
                <b-form-select :options='columns' v-model='size_by_group' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                <option value="raw_value">expression</option>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right" v-if="plot_type !== 'contour'">
                <label class="mb-0">Label</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for tooltip when hovering over a data marker'" icon="question-circle-fill"></b-icon>
                <b-form-select :options="columns" v-model='point_label' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Facet Row</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for row facets'" icon="question-circle-fill"></b-icon>
                <b-form-select :options="Object.keys(levels)" v-model='facet_row' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Facet Column</label>
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for column facets'" icon="question-circle-fill"></b-icon>
                <b-form-select :options="Object.keys(levels)" v-model='facet_col' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>

              <!-- 'advanced options' collapsable -->
              <b-form-group>
              <b-button block v-b-toggle.adv-collapse variant="outline-secondary">Advanced Options</b-button>
              </b-form-group>
              <b-collapse id="adv-collapse">
                <b-form-group label-align-sm="right" v-if="!(hide_x_labels || x_axis in levels)">
                  <label class="mb-0">X Tick Range</label>
                  <b-icon class="float-right" icon="question-circle-fill"
                    v-b-tooltip.hover="'Set X-axis range of plot.  Both numbers must be filled in.  Not applicable for categorical data series.'">
                  </b-icon>
                  <validation-observer slim>
                  <b-input-group class="flex-nowrap" size="sm" prepend="min" append="max">
                    <validation-provider name="x-min" vid="xmin_value" rules="required_if:xmax_value" v-slot="{ errors }">
                      <b-form-input class="pr-0" type="number" v-model="x_min" :state="errors[0] ? false : null">
                      </b-form-input>
                      <b-form-invalid-feedback>{{ errors[0] }}</b-form-invalid-feedback>
                    </validation-provider>
                    <span class="mt-2"> - </span>
                    <validation-provider name="x-max" vid="xmax_value" rules="required_if:xmin_value" v-slot="{ errors }">
                      <b-form-input class="pr-0" type="number" v-model="x_max" :state="errors[0] ? false : null">
                      </b-form-input>
                      <b-form-invalid-feedback>{{ errors[0] }}</b-form-invalid-feedback>
                    </validation-provider>
                  </b-input-group>
                  </validation-observer>
                </b-form-group>

                <b-form-group label-align-sm="right">
                  <label class="mb-0">X Axis Title</label>
                  <b-icon class="float-right" v-b-tooltip.hover="'Title of axis'" icon="question-circle-fill"></b-icon>
                  <b-form-input v-model="x_title"></b-form-input>
                </b-form-group>

                <!-- SAdkins TOOD: Try to have the invalid-feedback outside of b-input-group for just one error msg (also in x_axis) -->
                <b-form-group label-align-sm="right" v-if="!(hide_y_labels || y_axis in levels)">
                  <label class="mb-0">Y Tick Range</label>
                  <b-icon class="float-right" icon="question-circle-fill"
                    v-b-tooltip.hover="'Set Y-axis range of plot.  Both numbers must be filled in.  Not applicable for categorical data series.'">
                  </b-icon>
                  <validation-observer slim>
                  <b-input-group class="flex-nowrap" size="sm" prepend="min" append="max">
                    <validation-provider name="y-min" vid="ymin_value" rules="required_if:ymax_value" v-slot="{ errors }">
                      <b-form-input class="pr-0" type="number" v-model="y_min" :state="errors[0] ? false : null">
                      </b-form-input>
                      <b-form-invalid-feedback>{{ errors[0] }}</b-form-invalid-feedback>
                    </validation-provider>
                    <span class="mt-2"> - </span>
                    <validation-provider name="y-max" vid="ymax_value" rules="required_if:ymin_value" v-slot="{ errors }">
                      <b-form-input class="pr-0" type="number" v-model="y_max" :state="errors[0] ? false : null">
                      </b-form-input>
                      <b-form-invalid-feedback>{{ errors[0] }}</b-form-invalid-feedback>
                    </validation-provider>
                  </b-input-group>
                  </validation-observer>
                </b-form-group>

                <b-form-group label-align-sm="right">
                  <label class="mb-0">Y Axis Title</label>
                  <b-icon class="float-right" v-b-tooltip.hover="'Title of axis'" icon="question-circle-fill"></b-icon>
                  <b-form-input v-model="y_title"></b-form-input>
                </b-form-group>

                <b-form-group label-align-sm="right" v-if="['scatter', 'tsne_dynamic', 'tsne/umap_dynamic'].includes(plot_type)">
                  <label class="mb-0">Marker Size (sm <--> lg)</label>
                  <b-icon class="float-right" icon="question-circle-fill"
                    v-b-tooltip.hover="'Set a constant marker size.  If Size series is also selected, this will set the minimum marker size.'">
                  </b-icon>
                  <b-form-input type="range" v-model="marker_size" min="1" max="15" number></b-form-input>
                </b-form-group>

                <b-form-group label-align-sm="right" v-if="['scatter', 'violin'].includes(plot_type)">
                  <b-form-checkbox v-model='jitter'>
                    <label class="mb-0" id="jitter">Jitter points</label>
                  </b-form-checkbox>
                  <b-tooltip target="jitter" triggers="hover">
                    <span v-if="plot_type === 'scatter'">Check to convert scatter plot into strip plot</span>
                    <span v-if="plot_type === 'violin'">Check to convert violin plot into beeswarm plot</span>
                  </b-tooltip>
                </b-form-group>

                <b-form-group label-align-sm="right">
                  <b-form-checkbox v-model='hide_legend'>
                    <label class="mb-0" id="hide_legend">Hide Legend</label>
                  </b-form-checkbox>
                  <b-tooltip target="hide_legend" triggers="hover">
                    Check to not display legend in plot
                  </b-tooltip>
                </b-form-group>

                <!-- vertical lines -->
                <b-form-group label-align-sm="right" v-if="(!(x_axis in levels) && ['scatter'].includes(plot_type))">
                  <label class="mb-0">Vertical Line Information</label>
                  <b-icon class="float-right" icon="question-circle-fill"
                    v-b-tooltip.hover="'Add vertical lines to plot with various properties.  Not applicable if X-axis is categorical.'">
                  </b-icon>
                  <vertical-line v-for="(vl, index) in vlines" v-bind.sync="vl" v-bind:key="vl.id"></vertical-line>
                  <div>
                  <b-icon icon="plus-circle-fill" @click="addRow"></b-icon>
                  <b-icon icon="dash-circle-fill" @click="removeLast" v-if="vlines.length"></b-icon>
                  </div>
                </b-form-group>
              </b-collapse>

              <b-col>
                <b-button
                  @click='preview'
                  :disabled='x_axis === null'
                  class='btn-purple float-right'
                  size='sm'>
                    Preview Chart
                </b-button>
              </b-col>
            </b-form-group>
          </b-col>
        </b-row>
        <hr>
      </div>
    `,
    data() {
      return {
        x_axis: null,
        y_axis: "raw_value",
        z_axis: "raw_value",
        point_label: null,
        color_name: null,
        facet_row: null,
        facet_col: null,
        size_by_group: null, // Marker size is based on a group.
        marker_size: 3, // if size_by_group set, then this will be min marker size
        jitter: false,
        hide_x_labels: false,
        hide_y_labels: false,
        hide_legend: false,
        x_min: null,
        x_max: null,
        y_min: null,
        y_max: null,
        x_title: null,
        y_title: null,

        vlines: [],
        vline_counter: 1,
      };
    },
    computed: {
      ...Vuex.mapState([
        "dataset_id",
        "config",
        "columns",
        "levels", // Can use to determine categorical series
        "plot_type",
      ]),
    },
    created() {
      if ("x_axis" in this.config) this.x_axis = this.config.x_axis;
      if ("y_axis" in this.config && this.config.y_axis) {
        // we only want to change y_axis if there's a value other than null so
        // it defaults to "raw_value" set in data above
        this.y_axis = this.config.y_axis;
      }
      if ("z_axis" in this.config && this.config.z_axis) {
        // we only want to change z_axis if there's a value other than null so
        // it defaults to "raw_value" set in data above
        this.z_axis = this.config.z_axis;
      }
      if ("hide_x_labels" in this.config)
        this.hide_x_labels = this.config.hide_x_labels;
      if ("hide_y_labels" in this.config)
        this.hide_y_labels = this.config.hide_y_labels;
      if ("hide_legend" in this.config)
        this.hide_legend = this.config.hide_legend;
      if ("color_name" in this.config) this.color_name = this.config.color_name;
      if ("facet_row" in this.config) this.facet_row = this.config.facet_row;
      if ("facet_col" in this.config) this.facet_col = this.config.facet_col;
      if ("size_by_group" in this.config)
        this.size_by_group = this.config.size_by_group;
      if ("marker_size" in this.config && this.config.marker_size) {
        // we only want to change marker_size if there's a value other than null so
        // it defaults to the default size (3)
        this.marker_size = this.config.marker_size;
      }
      if ("jitter" in this.config && this.config.jitter) {
        // see 'marker_size' comment
        this.jitter = this.config.jitter;
      } else if (this.plot_type === 'violin') {
	      this.jitter = true;
      }
      if ("point_label" in this.config)
        this.point_label = this.config.point_label;
      if ("x_min" in this.config) this.x_min = this.config.x_min;
      if ("x_max" in this.config) this.x_max = this.config.x_max;
      if ("y_min" in this.config) this.y_min = this.config.y_min;
      if ("y_max" in this.config) this.y_max = this.config.y_max;
      if ("x_title" in this.config) this.x_title = this.config.x_title;
      if ("y_title" in this.config) this.y_title = this.config.y_title;
      if ("vlines" in this.config) this.vlines = this.config.vlines;
    },
    watch: {
      // Ensure a group is not used for two parameters
      x_axis(val) {
        if (this.y_axis === val) this.y_axis = null;
        if (this.z_axis === val) this.z_axis = null;
        if (this.size_by_group === val) this.size_by_group = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      y_axis(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.z_axis === val) this.z_axis = null;
        if (this.size_by_group === val) this.size_by_group = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      z_axis(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null; // TODO: 'Y' and 'Z' start as "raw_value"... need to prevent that
        if (this.size_by_group === val) this.size_by_group = null;
        // 'z' must be continuous and facets must be discrete so they will not overlap
      },
      color_name(val) {
        // Colors need to be cleared since the category is different.  New colors will be set after plot creation
        this.set_colors(null);
        this.set_color_palette(null);
        this.set_reverse_palette(false);
      },
      size_by_group(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.z_axis === val) this.z_axis = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      facet_row(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.size_by_group === val) this.size_by_group = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      facet_col(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.size_by_group === val) this.size_by_group = null;
        if (this.facet_row === val) this.facet_row = null;
      },
    },
    methods: {
      ...Vuex.mapActions([
        "fetch_plotly_data",
        "set_colors",
        "set_color_palette",
        "set_reverse_palette",
      ]),
      preview() {
        const {
          gene_symbol,
          analysis,
          colors,
          order,
          color_palette,
          reverse_palette,
        } = this.config;

        const config = {
          gene_symbol,
          analysis,
          colors,
          order,
          color_palette,
          reverse_palette,
          x_axis: this.x_axis,
          y_axis: this.y_axis,
          z_axis: this.z_axis,
          point_label: this.point_label,
          color_name: this.color_name,
          facet_row: this.facet_row,
          facet_col: this.facet_col,
          size_by_group: this.size_by_group,
          marker_size: this.marker_size,
          jitter: this.jitter,
          hide_x_labels: this.hide_x_labels,
          hide_y_labels: this.hide_y_labels,
          hide_legend: this.hide_legend,
          x_min: this.x_min,
          x_max: this.x_max,
          y_min: this.y_min,
          y_max: this.y_max,
          x_title: this.x_title,
          y_title: this.y_title,
          vlines: this.vlines,
        };

        const plot_type = this.plot_type;
        const dataset_id = this.dataset_id;
        this.fetch_plotly_data({ config, plot_type, dataset_id });
      },
      addRow() {
        // Add new 'vertical-line' component
        this.vlines.push({
          id: this.vline_counter++,
          vl_pos: null,
          vl_style: "solid",
        });
      },
      removeLast() {
        // Remove last 'vertical-line component
        this.vlines.pop();
      },
    },
    components: {
      verticalLine,
    },
  });

  const displayNameInput = Vue.component("display-name-input", {
    template: `
    <div id="display-name-input">
      <b-row class='mt-3'>
        <b-col>
        <h3>Display Name</h3>
        <b-form-group>
          <b-form-input
            type="text"
            :value='label'
            @input="update_display_name($event)"
          ></b-form-input>
        </b-form-group>
        </b-col>
      </b-row>
      <hr>
    </div>
    `,
    computed: {
      ...Vuex.mapState(["label"]),
    },
    methods: {
      ...Vuex.mapActions(["set_label"]),
      update_display_name(label) {
        this.set_label(label);
      },
    },
  });

  const saveDisplayBtn = Vue.component("save-display-btn", {
    template: `
    <b-form-group>
      <b-button @click='save' class='btn-purple float-right'>Save Display</b-button/>
    </b-form-group>
    `,
    props: ["display_id"],
    data() {
      return {
        variant: "",
      };
    },
    computed: {
      ...Vuex.mapState(["user", "dataset_id", "plot_type", "config", "label"]),
    },
    methods: {
      ...Vuex.mapActions(["update_display"]),
      async save() {
        const payload = {
          id: this.display_id,
          dataset_id: this.dataset_id,
          user_id: this.user.id,
          label: this.label,
          plot_type: this.plot_type,
          plotly_config: JSON.stringify({
            // depending on display type, this object will
            // have different properties
            ...this.config,
          }),
        };

        const res = await $.ajax({
          url: "./cgi/save_dataset_display.cgi",
          type: "POST",
          data: payload,
          dataType: "json",
        });

        if (res?.success) {
          if (this.display_id) {
            this.update_display(payload);
          }
          this.$router.push(`/dataset/${this.dataset_id}/displays`);
        } else {
          // not used
          this.variant = "danger";
        }
      },
    },
  });

  const geneSymbolInput = Vue.component("gene-symbol-input", {
    template: `
      <div id="gene-symbol-input">
        <b-row>
          <b-col>
            <b-icon class="float-right" v-b-tooltip.hover="'Enter gene to curate on.  Gene does not exist if autosearch box does not find it.'" icon="question-circle-fill"></b-icon>
            <h3>Gene Symbol</h3>
            <b-form-group>
              <b-form-input v-show='loading' disabled placeholder="Loading gene symbols..."></b-form-input>
              <vue-bootstrap-typeahead
                v-show="!loading"
                ref='gene_type_ahead'
                placeholder='Gene symbol'
                :value="config.gene_symbol"
                @hit='update_gene_symbol($event)'
                :data="gene_symbols"
                size='sm'
              />
              <slot name="svg_options"></slot>
              <slot name="button"></slot>
            </b-form-group>
          </b-col>
        </b-row>
        <hr>
      </div>
    `,
    props: ["analysis"],
    components: {
      VueBootstrapTypeahead,
    },
    data() {
      return {
        loading: false,
      };
    },
    watch: {
      analysis(analysis_id) {
        // when the analysis changes creating a tsne,
        // we want to fetch gene symbols for this h5ad
        const dataset_id = this.dataset_id;
        this.fetch_gene_symbols({ dataset_id, analysis_id });
      },
    },
    computed: {
      ...Vuex.mapState(["dataset_id", "config", "gene_symbols"]),
      is_gene_available() {
        return this.gene_symbols
          .map((gene) => gene.toLowerCase())
          .includes(this.config.gene_symbol.toLowerCase());
      },
    },
    async created() {
      const dataset_id = this.dataset_id;
      const analysis_id = this.analysis;
      this.loading = true;
      await this.fetch_gene_symbols({ dataset_id, analysis_id });
      this.loading = false;
    },
    async mounted() {
      // small hack to get around typeahead not allowing
      // a default value
      // https://github.com/alexurquhart/vue-bootstrap-typeahead/issues/22
      if (this.config.gene_symbol) {
        this.$refs.gene_type_ahead.inputValue = this.config.gene_symbol;
        // if there's a gene symbol we know that this is a saved
        // analysis and this gene exists
        this.$emit("gene-updated", true);
      }
    },
    methods: {
      ...Vuex.mapActions(["set_gene_symbol", "fetch_gene_symbols"]),
      update_gene_symbol(gene_symbol) {
        this.set_gene_symbol(gene_symbol);
        this.$emit("gene-updated", this.is_gene_available);
      },
    },
  });

  const displayOrder = Vue.component("display-order", {
    template: `
      <div id="display-order">
        <b-row>
          <b-col>
            <h4>Order</h4>
            <b-row v-for="dataseries in order" :key="dataseries.key">
            <b-col>
              <b-list-group class="mb-2">
              <h4>{{ dataseries.key }}</h4>
              <!-- plotly plots -->
              <draggable
                v-if="
                  !['tsne_static', 'umap_static', 'pca_static', 'svg'].includes(
                    plot_type
                  )
                "
                v-model="dataseries.value"
                @end="reorder_plotly_display"
              >
                <transition-group>
                  <b-list-group-item
                    v-for="elem in dataseries.value"
                    :key="elem"
                  >
                    {{ elem }}
                  </b-list-group-item>
                </transition-group>
              </draggable>
              <!-- scanpy plots -->
              <draggable
                v-else-if="
                  ['tsne_static', 'umap_static', 'pca_static'].includes(
                    plot_type
                  )
                "
                v-model="dataseries.value"
                @end="reorder_tsne_display"
              >
                <transition-group>
                  <b-list-group-item
                    v-for="elem in dataseries.value"
                    :key="elem"
                  >
                    {{ elem }}
                  </b-list-group-item>
                </transition-group>
              </draggable>
              </b-collapse>
              </b-list-group>
            </b-col>
            </b-row>
          </b-col>
        </b-row>
        <hr>
      </div>
    `,
    components: {
      draggable: vuedraggable,
    },
    data() {
      return {
        order: [],
      };
    },
    computed: {
      ...Vuex.mapState(["config", "plot_type", "dataset_id", "user"]),
    },
    created() {
      // Needed for initial display after first plotting preview
      this.get_order();

      this.unsubscribe = this.$store.subscribe((mutation, state) => {
        if (mutation.type === "set_order") {
          this.get_order();
        }
      });
    },
    beforeDestroy() {
      // If not present, subscriber will not stop even after component is destroyed
      this.unsubscribe();
    },
    methods: {
      ...Vuex.mapActions(["set_order", "fetch_plotly_data", "fetch_tsne_image"]),
      get_order() {
        const keys = Object.keys(this.config.order);
        const order = keys.map((key) => {
          return {
            key,
            value: [...this.config.order[key]],
          };
        });
        this.order = order;
      },
      reorder_plotly_display() {
        // Convert order from array of objects to a single object
        const order = this.order.reduce(
          (obj, item) => ((obj[item.key] = item.value), obj),
          {}
        );
        this.set_order(order);

        const config = this.config;
        const plot_type = this.plot_type;
        const dataset_id = this.dataset_id;

        this.fetch_plotly_data({ config, plot_type, dataset_id });
      },
      reorder_tsne_display() {
        // Convert order from array of objects to a single object
        const order = this.order.reduce(
          (obj, item) => ((obj[item.key] = item.value), obj),
          {}
        );
        this.set_order(order);
        const config = this.config;
        const plot_type = this.plot_type;
        const dataset_id = this.dataset_id;
        const analysis = config.analysis
        const analysis_owner_id = this.user.id;
        this.fetch_tsne_image({ config, plot_type, dataset_id, analysis, analysis_owner_id });
      },
    },
  });

  const displayColors = Vue.component("display-colors", {
    template: `
    <div id="display-colors">
      <b-row>
        <b-col>
          <h3>Color</h3>
          <b-row v-for="{name, color} in colors_array" :key='name'>
              <b-col>
                <label :for='name'>{{ name }}</label>
                <b-form-input
                  :key='name'
                  type="color"
                  :name="name"
                  :value='color'
                  @change='update_color(name, $event)'
                ></b-form-input>
                </b-col>
            </b-row>
        </b-col>
      </b-row>
      <hr>
    </div>
    `,
    data() {
      return {
        colors_array: [],
      };
    },
    computed: {
      ...Vuex.mapState(["config"]),
    },
    created() {
      // Needed for initial display after first plotting preview
      this.get_colors_array();

      this.unsubscribe = this.$store.subscribe((mutation, state) => {
        if (mutation.type === "set_colors") {
          this.get_colors_array();
        }
      });
    },
    beforeDestroy() {
      // If not present, subscriber will not stop even after component is destroyed
      this.unsubscribe();
    },
    methods: {
      ...Vuex.mapActions(["set_color"]),
      get_colors_array() {
        this.colors_array = Object.entries(this.config.colors).map(
          ([key, val]) => {
            return {
              name: key,
              color: val,
            };
          }
        );
      },
      update_color(name, color) {
        this.set_color({ name, color });
      },
    },
  });

  const displayPalettes = Vue.component("display-palettes", {
    template: `
    <div id="display-palettes">
      <b-row>
        <b-col>
          <h3>Color</h3>
          <b-row id="color-palette">
            <b-col>
              <b-form-group label-align-sm="right">
                <b-form-select size="sm" :options="options" v-model="palette">
                  <template slot="first">
                    <b-form-select-option :value="null" disabled>Please select an option</b-form-select-option>
                  </template>
                </b-form-select>
              </b-form-group>

              <b-form-group label-align-sm="right">
                <b-form-checkbox v-model='reverse_palette'>
                  Reverse colorscale
                </b-form-checkbox>
            </b-form-group>
            </b-col>
          </b-row>
        </b-col>
      </b-row>
      <hr>
    </div>
    `,
    data() {
      return {
        options: [
          {
            label: "Multi-color scales",
            options: [
              { value: "YlOrRd", text: "Yellow-Orange-Red" },
              { value: "Viridis", text: "Viridis" },
            ],
          },
          {
            label: "Single-color scales",
            options: [
              { value: "Greys", text: "Greyscale" },
              { value: "Blues", text: "Bluescale" },
              { value: "Purp", text: "Purplescale" },
            ],
          },
          {
            label: "Diverging Colorscales",
            options: [
              { value: "RdBu", text: "Red-Blue" },
              { value: "PiYG", text: "Pink-Yellow-Green" },
            ],
          },
        ],
        palette: null,
        reverse_palette: false,
      };
    },
    computed: {
      ...Vuex.mapState(["config", "plot_type", "dataset_id"]),
    },
    created() {
      // Needed for initial display after first plotting preview
      if ("color_palette" in this.config)
        this.palette = this.config.color_palette;
      if ("reverse_palette" in this.config)
        this.reverse_palette = this.config.reverse_palette;
    },
    watch: {
      palette(newval) {
        this.set_color_palette(newval);
        this.update_display();
      },
      reverse_palette(newval) {
        this.set_reverse_palette(newval);
        this.update_display();
      },
    },
    methods: {
      ...Vuex.mapActions([
        "set_color_palette",
        "set_reverse_palette",
        "fetch_plotly_data",
      ]),
      update_display() {
        const config = this.config;
        const plot_type = this.plot_type;
        const dataset_id = this.dataset_id;

        this.fetch_plotly_data({ config, plot_type, dataset_id });
      },
    },
  });

  const barDisplay = Vue.component("bar-display", {
    template: `
      <div id="bar-display">
        <gene-symbol-input
          v-model='config.gene_symbol'
          :analysis='config.analysis ? config.analysis.id : null'
          @gene-updated='is_gene_available = $event'
        ></gene-symbol-input>
        <transition name="fade" mode="out-in">
          <plotly-arguments v-if="config.gene_symbol !== '' && is_gene_available" />
        </transition>
        <slot name='chart'></slot>
        <transition name="fade" mode="out-in">
          <display-order
            v-if="'order' in config && Object.entries(config.order).length !== 0 && is_gene_available"
          ></display-order>
        </transition>
        <transition name="fade" mode="out-in">
          <display-colors
            v-if="is_there_data_to_save && is_gene_available && Object.entries(this.config.colors).length !== 0"
          ></display-colors>
        </transition>
        <transition name="fade" mode="out-in">
          <display-palettes
            v-if="is_there_data_to_save && is_gene_available && this.config.color_name !== (null || undefined) && Object.entries(this.config.colors).length === 0"
          ></display-palettes>
        </transition>
        <transition name="fade" mode="out-in">
          <display-name-input
            v-if='is_there_data_to_save && is_gene_available'
          ></display-name-input>
        </transition>
        <transition name="fade" mode="out-in">
          <save-display-btn
            v-if="is_there_data_to_save && is_gene_available"
            :display_id='display_id'
          ></save-display-btn>
        </transition>
      </div>
    `,
    props: ["display_id"],
    data() {
      return {
        is_gene_available: true,
      };
    },
    computed: {
      ...Vuex.mapState(["config"]),
      is_there_data_to_save() {
        return (
          "x_axis" in this.config &&
          "gene_symbol" in this.config &&
          this.config.gene_symbol !== ""
        );
      },
    },
    components: {
      PlotlyArguments,
      geneSymbolInput,
      displayNameInput,
      displayOrder,
      displayColors,
      displayPalettes,
      saveDisplayBtn,
    },
  });

  const lineDisplay = Vue.component("line-display", {
    extends: barDisplay,
  });

  const violinDisplay = Vue.component("violin-display", {
    extends: barDisplay,
  });

  const scatterDisplay = Vue.component("scatter-display", {
    extends: barDisplay,
  });

  const contourDisplay = Vue.component("contour-display", {
    extends: barDisplay,
  });

  const tsnePlotlyDisplay = Vue.component("tsne-plotly-display", {
    extends: scatterDisplay,
  });

  const plotlyDisplay = Vue.component("plotly-display", {
    template: `
      <div v-if="plot_type === 'bar'" id="plotly-display">
        <div>
          <transition name="fade" mode='out-in'>
            <bar-display
              v-if='!loading'
              :display_id='display_id'
            >
            </bar-display>
          </transition>
        </div>
      </div>
      <div v-else-if="plot_type === 'line'" id="plotly-display">
        <div>
          <transition name="fade" mode='out-in'>
            <line-display
              v-if='!loading'
              :display_id='display_id'
            >
            </line-display>
          </transition>
        </div>
      </div>
      <div v-else-if="plot_type === 'violin'" id="plotly-display">
        <div>
          <transition name="fade" mode='out-in'>
            <violin-display
              v-if='!loading'
              :display_id='display_id'
            >
            </violin-display>
          </transition>
        </div>
      </div>
      <div v-else-if='plot_type === "scatter"' id="plotly-display">
        <div>
          <transition name="fade" mode='out-in'>
            <scatter-display
              v-if='!loading'
              :display_id='display_id'
            >
            </scatter-display>
          </transition>
        </div>
      </div>
      <div v-else-if='plot_type === "contour"' id="plotly-display">
        <div>
          <transition name="fade" mode='out-in'>
            <contour-display
              v-if='!loading'
              :display_id='display_id'
            >
            </contour-display>
          </transition>
        </div>
      </div>
      <div v-else-if="['tsne_dynamic', 'tsne/umap_dynamic'].includes(plot_type)" id="plotly-display">
        <div>
          <transition name="fade" mode='out-in'>
            <tsne-plotly-display
              v-if='!loading'
              :display_id='display_id'
            >
            </tsne-plotly-display>
          </transition>
        </div>
      </div>
    </div>
    `,
    props: {
      display_id: {
        type: [String, null],
        default: null,
      },
    },
    components: {
      plotlyChart,
      barDisplay,
      lineDisplay,
      violinDisplay,
      scatterDisplay,
      contourDisplay,
      tsnePlotlyDisplay,
    },
    data() {
      return {
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState(["dataset_id", "plot_type", "config", "chart_data"]),
      is_creating_new_display() {
        return this.display_id === null;
      },
      is_there_data_to_draw() {
        return (
          Object.entries(this.chart_data).length !== 0 &&
          this.chart_data.constructor === Object
        );
      },
    },
    created() {
      this.fetch_h5ad_info({
        dataset_id: this.dataset_id,
        analysis: this.config.analysis,
      });
      if (!this.is_creating_new_display) {
        // if we are creating a new display, we do not
        // want to automatically generate a chart, and
        // wait for user to specify config options
        const config = this.config;
        const plot_type = this.plot_type;

        const dataset_id = this.dataset_id;
        this.fetch_plotly_data({ config, plot_type, dataset_id });
      }
    },
    methods: {
      ...Vuex.mapActions(["fetch_h5ad_info", "fetch_plotly_data"]),
      update_color({ name, color }) {
        const { data } = this.chart_data.plot_json;
        data
          .filter((el) => el.name === name)
          .forEach(({ marker }) => {
            marker.color = color;
          });

        const { colors } = this.config;
        colors[name] = color;

        // because vue wont detect these changes
        // we explitly reassign chart data with
        // new object
        this.chart_data = { ...this.chart_data };
      },
    },
  });

  const svgDisplay = Vue.component("svg-display", {
    template: `
      <div id="svg-display">
        <gene-symbol-input
        v-model='config.gene_symbol'
        @gene-updated='is_gene_available = $event'
        :analysis='config.analysis ? config.analysis.id : null'
        >
          <div slot="svg_options">
            <hr>
            <b-col class='mt-2'>
              <b-row>
                <label for='low'>Low Color</label>
                  <b-form-input
                    type="color"
                    name="low"
                    v-model='low_color'
                  >
                  </b-form-input>
              </b-row>
              <b-row>
                <label for="mid">Mid Color</label>
                <b-form-input
                  type="color"
                  name="mid"
                  v-model='mid_color'
                >
                </b-form-input>
              </b-row>
              <b-row>
                <label for="high">High Color</label>
                <b-form-input
                  type="color"
                  name="high"
                  v-model='high_color'
                >
                </b-form-input>
              </b-row>
            </b-col>
            </b-row>
            <b-row class='mt-3'>
              <b-col>
                <b-button
                  :disabled="config.gene_symbol === '' || !is_gene_available"
                  class='btn-purple float-right'
                  @click="preview"
                  >Preview SVG</b-button>
              </b-col>
            </b-row>
          </div>
        </gene-symbol-input>
        <transition name="fade" mode="out-in">
          <display-name-input
            v-if='is_there_data_to_draw && is_gene_available'
          ></display-name-input>
        </transition>
        <transition name="fade" mode="out-in">
          <save-display-btn
            v-if="is_there_data_to_draw && is_gene_available"
            :display_id='display_id'
          ></save-display-btn>
      </transition>
      </div>
    `,
    props: {
      display_id: {
        default: null,
      },
    },
    components: {
      geneSymbolInput,
      svgChart,
      displayNameInput,
      saveDisplayBtn,
    },
    data() {
      return {
        loading: false,
        is_gene_available: true,
      };
    },
    computed: {
      ...Vuex.mapState(["dataset_id", "config", "chart_data"]),
      low_color: {
        get() {
          return this.config.colors.low_color;
        },
        set(color) {
          this.set_color({ name: "low_color", color });
        },
      },
      mid_color: {
        get() {
          return this.config.colors.mid_color;
        },
        set(color) {
          this.set_color({ name: "mid_color", color });
        },
      },
      high_color: {
        get() {
          return this.config.colors.high_color;
        },
        set(color) {
          this.set_color({ name: "high_color", color });
        },
      },
      is_creating_new_display() {
        return this.display_id === null;
      },
      is_there_data_to_draw() {
        return (
          Object.entries(this.chart_data).length !== 0 &&
          this.chart_data.constructor === Object
        );
      },
      is_there_data_to_save() {
        return "colors" in this.config && "gene_symbol" in this.config;
      },
    },
    async created() {
      if (!this.is_creating_new_display) {
        // if we are creating a new display, we do not
        // want to automatically generate a chart, and
        // wait for user to specify config options
        const { gene_symbol } = this.config;
        const dataset_id = this.dataset_id;

        this.fetch_svg_data({ gene_symbol, dataset_id });
      }
    },
    methods: {
      ...Vuex.mapActions(["fetch_svg_data", "set_color"]),
      preview() {
        const { gene_symbol } = this.config;
        const dataset_id = this.dataset_id;

        this.fetch_svg_data({ gene_symbol, dataset_id });
      },
    },
  });

  const chooseAnalysis = Vue.component("choose-analysis", {
    template: `
      <b-card-body no-body
        class='elevation border-0 mt-5'>
      <b-row class='m-4'>
        <b-col>
          <h3>Choose Analyses</h3>
        </b-col>
        <b-col>
          <b-form id="choose-analysis">
            <b-form-radio-group
              v-model='private_or_public'
              buttons
              button-variant='light'
              :options="['Public', 'Private']"
            />
            <b-form-select
              :value='selected_analysis'
              @input='analysis_selected($event)'
              class="mt-2 mb-2 mr-sm-2 mb-sm-0"
              :disabled="private_or_public === null"
              :options='analyses'
            >
            </b-form-select>
          </b-form>
        </b-col>
    </b-row>
    </b-card-body>
    `,
    data() {
      return {
        private_or_public: null,
        selected_analysis: null,
        private: [],
        public: [],
      };
    },
    computed: {
      ...Vuex.mapState(["dataset_id", "config"]),
      analysis() {
        return this.private_or_public === "Public" ? this.public_labels.find(
          (el) => el.value === this.config.analysis_id
        ) : this.private_labels.find(
          (el) => el.value === this.config.analysis_id
        );
      },
      ana_private_or_public() {
        // If an analaysis id is passed,
        // check if its public or private
        return this.public
          .map((ana) => ana.id)
          .includes(this.config.analysis_id)
          ? "Public"
          : "Private";
      },
      private_labels() {
        return this.private.map((ana) => {
          return {
            value: ana.id,
            text: ana.label,
          };
        });
      },
      public_labels() {
        return this.public.map((ana) => {
          return {
            value: ana.id,
            text: ana.label,
          };
        });
      },
      analyses() {
        return this.private_or_public === "Public" ? this.public_labels : this.private_labels;
      },
    },
    async created() {
      const { data } = await axios.get(
        `./api/h5ad/${this.dataset_id}/analyses`
      );
      const public_analysis = data.public;
      const private_analysis = data.private;
      this.public = public_analysis;
      this.private = private_analysis;

      if (this.config.analysis_id) {
        this.private_or_public = this.ana_private_or_public;
        this.selected_analysis = this.analysis.value;
      }
    },
    methods: {
      ...Vuex.mapActions(["set_analysis_id"]),
      analysis_selected(analysis) {
        this.set_analysis_id(analysis);
      },
    },
  });

  const tsneArguments = Vue.component("tsne-arguments", {
    template: `
      <div id="tsne-arguments">
        <b-row>
          <b-col>
          <h3>Display Parameters</h3>
          <b-form-group>

              <b-form-group label-align-sm="right">
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for X axis'" icon="question-circle-fill"></b-icon>
                <label class="mb-0">X</label>
                <b-form-select :options='columns' v-model='x_axis' size='sm'>
                  <template slot="first">
                    <option :value=null></option>
                  </template>
                </b-form-select>
              </b-form-group>

              <b-form-group label-align-sm="right">
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for Y axis'" icon="question-circle-fill"></b-icon>
                <label class="mb-0">Y</label>
                <b-form-select :options='columns' v-model='y_axis' size='sm'>
                  <template slot="first">
                  <option value=null></option>
                  </template>
                </b-form-select>
              </b-form-group>

              <b-form-group label-align-sm="right">
                <b-form-checkbox v-model='show_colorized_legend'>
                  <label class="mb-0" id="show_colorized_legend">Show colorized legend</label>
                </b-form-checkbox>
                <b-tooltip target="show_colorized_legend" triggers="hover">
                  Check to enable ability to color by a category
                </b-tooltip>
              </b-form-group>

              <b-form-group v-show='show_colorized_legend' label-align-sm="right">
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for coloring plot'" icon="question-circle-fill"></b-icon>
                <label class="mb-0">Colorize legend by:</label>
                <b-form-select :options="Object.keys(levels)" v-model='colorize_legend_by' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>

              <b-form-group v-show='show_colorized_legend && colorize_legend_by' label-align-sm="right">
                <b-icon class="float-right" v-b-tooltip.hover="'Data series for creating separate plots by individual entities in a group'" icon="question-circle-fill"></b-icon>
                <label class="mb-0">Plot by group:</label>
                <b-form-select :options="Object.keys(levels)" v-model='plot_by_group' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
                <b-alert :show="plot_by_group && num_plot_by_group_levels >= 10" variant="warning">
                  Plot generation may be slow if the selected category has a large number of groups.
                </b-alert>
              </b-form-group>

              <b-form-group v-show='show_colorized_legend && plot_by_group' label-align-sm="right">
                <b-icon class="float-right" v-b-tooltip.hover="'Maximum number of plots per row.  If not provided, all plots will be on one row'" icon="question-circle-fill"></b-icon>
                <label class="mb-0">Max columns per row:</label>
                <b-form-input type="number" v-model='max_columns' number size='sm' min="1"></b-form-input>
              </b-form-group>

              <b-form-group v-show='show_colorized_legend && !plot_by_group' label-align-sm="right">
                <b-form-checkbox v-model='skip_gene_plot'>
                  <label class="mb-0" id="skip_gene_plot">Skip gene symbol plot</label>
                </b-form-checkbox>
                <b-tooltip target="skip_gene_plot" triggers="hover">
                  Check to skip the gene symbol plot
                </b-tooltip>
              </b-form-group>

              <b-form-group v-show='show_colorized_legend' label-align-sm="right">
                <b-form-checkbox v-model='horizontal_legend'>
                  <label class="mb-0" id="horizontal_legend">Place legend under plots</label>
                </b-form-checkbox>
                <b-tooltip target="horizontal_legend" triggers="hover">
                  Check to make horizontal legend along the bottom of the plotspace
                </b-tooltip>
              </b-form-group>

          </b-form-group>
          </b-col>
        </b-row>
        <hr>
      </div>
    `,
    components: {},
    data() {
      return {
        // Since the config may not have these values, create so they aren't undefined
        show_colorized_legend: false,
        horizontal_legend: false,
        skip_gene_plot: false,
        plot_by_group: null,
        max_columns: 4,
        colorize_legend_by: null,
      };
    },
    computed: {
      ...Vuex.mapState([
        "user",
        "dataset_id",
        "config",
        "columns",
        "plot_type",
        "image_data",
        "tsne_is_loading",
        "levels",
      ]),
      x_axis: {
        get() {
          return this.$store.state.config.x_axis;
        },
        set(value) {
          this.$store.commit("set_x_axis", value);
        },
      },
      y_axis: {
        get() {
          return this.$store.state.config.y_axis;
        },
        set(value) {
          this.$store.commit("set_y_axis", value);
        },
      },
      num_plot_by_group_levels() {
        if (this.plot_by_group) {
          return Object.keys(this.levels[this.plot_by_group]).length;
        }
        return -1;
      }
    },
    created() {
      if ("x_axis" in this.config) this.x_axis = this.config.x_axis;
      if ("y_axis" in this.config) this.y_axis = this.config.y_axis;
      if ("colorize_legend_by" in this.config)
        this.show_colorized_legend = true;
        this.colorize_legend_by = this.config.colorize_legend_by;
      if ("plot_by_group" in this.config)
        this.plot_by_group = this.config.plot_by_group;
      if ("max_columns" in this.config)
        this.max_columns = this.config.max_columns;
      if ("skip_gene_plot" in this.config)
        this.skip_gene_plot = this.config.skip_gene_plot;
      if ("horizontal_legend" in this.config)
        this.horizontal_legend = this.config.horizontal_legend;

      this.fetch_h5ad_info({
        dataset_id: this.dataset_id,
        analysis: this.config.analysis,
      });

      if (this.plot_params_ready()) {
        this.draw_image();
      }

    },
    watch: {
      show_colorized_legend(newval, oldval) {
        // if deselected, clear colorize legend select box
        if (newval !== true) {
          this.colorize_legend_by = null;
          this.skip_gene_plot = null;
          this.horizontal_legend = null;
          this.plot_by_group = null;
          this.max_columns = null;
        }
      },
      colorize_legend_by(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
            // Set order in config so "display-order" will render
            if (newval !== null && this.levels) {
              const colorize_key = this.colorize_legend_by;
              const order = {};
              order[colorize_key] = this.levels[colorize_key];

              if (this.plot_by_group !== null) {
                  // Add separately in case both are same dataseries group
                  const group_key = this.plot_by_group;
                  order[group_key] = this.levels[group_key];
              }

              // This is to prevent a bug where the levels have not been set yet when loading a display.
              if (oldval !== null) {
                this.$store.commit("set_order", order);
              }
            }
          this.draw_image();
        }
      },
      x_axis(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
          this.draw_image();
        }
      },
      y_axis(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
          this.draw_image();
        }
      },
      skip_gene_plot(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
          this.draw_image();
        }
      },
      horizontal_legend(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
          this.draw_image();
        }
      },
      plot_by_group(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
          // Plotting by group colors by gene symbol, so cannot skip gene plot
          if (newval !== null) this.skip_gene_plot = null;

          // Currently only works if colorize_legend_by is set
          if (this.colorize_legend_by !== null) {
            // Set order in config so "display-order" will render
            if (newval !== null && this.levels) {
              const group_key = this.plot_by_group;
              // Add separately in case both are same dataseries group
              const colorize_key = this.colorize_legend_by;
              const order = {};
              order[group_key] = this.levels[group_key];
              order[colorize_key] = this.levels[colorize_key];
              this.$store.commit("set_order", order);
            }
            this.draw_image();
          }
        }
      },
      max_columns(newval, oldval) {
        if (newval != oldval && this.plot_params_ready()) {
          this.draw_image();
        }
      },
    },
    methods: {
      ...Vuex.mapActions([
        "fetch_h5ad_info",
        "fetch_tsne_image",
        "set_image_data",
        "set_success",
        "set_message",
        "set_tsne_is_loading",
        "set_order",
      ]),
      plot_params_ready() {
        return this.x_axis && this.x_axis !== "null" && this.y_axis && this.y_axis !== "null";
      },
      draw_image() {
        const dataset_id = this.dataset_id;
        const analysis = this.config.analysis;
        const plot_type = this.plot_type;
        const analysis_owner_id = this.user.id;

        const config = {
          gene_symbol: this.config.gene_symbol,
          colorize_legend_by: this.colorize_legend_by,
          horizontal_legend: this.horizontal_legend,
          skip_gene_plot: this.skip_gene_plot,
          plot_by_group: this.plot_by_group,
          max_columns: this.max_columns,
          x_axis: this.x_axis,
          y_axis: this.y_axis,
          colors: this.colors,
          order: this.config.order,
          // helps stop caching issues
          timestamp: new Date().getTime(),
        };

        this.fetch_tsne_image({config, plot_type, dataset_id, analysis, analysis_owner_id})
      },
    },
  });

  const tsneDisplay = Vue.component("tsne-display", {
    template: `
      <div id="tsne-display">
        <gene-symbol-input
          v-model='config.gene_symbol'
          :analysis='config.analysis ? config.analysis.id : null'
          @gene-updated='is_gene_available = $event'
        ></gene-symbol-input>
        <transition name="fade" mode="out-in">
          <tsne-arguments v-if="config.gene_symbol !== '' && is_gene_available" />
        </transition>
        <transition name="fade" mode="out-in">
          <display-order
            v-if="
              'order' in config &&
              Object.entries(config.order).length !== 0 &&
              is_gene_available
            "
          ></display-order>
        </transition>
        <transition name="fade" mode="out-in">
        <display-name-input
          v-if='is_gene_available'
        ></display-name-input>
      </transition>
      <transition name="fade" mode="out-in">
        <save-display-btn
          v-if="is_gene_available"
          :display_id='display_id'
        ></save-display-btn>
      </transition>
      </div>
    `,
    props: ["display_id"],
    components: {
      chooseAnalysis,
      geneSymbolInput,
      tsneChart,
      displayNameInput,
      displayOrder,
      saveDisplayBtn,
    },
    data() {
      return {
        is_gene_available: true,
        show_tsne: false,
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState(["dataset_id", "config", "dataset_type", "analysis"]),
    },
  });

  const primaryConfig = Vue.component("primary-config", {
    template: `
      <div id="primary-config">
        <choose-display-type
        :analysis_id='config.analysis ? config.analysis.id : null'>
        </choose-display-type>
        <div v-if="is_type_plotly">
          <plotly-display></plotly-display>
        </div>
        <div v-else-if="is_type_svg">
          <svg-display></svg-display>
        </div>
        <div v-else-if="is_type_tsne">
          <tsne-display></tsne-display>
        </div>
      </div>
    `,
    components: {
      chooseDisplayType,
      plotlyDisplay,
      svgDisplay,
      tsneDisplay
    },
    computed: {
      ...Vuex.mapState(["config", "plot_type", "dataset_type"]),
      is_type_plotly() {
        return (
          this.plot_type === "bar" ||
          this.plot_type === "scatter" ||
          this.plot_type === "line" ||
          this.plot_type === "violin" ||
          this.plot_type === "contour" ||
          this.plot_type === "tsne/umap_dynamic"
        );
      },
      is_type_svg() {
        return this.plot_type === "svg";
      },
      is_type_tsne() {
        return (
          this.plot_type === "tsne_static" ||
          this.plot_type === "umap_static" ||
          this.plot_type === "pca_static" ||
          this.plot_type === "tsne"
        );
      },
    },
  });

  const chooseStoredAnalysis = Vue.component("choose-stored-analysis", {
    template: `
       <div id="choose-stored-analysis">
        <b-row>
          <b-col>
            <h3>Stored Analysis</h3>
            <b-form-select v-model="selected_analysis" class="mb-3" :disabled='loading'>
              <option v-if='loading' slot="first" :value="null">Loading stored analyses...</option>
              <option v-else :value="null">Please select stored analysis</option>
              <hr>
              <optgroup label="Public saved analysis">
                <option v-for="analysis in public_analyses" :value="{ id: analysis.id, type: analysis.type }"> {{ analysis.label }}</option>
                <option v-if="!public_analyses.length" disabled>No public analyses for this dataset</option>
                </optgroup>
                <hr>
                <optgroup label="Your saved analysis">
                <option v-for="analysis in private_analyses" :value="{ id: analysis.id, type: analysis.type }"> {{ analysis.label }} </option>
                <option v-if="!private_analyses.length" disabled>No saved analyses for this dataset</option>
              </optgroup>
            </b-form-select>
          </b-col>
        </b-row>
  </div>
    `,
    data() {
      return {
        selected_analysis: null,
        loading: true,
        public_analyses: [],
        private_analyses: [],
      };
    },
    async created() {
      this.loading = true;
      const { data } = await axios.get(
        `/./api/h5ad/${this.dataset_id}/analyses`
      );
      const { public: public_analyses, private: private_analyses } = data;

      this.public_analyses = public_analyses;
      this.private_analyses = private_analyses;
      this.loading = false;
    },
    watch: {
      selected_analysis(new_analysis, old_analysis) {
        this.set_analysis(new_analysis);
      },
    },
    methods: {
      ...Vuex.mapActions(["set_analysis"]),
    },
    computed: {
      ...Vuex.mapState(["dataset_id", "config"]),
    },
  });

  const storedAnalysisConfig = Vue.component("stored-analysis-config", {
    template: `
      <div id="stored-analysis-config">
        <choose-stored-analysis></choose-stored-analysis>
        <hr v-if="analysis">
        <primary-config v-if="analysis"></primary-config>
      </div>
    `,
    components: { chooseStoredAnalysis, primaryConfig },
    computed: {
      ...Vuex.mapState(["analysis"]),
    },
  });

  const configurationPanel = Vue.component("configuration-panel", {
    template: `
      <div id="configuration-panel">
        <b-row>
          <b-col>
          <h3>Dataset Type</h3>
            <b-form-group>
              <b-form-radio v-model="selected" value="primary">Primary Data</b-form-radio>
              <b-form-radio v-model="selected" value="analysis">Stored Analysis</b-form-radio>
            </b-form-group>
            </b-col>
        </b-row>
        <hr>
        <b-row>
          <b-col>
            <div v-if="selected == 'primary'">
              <primary-config></primary-config>
            </div>
            <div v-else>
              <stored-analysis-config v-if='dataset_id'></stored-analysis-config>
            </div>
          </b-col>
        </b-row>
      </div>
      `,
    components: {
      primaryConfig,
      storedAnalysisConfig,
    },
    data() {
      return {
        selected: "primary",
      };
    },
    computed: {
      ...Vuex.mapState(['dataset_type', 'dataset_id']),
    },
    watch: {
      selected(newValue) {
        this.set_dataset_type(newValue);
      },
    },
    created() {
      this.selected = this.dataset_type;
    },
    methods: {
      ...Vuex.mapActions(["set_dataset_type"]),
    },
  });

  const newDisplay = Vue.component("new-display", {
    template: `
      <b-container fluid id="new-display">
        <b-row>
          <b-col cols='2'>
            <configuration-panel></configuration-panel>
          </b-col>
          <b-col cols='10'>
            <plotly-chart v-if='is_type_plotly && is_there_data_to_draw' :data='chart_data' class="sticky-chart"></plotly-chart>
            <svg-chart
              v-if='is_type_svg && is_there_data_to_draw'
              :chart_data='chart_data'
              class="sticky-chart"
              :low_color="config.colors.low_color"
              :mid_color="config.colors.mid_color"
              :high_color="config.colors.high_color">
            </svg-chart>
            <tsne-chart
              v-if='is_type_tsne && gene_selected'
              display
              :gene_symbol='config.gene_symbol'
            ></tsne-chart>
          </b-col>
        </b-row>
      </b-container>
    `,
    components: {
      plotlyDisplay,
      plotlyChart,
      svgChart,
      tsneChart,
      svgDisplay,
      tsneDisplay,
      configurationPanel,
    },
    computed: {
      ...Vuex.mapState(["config", "plot_type", "chart_data", "analysis"]),
      is_type_plotly() {
        return (
          this.plot_type === "bar" ||
          this.plot_type === "scatter" ||
          this.plot_type === "line" ||
          this.plot_type === "violin" ||
          this.plot_type === "contour" ||
          this.plot_type === "tsne/umap_dynamic"
        );
      },
      is_type_svg() {
        return this.plot_type === "svg";
      },
      is_type_tsne() {
        // TODO: move to methods()
        return (
          this.plot_type === "tsne_static" ||
          this.plot_type === "umap_static" ||
          this.plot_type === "pca_static" ||
          this.plot_type === "tsne"
        );
      },
      is_there_data_to_draw() {
        return (
          Object.entries(this.chart_data).length !== 0 &&
          this.chart_data.constructor === Object
        );
      },
      gene_selected() {
        return this.config.gene_symbol;
      },
    },
  });

  const datasetDisplay = Vue.component("dataset-display", {
    template: `
      <b-container v-if='!loading' fluid id="dataset-display">
        <b-row>
          <b-col cols='2'>
            <plotly-display v-if="is_type_plotly" :display_id='display_id'></plotly-display>
            <svg-display v-else-if='is_type_svg' :display_id="display_id"></svg-display>
            <tsne-display v-else-if='is_type_tsne' :display_id='display_id'></tsne-display>
          </b-col>
          <b-col cols='10'>
            <plotly-chart
              v-if='is_type_plotly && is_there_data_to_draw'
              :data='chart_data'
              class="sticky-chart"
            ></plotly-chart>
            <svg-chart
              v-if='is_type_svg && is_there_data_to_draw'
              :chart_data='chart_data'
              class="sticky-chart"
              :low_color="config.colors.low_color"
              :mid_color="config.colors.mid_color"
              :high_color="config.colors.high_color">
            </svg-chart>
            <tsne-chart
              v-if='is_type_tsne'
            ></tsne-chart>
          </b-col>
        </b-row>
      </b-container>
    `,
    props: ["display_id"],
    components: {
      newDisplay,
      plotlyDisplay,
      tsneDisplay,
    },
    data() {
      return {
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState(["plot_type", "chart_data", "config", "user", "dataset_id"]),
      ...Vuex.mapGetters(["user_display"]),
      is_creating_new_display() {
        return this.display_id === "new";
      },
      is_type_plotly() {
        // handle legacy tsne dynamic plot option
        let plot_type = this.plot_type;
        if (plot_type === "tsne_dynamic") {
          plot_type = "tsne/umap_dynamic";
        }
        return (
          plot_type === "bar" ||
          plot_type === "scatter" ||
          plot_type === "line" ||
          plot_type === "violin" ||
          plot_type === "contour" ||
          plot_type === "tsne/umap_dynamic"
        );
      },
      is_type_svg() {
        return this.plot_type === "svg";
      },
      is_type_tsne() {
        return (
          this.plot_type === "tsne_static" ||
          this.plot_type === "umap_static" ||
          this.plot_type === "pca_static" ||
          this.plot_type === "tsne"
        );
      },
      is_there_data_to_draw() {
        return (
          Object.entries(this.chart_data).length !== 0 &&
          this.chart_data.constructor === Object
        );
      },
    },
    async created() {
      // If user displays not generated (such as refreshing "edit" route page, then generate first)
      if (! this.user_displays) {
        const user_id = this.user.id;
        const dataset_id = this.dataset_id;
        await this.fetch_user_displays({user_id, dataset_id});
      }
      const display_data = this.user_display(this.display_id);
      this.set_display_data(display_data);
    },
    methods: {
      ...Vuex.mapActions([
      "set_display_data",
      "fetch_user_displays",
    ]),
    },
  });

  const store = new Vuex.Store({
    state: {
      user: null,
      display_id: null,
      user_displays: [],
      owner_displays: [],
      default_display_id: 0,
      dataset_id: "",
      owner_id: null,
      config: {
        gene_symbol: "",
        analysis: null,
        colors: null, // Color mapping for dataseries groups
        color_palette: null, // Predefined swatches for continuous data
        reverse_palette: false,
        order: null,
        // Other properties will be set reactive (must add via Vue.set)
        // depending on the chart type
        // TODO: Branch off into its own module with states/getters/mutations/actions
      },
      gene_symbols: [],
      dataset_type: "primary",
      // why is analysis here too and within config?
      analysis: null,
      columns: [],
      levels: {},
      chart_data: {},
      is_public: false,
      label: "",
      title: "",
      plot_type: null,
      loading_chart: false,
      available_plot_types: {},
      plot_type_cancel_source: null,    // Cancel token for API call
      image_data: null,
      tsne_is_loading: false,
      success: 0,
      message: "",
    },

    getters: {
      is_user_owner(state) {
        return state.user.id === state.owner_id;
      },
      user_display(state) {
        return (display_id) =>
          state.user_displays.find((display) => display.id == display_id);
      },
      owner_display(state) {
        return (display_id) =>
          state.owner_displays.find((display) => display.id == display_id);
      },
      user_displays(state) {
        return state.user_displays.map((display) => {
          return {
            is_default: state.default_display_id == display.id,
            ...display,
          };
        });
      },
      owner_displays(state) {
        return state.owner_displays.map((display) => {
          return {
            is_default: state.default_display_id == display.id,
            ...display,
          };
        });
      },
    },

    mutations: {
      set_user(state, user) {
        state.user = user;
      },
      set_dataset_id(state, dataset_id) {
        state.dataset_id = dataset_id;
      },
      set_default_display_id(state, default_display_id) {
        state.default_display_id = default_display_id;
      },
      set_owner_id(state, owner_id) {
        state.owner_id = owner_id;
      },
      set_is_public(state, is_public) {
        state.is_public = is_public;
      },
      set_title(state, title) {
        state.title = title;
      },
      set_dataset_type(state, dataset_type) {
        state.dataset_type = dataset_type;
      },
      set_plot_type(state, plot_type) {
        // reset config, as different display types
        // has different configs
        state.config = {
          gene_symbol: "",
          analysis: state.config.analysis,
          colors: null,
        };

        if (
          plot_type === "bar" ||
          plot_type === "line" ||
          plot_type === "violin" ||
          plot_type === "scatter" ||
          plot_type === "contour" ||
          plot_type === "tsne/umap_dynamic"
        ) {
          Vue.set(state.config, "x_axis", null);
          Vue.set(state.config, "y_axis", null);
          Vue.set(state.config, "z_axis", null);
          Vue.set(state.config, "x_min", null);
          Vue.set(state.config, "y_min", null);
          Vue.set(state.config, "x_title", null);
          Vue.set(state.config, "y_title", null);
          Vue.set(state.config, "point_label", null);
          Vue.set(state.config, "hide_x_labels", false);
          Vue.set(state.config, "hide_y_labels", false);
          Vue.set(state.config, "hide_legend", false);
          Vue.set(state.config, "color_name", null);
          Vue.set(state.config, "facet_row", null);
          Vue.set(state.config, "facet_col", null);
          Vue.set(state.config, "size_by_group", null);
          Vue.set(state.config, "marker_size", null);
          Vue.set(state.config, "jitter", null);
          Vue.set(state.config, "vlines", []);
          Vue.set(state.config, "colors", {});
          Vue.set(state.config, "color_palette", null);
          Vue.set(state.config, "reverse_palette", false);
          Vue.set(state.config, "order", {});
        } else if (plot_type === "svg") {
          Vue.set(state.config, "colors", {
            // arbituary default colors (purple)
            low_color: "#e7d1d5",
            mid_color: null,
            high_color: "#401362",
          });
        } else if (
          plot_type === "tsne_static" ||
          plot_type === "umap_static" ||
          plot_type === "pca_static" ||
          plot_type === "tsne"
        ) {
          // tsne
          Vue.set(state.config, "x_axis", null);
          Vue.set(state.config, "y_axis", null);
          Vue.set(state.config, "colors", {});
          Vue.set(state.config, "order", {});
          Vue.set(state.config, "colorize_legend_by", null);
          Vue.set(state.config, "plot_by_group", null);
          Vue.set(state.config, "max_columns", null);
          Vue.set(state.config, "skip_gene_plot", false);
          Vue.set(state.config, "horizontal_legend", false);
        }
        state.plot_type = plot_type;
      },
      set_display_label(state, display_label) {
        state.display_label = display_label;
      },
      set_user_displays(state, user_displays) {
        state.user_displays = [...user_displays];
      },
      set_owner_displays(state, owner_displays) {
        state.owner_displays = [...owner_displays];
      },
      set_available_plot_types(state, available_plot_types) {
        state.available_plot_types = { ...available_plot_types };
      },
      set_plot_type_cancel_source (state, cancel_source) {
        state.cancel_source = cancel_source;
      },
      set_columns(state, columns) {
        state.columns = [...columns];
      },
      set_levels(state, levels) {
        state.levels = { ...levels };
      },
      reset_config(state) {
        state.config = { gene_symbol: "" };
      },
      set_config(state, config) {
        state.config = { ...config };
      },
      set_chart_data(state, data) {
        state.chart_data = { ...data };
      },
      set_x_axis(state, x_axis) {
        state.config.x_axis = x_axis;
      },
      set_y_axis(state, y_axis) {
        state.config.y_axis = y_axis;
      },
      set_z_axis(state, z_axis) {
        state.config.z_axis = z_axis;
      },
      set_x_min(state, x_min) {
        state.config.x_min = x_min;
      },
      set_y_min(state, y_min) {
        state.config.y_min = y_min;
      },
      set_x_max(state, x_max) {
        state.config.x_max = x_max;
      },
      set_y_max(state, y_max) {
        state.config.y_max = y_max;
      },
      set_x_title(state, x_title) {
        state.config.x_title = x_title;
      },
      set_y_title(state, y_title) {
        state.config.y_title = y_title;
      },
      set_vlines(state, vlines) {
        state.config.vlines = vlines;
      },
      set_point_label(state, point_label) {
        state.config.point_label = point_label;
      },
      set_hide_x_labels(state, hide) {
        state.config.hide_x_labels = hide;
      },
      set_hide_y_labels(state, hide) {
        state.config.hide_y_labels = hide;
      },
      set_hide_legend(state, hide) {
        state.config.hide_legend = hide;
      },
      set_color_name(state, color) {
        state.config.color_name = color;
      },
      set_facet_row(state, facet_row) {
        state.config.facet_row = facet_row;
      },
      set_facet_col(state, facet_col) {
        state.config.facet_col = facet_col;
      },
      set_size_by_group(state, size_by_group) {
        state.config.size_by_group = size_by_group;
      },
      set_marker_size(state, marker_size) {
        state.config.marker_size = marker_size;
      },
      set_jitter(state, jitter) {
        state.config.jitter = jitter;
      },
      set_order(state, order) {
        state.config.order = { ...order };
      },
      set_colors(state, colors) {
        if (typeof colors !== "string") state.config.colors = { ...colors };
      },
      set_color_palette(state, palette) {
        state.config.color_palette = palette;
      },
      set_reverse_palette(state, isReverse) {
        state.config.reverse_palette = isReverse;
      },
      set_color(state, { name, color }) {
        if (state.plot_type !== "svg") {
          // update plotly chart data
          const { data } = state.chart_data.plot_json;
          data
            .filter((el) => el.name === name)
            .forEach(({ marker }) => {
              marker.color = color;
            });
        }
        state.config.colors[name] = color;

        // because vue wont detect these changes
        // we explitly reassign chart data with
        // new object
        state.chart_data = { ...state.chart_data };
      },
      set_gene_symbol(state, gene_symbol) {
        state.config.gene_symbol = gene_symbol;
      },
      set_gene_symbols(state, gene_symbols) {
        state.gene_symbols = [...gene_symbols];
      },
      set_label(state, label) {
        state.label = label;
      },
      set_analysis_id(state, analysis) {
        state.config.analysis_id = analysis;
      },
      set_analysis(state, analysis) {
        // reset config, as different display types
        // has different configs
        state.config = {
          gene_symbol: "",
          colors: null,
        };
        state.plot_type = null;
        state.available_plot_types = {};

        state.analysis = analysis;
        state.config.analysis = analysis;
      },
      update_display(state, display) {
        const { id } = display;
        const i = state.user_displays.findIndex((el) => el.id == id);
        Vue.set(state.user_displays, i, display);
      },
      delete_display(state, { display_id }) {
        state.user_displays = [
          ...state.user_displays.filter((display) => display.id != display_id),
        ];
      },
      set_loading_chart(state, is_chart_loading) {
        state.loading_chart = is_chart_loading;
      },
      set_colorize_legend_by(state, legend_by) {
        state.config.colorize_legend_by = legend_by;
      },
      set_skip_gene_plot(state, skip_plot) {
        state.config.skip_gene_plot = skip_plot;
      },
      set_horizontal_legend(state, horizontal_legend) {
        state.config.horizontal_legend = horizontal_legend;
      },
      set_plot_by_group(state, group) {
        state.config.plot_by_group = group;
      },
      set_max_columns(state, max_cols) {
        state.config.max_columns = max_cols;
      },
      /*  Set earlier
      set_x_axis(state, x_axis) {
        state.config.x_axis = x_axis;
      },
      set_y_axis(state, y_axis) {
        state.config.y_axis = y_axis;
      },
      */
      set_image_data(state, image_data) {
        state.image_data = image_data;
      },
      set_success(state, success) {
        state.success = success;
      },
      set_message(state, message) {
        state.message = message;
      },
      set_tsne_is_loading(state, is_loading) {
        state.tsne_is_loading = is_loading;
      },
    },

    actions: {
      set_dataset_type({ commit }, dataset_type) {
        commit("set_dataset_type", dataset_type);
        // When display type is changed, we want to
        // reset these other settings so dependent
        // components
        commit("set_plot_type", null);
        commit("set_chart_data", {});
        commit("set_gene_symbol", "");
        commit("set_analysis", null);
      },
      update_display_label({ commit }, display_label) {
        commit("set_display_label", display_label);
      },
      set_dataset_id({ commit }, dataset_id) {
        commit("set_dataset_id", dataset_id);
      },
      set_owner_id({ commit }, owner_id) {
        commit("set_owner_id", owner_id);
      },
      set_is_public({ commit }, is_public) {
        commit("set_dataset_id", is_public);
      },
      set_title({ commit }, title) {
        commit("set_title", title);
      },
      set_tsne_is_loading({ commit }, is_loading) {
        commit("set_tsne_is_loading", is_loading);
      },
      set_dataset_info({ commit }, payload) {
        const { dataset_id, title, owner_id, is_public } = payload;
        commit("set_dataset_id", dataset_id);
        commit("set_owner_id", owner_id);
        commit("set_is_public", is_public);
        commit("set_title", title);
      },
      set_plot_type({ commit }, plot_type) {
        commit("set_plot_type", plot_type);
        // When display type is changed, we want to
        // reset these other settings so dependent
        // components
        commit("set_chart_data", {});
        commit("set_gene_symbol", "");
      },
      async fetch_dataset_info({ commit }, dataset_id) {
        commit("set_dataset_id", dataset_id);
        const { title, is_public, owner_id } = await $.ajax({
          url: "./cgi/get_dataset_info.cgi",
          type: "POST",
          data: { dataset_id },
          dataType: "json",
        });
        commit("set_owner_id", owner_id);
        commit("set_is_public", is_public);
        commit("set_title", title);
      },
      async fetch_user_displays({ commit }, { user_id, dataset_id }) {
        const displays = await $.ajax({
          url: "./cgi/get_dataset_displays.cgi",
          type: "POST",
          data: { user_id, dataset_id },
          dataType: "json",
        });
        // Filter out the multigene displays, which do not have the "gene_symbol" config property
        const curated_displays = displays.filter( display => display.plotly_config.hasOwnProperty('gene_symbol'));
        commit("set_user_displays", curated_displays);
      },
      async fetch_owner_displays({ commit }, { owner_id, dataset_id }) {
        const displays = await $.ajax({
          url: "./cgi/get_dataset_displays.cgi",
          type: "POST",
          data: { user_id: owner_id, dataset_id },
          dataType: "json",
        });
        // Filter out the multigene displays, which do not have the "gene_symbol" config property
        const curated_displays = displays.filter( display => display.plotly_config.hasOwnProperty('gene_symbol'));
        commit("set_owner_displays", curated_displays);
      },
      async fetch_available_plot_types(
        { commit, state },
        { user_id, session_id, dataset_id, analysis_id }
      ) {
          // Cancelling last axios call, if applicable
          if (state.cancel_source) {
            state.cancel_source.cancel('Newer "fetch_available_plot_types" call detected.');
          }

          // Create cancel token to cancel last request of this to prevent race condition
          // reference: https://github.com/axios/axios#cancellation
          const CancelToken = axios.CancelToken;
          const cancel_source = CancelToken.source()
          commit('set_plot_type_cancel_source', cancel_source);

        await axios.post(
          `/api/h5ad/${dataset_id}/availableDisplayTypes`,
          {
            user_id,
            session_id,
            dataset_id,
            analysis_id,
          },
          { cancelToken: cancel_source.token
          }).then((response) => {
            commit('set_available_plot_types', response.data);
	        }).catch((thrown) => {
            if (axios.isCancel(thrown)) {
                console.info('Request canceled:', thrown.message);
            } else {
                // handle error
                console.error(thrown);
            }
        });
      },
      async fetch_h5ad_info({ commit }, payload) {
        const { dataset_id, analysis } = payload;

        let data;
        if (analysis) {
          const response = await axios.get(
            `/api/h5ad/${dataset_id}?analysis_id=${analysis.id}`
          );
          data = response.data;
        } else {
          const response = await axios.get(`/api/h5ad/${dataset_id}`);
          data = response.data;
        }
        const { obs_columns, obs_levels } = data;
        commit("set_columns", obs_columns);
        commit("set_levels", obs_levels);
      },
      async fetch_plotly_data({ commit }, { config, plot_type, dataset_id }) {
        commit("set_loading_chart", true);
        const payload = { ...config, plot_type };
        const { data } = await axios.post(`/api/plot/${dataset_id}`, payload);
        commit("set_chart_data", data);

        const {
          plot_colors,
          plot_palette,
          reverse_palette,
          plot_order,
          x_axis,
          y_axis,
          z_axis,
          point_label,
          hide_x_labels,
          hide_y_labels,
          hide_legend,
          color_name,
          facet_row,
          facet_col,
          size_by_group,
          marker_size,
          jitter,
          x_min,
          y_min,
          x_max,
          y_max,
          x_title,
          y_title,
          vlines,
          success,
          message,
        } = data;

        commit("set_order", plot_order);
        commit("set_colors", plot_colors);
        commit("set_color_palette", plot_palette);
        commit("set_reverse_palette", reverse_palette);
        commit("set_x_axis", x_axis);
        commit("set_y_axis", y_axis);
        commit("set_z_axis", z_axis);
        commit("set_point_label", point_label);
        commit("set_hide_x_labels", hide_x_labels);
        commit("set_hide_y_labels", hide_y_labels);
        commit("set_hide_legend", hide_legend);
        commit("set_color_name", color_name);
        commit("set_facet_row", facet_row);
        commit("set_facet_col", facet_col);
        commit("set_size_by_group", size_by_group);
        commit("set_marker_size", marker_size);
        commit("set_jitter", jitter);
        commit("set_x_min", x_min);
        commit("set_y_min", y_min);
        commit("set_x_max", x_max);
        commit("set_y_max", y_max);
        commit("set_x_title", x_title);
        commit("set_y_title", y_title);
        commit("set_vlines", vlines);

        commit("set_success", success);
        commit("set_message", message);
        commit("set_loading_chart", false);
      },
      async fetch_tsne_image(
        { commit },
        { config, plot_type, dataset_id, analysis, analysis_owner_id }
      ) {
        commit("set_tsne_is_loading", true);
        const payload = { ...config, plot_type, analysis, analysis_owner_id };

        const { data } = await axios.post(`/api/plot/${dataset_id}/tsne`, payload);

        commit("set_x_axis", config.x_axis);
        commit("set_y_axis", config.y_axis);
        commit("set_colorize_legend_by", config.colorize_legend_by);
        commit("set_horizontal_legend", config.horizontal_legend);
        commit("set_skip_gene_plot", config.skip_gene_plot);
        commit("set_plot_by_group", config.plot_by_group);
        commit("set_max_columns", config.max_columns);
        commit("set_colors", config.colors);
        commit("set_order", config.order);

        commit("set_image_data", data.image);
        commit("set_success", data.success);
        commit("set_message", data.message);
        commit("set_tsne_is_loading", false);
      },
      set_index({ commit }, index) {
        commit("set_index", index);
      },
      set_order({ commit }, order) {
        commit("set_order", order);
        commit("set_levels", order);
      },
      set_color({ commit }, { name, color }) {
        commit("set_color", { name, color });
      },
      set_colors({ commit }, colors) {
        commit("set_colors", colors);
      },
      set_color_palette({ commit }, palette) {
        commit("set_color_palette", palette);
      },
      set_reverse_palette({ commit }, isReverse) {
        commit("set_reverse_palette", isReverse);
      },
      set_gene_symbol({ commit }, gene_symbol) {
        commit("set_gene_symbol", gene_symbol);
      },
      async fetch_gene_symbols({ commit }, { dataset_id, analysis_id }) {
        const base = `./api/h5ad/${dataset_id}/genes`;
        const query = analysis_id ? `?analysis=${analysis_id}` : "";

        const { data } = await axios.get(`${base}${query}`);
        commit("set_gene_symbols", data.gene_symbols);
      },
      set_label({ commit }, label) {
        commit("set_label", label);
      },
      set_display_data({ commit }, display) {
        let { label, plot_type, plotly_config: config } = display;

        commit("set_label", label);
        commit("set_plot_type", plot_type);
        commit("set_config", config);
      },
      reset({ commit }) {
        commit("set_label", "");
        commit("set_plot_type", null);
        commit("reset_config");
      },
      async fetch_svg_data({ commit }, { gene_symbol, dataset_id }) {
        const { data } = await axios.get(
          `/api/plot/${dataset_id}/svg?gene=${gene_symbol}`
        );
        commit("set_chart_data", data);
      },
      async fetch_default_display({ commit }, { user_id, dataset_id }) {
        const { default_display_id } = await $.ajax({
          url: "./cgi/get_default_display.cgi",
          type: "POST",
          data: { user_id, dataset_id },
          dataType: "json",
        });
        commit("set_default_display_id", default_display_id);
      },
      set_analysis_id({ commit }, analysis) {
        commit("set_analysis_id", analysis);
      },
      set_analysis({ commit }, analysis) {
        commit("set_analysis", analysis);
      },
      update_display({ commit }, display) {
        commit("update_display", display);
      },
      set_image_data({ commit }, image_data) {
        commit("set_image_data", image_data);
      },
      set_success({ commit }, success) {
        commit("set_success", success);
      },
      set_message({ commit }, message) {
        commit("set_message", message);
      },
      remove_display({ commit }, display_id) {
        commit("delete_display", display_id);
      },
      update_default_display_id({ commit }, { display_id }) {
        commit("set_default_display_id", display_id);
      },
    },
  });

  const routes = [
    {
      path: "/dataset/:dataset_id/displays",
      component: datasetCurator,
      props: true,
      children: [
        {
          path: "",
          name: "dashboard",
          component: datasetDisplays,
        },
        {
          path: "new",
          name: "new",
          component: newDisplay,
          beforeEnter(to, from, next) {
            // We want to reset our data that may have been loaded
            // from a previous display
            store.dispatch("reset");
            next();
          },
        },
        {
          path: ":display_id/edit",
          name: "edit",
          component: datasetDisplay,
          props: true,
        },
      ],
    }
  ];

  const router = new VueRouter({
    routes,
  });

  const app = new Vue({
    el: "#app",
    router,
    store,
    computed: {
      ...Vuex.mapState(["user"]),
    },
    created() {
      // We want to check for session when the curator app is first created
      sleep(500).then(() => {
        // If CURRENT_USER is defined at this point, add information as placeholder test
        if (CURRENT_USER) {
          this.$store.commit("set_user", CURRENT_USER);
        }
      });
    },
  });
};
