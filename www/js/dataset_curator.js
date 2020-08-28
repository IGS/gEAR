'use strict';

/*
 This script relies on the source having also included the
 common.js within this project (for login purposes)
*/

Vue.use(Vuex);
Vue.use(VueRouter);
Vue.use(bootstrapVue);

(() => {

  const datasetTitle = Vue.component('dataset-title', {
    template: `
      <b-row v-if='title' class='justify-content-md-center'>
        <div class="mt-5 col-12">
        <!--
          <h2 class='font-weight-light'>You are curating <span class='font-weight-bold'>{{ title }}</span></h2>
          <hr />
          -->
        </div>
      </b-row>
    `,
    computed: {
      ...Vuex.mapState(['dataset_id', 'title']),
    },
  });

  const addDisplayBtn = Vue.component('add-display-btn', {
    template: `
      <b-card
        @mouseover="isHovering = true"
        @mouseout="isHovering = false"
        @click="routeToAddDisplay"
        :class='classed'
      >
       <div class='mt-4 pt-5 card-body text-center'>
         <div class='display-4'>
           <i class="fas fa-plus"></i>
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
      ...Vuex.mapState(['dataset_id']),
      classed() {
        return {
          hovering: this.isHovering,
          'display-card': true,
          elevation: true,
          'border-0': true,
        };
      },
    },
    methods: {
      routeToAddDisplay() {
        this.$router.push(`displays/new`);
      },
    },
  });

  const tsneChart = Vue.component('tsne-chart', {
    template: `
     <div>
        <div v-if='is_loading' class='elevation border-0 mb-5 sticky-chart'>
          <div>
            <img class='ml-5 m-5' style="width:50px;height:50px;" src='../img/loading_search.gif'></img>
          </div>
        </div>
        <div v-if='plot_params_ready || preconfigured' ref='chart'>
          <img id='tsne_preview' class='img-fluid'></img>
        </div>
      </div>
    `,
    props: {
        data: {
            default: null,
        },
        preconfigured: {
            default: false,
        },
        display_data: {}
    },
    computed: {
      ...Vuex.mapState(['dataset_id', 'analysis', 'config', 'image_data', 'tsne_is_loading']),
      is_loading() {
          if (this.tsne_is_loading == true) {
              return true;
          } else {
              return false;
          }
      },
      plot_params_ready() {
          if (this.config.x_axis && this.config.y_axis) {
            return true;
          } else {
            return false;
          }
      }
    },
    created() {
      if (this.preconfigured) {
        const config = this.display_data.plotly_config;
        const dataset_id = this.dataset_id;
        this.fetch_tsne_image({config, dataset_id});
      }
    },
    watch: {
      image_data() {
        if (this.image_data) {
          $("#tsne_preview").attr('src', "data:image/png;base64," + this.image_data);
        }
      },
    },
    methods: {
      ...Vuex.mapActions(['fetch_tsne_image', ]),
      get_image_data() {
        // I don't believe this is used
      }
    },
  });

  const svgChart = Vue.component('svg-chart', {
    template: `
      <div>
        <div v-if='!display_data' class='elevation border-0 mb-5'>
          <b-card-body no-body
            class='elevation border-0 mb-5'>
            <div ref='chart'></div>
          </b-card-body>
        </div>
        <div v-else ref='chart' style='height:230px'></div>
      </div>
    `,
    props: {
      chart_data: {
        default: null,
      },
      low_color: {
        default: '',
      },
      high_color: {
        default: '',
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
        scoring_method: 'gene', // need to make this a toggle
      };
    },
    computed: {
      ...Vuex.mapState(['dataset_id']),
    },
    watch: {
      async svg(svg) {
        this.paths = svg.selectAll('path, circle');
        svg.select('svg').attr({
          width: '100%',
          height: this.display_data ? '200px' : '',
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
      high_color() {
        this.color_svg();
      },
      chart_data() {
        this.color_svg();
      },
    },
    async created() {
      const svg_path = `datasets_uploaded/${this.dataset_id}.svg`;
      Snap.load(svg_path, svg => {
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
        let chart_data, low_color, high_color;
        if (data) {
          const { plotly_config } = this.display_data;
          const { colors } = { ...plotly_config };
          low_color = colors.low_color;
          high_color = colors.high_color;
          chart_data = data;
        } else {
          chart_data = this.chart_data;
          low_color = this.low_color;
          high_color = this.high_color;
        }

        const score = chart_data.scores[this.scoring_method];
        const paths = this.paths;
        const { data: expression } = chart_data;
        if (
          this.scoring_method === 'gene' ||
          this.scoring_method === 'dataset'
        ) {
          const { min, max } = score;
          const color = d3
            .scaleLinear()
            .domain([min, max])
            .range([low_color, high_color]); // these colors should be stored in config

          const tissues = Object.keys(chart_data.data);

          paths.forEach(path => {
            const tissue = path.node.className.baseVal;
            if (tissues.includes(tissue)) {
              path.attr('fill', color(expression[tissue]));
            }
          });
        }
      },
    },
  });

  const plotlyChart = Vue.component('plotly-chart', {
    template: `
      <div>
        <div v-if='!display_data' class='elevation border-0 mb-5 sticky-chart'>
          <div v-if='!is_there_data'>
            <img class='ml-5 m-5' style="width:50px;height:50px;" src='../img/loading_search.gif'></img>
          </div>
          <div v-else ref='chart'>
            <img class='img-fluid' v-if='img' :src='imgData'></img>
          </div>
        </div>
        <div v-else ref='chart'>
          <div v-if='!imgData' class='col align-middle text-center mt-5 pt-4'>
            <p class=''>Loading</p>
          </div>
          <img v-else class='img-fluid' v-if='img' :src='imgData'></img>
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
        imgData: '',
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
    watch: {
      data() {
        this.draw_chart();
      },
    },
    mounted() {
      if (this.is_there_data) this.draw_chart();
    },
    computed: {
      ...Vuex.mapState(['dataset_id']),
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
        if (data) {
          this.loading = false;
          const { plot_json, plot_config } = data;
          if (this.img) {
            Plotly.toImage({ ...plot_json, plot_config }).then(url => {
              this.imgData = url;
            });
          } else {
            Plotly.newPlot(this.$refs.chart, { ...plot_json, plot_config });
          }
        } else {
          const { plot_json, plot_config } = this.data;
          if (this.img) {
            Plotly.toImage({ ...plot_json, plot_config }).then(url => {
              this.imgData = url;
            });
          } else {
            Plotly.newPlot(this.$refs.chart, { ...plot_json, plot_config });
          }
        }
      },
    },
  });

  const userDisplays = Vue.component('user-displays', {
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
                  <div
                    slot='header'
                  >
                  <p class='float-right'>
                  <b-badge
                    pill
                    class='purple'
                    variant="dark">{{ display_data.plotly_config.gene_symbol }}</b-badge>
                    <b-badge
                      pill
                      class='purple'
                      v-if='display_data.is_default' variant="dark">Default</b-badge></p>
                  </div>
                      <div v-if="display_data.plot_type === 'svg'">
                        <svg-chart :display_data='display_data'></svg-chart>
                      </div>
                      <div v-else-if="display_data.plot_type === 'tsne_static'">
                        <tsne-chart
                          :display_data='display_data'
                          :preconfigured=true
                        ></tsne-chart>
                      </div>
                      <plotly-chart v-else :display_data='display_data' :img='true'></plotly-chart>
                  <div
                    slot='footer'
                    class='text-right'
                  >
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
      ...Vuex.mapState(['user', 'dataset_id', 'default_display_id', 'config']),
      ...Vuex.mapGetters(['user_displays']),
      get_config_analysis_id() {
          if (typeof this.config.analysis == 'undefined' || this.config.analysis == null) {
              return null;
          } else {
              return this.config.analysis.id;
          }
      },
      user_id() {
        return this.user.id;
      },
    },
    methods: {
      ...Vuex.mapActions([
        'fetch_user_displays',
        'fetch_default_display',
        'remove_display',
        'update_default_display_id',
      ]),
      edit_display(display_id) {
        // this.$router.replace(`/dataset/${this.dataset_id}/displays/${display_id}/edit`)
        this.$router.push(`displays/${display_id}/edit`);
      },
      get_default_display() {
        const user_id = this.user_id;
        const dataset_id = this.dataset_id;

        return $.ajax({
          url: './cgi/get_default_display.cgi',
          type: 'POST',
          data: { user_id, dataset_id },
          dataType: 'json',
        });
      },
      async save_as_default(display_id) {
        const payload = {
          display_id,
          user_id: this.user.id,
          dataset_id: this.dataset_id,
        };
        await $.ajax({
          url: './cgi/save_default_display.cgi',
          type: 'POST',
          data: payload,
          dataType: 'json',
        });

        this.update_default_display_id({ display_id });
      },
      async delete_display(display_id) {
        const payload = {
          id: display_id,
          user_id: this.user_id,
        };

        const res = await $.ajax({
          url: './cgi/delete_dataset_display.cgi',
          type: 'POST',
          data: payload,
          dataType: 'json',
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
  const ownerDisplays = Vue.component('owner-displays', {
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
                <div
                  slot='header'
                >
                <p class='float-right'>
                <b-badge
                  pill
                  class='purple'
                  variant="dark">{{ display_data.plotly_config.gene_symbol }}</b-badge>
                  <b-badge
                    pill
                    class='purple'
                    v-if='display_data.is_default' variant="dark">Default</b-badge></p>
                </div>
                  <div v-if="display_data.plot_type === 'svg'">
                    <svg-chart :display_data='display_data'></svg-chart>
                  </div>
                  <div v-else-if="display_data.plot_type === 'tsne_static'">
                    <tsne-chart
                      :display_data='display_data'
                      :preconfigured=true
                    ></tsne-chart>
                    </div>
                    <div v-else>
                      <plotly-chart :display_data='display_data' :img='true'></plotly-chart>
                    </div>
                <div
                  slot='footer'
                  class='text-right'
                >
                  <b-button
                    size='sm'
                    variant='link'
                    style='color:#562a6f'
                    @click='save_as_default(display_data.id)'
                  >Make Default</b-button>
                </div>
            </b-card>
            </div>
        </div>
      </div>
    </transition>
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
        'user',
        'owner_id',
        'dataset_id',
        'default_display_id',
        'config'
      ]),
      ...Vuex.mapGetters(['is_user_owner', 'owner_displays']),

    },
    methods: {
      ...Vuex.mapActions([
        'fetch_owner_displays',
        'fetch_user_displays',
        'fetch_default_display',
        'update_default_display_id',
      ]),
      async save_as_default(display_id) {
        const payload = {
          display_id,
          user_id: this.user.id,
          dataset_id: this.dataset_id,
        };

        await $.ajax({
          url: './cgi/save_default_display.cgi',
          type: 'POST',
          data: payload,
          dataType: 'json',
        });

        this.update_default_display_id({ display_id });
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

  const datasetCurator = Vue.component('dataset-curator', {
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
    props: ['dataset_id', 'user'],
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
      ...Vuex.mapActions(['fetch_dataset_info']),
    },
  });

  const datasetDisplays = Vue.component('dataset-displays', {
    template: `
      <div>
        <user-displays></user-displays>
        <owner-displays v-if='owner_id && !is_user_owner'></owner-displays>
      </div>
    `,
    computed: {
      ...Vuex.mapState(['owner_id']),
      ...Vuex.mapGetters(['is_user_owner']),
    },
    components: {
      userDisplays,
      ownerDisplays,
    },
  });

  const chooseDisplayType = Vue.component('choose-display-type', {
    template: `
      <div>
        <b-row>
        <b-col>
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
    computed: {
      ...Vuex.mapState([
        'user',
        'plot_type',
        'dataset_id',
        'available_plot_types',
        'dataset_type',
      ]),
      display_options() {
        const display_options = Object.keys(this.available_plot_types)
          .filter(type => this.available_plot_types[type])
          .map((type, i) => type);

        if (
          this.dataset_type == 'primary' &&
          this.available_plot_types['tsne']
        ) {
          return display_options.filter(display => display !== 'tsne');
        } else {
          return display_options;
        }
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

      this.fetch_available_plot_types({ user_id, session_id, dataset_id });
    },
    methods: {
      ...Vuex.mapActions(['fetch_available_plot_types', 'set_plot_type']),
        update_plot_type(plot_type) {
        this.set_plot_type(plot_type);
        this.$emit('input', plot_type);
      },
    },
  });

  const PlotlyArguments = Vue.component('plotly-arguments', {
    template: `
      <div>
        <b-row>
          <b-col>
          <h3>Display Parameters</h3>
          <b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">X</label>
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
                <b-form-select :options='columns' v-model='y_axis' size='sm'>
                  <template slot="first">
                  <option value="raw_value">expression</option>
                  </template>
                </b-form-select>
                <b-form-checkbox v-model='hide_y_labels'>
                  Hide Y Tickmarks
                </b-form-checkbox>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Color</label>
                <b-form-select :options='columns' v-model='color_name' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                <option value="raw_value">expression</option>
                </b-form-select>
                <!--<b-link v-if='color' @click.stop='color_modal = !color_modal'><small class="text-muted pull-right primary">EDIT COLORS</small></b-link></p>
                    <b-modal v-model="color_modal">
                        <div v-for='level in color_levels'>
                          <label>{{ level }} </label>
                          <swatches colors="text-advanced" popover-to="left" ></swatches>
                        </div>
                    </b-modal>-->
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Label</label>
                <b-form-select :options="columns" v-model='point_label' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Facet Row</label>
                <b-form-select :options="columns" v-model='facet_row' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right">
                <label class="mb-0">Facet Column</label>
                <b-form-select :options='columns' v-model='facet_col' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>
              <b-form-group label-align-sm="right" v-if="plot_type === 'scatter'">
                <label class="mb-0">Marker Size (sm <--> lg)</label>
                <b-form-input type="range" v-model="marker_size" min="1" max="15" number></b-form-input>
              </b-form-group>
              <b-form-group label-align-sm="right" v-if="plot_type === 'scatter'">
                <label class="mb-0">Jitter (none <--> lots) </label>
                <b-form-input type="range" v-model="jitter" min="0" max="1" step="0.25" number></b-form-input>
              </b-form-group>
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
    components: {
        swatches: window['vue-swatches'],
    },
    data() {
      return {
        x_axis: null,
        y_axis: 'raw_value',
        point_label: null,
        color_name: null,
        facet_row: null,
        facet_col: null,
        marker_size: 3,
        jitter: 0,
        color_modal: false,
        hide_x_labels: false,
        hide_y_labels: false,
      };
    },
    computed: {
      ...Vuex.mapState([
        'dataset_id',
        'config',
        'columns',
        'levels',
        'plot_type',
      ]),
      color_levels() {
        return this.levels[this.color];
      },
    },
    created() {
      if ('x_axis' in this.config) this.x_axis = this.config.x_axis;
      if ('y_axis' in this.config && this.config.y_axis) {
        // we only want to change y_axis if there's a value other than null so
        // it defaults to "raw_value" set in data above
        this.y_axis = this.config.y_axis;
      }
      if ('hide_x_labels' in this.config)
        this.hide_x_labels = this.config.hide_x_labels;
      if ('hide_y_labels' in this.config)
        this.hide_y_labels = this.config.hide_y_labels;
      if ('color_name' in this.config) this.color_name = this.config.color_name;
      if ('facet_row' in this.config) this.facet_row = this.config.facet_row;
      if ('facet_col' in this.config) this.facet_col = this.config.facet_col;
      if ('marker_size' in this.config && this.config.marker_size) {
        // we only want to change marker_size if there's a value other than null so
        // it defaults to the default size (3)
        this.marker_size = this.config.marker_size;
      }
      if ('jitter' in this.config && this.config.jitter) {
        // see 'marker_size' comment
        this.jitter = this.config.jitter;
      }
      if ('point_label' in this.config) this.point_label = this.config.point_label;
    },
    watch: {
      x_axis(val) {
        if (this.y_axis === val) this.y_axis = null;
        if (this.point_label === val) this.point_label = null;
        if (this.color_name === val) this.color_name = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      y_axis(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.point_label === val) this.point_label = null;
        if (this.color_name === val) this.color_name = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      point_label(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.color_name === val) this.color_name = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      color_name(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.point_label === val) this.point_label = null;
        if (this.facet_row === val) this.facet_row = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      facet_row(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.point_label === val) this.point_label = null;
        if (this.color_name === val) this.color_name = null;
        if (this.facet_col === val) this.facet_col = null;
      },
      facet_col(val) {
        if (this.x_axis === val) this.x_axis = null;
        if (this.y_axis === val) this.y_axis = null;
        if (this.point_label === val) this.point_label = null;
        if (this.color_name === val) this.color_name = null;
        if (this.facet_row === val) this.facet_row = null;
      },
    },
    methods: {
      ...Vuex.mapActions(['fetch_plotly_data']),
      preview() {
        const { gene_symbol, analysis, colors } = this.config;

        const config = {
          gene_symbol,
          analysis,
          colors,
          x_axis: this.x_axis,
          y_axis: this.y_axis,
          point_label: this.point_label,
          color_name: this.color_name,
          facet_row: this.facet_row,
          facet_col: this.facet_col,
          marker_size: this.marker_size,
          jitter: this.jitter,
          hide_x_labels: this.hide_x_labels,
          hide_y_labels: this.hide_y_labels,
        };

        const plot_type = this.plot_type;
        const dataset_id = this.dataset_id;
        this.fetch_plotly_data({ config, plot_type, dataset_id });
      },
    },
  });

  const groupBy = Vue.component('group-by', {
    template: `
    <b-card-body no-body
      class='elevation border-0 mt-5'>
      <b-row class='m-4'>
        <b-col>
          <h3>Group By</h3>
          <p>Choose your columns to group by. When there are 3 selections, the display will subplot based on the first selection. The second selection is the column you are coloring, and the third option is what is on the x-axis. For line plots, the time_point column must be the last selection for the x-axis.</p>
        </b-col>
        <b-col>
          <b-form>
            <b-form-radio-group
              v-model='selected'
              buttons
              button-variant='light'
              :options="group_by_n"
            />
            <b-form-select
              v-for='(idx, i) in indices'
              :key='idx'
              class="mt-2 mb-2 mr-sm-2 mb-sm-0 group_by"
              :value="idx"
              @input="update_index"
              :options="columns"
              :disabled="idx === 'time_point' && plot_type === 'line'"
            >
              <option slot="first" :value="null">Choose index...</option>
            </b-form-select>
            <b-button
              @click='preview'
              variant="primary"
              class='mt-4 btn-purple float-right'>
                Preview Chart
            </b-button>
          </b-form>
        </b-col>
    </b-row>
    </b-card-body>
    `,
    data() {
      return {
        selected: 0,
      };
    },
    watch: {
      selected(new_value, old_value) {
        // when user changes selected, we need to reset
        // the index...
        this.set_index([]);
      },
    },
    computed: {
      ...Vuex.mapState([
        'dataset_id',
        'config',
        'columns',
        'levels',
        'plot_type',
      ]),
      group_by_n() {
        const [_, ...groups] = [...Array(this.columns.length + 1).keys()];
        return groups;
      },
      indices() {
        // TODO -- This code below is why after you group by, we cannot
        // go back up...fix this some way....

        // if we have an index, we want to show those selections
        let indices =
          this.config.index.length !== 0 ? this.config.index : this.columns;

        // force line plots to have time_point as the last axis
        if (this.plot_type === 'line') {
          indices = indices.filter(i => i !== 'time_point');
          indices.push('time_point');
          return indices.slice(indices.length - this.selected, indices.length);
        } else {
          return indices.slice(0, this.selected);
        }
      },
    },
    created() {
      // Make sure our select inputs do not exceed 3
      if (this.config.index.length !== 0) {
        this.selected =
          this.config.index.length <= 3 ? this.config.index.length : 3;
      } else {
        this.selected = this.columns.length <= 3 ? this.columns.length : 3;
      }
    },
    methods: {
      ...Vuex.mapActions(['set_index', 'fetch_plotly_data']),
      update_index() {
        const index = this.index_to_group_by();
        this.$emit('input', index);
      },
      index_to_group_by() {
        return [...document.querySelectorAll('.group_by')]
          .map(el => el.value)
          .filter(el => el);
      },
      preview() {
        const indices = this.index_to_group_by();
        this.set_index(indices);

        const { gene_symbol, index } = this.config;
        const config = { gene_symbol, index };
        const plot_type = this.plot_type;
        const dataset_id = this.dataset_id;

        this.fetch_plotly_data({ config, plot_type, dataset_id });
      },
    },
  });

  const displayNameInput = Vue.component('display-name-input', {
    template: `
    <div>
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
      ...Vuex.mapState(['label']),
    },
    methods: {
      ...Vuex.mapActions(['set_label']),
      update_display_name(label) {
        this.set_label(label);
      },
    },
  });

  const saveDisplayBtn = Vue.component('save-display-btn', {
    template: `
    <b-form-group>
      <b-button @click='save' class='btn-purple float-right'>Save Display</b-button/>
    </b-form-group>
        <!--<b-row class='mt-2'>
          <b-col></b-col>
          <b-col>
            <b-alert
              :show="dismissCountDown"
              dismissible
              fade
              :variant="variant"
              @dismissed="dismissCountDown=0"
              @dismiss-count-down="countDownChanged">
                <p v-if="variant=='success'">
                  Display Saved!
                </p>
                <p v-if="variant=='danger'">
                  Something went wrong...
                </p>
            </b-alert>
          </b-col>
        </b-row>-->
      `,
    props: ['display_id'],
    data() {
      return {
        dismissSecs: 3,
        variant: '',
        dismissCountDown: 0,
        showDismissibleAlert: false,
      };
    },
    computed: {
      ...Vuex.mapState(['user', 'dataset_id', 'plot_type', 'config', 'label']),
    },
    methods: {
      ...Vuex.mapActions(['update_display']),
      countDownChanged(dismissCountDown) {
        this.dismissCountDown = dismissCountDown;
      },
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
          url: './cgi/save_dataset_display.cgi',
          type: 'POST',
          data: payload,
          dataType: 'json',
        });

        if (res && res.success) {
          // this.variant = 'success';
          // this.dismissCountDown = this.dismissSecs;
          if (this.display_id) {
            this.update_display(payload);
          }
          this.$router.push(`/dataset/${this.dataset_id}/displays`);
        } else {
          this.variant = 'danger';
          this.dismissCountDown = this.dismissSecs;
        }
      },
    },
  });

  const geneSymbolInput = Vue.component('gene-symbol-input', {
    template: `
      <div>
        <b-row>
          <b-col>
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
    props: ['analysis'],
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
      ...Vuex.mapState(['dataset_id', 'config', 'gene_symbols']),
      is_gene_available() {
        return this.gene_symbols
          .map(gene => gene.toLowerCase())
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
        this.$emit('gene-updated', true);
      }
    },
    methods: {
      ...Vuex.mapActions(['set_gene_symbol', 'fetch_gene_symbols']),
      update_gene_symbol(gene_symbol) {
        this.set_gene_symbol(gene_symbol);
        this.$emit('gene-updated', this.is_gene_available);
      },
    },
  });

  // TODO: When new display parameters are updated, the default order is initially used (not passed in plotly_data API post)
  const displayOrder = Vue.component('display-order', {
    template: `
      <div>
        <b-row>
          <b-col>
            <h4>Order</h4>
            <b-row v-for="variable in order">
              <b-col>
              <b-list-group class="mb-2">
                <h4>{{ variable.key }}</h4>
                <draggable v-model="variable.value" @end="reorder_display">
                  <b-list-group-item v-for="val in variable.value" :key="val" class="draggable">{{ val }}</b-list-group-item>
                </draggable>
                <!--<b-list-group-item v-b-toggle='variable.key' class="d-flex fluid justify-content-between align-items-center">
                    <strong>{{ variable.key }}</strong>
                    <i class="fas fa-chevron-down"></i>
                </b-list-group-item>
                  <b-collapse :id="variable.key" :accordion="variable.key" role="tabpanel">
                    <draggable v-model="variable.value" @end="reorder_display">
                      <b-list-group-item v-for="val in variable.value" :key="val" class="draggable">{{ val }}</b-list-group-item>
                    </draggable>
                  </b-collapse>
                -->
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
        orderAccordian: {},
      };
    },
    computed: {
      ...Vuex.mapState(['levels', 'config', 'plot_type', 'dataset_id']),
    },
    created() {
      this.get_order();
    },
    methods: {
      ...Vuex.mapActions(['set_order', 'fetch_plotly_data']),
      get_order() {
        const keys = Object.keys(this.config.order);
        const order = keys.map(key => {
          return {
            key,
            value: [...this.config.order[key]],
          };
        });
        this.order = order;
      },
      reorder_display() {
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
    },
  });

  const displayColors = Vue.component('display-colors', {
    template: `
    <div>
      <b-row>
        <b-col>
          <h3>Color</h3>
          <b-row v-for="{name, color} in colors_array">
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
    computed: {
      ...Vuex.mapState(['config']),
      colors_array() {
        return Object.entries(this.config.colors).map(([key, val]) => {
          return {
            name: key,
            color: val,
          };
        });
      },
    },
    created() {
      // console.log(this.config.colors);
    },
    methods: {
      ...Vuex.mapActions(['set_color']),
      update_color(name, color) {
        this.set_color({ name, color });
      },
    },
  });

  const barDisplay = Vue.component('bar-display', {
    template: `
      <div>
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
    props: ['display_id'],
    data() {
      return {
        is_gene_available: true,
      };
    },
    computed: {
      ...Vuex.mapState(['config']),
      is_there_data_to_save() {
        return (
          'x_axis' in this.config &&
          // && 'colors' in this.config // && Object.entries(this.config.colors).length !== 0
          // && 'order' in this.config // && Object.entries(this.config.order).length !== 0
          'gene_symbol' in this.config &&
          this.config.gene_symbol !== ''
        );
      },
    },
    components: {
      PlotlyArguments,
      groupBy,
      geneSymbolInput,
      displayNameInput,
      displayOrder,
      displayColors,
      saveDisplayBtn,
    },
  });

  const lineDisplay = Vue.component('line-display', {
    extends: barDisplay,
  });

  const violinDisplay = Vue.component('violin-display', {
    extends: barDisplay,
  });

  const scatterDisplay = Vue.component('scatter-display', {
    extends: barDisplay,
  });

  const tsnePlotlyDisplay = Vue.component('tsne-plotly-display', {
    extends: scatterDisplay,
  });

  const plotlyDisplay = Vue.component('plotly-display', {
    template: `
      <div v-if="plot_type === 'bar'">
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
      <div v-else-if="plot_type === 'line'">
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
      <div v-else-if="plot_type === 'violin'">
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
      <div v-else-if='plot_type === "scatter"'>
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
      <div v-else-if='plot_type === "tsne_dynamic"'>
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
      tsnePlotlyDisplay
    },
    data() {
      return {
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState(['dataset_id', 'plot_type', 'config', 'chart_data']),
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
      ...Vuex.mapActions(['fetch_h5ad_info', 'fetch_plotly_data']),
      update_color({ name, color }) {
        const { data } = this.chart_data.plot_json;
        data
          .filter(el => el.name === name)
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

  const svgDisplay = Vue.component('svg-display', {
    template: `
      <div>
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
      ...Vuex.mapState(['dataset_id', 'config', 'chart_data']),
      low_color: {
        get() {
          return this.config.colors.low_color;
        },
        set(color) {
          this.set_color({ name: 'low_color', color });
        },
      },
      high_color: {
        get() {
          return this.config.colors.high_color;
        },
        set(color) {
          this.set_color({ name: 'high_color', color });
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
        return 'colors' in this.config && 'gene_symbol' in this.config;
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
      ...Vuex.mapActions(['fetch_svg_data', 'set_color']),
      preview() {
        const { gene_symbol } = this.config;
        const dataset_id = this.dataset_id;

        this.fetch_svg_data({ gene_symbol, dataset_id });
      },
    },
  });

  const chooseAnalysis = Vue.component('choose-analysis', {
    template: `
      <b-card-body no-body
        class='elevation border-0 mt-5'>
      <b-row class='m-4'>
        <b-col>
          <h3>Choose Analyses</h3>
        </b-col>
        <b-col>
          <b-form>
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
      ...Vuex.mapState(['dataset_id', 'config']),
      analysis() {
        if (this.private_or_public === 'Public') {
          return this.public_labels.find(
            el => el.value === this.config.analysis_id
          );
        } else {
          return this.private_labels.find(
            el => el.value === this.config.analysis_id
          );
        }
      },
      ana_private_or_public() {
        // If an analaysis id is passed,
        // check if its public or private
        return this.public.map(ana => ana.id).includes(this.config.analysis_id)
          ? 'Public'
          : 'Private';
      },
      private_labels() {
        return this.private.map(ana => {
          return {
            value: ana.id,
            text: ana.label,
          };
        });
      },
      public_labels() {
        return this.public.map(ana => {
          return {
            value: ana.id,
            text: ana.label,
          };
        });
      },
      analyses() {
        if (this.private_or_public === 'Public') {
          return this.public_labels;
        } else {
          return this.private_labels;
        }
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
      ...Vuex.mapActions(['set_analysis_id']),
      analysis_selected(analysis) {
        this.set_analysis_id(analysis);
      },
    },
  });

  const tsneArguments = Vue.component('tsne-arguments', {
    template: `
      <div>
        <b-row>
          <b-col>
          <h3>Display Parameters</h3>
          <b-form-group>

              <b-form-group label-align-sm="right">
                <label class="mb-0">X</label>
                <b-form-select :options='columns' v-model='x_axis' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>

              <b-form-group label-align-sm="right">
                <label class="mb-0">Y</label>
                <b-form-select :options='columns' v-model='y_axis' size='sm'>
                  <template slot="first">
                  <option value="raw_value">expression</option>
                  </template>
                </b-form-select>
              </b-form-group>

              <b-form-group label-align-sm="right">
                <b-form-checkbox v-model='show_colorized_legend'>
                  Show colorized legend
                </b-form-checkbox>
              </b-form-group>

              <b-form-group v-if='show_colorized_legend' label-align-sm="right">
                <label class="mb-0">Colorize legend by:</label>
                <b-form-select :options="columns" v-model='colorize_legend_by' size='sm'>
                  <template slot="first">
                    <option :value="null"></option>
                  </template>
                </b-form-select>
              </b-form-group>

            </b-form-group>
          </b-col>
        </b-row>
        <hr>
      </div>
    `,
    components: {

    },
    data() {
      return {
          show_colorized_legend: false,
      };
    },
    computed: {
        ...Vuex.mapState(['dataset_id','config','columns','plot_type','image_data','tsne_is_loading']),
        colorize_legend_by: {
            get () {
                return this.$store.state.config.colorize_legend_by;
            },
            set (value) {
                this.$store.commit('set_colorize_legend_by', value);
            }
        },
        x_axis: {
          get () {
              return this.$store.state.config.x_axis;
          },
          set (value) {
              this.$store.commit('set_x_axis', value);
          }
        },
        y_axis: {
          get () {
              return this.$store.state.config.y_axis;
          },
          set (value) {
              this.$store.commit('set_y_axis', value);
          }
        }
    },
    created() {
        if ('x_axis' in this.config) this.x_axis = this.config.x_axis;
        if ('y_axis' in this.config) this.y_axis = this.config.y_axis;
        if ('colorize_legend_by' in this.config) this.colorize_legend_by = this.config.colorize_legend_by;

        this.fetch_h5ad_info({
            dataset_id: this.dataset_id,
            analysis: this.config.analysis,
        });
    },
    watch: {
        show_colorized_legend: function (val, oldval) {
            // if deselected, clear colorize legend select box
            if (this.show_colorized_legend === false) this.colorize_legend_by = null;
        },
        colorize_legend_by: function (newval, oldval) {
          if (newval != oldval && this.plot_params_ready()) {
            this.draw_image();
          }
        },
        x_axis: function (newval, oldval) {
          if (newval != oldval && this.plot_params_ready()) {
            this.draw_image();
          }
        },
        y_axis: function (newval, oldval) {
          if (newval != oldval && this.plot_params_ready()) {
            this.draw_image();
          }
        },
    },
    methods: {
      ...Vuex.mapActions(['fetch_h5ad_info', 'fetch_plotly_data', 'set_image_data', 'set_tsne_is_loading']),
      plot_params_ready() {
        if (this.x_axis && this.y_axis) {
          return true;
        } else {
          return false;
        }
      },
      get_image_data(gene_symbol) {
        // then craziness: https://stackoverflow.com/a/48980526
        // shift this out when the fetch_tsne_image method is done in Vuex
        return axios.get(`/api/plot/${this.dataset_id}/tsne`, {
            params: {
                gene: this.config.gene_symbol,
                analysis: this.analysis_id,
                colorize_by: this.colorize_legend_by,
                x_axis: this.x_axis,
                y_axis: this.y_axis,
                analysis_owner_id: this.user_id,
                colors: this.colors,
                // helps stop caching issues
                timestamp: new Date().getTime()
            }
        }).then(
          response => {return response.data}
        )
      },
      draw_image() {
        this.set_tsne_is_loading(true);
        this.get_image_data(this.config.gene_symbol).then(
          data => {
            this.set_tsne_is_loading(false);
            this.set_image_data(data);
          }
        );
      }
    },
  });

  const tsneDisplay = Vue.component('tsne-display', {
    template: `
      <div>
        <gene-symbol-input
          v-model='config.gene_symbol'
          :analysis='config.analysis ? config.analysis.id : null'
          @gene-updated='is_gene_available = $event'
        ></gene-symbol-input>
        <transition name="fade" mode="out-in">
          <tsne-arguments v-if="config.gene_symbol !== '' && is_gene_available" />
        </transition>

        <!-- We don't need a preview button, we can immediately show after gene selection
        <transition name="fade" mode="out-in">
              <b-button
                        v-if='is_gene_available'
                        @click="show_tsne = true"
                        class='mt-4 btn-purple'>
                        Preview tSNE
              </b-button>
        </transition>
        -->
        <transition name="fade" mode="out-in">
        <display-name-input
          v-if='is_gene_available'
        ></display-name-input>
      </transition>
      <!--<transition name="fade" mode="out-in">
        <display-colors
          v-if="is_gene_available && config.colorize_legend_by"
        ></display-colors>
      </transition>-->
      <transition name="fade" mode="out-in">
        <save-display-btn
          v-if="is_gene_available"
          :display_id='display_id'
        ></save-display-btn>
      </transition>
      </div>
    `,
    props: ['display_id'],
    components: {
      chooseAnalysis,
      geneSymbolInput,
      tsneChart,
      displayNameInput,
      displayColors,
      saveDisplayBtn,
    },
    data() {
      return {
        is_gene_available: false,
        show_tsne: false,
        loading: false,
      };
    },
    computed: {
      ...Vuex.mapState(['dataset_id', 'config', 'dataset_type', 'analysis']),
    },
  });

  const primaryConfig = Vue.component('primary-config', {
    template: `
      <div>
        <choose-display-type></choose-display-type>
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
    },
    computed: {
      ...Vuex.mapState(['config', 'plot_type', 'dataset_type']),
      is_type_plotly() {
        return (
          this.plot_type == 'bar' ||
          this.plot_type == 'scatter' ||
          this.plot_type == 'line' ||
          this.plot_type == 'violin' ||
          this.plot_type == 'tsne_dynamic'
        );
      },
      is_type_svg() {
        return this.plot_type === 'svg';
      },
      is_type_tsne() {
        return (
            this.plot_type === 'tsne_static' ||
            this.plot_type === 'tsne'
        );
      },
    },
  });

  const chooseStoredAnalysis = Vue.component('choose-stored-analysis', {
    template: `
       <div>
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
      ...Vuex.mapActions(['set_analysis']),
    },
    computed: {
      ...Vuex.mapState(['dataset_id', 'config']),
    },
  });

  const storedAnalysisConfig = Vue.component('stored-analysis-config', {
    template: `
      <div>
        <choose-stored-analysis></choose-stored-analysis>
        <hr v-if="analysis">
        <primary-config v-if="analysis"></primary-config>
      </div>
    `,
    components: { chooseStoredAnalysis, primaryConfig },
    computed: {
      ...Vuex.mapState(['analysis']),
    },
  });

  const configurationPanel = Vue.component('configuration-panel', {
    template: `
      <div>
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
              <stored-analysis-config></stored-analysis-config>
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
        selected: 'primary',
      };
    },
    computed: {
      ...Vuex.mapState(['dataset_type']),
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
      ...Vuex.mapActions(['set_dataset_type']),
    },
  });

  const newDisplay = Vue.component('new-display', {
    template: `
      <b-container fluid>
        <b-row>
          <b-col cols='2'>
            <configuration-panel></configuration-panel>
          </b-col>
          <b-col cols='10'>
            <plotly-chart v-if='is_type_plotly && is_there_data_to_draw' :data='chart_data' class="sticky-chart"></plotly-chart>
            <svg-chart
              v-else-if='is_type_svg && is_there_data_to_draw'
              :chart_data='chart_data'
              class="sticky-chart"
              :low_color="config.colors.low_color"
              :high_color="config.colors.high_color"/>
            <tsne-chart
              v-else-if='is_type_tsne && gene_selected'
              display
              :gene_symbol='config.gene_symbol'
            ></tsne-chart>
          </b-col>
        </b-row>
      </b-container>
    `,
    components: {
      chooseDisplayType,
      plotlyDisplay,
      plotlyChart,
      svgChart,
      tsneChart,
      svgDisplay,
      tsneDisplay,
      configurationPanel,
    },
    computed: {
        ...Vuex.mapState(['config', 'plot_type', 'chart_data', 'analysis']),
      is_type_plotly() {
        return (
          this.plot_type == 'bar' ||
          this.plot_type == 'scatter' ||
          this.plot_type == 'line' ||
          this.plot_type == 'violin' ||
          this.plot_type == 'tsne_dynamic'
        );
      },
      is_type_svg() {
        return this.plot_type === 'svg';
      },
      is_type_tsne() {
        return this.plot_type === 'tsne_static';
      },
      is_there_data_to_draw() {
        return (
          Object.entries(this.chart_data).length !== 0 &&
          this.chart_data.constructor === Object
        );
      },
      gene_selected() {
        return this.config.gene_symbol;
      }
    },
  });

  const datasetDisplay = Vue.component('dataset-display', {
    template: `
      <b-container v-if='!loading' fluid>
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
              v-else-if='is_type_svg && is_there_data_to_draw'
              :chart_data='chart_data'
              class="sticky-chart"
              :low_color="config.colors.low_color"
              :high_color="config.colors.high_color"/>
            <tsne-chart
              v-else-if='is_type_tsne'
              :analysis_id='config.analysis.id'
              :gene_symbol='config.gene_symbol'
            ></tsne-chart>
          </b-col>
        </b-row>
      </b-container>
    `,
    props: ['display_id'],
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
      ...Vuex.mapState(['plot_type', 'chart_data', 'config']),
      ...Vuex.mapGetters(['user_display']),
      is_creating_new_display() {
        return this.display_id === 'new';
      },
      is_type_plotly() {
        return (
          this.plot_type === 'bar' ||
          this.plot_type === 'scatter' ||
          this.plot_type === 'line' ||
          this.plot_type === 'violin' ||
          this.plot_type === 'tsne_dynamic'
        );
      },
      is_type_svg() {
        return this.plot_type === 'svg';
      },
      is_type_tsne() {
        return this.plot_type === 'tsne_static';
      },
      is_there_data_to_draw() {
        return (
          Object.entries(this.chart_data).length !== 0 &&
          this.chart_data.constructor === Object
        );
      },
    },
    created() {
      const display_data = this.user_display(this.display_id);
      this.set_display_data(display_data);

      // this.loading = true;

      // if (this.display_id !== 'new') {
      //   // if user is not creating a new display
      //   // fetch data for the display they are editing
      //   const { data } = await this.fetch_display();
      //   this.label = data.label;
      //   this.plot_type = data.plot_type;
      //   this.config = JSON.parse(data.plotly_config);
      // }
      // this.loading = false;
    },
    methods: {
      ...Vuex.mapActions(['set_display_data']),
      // initial_label() {
      //   const today = new Date();
      //   const date = `${today.getFullYear()}-${today.getMonth()+1}-${today.getDate()}`;
      //   const time = `${today.getHours()}:${today.getMinutes()}:${today.getSeconds()}`;
      //   return `${date} ${time}`;
      // },
    },
  });

  const store = new Vuex.Store({
    state: {
      user: null,
      display_id: null,
      user_displays: [],
      owner_displays: [],
      default_display_id: 0,
      dataset_id: '',
      owner_id: null,
      config: {
        gene_symbol: '',
        analysis: null,
        colors: null,
        // Other properties will be set reactive (must add via Vue.set)
        // depending on the chart type
      },
      gene_symbols: [],
      dataset_type: 'primary',
      // why is analysis here too and within config?
      analysis: null,
      columns: [],
      levels: {},
      chart_data: {},
      is_public: false,
      label: '',
      title: '',
      plot_type: null,
      loading_chart: false,
      available_plot_types: {},
      image_data: null,
      tsne_is_loading: false
    },
    getters: {
      is_user_owner(state) {
        return state.user.id === state.owner_id;
      },
      user_display(state) {
        return display_id =>
          state.user_displays.find(display => display.id == display_id);
      },
      owner_display(state) {
        return display_id =>
          state.owner_displays.find(display => display.id == display_id);
      },
      user_displays(state) {
        return state.user_displays.map(display => {
          return {
            is_default: state.default_display_id == display.id ? true : false,
            ...display,
          };
        });
      },
      owner_displays(state) {
        return state.owner_displays.map(display => {
          return {
            is_default: state.default_display_id == display.id ? true : false,
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
          gene_symbol: '',
          analysis: state.config.analysis,
          colors: null
        };

        if (
          plot_type === 'bar' ||
          plot_type === 'line' ||
          plot_type === 'violin' ||
          plot_type === 'scatter' ||
          plot_type === 'tsne_dynamic'
        ) {
          Vue.set(state.config, 'x_axis', null);
          Vue.set(state.config, 'y_axis', null);
          Vue.set(state.config, 'point_label', null);
          Vue.set(state.config, 'hide_x_labels', false);
          Vue.set(state.config, 'hide_y_labels', false);
          Vue.set(state.config, 'color_name', null);
          Vue.set(state.config, 'facet_row', null);
          Vue.set(state.config, 'facet_col', null);
          Vue.set(state.config, 'marker_size', null);
          Vue.set(state.config, 'jitter', null);
          Vue.set(state.config, 'colors', {});
          Vue.set(state.config, 'order', {});
        } else if (plot_type === 'svg') {
          Vue.set(state.config, 'colors', {
            // arbituary default colors (purple)
            low_color: '#e7d1d5',
            high_color: '#401362',
          });
        } else if (plot_type === 'tsne_static') {
          // tsne
          Vue.set(state.config, 'x_axis', null);
          Vue.set(state.config, 'y_axis', null);
          Vue.set(state.config, 'colors', {});
          Vue.set(state.config, 'colorize_legend_by', null);
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
      set_columns(state, columns) {
        state.columns = [...columns];
      },
      set_levels(state, levels) {
        state.levels = { ...levels };
      },
      reset_config(state) {
        state.config = { gene_symbol: '' };
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
      set_point_label(state,point_label) {
        state.config.point_label = point_label;
      },
      set_hide_x_labels(state, hide) {
        state.config.hide_x_labels = hide;
      },
      set_hide_y_labels(state, hide) {
        state.config.hide_y_labels = hide;
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
        if (typeof colors !== 'string') state.config.colors = { ...colors };
      },
      set_color(state, { name, color }) {
        if (state.plot_type !== 'svg') {
          // update plotly chart data
          const { data } = state.chart_data.plot_json;
          data
            .filter(el => el.name === name)
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
          gene_symbol: '',
          colors: null
        };
        state.plot_type = null;
        state.available_plot_types = {};

        state.analysis = analysis;
        state.config.analysis = analysis;
      },
      update_display(state, display) {
        const { id } = display;
        const i = state.user_displays.findIndex(el => el.id == id);
        Vue.set(state.user_displays, i, display);
      },
      delete_display(state, { display_id }) {
        state.user_displays = [
          ...state.user_displays.filter(display => display.id != display_id),
        ];
      },
      set_loading_chart(state, is_chart_loading) {
        state.loading_chart = is_chart_loading;
      },
      set_colorize_legend_by(state, legend_by) {
          state.config.colorize_legend_by = legend_by;
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
      set_tsne_is_loading(state, is_loading) {
        state.tsne_is_loading = is_loading;
      }
    },

    actions: {
      set_dataset_type({ commit }, dataset_type) {
        commit('set_dataset_type', dataset_type);
        // When display type is changed, we want to
        // reset these other settings so dependent
        // components
        commit('set_plot_type', null);
        commit('set_chart_data', {});
        commit('set_gene_symbol', '');
        commit('set_analysis', null);
      },
      update_display_label({ commit }, display_label) {
        commit('set_display_label', display_label);
      },
      set_dataset_id({ commit }, dataset_id) {
        commit('set_dataset_id', dataset_id);
      },
      set_owner_id({ commit }, owner_id) {
        commit('set_owner_id', owner_id);
      },
      set_is_public({ commit }, is_public) {
        commit('set_dataset_id', is_public);
      },
      set_title({ commit }, title) {
        commit('set_title', title);
      },
      set_tsne_is_loading({ commit }, is_loading) {
        commit('set_tsne_is_loading', is_loading);
      },
      set_dataset_info({ commit }, payload) {
        const { dataset_id, title, owner_id, is_public } = payload;
        commit('set_dataset_id', dataset_id);
        commit('set_owner_id', owner_id);
        commit('set_is_public', is_public);
        commit('set_title', title);
      },
        set_plot_type({ commit }, plot_type) {
        commit('set_plot_type', plot_type);
        // When display type is changed, we want to
        // reset these other settings so dependent
        // components
        commit('set_chart_data', {});
        commit('set_gene_symbol', '');
      },
      async fetch_dataset_info({ commit }, dataset_id) {
        commit('set_dataset_id', dataset_id);
        const { title, is_public, owner_id } = await $.ajax({
          url: './cgi/get_dataset_info.cgi',
          type: 'POST',
          data: { dataset_id },
          dataType: 'json',
        });
        commit('set_owner_id', owner_id);
        commit('set_is_public', is_public);
        commit('set_title', title);
      },
      async fetch_user_displays({ commit }, { user_id, dataset_id }) {
        let displays = await $.ajax({
          url: './cgi/get_dataset_displays.cgi',
          type: 'POST',
          data: { user_id, dataset_id },
          dataType: 'json',
        });
        commit('set_user_displays', displays);
      },
      async fetch_owner_displays({ commit }, { owner_id, dataset_id }) {
        const displays = await $.ajax({
          url: './cgi/get_dataset_displays.cgi',
          type: 'POST',
          data: { user_id: owner_id, dataset_id },
          dataType: 'json',
        });

        commit('set_owner_displays', displays);
      },
      async fetch_available_plot_types(
        { commit, state },
        { user_id, session_id, dataset_id }
      ) {
        const { config } = state;
        const analysis_id = config.analysis ? config.analysis.id : null;
        const { data: available_plot_types } = await axios.post(
          `/api/h5ad/${dataset_id}/availableDisplayTypes`,
          {
            user_id,
            session_id,
            dataset_id,
            analysis_id,
          }
        );
        commit('set_available_plot_types', available_plot_types);
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
        commit('set_columns', obs_columns);
        commit('set_levels', obs_levels);
      },
      async fetch_plotly_data({ commit }, { config, plot_type, dataset_id }) {
        commit('set_loading_chart', true);
        const payload = { ...config, plot_type };
        const { gene_symbol } = config;

        const { data } = await axios.post(`/api/plot/${dataset_id}`, payload);
        commit('set_chart_data', data);

        const {
          plot_colors,
          plot_order,
          x_axis,
          y_axis,
          point_label,
          hide_x_labels,
          hide_y_labels,
          color_name,
          facet_row,
          facet_col,
          marker_size,
          jitter,
        } = data;

        commit('set_order', plot_order);
        commit('set_colors', plot_colors);
        commit('set_x_axis', x_axis);
        commit('set_y_axis', y_axis);
        commit('set_point_label', point_label);
        commit('set_hide_x_labels', hide_x_labels);
        commit('set_hide_y_labels', hide_y_labels);
        commit('set_color_name', color_name);
        commit('set_facet_row', facet_row);
        commit('set_facet_col', facet_col);
        commit('set_marker_size', marker_size);
        commit('set_jitter', jitter);

        commit('set_loading_chart', false);
      },
      async fetch_tsne_image({ commit }, {config, dataset_id}) {
        commit('set_tsne_is_loading', true);

        const {data} = await axios.get(`/api/plot/${dataset_id}/tsne`, {
          params: {
              gene: config.gene_symbol,
              analysis: this.analysis_id,
              colorize_by: config.colorize_legend_by,
              x_axis: config.x_axis,
              y_axis: config.y_axis,
              analysis_owner_id: this.user_id,
              colors: config.colors,
              // helps stop caching issues
              timestamp: new Date().getTime()
          }
        });

        commit('set_image_data', data);
        commit('set_tsne_is_loading', false);
      },
      set_index({ commit }, index) {
        commit('set_index', index);
      },
      set_order({ commit }, order) {
        commit('set_order', order);
        commit('set_levels', order);
      },
      set_color({ commit }, { name, color }) {
        commit('set_color', { name, color });
      },
      set_gene_symbol({ commit }, gene_symbol) {
        commit('set_gene_symbol', gene_symbol);
      },
      async fetch_gene_symbols({ commit }, { dataset_id, analysis_id }) {
        const base = `./api/h5ad/${dataset_id}/genes`;
        const query = analysis_id ? `?analysis=${analysis_id}` : '';

        const { data } = await axios.get(`${base}${query}`);
        commit('set_gene_symbols', data.gene_symbols);
      },
      set_label({ commit }, label) {
        commit('set_label', label);
      },
      set_display_data({ commit }, display) {
        let { label, plot_type, plotly_config: config } = display;

        commit('set_label', label);
        commit('set_plot_type', plot_type);
        commit('set_config', config);
      },
      reset({ commit }) {
        commit('set_label', '');
        commit('set_plot_type', null);
        commit('reset_config');
      },
      async fetch_svg_data({ commit }, { gene_symbol, dataset_id }) {
        const { data } = await axios.get(
          `/api/plot/${dataset_id}/svg?gene=${gene_symbol}`
        );
        commit('set_chart_data', data);
      },
      async fetch_default_display({ commit }, { user_id, dataset_id }) {
        const { default_display_id } = await $.ajax({
          url: './cgi/get_default_display.cgi',
          type: 'POST',
          data: { user_id, dataset_id },
          dataType: 'json',
        });
        commit('set_default_display_id', default_display_id);
      },
      set_analysis_id({ commit }, analysis) {
        commit('set_analysis_id', analysis);
      },
      set_analysis({ commit }, analysis) {
        commit('set_analysis', analysis);
      },
      update_display({ commit }, display) {
        commit('update_display', display);
      },
      set_image_data({ commit }, image_data) {
        commit('set_image_data', image_data)
      },
      remove_display({ commit }, display_id) {
        commit('delete_display', display_id);
      },
      update_default_display_id({ commit }, { display_id }) {
        commit('set_default_display_id', display_id);
      },
    },
  });

  const routes = [
    {
      path: '/dataset/:dataset_id/displays/',
      component: datasetCurator,
      props: true,
      children: [
        {
          name: 'dashboard',
          path: '/',
          component: datasetDisplays,
        },
        {
          path: 'new',
          name: 'new',
          component: newDisplay,
          beforeEnter(to, from, next) {
            // We want to reset our data that may have been loaded
            // from a previous display
            store.dispatch('reset');
            next();
          },
        },
        {
          path: ':display_id/edit',
          component: datasetDisplay,
          props: true,
        },
      ],
    },
  ];

  const router = new VueRouter({
    routes,
  });

  const app = new Vue({
    el: '#app',
    router,
    store,
    computed: {
      ...Vuex.mapState(['user']),
    },
    created() {
      // check if the user is already logged in
      // SADKINS - 7/10/2020 - this was originally a Vue component but I removed that due to redundancy
      check_for_login();

      // We want to check for session when the curator app is first created
      session_id = Cookies.get('gear_session_id');
      sleep(500).then(() => {
        // If CURRENT_USER is defined at this point, add information as placeholder test
        if (CURRENT_USER) {
          this.$store.commit('set_user', CURRENT_USER);
        }
      })
    },
  });
})();
