/*
SAdkins note - The styling on this page differs from other gEAR JS pages
I was trying to follow the recommended Javascript style in using camelCase
in addition using the StandardJS style linter ("semistandrd" variant)
with some refactoring done with the P42+ VSCode extension (to modernize code)
However, code inherited from common.js is still in snake_case rather than camelCase
*/

'use strict';
/* global $, axios, Plotly, CURRENT_USER, session_id, check_for_login */

let numObs = 0; // dummy value to initialize with

let obsFilters = {};
let sortCategories = {"facet": null, "axis": null, "groupby": null}; // Control sorting order of subplots and axis categories
let genesFilters = [];

let plotConfig = {};  // Plot config that is passed to API or stored in DB

let selectedGenes = null;  // Genes selected in a plot (e.g. volcano with lasso tool)

let datasetId = null;
let defaultDisplayId = null;
let displayId = null;
let obsLevels = null;
let geneSymbols = null;

let datasetTree = new DatasetTree({treeDiv: '#dataset_tree'});
let geneCartTree = new GeneCartTree({treeDiv: '#gene_cart_tree'});

const plotTypes = ['dotplot', 'heatmap', 'mg_violin', 'quadrant', 'volcano'];

const dotplotOptsIds = ["#obs_groupby_container"]; //, "#obs_facet_container",];
const heatmapOptsIds = ["#heatmap_options_container", "#adv_heatmap_opts"];//, "#obs_facet_container",];
const quadrantOptsIds = ["#quadrant_options_container", "#de_test_container", "#include_zero_foldchange_container"];
const violinOptsIds = ["#obs_groupby_container", "#adv_violin_opts"]; //, "#obs_facet_container",];
const volcanoOptsIds = ["#volcano_options_container", "#de_test_container", "#adjusted_pvals_checkbox_container", "#annot_nonsig_checkbox_container"];

// color palettes
const continuous_palettes = [
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
];
const discrete_palettes = ["alphabet", "vivid", "light24", "dark24"];

// Async to ensure data is fetched before proceeding.
// This self-invoking function loads the initial state of the page.
(async () => {
  // Hide further configs until a dataset is chosen.
  // Changing the dataset will start triggering these to show
  $('#plot_type_container').hide();
  $('#advanced_options_container').hide();
  $('#gene_container').hide();

  // check if the user is already logged in
  await check_for_login();

  $('#pre_dataset_spinner').show();
  // Load gene carts and datasets before the dropdown appears
  await reloadTrees ();
  $('#pre_dataset_spinner').hide();

  // Initialize plot types
   $('#plot_type_select').select2({
    placeholder: 'Choose how to plot',
    width: '25%'
  });

  // If brought here by the "gene search results" page, curate on the dataset ID that referred us
  const linkedDatasetId = getUrlParameter("dataset_id");
  if (linkedDatasetId) {
    $('#dataset').val(linkedDatasetId);
    $('#dataset').text(datasetTree.treeData.find(e => e.dataset_id === linkedDatasetId).text);
    $('#dataset').trigger('change');
  }

  // Create observer to watch if user changes (ie. successful login does not refresh page)
  // See: https://developer.mozilla.org/en-US/docs/Web/API/MutationObserver

  // Select the node that will be observed for mutations
  const targetNode = document.getElementById('loggedin_controls');
  // Create an observer instance linked to the callback function
  const observer = new MutationObserver(reloadTrees);
  // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
  // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
  observer.observe(targetNode, { attributes: true });
})();

// Call API to return plot JSON data
async function getData (datasetId, payload) {
  try {
    return await axios.post(`/api/plot/${datasetId}/mg_dash`, {
      ...payload
    })
  } catch (e) {

    const message = "There was an error in making this plot. Please contact the gEAR team using the 'Contact' button at the top of the page.";
    const success = -1;
    const data = {message, success};
    return {data};
  }
}

// Call API to return a list of the dataset's gene symbols
async function fetchGeneSymbols (payload) {
  const { datasetId, analysis } = payload;
  const base = `./api/h5ad/${datasetId}/genes`;
  const query = analysis ? `?analysis=${analysis.id}` : '';

  const { data } = await axios.get(`${base}${query}`);
  return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
}

// Call API to return all observations
async function fetchH5adInfo (payload) {
  const { datasetId, analysis } = payload;
  const base = `./api/h5ad/${datasetId}`;
  const query = analysis ? `?analysis=${analysis.id}` : '';
  const { data } = await axios.get(`${base}${query}`);
  return data;
}

async function loadDatasets () {
  await $.ajax({
    type: 'POST',
    url: './cgi/get_h5ad_dataset_list.cgi',
    data: {
      session_id
    },
    dataType: 'json',
    success(data) {
      let counter = 0;

      // Populate select box with dataset information owned by the user
      let userDatasets = [];
      if (data.user.datasets.length > 0) {
        // User has some profiles
        $.each(data.user.datasets, (_i, item) => {
          userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
        });
      }
      // Next, add datasets shared with the user
      let sharedDatasets = [];
      if (data.shared_with_user.datasets.length > 0) {
        $.each(data.shared_with_user.datasets, (_i, item) => {
          sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
        });
      }
      // Now, add public datasets
      let domainDatasets = [];
      if (data.public.datasets.length > 0) {
        $.each(data.public.datasets, (_i, item) => {
          domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
        });
      }

      datasetTree.userDatasets = userDatasets;
      datasetTree.sharedDatasets = sharedDatasets;
      datasetTree.domainDatasets = domainDatasets;
      datasetTree.generateTree();
    },
    error(_xhr, _status, msg) {
      console.error(`Failed to load dataset list because msg: ${msg}`);
    }
  });
}

// Draw plotly chart to image
async function drawPreviewImage (display) {
  // check if config has been stringified
  let config;
  config = typeof display.plotly_config === 'string' ? JSON.parse(display.plotly_config) : display.plotly_config;

  const { data } = await getData(datasetId, config);

  // If there was an error in the plot, put alert up
  if ( data.success < 1 ) {
    $(`#modal-display-img-${display.id} + .js-plot-error`).show();
    $(`#modal-display-img-${display.id} + .js-plot-error`).html(data.message);
    $(`#modal-display-${display.id}-loading`).hide();
    return;
  }

  const { plot_json: plotlyJson, plot_config: plotlyConfig } = data;
  Plotly.toImage(
    { ...plotlyJson, plotlyConfig },
    { height: 500, width: 500 }
  ).then(url => {
    $(`#modal-display-img-${display.id}`).prop('src', url);
  }).then(
    () => { $(`#modal-display-${display.id}-loading`).hide(); }
  );
}

// Draw plotly chart in HTML
function drawChart (data, datasetId) {
  const targetDiv = `dataset_${datasetId}_h5ad`;
  const parentDiv = `dataset_${datasetId}`;
  const { plot_json: plotlyJson, plot_config: plotlyConfig, message, success } = data;

  // Since default plots are now added after dataset selection, wipe the plot when a new one needs to be drawn
  $(`#${targetDiv}`).empty()

  // If there was an error in the plot, put alert up
  if ( success < 1 || !plotlyJson.layout) {
    $(`#${parentDiv} .js-plot-error`).show();
    $(`#${parentDiv} .js-plot-error`).html(message);
    return;
  }

  const layoutMods = {
    //height: targetDiv.clientHeight,
    //width: targetDiv.clientWidth,
  };

  // NOTE: This will definitely affect the layout on the gene search results page
  // if the "height" style for the container in CSS is removed.
  if ($('#plot_type_select').select2('data')[0].id === 'heatmap') {
    if (genesFilters.length > 50) {
      layoutMods.height = genesFilters.length * 10;
    }
  } else if (['quadrant', 'volcano'].includes($('#plot_type_select').select2('data')[0].id)) {
    layoutMods.height = 800;
    //layoutMods.width = 1080; // If window is not wide enough, the plotly option icons will overlap contents on the right
  } else if ($('#plot_type_select').select2('data')[0].id === "mg_violin" && $("#stacked_violin").is(":checked")){
    layoutMods.height = 800;
  }

  // Overwrite plot layout and config values with custom ones from display
  const layout = {
    ...plotlyJson.layout,
    ...layoutMods
  };

  const configMods = {
    responsive: true
  };

  const config = {
    ...plotlyConfig,
    ...configMods
  };
  Plotly.newPlot(targetDiv, plotlyJson.data, layout, config);

  // Show any warnings from the API call
  if (message && success > 1) {
    $(`#${parentDiv} .js-plot-warning`).show();
    $(`#${parentDiv} .js-plot-warning`).html(`<ul>${message}</ul>`);
  }

  // If plot data is selected, create the right-column table and do other misc things
  $(`#dataset_${datasetId}_h5ad`).on("plotly_selected", (_e, data) => {

    if (!(['volcano', 'quadrant'].includes(plotConfig.plot_type))) {
      return;
    }

    $("#selected_genes_btn").prop("disabled", false);

    // Note: the jQuery implementation has slightly different arguments than what is in the plotlyJS implementation
    // We want 'data', which returns the eventData PlotlyJS events normally return
    selectedGenes = [];

    data.points.forEach((pt) => {
      selectedGenes.push({
        gene_symbol: pt.data.text[pt.pointNumber],
      });
    });

    // Sort in alphabetical order
    selectedGenes.sort();

    const template = $.templates("#selected_genes_tmpl");
    const htmlOutput = template.render(selectedGenes);
    $("#selected_genes_list").html(htmlOutput);

  });
}

// Submit API request and draw the HTML
async function draw (datasetId, payload) {
  const {data } = await getData(datasetId, payload);
  drawChart(data, datasetId);
}

// If user changes, update genecart/profile trees
async function reloadTrees(){

  // Update dataset and genecart trees in parallel
  // Works if they were not populated or previously populated
  await Promise.all([loadDatasets(), loadGeneCarts()]);
}

// Render the gene-selection dropdown menu
function createGeneDropdown (genes) {
  const tmpl = $.templates('#gene_dropdown_tmpl');
  const html = tmpl.render({ genes });
  $('#gene_dropdown_container').html(html);
  $('#gene_dropdown').select2({
    placeholder: 'To search, click to select or start typing some gene names',
    allowClear: true,
    width: 'resolve'
  });
  $('#gene_spinner').hide();
}

// Render the observation groupby field HTML
function createObsGroupbyField (obsLevels) {
  const tmpl = $.templates('#obs_groupby_tmpl'); // Get compiled template using jQuery selector for the script block
  const html = tmpl.render(obsLevels); // Render template using data - as HTML string
  $('#obs_groupby_container').html(html); // Insert HTML string into DOM
}

// Render the observation filter dropdowns
function createObsFilterDropdowns (obsLevels) {
  const tmpl = $.templates('#obs_filters_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_filters_container').html(html);
  $('select.js-obs-levels').select2({
    placeholder: 'Start typing to include groups from this category. Click "All" to use all groups',
    allowClear: true,
    width: 'resolve'
  });
}

// Render the observation facet field
function createObsFacetField (obsLevels) {
  const tmpl = $.templates('#obs_facet_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_facet_container').html(html);
}

// Render the observation axis sort field
function createObsSortField (obsLevels) {
  const tmpl = $.templates('#obs_axis_sort_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_axis_sort_container').html(html);
}

// Render the observation axis sort field
function createObsClusterBarField (obsLevels) {
  const tmpl = $.templates('#obs_clusterbar_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_clusterbar_container').html(html);
}

// Render the sortable list for the chosen category
function createObsSortable (obsLevel, scope) {
  const propData = $(`#${obsLevel}_dropdown`).select2('data');
  const sortData = propData.map((elem) => elem.id);
  const tmpl = $.templates('#obs_sortable_tmpl');
  const html = tmpl.render({ sortData });
  $(`#${scope}_sortable`).html(html);
  $(`#${scope}_sortable`).sortable();

  // Store category name to retrieve later.
  sortCategories[scope] = obsLevel;

  $(`#${scope}_order_label`).show();
}

// Render dropdowns specific to the dot plot
function createDotplotDropdowns (obsLevels) {
  createObsGroupbyField (obsLevels);
  //createObsFacetField//(obsLevels);
}

// Render dropdowns specific to the heatmap plot
function createHeatmapDropdowns (obsLevels) {
  //createObsFacetField(obsLevels);
  createObsSortField(obsLevels);
  createObsClusterBarField(obsLevels);

  // Initialize differential expression test dropdown
  $('#distance_select').select2({
    placeholder: 'Choose distance metric',
    width: '25%'
  });

  $('#cluster_obs_warning').hide();
}

// Render dropdowns specific to the quadrant plot
function createQuadrantDropdowns (obsLevels) {
  const tmpl = $.templates('#select_conditions_tmpl');
  const html = tmpl.render(obsLevels);
  $('#quadrant_compare1_condition').html(html);
  $('#quadrant_compare1_condition').select2({
    placeholder: 'Select the first query condition.',
    width: '25%'
  });
  $('#quadrant_compare2_condition').html(html);
  $('#quadrant_compare2_condition').select2({
    placeholder: 'Select the second query condition.',
    width: '25%'
  });
  $('#quadrant_ref_condition').html(html);
  $('#quadrant_ref_condition').select2({
    placeholder: 'Select the reference condition.',
    width: '25%'
  });

  // Initialize differential expression test dropdown
  $('#de_test_select').select2({
    placeholder: 'Choose DE testing algorithm',
    width: '25%'
  });
}

// Render dropdowns specific to the violin plot
function createViolinDropdowns (obsLevels) {
  createObsGroupbyField (obsLevels);
  //createObsFacetField(obsLevels);
}

// Render dropdowns specific to the volcano plot
function createVolcanoDropdowns (obsLevels) {
  const tmpl = $.templates('#select_conditions_tmpl');
  const html = tmpl.render(obsLevels);
  $('#volcano_query_condition').html(html);
  $('#volcano_query_condition').select2({
    placeholder: 'Select the query condition.',
    width: '25%'
  });
  $('#volcano_ref_condition').html(html);
  $('#volcano_ref_condition').select2({
    placeholder: 'Select the reference condition.',
    width: '25%'
  });

  // Initialize differential expression test dropdown
  $('#de_test_select').select2({
    placeholder: 'Choose DE testing algorithm',
    width: '25%'
  });

}

function curateObservations (obsLevels) {
  // Delete useless filters
  for (const property in obsLevels) {
    if (property === 'color' || property.endsWith('colors')) {
      delete obsLevels[property];
    } else if (obsLevels[property].length === 1) {
      delete obsLevels[property];
    }
  }
  return obsLevels;
}

// Generate a list of saved plot displays the user has access in viewing.
async function loadSavedDisplays (datasetId, defaultDisplayId=null) {
  const datasetData = await fetchDatasetInfo(datasetId);
  const { owner_id: ownerId } = datasetData;
  const userDisplays = await fetchUserDisplays(CURRENT_USER.id, datasetId);
  // Do not duplicate user displays in the owner display area as it can cause HTML element issues
  const ownerDisplays = CURRENT_USER.id === ownerId ? [] : await fetchOwnerDisplays(ownerId, datasetId);

  // Filter displays to those only with multigene plot types
  const mgUserDisplays = userDisplays.filter(d => plotTypes.includes(d.plot_type));
  const mgOwnerDisplays = ownerDisplays.filter(d => plotTypes.includes(d.plot_type));

  //
  mgUserDisplays.forEach(display => {
    display.is_default = false;
    if (defaultDisplayId && display.id === defaultDisplayId ) {
      display.is_default = true;
    }
    drawPreviewImage(display);
  });
  mgOwnerDisplays.forEach(display => {
    display.is_default = false;
    if (defaultDisplayId && display.id === defaultDisplayId ) {
      display.is_default = true;
    }
    drawPreviewImage(display);
  });

  const displaysTmpl = $.templates('#saved_display_modal_tmpl');
  const displaysHtml = displaysTmpl.render({
    dataset_id: datasetId,
    user_displays: mgUserDisplays,
    owner_displays: mgOwnerDisplays
  });

  $('#saved_display_modal').html(displaysHtml);
}

// Populate the HTML config options based on what was in the plot
function loadDisplayConfigHtml (plotConfig) {
  // NOTE: The calling function also clicks "#reset_opts", so the options are rendered already
  // Populate filter-by dropdowns
  obsFilters = plotConfig.obs_filters;
  for (const property in obsFilters) {
    $(`#${property}_dropdown`).val(obsFilters[property]);
    $(`#${property}_dropdown`).trigger('change');
  }

  // Populate plot type-specific dropdowns and checkbox options
  switch ($('#plot_type_select').val()) {
    case 'dotplot':
      $(`#${plotConfig.groupby_filter}_groupby`).prop('checked', true).click();
      $(`#${plotConfig.facet_col}_facet`).prop('checked', true).click();
      break;
    case 'heatmap':
      $(`#${plotConfig.groupby_filter}_groupby`).prop('checked', true).click();
      $(`#${plotConfig.axis_sort_col}_sort`).prop('checked', true).click();
      $(`#${plotConfig.facet_col}_facet`).prop('checked', true).click();
      for (const field in plotConfig.clusterbar_fields) {
        $(`#${field}_clusterbar`).prop('checked', true).click();
      }
      $('#matrixplot').prop('checked', plotConfig.matrixplot);
      $('#cluster_obs').prop('checked', plotConfig.cluster_obs);
      $('#cluster_genes').prop('checked', plotConfig.cluster_obs);
      $('#flip_axes').prop('checked', plotConfig.flip_axes);
      $('#distance_select').val(plotConfig.distance_metric);
      $('#distance_select').trigger('change');
      break;
    case 'mg_violin':
      $(`#${plotConfig.groupby_filter}_groupby`).prop('checked', true).click();
      $(`#${plotConfig.facet_col}_facet`).prop('checked', true).click();
      $('#stacked_violin').prop('checked', plotConfig.stacked_violin);
      $('#violin_add_points').prop('checked', plotConfig.violin_add_points);
      break;
    case 'quadrant':
      $('#include_zero_foldchange').prop('checked', plotConfig.include_zero_fc);
      $("#quadrant_foldchange_cutoff").val(plotConfig.fold_change_cutoff);
      $("#quadrant_fdr_cutoff").val(plotConfig.fdr_cutoff);
      $('#quadrant_compare1_condition').val(plotConfig.compare1_condition);
      $('#quadrant_compare1_condition').trigger('change');
      $('#quadrant_compare2_condition').val(plotConfig.compare2_condition);
      $('#quadrant_compare2_condition').trigger('change');
      $('#quadrant_ref_condition').val(plotConfig.ref_condition);
      $('#quadrant_ref_condition').trigger('change');

      $('#de_test_select').val(plotConfig.de_test_algo);
      $('#de_test_select').trigger('change');
      break;
    default:
      // volcano
      $('#adj_pvals').prop('checked', plotConfig.adj_pvals);
      $('#annot_nonsig').prop('checked', plotConfig.annotate_nonsignificant)
      $('#volcano_query_condition').val(plotConfig.query_condition);
      $('#volcano_query_condition').trigger('change');
      $('#volcano_ref_condition').val(plotConfig.ref_condition);
      $('#volcano_ref_condition').trigger('change');
      $('#de_test_select').val(plotConfig.de_test_algo);
      $('#de_test_select').trigger('change');
  }
}

// Load all saved gene carts for the current user
function loadGeneCarts () {
  const d = new $.Deferred();

  if (!session_id) {
    // User is not logged in. Hide gene carts container
    $('#gene_cart_container').hide();
    d.resolve();
  } else {
    $.ajax({
      url: './cgi/get_user_gene_carts.cgi',
      type: 'post',
      data: { session_id },
      dataType: 'json',
      success(data, _textStatus, _jqXHR) { // source https://stackoverflow.com/a/20915207/2900840
        const userGeneCarts = [];

        if (data.gene_carts.length > 0) {
          // User has some profiles
          $.each(data.gene_carts, (_i, item) => {
            userGeneCarts.push({ value: item.id, text: item.label });
          });

          // No domain gene carts yet
          geneCartTree.userGeneCarts = userGeneCarts;
          geneCartTree.generateTree();

          $('#gene_cart_container').show();

        } else {
          $('#gene_cart_container').hide();
        }
        d.resolve();
      },
      error(_jqXHR, _errorThrown) {
        // display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        d.fail();
      }
    });
  }
  d.promise();
}

function saveGeneCart () {
  // must have access to USER_SESSION_ID
  const gc = new GeneCart({
    session_id: CURRENT_USER.session_id
    , label: $("#gene_cart_name").val()
    , gctype: "unweighted-list"
    , organism_id:  $("#dataset").data('organism-id')
    , is_public: 0
  });

  selectedGenes.forEach((pt) => {
    const gene = new Gene({
      id: pt.gene_id,
      gene_symbol: pt.gene_symbol,
    });
    gc.add_gene(gene);
  });

  gc.save(updateUIAfterGeneCartSaveSuccess, updateUIAfterGeneCartSaveFailure);
}

function updateUIAfterGeneCartSaveSuccess(gc) {
  $('#saved_gene_cart_confirmation').text('Gene cart successfully saved!');
  $('#saved_gene_cart_confirmation').addClass('text-success');
  $('#saved_gene_cart_confirmation').show();
}

function updateUIAfterGeneCartSaveFailure(gc) {
  $('#saved_gene_cart_confirmation').text('Issue with saving gene cart.');
  $('#saved_gene_cart_confirmation').addClass('text-danger');
  $('#saved_gene_cart_confirmation').show();
}

function fetchDatasetInfo (datasetId) {
  return $.ajax({
    url: './cgi/get_dataset_info.cgi',
    type: 'POST',
    data: { dataset_id: datasetId },
    dataType: 'json'
  });
}

function fetchUserDisplays (userId, datasetId) {
  return $.ajax({
    url: './cgi/get_dataset_displays.cgi',
    type: 'POST',
    data: { user_id: userId, dataset_id: datasetId },
    dataType: 'json'
  });
}
function fetchOwnerDisplays (ownerId, datasetId) {
  return $.ajax({
    url: './cgi/get_dataset_displays.cgi',
    type: 'POST',
    data: { user_id: ownerId, dataset_id: datasetId },
    dataType: 'json'
  });
}

function getDefaultDisplay (datasetId) {
  return $.ajax({
    url: './cgi/get_default_display.cgi',
    type: 'POST',
    data: {
      user_id: CURRENT_USER.id,
      dataset_id: datasetId,
      is_multigene: 1
    },
    dataType: 'json'
  });
}

$('#dataset').change(async function () {
  datasetId = $('#dataset').val();
  displayId = null;

  $('#post_dataset_spinner').show();

  // Obtain default display ID for this dataset
  const {default_display_id: defaultDisplayId} = await getDefaultDisplay(datasetId);

  // Populate saved displays modal
  loadSavedDisplays(datasetId, defaultDisplayId);

  $('#load_saved_plots').show();
  $('#plot_type_container').show();

  $('#gene_container').show();
  $('#gene_spinner').show();
  $('#genes_not_found').hide();
  $('#too_many_genes_warning').hide();
  $('#gene_cart').text('Choose gene cart'); // Reset gene cart text

  // Create promises to get genes and observations for this dataset
  const geneSymbolsPromise = fetchGeneSymbols({ datasetId, undefined })
  const h5adPromise =  fetchH5adInfo({ datasetId, undefined });

  // Execute both in parallel
  geneSymbols = await geneSymbolsPromise;
  const data = await h5adPromise;

  createGeneDropdown(geneSymbols);  // gene_spinner hidden here

  $('#post_dataset_spinner').hide();

  // Cannot cluster with just one gene (because function is only available
  // in dash.clustergram which requires 2 or more genes in plot)
  // Adding a gene will trigger a change to enable the property
  $("#cluster_obs").prop("disabled", true);
  $("#cluster_genes").prop("disabled", true);

  // Get categorical observations for this dataset
  obsLevels = curateObservations(data.obs_levels);
  numObs = data.num_obs;

  $('#create_plot').show();
  $("#create_plot").prop("disabled", true);
  $('#reset_opts').show();

  // If a plot type was already selected,
  // reset the options so configs are populated for the current dataset
  if ($('#plot_type_select').val() ) {
    $("#create_plot").prop("disabled", false);
    $('#reset_opts').click();
  }

  /* SAdkins - Commenting out for the time being.
  Am annoyed that I have to wait for the plot to load (when I don't want to use this plot)

  // Load an initial plot (just to populate the plot space)
  if (defaultDisplayId) {
    // Easiest way to do this.  Also populates the conditions
    $(`#${defaultDisplayId}_load`).click();
  } else {
    if (! obsLevels) {
      return;
    }
    // Use the first field in our categorical observations for the volcano
    // but "cluster" takes precedence
    let field = Object.keys(obsLevels)[0];
    if ("cluster" in obsLevels) {
       field = "cluster";
    }
    // Cannot make volcano if this field does not have two conditions
    // Instead of trying to find ways to make the plot, just cut our losses
    if (obsLevels[field].length < 2) {
      return;
    }
    const loadPlotConfig = {
      plot_type: 'volcano'
      , obs_filters: obsLevels // Just keep everything
      , query_condition: `${field};-;${obsLevels[field][0]}`
      , ref_condition: `${field};-;${obsLevels[field][1]}`
      , use_adj_pvals: true
    };
    plotConfig = loadPlotConfig;

    // Draw the updated chart
    $('#dataset_spinner').show();
    const plotTemplate = $.templates('#dataset_plot_tmpl');
    const plotHtml = plotTemplate.render({ dataset_id: datasetId });
    $('#dataset_plot').html(plotHtml);
    //NOTE: Height will not change since select2 element is not updated
    await draw(datasetId, loadPlotConfig);
    $('#dataset_spinner').hide();

    // Show plot options and disable selected genes button (since genes are not selected anymore)
    $('#post_plot_options').show();
    $("#selected_genes_btn").prop("disabled", true);

    // Set plot config so that
    $('#plot_type_select').val(plotConfig.plot_type);
    $('#plot_type_select').trigger('change');
  }
  */
});

$("#save_gene_cart").on("click", () => {
  $("#save_gene_cart").prop("disabled", true);

  if (CURRENT_USER) {
    saveGeneCart();
  } else {
    window.alert("You must be signed in to do that.");
  }
  $("#save_gene_cart").prop('disabled', false);
});

$("#download_plot").on("click", () => {
  Plotly.downloadImage(
    `dataset_${datasetId}_h5ad`, {
      width: $('#plot_download_width').val()
      , height: $('#plot_download_height').val()
      , scale: $('#plot_download_scale').val()
    }
  )
});

// Load user's gene carts
$('#gene_cart').change(function () {
  let geneCartId = $(this).val();
  const params = { session_id, gene_cart_id: geneCartId };
  const d = new $.Deferred(); // Causes editable to wait until results are returned
  // User is not logged in
  if (!session_id) {
    d.resolve();
  } else {
    // User is logged in

    // Get the gene cart members and populate the gene symbol search bar
    $.ajax({
      url: './cgi/get_gene_cart_members.cgi',
      type: 'post',
      data: params,
      success(data, _newValue, _oldValue) {
        if (data.success === 1) {
          // Append gene symbols to search bar
          const geneCartSymbols = [];

          // format gene symbols into search string
          $.each(data.gene_symbols, (_i, item) => {
            geneCartSymbols.push(item.label);
          });

          const geneCartSymbolsLowerCase = geneCartSymbols.map(x => x.toLowerCase());
          const geneSymbolsLowerCase = geneSymbols.map(x => x.toLowerCase());

          // Get genes from gene cart that are present in dataset's genes.  Preserve casing of dataset's genes.
          const intersection = geneSymbols.filter(x => geneCartSymbolsLowerCase.includes(x.toLowerCase()));

          // Get genes from gene cart that are not in dataset's genes
          const difference = geneCartSymbols.filter(x => !geneSymbolsLowerCase.includes(x.toLowerCase()));

          $('#gene_dropdown').val(intersection);
          $('#gene_dropdown').trigger('change');

          if (difference.length > 0) {
            const differenceString = difference.join(', ');
            $('#genes_not_found').text(`The following gene cart genes were not found in this dataset: ${differenceString}`);
            $('#genes_not_found').show();
          } else {
            $('#genes_not_found').hide();
          }
        }
        d.resolve();
      }
    });
  }
  return d.promise();
});

// Some options are specific to certain plot types
$('#plot_type_select').change(() => {
  $('#reset_opts').click();  // Reset all options
  $('#selected_genes_field').hide();
  $("#create_plot").prop("disabled", false);
  $('#advanced_options_container').show();

  dotplotOptsIds.forEach(id => {
    $(id).hide();
  })
  heatmapOptsIds.forEach(id => {
    $(id).hide();
  });
  quadrantOptsIds.forEach(id => {
    $(id).hide();
  })
  violinOptsIds.forEach(id => {
    $(id).hide();
  });
  volcanoOptsIds.forEach(id => {
    $(id).hide();
  });

  switch ($('#plot_type_select').val()) {
    case 'dotplot':
      dotplotOptsIds.forEach(id => {
        $(id).show();
      })
      $("#gene_selection_help").text("Choose the genes to include in plot.");
      $('#advanced_options_container').hide();  // no advanced opts yet for this plottype
      break;
    case 'heatmap':
      heatmapOptsIds.forEach(id => {
        $(id).show();
      });
      $("#gene_selection_help").text("Choose the genes to include in plot.");
      break;
    case 'mg_violin':
      violinOptsIds.forEach(id => {
        $(id).show();
      });
      $("#gene_selection_help").text("Choose the genes to include in plot.");
      break;
    case 'quadrant':
      quadrantOptsIds.forEach(id => {
        $(id).show();
      })
      $("#gene_selection_help").text("OPTIONAL: Gene selection is optional for this plot type. Selected genes are annotated in the plot.");
      break;
    default:
      // volcano
      volcanoOptsIds.forEach(id => {
        $(id).show();
      });
      $("#gene_selection_help").text("OPTIONAL: Gene selection is optional for this plot type. Selected genes are annotated in the plot.");
    }
});

$(document).on('change', '#cluster_obs', () => {
  $('#cluster_obs_warning').hide();
  if ($('#cluster_obs').is(':checked') && numObs >= 1000){
    $('#cluster_obs_warning').show();
  }
});

$(document).on('change', '#gene_dropdown', () => {
  const genesFilters = $('#gene_dropdown').select2('data').map((elem) => elem.id);
  // Show warning if too many genes are entered
  $("#too_many_genes_warning").hide();
  if (genesFilters.length > 10) {
    $("#too_many_genes_warning").text(`There are currently ${genesFilters.length} genes to be plotted. Be aware that with some plots, a high number of genes can make the plot congested or unreadable.`);
    $("#too_many_genes_warning").show();
  }

  // Cannot cluster columns with just one gene (because function is only available
  // in dash.clustergram which requires 2 or more genes in plot)
  if (genesFilters.length > 1) {
    $("#cluster_obs").prop("disabled", false);
    $("#cluster_genes").prop("disabled", false);
  } else {
    $("#cluster_obs").prop("disabled", true);
    $("#cluster_obs").prop("checked", false);
    $("#cluster_genes").prop("disabled", true);
    $("#cluster_genes").prop("checked", false);
  }
});

// When a column is chosen, populate the sortable list
$(document).on('change', 'input[name="obs_facet"]', () => {
  const obsLevel = $('input[name="obs_facet"]:checked').val();
  createObsSortable(obsLevel, "facet");
});

$(document).on('change', 'input[name="obs_sort"]', () => {
  const obsLevel = $('input[name="obs_sort"]:checked').val();
  createObsSortable(obsLevel, "axis");
});

// "obs_groupby" stuff controls the plot axis content like "obs_sort"
$(document).on('change', 'input[name="obs_groupby"]', () => {
  const obsLevel = $('input[name="obs_groupby"]:checked').val();
  createObsSortable(obsLevel, "groupby");
});

$(document).on('click', '#create_plot', async () => {

  // Reset plot errors and warnings for both plots
  $('.js-plot-error').empty().hide();
  $('.js-plot-warning').empty().hide();

  plotConfig = {};

  const plotType = $('#plot_type_select').select2('data')[0].id;
  plotConfig.plot_type = plotType;

  // Update filters based on selection
  obsFilters = {};
  for (const property in obsLevels) {
    const propData = $(`#${property}_dropdown`).select2('data');
    obsFilters[property] = propData.map((elem) => elem.id);

    // If no groups for an observation are selected, delete filter
    if (!obsFilters[property].length) {
      delete obsFilters[property];
    }
  }
  plotConfig.obs_filters = obsFilters;

  if (!plotType) {
    window.alert('Please select a plot type.');
    return;
  }

  plotConfig.gene_symbols = genesFilters = $('#gene_dropdown').select2('data').map((elem) => elem.id);

  if (!["volcano", "quadrant"].includes(plotType)) {
    if (Object.keys(obsFilters).length) {
      window.alert('At least one observation must have categories filtered.');
      return;
    }
  }

  let sortOrder = {};
  let categoriesUsed = [];
  if (sortCategories.axis || sortCategories.groupby) {
    // Grab the sorted order of the list and convert to array
    if (sortCategories.axis) {
      sortOrder[sortCategories.axis] = $('#axis_sortable').sortable("toArray", {attribute:"value"});
      categoriesUsed.push(sortCategories.axis);
    }
    if (sortCategories.groupby) {
      sortOrder[sortCategories.groupby] = $('#groupby_sortable').sortable("toArray", {attribute:"value"});
      categoriesUsed.push(sortCategories.groupby);
    }

  }
  // This should be rare, but just use the axis order if both are the same category
  if (sortCategories.facet && !(cats_used.includes(sortCategories.facet))) {
    sortOrder[sortCategories.facet] = $('#facet_sortable').sortable("toArray", {attribute:"value"});
  }
  plotConfig.sort_order = sortOrder;

  // Add specific plotConfig options depending on plot type
  switch (plotType) {
    case 'dotplot':
      if ((plotConfig.gene_symbols).length < 1) {
        window.alert('At least one gene must be provided.');
        return;
      }
      plotConfig.groupby_filter = $('input[name="obs_groupby"]:checked').val();
      if (!plotConfig.groupby_filter) {
        window.alert("Must select a groupby filter for dot plots.");
        return;
      }
      plotConfig.groupby_filter = $('input[name="obs_groupby"]:checked').val();
      plotConfig.facet_col = $('input[name="obs_facet"]:checked').val();
      break;
    case 'heatmap':
      if ((plotConfig.gene_symbols).length < 2) {
        window.alert("Must select at least 2 genes to generate a heatmap");
        return;
      }
      plotConfig.clusterbar_fields = [];
      $('input[name="obs_clusterbar"]:checked').each( (idx, elem) => {
        plotConfig.clusterbar_fields.push($(elem).val());
      });
      plotConfig.matrixplot = $('#matrixplot').is(':checked');
      plotConfig.axis_sort_col = $('input[name="obs_sort"]:checked').val();
      plotConfig.facet_col = $('input[name="obs_facet"]:checked').val();
      if (plotConfig.matrixplot && !plotConfig.axis_sort_col) {
        window.alert("Must choose a 'Sort observations on axis' category to aggregate means for the matrixplot.");
        return
      }
      plotConfig.cluster_obs = $('#cluster_obs').is(':checked');
      plotConfig.cluster_genes = $('#cluster_genes').is(':checked');
      plotConfig.flip_axes = $('#flip_axes').is(':checked');
      plotConfig.distance_metric = $('#distance_select').select2('data')[0].id;
      break;
    case 'mg_violin':
      if ((plotConfig.gene_symbols).length < 1) {
        window.alert('At least one gene must be provided.');
        return;
      }
      plotConfig.groupby_filter = $('input[name="obs_groupby"]:checked').val();
      if (!plotConfig.groupby_filter) {
        window.alert("Must select a groupby filter for violin plots.");
        return;
      }
      plotConfig.facet_col = $('input[name="obs_facet"]:checked').val();
      plotConfig.stacked_violin = $('#stacked_violin').is(':checked');
      plotConfig.violin_add_points = $('#violin_add_points').is(':checked');
      break;
    case 'quadrant':
      plotConfig.include_zero_fc = $('#include_zero_foldchange').is(':checked');
      plotConfig.fold_change_cutoff = Number($("#quadrant_foldchange_cutoff").val());
      plotConfig.fdr_cutoff = Number($("#quadrant_fdr_cutoff").val());
      plotConfig.de_test_algo = $('#de_test_select').select2('data')[0].id;
      if (! plotConfig.de_test_algo) {
        window.alert('Must select a DE statistical test.');
        return;
      }
      plotConfig.compare1_condition = $('#quadrant_compare1_condition').select2('data')[0].id;
      plotConfig.compare2_condition = $('#quadrant_compare2_condition').select2('data')[0].id;
      plotConfig.ref_condition = $('#quadrant_ref_condition').select2('data')[0].id;
      break;
    default:
      // volcano
      plotConfig.adjust_pvals = $('#adj_pvals').is(':checked');
      plotConfig.annotate_nonsignificant = $('#annot_nonsig').is(':checked');
      plotConfig.de_test_algo = $('#de_test_select').select2('data')[0].id;
      if (! plotConfig.de_test_algo) {
        window.alert('Must select a DE statistical test.');
        return;
      }
      plotConfig.query_condition = $('#volcano_query_condition').select2('data')[0].id;
      plotConfig.ref_condition = $('#volcano_ref_condition').select2('data')[0].id;
      // Validation related to the conditions
      if (!(plotConfig.query_condition && plotConfig.ref_condition)) {
        window.alert('Both comparision conditions must be selected to generate a volcano plot.');
        return;
      }
      const conditionKey = plotConfig.query_condition.split(';-;')[0];
      if (plotConfig.query_condition.split(';-;')[0] !== plotConfig.ref_condition.split(';-;')[0]) {
        window.alert('Please choose 2 conditions from the same observation group.');
        return;
      }

      if (plotConfig.query_condition.split(';-;')[1] === plotConfig.ref_condition.split(';-;')[1]) {
        window.alert('Please choose 2 different conditions.');
        return;
      }

      // If condition category was filtered, the selected groups must be present
      if (conditionKey in obsFilters
        && !(obsFilters[conditionKey].includes(plotConfig.query_condition.split(';-;')[1])
        && obsFilters[conditionKey].includes(plotConfig.ref_condition.split(';-;')[1]))) {
        window.alert('Condition observation is found in filters list, but one or both condition groups is filtered out. Please adjust.');
        return;
      }
  }

  // Render dataset plot HTML
  const plotTemplate = $.templates('#dataset_plot_tmpl');
  const plotHtml = plotTemplate.render({ dataset_id: datasetId });
  $('#dataset_plot').html(plotHtml);

  // Draw the updated chart
  $('#dataset_spinner').show();
  await draw(datasetId, plotConfig);
  $('#dataset_spinner').hide();

  // Show plot options and disable selected genes button (since genes are not selected anymore)
  $('#post_plot_options').show();
  $("#selected_genes_btn").prop("disabled", true);
});

// If "all" button is clicked, populate dropdown with all groups in this observation
$(document).on('click', '.all', function () {
  const id = this.id;
  const group = id.replace('_all', '');

  $(`#${group}_dropdown`).val(obsLevels[group]);
  $(`#${group}_dropdown`).trigger('change'); // This actually triggers select2 to show the dropdown vals
});

// Reset observation filters choices to be empty
$(document).on('click', '#reset_opts', async function () {
  $('#options_container').show();
  $('#options_spinner').show();

  // Reset sorting order
  sortCategories = {"facet": null, "axis": null, "groupby": null};
  $(`#axis_order_label`).hide();
  $(`#facet_order_label`).hide();

  // Update fields dependent on dataset observations
  createObsFilterDropdowns(obsLevels);
  switch ($('#plot_type_select').val()) {
    case 'dotplot':
      createDotplotDropdowns(obsLevels);
      break;
    case 'heatmap':
      createHeatmapDropdowns(obsLevels);
      break;
    case 'mg_violin':
      createViolinDropdowns(obsLevels);
      break;
    case 'quadrant':
      createQuadrantDropdowns(obsLevels);
      break;
    default:
      // volcano
      createVolcanoDropdowns(obsLevels);
  }

  $('.all').click();  // Include all groups for every category (filter nothing)

  $('#options_spinner').hide();
});

// If advanced options collapsable is clicked, toggle arrow b/t up and down
$(document).on('click', '#advanced_options_button', () => {
  if ($('#advanced_options_arrow').hasClass('fa-caret-down')) {
    $('#advanced_options_arrow').removeClass('fa-caret-down');
    $('#advanced_options_arrow').addClass('fa-caret-up');
  } else {
    $('#advanced_options_arrow').addClass('fa-caret-down');
    $('#advanced_options_arrow').removeClass('fa-caret-up');
  }
});

// Save plot
$(document).on('click', '#save_display_btn', async function () {
  $('#saved_plot_confirmation').hide();
  $('#saved_plot_confirmation').removeClass('text-success');
  $('#saved_plot_confirmation').removeClass('text-danger');

  const plotType = $('#plot_type_select').select2('data')[0].id;

  const payload = {
    id: null, // Want to save as a new display
    dataset_id: datasetId,
    user_id: CURRENT_USER.id,
    label: $('#display_name').val(),
    plot_type: plotType,
    plotly_config: JSON.stringify({
      // depending on display type, this object will
      // have different properties
      ...plotConfig
    })
  };

  const res = await $.ajax({
    url: './cgi/save_dataset_display.cgi',
    type: 'POST',
    data: payload,
    dataType: 'json'
  });

  if (res?.success) {

    let msg = 'Plot successfully saved'

    if ($("#save_as_default_check").is(':checked') && res.display_id) {
      displayId = res.display_id;
      const res2 = await $.ajax({
        url: './cgi/save_default_display.cgi',
        type: 'POST',
        data: {
          user_id: CURRENT_USER.id,
          dataset_id: datasetId,
          display_id: displayId,
          is_multigene: 1
        },
        dataType: 'json'
      });

      if (res2?.success) {
          // Swap current default buttons
        $('.js-current-default')
          .prop('disabled', false)
          .addClass('js-save-default')
          .addClass('btn-purple')
          .removeClass('btn-secondary')
          .removeClass('js-current-default')
          .text("Make Default");
        $(`#${displayId}_default`)
          .prop('disabled', true)
          .removeClass('js-save-default')
          .removeClass('btn-purple')
          .addClass('btn-secondary')
          .addClass('js-current-default')
          .text("Default");
          msg += " and was set as the default display.";
      } else {
          msg += " but there was an issue saving as the default display.";
      }
    } else {
      msg += " but not set as default display."
    }

    $('#saved_plot_confirmation').text(msg);
    $('#saved_plot_confirmation').addClass('text-success');
    $('#saved_plot_confirmation').show();
  } else {
    $('#saved_plot_confirmation').text('There was an issue saving the plot');
    $('#saved_plot_confirmation').addClass('text-danger');
    $('#saved_plot_confirmation').show();
  }

  // Update saved displays modal so new plot is included
  loadSavedDisplays(datasetId);
});

// Load display information back into the curator page
$(document).on('click', '.js-load-display', async function () {
  const id = this.id;
  displayId = id.replace('_load', '');

  const display = await $.ajax({
    url: './cgi/get_dataset_display.cgi',
    type: 'POST',
    data: { display_id: displayId },
    dataType: 'json'
  });

  plotConfig = display.plotly_config;

  // Load plot type
  $('#plot_type_select').val(display.plot_type);
  $('#plot_type_select').trigger('change');

  // Load gene symbols
  geneSymbols = plotConfig.gene_symbols;
  $('#gene_dropdown').val(geneSymbols);
  $('#gene_dropdown').trigger('change');

  // Load config options
  loadDisplayConfigHtml(plotConfig);

  // Hide modal box
  $('#load_plots_modal').modal('hide');

  // Draw the updated chart
  $('#dataset_spinner').show();
  const plotTemplate = $.templates('#dataset_plot_tmpl');
  const plotHtml = plotTemplate.render({ dataset_id: datasetId });
  $('#dataset_plot').html(plotHtml);
  await draw(datasetId, plotConfig);
  $('#dataset_spinner').hide();

  // Show plot options
  $('#post_plot_options').show();
  $("#selected_genes_btn").prop("disabled", true);

});

// Delete user display
$(document).on('click', '.js-delete-display', async function () {
  $('#delete_display_confirmation').hide();
  $('#delete_display_confirmation').removeClass('alert-success');
  $('#delete_display_confirmation').removeClass('alert-danger');

  const id = this.id;
  const displayId = id.replace('_delete', '');

  const res = await $.ajax({
    url: './cgi/delete_dataset_display.cgi',
    type: 'POST',
    data: { id: displayId, user_id: CURRENT_USER.id },
    dataType: 'json'
  });

  $('#load_plots_modal').modal('hide');
  $('#delete_display_confirmation').show();
  if (res?.success) {
    $('#delete_display_confirmation_text').text('Display was successfully deleted.');
    $('#delete_display_confirmation').addClass('alert-success');
  } else {
    $('#delete_display_confirmation_text').text('There was an issue deleting the saved display.');
    $('#delete_display_confirmation').addClass('alert-danger');
  }

  // Update saved displays, now that display has been deleted
  loadSavedDisplays(datasetId);
});

// Save this particular display as the user's default display
$(document).on('click', '.js-save-default', async function () {
  $('#saved_default_confirmation').hide();
  $('#saved_default_confirmation').removeClass('alert-success');
  $('#saved_default_confirmation').removeClass('alert-danger');

  const id = this.id;
  const displayId = id.replace('_default', '');

  const res = await $.ajax({
    url: './cgi/save_default_display.cgi',
    type: 'POST',
    data: {
      user_id: CURRENT_USER.id,
      dataset_id: datasetId,
      display_id: displayId,
      is_multigene: 1
    },
    dataType: 'json'
  });

  // Swap current default buttons
  $('.js-current-default')
    .prop('disabled', false)
    .addClass('js-save-default')
    .addClass('btn-purple')
    .removeClass('btn-secondary')
    .removeClass('js-current-default')
    .text("Make Default")
  $(`#${displayId}_default`)
    .prop('disabled', true)
    .removeClass('js-save-default')
    .removeClass('btn-purple')
    .addClass('btn-secondary')
    .addClass('js-current-default')
    .text("Default")

  $('#saved_default_confirmation').show();
  if (res?.success) {
    $('#saved_default_confirmation_text').text('Display successfully saved as your new default.');
    $('#saved_default_confirmation').addClass('alert-success');
  } else {
    $('#saved_default_confirmation_text').text('There was an issue setting the default display.');
    $('#saved_default_confirmation').addClass('alert-danger');
  }
});
