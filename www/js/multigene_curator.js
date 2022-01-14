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
let sortCategories = {"primary": null, "secondary": null}; // Control sorting order of chosen categories
let genesFilter = [];

let plotConfig = {};  // Plot config that is passed to API or stored in DB

let selectedGenes = null;  // Genes selected in a plot (e.g. volcano with lasso tool)

let datasetId = null;
let defaultDisplayId = null;
let displayId = null;
let obsLevels = null;
let geneSymbols = null;

const datasetTree = new DatasetTree({treeDiv: '#dataset_tree'});
const geneCartTree = new GeneCartTree({treeDiv: '#gene_cart_tree'});

const plotTypes = ['dotplot', 'heatmap', 'mg_violin', 'quadrant', 'volcano'];

const dotplotOptsIds = ["#obs_primary_container", "#obs_secondary_container"];
const heatmapOptsIds = ["#heatmap_options_container", "#obs_primary_container", "#obs_secondary_container"];
const quadrantOptsIds = ["#quadrant_options_container", "#de_test_container"];
const violinOptsIds = ["#violin_options_container", "#obs_primary_container", "#obs_secondary_container"];
const volcanoOptsIds = ["#volcano_options_container", "#de_test_container"];

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
window.onload=() => {
  // Hide further configs until a dataset is chosen.
  // Changing the dataset will start triggering these to show
  $('#plot_type_container').hide();
  $('#gene_container').hide();

  // check if the user is already logged in
  check_for_login();

  // Load gene carts and datasets before the dropdown appears
  reloadTrees ();

  // Initialize plot types
   $('#plot_type_select').select2({
    placeholder: 'Choose how to plot',
    width: '25%'
  });

  // If brought here by the "gene search results" page, curate on the dataset ID that referred us
  const linkedDatasetId = getUrlParameter("dataset_id");
  if (linkedDatasetId) {
    $('#dataset').val(linkedDatasetId);
    try {
        $('#dataset').text(dataset_tree.treeData.find(e => e.dataset_id === linkedDatasetId).text);
        $("#dataset").trigger("change");
    } catch {
        console.error(`Dataset id ${linkedDatasetId} was not returned as a public/private/shared dataset`);
    }
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
};

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
  $('#pre_dataset_spinner').show();
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
      const userDatasets = [];
      if (data.user.datasets.length > 0) {
        // User has some profiles
        $.each(data.user.datasets, (_i, item) => {
          if (item) {
            userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
          }
        });
      }
      // Next, add datasets shared with the user
      const sharedDatasets = [];
      if (data.shared_with_user.datasets.length > 0) {
        $.each(data.shared_with_user.datasets, (_i, item) => {
          if (item) {
            sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
          }
        });
      }
      // Now, add public datasets
      const domainDatasets = [];
      if (data.public.datasets.length > 0) {
        $.each(data.public.datasets, (_i, item) => {
          if (item) {
            domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
          }
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
  $('#pre_dataset_spinner').hide();
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

  // Make some complex edits to the plotly layout
  const plotType = $('#plot_type_select').select2('data')[0].id;
  if (plotType === 'heatmap') {
    setHeatmapHeightBasedOnGenes(plotlyJson.layout, genesFilter);
  } else if (plotType=== "mg_violin" && $("#stacked_violin").is(":checked")){
    adjustStackedViolinHeight(plotlyJson.layout);
  }

  const configMods = {
    responsive: true
  };

  const config = {
    ...plotlyConfig,
    ...configMods
  };
  Plotly.newPlot(targetDiv, plotlyJson.data, plotlyJson.layout, config);

  // Update plot with custom plot config stuff stored in plot_display_config.js
  const curator_conf = post_plotly_config.curator;
  for (const idx in curator_conf) {
    const conf = curator_conf[idx];
    // Get config (data and/or layout info) for the plot type chosen, if it exists
    if (conf.plot_type == plotType) {
      const update_data = "data" in conf ? conf.data : {};
      const update_layout = "layout" in conf ? conf.layout : {};
      Plotly.update(targetDiv, update_data, update_layout)
    }
  }


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

// Render the observation primary field HTML
function createObsPrimaryField (obsLevels) {
  const tmpl = $.templates('#obs_primary_tmpl'); // Get compiled template using jQuery selector for the script block
  const html = tmpl.render(obsLevels); // Render template using data - as HTML string
  $('#obs_primary_container').html(html); // Insert HTML string into DOM
}

// Render the observation secondary field
function createObsSecondaryField (obsLevels) {
  const tmpl = $.templates('#obs_secondary_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_secondary_container').html(html);
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

// Render clusterbars for heatmap
function createObsClusterBarField (obsLevels) {
  const tmpl = $.templates('#obs_clusterbar_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_clusterbar_container').html(html);
}

// Render the sortable list for the chosen category
function createObsSortable (obsLevel, scope) {
  const escapedObsLevel = $.escapeSelector(obsLevel);
  const propData = $(`#${escapedObsLevel}_dropdown`).select2('data');
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
  createObsPrimaryField (obsLevels);
  $('#none_primary_group').hide();  // Primary is required
  createObsSecondaryField(obsLevels);
}

// Render dropdowns specific to the heatmap plot
function createHeatmapDropdowns (obsLevels) {
  createObsPrimaryField(obsLevels);
  createObsSecondaryField(obsLevels);
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
  // Filter only categories with at least 3 conditions
  const goodObsLevels = {}
  for (const category in obsLevels) {
    if (obsLevels[category].length >= 3) {
      goodObsLevels[category] = obsLevels[category]
    }
  }

  const tmpl = $.templates('#select_conditions_tmpl');
  const html = tmpl.render(goodObsLevels);
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
  createObsPrimaryField (obsLevels);
  $('#none_primary_group').hide();  // Primary is required
  createObsSecondaryField(obsLevels);
}

// Render dropdowns specific to the volcano plot
function createVolcanoDropdowns (obsLevels) {
  // Filter only categories wtih at least 2 conditions
  const goodObsLevels = {}
  for (const category in obsLevels) {
    if (obsLevels[category].length >= 2) {
      goodObsLevels[category] = obsLevels[category]
    }
  }

  const tmpl = $.templates('#select_conditions_tmpl');
  const html = tmpl.render(goodObsLevels);
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
    const escapedProperty = $.escapeSelector(property);
    $(`#${escapedProperty}_dropdown`).val(obsFilters[property]);
    $(`#${escapedProperty}_dropdown`).trigger('change');
  }

  const escapedPrimary = $.escapeSelector(plotConfig.primary_col);
  const escapedSecondary = $.escapeSelector(plotConfig.secondary_col);

  $(`#${escapedPrimary}_primary`).prop('checked', true).click();
  $(`#${escapedSecondary}_secondary`).prop('checked', true).click();

  $('#plot_title').val(plotConfig.plot_title);

  // Populate plot type-specific dropdowns and checkbox options
  switch ($('#plot_type_select').val()) {
    case 'dotplot':
      break;
    case 'heatmap':
      for (const field in plotConfig.clusterbar_fields) {
        const escapedField = $.escapeSelector(field);
        $(`#${escapedField}_clusterbar`).prop('checked', true).click();
      }
      $('#matrixplot').prop('checked', plotConfig.matrixplot);
      $('#cluster_obs').prop('checked', plotConfig.cluster_obs);
      $('#cluster_genes').prop('checked', plotConfig.cluster_obs);
      $('#flip_axes').prop('checked', plotConfig.flip_axes);
      $('#distance_select').val(plotConfig.distance_metric);
      $('#distance_select').trigger('change');
      break;
    case 'mg_violin':
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
        const domainGeneCarts = [];
        const userGeneCarts = [];

        if (data.gene_carts.length > 0) {
          // User has some profiles
          $.each(data.gene_carts, (_i, item) => {
            userGeneCarts.push({ value: item.id, text: item.label });
          });
        }

          // Now, add public gene carts
        if (data.domain_gene_carts.length > 0) {
          $.each(data.domain_gene_carts, (_i, item) => {
            domainGeneCarts.push({ value: item.id, text: item.label });
          });
        }

        geneCartTree.domainGeneCarts = domainGeneCarts;
        geneCartTree.userGeneCarts = userGeneCarts;
        geneCartTree.generateTree();

        if (userGeneCarts || domainGeneCarts) {
          $('#gene_cart_container').show();
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
  $('#gene_cart_clear').click();

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

  // Determine if at least one category has at least two groups (for volcanos) or at least three groups (for quadrants)
  // If condition is not met, disable the plot option
  let hasTwoGroupCategory = false;
  let hasThreeGroupCategory = false;

  for (const category in obsLevels) {
    if (obsLevels[category].length >= 3) {
      hasThreeGroupCategory = true;
    } else if (obsLevels[category].length >= 2) {
      hasTwoGroupCategory = true;
    }
    if (hasTwoGroupCategory && hasThreeGroupCategory) {
      break;
    }
  }
  $('#volcano_opt').prop("disabled", false);
  if (!hasTwoGroupCategory) {
    $('#volcano_opt').prop("disabled", true);
  }
  $('#quadrant_opt').prop("disabled", false);
  if (!hasThreeGroupCategory) {
    $('#quadrant_opt').prop("disabled", true);
  }

  $('#create_plot').show();
  $("#create_plot").prop("disabled", true);
  $('#reset_opts').show();

  // If a plot type was already selected,
  // reset the options so configs are populated for the current dataset
  if ($('#plot_type_select').val() ) {
    $("#create_plot").prop("disabled", false);
    $('#reset_opts').click();
  }
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
  const geneCartId = $(this).val();
  const params = { session_id, gene_cart_id: geneCartId };
  const d = new $.Deferred(); // Causes editable to wait until results are returned
  // User is not logged in
  if (!session_id) {
    d.resolve();
  } else {
    // User is logged in

    $('#gene_spinner').show();
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
    $('#gene_spinner').hide();
  }
  return d.promise();
});

// Some options are specific to certain plot types
$('#plot_type_select').change(() => {
  $('#reset_opts').click();  // Reset all options
  $('#selected_genes_field').hide();
  $("#create_plot").prop("disabled", false);

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
  const genesFilter = $('#gene_dropdown').select2('data').map((elem) => elem.id);
  // Show warning if too many genes are entered
  $("#too_many_genes_warning").hide();
  if (genesFilter.length > 10) {
    $("#too_many_genes_warning").text(`There are currently ${genesFilter.length} genes to be plotted. Be aware that with some plots, a high number of genes can make the plot congested or unreadable.`);
    $("#too_many_genes_warning").show();
  }

  // Cannot cluster columns with just one gene (because function is only available
  // in dash.clustergram which requires 2 or more genes in plot)
  if (genesFilter.length > 1) {
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
$(document).on('change', 'input[name="obs_primary"]', () => {
  const obsLevel = $('input[name="obs_primary"]:checked').val();
  if (obsLevel !== "none") {
    createObsSortable(obsLevel, "primary");
  } else {
    $('#primary_sortable').empty();
  }
});

$(document).on('change', 'input[name="obs_secondary"]', (e) => {
  const obsLevel = $('input[name="obs_secondary"]:checked').val();
  if ( obsLevel !== "none") {
    createObsSortable(obsLevel, "secondary");
  } else {
    $('#secondary_sortable').empty();
  }
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
    const escapedProperty = $.escapeSelector(property);
    const propData = $(`#${escapedProperty}_dropdown`).select2('data');
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

  plotConfig.gene_symbols = genesFilter = $('#gene_dropdown').select2('data').map((elem) => elem.id);

  const sortOrder = {};
  const categoriesUsed = [];
  // Grab the sorted order of the list and convert to array
  if (sortCategories.primary) {
    sortOrder[sortCategories.primary] = $('#primary_sortable').sortable("toArray", {attribute:"value"});
    categoriesUsed.push(sortCategories.primary);
  }
  // This should be rare, but just use the primary order if both are the same category
  if (sortCategories.secondary && !(categoriesUsed.includes(sortCategories.secondary))) {
    sortOrder[sortCategories.secondary] = $('#secondary_sortable').sortable("toArray", {attribute:"value"});
  }
  plotConfig.sort_order = sortOrder;

  if ($('input[name="obs_primary"]:checked').val() !== "none"){
    plotConfig.primary_col = $('input[name="obs_primary"]:checked').val();
  }
  if ($('input[name="obs_secondary"]:checked').val() !== "none"){
    plotConfig.secondary_col = $('input[name="obs_secondary"]:checked').val();
  }

  plotConfig.plot_title = $('#plot_title').val();

  // Add specific plotConfig options depending on plot type
  switch (plotType) {
    case 'dotplot':
      if ((plotConfig.gene_symbols).length < 1) {
        window.alert('At least one gene must be provided.');
        return;
      }
      if (!plotConfig.primary_col) {
        window.alert("Must select at least a primary category to aggregate groups by for dot plots.");
        return;
      }
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
      if (plotConfig.matrixplot && !plotConfig.primary_col) {
        window.alert("Must choose at least a primary category to aggregate means for the matrixplot.");
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
      if (!plotConfig.primary_col) {
        window.alert("Must select at least a primary category to aggregate groups by for violin plots.");
        return;
      }
      plotConfig.stacked_violin = $('#stacked_violin').is(':checked');
      plotConfig.violin_add_points = $('#violin_add_points').is(':checked');
      break;
    case 'quadrant':
      plotConfig.include_zero_fc = $('#include_zero_foldchange').is(':checked');
      plotConfig.fold_change_cutoff = Number($("#quadrant_foldchange_cutoff").val());
      plotConfig.fdr_cutoff = Number($("#quadrant_fdr_cutoff").val());
      if (! $('#de_test_select').select2('data')[0].id) {
        window.alert('Must select a DE statistical test.');
        return;
      }
      plotConfig.de_test_algo = $('#de_test_select').select2('data')[0].id;
      plotConfig.compare1_condition = $('#quadrant_compare1_condition').select2('data')[0].id;
      plotConfig.compare2_condition = $('#quadrant_compare2_condition').select2('data')[0].id;
      plotConfig.ref_condition = $('#quadrant_ref_condition').select2('data')[0].id;
      break;
    default:
      // volcano
      plotConfig.adj_pvals = $('#adj_pvals').is(':checked');
      plotConfig.annotate_nonsignificant = $('#annot_nonsig').is(':checked');
      if (! $('#de_test_select').select2('data')[0].id) {
        window.alert('Must select a DE statistical test.');
        return;
      }
      plotConfig.de_test_algo = $('#de_test_select').select2('data')[0].id;
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
$(document).on('click', '.js-all', function () {
  const id = this.id;
  const group = id.replace('_all', '');
  const escapedGroup = $.escapeSelector(group);

  $(`#${escapedGroup}_dropdown`).val(obsLevels[group]);
  $(`#${escapedGroup}_dropdown`).trigger('change'); // This actually triggers select2 to show the dropdown vals
});

// Clear gene cart
$(document).on('click', '#gene_cart_clear', () => {
  $('#gene_dropdown').val('');  // Clear genes
  $('#gene_dropdown').trigger('change');
  $('#gene_cart').text('Choose gene cart'); // Reset gene cart text
  $('#gene_cart').val('');
});

// Reset observation filters choices to be empty
$(document).on('click', '#reset_opts', async function () {
  $('#options_container').show();
  $('#options_spinner').show();

  // Reset sorting order
  sortCategories = {"primary": null, "secondary": null};
  $(`#primary_order_label`).hide();
  $(`#secondary_order_label`).hide();

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

  $('.js-all').click();  // Include all groups for every category (filter nothing)

  $('#options_spinner').hide();
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
  } else {
    $('#saved_plot_confirmation').text('There was an issue saving the plot');
    $('#saved_plot_confirmation').addClass('text-danger');
  }
  $('#saved_plot_confirmation').show();

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
  $('#plot_spinner').show();
  const plotTemplate = $.templates('#dataset_plot_tmpl');
  const plotHtml = plotTemplate.render({ dataset_id: datasetId });
  $('#dataset_plot').html(plotHtml);
  await draw(datasetId, plotConfig);
  $('#plot_spinner').hide();

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
