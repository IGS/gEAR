/*
SAdkins note - The styling on this page differs from other gEAR JS pages
I was trying to follow the recommended Javascript style in using camelCase
in addition to using the StandardJS style ("semistandrd" variant)
However, code inherited from common.js is still in snake_case rather than camelCase
*/

'use strict';
/* global $, axios, Plotly, CURRENT_USER, session_id, check_for_login */

let obsFilters = {};
let genesFilters = [];
let supplementaryGenesFilters = [];

let plotConfig = {};

let datasetId = null;
let displayId = null;
let obsLevels = null;
let geneSymbols = null;

// Async to ensure data is fetched before proceeding
(async () => {
  // check if the user is already logged in
  await check_for_login();

  // Initialize tooltips
  $(function () {
    $('[data-toggle="tooltip"]').tooltip();
  });

  // Initialize datasets available to the user
  $('#dataset_select').select2({
    placeholder: 'To search, click to select or start typing a dataset name'
  });
  await populateDatasets();

  // Initialize plot types
  $('#plot_type_select').select2({
    placeholder: 'Choose how to plot'
  });

  // Hide further configs until a dataset is chosen
  $('#plot_type_container').hide();
  $('#advanced_options_container').hide();
  $('#gene_container').hide();
})();

// Call API to return plot JSON data
async function getData (datasetId, payload) {
  return await axios.post(`/api/plot/${datasetId}/mg_dash`, {
    ...payload
  });
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

async function populateDatasets () {
  $.ajax({
    type: 'POST',
    url: './cgi/get_h5ad_dataset_list.cgi',
    data: {
      session_id
    },
    dataType: 'json',
    success: function (data) {
      // Populate select box with dataset information owned by the user
      if (data.user.datasets.length > 0) {
        const userDatasetListTmpl = $.templates('#dataset_list_tmpl');
        const userDatasetListHtml = userDatasetListTmpl.render(
          data.user.datasets
        );
        $('#dataset_ids_user').html(userDatasetListHtml);
      } else {
        $('#dataset_id .user-initial').html('Not logged in');
      }

      // Next, add datasets shared with the user
      if (data.shared_with_user.datasets.length > 0) {
        const sharedWithUserDatasetListTmpl = $.templates(
          '#dataset_list_tmpl'
        );
        const sharedWithUserDatasetListHtml = sharedWithUserDatasetListTmpl.render(
          data.shared_with_user.datasets
        );
        $('#dataset_ids_shared_with_user').html(
          sharedWithUserDatasetListHtml
        );
      }

      // Now, add public datasets
      if (data.public.datasets.length > 0) {
        const publicDatasetListTmpl = $.templates('#dataset_list_tmpl');
        const publicDatasetListHtml = publicDatasetListTmpl.render(
          data.public.datasets
        );
        $('#dataset_ids_public').html(publicDatasetListHtml);
      }
    },
    error: function (xhr, status, msg) {
      console.error('Failed to load dataset list because msg: ' + msg);
    }
  });
}

// Draw plotly chart to image
async function drawPreviewImage (display) {
  // check if config has been stringified
  let config;
  if (typeof display.plotly_config === 'string') {
    config = JSON.parse(display.plotly_config);
  } else {
    config = display.plotly_config;
  }

  const { data } = await getData(datasetId, config);
  const { plot_json: plotJson, plot_config: plotConfig } = data;
  Plotly.toImage(
    { ...plotJson, plotConfig },
    { height: 500, width: 500 }
  ).then(url => {
    $(`#modal-display-img-${display.id}`).prop('src', url);
  }).then(
    () => { $(`#modal-display-${display.id}-loading`).hide(); }
  );
}

// Draw plotly chart in HTML
function drawChart (data, datasetId, supplementary = false) {
  const targetDiv = supplementary ? `dataset_${datasetId}_secondary` : `dataset_${datasetId}_h5ad`;
  const { plot_json: plotJson, plot_config: plotConfig, message, success} = data;

  const layoutMods = {
    height: targetDiv.clientHeight,
    width: targetDiv.clientWidth
  };

  // NOTE: This will definitely affect the layout on the gene search results page
  // if the "height" style for the container in CSS is removed.
  if (!supplementary && $('#plot_type_select').select2('data')[0].id === 'heatmap') {
    if (genesFilters.length > 50) {
      layoutMods.height = genesFilters.length * 10;
    }
  }

  // Overwrite plot layout and config values with custom ones from display
  const layout = {
    ...plotJson.layout,
    ...layoutMods
  };

  const configMods = {
    responsive: false
  };

  const config = {
    ...plotConfig,
    ...configMods
  };
  Plotly.newPlot(targetDiv, plotJson.data, layout, config);

  if (message) {
    if (success < 1) {
      $(targetDiv + '.js-plot-error').text(message);
      $(targetDiv + '.js-plot-error').show();
    } else {
      $(targetDiv + '.js-plot-warning').text(message);
      $(targetDiv + '.js-plot-warning').show();
    }

  }
}

// Submit API request and draw the HTML
async function draw (datasetId, payload, supplementary = false) {
  const {
    data
  } = await getData(datasetId, payload);
  drawChart(data, datasetId, supplementary);
}

function createGeneDropdown (genes) {
  const tmpl = $.templates('#gene_dropdown_tmpl');
  const data = { genes: genes };
  const html = tmpl.render(data);
  $('#gene_dropdown_container').html(html);
  $('#gene_dropdown').select2({
    placeholder: 'To search, click to select or start typing some gene names',
    allowClear: true
  });
}

// Render the observation groupby field HTML
function createObsGroupbyField (obsLevels) {
  const tmpl = $.templates('#obs_groupby_tmpl'); // Get compiled template using jQuery selector for the script block
  const html = tmpl.render(obsLevels); // Render template using data - as HTML string
  $('#obs_groupby_container').html(html); // Insert HTML string into DOM
}

// Render the observation filter dropdowns
function createObsDropdowns (obsLevels) {
  const tmpl = $.templates('#obs_dropdowns_tmpl');
  const html = tmpl.render(obsLevels);
  $('#obs_dropdowns_container').html(html);
  $('select.js-obs-levels').select2({
    placeholder: 'Start typing to filter categories. Click "All" to use all categories',
    allowClear: true
  });
}

// Render the volcano condition selection dropdown
function createVolcanoDropdowns (obsLevels) {
  const tmpl = $.templates('#volcano_options_tmpl');
  const html = tmpl.render(obsLevels);
  $('#volcano_condition1').html(html);
  $('#volcano_condition1').select2({
    placeholder: 'Select the first condition to compare with.'
  });
  $('#volcano_condition2').html(html);
  $('#volcano_condition2').select2({
    placeholder: 'Select the second condition to compare with.'
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
async function loadSavedDisplays (datasetId) {
  const datasetData = await fetchDatasetInfo(datasetId);
  const { owner_id: ownerId } = datasetData;
  const userDisplays = await fetchUserDisplays(CURRENT_USER.id, datasetId);
  const ownerDisplays = await fetchOwnerDisplays(ownerId, datasetId);

  // Filter displays to those only with multigene plot types
  const mgUserDisplays = userDisplays.filter(d => ['heatmap', 'mg_violin', 'volcano'].includes(d.plot_type));
  const mgOwnerDisplays = ownerDisplays.filter(d => ['heatmap', 'mg_violin', 'volcano'].includes(d.plot_type));

  // TODO: Determine default display and hide "make_default" button

  mgUserDisplays.forEach(display => {
    drawPreviewImage(display);
  });
  mgOwnerDisplays.forEach(display => {
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
// Create a few of the dropdown options
  createObsGroupbyField(obsLevels);
  createObsDropdowns(obsLevels);
  createVolcanoDropdowns(obsLevels);

  // Populate any checkboxes and radio buttons
  $('#cluster_cols').prop('checked', plotConfig.cluster_cols);
  $('#adj_pvals').prop('checked', plotConfig.adj_pvals);
  $(`#${plotConfig.groupby_filter}_groupby`).prop('checked', true).click();

  // Populate filter-by dropdowns
  obsFilters = plotConfig.obs_filters;
  for (const property in obsFilters) {
    $(`#${property}_dropdown`).val(obsFilters[property]);
    $(`#${property}_dropdown`).trigger('change');
  }

  // Populate volcano conditions.
  $('#volcano_condition1').val(plotConfig.condition1);
  $('#volcano_condition1').trigger('change');
  $('#volcano_condition2').val(plotConfig.condition2);
  $('#volcano_condition2').trigger('change');
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
      success: function (data, textStatus, jqXHR) { // source https://stackoverflow.com/a/20915207/2900840
        const userGeneCarts = [];

        if (data.gene_carts.length > 0) {
          // User has some profiles
          $.each(data.gene_carts, function (i, item) {
            userGeneCarts.push({ value: item.id, text: item.label });
          });
          $('#gene_cart_container').show();
          const tmpl = $.templates('#gene_cart_tmpl');
          const html = tmpl.render(userGeneCarts);
          $('#selected_gene_cart').html(html);
          // Add blank option so select2 will show placeholder
          $('#selected_gene_cart').prepend("<option value=''></option>").val('');
          $('#selected_gene_cart').select2({
            placeholder: 'Select a preloaded set of genes.',
            allowClear: true
          });
        } else {
          $('#gene_cart_container').hide();
        }
        d.resolve();
      },
      error: function (jqXHR, textStatus, errorThrown) {
        // display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        d.fail();
      }
    });
  }
  d.promise();
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

$('#dataset_select').change(async function () {
  datasetId = $('#dataset_select').select2('data')[0].id;
  displayId = null;

  // Populate saved displays modal
  loadSavedDisplays(datasetId);

  $('#load_saved_plots').show();
  $('#plot_type_container').show();
  $('#options_spinner').show();
  $('#advanced_options_container').show();
  $('#gene_container').show();

  // Get genes for this dataset
  geneSymbols = await fetchGeneSymbols({ datasetId, undefined });
  createGeneDropdown(geneSymbols);
  loadGeneCarts();
  $('#genes_not_found').hide();

  // Get categorical observations for this dataset
  const data = await fetchH5adInfo({ datasetId, undefined });
  obsLevels = curateObservations(data.obs_levels);
  createObsGroupbyField(obsLevels);
  createObsDropdowns(obsLevels);
  createVolcanoDropdowns(obsLevels);

  // Ensure genes and observation columns dropdown toooltip shows
  $(function () {
    $('[data-toggle="tooltip"]').tooltip();
  });

  $('#options_spinner').hide();
  $('#update_plot').show();
  $('#reset_obs').show();
});

// Load user's gene carts
$('#selected_gene_cart').change(function () {
  const geneCartId = $('#selected_gene_cart').select2('data')[0].id;
  const params = { session_id: session_id, gene_cart_id: geneCartId };
  const d = new $.Deferred(); // Causes editable to wait until results are returned
  // User is not logged in
  if (!session_id) {
    d.resolve();
  } else {
    // User is logged in
    $('#search_gene_symbol').prop('disabled', true);
    $('#selected_gene_cart_loading_c').show();

    // Get the gene cart members and populate the gene symbol search bar
    $.ajax({
      url: './cgi/get_gene_cart_members.cgi',
      type: 'post',
      data: params,
      success: function (data, newValue, oldValue) {
        if (data.success === 1) {
          // Append gene symbols to search bar
          const geneCartSymbols = [];

          // format gene symbols into search string
          $.each(data.gene_symbols, function (i, item) {
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
            $('#genes_not_found').text('The following gene cart genes were not found in this dataset: ' + differenceString);
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

// Cannot groupby if observations are going to be clustered
$('#cluster_cols').change(function () {
  if (this.checked) {
    $('.js-obs-groupby').prop('checked', false);
    $('#obs_groupby_container').hide();
  } else {
    $('#obs_groupby_container').show();
  }
});

// Some options are specific to certain plot types
$('#plot_type_select').change(function () {
  switch ($('#plot_type_select').val()) {
    case 'heatmap':
      $('#cluster_cols_checkbox_container').show();
      $('#obs_groupby_container').show();
      $('#volcano_options_container').hide();
      $('adjusted_pvals_checkbox_container').hide();
      break;
    case 'mg_violin':
      $('#cluster_cols').prop('checked', false);
      $('#cluster_cols_checkbox_container').hide();
      $('#obs_groupby_container').show();
      $('#volcano_options_container').hide();
      $('adjusted_pvals_checkbox_container').hide();
      break;
    default:
      // volcano
      $('#cluster_cols').prop('checked', false);
      $('#cluster_cols_checkbox_container').hide();
      $('.js-obs-groupby').prop('checked', false);
      $('#obs_groupby_container').hide();
      $('#volcano_options_container').show();
      $('adjusted_pvals_checkbox_container').show();
  }
});

$(document).on('click', '#update_plot', async function () {
  // Remove supplementary plot and reset its genes filter
  if (supplementaryGenesFilters.length) {
    supplementaryGenesFilters = [];
    $('#supplementary_plot').remove();
  }

  // Render dataset plot HTML
  const plotTemplate = $.templates('#dataset_plot_tmpl');
  const plotHtml = plotTemplate.render({ dataset_id: datasetId });
  $('#dataset_plot').html(plotHtml);

  const plotType = $('#plot_type_select').select2('data')[0].id;

  const groupbyFilter = $('input[name="obs_groupby"]:checked').val();

  // Update filters based on selection
  obsFilters = {};
  for (const property in obsLevels) {
    const propData = $(`#${property}_dropdown`).select2('data');
    obsFilters[property] = propData.map(function (elem) {
      return elem.id;
    });

    // If no groups for an observation are selected, delete filter
    if (!obsFilters[property].length) {
      delete obsFilters[property];
    }
  }

  if (!plotType) {
    window.alert('Please select a plot type.');
    return;
  }

  if (!(plotType === 'volcano' || Object.keys(obsFilters).length)) {
    window.alert('At least one observation must have categories filtered.');
    return;
  }

  genesFilters = $('#gene_dropdown').select2('data').map(function (elem) {
    return elem.id;
  });

  if (!(plotType === 'volcano' || genesFilters.length)) {
    window.alert('At least one gene must be provided.');
    return;
  }

  const clusterCols = $('#cluster_cols').is(':checked');
  const adjustPvals = $('#adj_pvals').is(':checked');

  const condition1 = $('#volcano_condition1').select2('data')[0].id;
  const condition2 = $('#volcano_condition2').select2('data')[0].id;

  if (plotType === 'volcano' && !(condition1 && condition2)) {
    window.alert('Both comparision conditions must be selected to generate a volcano plot.');
    return;
  }

  // Validation related to the conditions
  if (condition1 && condition2) {
    const conditionKey = condition1.split(';-;')[0];
    if (condition1.split(';-;')[0] !== condition2.split(';-;')[0]) {
      window.alert('Please choose 2 conditions from the same observation group.');
      return;
    }

    if (condition1.split(';-;')[1] === condition2.split(';-;')[1]) {
      window.alert('Please choose 2 different conditions.');
      return;
    }

    // If condition category was filtered, the selected groups must be present
    if (conditionKey in obsFilters) {
      if (!(obsFilters[conditionKey].includes(condition1.split(';-;')[1]) &&
        obsFilters[conditionKey].includes(condition2.split(';-;')[1]))) {
        window.alert('Condition observation is found in filters list, but one or both condition groups is filtered out. Please adjust.');
        return;
      }
    }
  }

  plotConfig = {
    groupby_filter: groupbyFilter,
    plot_type: plotType,
    gene_symbols: genesFilters,
    obs_filters: obsFilters,
    cluster_cols: clusterCols,
    adj_pvals: adjustPvals,
    condition1: condition1,
    condition2: condition2
  };

  // Draw the updated chart
  $('#dataset_spinner').show();
  await draw(datasetId, plotConfig);
  $('#dataset_spinner').hide();

  const saveTemplate = $.templates('#save_plot_tmpl');
  const saveHtml = saveTemplate.render({ plot_type: plotType });
  $('#save_functions').html(saveHtml);
});

// If "all" button is clicked, populate dropdown with all groups in this observation
$(document).on('click', '.all', function () {
  const id = this.id;
  const group = id.replace('_all', '');

  $(`#${group}_dropdown`).val(obsLevels[group]);
  $(`#${group}_dropdown`).trigger('change'); // This actually triggers select2 to show the dropdown vals
});

// If gene is clicked in plot, display supplementary violin plot
$(document).on('click', 'g.y5tick text a', async function () {
  const gene = $(this).text();

  // Add or remove gene depending on if it is already in array
  const index = supplementaryGenesFilters.indexOf(gene);
  if (index === -1) {
    supplementaryGenesFilters.push(gene);
    $(this).parent().css('fill', 'crimson');
  } else {
    supplementaryGenesFilters.splice(index, 1);
    $(this).parent().css('fill', 'rgb(42, 63, 95)'); // original default fill color
  }

  // Render supplementary plot HTML
  if (supplementaryGenesFilters.length) {
    // Render supplementary plot HTML
    const plotTemplate = $.templates('#supplementary_plot_tmpl');
    const plotHtml = plotTemplate.render({ dataset_id: datasetId });
    $('#supplementary_plot').html(plotHtml);

    // Draw the supplementary chart
    const groupbyFilter = $('input[name="obs_groupby"]:checked').val();
    plotConfig = {
      groupby_filter: groupbyFilter,
      plot_type: 'mg_violin',
      gene_symbols: supplementaryGenesFilters,
      obs_filters: obsFilters
    };
    $('#supplementary_spinner').show();
    await draw(datasetId, plotConfig, true);
    $('#supplementary_spinner').hide();
  }
});

// Reset observation filters choices to be empty
$(document).on('click', '#reset_obs', async function () {
  // Get categorical observations for this dataset
  const data = await fetchH5adInfo({ datasetId, undefined });
  obsLevels = curateObservations(data.obs_levels);

  // Update fields dependent on dataset observations
  createObsGroupbyField(obsLevels);
  createObsDropdowns(obsLevels);
  createVolcanoDropdowns(obsLevels);
});

// If advanced options collapsable is clicked, toggle arrow b/t up and down
$(document).on('click', '#advanced_options_button', function () {
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
    id: displayId,
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

  if (res && res.success) {
    $('#saved_plot_confirmation').text('Plot successfully saved');
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
$(document).on('click', '.js-edit-display', async function () {
  const id = this.id;
  displayId = id.replace('_edit', '');

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
  if (res && res.success) {
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
  const displayId = id.replace('_save_default', '');

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

  $('#saved_default_confirmation').show();
  if (res && res.success) {
    $('#saved_default_confirmation_text').text('Display successfully saved as your new default.');
    $('#saved_default_confirmation').addClass('alert-success');
  } else {
    $('#saved_default_confirmation_text').text('There was an issue setting the default display.');
    $('#saved_default_confirmation').addClass('alert-danger');
  }
});
