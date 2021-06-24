"use strict";

var obsFilters = {};
var genesFilters = [];
var supplementaryGenesFilters = [];

var datasetId = null;
var obsLevels = null;
var geneSymbols = null;    //TODO: get encoded string from POST and decode to array to pre-populate


// Async to ensure data is fetched before proceeding
(async () => {

  // check if the user is already logged in
  await check_for_login();
  session_id = Cookies.get("gear_session_id");

  $(function () {
    $('[data-toggle="tooltip"]').tooltip()
  })

  $('#dataset_select').select2({
    placeholder: 'To search, click to select or start typing a dataset name',
  });
  await populateDatasets();

  // Initialize plot types
  $('#plot_type_select').select2({
    placeholder: 'Choose how to plot',
    allowClear: true,
  });

  // Hide further configs until a dataset is chosen
  $('#plot_type_container').hide();
  $('#advanced_options_container').hide();
  $('#gene_container').hide();

  $("#dataset_select").change(async function () {
    datasetId = $("#dataset_select").select2('data')[0].id;

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
    obsLevels = curateObservations(data['obs_levels']);
    createObsGroupbyField(obsLevels);
    createObsDropdowns(obsLevels);

    // Ensure genes and observation columns dropdown toooltip shows
    $(function () {
      $('[data-toggle="tooltip"]').tooltip()
    })

    $('#options_spinner').hide();
    $("#update_plot").prop("disabled", false);
    $("#reset_obs").prop("disabled", false);

  })

  $(document).on('click', "#update_plot", function () {

    // Remove supplementary plot and reset its genes filter
    if (supplementaryGenesFilters.length) {
      supplementaryGenesFilters = []
      $('#supplementary_plot').remove()
    }

    // Render dataset plot HTML
    const plotTemplate = $.templates("#dataset_plot_tmpl");
    const plotHtml = plotTemplate.render({ dataset_id: datasetId });
    $('#dataset_plot').html(plotHtml);

    $('#dataset_spinner').show();

    var plotType = $("#plot_type_select").select2('data')[0].id;

    var groupbyFilter = $('input[name="obs_groupby"]:checked').val();

    // Update filters based on selection
    obsFilters = {};
    for (var property in obsLevels) {
      var propData = $(`#${property}_dropdown`).select2('data');
      obsFilters[property] = propData.map(function (elem) {
        return elem.id;
      });

      // If no groups for an observation are selected, delete filter
      if (!obsFilters[property].length) {
        delete obsFilters[property];
      }
    }

    if (!Object.keys(obsFilters).length) {
      alert("At least one observation must have categories filtered.");
      return;
    }

    genesFilters = $("#gene_dropdown").select2('data').map(function (elem) {
      return elem.id;
    })

    if (!genesFilters.length) {
      alert("At least one gene must be provided.");
      return;
    }

    const clusterCols = $('#cluster_cols').is(':checked');
    const payload = { clusterCols, groupbyFilter };

    // Draw the updated chart
    draw(datasetId, plotType, genesFilters, obsFilters, payload);
    $('#dataset_spinner').hide();

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
      case "heatmap":
        $('#obs_checkbox_container').show();
        $('#obs_groupby_container').show();
        break;
      case "mg_violin":
        $('#cluster_cols').prop('checked', false);
        $('#obs_checkbox_container').hide();
        $('#obs_groupby_container').show();
        break;
      default:
        $('#cluster_cols').prop('checked', false);
        $('#obs_checkbox_container').hide();
        $('.js-obs-groupby').prop('checked', false);
        $('#obs_groupby_container').hide();
    }

  });

  // If "all" button is clicked, populate dropdown with all groups in this observation
  $(document).on('click', ".all", function () {
    var id = this.id;
    var group = id.replace("_all", "");

    $(`#${group}_dropdown`).val(obsLevels[group]);
    $(`#${group}_dropdown`).trigger('change');  // This actually triggers select2 to show the dropdown vals
  });

  // If gene is clicked in plot display supplementary violin plot
  $(document).on('click', "g.y5tick text a", function () {
    var gene = $(this).text();

    // Add or remove gene depending on if it is already in array
    const index = supplementaryGenesFilters.indexOf(gene);
    if (index === -1) {
      supplementaryGenesFilters.push(gene);
      $(this).parent().css("fill", "crimson");
    } else {
      supplementaryGenesFilters.splice(index, 1);
      $(this).parent().css("fill", "rgb(42, 63, 95)"); // original default fill color
    }

    // Render supplementary plot HTML
    if (supplementaryGenesFilters.length) {
      // Render supplementary plot HTML
      const plotTemplate = $.templates("#supplementary_plot_tmpl");
      const plotHtml = plotTemplate.render({ dataset_id: datasetId });
      $('#supplementary_plot').html(plotHtml);
      $('#supplementary_spinner').show();

      // Draw the supplementary chart
      var groupbyFilter = $('input[name="obs_groupby"]:checked').val();
      var payload = { groupbyFilter }
      draw(datasetId, "mg_violin", supplementaryGenesFilters, obsFilters, payload, true);
      $('#supplementary_spinner').hide();

    }
  });

  $(document).on('click', "#reset_obs", async function () {
    // Get categorical observations for this dataset
    const data = await fetchH5adInfo({ datasetId, undefined });
    obsLevels = curateObservations(data.obsLevels);
    createObsGroupbyField(obsLevels);
    createObsDropdowns(obsLevels);

    // Update list of observations to groupby by
    createObsGroupbyField(obsLevels);
    createObsDropdowns(obsLevels);
  });

  // If advanced options collapsable is clicked, toggle arrow b/t up and down
  $(document).on('click', "#advanced_options_button", function () {
    if ($("#advanced_options_arrow").hasClass("fa-caret-down")) {
      $("#advanced_options_arrow").removeClass("fa-caret-down");
      $("#advanced_options_arrow").addClass("fa-caret-up");
    } else {
      $("#advanced_options_arrow").addClass("fa-caret-down");
      $("#advanced_options_arrow").removeClass("fa-caret-up");
    }
  });

  // Load user's gene carts
  $('#selected_gene_cart').change(function () {
    session_id = Cookies.get('gear_session_id');
    var geneCartId = $("#selected_gene_cart").select2('data')[0].id;
    var params = { "session_id": session_id, "gene_cart_id": geneCartId };
    var d = new $.Deferred(); //Causes editable to wait until results are returned
    //User is not logged in
    if (!session_id) {
      d.resolve();
    } else {
      //User is logged in
      $("#search_gene_symbol").prop("disabled", true);
      $("#selected_gene_cart_loading_c").show();

      //Get the gene cart members and populate the gene symbol search bar
      $.ajax({
        url: './cgi/get_gene_cart_members.cgi',
        type: 'post',
        data: params,
        success: function (data, newValue, oldValue) {
          if (data['success'] == 1) {
            // Append gene symbols to search bar
            var geneCartSymbols = [];

            //format gene symbols into search string
            $.each(data['gene_symbols'], function (i, item) {
              geneCartSymbols.push(item['label']);
            });

            let geneCartSymbolsLowerCase = geneCartSymbols.map(x => x.toLowerCase());
            let geneSymbolsLowerCase = geneSymbols.map(x => x.toLowerCase());

            // Get genes from gene cart that are present in dataset's genes.  Preserve casing of dataset's genes.
            let intersection = geneSymbols.filter(x => geneCartSymbolsLowerCase.includes(x.toLowerCase()));
            let uniqIntersection = [... new Set(intersection)];  // Dataset may have a gene repeated in it, so resolve this.

            // Get genes from gene cart that are not in dataset's genes
            let difference = geneCartSymbols.filter(x => !geneSymbolsLowerCase.includes(x.toLowerCase()));

            $('#gene_dropdown').val(uniqIntersection);
            $('#gene_dropdown').trigger("change");

            if (difference.length > 0) {
              let differenceString = difference.join(", ");
              $("#genes_not_found").text("The following gene cart genes were not found in this dataset: " + differenceString);
              $("#genes_not_found").show();
            } else {
              $("#genes_not_found").hide();
            }
          }
          d.resolve();

        }
      });
    }
    return d.promise();
  });

})();


async function getData(datasetId, plotType, genesFilters, obsFilters, payload) {
  return await axios.post(`/api/plot/${datasetId}/mg_dash`, {
    ...payload,
    plot_type: plotType,
    gene_symbols: genesFilters,
    obs_filters: obsFilters,
  });
}

async function fetchGeneSymbols(payload) {
  const { datasetId, analysis } = payload;
  const base = `./api/h5ad/${datasetId}/genes`;
  const query = analysis ? `?analysis=${analysis.id}` : "";

  const { data } = await axios.get(`${base}${query}`);
  return data['gene_symbols'];
}

async function fetchH5adInfo(payload) {
  const { datasetId, analysis } = payload;
  const base = `./api/h5ad/${datasetId}`;
  const query = analysis ? `?analysis=${analysis.id}` : "";
  const { data } = await axios.get(`${base}${query}`);
  return data;
}

async function populateDatasets() {
  $.ajax({
    type: "POST",
    url: "./cgi/get_h5ad_dataset_list.cgi",
    data: {
      session_id: undefined,  // TODO: define with user login info
    },
    dataType: "json",
    success: function (data) {
      // Populate select box with dataset information owned by the user
      if (data["user"]["datasets"].length > 0) {
        var userDatasetListTmpl = $.templates("#dataset_list_tmpl");
        var userDatasetListHtml = userDatasetListTmpl.render(
          data["user"]["datasets"]
        );
        $("#dataset_ids_user").html(userDatasetListHtml);
      } else {
        $("#dataset_id .user-initial").html("Not logged in");
      }

      // Next, add datasets shared with the user
      if (data["shared_with_user"]["datasets"].length > 0) {
        var sharedWithUserDatasetListTmpl = $.templates(
          "#dataset_list_tmpl"
        );
        var sharedWithUserDatasetListHtml = sharedWithUserDatasetListTmpl.render(
          data["shared_with_user"]["datasets"]
        );
        $("#dataset_ids_shared_with_user").html(
          sharedWithUserDatasetListHtml
        );
      }

      // Now, add public datasets
      if (data["public"]["datasets"].length > 0) {
        var publicDatasetListTmpl = $.templates("#dataset_list_tmpl");
        var publicDatasetListHtml = publicDatasetListTmpl.render(
          data["public"]["datasets"]
        );
        $("#dataset_ids_public").html(publicDatasetListHtml);
      }

    },
    error: function (xhr, status, msg) {
      console.error("Failed to load dataset list because msg: " + msg);
    },
  });
}

// Draw plotly chart in HTML
function drawChart(data, datasetId, supplementary = false) {
  const targetDiv = supplementary ? `dataset_${datasetId}_secondary` : `dataset_${datasetId}_h5ad`;
  const {'plot_json': plotJson, 'plot_config': plotConfig} = data;

  var layoutMods = {
    height: targetDiv.clientHeight,
    width: targetDiv.clientWidth,
  };

  // Overwrite plot layout and config values with custom ones from display
  var layout = {
    ...plotJson.layout,
    ...layoutMods,
  };

  var configMods = {
    responsive: false,
  };

  const config = {
    ...plotConfig,
    ...configMods,
  };
  Plotly.newPlot(targetDiv, plotJson.data, layout, config);
}

// Submit API request and draw the HTML
async function draw(datasetId, plotType, genesFilters, obsFilters, payload, supplementary = false) {
  const {
    data
  } = await getData(datasetId, plotType, genesFilters, obsFilters, payload);
  drawChart(data, datasetId, supplementary);
}

function createGeneDropdown(genes) {
  var tmpl = $.templates("#gene_dropdown_tmpl");
  var data = { genes: genes }
  var html = tmpl.render(data);
  $("#gene_dropdown_container").html(html);
  $('#gene_dropdown').select2({
    placeholder: 'To search, click to select or start typing some gene names',
    allowClear: true
  });
}

function createObsGroupbyField(obsLevels) {
  var tmpl = $.templates("#obs_groupby_tmpl");  // Get compiled template using jQuery selector for the script block
  var html = tmpl.render(obsLevels);                 // Render template using data - as HTML string
  $("#obs_groupby_container").html(html);       // Insert HTML string into DOM
}

function createObsDropdowns(obsLevels) {
  var tmpl = $.templates("#obs_dropdowns_tmpl");
  var html = tmpl.render(obsLevels);
  $("#obs_dropdowns_container").html(html);
  $('select.js-obs-levels').select2({
    placeholder: 'Start typing to filter categories. Click "All" to use all categories',
    allowClear: true
  });
}

function curateObservations(obsLevels) {
  // Delete useless filters
  for (var property in obsLevels) {
    if (property === "color" || property.endsWith("colors")) {
      delete obsLevels[property];
    } else if (obsLevels[property].length === 1) {
      delete obsLevels[property];
    }
  }
  return obsLevels
}

function loadGeneCarts() {
  var d = new $.Deferred();
  session_id = Cookies.get('gear_session_id');

  if (!session_id) {
    //User is not logged in. Hide gene carts container
    $("#gene_cart_container").hide();
    d.resolve();
  } else {
    $.ajax({
      url: './cgi/get_user_gene_carts.cgi',
      type: 'post',
      data: { 'session_id': session_id },
      dataType: 'json',
      success: function (data, textStatus, jqXHR) { //source https://stackoverflow.com/a/20915207/2900840
        var userGeneCarts = [];

        if (data['gene_carts'].length > 0) {
          //User has some profiles
          $.each(data['gene_carts'], function (i, item) {
            userGeneCarts.push({ value: item['id'], text: item['label'] });

          });
          $("#gene_cart_container").show();
          var tmpl = $.templates("#gene_cart_tmpl");
          var html = tmpl.render(userGeneCarts);
          $("#selected_gene_cart").html(html);
          // Add blank option so select2 will show placeholder
          $("#selected_gene_cart").prepend("<option value=''></option>").val('');
          $('#selected_gene_cart').select2({
            placeholder: 'Select a preloaded set of genes.',
            allowClear: true
          });
        } else {
          $("#gene_cart_container").hide();
        }
        d.resolve();
      },
      error: function (jqXHR, textStatus, errorThrown) {
        //display_error_bar(jqXHR.status + ' ' + errorThrown.name);
        d.fail();
      }
    });
  }
  d.promise();
}
