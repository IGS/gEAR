// NOTE: - SAdkins - 11/3/21 - Refactored code using P42 VSCode extension https://p42.ai/documentation/code-action/

let plot_data = null;
let selected_data = null;

let condition_x = null;
let condition_y = null;

// SAdkins - 2/15/21 - This is a list of datasets already log10-transformed where if selected will use log10 as the default dropdown option
// This is meant to be a short-term solution until more people specify their data is transformed via the metadata
const log10_transformed_datasets = [
"320ca057-0119-4f32-8397-7761ea084ed1"
, "df726e89-b7ac-d798-83bf-2bd69d7f3b52"
, "bad48d04-db27-26bc-2324-e88506f751fd"
, "dbd715bf-778a-4923-6fe7-c587987cdb00"
, "c8d99d13-394f-a87f-5d3a-395968fdb619"
, "bee735e5-d180-332c-7892-dd751dd76bb8"
, "c4f16a12-9e98-47be-4335-b8321282919e"
, "17a07bf4-b41a-d9c3-9aa7-b4729390f57a"
, "6a0a2bca-0f86-59d0-4e3d-4457be3a71ff"
, "39e01b71-415f-afa7-0c64-f0e996be0fb7"
, "6482c608-a6bd-d8b1-6bc1-5b53c34ed61c"
, "0c5a4c18-c2a9-930c-6e52-ef411f54eb67"
, "3c02d449-61ab-4bcd-f100-5f5937b1794e"
, "23e3797f-3016-8142-cbe8-69b03131ad95"
, "b16eeb8d-d68e-c7c9-9dc9-a3f4821e9192"
, "b96f448a-315d-549d-6e8a-83cdf1ce1b5c"
, "b0420910-a0fa-e920-152d-420b6275d3af"
, "f1ce4e63-3577-8020-8307-e88f1fb98953"
, "2f79f784-f7f7-7dc3-9b3e-4c87a4346d91"
, "c32835d3-cac4-bb0e-a90a-0b41dec6617a"
, "fbe1296e-572c-d388-b9d1-6e2a6bf10b0a"
, "1b12dde9-1762-7564-8fbd-1b07b750505f"
, "a2dd9f06-5223-0779-8dfc-8dce7a3897e1"
, "f7de7db2-b4cb-ebe3-7f1f-b278f46f1a7f"
, "e34fa5c6-1083-cacb-eedf-23f59f2e005f"
, "0c5fb6b0-31ab-6bfc-075d-76756ccd56b4"
, "a183b2e6-ab38-458a-52a6-5eb014d073da"
, "c4f16a12-9e98-47be-4335-b8321282919e"
, "2a25e445-2776-8913-076f-9a147a43e8b4"
, "2786d849-f11c-2de6-b22e-12c940aafe07"
, "2e3423b3-74db-d436-8357-abb3031d47e9"
, "4cb2ac62-c283-86a9-83cb-2c1b381948f2"
, "d0659d69-1a33-8b84-252c-f7ded46aa3d6"
, "cee5325d-434f-fefe-d2e6-e0be39421951"
, "34f8f131-8158-db83-7df9-db9003797dff"
, "7ddb4965-e710-faf7-ee26-4ce95d7602a8"
, "f122cac5-c79f-8ea2-166e-42415916db11"
, "173ab634-a2b1-87bc-f1ef-d288de0bcd1a"
];

// TODO: Have mechanism to convert non-categorical column to categorical if it was erroneously added as numerical

const dataset_tree = new DatasetTree({treeDiv: '#dataset_tree'});

window.onload = () => {
  // check if the user is already logged in
  check_for_login();
  session_id = Cookies.get("gear_session_id");

  $(".btn-apply-filter").on("click", () => {
    $(".initial_instructions").hide();
    $("#myChart").html("");
    $("#error_loading_c").hide();
    $("#plot_loading").show();
    load_comparison_graph();
  });

  /***** gene cart stuff *****/
  $("#create_gene_cart").on("click", () => {
    $("#create_gene_cart_dialog").show("fade");
  });
  $("#cancel_save_gene_cart").on("click", () => {
    $("#create_gene_cart_dialog").hide("fade");
    $("#gene_cart_name").val("");
  });
  $("#gene_cart_name").on("input", function () {
    if ($(this).val() == "") {
      $("#save_gene_cart").prop("disabled", true);
    } else {
      $("#save_gene_cart").prop("disabled", false);
    }
  });
  $("#save_gene_cart").on("click", () => {
    $("#save_gene_cart").prop("disabled", true);

    if (CURRENT_USER) {
      save_gene_cart();
    } else {
      alert("You must be signed in to do that.");
    }
  });
  /***** end gene cart stuff *****/

  $("#dataset_id").on("change", () => {
      populate_condition_selection_control();
      // Change the default if the dataset is already log10 transformed
      if (log10_transformed_datasets.includes($('#dataset_id').val())) {
        $('#log_base').val('10');
      } else {
        $('#log_base').val('2');
      }
    });

  $("#statistical_test").on("change", () => {
    if ($("#statistical_test").val()) {
      $("#test_pval_cutoff").prop("disabled", false);
    } else {
      $("#test_pval_cutoff").prop("disabled", true);
    }
  });

  // Create observer to watch if user changes (ie. successful login does not refresh page)
  // See: https://developer.mozilla.org/en-US/docs/Web/API/MutationObserver

  // Select the node that will be observed for mutations
  const target_node = document.getElementById('loggedin_controls');
  // Create an observer instance linked to the callback function
  const observer = new MutationObserver(populate_dataset_selection_controls);
  // For the "config" settings, do not monitor the subtree of nodes as that will trigger the callback multiple times.
  // Just seeing #loggedin_controls go from hidden (not logged in) to shown (logged in) is enough to trigger.
  observer.observe(target_node, { attributes: true });
};

function download_selected_genes() {
  // Builds a file in memory for the user to download.  Completely client-side.
  // plot_data contains three keys: x, y and symbols
  // build the file string from this

  const x_label = $('#x_label').val().length ? $('#x_label').val() : "x-condition";
  const y_label = $('#y_label').val().length ? $('#y_label').val() : "y-condition";

  const file_contents =
      $("#log_base").val() == "raw"
    ? "gene_symbol\tp-value\tfold change\t"
  + x_label + "\t"
  + y_label + "\n"
    : "gene_symbol\tp-value\tfold change\t"
  + x_label + " (log" + $("#log_base").val() +")\t"
  + y_label + " (log" + $("#log_base").val() +")\n";

  selected_data.points.forEach((pt) => {
    // Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
    file_contents +=
      `${pt.data.id[pt.pointNumber]}\t`
      + ($("#statistical_test").val() ? pt.data.pvals[pt.pointNumber].toExponential(2) : "NA") + "\t"
      + pt.data.foldchange[pt.pointNumber].toFixed(1) + "\t"
      + pt.x.toFixed(1) + "\t"
      + pt.y.toFixed(1) + "\n";
  });

  const element = document.createElement("a");
  element.setAttribute(
    "href",
    `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(file_contents)}`
  );
  element.setAttribute("download", "selected_genes.tsv");
  element.style.display = "none";
  document.body.appendChild(element);
  element.click();
  document.body.removeChild(element);
}

// Build list of selected categories and groups
function update_selected_conditions() {
  const condition = {};
  // Create current object of which groups are to be included (checked)
  $('#conditions_accordion').find('.js-group-check').each(function(){
    const id = $(this).data("group");

    const category = id.split(';-;')[0];
    const group = id.split(';-;')[1];

    if (Object.keys(condition).indexOf(category) === -1) {
      condition[category] = [];
    }

    if ($(this).prop("checked")) {
      condition[category].push(group);
    }
  });
  return condition;
}

function load_comparison_graph() {
  // Save current state of active condition tab
  if ($("#condition_x_tab").hasClass("active")) {
    condition_x = update_selected_conditions();
  } else {
    // on #condition_y_tab
    condition_y = update_selected_conditions();
  }

  const dataset_id = $("#dataset_id").val();
  const dataset_text = $("#dataset_id").text();
  const sanitized_condition_x = sanitize_condition(condition_x);
  const sanitized_condition_y = sanitize_condition(condition_y);
  const condition_x_string = JSON.stringify(sanitized_condition_x);
  const condition_y_string = JSON.stringify(sanitized_condition_y);

  // empty error message, so that user/helper won't get confused
  $("#ticket_error_msg").empty();
  $("#error_loading_c").hide();
  $("#genes_not_found").empty().hide();

  $.ajax({
    url: "./cgi/get_dataset_comparison.cgi",
    type: "POST",
    data: {
      dataset_id,
      condition_x: condition_x_string,
      condition_y: condition_y_string,
      fold_change_cutoff: $("#fold_change_cutoff").val(),
      std_dev_num_cutoff: $("#std_dev_num_cutoff").val(),
      log_transformation: $("#log_base").val(),
      statistical_test: $("#statistical_test").val(),
    },
    dataType: "json",
    success(data, textStatus, jqXHR) {
      if (data.success == 1) {
        $("#fold_change_std_dev").html(data.fold_change_std_dev);
        plot_data = data;
        plot_data_to_graph(data);
      } else {
        // Handle graphing failures
        $("#plot_loading").hide();
        $("#ticket_dataset_id").text(dataset_id);
        $("#ticket_dataset_text").text(dataset_text);
        $("#ticket_datasetx_condition").text(condition_x_string);
        $("#ticket_datasety_condition").text(condition_y_string);
        $("#ticket_error_msg").html(data.error);
        $("#error_loading_c").show();
      }
    },
    error(jqXHR, textStatus, errorThrown) {
      // Handle graphing failures
      $("#plot_loading").hide();
      $("#ticket_dataset_id").text(dataset_id);
      $("#ticket_dataset_text").text(dataset_text);
      $("#ticket_datasetx_condition").text(condition_x_string);
      $("#ticket_datasety_condition").text(condition_y_string);
      $("#error_loading_c").show();
    },
  });
}

async function fetch_h5ad_observations (dataset_id) {
  const base = `./api/h5ad/${dataset_id}`;
  const { data } = await axios.get(base);
  return data;
}

async function populate_condition_selection_control() {
  const dataset_id = $("#dataset_id").val();
  $("#conditions_accordion").html("<p>Loading ... </p>");
  const obs_data = await fetch_h5ad_observations(dataset_id);
  const cat_obs = obs_data.obs_levels;  // cat->groups
  const all_obs = obs_data.obs_columns; // Array
  const noncat_obs = Object.values(all_obs).filter(x => Object.keys(cat_obs).indexOf(x) === -1);  // Array

  // Render templates
  const selector_tmpl = $.templates("#dataset_condition_options");
  const selector_html = selector_tmpl.render(cat_obs);
  $("#conditions_accordion").html(selector_html);

  if (noncat_obs.length) {
    const noncat_tmpl = $.templates("#non_categories_list");
    const noncat_html = noncat_tmpl.render({noncat_obs});
    $("#noncats").html(noncat_html);
  }

  // Since we want the default state of the category groups to be unchecked,
  // we do not set an initial state of stored checked conditions

  if (obs_data.has_replicates == 1) {
    $("#statistical_test_label").html("");
    $("#statistical_test").attr("disabled", false);
  } else {
    $("#statistical_test_label").html(
      "Not applicable since this dataset has no replicates"
    );
    $("#statistical_test").attr("disabled", true);
  }
}

$('#condition_x_tab').on('shown.bs.tab', (e) => {

  condition_y = update_selected_conditions();

  // Load condition_x stuff
  $('.js-cat-check').prop("checked", false);
  $('.js-group-check').prop("checked", false);
  for (const cat in condition_x) {
    for (const elem of condition_x[cat]) {
      $(`input[data-group="${cat};-;${elem}"]`).prop("checked", true);
    }
    $('.js-group-check').change();  // trigger so the cat checkbox matches up
  }

})

$('#condition_y_tab').on('shown.bs.tab', (e) => {

  condition_x = update_selected_conditions();

  // Load condition_y stuff
  $('.js-cat-check').prop("checked", false);
  $('.js-group-check').prop("checked", false);
  for (const cat in condition_y) {
    for (const elem of condition_y[cat]) {
      $(`input[data-group="${cat};-;${elem}"]`).prop("checked", true);
    }
    $('.js-group-check').change();  // trigger so the cat checkbox matches up
  }
})

$(document).on('change', '.js-cat-check', function (e) {
  // If turned on, check all group boxes
  // If turned off, uncheck all group boxes
  const checked = $(this).prop("checked");

  const id = this.id;
  const category = id.replace('_check', '');
  const escapedCategory = $.escapeSelector(category);
  const category_collaspable = $(`#${escapedCategory}_body`);

  category_collaspable.find('input[type="checkbox"]').prop({checked});

  // Expand collaspable since category was focused on
  category_collaspable.collapse('show');

  // Update the selected conditons div and the "axis label" input boxes
  const template = $.templates("#selected_conditions_list");

  if ($("#condition_x_tab").hasClass("active")) {
    const curr_condition_x = update_selected_conditions();
    const sanitized_condition_x = sanitize_condition(curr_condition_x);
    const htmlOutput = template.render(sanitized_condition_x);
    $('#selected_x_condition').html(htmlOutput);
    $('#x_label').val(JSON.stringify(sanitized_condition_x).replace("{", "").replace("}", ""));

  } else {
    // on #condition_y_tab
    const curr_condition_y = update_selected_conditions();
    const sanitized_condition_y = sanitize_condition(curr_condition_y);
    const htmlOutput = template.render(sanitized_condition_y);
    $('#selected_y_condition').html(htmlOutput);
    $('#y_label').val(JSON.stringify(sanitized_condition_y).replace("{", "").replace("}", ""));
  }
})

$(document).on('click', '.js-cat-collapse', function (e) {
  // If category was clicked, then toggle collapsable element
  // Controlling via JS instead of "data-target" since we may need to escape CSS selectors
  const id = this.id;
  const category = id.replace('_collapse', '');
  const escapedCategory = $.escapeSelector(category);
  const category_collaspable = $(`#${escapedCategory}_body`);
  category_collaspable.collapse('toggle');
})

$(document).on('change', '.js-group-check', function(e) {
  // https://css-tricks.com/indeterminate-checkboxes/
  // After changing checkbox status, check siblings
  // and determine if category checkbox should be
  // checked, not checked, or indeterminate

  const checked = $(this).prop("checked");

  // Get category name out of the checkbox ID
  const id = $(this).data("group");
  const category = id.split(';-;')[0]
  const escapedCategory = $.escapeSelector(category);
  const category_header = $(`#${escapedCategory}_check`);
  const category_collaspable = $(`#${escapedCategory}_body`);

  // Get checked status of all other checkboxes in this category
  // If there is a combination of checked/unchecked the "each" loop breaks early
  let all = true;
  $(category_collaspable).find('input[type="checkbox"]').each(function(){
    const return_value = all = ($(this).prop("checked") === checked);
    return return_value;
  });

  if (all) {
    // All group checkboxes are the same as the category checkbox
    category_header.prop({
      "indeterminate": false,
      "checked": checked
    });
  } else {
    // All group checkbox states are mixed.  Category checkbox is indeterminate and unchecked
    category_header.prop({
      "indeterminate": true,
      "checked": false
    });
  }

  // Update the selected conditons div and the "axis label" input boxes
  const template = $.templates("#selected_conditions_list");
  if ($("#condition_x_tab").hasClass("active")) {
    const curr_condition_x = update_selected_conditions();
    const sanitized_condition_x = sanitize_condition(curr_condition_x);
    const htmlOutput = template.render(sanitized_condition_x);
    $('#selected_x_condition').html(htmlOutput);
    $('#x_label').val(JSON.stringify(sanitized_condition_x).replace("{", "").replace("}", ""));

  } else {
    // on #condition_y_tab
    const curr_condition_y = update_selected_conditions();
    const sanitized_condition_y = sanitize_condition(curr_condition_y);
    const htmlOutput = template.render(sanitized_condition_y);
    $('#selected_y_condition').html(htmlOutput);
    $('#y_label').val(JSON.stringify(sanitized_condition_y).replace("{", "").replace("}", ""));
  }
})

function sanitize_condition(condition) {
  const sanitized_condition = {}
  for (const property in condition) {
    // If no groups for an observation are selected, delete filter
    if (condition[property].length) {
      sanitized_condition[property] = condition[property];
    }
  }
  return sanitized_condition;
}

async function populate_dataset_selection_controls() {
  const dataset_id = getUrlParameter("dataset_id");
  $('#pre_dataset_spinner').show();
  await $.ajax({
    type: "POST",
    url: "./cgi/get_h5ad_dataset_list.cgi",
    data: {
      session_id: CURRENT_USER.session_id,
      for_page: "compare_dataset",
      include_dataset_id: dataset_id,
    },
    dataType: "json",
    success(data) {
      let counter = 0
      // Populate select box with dataset information owned by the user
      const user_datasets = [];
      if (data.user.datasets.length > 0) {
        // User has some profiles
        $.each(data.user.datasets, (_i, item) => {
          user_datasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
        });
      }
      // Next, add datasets shared with the user
      const shared_datasets = [];
      if (data.shared_with_user.datasets.length > 0) {
        // User has some profiles
        $.each(data.shared_with_user.datasets, (_i, item) => {
          shared_datasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id  });
        });
      }
      // Now, add public datasets
      const domain_datasets = [];
      if (data.public.datasets.length > 0) {
        // User has some profiles
        $.each(data.public.datasets, (_i, item) => {
          domain_datasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id  });
        });
      }

      dataset_tree.userDatasets = user_datasets;
      dataset_tree.sharedDatasets = shared_datasets;
      dataset_tree.domainDatasets = domain_datasets;
      dataset_tree.generateTree();

      // was there a requested dataset ID already?
      if (dataset_id !== undefined) {
        $("#dataset_id").val(dataset_id);
        $('#dataset_id').text(dataset_tree.treeData.find(e => e.dataset_id === dataset_id).text);
        $("#dataset_id").trigger("change");
      }
    },
    error(xhr, status, msg) {
      report_error(`Failed to load dataset list because msg: ${msg}`);
    },
  });
  $('#pre_dataset_spinner').hide();
}

function plot_data_to_graph(data) {
  $("#tbl_selected_genes").hide();
  $("#selection_methods_c").show();

  const point_labels = [];
  let perform_ranking = false;

  if ($("#statistical_test").val()) {
    perform_ranking = true;
  }

  let plotdata = null;

  if (!perform_ranking) {
    for (i = 0; i < data.symbols.length; i++) {
      point_labels.push(`Gene symbol: ${data.symbols[i]}`);
    }

    plotdata = [
      {
        id: data.symbols,
        pvals: data.pvals_adj,
        x: data.x,
        y: data.y,
        foldchange: data.fold_changes,
        mode: "markers",
        type: "scatter",
        text: point_labels,
        marker: {
          color: "#2F103E",
          size: 4,
        },
      },
    ];
  } else {
    const pval_cutoff = parseFloat($("#test_pval_cutoff").val());
    const passing = { x: [], y: [], labels: [], id: [], pvals: [], foldchange: []};
    const failing = { x: [], y: [], labels: [], id: [], pvals: [], foldchange: []};

    for (i = 0; i < data.x.length; i++) {
      // pvals_adj array consist of 1-element arrays, so let's flatten to prevent potential issues
      // Caused by rank_genes_groups output (1 inner array per query comparison group)
      data.pvals_adj = data.pvals_adj.flat();

      const this_pval = parseFloat(data.pvals_adj[i]);

      if ((this_pval <= pval_cutoff)) {
        // good scoring match
        passing.x.push(data.x[i]);
        passing.y.push(data.y[i]);
        passing.foldchange.push(data.fold_changes[i]);
        passing.labels.push(
          "Gene symbol: " +
            data.symbols[i] +
            "   P-value: " +
            this_pval.toPrecision(6)
        );
        passing.id.push(data.symbols[i]);
        passing.pvals.push(data.pvals_adj[i]);
      } else {
        // this one didn't pass the p-value cutoff
        failing.x.push(data.x[i]);
        failing.y.push(data.y[i]);
        failing.foldchange.push(data.fold_changes[i]);
        failing.labels.push(
          "Gene symbol: " +
            data.symbols[i] +
            "   P-value: " +
            this_pval.toPrecision(6)
        );
        failing.id.push(data.symbols[i]);
        failing.pvals.push(data.pvals_adj[i]);
      }
    }

    plotdata = $("input[name='stat_action']:checked").val() == "colorize" ? [
      {
        id: passing.id,
        pvals: passing.pvals,
        x: passing.x,
        y: passing.y,
        foldchange: passing.foldchange,
        mode: "markers",
        name: "Passed cutoff",
        type: "scatter",
        text: passing.labels,
        marker: {
          color: "#FF0000",
          size: 4,
        },
      },
      {
        id: failing.id,
        pvals: failing.pvals,
        x: failing.x,
        y: failing.y,
        foldchange: failing.foldchange,
        mode: "markers",
        name: "Did not pass cutoff",
        type: "scatter",
        text: failing.labels,
        marker: {
          color: "#A1A1A1",
          size: 4,
        },
      },
    ] : [
      {
        id: passing.id,
        pvals: passing.pvals,
        x: passing.x,
        y: passing.y,
        foldchange: passing.foldchange,
        mode: "markers",
        type: "scatter",
        text: passing.labels,
        marker: {
          color: "#2F103E",
          size: 4,
        },
      },
    ];
  }

  const layout = {
    title: $("#dataset_id").text(),
    xaxis: {
      title: $('#x_label').val().length ? $('#x_label').val() : JSON.stringify(data.condition_x_idx),
      type: "",
    },
    yaxis: {
      title: $('#y_label').val().length ? $('#y_label').val() : JSON.stringify(data.condition_y_idx),
      type: "",
    },
    annotations: [],
    margin: { t: 40 },
    hovermode: "closest",
    dragmode: "select",
  };

  // Take genes to search for and highlight their datapoint in the plot
  const genes_not_found = [];
  if ($('#highlighted_genes').val()) {
    const searched_genes = $('#highlighted_genes').val().replace(/\s/g, "").split(",");
    searched_genes.forEach((gene) => {
      let found = false;
      plots:
      for (i = 0; i < plotdata.length; i++) {
        genes:
        for (j = 0; j < plotdata[i].id.length; j++) {
          if (gene.toLowerCase() === plotdata[i].id[j].toLowerCase() ) {
            // If gene is found add an annotation arrow
            layout.annotations.push({
              xref: "x",
              yref: "y",
              x: plotdata[i].x[j],
              y: plotdata[i].y[j],
              text:plotdata[i].id[j],
              font: {
                color: "crimson",
              },
              showarrow: true,
              arrowcolor: "crimson",

            });
            found = true;
            break plots;
          }
        }
      }
      if (! found)
        genes_not_found.push(gene);
    });
  }


  $("#plot_loading").hide();
  const graphDiv = document.getElementById("myChart");
  Plotly.newPlot(graphDiv, plotdata, layout, { showLink: false });
  $("#selected_label").hide();
  $("#controls_label").show();

  // If searched-for genes were not found, display under plot
  if (genes_not_found.length) {
    const genes_not_found_str = genes_not_found.join(", ");
    $("#genes_not_found").text(`Searched genes not found: ${genes_not_found_str}`);
    $("#genes_not_found").show();
  } else {
    $("#genes_not_found").hide();
  }

  // If plot data is selected, create the right-column table and do other misc things
  graphDiv.on("plotly_selected", (eventData) => {
    selected_data = eventData;
    selected_gene_data = [];

    eventData.points.forEach((pt) => {
      // Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
      // Each trace has its own "pointNumber" ids so gene symbols and pvalues needed to be passed in for each plotdata trace
      selected_gene_data.push({
        gene_symbol: pt.data.id[pt.pointNumber],
        pvals: $("#statistical_test").val() ? pt.data.pvals[pt.pointNumber].toExponential(2) : "NA",
        foldchange: pt.data.foldchange[pt.pointNumber].toFixed(1),
      });
    });

    // Sort by adjusted p-value in descending order either by fold change or p-values
    selected_gene_data.sort((a, b) => b.foldchange - a.foldchange);
    if ($("#statistical_test").val())
      selected_gene_data.sort((a, b) => a.pvals - b.pvals);

    const template = $.templates("#selected_genes_tmpl");
    const htmlOutput = template.render(selected_gene_data);
    $("#selected_genes_c").html(htmlOutput);

    // Highlight table rows that match searched genes
    if ($('#highlighted_genes').val()) {
      const searched_genes = $('#highlighted_genes').val().replace(/\s/g, "").split(",");
      // Select the first column (gene_symbols) in each row
      $("#selected_genes_c tr td:first-child").each(function() {
        const table_gene = $(this).text();
        searched_genes.forEach((gene) => {
          if (gene.toLowerCase() === table_gene.toLowerCase() ) {
            $(this).parent().addClass("table-success");
          }
        });
      })
    }

    // toggle visibilities
    $(".selection_instructions").hide();
    $("#saved_gene_cart_info_c").hide();
    $("#tbl_selected_genes").show();

    $("#controls_label").hide();
    $("#selected_label").show();

    if ($("#log_base").val() == "raw") {
      $("#tbl_selected_genes_transformation_row").hide();
    } else {
      $("#table_transformation_label").text(`Log${$("#log_base").val()}`);
      $("#tbl_selected_genes_transformation_row").show();
    }
  });

  window.onresize = () => {
    Plotly.Plots.resize(graphDiv);
  };
}

function save_gene_cart() {
  // must have access to USER_SESSION_ID
  const gc = new GeneCart({
      session_id: CURRENT_USER.session_id,
      label: $("#gene_cart_name").val(),
      gctype: 'unweighted-list',
      organism_id: $("#dataset_").data('organism-id'),
      is_public: 0
  });

  selected_data.points.forEach((pt) => {
    const gene = new Gene({
      id: plot_data.gene_ids[pt.pointNumber],
      gene_symbol: plot_data.symbols[pt.pointNumber],
    });
    gc.add_gene(gene);
  });

  gc.save(update_ui_after_gene_cart_save_success, update_ui_after_gene_cart_save_failure);
}

// Sort selected gene table (using already generated table data)
// Taken from https://www.w3schools.com/howto/howto_js_sort_table.asp
function sortTable(n) {
  let table;
  let rows;
  let switching;
  let i;
  let x;
  let y;
  let shouldSwitch;
  let dir;
  let switchcount = 0;
  table = document.getElementById("tbl_selected_genes");

  switching = true;
  // Set the sorting direction to ascending:
  dir = "asc";
  /* Make a loop that will continue until
    no switching has been done: */
  while (switching) {
    // Start by saying: no switching is done:
    switching = false;
    rows = table.rows;
    /* Loop through all table rows (except the
      first, which contains table headers): */
    for (i = 1; i < rows.length - 1; i++) {
      // Start by saying there should be no switching:
      shouldSwitch = false;
      /* Get the two elements you want to compare,
        one from current row and one from the next: */
      x = rows[i].getElementsByTagName("td")[n];
      y = rows[i + 1].getElementsByTagName("td")[n];
      /* Check if the two rows should switch place,
        based on the direction, asc or desc: */
      if (dir == "asc") {
        // First column is gene_symbol... rest are numbers
        if (n === 0 && x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
          // If so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
        if (Number(x.innerHTML) > Number(y.innerHTML)) {
          shouldSwitch = true;
          break;
        }
      } else if (dir == "desc") {
        if (n === 0 && x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
          // If so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
        if (Number(x.innerHTML) < Number(y.innerHTML)) {
          shouldSwitch = true;
          break;
        }
      }
    }
    if (shouldSwitch) {
      /* If a switch has been marked, make the switch
        and mark that a switch has been done: */
      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
      switching = true;
      // Each time a switch is done, increase this count by 1:
      switchcount++;
    } else {
      /* If no switching has been done AND the direction is "asc",
        set the direction to "desc" and run the while loop again. */
      if (switchcount == 0 && dir == "asc") {
        dir = "desc";
        switching = true;
      }
    }
  }
}

function update_ui_after_gene_cart_save_success(gc) {
  $("#create_gene_cart_dialog").hide("fade");
  $("#saved_gene_cart_info_c > h3").html(`Cart: ${gc.label}`);
  $("#gene_cart_member_count").html(gc.genes.length);
  $("#saved_gene_cart_info_c").show();
}

function update_ui_after_gene_cart_save_failure(gc) {
  $("#create_gene_cart_dialog").hide("fade");
  $("#saved_gene_cart_info_c > h3").html("There was an issue saving the gene cart.");
  $("#saved_gene_cart_info_c").show();
}
