// TODO
//  - Make existing plots disappear right when the user hits 'Plot' to redraw

plot_data = null;
selected_data = null;

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
];

window.onload = function () {
  // check if the user is already logged in
  check_for_login();
  session_id = Cookies.get("gear_session_id");

  $(".btn-apply-filter").on("click", function () {
    $(".initial_instructions").hide();
    $("#myChart").html("");
    $("#error_loading_c").hide();
    $("#plot_loading").show();
    load_comparison_graph();
  });

  /***** gene cart stuff *****/
  $("#create_gene_cart").on("click", function () {
    $("#create_gene_cart_dialog").show("fade");
  });
  $("#cancel_save_gene_cart").on("click", function () {
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
  $("#save_gene_cart").on("click", function () {
    $("#save_gene_cart").prop("disabled", true);

    if (CURRENT_USER) {
      save_gene_cart();
    } else {
      alert("You must be signed in to do that.");
    }
  });
  /***** end gene cart stuff *****/

  $("#dataset_id").on("change", function () {
    populate_condition_selection_control();
    // Change the default if the dataset is already log10 transformed
    if (log10_transformed_datasets.includes($('#dataset_id').val())) {
        $('#log_base').val('10');
    } else {
        $('#log_base').val('2');
    }
  });

  $("#statistical_test").on("change", function () {
    if ($("#statistical_test").val()) {
      $("#test_pval_cutoff").prop("disabled", false);
    } else {
      $("#test_pval_cutoff").prop("disabled", true);
    }
  });

  // initially disable the condition selectors
  $("#dataset1_conditions").attr("disabled", "disabled");
  $("#dataset2_conditions").attr("disabled", "disabled");
};

function download_selected_genes() {
  // Builds a file in memory for the user to download.  Completely client-side.
  // plot_data contains three keys: x, y and symbols
  // build the file string from this
  if ($("#log_base").val() == "raw") {
    file_contents =
      "gene_symbol\tp-value\tfold change\t"
      + $("#dataset1_conditions").val() + "\t"
      + $("#dataset2_conditions").val() + "\n";
  } else {
    file_contents =
      "gene_symbol\tp-value\tfold change\t"
      + $("#dataset1_conditions").val() + " (log" + $("#log_base").val() +")\t"
      + $("#dataset2_conditions").val() + " (log" + $("#log_base").val() +")\n"
  }

  selected_data.points.forEach(function (pt) {
    // Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
    file_contents +=
      pt.data.id[pt.pointNumber] + "\t"
      + ($("#statistical_test").val() ? pt.data.pvals[pt.pointNumber].toExponential(2) : "NA") + "\t"
      + pt.data.foldchange[pt.pointNumber].toFixed(1) + "\t"
      + pt.x.toFixed(1) + "\t"
      + pt.y.toFixed(1) + "\n";
  });

  var element = document.createElement("a");
  element.setAttribute(
    "href",
    "data:text/tab-separated-values;charset=utf-8," + encodeURIComponent(file_contents)
  );
  element.setAttribute("download", "selected_genes.tsv");
  element.style.display = "none";
  document.body.appendChild(element);
  element.click();
  document.body.removeChild(element);
}

function load_comparison_graph() {
  $.ajax({
    url: "./cgi/get_dataset_comparison.cgi",
    type: "POST",
    data: {
      dataset1_id: $("#dataset_id").val(),
      dataset1_condition: $("#dataset1_conditions").val(),
      dataset2_condition: $("#dataset2_conditions").val(),
      fold_change_cutoff: $("#fold_change_cutoff").val(),
      std_dev_num_cutoff: $("#std_dev_num_cutoff").val(),
      log_transformation: $("#log_base").val(),
      statistical_test: $("#statistical_test").val(),
    },
    dataType: "json",
    success: function (data, textStatus, jqXHR) {
      if (data["success"] == 1) {
        $("#fold_change_std_dev").html(data["fold_change_std_dev"]);
        plot_data = data;
        plot_data_to_graph(data);
      } else {
        console.log("CGI unsuccessful");
        // Handle graphing failures
        $("#plot_loading").hide();
        $("#ticket_datasetx_id").text($("#dataset_id").val());
        $("#ticket_datasetx_condition").text($("#dataset1_conditions").val());
        $("#ticket_datasety_id").text($("#dataset_id").val());
        $("#ticket_datasety_condition").text($("#dataset2_conditions").val());
        $("#ticket_error_msg").text(data["error"]);
        $("#error_loading_c").show();
      }
    },
    error: function (jqXHR, textStatus, errorThrown) {
      // Handle graphing failures
      $("#plot_loading").hide();
      $("#ticket_datasetx_id").text($("#dataset_id").val());
      $("#ticket_datasetx_condition").text($("#dataset1_conditions").val());
      $("#ticket_datasety_id").text($("#dataset_id").val());
      $("#ticket_datasety_condition").text($("#dataset2_conditions").val());
      $("#error_loading_c").show();
    },
  });
}

function populate_condition_selection_control() {
  dataset_id = $("#dataset_id").val();
  $("#dataset1_conditions").attr("disabled", "disabled");
  $("#dataset2_conditions").attr("disabled", "disabled");
  $("#dataset1_conditions").html("<option>Loading ... </option>");
  $("#dataset2_conditions").html("<option>Loading ... </option>");

  $.ajax({
    url: "./cgi/get_condition_list.cgi",
    type: "POST",
    data: { dataset_id: dataset_id },
    dataType: "json",
    success: function (data, textStatus, jqXHR) {
      if (data["success"] === 0) {
        display_error_bar(data["error"]);
      } else if (data["success"] === -1) {
        $("#dataset1_conditions > option")
          .html("This dataset is not ready to compare.")
          .attr("selected", "selected");
      } else {
        var selectorTmpl = $.templates("#dataset_condition_options");
        var selectorHtml = selectorTmpl.render(data["conditions"]);
        $("#dataset1_conditions").html(selectorHtml);
        $("#dataset2_conditions").html(selectorHtml);

        $("#dataset1_conditions").removeAttr("disabled");
        $("#dataset2_conditions").removeAttr("disabled");

        if (data["has_replicates"] == 1) {
          $("#statistical_test_label").html("");
          $("#statistical_test").attr("disabled", false);
        } else {
          $("#statistical_test_label").html(
            "Not applicable since this dataset has no replicates"
          );
          $("#statistical_test").attr("disabled", true);
        }
      }
    },
    error: function (jqXHR, textStatus, errorThrown) {
      //alert("Failure!  Status: (" + textStatus + ") Error: (" + errorThrown + ")");
      console.log("textStatus= ", textStatus);
      console.log("errorThrown= ", errorThrown);
      display_error_bar(jqXHR.status + " " + errorThrown.name);
    },
  });
}

function populate_dataset_selection_controls() {
  $.ajax({
    type: "POST",
    url: "./cgi/get_h5ad_dataset_list.cgi",
    data: {
      session_id: CURRENT_USER.session_id,
      for_page: "compare_dataset",
      include_dataset_id: getUrlParameter("dataset_id"),
    },
    dataType: "json",
    success: function (data) {
      if (data["user"]["datasets"].length > 0) {
        var user_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
        var user_dataset_list_html = user_dataset_list_tmpl.render(
          data["user"]["datasets"]
        );
        $("#dataset_ids_user").html(user_dataset_list_html);
      } else {
        $("#dataset_id .user_initial").html("Not logged in");
      }

      if (data["shared_with_user"]["datasets"].length > 0) {
        var shared_with_user_dataset_list_tmpl = $.templates(
          "#dataset_list_tmpl"
        );
        var shared_with_user_dataset_list_html = shared_with_user_dataset_list_tmpl.render(
          data["shared_with_user"]["datasets"]
        );
        $("#dataset_ids_shared_with_user").html(
          shared_with_user_dataset_list_html
        );
      }

      if (data["public"]["datasets"].length > 0) {
        var public_dataset_list_tmpl = $.templates("#dataset_list_tmpl");
        var public_dataset_list_html = public_dataset_list_tmpl.render(
          data["public"]["datasets"]
        );
        $("#dataset_ids_public").html(public_dataset_list_html);
      }

      // was there a requested dataset ID already?
      var dataset_id = getUrlParameter("dataset_id");
      if (dataset_id !== undefined) {
        $("#dataset_id").val(dataset_id);
        $("#dataset_id").trigger("change");
      }
    },
    error: function (xhr, status, msg) {
      report_error("Failed to load dataset list because msg: " + msg);
    },
  });
}

function plot_data_to_graph(data) {
  $("#tbl_selected_genes").hide();
  $("#selection_methods_c").show();

  var point_labels = [];
  var perform_ranking = false;

  if ($("#statistical_test").val()) {
    perform_ranking = true;
  }

  var plotdata = null;

  if (perform_ranking) {
    var pval_cutoff = parseFloat($("#test_pval_cutoff").val());
    var passing = { x: [], y: [], labels: [], id: [], pvals: [], foldchange: []};
    var failing = { x: [], y: [], labels: [], id: [], pvals: [], foldchange: []};

    for (i = 0; i < data["x"].length; i++) {
      // pvals_adj array consist of 1-element arrays, so let's flatten to prevent potential issues
      // Probably happened when pulling from AnnData object in get_dataset_comparison.cgi
      data["pvals_adj"] = data["pvals_adj"].flat();

      var this_pval = parseFloat(data["pvals_adj"][i]);

      if (this_pval <= pval_cutoff) {
        // good scoring match
        passing["x"].push(data["x"][i]);
        passing["y"].push(data["y"][i]);
        passing["foldchange"].push(data["fold_changes"][i]);
        passing["labels"].push(
          "Gene symbol: " +
            data["symbols"][i] +
            "   P-value: " +
            this_pval.toPrecision(6)
        );
        passing["id"].push(data["symbols"][i]);
        passing["pvals"].push(data["pvals_adj"][i]);
      } else {
        // this one didn't pass the p-value cutoff
        failing["x"].push(data["x"][i]);
        failing["y"].push(data["y"][i]);
        failing["foldchange"].push(data["fold_changes"][i]);
        failing["labels"].push(
          "Gene symbol: " +
            data["symbols"][i] +
            "   P-value: " +
            this_pval.toPrecision(6)
        );
        failing["id"].push(data["symbols"][i]);
        failing["pvals"].push(data["pvals_adj"][i]);
      }
    }

    if ($("input[name='stat_action']:checked").val() == "colorize") {
      plotdata = [
        {
          id: passing["id"],
          pvals: passing["pvals"],
          x: passing["x"],
          y: passing["y"],
          foldchange: passing["foldchange"],
          mode: "markers",
          name: "Passed cutoff",
          type: "scatter",
          text: passing["labels"],
          marker: {
            color: "#FF0000",
            size: 4,
          },
        },
        {
          id: failing["id"],
          pvals: failing["pvals"],
          x: failing["x"],
          y: failing["y"],
          foldchange: failing["foldchange"],
          mode: "markers",
          name: "Did not pass cutoff",
          type: "scatter",
          text: failing["labels"],
          marker: {
            color: "#A1A1A1",
            size: 4,
          },
        },
      ];
    } else {
      plotdata = [
        {
          id: passing["id"],
          pvals: passing["pvals"],
          x: passing["x"],
          y: passing["y"],
          foldchange: passing["foldchange"],
          mode: "markers",
          type: "scatter",
          text: passing["labels"],
          marker: {
            color: "#2F103E",
            size: 4,
          },
        },
      ];
    }
  } else {
    for (i = 0; i < data["symbols"].length; i++) {
      point_labels.push("Gene symbol: " + data["symbols"][i]);
    }

    plotdata = [
      {
        id: data["symbols"],
        pvals: data["pvals_adj"],
        x: data["x"],
        y: data["y"],
        foldchange: data["fold_changes"],
        mode: "markers",
        type: "scatter",
        text: point_labels,
        marker: {
          color: "#2F103E",
          size: 4,
        },
      },
    ];
  }

  var layout = {
    title: $("#dataset_id option:selected").text(),
    xaxis: {
      title: $("#dataset1_conditions option:selected").text(),
      type: "",
    },
    yaxis: {
      title: $("#dataset2_conditions option:selected").text(),
      type: "",
    },
    annotations: [],
    margin: { t: 40 },
    hovermode: "closest",
    dragmode: "select",
  };

  // Take genes to search for and highlight their datapoint in the plot
  var genes_not_found = [];
  if ($('#highlighted_genes').val()) {
    var searched_genes = $('#highlighted_genes').val().replace(/\s/g, "").split(",");
    searched_genes.forEach((gene) => {
      var found = false;
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
  var graphDiv = document.getElementById("myChart");
  Plotly.newPlot(graphDiv, plotdata, layout, { showLink: false });
  $("#selected_label").hide();
  $("#controls_label").show();

  // If searched-for genes were not found, display under plot
  if (genes_not_found.length) {
    var genes_not_found_str = genes_not_found.join(", ");
    $("#genes_not_found").text("Searched genes not found: " + genes_not_found_str);
    $("#genes_not_found").show();
  } else {
    $("#genes_not_found").hide();
  }

  // If plot data is selected, create the right-column table and do other misc things
  graphDiv.on("plotly_selected", function (eventData) {
    selected_data = eventData;
    selected_gene_data = [];

    eventData.points.forEach(function (pt) {
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

    var template = $.templates("#selected_genes_tmpl");
    var htmlOutput = template.render(selected_gene_data);
    $("#selected_genes_c").html(htmlOutput);

    // Highlight table rows that match searched genes
    if ($('#highlighted_genes').val()) {
      var searched_genes = $('#highlighted_genes').val().replace(/\s/g, "").split(",");
      // Select the first column (gene_symbols) in each row
      $("#selected_genes_c tr td:first-child").each(function() {
        var table_gene = $(this).text();
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
      $("#table_transformation_label").text("Log" + $("#log_base").val());
      $("#tbl_selected_genes_transformation_row").show();
    }
  });

  window.onresize = function () {
    Plotly.Plots.resize(graphDiv);
  };
}

function save_gene_cart() {
  // must have access to USER_SESSION_ID
  var gc = new GeneCart({
    session_id: CURRENT_USER.session_id,
    label: $("#gene_cart_name").val(),
  });

  selected_data.points.forEach(function (pt) {
    var gene = new Gene({
      id: plot_data["gene_ids"][pt.pointNumber],
      gene_symbol: plot_data["symbols"][pt.pointNumber],
    });
    gc.add_gene(gene);
  });

  gc.save(update_ui_after_gene_cart_save);
}

// Sort selected gene table (using already generated table data)
// Taken from https://www.w3schools.com/howto/howto_js_sort_table.asp
function sortTable(n) {
  var table,
    rows,
    switching,
    i,
    x,
    y,
    shouldSwitch,
    dir,
    switchcount = 0;
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
        if (n === 0) {
          if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
            // If so, mark as a switch and break the loop:
            shouldSwitch = true;
            break;
          }
        }
        if (Number(x.innerHTML) > Number(y.innerHTML)) {
          shouldSwitch = true;
          break;
        }
      } else if (dir == "desc") {
        if (n === 0) {
          if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
            // If so, mark as a switch and break the loop:
            shouldSwitch = true;
            break;
          }
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

function update_ui_after_gene_cart_save(gc) {
  $("#create_gene_cart_dialog").hide("fade");
  $("#saved_gene_cart_info_c > h3").html("Cart: " + gc.label);
  $("#gene_cart_member_count").html(gc.genes.length);
  $("#saved_gene_cart_info_c").show();
}
