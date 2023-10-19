'use strict';

// SAdkins - 2/15/21 - This is a list of datasets already log10-transformed where if selected will use log10 as the default dropdown option
// This is meant to be a short-term solution until more people specify their data is transformed via the metadata
const LOG10_TRANSFORMED_DATASETS = [
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
, '48bab518-439e-4a17-b868-6b225abf2c73'
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
, "80eadbe6-49ac-8eaf-f2fb-e07706cf117b"
];

let sessionId;
let colorblindMode;
let facetWidget;
let datasetId;
let organismId;	// Used for saving as gene cart
let compareData;;
let selectedData;

const datasetTree = new DatasetTree({
    element: document.getElementById("dataset_tree")
    , searchElement: document.getElementById("dataset_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current_dataset_c").classList.remove("is-hidden");
        document.getElementById("current_dataset").textContent = e.node.title;
        document.getElementById("current_dataset_post").textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        datasetId = newDatasetId;

        // Click to get to next step
        document.getElementById("condition_compare_s").click();

        // Clear "success/failure" icons
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.classList.add("is-hidden");
        }
		for (const elt of document.getElementsByClassName("js-step-failure")) {
			elt.classList.add("is-hidden");
		}

		const compareSeriesElt = document.getElementById("compare_series");
		compareSeriesElt.parentElement.classList.add("is-loading");

		// Create facet widget, which will refresh filters
		facetWidget = await createFacetWidget(sessionId, datasetId, null, {});
		document.getElementById("facet_content").classList.remove("is-hidden");
		document.getElementById("selected_facets").classList.remove("is-hidden");

		// Update compare series options
		const catColumns = facetWidget.aggregations.map((agg) => agg.name);
		updateSeriesOptions("js-compare", catColumns);

		compareSeriesElt.parentElement.classList.remove("is-loading");

    })
});

const createFacetWidget = async (sessionId, datasetId, analysisId, filters) => {
    document.getElementById("selected_facets_loader").classList.remove("is-hidden")

    const {aggregations, total_count:totalCount} = await fetchAggregations(sessionId, datasetId, analysisId, filters);
    document.getElementById("num_selected").textContent = totalCount;


    const facetWidget = new FacetWidget({
        aggregations,
        filters,
        onFilterChange: async (filters) => {
            if (filters) {
                try {
                    const {aggregations, total_count:totalCount} = await fetchAggregations(sessionId, datasetId, analysisId, filters);
                    facetWidget.updateAggregations(aggregations);
                    document.getElementById("num_selected").textContent = totalCount;
                } catch (error) {
                    logErrorInConsole(error);
                }
            } else {
                // Save an extra API call
                facetWidget.updateAggregations(facetWidget.aggregations);
            }
        },
		filterHeaderExtraClasses:"has-background-white"
    });
    document.getElementById("selected_facets_loader").classList.add("is-hidden")
    return facetWidget;
}

/* Creates a Toast-style message in the upper-corner of the screen. */
const createToast = (msg, levelClass="is-danger") => {
    const template = `
    <div class="notification js-toast ${levelClass} animate__animated animate__fadeInUp animate__faster">
        <button class="delete"></button>
        ${msg}
    </div>
    `
    const html = generateElements(template);

    const numToasts = document.querySelectorAll(".js-toast.notification").length;

    if (document.querySelector(".js-toast.notification")) {
        // If .js-toast notifications are present, append under final notification
        // This is to prevent overlapping toast notifications
        document.querySelector(".js-toast.notification:last-of-type").insertAdjacentElement("afterend", html);
        // Position new toast under previous toast with CSS
        html.style.setProperty("top", `${(numToasts * 70) + 30}px`);
    } else {
        // Otherwise prepend to top of main content
        document.getElementById("main_c").prepend(html);
    }

    // This should get the newly added notification since it is now the first
    document.querySelector(".js-toast.notification .delete").addEventListener("click", (event) => {
        const notification = event.target.closet(".js-toast.notification");
        notification.remove(notification);
    });

    // For a success message, remove it after 3 seconds
    if (levelClass === "is-success") {
        const notification = document.querySelector(".js-toast.notification:last-of-type");
        notification.classList.remove("animate__fadeInUp");
        notification.classList.remove("animate__faster");
        notification.classList.add("animate__fadeOutDown");
        notification.classList.add("animate__slower");
    }
}

const fetchAggregations = async (session_id, dataset_id, analysis_id, filters) => {
    const payload = {session_id, dataset_id, analysis_id, filters};
    try {
        const {data} = await axios.post(`/api/h5ad/${dataset_id}/aggregations`, payload);
        if (data.hasOwnProperty("success") && data.success < 1) {
            throw new Error(data?.message || "Could not fetch number of observations for this dataset. Please contact the gEAR team.");
        }
        const {aggregations, total_count} = data;
        return {aggregations, total_count};
    } catch (error) {
        logErrorInConsole(error);
    }
}

const fetchDatasetComparison = async (dataset_id, obs_filters, condition_x, condition_y, fold_change_cutoff, std_dev_num_cutoff, log_transformation, statistical_test) => {
	const payload = {dataset_id, obs_filters, condition_x, condition_y, fold_change_cutoff, std_dev_num_cutoff, log_transformation, statistical_test};
	try {
		const {data} = await axios.post("cgi/get_dataset_comparison.cgi", convertToFormData(payload));
		return data;
	} catch (error) {
		logErrorInConsole(error);
		const msg = "Could not fetch dataset comparison. Please contact the gEAR team."
		createToast(msg);
		throw new Error(msg);
	}
}

const fetchDatasets = async (session_id) => {
    const payload = {session_id}
    try {
        const {data} = await axios.post("cgi/get_h5ad_dataset_list.cgi", convertToFormData(payload));
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch datasets. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

const getComparisons = async (event) => {

	// set loading icon
	event.target.classList.add("is-loading");

	const filters = JSON.stringify(facetWidget.filters);

	const compareSeries = document.getElementById("compare_series").value
	const xSeries = document.getElementById("compare_x_series").value
	const ySeries = document.getElementById("compare_y_series").value

	const xCondition = `${compareSeries};-;${xSeries}`;
	const yCondition = `${compareSeries};-;${ySeries}`;

	const foldChangeCutoff = document.getElementById("fc_cutoff").value;
	const stdDevNumCutoff = document.getElementById("standard_deviation").value;
	const logTransformation = document.getElementById("log_base").value;
	const statisticalTest = document.getElementById("statistical_test").value;

	try {
		const data = await fetchDatasetComparison(datasetId, filters, xCondition, yCondition, foldChangeCutoff, stdDevNumCutoff, logTransformation, statisticalTest);
		if (data?.success < 1) {
			throw new Error(data?.message || "Could not fetch dataset comparison. Please contact the gEAR team.");
		}
		compareData = data;
		plotDataToGraph(compareData);

		// Hide this view
		document.getElementById("content_c").classList.add("is-hidden");
		// Generate and display "post-plotting" view/container
		document.getElementById("post_plot_content_c").classList.remove("is-hidden");

	} catch (error) {
		console.error(error);
		handleGetComparisonError(datasetId, xSeries, ySeries);
	} finally {
		event.target.classList.remove("is-loading");
	}

	// When a plot configuration ID is selected, populate the plot configuration post textbox
	const plotConfigElts = ["statistical_test", "pval_cutoff", "cutoff_filter_action", "log_base", "fc_cutoff", "standard_deviation"];
	for (const elt of plotConfigElts) {
		// if value is empty, set to "None", or if disabled, set to "N/A"
		let value = document.getElementById(elt).disabled ? "N/A" : document.getElementById(elt).value || "None"

		// Append extra flavor text
		if (elt == "log_base" && !(value === "raw" )) {
			value = `log${value}`
		}

		if (elt == "standard_deviation") {
			value = `±${value}`
		}

		document.getElementById(`${elt}_post`).textContent = value;
	}

}

const getSeriesItems = (series) => {
	return facetWidget.aggregations.find((agg) => agg.name === series).items;
}

const getSeriesNames = (seriesItems) => {
	return seriesItems.map((item) => item.name);
}

const handleGetComparisonError = (datasetID, conditionX, conditionY) => {
	const msg = `Could not fetch dataset comparison. Please contact the gEAR team.`;
	createToast(msg);
	console.error(msg);
}

/* Transform and load dataset data into a "tree" format */
const loadDatasetTree = async () => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const datasetData = await fetchDatasets(sessionId);

        let counter = 0;

        // Populate select box with dataset information owned by the user
        if (datasetData.user.datasets.length > 0) {
            // User has some profiles
            for (const item of datasetData.user.datasets) {
                if (item) {
                    userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        // Next, add datasets shared with the user
        if (datasetData.shared_with_user.datasets.length > 0) {
            for (const item of datasetData.shared_with_user.datasets) {
                if (item) {
                    sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        // Now, add public datasets
        if (datasetData.public.datasets.length > 0) {
            for (const item of datasetData.public.datasets) {
                if (item) {
                    domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
                }
            };
        }
        datasetTree.userDatasets = userDatasets;
        datasetTree.sharedDatasets = sharedDatasets;
        datasetTree.domainDatasets = domainDatasets;
        datasetTree.generateTree();
    } catch (error) {
        document.getElementById("dataset_s_failed").classList.remove("is-hidden");
    }

}

const plotDataToGraph = (data) => {

	const statisticalTest = document.getElementById("statistical_test").value;

	const pointLabels = [];
	const performRanking = statisticalTest ? true : false;

	const plotData = [];

	if (performRanking) {
		const pValCutoff = document.getElementById("pval_cutoff").value;
		const pvalCutoff = parseFloat(pValCutoff);
		const passing = { x: [], y: [], labels: [], id: [], pvals: [], foldchange: []};
		const failing = { x: [], y: [], labels: [], id: [], pvals: [], foldchange: []};

		data.x.forEach((trace, i) => {
			// pvals_adj array consist of 1-element arrays, so let's flatten to prevent potential issues
			// Caused by rank_genes_groups output (1 inner array per query comparison group)
			data.pvals_adj = data.pvals_adj.flat();

			const thisPval = parseFloat(data.pvals_adj[i]);

			const arrayToPushInto = (thisPval <= pvalCutoff) ? passing : failing;

			arrayToPushInto.x.push(trace);
			arrayToPushInto.y.push(data.y[i]);
			arrayToPushInto.foldchange.push(data.fold_changes[i]);
			arrayToPushInto.labels.push(
			"Gene symbol: " +
				data.symbols[i] +
				"   P-value: " +
				thisPval.toPrecision(6)
			);
			arrayToPushInto.id.push(data.symbols[i]);
			arrayToPushInto.pvals.push(data.pvals_adj[i]);

		});

		const passColor = CURRENT_USER.colorblind_mode ? 'rgb(0, 34, 78)' : "#FF0000";
		const failColor = CURRENT_USER.colorblind_mode ? 'rgb(254, 232, 56)' : "#A1A1A1";

		const statAction = document.getElementById("cutoff_filter_action").value;
		if (statAction === "colorize") {
			const passingObj = {
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
					color: passColor,
					size: 4,
					},
				}
			const failingObj = {
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
					color: failColor,
					size: 4,
					},
				}
			plotData.push(passingObj);
			plotData.push(failingObj);
		} else {

			const passingObj = {
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
			}
			plotData.push(passingObj);
		}

	} else {
		for (const gene of data.symbols) {
			pointLabels.push(`Gene symbol: ${gene}`);
		}

		const dataObj = {
				id: data.symbols,
				pvals: data.pvals_adj,
				x: data.x,
				y: data.y,
				foldchange: data.fold_changes,
				mode: "markers",
				type: "scatter",
				text: pointLabels,
				marker: {
				color: "#2F103E",
				size: 4,
				},
			}
		plotData.push(dataObj);
	}

	const layout = {
		title: "Dataset Comparison",
		xaxis: {
			title: data.condition_x
		},
		yaxis: {
			title: data.condition_y,
		},
		annotations: [],
		hovermode: "closest",
		dragmode: "select",
	};

	/*
	// Take genes to search for and highlight their datapoint in the plot

	const annotationColor = CURRENT_USER.colorblind_mode ? 'rgb(125, 124, 118)' : "crimson";

	const highlightedGenes = document.getElementById("highlighted_genes").value;
	const genesNotFound = [];
	if (highlightedGenes) {
		const searchedGenes = highlightedGenes.replace(/,?\s/g, ",").split(",");
		searchedGenes.forEach((gene) => {
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
						color: annotationColor,
					},
					showarrow: true,
					arrowcolor: annotationColor,

					});
					found = true;
					break plots;
				}
			}
		}
		if (! found)
			genesNotFound.push(gene);
		});
	}
	*/


	const config = {
		editable: true, // allows user to edit plot title, axis labels, and legend
		showLink: false
	}

	const plotContainer = document.getElementById("plot_container");
	plotContainer.replaceChildren();    // erase plot

	// NOTE: Plot initially is created to a default width but is responsive.
	// Noticed container within our "column" will make full-width go beyond the screen
	const divElt = generateElements('<div class="container is-max-desktop" id="plotly_preview"></div>');
	plotContainer.append(divElt);
	Plotly.purge("plotly_preview"); // clear old Plotly plots

	Plotly.newPlot("plotly_preview", plotData, layout, config);


	// If searched-for genes were not found, display under plot
	/*
	if (genesNotFound.length) {
		const genesNotFoundString = genesNotFound.join(", ");
		$("#genes_not_found").text(`Searched genes not found: ${genesNotFoundString}`);
		$("#genes_not_found").show();
	} else {
		$("#genes_not_found").hide();
	}
	*/

	const plotlyPreview = document.getElementById("plotly_preview");

	// If plot data is selected, create the right-column table and do other misc things
	plotlyPreview.on("plotly_selected", (eventData) => {
		selectedData = eventData;
		const selectedGeneData = [];

		eventData.points.forEach((pt) => {
			// Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
			// Each trace has its own "pointNumber" ids so gene symbols and pvalues needed to be passed in for each plotdata trace
			selectedGeneData.push({
				gene_symbol: pt.data.id[pt.pointNumber],
				pvals: statisticalTest ? pt.data.pvals[pt.pointNumber].toExponential(2) : "NA",
				foldchange: pt.data.foldchange[pt.pointNumber].toFixed(1),
			});
		});

		// Sort by adjusted p-value in descending order either by fold change or p-values
		selectedGeneData.sort((a, b) => b.foldchange - a.foldchange);
		if (statisticalTest)
			selectedGeneData.sort((a, b) => a.pvals - b.pvals);

		/*
		const template = $.templates("#selected_genes_tmpl");
		const htmlOutput = template.render(selectedGeneData);
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
		*/

	});

	/*window.onresize = () => {
		Plotly.Plots.resize(graphDiv);
	};*/

	const plotlyNote = generateElements(`
		<div id="editable_tip" class="message is-info is-light">
			<div class="message-body">
	  			<p><strong>Tip:</strong> You can click the plot title or axis labels to edit them. Hit Enter to apply edit.</p>
			</div>
  		</div>`);
	plotlyPreview.append(plotlyNote);
}

const populatePostCompareBox = (scope, series, group) => {
	// Find box
	const boxElt = document.querySelector(`#${scope}_post_c .notification`);
	boxElt.replaceChildren();
	// Add series as mini-subtitle and group as tag
	const seriesElt = document.createElement("div");
	seriesElt.classList.add("has-text-weight-semibold", "mb-2");
	seriesElt.textContent = series;
	boxElt.append(seriesElt);
	const groupElt = document.createElement("span");
	groupElt.classList.add("tag", "is-light", "is-primary");
	groupElt.textContent = group;
	boxElt.append(groupElt);

}

const sanitizeCondition = (condition) => {
	const sanitized_condition = {}
	for (const property in condition) {
		// If no groups for an observation are selected, delete filter
		if (condition[property].length) {
		sanitized_condition[property] = condition[property];
		}
	}
	return sanitized_condition;
}

// For plotting options, populate select menus with category groups
const updateSeriesOptions = (classSelector, seriesArray) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add categories
        for (const group of seriesArray.sort()) {

            const option = document.createElement("option");
            option.textContent = group;
            option.value = group;
            elt.append(option);
        }
    }
}

// For a given categorical series (e.g. "celltype"), update the "group" options
const updateGroupOptions = (classSelector, groupsArray) => {

    for (const elt of document.getElementsByClassName(classSelector)) {
        elt.replaceChildren();

        // Append empty placeholder element
        const firstOption = document.createElement("option");
        elt.append(firstOption);

        // Add categories
        for (const group of groupsArray.sort()) {
            const option = document.createElement("option");
            option.textContent = group;
            option.value = group;
            elt.append(option);
        }
    }
}

const validateCompareGroups = (event) => {
	const compareGroups = [...document.getElementsByClassName("js-compare-groups")].map((elt) => elt.value);
	// Filter out empty values and duplicates
	const uniqueGroups = [...new Set(compareGroups)].filter(x => x);
	// Get all unselected groups
	const series = document.querySelector(".js-compare").value;
	const seriesItems = getSeriesItems(series);
	const seriesNames = getSeriesNames(seriesItems);
	const unselectedGroups = seriesNames.filter((group) => !uniqueGroups.includes(group));

	for (const innerClassElt of document.getElementsByClassName("js-compare-groups")) {
		// enable all unselected groups
		for (const group of unselectedGroups) {
			const opt = innerClassElt.querySelector(`option[value="${group}"]`);
			opt.removeAttribute("disabled");
		}
		// disable unique groups in other compare groups
		for (const group of uniqueGroups) {
			if (innerClassElt.id !== event.target.id) {
				const opt = innerClassElt.querySelector(`option[value="${group}"]`);
				opt.setAttribute("disabled", "disabled");
			}
		}
	}
}

const validatePlotRequirements = (event) => {
    const elt = event.target;
    // Reset "status" classes
    elt.classList.remove("is-success", "is-danger");
    if (elt.value) {
        elt.parentElement.classList.remove("is-danger");
        elt.parentElement.classList.add("is-success");

        const validationElts = document.getElementsByClassName("js-plot-req");

        // If every validation param has been filled out, it's OK to plot
        // NOTE: need to ensure pre- and post- param elements are filled before this function is called
        if ([...validationElts].every(element => element.value)) {
            for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
                plotBtn.disabled = false;
            }
            document.getElementById("condition_compare_s_success").classList.remove("is-hidden");
            document.getElementById("condition_compare_s_failed").classList.add("is-hidden");
        }
        return;
    }

    // Required paramater has no value. Indicate it and disable plot buttons
    elt.parentElement.classList.add("is-danger");
    elt.parentElement.classList.remove("is-success");

    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
        plotBtn.disabled = true;
    }
    document.getElementById("condition_compare_s_success").classList.add("is-hidden");
    document.getElementById("condition_compare_s_failed").classList.remove("is-hidden");
}

/* --- Event listeners --- */

document.getElementById("statistical_test").addEventListener("change", (event) => {
	const pvalCutoff = document.getElementById("pval_cutoff");
	const cutoffFilterAction = document.getElementById("cutoff_filter_action");
	pvalCutoff.disabled = event.target.value ? false : true;
	cutoffFilterAction.disabled = event.target.value ? false : true;
});

const validationElts = document.getElementsByClassName("js-plot-req");
for (const classElt of validationElts ) {
	classElt.addEventListener("change", validatePlotRequirements);
}

// When compare series changes, update the compare groups
for (const classElt of document.getElementsByClassName("js-compare")) {
	classElt.addEventListener("change", async (event) => {
		const compareSeries = event.target.value;
		const seriesItems = getSeriesItems(compareSeries);
		const seriesNames = getSeriesNames(seriesItems);

		updateGroupOptions("js-compare-x", seriesNames);
		updateGroupOptions("js-compare-y", seriesNames);
	})
}

// When compare groups change, prevent the same group from being selected in the other compare groups
for (const classElt of document.getElementsByClassName("js-compare-groups")) {
	classElt.addEventListener("change", validateCompareGroups);
}

for (const classElt of document.getElementsByClassName("js-compare-x")) {
	classElt.addEventListener("change", (event) => {
		const compareX = event.target.value;
		const compareSeries = document.getElementById("compare_series").value;
		populatePostCompareBox("compare_x", compareSeries, compareX);
	})
}

for (const classElt of document.getElementsByClassName("js-compare-y")) {
	classElt.addEventListener("change", (event) => {
		const compareX = event.target.value;
		const compareSeries = document.getElementById("compare_series").value;
		populatePostCompareBox("compare_y", compareSeries, compareX);
	})
}

for (const classElt of document.getElementsByClassName("js-plot-btn")) {
	classElt.addEventListener("click", getComparisons);
}

document.getElementById("edit_params").addEventListener("click", (event) => {
    event.target.classList.add("is-loading");
    // Hide this view
    document.getElementById("content_c").classList.remove("is-hidden");
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").classList.add("is-hidden");

    event.target.classList.remove("is-loading");
})

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// Update with current page info
	document.querySelector("#header_bar .navbar-item").textContent = "Comparison Tool";
	for (const elt of document.querySelectorAll("#primary_nav > aside > ul > li > span > span > a")) {
		elt.classList.remove("is-active");
	}

	document.querySelector("a[tool='compare'").classList.add("is-active");

    sessionId = CURRENT_USER.session_id;
    colorblindMode = CURRENT_USER.colorblind_mode || false;
    Cookies.set('gear_session_id', sessionId, { expires: 7 });

	try {
		await loadDatasetTree()
        // If brought here by the "gene search results" page, curate on the dataset ID that referred us
        const urlParams = new URLSearchParams(window.location.search);
        if (urlParams.has("dataset_id")) {
            const linkedDatasetId = urlParams.get("dataset_id");
            try {
                // find DatasetTree node and trigger "activate"
                const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
                foundNode.setActive(true);
                datasetId = linkedDatasetId;
            } catch (error) {
                createToast(`Dataset id ${linkedDatasetId} was not found as a public/private/shared dataset`);
                throw new Error(error);
            }
        }
	} catch (error) {
		logErrorInConsole(error);
	}

};