'use strict';

import { apiCallsMixin, createToast, disableAndHideElement, getCurrentUser, logErrorInConsole, registerPageSpecificLoginUIUpdates } from './common.v2.js';
import { Gene, WeightedGene } from "./classes/gene.js";
import { GeneCart, WeightedGeneCart } from "./classes/genecart.v2.js";
import { DatasetTree } from "./classes/tree.js";
import { fetchGeneCartData, geneCollectionState } from '../include/gene-collection-selector/gene-collection-selector.js';

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
let facetWidget;	// stores aggregation data
let datasetId;
let organismId;	// Used for saving as gene cart
let compareData;;
let selectedGeneData;
let manuallyEnteredGenes = new Set();

// imported from gene-collection-selector.js
// let geneCollectionState.selectedGenes = new Set();

// Storing user's plot text edits, so they can be restored if user replots
let titleText = null;
let xaxisText = null;
let yaxisText = null;

const datasetTree = new DatasetTree({
    element: document.getElementById("dataset-tree")
    , searchElement: document.getElementById("dataset-query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current-dataset-c").classList.remove("is-hidden");
        document.getElementById("current-dataset").textContent = e.node.title;
        document.getElementById("current-dataset-post").textContent = e.node.title;

        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        datasetId = newDatasetId;

        // Hide tree by clicking the "toggle" button
        document.getElementById("btn-toggle-dataset-tree").click();

        // Click to get to next step
        document.getElementById("condition-compare-s").click();

        // Clear "success/failure" icons
        for (const elt of document.getElementsByClassName("js-step-success")) {
            elt.classList.add("is-hidden");
        }
		for (const elt of document.getElementsByClassName("js-step-failure")) {
			elt.classList.add("is-hidden");
		}

        // collapse tree
        e.node.tree.expandAll(false);

        document.getElementById("dataset-s-success").classList.remove("is-hidden");
        document.getElementById("dataset-s-failed").classList.add("is-hidden");

		const compareSeriesElt = document.getElementById("compare-series");
		compareSeriesElt.parentElement.classList.add("is-loading");

		// Clear selected gene tags
		document.getElementById("gene-tags").replaceChildren();

		// Clear compare groups
		for (const classElt of document.getElementsByClassName("js-compare-groups")) {
			classElt.replaceChildren();
		}

		// Update compare series options
		try {
			facetWidget = await createFacetWidget(datasetId, null, {}); // Initial fetching of categorical columns
			const catColumns = facetWidget.aggregations.map((agg) => agg.name);
			updateSeriesOptions("js-compare", catColumns);
		} catch (error) {
			return
		} finally {
			compareSeriesElt.parentElement.classList.remove("is-loading");
		}

    })
});

/**
 * Activates a dataset in the dataset tree based on a URL parameter.
 *
 * This function checks if the specified URL parameter exists, optionally fetches additional
 * dataset information using a provided function, and then activates the corresponding dataset
 * node in the dataset tree. If the dataset cannot be found or accessed, a toast notification
 * is displayed and an error is thrown.
 *
 * @async
 * @param {string} paramName - The name of the URL parameter to look for.
 * @param {function} [fetchInfoFn] - Optional async function to fetch dataset info using the parameter value.
 *        Should return a Promise that resolves to an array of objects containing a `dataset_id` property.
 * @throws {Error} If the dataset cannot be accessed or found in the dataset tree.
 */
const activateDatasetFromParam = async (paramName, fetchInfoFn) => {
    if (!urlParams.has(paramName)) {
        return;
    }
    const paramValue = urlParams.get(paramName);
    let linkedDatasetId;
    try {
        if (fetchInfoFn) {
            const data = await fetchInfoFn(paramValue);
            linkedDatasetId = data.datasets[0].id;
            if (!linkedDatasetId) {
                throw new Error(`Accessible dataset for ${paramName} ${paramValue} was not found`);
            }
        } else {
            linkedDatasetId = paramValue;
        }
    } catch (error) {
        createToast(error.message);
        throw new Error(error);
    }

    try {
        // find DatasetTree node and trigger "activate"
        const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
        foundNode.setActive(true, {focusTree:true});
        datasetTree.tree.setActiveNode(foundNode);
        datasetTree.selectCallback({node: foundNode});  // manually trigger the "activate" event.
        datasetId = linkedDatasetId;
    } catch (error) {
        createToast(`Dataset id ${linkedDatasetId} was not found as a public/private/shared dataset`);
        throw new Error(error);
    }
}

const adjustGeneTableLabels = () => {
    const geneFoldchanges = document.getElementById("tbl-gene-foldchanges");
	geneFoldchanges.replaceChildren();
	const logBase = document.getElementById("log-base").value;

	const spanIcon = document.createElement("span");
	spanIcon.classList.add("icon");
	const i = document.createElement("i");
	i.classList.add("mdi", "mdi-sort-numeric-ascending");
	i.setAttribute("aria-hidden", "true");
	spanIcon.appendChild(i);
	geneFoldchanges.appendChild(spanIcon);

	if (logBase === "raw") {
		geneFoldchanges.prepend("Fold Change ");
		return;
	}
	geneFoldchanges.prepend(`Log${logBase} Fold Change `);
}

const appendGeneTagButton = (geneTagElt) => {
    // Add delete button
    const deleteBtnElt = document.createElement("button");
    deleteBtnElt.classList.add("delete", "is-small");
    geneTagElt.appendChild(deleteBtnElt);
    deleteBtnElt.addEventListener("click", (event) => {
        // Remove gene from geneCollectionState.selectedGenes
        const gene = event.target.parentNode.textContent;
		geneCollectionState.selectedGenes.delete(gene);
		event.target.parentNode.remove();

		// Remove gene from manually entered genes textbox
		manuallyEnteredGenes.delete(gene);
		document.getElementById("genes-manually-entered").value = Array.from(manuallyEnteredGenes).join(" ");

		// Update graph
		updatePlotAnnotations(Array.from(geneCollectionState.selectedGenes).sort());

        // Remove checkmark from gene lists dropdown
        const geneListLabel = document.querySelector(`#dropdown-content-genes .gene-item-label[text="${gene}"]`);
        if (!geneListLabel) {
            return;
        }
        const geneListElt = geneListLabel.parentElement;
        const geneListI = geneListElt.querySelector("i.toggler")
        if (!geneListI.classList.contains("mdi-check")) {
            return;
        }
        geneListI.classList.replace("mdi-check", "mdi-plus");
        geneListI.classList.replace("gene-list-item-remove", "gene-list-item-add");
        geneListElt.classList.remove("is-selected");


    });
}

const chooseGenes = (event) => {
    // Triggered when a gene is selected

    // Delete existing tags
    const geneTagsElt = document.getElementById("gene-tags");
    geneTagsElt.replaceChildren();

	if (geneCollectionState.selectedGenes.size == 0) return;  // Do not trigger after initial population

    // Update list of gene tags
	const sortedGenes = Array.from(geneCollectionState.selectedGenes).sort();
    for (const opt in sortedGenes) {
        const geneTagElt = document.createElement("span");
        geneTagElt.classList.add("tag", "is-primary", "mx-1");
        geneTagElt.textContent = sortedGenes[opt];
        appendGeneTagButton(geneTagElt);
        geneTagsElt.appendChild(geneTagElt);
    }

    document.getElementById("gene-tags-c").classList.remove("is-hidden");

    // If more than 10 tags, hide the rest and add a "show more" button
    if (geneCollectionState.selectedGenes.size > 10) {
        const geneTags = geneTagsElt.querySelectorAll("span.tag");
        for (let i = 10; i < geneTags.length; i++) {
            geneTags[i].classList.add("is-hidden");
        }
        // Add show more button
        const showMoreBtnElt = document.createElement("button");
        showMoreBtnElt.classList.add("tag", "button", "is-small", "is-primary", "is-light");
        const numToDisplay = geneCollectionState.selectedGenes.size - 10;
        showMoreBtnElt.textContent = `+${numToDisplay} more`;
        showMoreBtnElt.addEventListener("click", (event) => {
            const geneTags = geneTagsElt.querySelectorAll("span.tag");
            for (let i = 10; i < geneTags.length; i++) {
                geneTags[i].classList.remove("is-hidden");
            }
            event.target.remove();
        });
        geneTagsElt.appendChild(showMoreBtnElt);
    }

	updatePlotAnnotations(sortedGenes);

}

const clearGenes = (event) => {
    document.getElementById("clear-genes-btn").classList.add("is-loading");
	document.getElementById("gene-tags").replaceChildren();
	geneCollectionState.selectedGenes.clear();
	document.getElementById("dropdown-gene-list-cancel").click();	// clear the dropdown
	// Remove gene from manually entered genes textbox
	manuallyEnteredGenes.clear();
	document.getElementById("genes-manually-entered").value = "";

	updatePlotAnnotations([]);
    document.getElementById("clear-genes-btn").classList.remove("is-loading");
}

const createFacetWidget = async (datasetId, analysisId, filters) => {
    document.getElementById("selected-facets-loader").classList.remove("is-hidden")
	document.getElementById("facet-content").classList.add("is-hidden");
	document.getElementById("selected-facets").classList.add("is-hidden");

	try {
    	const {aggregations, total_count:totalCount} = await fetchAggregations(datasetId, analysisId, filters);
    	document.getElementById("num-selected").textContent = totalCount;
	} catch (error) {
		document.getElementById("num-selected").textContent = "0";
		throw error;
	}

    const facetWidget = new FacetWidget({
        aggregations,
        filters,
        onFilterChange: async (filters) => {
            if (filters) {
                try {
                    const {aggregations, total_count:totalCount} = await fetchAggregations(datasetId, analysisId, filters);
                    facetWidget.updateAggregations(aggregations);
                    document.getElementById("num-selected").textContent = totalCount;
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
    document.getElementById("selected-facets-loader").classList.add("is-hidden")
	document.getElementById("facet-content").classList.remove("is-hidden");
	document.getElementById("selected-facets").classList.remove("is-hidden");
    return facetWidget;
}

/**
 * Updates the UI to display tags for the selected filters.
 * Converts selected filters into a "category:group" format and creates
 * corresponding tags to display in the "selected-filters-tags" element.
 * If no filters are selected, the "selected-filters-all" element is shown instead.
 *
 * @param {Object} filters - An object representing the selected filters.
 *                            Each key is a filter category, and its value is an array of selected groups.
 */
const createFilterTags = (filters) => {
	const selectedFiltersTags = document.getElementById('selected-filters-tags');
	const allSelectedTag = document.getElementById("selected-filters-all");

	// For each selected filter, convert into a "category:group" format
	const filterGroups = [];
	for (const filter in filters) {
		for (const group of filters[filter]) {
			filterGroups.push(`${filter}:${group}`);
		}
	}
	selectedFiltersTags.replaceChildren();
	allSelectedTag.classList.remove("is-hidden");

	if (filterGroups.length) {
		allSelectedTag.classList.add("is-hidden");
		// Create a tag for each selected filter
		for (const filterTag of filterGroups) {
			const tag = document.createElement("span");
			tag.classList.add("tag", "is-dark", "is-rounded", "is-small", "mx-1");
			tag.textContent = filterTag;
			selectedFiltersTags.appendChild(tag);
		}
	}
}

const downloadSelectedGenes = (event) => {
    event.preventDefault();

	// Builds a file in memory for the user to download.  Completely client-side.
	// plot_data contains three keys: x, y and symbols
	// build the file string from this

    // Adjust headers to the plot type
	const xLabel = JSON.stringify([...document.querySelectorAll("#compare-x input:checked")].map((elt) => elt.value));
	const yLabel = JSON.stringify([...document.querySelectorAll("#compare-y input:checked")].map((elt) => elt.value));

	const logBase = document.getElementById("log-base").value;

	let fileContents =
		logBase === "raw"
		? "gene_symbol\tp-value\traw fold change\t"
		+ xLabel + "\t"
		+ yLabel + "\n"
		: "gene_symbol\tp-value\traw fold change\t"
		+ xLabel + " (log" + logBase +")\t"
		+ yLabel + " (log" + logBase +")\n";


	selectedGeneData.forEach((gene) => {
		// Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
		fileContents +=
			`${gene.gene_symbol}\t`
			+ `${gene.pval}\t`
			+ `${gene.foldchange}\t`
			+ `${gene.x}\t`
			+ `${gene.y}\n`;
	});

	const element = document.createElement("a");
	element.setAttribute(
		"href",
		`data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(fileContents)}`
	);
	element.setAttribute("download", "geneCollectionState.selectedGenes.tsv");
	element.style.display = "none";
	document.body.appendChild(element);
	element.click();
	document.body.removeChild(element);
}

const fetchAggregations = async (datasetId, analysisId, filters) => {
    try {
        const data = await apiCallsMixin.fetchAggregations(datasetId, analysisId, filters)
        if (data.hasOwnProperty("success") && data.success < 1) {
            throw new Error(data?.message || "Could not fetch number of observations for this dataset. Please contact the gEAR team.");
        }
        const {aggregations, total_count} = data;
        return {aggregations, total_count};
    } catch (error) {
        logErrorInConsole(error);
		createToast(error.message);
		throw error
    }
}

const fetchDatasetComparison = async (datasetId, filters, compareKey, conditionX, conditionY, foldChangeCutoff, stDevNumCutoff, logBase, statisticalTestAction) => {
	try {
		return await apiCallsMixin.fetchDatasetComparison(datasetId, filters, compareKey, conditionX, conditionY, foldChangeCutoff, stDevNumCutoff, logBase, statisticalTestAction);
	} catch (error) {
		const msg = "Could not fetch dataset comparison. Please contact the gEAR team."
		throw new Error(msg);
	}
}

const getComparisons = async (event) => {

	// set loading icon
	event.target.classList.add("is-loading");

	const filters = JSON.stringify(facetWidget.filters);
	createFilterTags(facetWidget.filters);

	const compareSeries = document.getElementById("compare-series").value

	// Get all checked x and y series
	const checkedX = JSON.stringify([...document.querySelectorAll("#compare-x input:checked")].map((elt) => elt.value));
	const checkedY = JSON.stringify([...document.querySelectorAll("#compare-y input:checked")].map((elt) => elt.value));

	const foldChangeCutoff = document.getElementById("fc-cutoff").value;
	const stdDevNumCutoff = document.getElementById("standard-deviation").value;
	const logTransformation = document.getElementById("log-base").value;
	const statisticalTest = document.getElementById("statistical-test").value;

	try {
		const data = await fetchDatasetComparison(datasetId, filters, compareSeries, checkedX, checkedY, foldChangeCutoff, stdDevNumCutoff, logTransformation, statisticalTest);
		if (data?.success < 1) {
			throw new Error(data?.message || "Could not fetch dataset comparison. Please contact the gEAR team.");
		}
		compareData = data;
		plotDataToGraph(compareData);

		// If any genes selected, update plot annotations (since plot was previously purged)
		const sortedGenes = Array.from(geneCollectionState.selectedGenes).sort();
		updatePlotAnnotations(sortedGenes);

		// Hide this view
		document.getElementById("content-c").classList.add("is-hidden");
		// Generate and display "post-plotting" view/container
		document.getElementById("post-plot-content-c").classList.remove("is-hidden");

	} catch (error) {
		console.error(error);
		handleGetComparisonError(datasetId, checkedX, checkedY);
	} finally {
		event.target.classList.remove("is-loading");
	}

	// When a plot configuration ID is selected, populate the plot configuration post textbox
	const plotConfigElts = ["statistical-test", "pval-cutoff", "cutoff-filter-action", "log-base", "fc-cutoff", "standard-deviation"];
	for (const elt of plotConfigElts) {
		// if value is empty, set to "None", or if disabled, set to "N/A"
		let value = document.getElementById(elt).disabled ? "N/A" : document.getElementById(elt).value || "None"

		// Append extra flavor text
		if (elt == "log-base" && !(value === "raw" )) {
			value = `log${value}`
		}

		if (elt == "standard-deviation" && !(value === "0" )) {
			value = `±${value}`
		}

		document.getElementById(`${elt}-post`).textContent = value;
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

/**
 * Highlights rows in a gene table based on the searched genes.
 *
 * @param {string[]} searchedGenes - An array of genes to search for.
 */
const highlightTableGenes = (searchedGenes=[]) => {
	const geneTableBody = document.getElementById("gene-table-body");
	// Select the first column (gene_symbols) in each row
	for (const row of geneTableBody.children) {
		// clear any previous highlighting
		row.classList.remove("has-background-success");
		const tableGene = row.children[0].textContent;
		for (const gene of searchedGenes) {
			if (gene.toLowerCase() === tableGene.toLowerCase() ) {
				row.classList.add("has-background-success");
			}
		};
	}
}

/* Transform and load dataset data into a "tree" format */
const loadDatasetTree = async (shareId) => {
    const userDatasets = [];
    const sharedDatasets = [];
    const domainDatasets = [];
    try {
        const datasetData = await apiCallsMixin.fetchAllDatasets(shareId);

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
		createToast("Could not fetch datasets. Please contact the gEAR team.");
        document.getElementById("dataset-s-success").classList.add("is-hidden");
        document.getElementById("dataset-s-failed").classList.remove("is-hidden");
    }
}

const plotDataToGraph = (data) => {

	const statisticalTest = document.getElementById("statistical-test").value;

	const pointLabels = [];
	const performRanking = statisticalTest ? true : false;

	const plotData = [];

	if (performRanking) {
		const pValCutoff = document.getElementById("pval-cutoff").value;
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

		const passColor = getCurrentUser().colorblind_mode ? 'rgb(0, 34, 78)' : "#FF0000";
		const failColor = getCurrentUser().colorblind_mode ? 'rgb(254, 232, 56)' : "#A1A1A1";

		const statAction = document.getElementById("cutoff-filter-action").value;
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
						color: new Array(passing.x.length).fill(passColor, 0, passing.x.length),
						size: 4,
					},
				}
			// store original marker color as a deep copy
			passingObj.marker.origColor = JSON.parse(JSON.stringify(passingObj.marker.color));

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
						color: new Array(failing.x.length).fill(failColor, 0, failing.x.length),
						size: 4,
					},
				}
			// store original marker color as a deep copy
			failingObj.marker.origColor = JSON.parse(JSON.stringify(failingObj.marker.color));

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
					color: new Array(passing.x.length).fill("#000000" , 0, passing.x.length),
					size: 4,
				},
			}
			// store original marker color as a deep copy
			passingObj.marker.origColor = JSON.parse(JSON.stringify(passingObj.marker.color));

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
				color: new Array(data.x.length).fill("#000000", 0, data.x.length),
				size: 4,
			},
		}
		// store original marker color as a deep copy
		dataObj.marker.origColor = JSON.parse(JSON.stringify(dataObj.marker.color));

		plotData.push(dataObj);
	}

	const layout = {
		title: titleText || "Dataset Comparison",
		xaxis: {
			title: xaxisText || data.condition_x.join(", "),
		},
		yaxis: {
			title: yaxisText || data.condition_y.join(", "),
		},
		annotations: [],
		hovermode: "closest",
		dragmode: "select",
	};


	const config = {
		editable: true, // allows user to edit plot title, axis labels, and legend
		showLink: false
	}

	const plotContainer = document.getElementById("plot-container");
	plotContainer.replaceChildren();    // erase plot

	// NOTE: Plot initially is created to a default width but is responsive.
	// Noticed container within our "column" will make full-width go beyond the screen

	const plotlyPreview = document.createElement("div");
	plotlyPreview.id = "plotly-preview";
	plotlyPreview.classList.add("container", "is-max-desktop");
	plotContainer.append(plotlyPreview);

	Plotly.purge("plotly-preview"); // clear old Plotly plots

	Plotly.newPlot("plotly-preview", plotData, layout, config);

	// Hide table when plot is first loaded
	document.getElementById("tbl-selected-genes").classList.add("is-hidden");

	// If plot data is selected, create the right-column table and do other misc things
	plotlyPreview.on("plotly_selected", (eventData) => {

		// Hide selected genes table and disable unweighted radio button if no genes are selected
		document.getElementById("tbl-selected-genes").classList.add("is-hidden");
		document.getElementById("download-selected-genes-btn").classList.add("is-hidden");
		document.querySelector("input[name='genecart_type'][value='unweighted']").disabled = true;
		document.querySelector("input[name='genecart_type'][value='unweighted']").parentElement.setAttribute("disabled", "disabled");
		// click the weighted radio button
		document.querySelector("input[name='genecart_type'][value='weighted']").click();

		if (eventData?.points.length) {
			document.getElementById("tbl-selected-genes").classList.remove("is-hidden");
			document.getElementById("download-selected-genes-btn").classList.remove("is-hidden");
			document.querySelector("input[name='genecart_type'][value='unweighted']").disabled = false;
			document.querySelector("input[name='genecart_type'][value='unweighted']").parentElement.removeAttribute("disabled");
			// click the unweighted radio button
			document.querySelector("input[name='genecart_type'][value='unweighted']").click();

			adjustGeneTableLabels();
			populateGeneTable(eventData);
		}

		// Get genes from gene tags
		const geneTags = document.querySelectorAll("#gene-tags span.tag");
		const searchedGenes = [];
		for (const tag of geneTags) {
			searchedGenes.push(tag.textContent);
		}

		// Highlight table rows that match searched genes
		if (searchedGenes) {
			highlightTableGenes(searchedGenes);
		}
	});

	// Handler for when plot text is edited
	plotlyPreview.on("plotly_relayout", (eventData) => {
		// If plot title, x-axis or y-axis is edited, save the text
		if (eventData?.["title.text"]) {
			titleText = eventData["title.text"];
		}
		if (eventData?.["xaxis.title.text"]) {
			xaxisText = eventData["xaxis.title.text"];
		}
		if (eventData?.["yaxis.title.text"]) {
			yaxisText = eventData["yaxis.title.text"];
		}
	});

	const plotlyNote = document.createElement("div");
	plotlyNote.id = "tip-on-editing";
	plotlyNote.classList.add("notification", "content", "is-info", "is-light");
	plotlyNote.innerHTML = `<p><strong>Tip:</strong> Use the Plotly box and lasso select tools (upper-right) to select genes to view as a table. Double-clicking clears any selections.</p>

	<p>You can also click the plot title or axis labels to edit them. Hit Enter to apply edit.</p>`;
	plotlyPreview.append(plotlyNote);
}


const populateGeneTable = (data) => {
	const statisticalTest = document.getElementById("statistical-test").value;

	selectedGeneData = [];

	data.points.forEach((pt) => {
		// Some warnings on using toFixed() here: https://stackoverflow.com/a/12698296/1368079
		// Each trace has its own "pointNumber" ids so gene symbols and pvalues needed to be passed in for each plotdata trace
		selectedGeneData.push({
			gene_symbol: pt.data.id[pt.pointNumber],
			pval: statisticalTest ? pt.data.pvals[pt.pointNumber].toExponential(2) : "NA",
			foldchange: pt.data.foldchange[pt.pointNumber].toFixed(1),
			x: pt.data.x[pt.pointNumber].toFixed(1),
			y: pt.data.y[pt.pointNumber].toFixed(1)
		});
	});

	// Sort by adjusted p-value in descending order either by fold change or p-values
	selectedGeneData.sort((a, b) => b.foldchange - a.foldchange);
	if (statisticalTest)
		selectedGeneData.sort((a, b) => a.pval - b.pval);


    const geneTableBody = document.getElementById("gene-table-body");
    geneTableBody.replaceChildren();

    for (const gene of selectedGeneData) {
        const row = document.createElement("tr");
        row.innerHTML = `<td>${gene.gene_symbol}</td><td>${gene.pval}</td><td>${gene.foldchange}</td>`;
        geneTableBody.appendChild(row);
    }

	// If not statistical test, hide p-value column (deleting can cause issues with subsequent calls to this function)
	const pvalColumn = document.getElementById("tbl-gene-pvalues");
	pvalColumn.classList.remove("is-hidden");
	for (const pvalCell of document.querySelectorAll("#tbl-selected-genes tbody tr td:nth-child(2)")) {
		pvalCell.classList.remove("is-hidden");
	}
	if (!statisticalTest) {
		pvalColumn.classList.add("is-hidden");
		for (const pvalCell of document.querySelectorAll("#tbl-selected-genes tbody tr td:nth-child(2)")) {
			pvalCell.classList.add("is-hidden");
		}
	}
	// Should be sorted by logFC now

}

const populatePostCompareBox = (scope, series, groups) => {
	// Find box
	const boxElt = document.querySelector(`#${scope}-post-c .notification`);
	boxElt.replaceChildren();

	// Add series as mini-subtitle and group as tag
	const seriesElt = document.createElement("div");
	seriesElt.classList.add("has-text-weight-semibold", "mb-2");
	seriesElt.textContent = series;
	boxElt.append(seriesElt);

	const tagsElt = document.createElement("div");
	tagsElt.classList.add("tags");
	boxElt.append(tagsElt);

	for (const group of groups) {

		const groupElt = document.createElement("span");
		groupElt.classList.add("tag", "is-dark", "is-rounded");
		groupElt.textContent = group;
		tagsElt.append(groupElt);
	}
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

const saveGeneCart = () => {
    // must have access to USER_SESSION_ID
    const gc = new GeneCart({
        session_id: sessionId
        , label: document.getElementById("new-genecart-label").value
        , gctype: "unweighted-list"
        , organism_id:  organismId
        , is_public: 0
    });

    for (const sg of selectedGeneData) {
        const gene = new Gene({
            id: sg.ensembl_id, // Ensembl ID stored in "customdata" property
            gene_symbol: sg.gene_symbol,
        });
        gc.addGene(gene);
    }

    gc.save(updateUIAfterGeneCartSaveSuccess, updateUIAfterGeneCartSaveFailure);
}

const saveWeightedGeneCart = () => {

	// must have access to USER_SESSION_ID


	// Saving raw FC by default so it is easy to transform weight as needed
	const weightLabels = ["FC"];

	const gc = new WeightedGeneCart({
		session_id: sessionId
		, label:  document.getElementById("new-genecart-label").value
		, gctype: 'weighted-list'
		, organism_id: organismId
		, is_public: 0
	}, weightLabels);

	compareData.gene_ids.forEach((gene_id, i) => {
		const weights = [compareData.fold_changes[i]];

		const gene = new WeightedGene({
			id: gene_id,
			gene_symbol: compareData.symbols[i]
		}, weights);
		gc.addGene(gene);
	});

	gc.save(updateUIAfterGeneCartSaveSuccess, updateUIAfterGeneCartSaveFailure);
}

// Taken from https://www.w3schools.com/howto/howto_js_sort_table.asp
const sortGeneTable = (mode) => {
	let table;
	let rows;
	let switching;
	let i;
	let x;
	let y;
	let shouldSwitch;
	let dir;
	let switchcount = 0;
	table = document.getElementById("tbl-selected-genes");

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
            x = rows[i].getElementsByTagName("td")[mode];
            y = rows[i + 1].getElementsByTagName("td")[mode];
            /* Check if the two rows should switch place,
                based on the direction, asc or desc: */
            if (dir == "asc") {
                // First column is gene_symbol... rest are numbers
                if (mode === 0 && x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                    // If so, mark as a switch and break the loop:
                    shouldSwitch = true;
                    break;
                }
                if (Number(x.innerHTML) > Number(y.innerHTML)) {
                    shouldSwitch = true;
                    break;
                }
            } else if (dir == "desc") {
                if (mode === 0 && x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
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

    // Reset other sort icons to "ascending" state, to show what direction they will sort when clicked
    const otherTblHeaders = document.querySelectorAll(`.js-tbl-gene-header:not(:nth-child(${mode + 1}))`);
    for (const tblHeader of otherTblHeaders) {
        const currIcon = tblHeader.querySelector("i");
        if (mode == 0) {
            currIcon.classList.remove("mdi-sort-alphabetical-descending");
            currIcon.classList.add("mdi-sort-alphabetical-ascending");
        } else {
            currIcon.classList.remove("mdi-sort-numeric-descending");
            currIcon.classList.add("mdi-sort-numeric-ascending");
        }
    }

    // toggle the mdi icons between ascending / descending
    // icon needs to reflect the current state of the sort
    const selectedTblHeader = document.querySelector(`.js-tbl-gene-header:nth-child(${mode + 1})`);
    const currIcon = selectedTblHeader.querySelector("i");
    if (dir == "asc") {
        if (mode == 0) {
            currIcon.classList.remove("mdi-sort-alphabetical-descending");
            currIcon.classList.add("mdi-sort-alphabetical-ascending");
        } else {
            currIcon.classList.remove("mdi-sort-numeric-descending");
            currIcon.classList.add("mdi-sort-numeric-ascending");
        }
    } else {
        if (mode == 0) {
            currIcon.classList.remove("mdi-sort-alphabetical-ascending");
            currIcon.classList.add("mdi-sort-alphabetical-descending");
        } else {
            currIcon.classList.remove("mdi-sort-numeric-ascending");
            currIcon.classList.add("mdi-sort-numeric-descending");
        }
    }
}

// For a given categorical series (e.g. "celltype"), add checkboxes for each category
const updateGroupOptions = (selectorId, groupsArray, series) => {

	const elt = document.getElementById(selectorId);
	elt.classList.remove("is-hidden");

	// Add categories
	for (const group of groupsArray.sort()) {

		const checkbox = document.createElement("input");
		checkbox.type = "checkbox";
		checkbox.id = `${selectorId}-${group}`;
		checkbox.name = group;
		checkbox.value = group;

		const label = document.createElement("label");
		label.classList.add("checkbox");
		label.htmlFor = `${selectorId}-${group}`;
		label.textContent = ` ${group}`;
		label.prepend(checkbox);

		// If group has aggregations count of 0 (no data after filtering), disable checkbox and label
		checkbox.disabled = false;
		label.removeAttribute("disabled");
		if (facetWidget.aggregations.find((agg) => agg.name === series).items.find((item) => item.name === group).count === 0) {
			checkbox.disabled = true;
			label.setAttribute("disabled", "disabled");
		}


		// create .control div to ensure checkboxs are vertically aligned
		const control = document.createElement("div");
		control.classList.add("control", "m-1");
		control.append(label);

		elt.append(control);
	}
}

// Update the plotly graph with the selected genes
const updatePlotAnnotations = (genes) => {
	// Take genes to search for and highlight their datapoint in the plot

	const plotlyPreview = document.getElementById("plotly-preview");
	const plotData = plotlyPreview.data;
	const layout = plotlyPreview.layout;

	const annotationColor = getCurrentUser().colorblind_mode ? "orange" : "cyan";

	layout.annotations = [];

	// Reset all trace colors
	for (const trace of plotData) {
		trace.id.forEach((element, i) => {
			trace.marker.color[i] = trace.marker.origColor[i];
		});
	}

	genes.forEach((gene) => {
		for (const trace of plotData) {
			trace.id.forEach((element, i) => {
				if (gene.toLowerCase() !== element.toLowerCase() ) {
					return;
				}

				// If gene is found add an annotation arrow
				layout.annotations.push({
					xref: "x",
					yref: "y",
					x: trace.x[i],
					y: trace.y[i],
					text:element,
					bgcolor: annotationColor,
					showarrow: true,
					arrowcolor: annotationColor,
					opacity: 0.8,

				});

				// change trace dot to match annotation
				trace.marker.color[i] = annotationColor;

			});
		}
	});

	// If no annotations, add warning text that all genes were filtered out
	if (!layout.annotations.length && genes.length) {
		layout.annotations.push({
			xref: "paper",
			yref: "paper",
			x: 0,
			y: 1,
			text: "No selected genes were found in this plot.",
			bgcolor: "lightyellow",
			showarrow: false,
			opacity: 0.8,
		});
	}

	// update the Plotly layout
	Plotly.relayout(plotlyPreview, layout);

	// Update table highlighting
	highlightTableGenes(genes);

}

// For plotting options, populate select menus with category groups
const updateSeriesOptions = (classSelector, seriesArray) => {

	if (!seriesArray.length) {
		createToast("No categorical data series found for this dataset. Please choose another dataset or contact the gEAR team.");
		return;
	}

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

const updateUIAfterGeneCartSaveSuccess = (gc) => {
	createToast("Gene cart saved successfully", "is-success");
}

const updateUIAfterGeneCartSaveFailure = (gc, message) => {
    createToast(message);
}

const validateCompareGroups = (event) => {

	const group = event.target.value;
	const checked = event.target.checked;

	for (const innerClassElt of document.getElementsByClassName("js-compare-groups")) {
		// BUG: Checking via label click enables X compare group
		if (checked) {
			// disable unique groups in other compare groups
			if (innerClassElt.closest(".js-compare-groups").id !== event.target.closest(".js-compare-groups").id) {
				const checkbox = innerClassElt.querySelector(`input[value="${group}"]`);
				checkbox.setAttribute("disabled", "disabled");
				checkbox.parentElement.setAttribute("disabled", "disabled");
			}
		} else {
			const checkbox = innerClassElt.querySelector(`input[value="${group}"]`);
			checkbox.removeAttribute("disabled");
			checkbox.parentElement.removeAttribute("disabled");
		}
	}
}

const validatePlotRequirements = (event) => {
    const elt = event.target;
    // Reset "status" classes
    elt.classList.remove("is-success", "is-danger");

    for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
        plotBtn.disabled = true;
    }

	// We need at least one compare-x and one compare-y checkbox checked
	const checkedX = [...document.querySelectorAll(".js-compare-x input:checked")].map((elt) => elt.value);
	const checkedY = [...document.querySelectorAll(".js-compare-y input:checked")].map((elt) => elt.value);

	if (checkedX.length && checkedY.length) {
		// Enable plot button
		for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
			plotBtn.disabled = false;

			document.getElementById("condition-compare-s-failed").classList.add("is-hidden");
		}
		return;
	}

    document.getElementById("condition-compare-s-success").classList.add("is-hidden");
}

/* --- Event listeners --- */

document.getElementById("statistical-test").addEventListener("change", (event) => {
	const pvalCutoff = document.getElementById("pval-cutoff");
	const cutoffFilterAction = document.getElementById("cutoff-filter-action");
	pvalCutoff.disabled = event.target.value ? false : true;
	cutoffFilterAction.disabled = event.target.value ? false : true;
});

// When compare series changes, update the compare groups
for (const classElt of document.getElementsByClassName("js-compare")) {
	const compareSeriesNotification = document.getElementById("select-compare-series-notification");
	classElt.addEventListener("change", async (event) => {

		// Disable plot button until conditions are met
		for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
			plotBtn.disabled = true;
		}

		// Hide and clear compare groups
		for (const classElt of document.getElementsByClassName("js-compare-groups")) {
			classElt.parentElement.classList.add("is-hidden");
			classElt.replaceChildren();
		}

		const compareSeries = event.target.value;
		compareSeriesNotification.classList.remove("is-hidden", "is-danger");
		compareSeriesNotification.classList.add("is-warning");
		compareSeriesNotification.textContent = "Please select a series to compare to choose X and Y conditions";

		if (!compareSeries) return;

		const seriesItems = getSeriesItems(compareSeries);
		const seriesNames = getSeriesNames(seriesItems);

		// at least 2 of the series items must have 1+ aggregation total, or else we can't compare
		const seriesAggCountsFiltered = seriesItems.filter((item) => item.count > 0);
		if (seriesAggCountsFiltered.length < 2) {
			compareSeriesNotification.classList.remove("is-warning");
			compareSeriesNotification.classList.add("is-danger");
			compareSeriesNotification.textContent = `At least 2 groups within ${compareSeries} must each have one or more observations (after filtering) to compare`;
			return;
		}

		compareSeriesNotification.classList.add("is-hidden");

		// Reset facet groups
		facetWidget = await createFacetWidget(datasetId, null, {}); // will also remove any existing facet widget

		updateGroupOptions("compare-x", seriesNames, compareSeries);
		updateGroupOptions("compare-y", seriesNames, compareSeries);

		// Show compare groups since things have validated.
		for (const classElt of document.getElementsByClassName("js-compare-groups")) {
			classElt.parentElement.classList.remove("is-hidden");
		}

		// Hide the chosen group's facet element
		const facetElt = document.getElementById(`filter-${compareSeries}`);
		facetElt.classList.add("is-hidden");

	})
}

// When compare groups change, prevent the same group from being selected in the other compare groups
for (const classElt of document.getElementsByClassName("js-compare-groups")) {
	classElt.addEventListener("change", validateCompareGroups);
	classElt.addEventListener("change", validatePlotRequirements);
}

for (const classElt of document.getElementsByClassName("js-compare-x")) {
	classElt.addEventListener("change", (event) => {
	// We need at least one compare-x and one compare-y checkbox checked
		const checkedX = [...document.querySelectorAll(".js-compare-x input:checked")].map((elt) => elt.value);
		const compareSeries = document.getElementById("compare-series").value;
		populatePostCompareBox("compare-x", compareSeries, checkedX);
	})
}

for (const classElt of document.getElementsByClassName("js-compare-y")) {
	classElt.addEventListener("change", (event) => {
		const checkedY = [...document.querySelectorAll(".js-compare-y input:checked")].map((elt) => elt.value);
		const compareSeries = document.getElementById("compare-series").value;
		populatePostCompareBox("compare-y", compareSeries, checkedY);
	})
}

for (const classElt of document.getElementsByClassName("js-plot-btn")) {
	classElt.addEventListener("click", getComparisons);
}

document.getElementById("edit-params").addEventListener("click", (event) => {
    event.target.classList.add("is-loading");
	// Clear existing plot annotations
	clearGenes();
    // Hide this view
    document.getElementById("content-c").classList.remove("is-hidden");
    // Generate and display "post-plotting" view/container
    document.getElementById("post-plot-content-c").classList.add("is-hidden");

    event.target.classList.remove("is-loading");
})

document.getElementById("clear-genes-btn").addEventListener("click", clearGenes);

document.getElementById("new-genecart-label").addEventListener("input", (event) => {
    const saveBtn = document.getElementById("save-genecart-btn");
    saveBtn.disabled = event.target.value ? false : true;
});

document.getElementById("save-genecart-btn").addEventListener("click", (event) => {
    event.preventDefault();
    event.target.classList.add("is-loading");
    // get value of genecart radio button group
    const geneCartName = document.querySelector("input[name='genecart_type']:checked").value;
    if (getCurrentUser()) {
        if (geneCartName === "unweighted") {
            saveGeneCart();
        } else {
            saveWeightedGeneCart();
        }
    }
    event.target.classList.remove("is-loading");
});

// handle when the dropdown-gene-list-search-input input box is changed
document.getElementById('genes-manually-entered').addEventListener('change', (event) => {
    const searchTermString = event.target.value;
    const newManuallyEnteredGenes = searchTermString.length > 0 ? new Set(searchTermString.split(/[ ,]+/)) : new Set();

    // Remove genes that have been deleted from the geneCollectionState.selectedGenes set
    for (const gene of manuallyEnteredGenes) {
        if (!newManuallyEnteredGenes.has(gene)) {
            geneCollectionState.selectedGenes.delete(gene);
        }
    }

    // Add new genes to the geneCollectionState.selectedGenes set
    for (const gene of newManuallyEnteredGenes) {
        geneCollectionState.selectedGenes.add(gene);
    }

    manuallyEnteredGenes = newManuallyEnteredGenes;
    chooseGenes(null);
});

// Use the dataset input selector button to toggle the dataset selection div
document.getElementById("btn-toggle-dataset-tree").addEventListener("click", (event) => {
    // Toggle the dataset selection div
    const selectionDiv = document.getElementById("dataset-selection-c");
    if (selectionDiv.classList.contains("is-hidden")) {
        selectionDiv.classList.remove("is-hidden");
        event.target.textContent = "Collapse dataset selection tool";
    } else {
        selectionDiv.classList.add("is-hidden");
        event.target.textContent = "Expand dataset selection tool";
    }
});

document.getElementById('dropdown-gene-list-proceed').addEventListener('click', chooseGenes);

document.getElementById("download-selected-genes-btn").addEventListener("click", downloadSelectedGenes);

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// Update with current page info
	document.getElementById("page-header-label").textContent = "Comparison Tool";
    sessionId = getCurrentUser().session_id;

	if (! sessionId ) {
		// TODO: Add master override to prevent other triggers from enabling saving
        createToast("Not logged in so saving gene carts is disabled.", "is-warning");
        disableAndHideElement(document.getElementById("save-genecart-btn"));
    }


	try {
        // If brought here by the "gene search results" page, curate on the dataset ID that referred us
        const urlParams = new URLSearchParams(window.location.search);
        const shareId = urlParams.get("share_id");

		await Promise.all([
			loadDatasetTree(shareId),
			fetchGeneCartData()
		]);
        // Usage inside handlePageSpecificLoginUIUpdates
        if (urlParams.has("share_id")) {
            return await activateDatasetFromParam("share_id", async (shareId) =>
                await apiCallsMixin.fetchDatasetListInfo({permalink_share_id: shareId})
            );
        } else if (urlParams.has("dataset_id")) {
    		// Legacy support for dataset_id

            await activateDatasetFromParam("dataset_id");
        }
	} catch (error) {
		logErrorInConsole(error);
	}

};

registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);