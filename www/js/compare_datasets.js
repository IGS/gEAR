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
let facetWidget;
let datasetId;
let organismId;	// Used for saving as gene cart
let compareData;;
let selectedGeneData;
let geneSelect;

// Storing user's plot text edits, so they can be restored if user replots
let titleText = null;
let xaxisText = null;
let yaxisText = null;

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

		// Clear selected gene tags
		document.getElementById("gene_tags").replaceChildren();

		// Clear compare groups
		for (const classElt of document.getElementsByClassName("js-compare-groups")) {
			classElt.replaceChildren();
		}

    	// Creates gene select instance that allows for multiple selection
		geneSelect = createGeneSelectInstance("gene_select", geneSelect);
		// Populate gene select element
		await geneSelectUpdate()


		// Create facet widget, which will refresh filters
		facetWidget = await createFacetWidget(datasetId, null, {});
		document.getElementById("facet_content").classList.remove("is-hidden");
		document.getElementById("selected_facets").classList.remove("is-hidden");

		// Update compare series options
		const catColumns = facetWidget.aggregations.map((agg) => agg.name);
		updateSeriesOptions("js-compare", catColumns);

		compareSeriesElt.parentElement.classList.remove("is-loading");

    })
});

const geneCartTree = new GeneCartTree({
    element: document.getElementById("genecart_tree")
    , searchElement: document.getElementById("genecart_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "genecart") {
            return;
        }

        // Get gene symbols from gene cart
        const geneCartId = e.node.data.orig_id;
        const geneCartMembers = await fetchGeneCartMembers(geneCartId);
        const geneCartSymbols = geneCartMembers.map((item) => item.label);

        // Normalize gene symbols to lowercase
        const geneSelectSymbols = geneSelect.data.map((opt) => opt.value);
        const geneCartSymbolsLowerCase = geneCartSymbols.map((x) => x.toLowerCase());

        const geneSelectedOptions = geneSelect.selectedOptions.map((opt) => opt.data.value);

        // Get genes from gene cart that are present in dataset's genes.  Preserve casing of dataset's genes.
        const geneCartIntersection = geneSelectSymbols.filter((x) => geneCartSymbolsLowerCase.includes(x.toLowerCase()));
        // Add in already selected genes (union)
        const geneSelectIntersection = [...new Set(geneCartIntersection.concat(geneSelectedOptions))];

        // change all options to be unselected
        const origSelect = document.getElementById("gene_select");
        for (const opt of origSelect.options) {
            opt.removeAttribute("selected");
        }

        // Assign intersection genes to geneSelect "selected" options
        for (const gene of geneSelectIntersection) {
            const opt = origSelect.querySelector(`option[value="${gene}"]`);
            try {
                opt.setAttribute("selected", "selected");
            } catch (error) {
                // sanity check
                const msg = `Could not add gene ${gene} to gene select.`;
                console.warn(msg);
            }
        }

        geneSelect.update();
        trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to load tags
    })
});

const adjustGeneTableLabels = () => {
    const geneFoldchanges = document.getElementById("tbl_gene_foldchanges");
	const log_base = document.getElementById("log_base").value;

	const spanIcon = document.createElement("span");
	spanIcon.classList.add("icon");
	const i = document.createElement("i");
	i.classList.add("mdi", "mdi-sort-numeric-ascending");
	i.setAttribute("aria-hidden", "true");
	spanIcon.appendChild(i);
	geneFoldchanges.appendChild(spanIcon);

	if (log_base === "raw") {
		geneFoldchanges.prepend("Fold Change ");
		return;
	}
	geneFoldchanges.prepend(`Log${log_base} Fold Change `);
}

const appendGeneTagButton = (geneTagElt) => {
    // Add delete button
    const deleteBtnElt = document.createElement("button");
    deleteBtnElt.classList.add("delete", "is-small");
    geneTagElt.appendChild(deleteBtnElt);
    deleteBtnElt.addEventListener("click", (event) => {
        // Remove gene from geneSelect
        const gene = event.target.parentNode.textContent;
        const geneSelectElt = document.getElementById("gene_select");
        geneSelectElt.querySelector(`option[value="${gene}"]`).removeAttribute("selected");

        geneSelect.update();
        trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to load tags
    });

    // ? Should i add ellipses for too many genes? Should I make the box collapsable?
}

const chooseGene = (event) => {
    // Triggered when a gene is selected

    // Delete existing tags
    const geneTagsElt = document.getElementById("gene_tags");
    geneTagsElt.replaceChildren();

    if (!geneSelect.selectedOptions.length) return;   // Do not trigger after initial population

    // Update list of gene tags
    const sortedGenes = geneSelect.selectedOptions.map((opt) => opt.data.value).sort();
    for (const opt in sortedGenes) {
        const geneTagElt = document.createElement("span");
        geneTagElt.classList.add("tag", "is-primary", "mx-1");
        geneTagElt.textContent = sortedGenes[opt];
        appendGeneTagButton(geneTagElt);
        geneTagsElt.appendChild(geneTagElt);
    }

    document.getElementById("gene_tags_c").classList.remove("is-hidden");
    if (!geneSelect.selectedOptions.length) {
        document.getElementById("gene_tags_c").classList.add("is-hidden");
    }

    // If more than 10 tags, hide the rest and add a "show more" button
    if (geneSelect.selectedOptions.length > 10) {
        const geneTags = geneTagsElt.querySelectorAll("span.tag");
        for (let i = 10; i < geneTags.length; i++) {
            geneTags[i].classList.add("is-hidden");
        }
        // Add show more button
        const showMoreBtnElt = document.createElement("button");
        showMoreBtnElt.classList.add("tag", "button", "is-small", "is-primary", "is-light");
        const numToDisplay = geneSelect.selectedOptions.length - 10;
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
    document.getElementById("clear_genes_btn").classList.add("is-loading");
    geneSelect.clear();
	updatePlotAnnotations([]);
    document.getElementById("clear_genes_btn").classList.remove("is-loading");
}

const createFacetWidget = async (datasetId, analysisId, filters) => {
    document.getElementById("selected_facets_loader").classList.remove("is-hidden")

    const {aggregations, total_count:totalCount} = await fetchAggregations(datasetId, analysisId, filters);
    document.getElementById("num_selected").textContent = totalCount;


    const facetWidget = new FacetWidget({
        aggregations,
        filters,
        onFilterChange: async (filters) => {
            if (filters) {
                try {
                    const {aggregations, total_count:totalCount} = await fetchAggregations(datasetId, analysisId, filters);
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

const createGeneSelectInstance = (idSelector, geneSelect=null) => {
    // NOTE: Updating the list of genes can be memory-intensive if there are a lot of genes
    // and (I've noticed) if multiple select2 elements for genes are present.

    // If object exists, just update it with the revised data and return
    if (geneSelect) {
        geneSelect.update();
        return geneSelect;
    }

    return NiceSelect.bind(document.getElementById(idSelector), {
        placeholder: 'To search, start typing a gene name',
        searchtext: 'To search, start typing a gene name',
        searchable: true,
        allowClear: true,
    });
}

const downloadSelectedGenes = (event) => {
    event.preventDefault();

	// Builds a file in memory for the user to download.  Completely client-side.
	// plot_data contains three keys: x, y and symbols
	// build the file string from this

    // Adjust headers to the plot type
	const xLabel = JSON.stringify([...document.querySelectorAll("#compare_x input:checked")].map((elt) => elt.value));
	const yLabel = JSON.stringify([...document.querySelectorAll("#compare_y input:checked")].map((elt) => elt.value));

	const logBase = document.getElementById("log_base").value;

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
	element.setAttribute("download", "selected_genes.tsv");
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

const fetchDatasets = async () => {
    try {
        return await apiCallsMixin.fetchDatasets();
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch datasets. Please contact the gEAR team."
        createToast(msg);
        throw new Error(msg);
    }
}

/* Fetch gene collection members */
const fetchGeneCartMembers = async (geneCartId) => {
    try {
        const {gene_symbols, success} = await apiCallsMixin.fetchGeneCartMembers(geneCartId);
        if (!success) {
            throw new Error("Could not fetch gene collection members. You can still enter genes manually.");
        }
        return gene_symbols;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene collection members. You can still enter genes manually.";
        createToast(msg);
        throw new Error(msg);
    }
}

/* Fetch gene collections */
const fetchGeneCarts = async () => {
    const cartType = "unweighted-list";
    try {
        return await apiCallsMixin.fetchGeneCarts(cartType);
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene collections. You can still enter genes manually.";
        createToast(msg);
        throw new Error(msg);
    }
}

const fetchGeneSymbols = async (datasetId, analysisId) => {
	try {
		const data = await apiCallsMixin.fetchGeneSymbols(datasetId, analysisId);
		return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
	} catch (error) {
		logErrorInConsole(error);
		const msg = "Could not fetch gene symbols for this dataset. Please contact the gEAR team."
		createToast(msg);
		return [];
	}
}

const geneSelectUpdate = async (analysisId=null) => {
    // Populate gene select element
    try {
        const geneSymbols = await fetchGeneSymbols(datasetId, analysisId);
        updateGeneOptions(geneSymbols); // Come from curator specific code
    } catch (error) {
		logErrorInConsole(error);
	}
}

const getComparisons = async (event) => {

	// set loading icon
	event.target.classList.add("is-loading");

	const filters = JSON.stringify(facetWidget.filters);

	const compareSeries = document.getElementById("compare_series").value

	// Get all checked x and y series
	const checkedX = JSON.stringify([...document.querySelectorAll("#compare_x input:checked")].map((elt) => elt.value));
	const checkedY = JSON.stringify([...document.querySelectorAll("#compare_y input:checked")].map((elt) => elt.value));

	const foldChangeCutoff = document.getElementById("fc_cutoff").value;
	const stdDevNumCutoff = document.getElementById("standard_deviation").value;
	const logTransformation = document.getElementById("log_base").value;
	const statisticalTest = document.getElementById("statistical_test").value;

	try {
		const data = await fetchDatasetComparison(datasetId, filters, compareSeries, checkedX, checkedY, foldChangeCutoff, stdDevNumCutoff, logTransformation, statisticalTest);
		if (data?.success < 1) {
			throw new Error(data?.message || "Could not fetch dataset comparison. Please contact the gEAR team.");
		}
		compareData = data;
		plotDataToGraph(compareData);

		// If any genes selected, update plot annotations (since plot was previously purged)
		const sortedGenes = geneSelect.selectedOptions.map((opt) => opt.data.value).sort();
		updatePlotAnnotations(sortedGenes);

        // Show button to add genes to gene cart
        document.getElementById("gene_cart_btn_c").classList.remove("is-hidden");

		// Hide this view
		document.getElementById("content_c").classList.add("is-hidden");
		// Generate and display "post-plotting" view/container
		document.getElementById("post_plot_content_c").classList.remove("is-hidden");

	} catch (error) {
		console.error(error);
		handleGetComparisonError(datasetId, checkedX, checkedY);
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

		if (elt == "standard_deviation" && !(value === "0" )) {
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
        const datasetData = await fetchDatasets();

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

/* Transform and load gene collection data into a "tree" format */
const loadGeneCarts = async () => {
    try {
        const geneCartData = await fetchGeneCarts();
        const carts = {};
        const cartTypes = ['domain', 'user', 'group', 'shared', 'public'];
        let cartsFound = false;

        // Loop through the different types of gene collections and add them to the carts object
        for (const ctype of cartTypes) {
            carts[ctype] = [];

            if (geneCartData[`${ctype}_carts`].length > 0) {
                cartsFound = true;

                for (const item of geneCartData[`${ctype}_carts`]) {
					const fullLabel = `${item.label} (${item.num_genes} genes)`; // Add number of genes to label
                    carts[ctype].push({value: item.id, text: fullLabel });
                };
            }
        }

        geneCartTree.domainGeneCarts = carts.domain;
        geneCartTree.userGeneCarts = carts.user;
        geneCartTree.groupGeneCarts = carts.group;
        geneCartTree.sharedGeneCarts = carts.shared;
        geneCartTree.publicGeneCarts = carts.public;
        geneCartTree.generateTree();
        /*if (!cartsFound ) {
            // ? Put some warning if carts not found
            $('#gene_cart_container').show();
        }*/

    } catch (error) {
        document.getElementById("gene_s_failed").classList.remove("is-hidden");
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

	const plotContainer = document.getElementById("plot_container");
	plotContainer.replaceChildren();    // erase plot

	// NOTE: Plot initially is created to a default width but is responsive.
	// Noticed container within our "column" will make full-width go beyond the screen

	const plotlyPreview = document.createElement("div");
	plotlyPreview.id = "plotly_preview";
	plotlyPreview.classList.add("container", "is-max-desktop");
	plotContainer.append(plotlyPreview);

	Plotly.purge("plotly_preview"); // clear old Plotly plots

	Plotly.newPlot("plotly_preview", plotData, layout, config);

	// If plot data is selected, create the right-column table and do other misc things
	plotlyPreview.on("plotly_selected", (eventData) => {

		// Hide selected genes table and disable unweighted radio button if no genes are selected
		document.getElementById("tbl_selected_genes").classList.add("is-hidden");
		document.getElementById("download_selected_genes_btn").classList.add("is-hidden");
		document.querySelector("input[name='genecart_type'][value='unweighted']").disabled = true;
		document.querySelector("input[name='genecart_type'][value='unweighted']").parentElement.setAttribute("disabled", "disabled");

		if (eventData?.points.length) {
			document.getElementById("tbl_selected_genes").classList.remove("is-hidden");
			document.getElementById("download_selected_genes_btn").classList.remove("is-hidden");
			document.querySelector("input[name='genecart_type'][value='unweighted']").disabled = false;
			document.querySelector("input[name='genecart_type'][value='unweighted']").parentElement.removeAttribute("disabled");

			adjustGeneTableLabels();
			populateGeneTable(eventData);
		}

		// Get genes from gene tags
		const geneTags = document.querySelectorAll("#gene_tags span.tag");
		const searchedGenes = [];
		for (const tag of geneTags) {
			searchedGenes.push(tag.textContent);
		}

		// Highlight table rows that match searched genes
		if (searchedGenes) {
			const geneTableBody = document.getElementById("gene_table_body");
			// Select the first column (gene_symbols) in each row
			for (const row of geneTableBody.children) {
				const tableGene = row.children[0].textContent;
				for (const gene of searchedGenes) {
					if (gene.toLowerCase() === tableGene.toLowerCase() ) {
						row.classList.add("has-background-success-light");
					}
				};
			}
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
	plotlyNote.id = "tip_on_editing";
	plotlyNote.classList.add("notification", "content", "is-info", "is-light");
	plotlyNote.innerHTML = `<p><strong>Tip:</strong> Use the Plotly box and lasso select tools (upper-right) to select genes to view as a table.</p>

	<p>You can also click the plot title or axis labels to edit them. Hit Enter to apply edit.</p>`;
	plotlyPreview.append(plotlyNote);
}


const populateGeneTable = (data) => {
	const statisticalTest = document.getElementById("statistical_test").value;

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


    const geneTableBody = document.getElementById("gene_table_body");
    geneTableBody.replaceChildren();

    for (const gene of selectedGeneData) {
        const row = document.createElement("tr");
        row.innerHTML = `<td>${gene.gene_symbol}</td><td>${gene.pval}</td><td>${gene.foldchange}</td>`;
        geneTableBody.appendChild(row);
    }

	// If not statistical test, delete p-value column
	if (!statisticalTest) {
		const pvalColumn = document.querySelector("#tbl_selected_genes thead tr th:nth-child(2)");
		pvalColumn.remove();
		for (const pvalCell of document.querySelectorAll("#tbl_selected_genes tbody tr td:nth-child(2)")) {
			pvalCell.remove();
		}
	}
	// Should be sorted by logFC now

}

const populatePostCompareBox = (scope, series, groups) => {
	// Find box
	const boxElt = document.querySelector(`#${scope}_post_c .notification`);
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
        , label: document.getElementById("new_genecart_label").value
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
		, label:  document.getElementById("new_genecart_label").value
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

const updateGeneOptions = (geneSymbols) => {

    const geneSelectElt = document.getElementById("gene_select");
    geneSelectElt.replaceChildren();

	geneSelectElt.parentElement.classList.add("is-loading");

    // Append empty placeholder element
    const firstOption = document.createElement("option");
    firstOption.textContent = "Please select a gene";
    geneSelectElt.append(firstOption);

    for (const gene of geneSymbols.sort()) {
        const option = document.createElement("option");
        option.textContent = gene;
        option.value = gene;
        geneSelectElt.append(option);
    }

    // Update the nice-select2 element to reflect this.
    // This function is always called in the 1st view, so only update that
    geneSelect.update();

	geneSelectElt.parentElement.classList.remove("is-loading");

}

// For a given categorical series (e.g. "celltype"), add checkboxes for each category
const updateGroupOptions = (selectorId, groupsArray, series) => {

	const elt = document.getElementById(selectorId);
	elt.classList.remove("is-hidden");

	// Add categories
	for (const group of groupsArray.sort()) {

		const checkbox = document.createElement("input");
		checkbox.type = "checkbox";
		checkbox.id = `${selectorId}_${group}`;
		checkbox.name = group;
		checkbox.value = group;

		const label = document.createElement("label");
		label.classList.add("checkbox");
		label.htmlFor = `${selectorId}_${group}`;
		label.textContent = ` ${group}`;
		label.prepend(checkbox);

		// If group has aggregations count of 0 (no data after filtering), disable checkbox and label
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

	const plotlyPreview = document.getElementById("plotly_preview");
	const plotData = plotlyPreview.data;
	const layout = plotlyPreview.layout;

	const annotationColor = CURRENT_USER.colorblind_mode ? "orange" : "cyan";

	layout.annotations = [];

	// Reset all trace colors
	for (const trace of plotData) {
		trace.marker.color = trace.marker.origColor;
	}

	genes.forEach((gene) => {
		let found = false;
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

				found = true;
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

const updateUIAfterGeneCartSaveSuccess = (gc) => {
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

			document.getElementById("condition_compare_s_failed").classList.add("is-hidden");
		}
		return;
	}

    document.getElementById("condition_compare_s_success").classList.add("is-hidden");
}

/* --- Event listeners --- */

document.getElementById("statistical_test").addEventListener("change", (event) => {
	const pvalCutoff = document.getElementById("pval_cutoff");
	const cutoffFilterAction = document.getElementById("cutoff_filter_action");
	pvalCutoff.disabled = event.target.value ? false : true;
	cutoffFilterAction.disabled = event.target.value ? false : true;
});

// When compare series changes, update the compare groups
for (const classElt of document.getElementsByClassName("js-compare")) {
	const compareSeriesNotification = document.getElementById("select_compare_series_notification");
	classElt.addEventListener("change", async (event) => {
		const compareSeries = event.target.value;
		compareSeriesNotification.classList.remove("is-hidden", "is-danger");
		compareSeriesNotification.classList.add("is-warning");
		compareSeriesNotification.textContent = "Please select a series to compare first";

		for (const classElt of document.getElementsByClassName("js-compare-groups")) {
			classElt.classList.add("is-hidden");
			classElt.replaceChildren();
		}

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

		updateGroupOptions("compare_x", seriesNames, compareSeries);
		updateGroupOptions("compare_y", seriesNames, compareSeries);
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
		const compareSeries = document.getElementById("compare_series").value;
		populatePostCompareBox("compare_x", compareSeries, checkedX);
	})
}

for (const classElt of document.getElementsByClassName("js-compare-y")) {
	classElt.addEventListener("change", (event) => {
		const checkedY = [...document.querySelectorAll(".js-compare-y input:checked")].map((elt) => elt.value);
		const compareSeries = document.getElementById("compare_series").value;
		populatePostCompareBox("compare_y", compareSeries, checkedY);
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

document.getElementById("clear_genes_btn").addEventListener("click", clearGenes);

const geneSelectElts = document.querySelectorAll("select.js-gene-select");
for (const geneSelectElt of geneSelectElts) {
    geneSelectElt.addEventListener("change", chooseGene);
}

// code from Bulma documentation to handle modals
document.getElementById("gene_cart_btn").addEventListener("click", ($trigger) => {
    const closestButton = $trigger.target.closest(".button");
    const modal = closestButton.dataset.target;
    const $target = document.getElementById(modal);
    openModal($target);

});

document.getElementById("new_genecart_label").addEventListener("input", (event) => {
    const saveBtn = document.getElementById("save_genecart_btn");
    saveBtn.disabled = event.target.value ? false : true;
});

document.getElementById("save_genecart_btn").addEventListener("click", (event) => {
    event.preventDefault();
    event.target.classList.add("is-loading");
    // get value of genecart radio button group
    const geneCartName = document.querySelector("input[name='genecart_type']:checked").value;
    if (CURRENT_USER) {
        if (geneCartName === "unweighted") {
            saveGeneCart();
        } else {
            saveWeightedGeneCart();
        }
    }
    event.target.classList.remove("is-loading");
});

document.getElementById("download_selected_genes_btn").addEventListener("click", downloadSelectedGenes);

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// Update with current page info
	document.getElementById("page-header-label").textContent = "Comparison Tool";
	for (const elt of document.querySelectorAll("#primary_nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}

	document.querySelector("a[tool='compare'").classList.add("is-active");

    sessionId = CURRENT_USER.session_id;

	if (! sessionId ) {
		// TODO: Add master override to prevent other triggers from enabling saving
        createToast("Not logged in so saving gene carts is disabled.");
        document.getElementById("gene_cart_btn").disabled = true;
    }


	try {
		await Promise.all([
			loadDatasetTree(),
			loadGeneCarts()
		]);
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