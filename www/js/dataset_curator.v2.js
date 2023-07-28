// I use camelCase for my variable/function names to adhere to JS style standards
// Exception being functions that do fetch calls, so we can use JS destructuring on the payload

'use strict';

//TODO  - Create "error" reporting in each step

let plotConfig = {};  // Plot config that is passed to API or stored in DB

let datasetId = null;
let organismId = null;

let analysisSelect = null;
let plotSelect = null;
let geneSelect = null;

const userId = 622;     // ! It's me
const sessionId = "ee95e48d-c512-4083-bf05-ca9f65e2c12a"    // ! It's me
const session_id = sessionId;

Cookies.set('gear_session_id', sessionId, { expires: 7 });

const datasetTree = new DatasetTree({
    element: document.getElementById("dataset_tree")
    , searchElement: document.getElementById("dataset_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "dataset") {
            return;
        }
        document.getElementById("current_dataset_c").style.display = "";
        document.getElementById("current_dataset").textContent = e.node.title;
        const newDatasetId = e.node.data.dataset_id;
        organismId = e.node.data.organism_id;

        // We don't want to needless run this if the same dataset was clicked
        if (newDatasetId === datasetId) {
            return;
        }

        // Click to get to next step
        document.getElementById("load_plot_s").click();
        document.getElementById('new_display').disabled = false;

        // Clear "current <whatever>" text
        document.getElementById("current_gene_c").style.display = "none";

        // Clear (and update) options within nice-select2 structure.
        analysisSelect.clear();
        geneSelect.clear(); // BUG: Figure out why this is triggering twice (check function in Github)

        datasetId = newDatasetId;
        // Fetch dataset information
        const {owner_id: ownerId} = await fetchDatasetInfo(datasetId);

        // displays
        const userDisplays = await fetchDatasetDisplays(userId, datasetId);
        const ownerDisplays = userId === ownerId ? [] : await fetchDatasetDisplays(ownerId, datasetId);
        const defaultDisplayId = await fetchDefaultDisplay(userId, datasetId);
        renderDisplayCards(userDisplays, ownerDisplays, defaultDisplayId);

        // analyses
        const {public_analyses: publicAnalyses, private_analyses: privateAnalyses} = await fetchAnalyses(datasetId);
        updateAnalysesOptions(privateAnalyses, publicAnalyses);

        // plot types
        const availablePlotTypes = await fetchAvailablePlotTypes(userId, sessionId, datasetId, undefined);
        for (const plotType in availablePlotTypes) {
            const isAllowed = availablePlotTypes[plotType];
            if (plotType === "tsne/umap_dynamic") {
                document.getElementById("tsne_dyna_opt").disabled = !isAllowed;
            } else {
                document.getElementById(`${plotType}_opt`).disabled = !isAllowed;
            }
        }
    })
});

const chooseAnalysis = async () => {

    const analysisId = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;

    // Populate gene select element
    const geneSymbols = await fetchGeneSymbols(datasetId, analysisId);
    updateGeneSymbolOptions(geneSymbols);
}

/* Display has been chosen, so display analysis and plot type options */
const chooseDisplay = async () => {
    document.getElementById("analysis_type_select").disabled = false;
    analysisSelect.update();

    document.getElementById("plot_type_select").disabled = false;
    plotSelect.update();

    // Populate gene select element
    const geneSymbols = await fetchGeneSymbols(datasetId, undefined);
    updateGeneSymbolOptions(geneSymbols);

    document.getElementById("plot_type_s").click();
}

const chooseGene = () => {
    if (!geneSelect.selectedOptions.length) return;   // Do not trigger after initial population

    document.getElementById("current_gene_c").style.display = "";
    document.getElementById("current_gene").textContent = getSelect2Value(geneSelect);

    document.getElementById("plot_options_s").click();
}

const choosePlot = () => {
    if (!plotSelect.selectedOptions.length) return;   // Do not trigger after setting disable/enable on options

    includePlotParamOptions();
    document.getElementById("gene_s").click();
}

const cloneDisplay = async (display) => {
    // Open the next step
    await chooseDisplay();

    // Read clone config to populate analysis, plot type, gnee and plot-specific options
    // TODO: Complete
    plotConfig = display.plotConfig;
    document.getElementById("gene_select").value = plotConfig.gene_symbol
    geneSelect.update();

}

const createAnalysisSelectInstance = () => {
    analysisSelect = NiceSelect.bind(document.getElementById("analysis_type_select"), {
        placeholder: 'Select an analysis.',
        allowClear: true,
    });
}

const createGeneSelectInstance = () => {
    geneSelect = NiceSelect.bind(document.getElementById("gene_select"), {
        placeholder: 'To search, start typing a gene name',
        searchtext: 'To search, start typing a gene name',
        searchable: true,
        allowClear: true,
    });

    // BUG: Figure out "scrollbar"
}

const createPlotSelectInstance = () => {
    // Initialize plot types
    plotSelect = NiceSelect.bind(document.getElementById("plot_type_select"), {
        placeholder: 'Choose how to plot',
        minimumResultsForSearch: -1
    });
}

const deleteDisplay = async(user_id, display_id) => {
    const payload = {user_id, display_id};
    return await axios.post("/cgi/delete_dataset_display.cgi", payload);
}

const fetchAnalyses = async (datasetId) => {
    const { data } = await axios.get(`./api/h5ad/${datasetId}/analyses`);

    const { public: public_analyses, private: private_analyses } = data;
    return {public_analyses, private_analyses};

}

const fetchAvailablePlotTypes = async (user_id, session_id, dataset_id, analysis_id) => {
    const payload = {user_id, session_id, dataset_id, analysis_id};
    const {data} = await axios.post(`/api/h5ad/${dataset_id}/availableDisplayTypes`, payload);
    return data
}

const fetchDatasetDisplay = async (display_id) => {
    const payload = {display_id};
    // POST due to payload variables being sensitive
    const {data} = await axios.post("/cgi/get_dataset_display.cgi", convertToFormData(payload));
    return data
}

const fetchDatasetDisplayImage = async (dataset_id, display_id) => {
    const payload = {dataset_id, display_id};
    // POST due to payload variables being sensitive
    const {data} = await axios.post("/cgi/get_dataset_display_image.cgi", convertToFormData(payload));
    return data
}

const fetchDatasetDisplays = async (user_id, dataset_id) => {
    const payload = {user_id, dataset_id};
    // POST due to payload variables being sensitive
    const {data} = await axios.post("/cgi/get_dataset_displays.cgi", convertToFormData(payload));
    return data;
}

const fetchDatasetInfo = async (dataset_id) => {
    const payload = {dataset_id};
    const {data} = await axios.post("/cgi/get_dataset_info.cgi", convertToFormData(payload));
    const {title, is_public, owner_id} = data;
    return {title, is_public, owner_id};
}

const fetchDatasets = async (session_id) => {
    const payload = {session_id}
    const {data} = await axios.post("cgi/get_h5ad_dataset_list.cgi", convertToFormData(payload));
    return data;
}

const fetchDefaultDisplay = async (user_id, dataset_id) => {
    const payload = {user_id, dataset_id};
    // POST due to payload variables being sensitive
    const {data} =  await axios.post("/cgi/get_default_display.cgi", convertToFormData(payload));
    const {default_display_id} = data;
    return default_display_id;
}

const fetchGeneSymbols = async (datasetId, analysisId) => {
    let url = `./api/h5ad/${datasetId}/genes`;
    if (analysisId) url += `?analysis_id=${analysisId}`;

    const { data } = await axios.get(url);
    return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
}

const fetchH5adInfo = async (datasetId, analysisId) => {
    let url = `/api/h5ad/${datasetId}`
    if (analysisId) url += `?analysis_id=${analysisId}`;
    const {data} = await axios.get(url);
    const { obs_columns, obs_levels } = data;
    return { obs_columns, obs_levels };
}

const fetchPlotlyData = async (plotConfig, datasetId, plot_type, colorblind_mode)  => {
    // TODO: Set loading
    const payload = { ...plotConfig,  plot_type, colorblind_mode };
    const { data } = await axios.post(`/api/plot/${datasetId}`, payload);
    // TODO: Set plot
}

const fetchSvgData = async (geneSymbol, datasetId) => {
    const { data } = await axios.get(`/api/plot/${datasetId}/svg?gene=${geneSymbol}`);
    //TODO: Set plot
};

const fetchTsneImage = async (plotConfig, datasetId, plot_type, analysis, analysis_owner_id, colorblind_mode) => {
    // TODO: Set loading
    const payload = { ...plotConfig, plot_type, analysis, analysis_owner_id, colorblind_mode };
    const { data } = await axios.post(`/api/plot/${datasetId}/tsne`, payload);
    // TODO: set plot
}

const getSelect2Value = (select) => {
    // Get value from select2 element
    return select.selectedOptions[0].data.value;
}

const includePlotParamOptions = async () => {
    const plotType = getSelect2Value(plotSelect);

    const plotlyPlots = ["bar", "line", "scatter", "tsne_dyna", "violin"];
    const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];

    const plotOptionsElt = document.getElementById("plot_options_collapsable");
    plotOptionsElt.replaceChildren();

    // include plotting backend options
    if (plotlyPlots.includes(plotType)) {
        const response = await fetch("../include/plot_types/single_gene_plotly.html", {cache: "reload"});
        const body = await response.text();
        plotOptionsElt.innerHTML = body;

        // populate advanced options for specific plot types
        const plotSpecificOptionsElt = document.getElementById("plot_specific_options");
        if (["scatter", "tsne_dyna"].includes(plotType)) {
            const response = await fetch("../include/plot_types/advanced_scatter.html", {cache: "reload"});
            const body = await response.text();
            plotSpecificOptionsElt.innerHTML = body;
        } else if (plotType == "violin") {
            const response = await fetch("../include/plot_types/advanced_violin.html", {cache: "reload"});
            const body = await response.text();
            plotSpecificOptionsElt.innerHTML = body;
        }

        setupPlotlyOptions();
        return;
    }
    if (scanpyPlots.includes(plotType)) {
        const response = await fetch("../include/plot_types/tsne_static.html", {cache: "reload"});
        const body = await response.text();
        plotOptionsElt.innerHTML = body;

        setupScanpyOptions();
        return;
    }
    if (plotType === "svg") {
        const response = await fetch("../include/plot_types/svg.html", {cache: "reload"});
        const body = await response.text();
        plotOptionsElt.innerHTML = body;
        return;
    }

    console.warn(`Plot type ${plotType} not recognized.`)

}

/* Transform and load dataset data into a "tree" format */
const loadDatasetTree = async () => {

    const datasetData = await fetchDatasets(sessionId);

    let counter = 0;

    // Populate select box with dataset information owned by the user
    const userDatasets = [];
    if (datasetData.user.datasets.length > 0) {
        // User has some profiles
        datasetData.user.datasets.forEach((item) => {
            if (item) {
                userDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
            }
        });
    }
    // Next, add datasets shared with the user
    const sharedDatasets = [];
    if (datasetData.shared_with_user.datasets.length > 0) {
        datasetData.shared_with_user.datasets.forEach((item) => {
            if (item) {
                sharedDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
            }
        });
    }
    // Now, add public datasets
    const domainDatasets = [];
    if (datasetData.public.datasets.length > 0) {
        datasetData.public.datasets.forEach((item) => {
            if (item) {
                domainDatasets.push({ value: counter++, text: item.title, dataset_id : item.id, organism_id: item.organism_id });
            }
        });
    }

    datasetTree.userDatasets = userDatasets;
    datasetTree.sharedDatasets = sharedDatasets;
    datasetTree.domainDatasets = domainDatasets;
    datasetTree.generateTree();
}

const renderDisplayCards = async (userDisplays, ownerDisplays, defaultDisplayId) => {
    const userDisplaysElt = document.getElementById("user_displays");
    const ownerDisplaysElt = document.getElementById("owner_displays");

    // Empty existing displays
    userDisplaysElt.replaceChildren();
    ownerDisplaysElt.replaceChildren();

    for (const display of userDisplays) {
        const displayUrl = await fetchDatasetDisplayImage(datasetId, display.id);
        const geneSymbol = display.plotly_config.gene_symbol;

        const label = display.label || "";

        // TODO - Get footer button styles correct
        const template = `
                    <div id="${display.id}_display" class="column is-one-quarter">
                        <div class="box card has-background-primary-light has-text-primary">
                            <header class="card-header">
                                <p class="card-header-title has-text-black">${label}</p>
                            </header>
                            <div class="card-image">
                                <figure class="image is-4by3">
                                    <img src="${displayUrl} " alt="Saved display">
                                </figure>
                            </div>
                            <div class="card-content">
                                <p class="subtitle">Gene: ${geneSymbol}</p>
                            </div>
                            <footer class="card-footer">
                                <button class="js-display-default card-footer-item is-link" id="${display.id}_default">Set as Default</button>
                                <button class="card-footer-item button is-link" id="${display.id}_clone">Clone</button>
                                <button class="card-footer-item button is-link" id="${display.id}_delete">Delete</button>
                            </footer>
                        </div>
                    </div>`;

        const htmlCollection = generateElements(template);
        userDisplaysElt.append(htmlCollection);

        // Edit default properties if this is the default display
        const defaultElt = document.getElementById(`${display.id}_default`);
        if (display.id === defaultDisplayId) {
            defaultElt.textContent = "Default";
            defaultElt.disabled = true;
        }

        defaultElt.addEventListener("click", (event) => saveDefaultDisplay(display.id));
        document.getElementById(`${display.id}_clone`).addEventListener("click", (event) => cloneDisplay(display));
        document.getElementById(`${display.id}_delete`).addEventListener("click", (event) => deleteDisplay(userId, display.id));

    }

    for (const display of ownerDisplays) {
        const displayUrl = await fetchDatasetDisplayImage(datasetId, display.id);
        const geneSymbol = display.plotly_config.gene_symbol;

        const label = display.label || `Unnamed ${display.plot_type} display`;

        const template = `
                    <div id="${display.id}_display" class="column is-one-quarter">
                        <div class="box card has-background-primary-light has-text-primary">
                            <header class="card-header">
                                <p class="card-header-title has-text-black">${label}</p>
                            </header>
                            <div class="card-image">
                                <figure class="image is-4by3">
                                    <img src="${displayUrl} " alt="Saved display">
                                </figure>
                            </div>
                            <div class="card-content">
                                <p class="subtitle is-6">Gene: ${geneSymbol}</p>
                            </div>
                            <footer class="card-footer">
                                <button class="js-display-default card-footer-item button is-link" id="${display.id}_default">Set as Default</button>
                                <button class="card-footer-item button is-link" id="${display.id}_clone">Clone</button>
                            </footer>
                        </div>
                    </div>`;

        const htmlCollection = generateElements(template);
        ownerDisplaysElt.append(htmlCollection);

        // Edit default properties if this is the default display
        const defaultElt = document.getElementById(`${display.id}_default`);
        if (display.id === defaultDisplayId) {
            defaultElt.textContent = "Default";
            defaultElt.disabled = true;
        }

        defaultElt.addEventListener("click", (event) => saveDefaultDisplay(display.id));
        document.getElementById(`${display.id}_clone`).addEventListener("click", (event) => cloneDisplay(display));

    }
}

const saveDatasetDisplay = async(displayId, dataset_id, user_id, label, plot_type, plotConfig) => {
    // NOTE: Saving all displays as new displays (clone) instead of overwriting. User can always delete excess displays
    const payload = {
        id: displayId,
        dataset_id,
        user_id,
        label,
        plot_type,
        plotly_config: JSON.stringify({
        // depending on display type, this object will
        // have different properties
        ...plotConfig,
        }),
    };

    return await axios.post("/cgi/save_dataset_display.cgi", convertToFormData(payload));
}

const saveDefaultDisplay = async (displayId) => {
    const payload = {display_id: displayId};
    await axios.post("/cgi/save_default_display.cgi", convertToFormData(payload))
    .catch(function (error) {
        if (error.response) {
            // The request was made and the server responded with a status code
            // that falls out of the range of 2xx
            console.error(error.response.data);
            console.error(error.response.status);
            console.error(error.response.headers);
        } else if (error.request) {
            // The request was made but no response was received
            // `error.request` is an instance of XMLHttpRequest in the browser and an instance of
            // http.ClientRequest in node.js
            console.error(error.request);
        } else {
            // Something happened in setting up the request that triggered an Error
            console.error('Error', error.message);
        }
        console.error(error.config);
    });

    //Update labels of displays... this becomes "Default", others become "Make Default"
    for (const elt of document.getElementsByClassName("js-display-default")) {
        elt.disabled = false;
        elt.textContent = "Set as Default";
    }

    const currentDefaultElt = document.getElementById(`${displayId}_default`);
    currentDefaultElt.disabled = true;
    currentDefaultElt.textContent = "Default";
}

/* Set up the plotly-based plot options, such as "select" elements, events, etc. */
const setupPlotlyOptions = async () => {

    const analysisId = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const plotType = getSelect2Value(plotSelect);
    const {obs_columns: allColumns, obs_levels: levels} = await fetchH5adInfo(datasetId, analysisId);
    const catColumns = Object.keys(levels);

    updateSeriesOptions("x_axis_series", allColumns, true);
    updateSeriesOptions("y_axis_series", allColumns, true, "raw_value");
    updateSeriesOptions("color_series", allColumns, true);
    updateSeriesOptions("label_series", allColumns, true);
    updateSeriesOptions("facet_row_series", catColumns, false);
    updateSeriesOptions("facet_col_series", catColumns, false);

    if (["scatter", "tsne_dyna"].includes(plotType)) {
        const difference = (arr1, arr2) => arr1.filter(x => !arr2.includes(x))
        const continuousColumns = difference(allColumns, catColumns);

        updateSeriesOptions("marker_size_series", continuousColumns, true);
    }

    // TODO: Set up validation checkers
}

/* Set up the scanpy-based plot options, such as "select" elements, events, etc. */
const setupScanpyOptions = async () => {
    const analysisId = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const plotType = getSelect2Value(plotType);
    const {allColumns, catColumns} = await fetchH5adInfo(datasetId, analysisId);

    let xDefaultOption = null;
    let yDefaultOption = null;

    if (analysisId) {
        switch (plotType) {
            case "pca_static":
                xDefaultOption = "X_pca_1";
                yDefaultOption = "X_pca_2";
                break;
            case "tsne_static":
                xDefaultOption = "X_tsne_1";
                yDefaultOption = "X_tsne_2";
                break;
            case "umap_static":
                xDefaultOption = "X_umap_1";
                yDefaultOption = "X_umap_2";
                break;
        }
    }

    updateSeriesOptions("x_axis_series", allColumns, true, xDefaultOption);
    updateSeriesOptions("y_axis_series", allColumns, true, yDefaultOption);
    updateSeriesOptions("plot_by_series", catColumns, false);

    // TODO: set up validation checkers
}

const updateAnalysesOptions = (privateAnalyses, publicAnalyses) => {
    const privateAnalysesElt = document.getElementById("private_analyses");
    const publicAnalysesElt = document.getElementById("public_analyses");

    // Empty the old optgroups
    privateAnalysesElt.replaceChildren();
    publicAnalysesElt.replaceChildren();

    // Show message that no analyses are present if none exist
    if (!privateAnalyses.length) {
        const option = document.createElement("option");
        option.text = "You have no saved analyses for this dataset.";
        option.disabled = true;
        privateAnalysesElt.append(option);
    }

    if (!publicAnalyses.length) {
        const option = document.createElement("option");
        option.text = "There are no public analyses for this dataset.";
        option.disabled = true;
        publicAnalysesElt.append(option);
    }

    // Load each analysis as an option
    for (const analysis of privateAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.dataset.id = analysis.id;
        option.dataset.type = analysis.type;
        privateAnalysesElt.append(option);
    }
    for (const analysis of publicAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.dataset.id = analysis.id;
        option.dataset.type = analysis.type;
        publicAnalysesElt.append(option);
    }

    // Update seel
    analysisSelect.update();
}

const updateGeneSymbolOptions = (geneSymbols) => {
    const geneSelectElt = document.getElementById("gene_select");

    geneSelectElt.replaceChildren();

    // Append empty placeholder element
    const firstOption = document.createElement("option");
    geneSelectElt.append(firstOption);

    for (const gene of geneSymbols.sort()) {
        const option = document.createElement("option");
        option.textContent = gene;
        option.value = gene;
        geneSelectElt.append(option);
    }

    // Update the nice-select2 element to reflect this.
    geneSelect.update();
}

// For plotting options, populate select menus with category groups
const updateSeriesOptions = (element, seriesArray, addExpression, defaultOption) => {
    const elt = document.getElementById(element);

    elt.replaceChildren();

    // Append empty placeholder element
    const firstOption = document.createElement("option");
    elt.append(firstOption);

    // Add an expression option (since expression is not in the categories)
    if (addExpression) {
        const expression = document.createElement("option");
        expression.textContent = "Expression";
        expression.value = "raw_value";
        if ("raw_value" === defaultOption) {
            expression.selected = true;
        }
    }

    // Add categories
    for (const group of seriesArray.sort()) {
        const option = document.createElement("option");
        option.textContent = group;
        option.value = group;
        elt.append(option);
        if (group === defaultOption) {
            option.selected = true;
        }
    }

}

window.onload = () => {

    loadDatasetTree();

    // If brought here by the "gene search results" page, curate on the dataset ID that referred us
    /*const linkedDatasetId = getUrlParameter("dataset_id");
    if (linkedDatasetId) {
        try {
            // find DatasetTree node and trigger "activate"
            const foundNode = datasetTree.findFirst(e => e.data.dataset_id === linkedDatasetId);
            foundNode.setActive(true);
            datasetId = linkedDatasetId;
            $("#dataset").val(linkedDatasetId);
        } catch {
            console.error(`Dataset id ${linkedDatasetId} was not returned as a public/private/shared dataset`);
        }
    }*/

    createAnalysisSelectInstance();
    createPlotSelectInstance();
    createGeneSelectInstance();
};

document.getElementById("new_display").addEventListener("click", chooseDisplay);
document.getElementById("analysis_type_select").addEventListener("change", chooseAnalysis);
document.getElementById("plot_type_select").addEventListener("change", choosePlot);
document.getElementById("gene_select").addEventListener("change", chooseGene);

