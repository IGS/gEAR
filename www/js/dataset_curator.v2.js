// I use camelCase for my variable/function names to adhere to JS style standards
// Exception being functions that do fetch calls, so we can use JS destructuring on the payload

'use strict';

const session_id = 'ee95e48d-c512-4083-bf05-ca9f65e2c12a';

//TODO - Create "error" reporting in each step

let numObs = 0; // dummy value to initialize with

let obsFilters = {};

let plotConfig = {};  // Plot config that is passed to API or stored in DB
let plotData = null;    // Plotly "data" JSON

let datasetId = null;
let organismId = null;
let displayId = null;
let obsLevels = null;
let obsNotUsed = null;
let geneSymbol = null;

const userId = 622;     // ! It's me
const sessionId = "ee95e48d-c512-4083-bf05-ca9f65e2c12a"
Cookies.set('gear_session_id', sessionId, { expires: 7 });

const datasetTree = new DatasetTree({
    element: document.getElementById("dataset_tree")
    , searchElement: document.getElementById("dataset_query")
    , selectCallback: (async (e) => {
        if (e.node.type === "dataset") {

            document.getElementById("current_dataset_c").style.display = "";
            document.getElementById("current_dataset").textContent = e.node.title
            const newDatasetId = e.node.data.dataset_id;
            organismId = e.node.data.organism_id;

            // We don't want to needless run this if the same dataset was clicked
            if (newDatasetId !== datasetId) {
                datasetId = newDatasetId;
                // Fetch dataset information, such as analysis, genes, and AnnData stuff
                const {owner_id: ownerId} = await fetchDatasetInfo(datasetId);
                const userDisplays = await fetchDatasetDisplays(userId, datasetId);
                const ownerDisplays = userId === ownerId ? [] : await fetchDatasetDisplays(ownerId, datasetId);

                const analysesPromise = fetchAnalyses(datasetId);

                const defaultDisplayPromise = fetchDefaultDisplay(userId, datasetId);
                const geneSymbolsPromise = fetchGeneSymbols(datasetId, undefined);

                const {public_analyses: publicAnalyses, private_analyses: privateAnalyses} = await analysesPromise;

                const defaultDisplayId = await defaultDisplayPromise.data;

                for (const display of userDisplays) {
                    const displayUrl = await fetchDatasetDisplayImage(datasetId, display.id);
                    const geneSymbol = display.plotly_config.gene_symbol;

                    const label = display.label || "";

                    const defaultElt = (display.id === defaultDisplayId)
                        ? '<a href="#" class="card-footer-item">Default</a>'
                        : '<a href="#" class="card-footer-item">Set as Default</a>';

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
                                ${defaultElt}
                                <a href="#" class="card-footer-item">Clone</a>
                                <a href="#" class="card-footer-item">Delete</a>
                            </footer>
                        </div>
                    </div>`

                    const parent = document.getElementById("user_displays");
                    const htmlCollection = generateElements(template);
                    parent.append(htmlCollection)
                }

                for (const display of ownerDisplays) {
                    const displayUrl = await fetchDatasetDisplayImage(datasetId, display.id);
                    const geneSymbol = display.plotly_config.gene_symbol;

                    const label = display.label || "";

                    const defaultElt = (display.id === defaultDisplayId)
                        ? '<a href="#" class="card-footer-item">Default</a>'
                        : '<a href="#" class="card-footer-item">Set as Default</a>';

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
                                ${defaultElt}
                                <a href="#" class="card-footer-item">Clone</a>
                            </footer>
                        </div>
                    </div>`

                    const parent = document.getElementById("owner_displays");
                    const htmlCollection = generateElements(template);
                    parent.append(htmlCollection)
                }


                // Populate gene select element
                const geneSymbols = await geneSymbolsPromise;
                updateGeneSymbolOptions(geneSymbols);

                // Click to get to next step
                document.getElementById("load_plot_s").click();

            }
        }
    })
});

let plotSelect = null;
let geneSelect = null;

const createGeneSelectInstance = () => {
    geneSelect = NiceSelect.bind(document.getElementById("gene_select"), {
        placeholder: 'To search, start typing a gene name',
        searchable: true,
        allowClear: true,
        width: "resolve"
    });
}

const deleteDisplay = async(user_id, display_id) => {
    const payload = {user_id, display_id};
    return await axios.post("/cgi/delete_dataset_display.cgi", payload);
    // TODO: Update displays
}

const fetchAnalyses = async (datasetId) => {
    const { data } = await axios.get(`./api/h5ad/${datasetId}/analyses`);

    const { public: public_analyses, private: private_analyses } = data;
    return {public_analyses, private_analyses};

}

const fetchAvailablePlotTypes = async (user_id, session_id, dataset_id, analysis_id) => {
    const payload = {user_id, session_id, dataset_id, analysis_id};
    const res = await axios.post(`/api/h5ad/${dataset_id}/availableDisplayTypes`, payload);
    if (!res.success) {
        throw new Error(datasetId.message);
    }
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
    return data
}

const fetchGeneSymbols = async (datasetId, analysisId) => {
    const url = `./api/h5ad/${datasetId}/genes`;
    if (analysisId) url += `?analysis_id=${analysisId}`;

    const { data } = await axios.get(url);
    return [...new Set(data.gene_symbols)]; // Dataset may have a gene repeated in it, so resolve this.
}

const fetchH5adInfo = async (datasetId, analysisId) => {
    const url = `/api/h5ad/${datasetId}`
    if (analysisId) url += `?analysis_id=${analysisId}`;
    const {data} = await axios.get(url);
    const { obs_columns, obs_levels } = data;
    return {obs_columns, obs_levels}
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

const getSelectedValue = (select) => {
    // Get value from select2 element
    return select.selectedOptions[0].data.value;

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

const saveDatasetDisplay = async(displayId, dataset_id, user_id, label, plot_type, plotConfig) => {
    const payload = {
        id: displayId,
        dataset_id,
        user_id,
        label,
        plot_type,
        plotly_config: JSON.stringify({
        // depending on display type, this object will
        // have different properties
        ...this.config,
        }),
    };

    const res = await axios.post("/cgi/save_dataset_display.cgi", payload)
    if (res?.success) {
        if (displayId) {
            // TODO: update current display
        }
    }
}

const saveDefaultDisplay = async (display_id) => {
    const payload = {display_id};
    return await axios.post("/cgi/save_default_display.cgi", payload);
    // TODO: Update labels of displays
}

const updateGeneSymbolOptions= (geneSymbols) => {
    const geneOptionElts = geneSymbols.sort().map(cat => `<option>${cat}</option>`);
    if (geneSelect instanceof Object) {
        geneSelect.destroy();
    }
    const geneSelectTemplate = `
                <select id="gene_select">
                <option></option>
                ${geneOptionElts.join("")}
                <select>
                `;
    const geneSelectHtmlCollection = generateElements(geneSelectTemplate);
    document.getElementById("gene_collapsable").replaceChildren(geneSelectHtmlCollection);
    createGeneSelectInstance();
}

window.onload = () => {

    loadDatasetTree();

    // TODO: If url param "dataset_id"=<dataset_id>, then pre-populate

    // TODO: Figure out "scrollbar"
    // Initialize plot types
    plotSelect = NiceSelect.bind(document.getElementById("plot_type_select"), {
        placeholder: 'Choose how to plot',
        width: '25%',
        minimumResultsForSearch: -1
    })

};


