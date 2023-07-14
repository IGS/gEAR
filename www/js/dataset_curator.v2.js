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
let displayId = null;
let obsLevels = null;
let obsNotUsed = null;
let geneSymbol = null;

let datasetTree = null;

let plotSelect = null;
let geneSelect = null;

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
    return await axios.post("/cgi/get_dataset_display.cgi", payload);
}

const fetchDatsetDisplays = async (user_id, dataset_id) => {
    const payload = {user_id, dataset_id};
    // POST due to payload variables being sensitive
    return await axios.post("/cgi/get_dataset_displays.cgi", payload);
}

const fetchDatasetInfo = async (dataset_id) => {
    const payload = {dataset_id};
    const {title, is_public, owner_id} = await axios.post("/cgi/get_dataset_info.cgi", payload);
    return {title, is_public, owner_id};
}

const fetchDatasets = async (session_id) => {
    const payload = {session_id}
    const {data} = await axios.post("cgi/get_h5ad_dataset_list.cgi", convertToFormData(payload));

    // Populate select box with dataset information owned by the user
    const userDatasets = [];
    if (data.user.datasets.length > 0) {
        // User has some profiles
        data.user.datasets.forEach((item) => {
            if (item) {
                userDatasets.push({  title: item.title, type:"show", dataset_id: item.id, organism_id: item.organism_id });
            }
        });
    }
    // Next, add datasets shared with the user
    const sharedDatasets = [];
    if (data.shared_with_user.datasets.length > 0) {
        data.shared_with_user.datasets.forEach((item) => {
            if (item) {
                sharedDatasets.push({  title: item.title, type:"show", dataset_id: item.id, organism_id: item.organism_id });
            }
        });
    }
    // Now, add public datasets
    const domainDatasets = [];
    if (data.public.datasets.length > 0) {
        data.public.datasets.forEach((item) => {
            if (item) {
                domainDatasets.push({  title: item.title, type:"show", dataset_id: item.id, organism_id: item.organism_id });
            }
        });
    }

    return [
        {
            title: `Public datasets (${domainDatasets.length})`,
            type: "folder",
            children: domainDatasets
        },
        {
            title: `Shared Datasets (${sharedDatasets.length})`,
            type: "folder",
            children: sharedDatasets
        },
        {
            title: `Your Datasets (${userDatasets.length})`,
            type: "folder",
            children: userDatasets
        },
    ];
}

const fetchDefaultDisplay = async (user_id, dataset_id) => {
    const payload = {user_id, dataset_id};
    // POST due to payload variables being sensitive
    return await axios.post("/cgi/get_default_display.cgi", payload);
}

const fetchGeneSymbols = async (datasetId, analysisId) => {
    const url = `./api/h5ad/${datasetId}/genes`;
    if (analysisId) url += `?analysis_id=${analysisId}`;

    const { data } = await axios.get(url);
    return data.gene_symbols

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

const getSelectedPlotType = () => {
    // Get value from select2 element
    return plotSelect.selectedOptions[0].data.value;

}

const loadDatasetTree = async () => {

    // Change icons
    // https://mar10.github.io/wunderbaum/api/variables/common.iconMap.html


    datasetTree = new mar10.Wunderbaum({
        // https://mar10.github.io/wunderbaum/api/interfaces/wb_options.WunderbaumOptions.htm
        element: document.getElementById("dataset_tree"),
        source: await fetchDatasets(session_id),
        filter: {
            // https://github.com/mar10/wunderbaum/blob/c393a7a4e01a8c6bb9084002d3a70b9e18bcc974/src/wb_ext_filter.ts
            connectInput: "input#dataset_query",
            autoExpand: true,
            mode: "hide",
        },
        types: {
            folder: {icon: "mdi mdi-folder"},
            show: {icon: "mdi mdi-card-account-details-outline"}
        },
        autoCollapse: true,
        debugLevel: 2,
        //event handlers
        init: (e) => {
            // After tree is created and data is loaded
            e.tree.setFocus();
        },
        render: (e) => {
            // e.node was rendered. We may now modify the markup...
        },
    });
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

window.onload = () => {

    loadDatasetTree();

    // Initialize plot types
    plotSelect = NiceSelect.bind(document.getElementById("plot_type_select"), {
        placeholder: 'Choose how to plot',
        width: '25%',
        minimumResultsForSearch: -1
    })

    geneSelect = NiceSelect.bind(document.getElementById("gene_select"), {
        placeholder: 'To search, start typing a gene name',
        allowClear: true,
        width: 'resolve'
    });

};