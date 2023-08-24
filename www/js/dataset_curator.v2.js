// I use camelCase for my variable/function names to adhere to JS style standards
// Exception being functions that do fetch calls, so we can use JS destructuring on the payload

'use strict';

//TODO  - Create "error" reporting in each step

let plotConfig = {};  // Plot config that is passed to API or stored in DB

let datasetId = null;
let organismId = null;
let analysisObj = null;

let analysisSelect = null;
let plotSelect = null;
let geneSelect = null;

const userId = 622;     // ! It's me
const sessionId = "ee95e48d-c512-4083-bf05-ca9f65e2c12a"    // ! It's me
const session_id = sessionId;
const colorblindMode = false;

Cookies.set('gear_session_id', sessionId, { expires: 7 });

const plotlyPlots = ["bar", "line", "scatter", "tsne_dyna", "violin"];
const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];

const plotlyElt2Prop = {
    "x_axis_series":"x_axis"
    , "y_axis_series":"y_axis"
    , "label_series":"point_label"
    , "hide_x_ticks":"hide_x_labels"
    , "hide_y_ticks":"hide_y_labels"
    , "color_series":"color_name"
    , "facet_row_series":"facet_row"
    , "facet_col_series":"facet_col"
    , "x_axis_title":"x_title"
    , "y_axis_title":"y_title"
    , "xmin_value":"x_min"
    , "ymin_value":"y_min"
    , "xmax_value":"x_max"
    , "ymax_value":"y_max"
    , "hide_legend":"hide_legend"
    , "add_jitter":"jitter"
    , "marker_size_series":"size_by_group"
    , "marker_size":"markersize"
}
/*
    'vlines_body":"vlines"    // This is a special case
    color_map = req.get('colors')
    palette = req.get('color_palette')
    reverse_palette = req.get('reverse_palette')
    order = req.get('order', {})
*/

const scanpyElt2Prop = {
    "x_axis_series":"x_axis"
    , "y_axis_series":"y_axis"
    , "colorize_legend_by":"colorize_legend_by"
    , "plot_by_series":"plot_by_group"
    , "max_columns":"max_columns"
    , "skip_gene_plot":"skip_gene_plot"
    , "horizontal_legend":"horizontal_legend"
}
/*
        order = req.get('order', {})
*/

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
        document.getElementById('new_display').classList.add("is-loading");

        // Clear "current <whatever>" text
        document.getElementById("current_gene_c").style.display = "none";
        document.getElementById("current_analysis_c").style.display = "none";
        document.getElementById("current_plot_type_c").style.display = "none";


        datasetId = newDatasetId;

        // Clear (and update) options within nice-select2 structure.
        analysisObj = null;
        analysisSelect.clear(); // BUG: Figure out why this is triggering gene-population twice (check function in Github)
        geneSelect.clear();
        plotSelect.clear();

        // Fetch dataset information
        const {owner_id: ownerId} = await fetchDatasetInfo(datasetId);

        // displays
        const userDisplays = await fetchDatasetDisplays(userId, datasetId);
        const ownerDisplays = userId === ownerId ? [] : await fetchDatasetDisplays(ownerId, datasetId);
        const defaultDisplayId = await fetchDefaultDisplay(userId, datasetId);
        renderDisplayCards(userDisplays, ownerDisplays, defaultDisplayId);
        document.getElementById('new_display').classList.remove("is-loading");
        document.getElementById('new_display').disabled = false;

    })
});

const capitalized = (string) => string.charAt(0).toUpperCase() + string.slice(1)

const chooseAnalysis = async () => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const analysisText = analysisId || "Primary Analysis";

    // Display current selected analysis
    document.getElementById("current_analysis_c").style.display = "";
    document.getElementById("current_analysis").textContent = analysisText;

    // NOTE: For now, we can just pass analysis id only to tSNE and be fine
    // Any private dataset will belong to our user. Any public datasets can be found by the API "get_analysis" code.
    if (analysisId) {
        analysisObj = {id: analysisId};
    }

    if (analysisSelect.data.length < 5) return; // Have not retrieved analyses from API yet

    // Populate gene select element
    const geneSymbols = await fetchGeneSymbols(datasetId, analysisId);
    updateGeneSymbolOptions(geneSymbols);
}

/* Display has been chosen, so display analysis and plot type options */
const chooseDisplay = async () => {
    document.getElementById('new_display').classList.add("is-loading");
    document.getElementById("analysis_select").disabled = false;

    document.getElementById("plot_type_select").disabled = false;

    // analyses
    const {public_analyses: publicAnalyses, private_analyses: privateAnalyses} = await fetchAnalyses(datasetId);
    updateAnalysesOptions(privateAnalyses, publicAnalyses);
    analysisSelect.update();

    // plot types
    // NOTE: Believe this triggers the plotSelect "change" element
    const availablePlotTypes = await fetchAvailablePlotTypes(userId, sessionId, datasetId, undefined);
    for (const plotType in availablePlotTypes) {
        const isAllowed = availablePlotTypes[plotType];
        if (plotType === "tsne/umap_dynamic") {
            document.getElementById("tsne_dyna_opt").disabled = !isAllowed;
        } else {
            document.getElementById(`${plotType}_opt`).disabled = !isAllowed;
        }
    }
    plotSelect.update();
    document.getElementById('new_display').classList.remove("is-loading");


    // Populate gene select element
    const geneSymbols = await fetchGeneSymbols(datasetId, undefined);
    updateGeneSymbolOptions(geneSymbols);

    document.getElementById("plot_type_s").click();
}

const chooseGene = () => {
    if (!geneSelect.selectedOptions.length) return;   // Do not trigger after initial population

    // Display current selected gene
    document.getElementById("current_gene_c").style.display = "";
    document.getElementById("current_gene").textContent = getSelect2Value(geneSelect);

    document.getElementById("plot_options_s").click();
}

const choosePlotType = () => {
    plotConfig = {};    // Reset the plot config parameters
    if (!plotSelect.selectedOptions.length) return;   // Do not trigger after setting disable/enable on options

    // Do not display if default opt is chosen
    const plotType = getSelect2Value(plotSelect)
    if (plotType === "nope") return;

    // Display current selected plot type
    document.getElementById("current_plot_type_c").style.display = "";
    document.getElementById("current_plot_type").textContent = plotType;

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

const colorSVG = (chartData, plotConfig) => {
    // I found adding the mid color for the colorblind mode  skews the whole scheme towards the high color
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : plotConfig["low_color"];
    const midColor = colorblindMode ? null : plotConfig["mid_color"];
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : plotConfig["high_color"];

    // for those fields which have no reading, a specific value is sometimes put in instead
    // These are colored a neutral color
    const NA_FIELD_PLACEHOLDER = -0.012345679104328156;
    const NA_FIELD_COLOR = '#808080';

    //const scoreMethod = document.getElementById("scoring_method").value;
    const score = chartData.scores["gene"]
    const { min, max } = score;
    let color = null;
    // are we doing a three- or two-color gradient?
    if (midColor) {
        if (min >= 0) {
            // All values greater than 0, do right side of three-color
            color = d3
                .scaleLinear()
                .domain([min, max])
                .range([midColor, highColor]);
        } else if (max <= 0) {
            // All values under 0, do left side of three-color
            color = d3
                .scaleLinear()
                .domain([min, max])
                .range([lowColor, midColor]);
        } else {
            // We have a good value range, do the three-color
            color = d3
                .scaleLinear()
                .domain([min, 0, max])
                .range([lowColor, midColor, highColor]);
        }
    } else {
        color = d3
            .scaleLinear()
            .domain([min, max])
            .range([lowColor, highColor]);
    }


    // Load SVG file and set up the window
    const svg = document.getElementById("plot_container");
    const snap = Snap(svg);
    const svg_path = `datasets_uploaded/${datasetId}.svg`;
    Snap.load(svg_path, async (path) => {
        await snap.append(path)

        snap.select("svg").attr({
            width: "100%",
        });

        // Fill in tissue classes with the expression colors
        const {data: expression} = chartData;
        const tissues = Object.keys(chartData.data);   // dataframe
        const paths = Snap.selectAll("path, circle");
        paths.forEach((path) => {
            const tissue = path.node.className.baseVal;
            if (tissues.includes(tissue)) {
                if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                    path.attr('fill', NA_FIELD_COLOR);
                } else {
                    path.attr('fill', color(expression[tissue]));
                }
            }
        });

        // TODO: Potentially replicate some of the features in display.js like log-transforms and tooltips
    });

}

const createAnalysisSelectInstance = () => {
    analysisSelect = NiceSelect.bind(document.getElementById("analysis_select"), {
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

const createPlot = async () => {

    const plotType = getSelect2Value(plotSelect);
    const geneSymbol = getSelect2Value(geneSelect);
    plotConfig["gene_symbol"] = geneSymbol;

    const plotBtn = document.getElementById("plot_btn");
    const plotContainer = document.getElementById("plot_container");
    const paramsContainer = document.getElementById("params_container");

    plotContainer.replaceChildren();    // erase plot
    paramsContainer.replaceChildren();    // erase params

    // Set loading
    plotBtn.classList.add("is-loading");

    // Call API route by plot type
    if (plotlyPlots.includes(plotType)) {
        for (const elt in plotlyElt2Prop) {
            // Some elements are only present in certain plot type configurations
            if (document.getElementById(elt)) {
                plotConfig[plotlyElt2Prop[elt]] = document.getElementById(elt).value || undefined;
                if (document.getElementById(elt)?.type == "checkbox") {
                    plotConfig[plotlyElt2Prop[elt]] = document.getElementById(elt).checked;
                }
            }
        }
        // Get order
        plotConfig["order"] = {};
        for (const param of ["x_axis", "y_axis", "faceet_row", "facet_col"]) {
            if (document.getElementById(`${param}_order_list`)) {
                const series = document.getElementById(`${param}_series`).value;
                const serialized = sortable(`#${param}_order_list`, 'serialize')[0].items;
                // Sort by "sortable" index position
                plotConfig["order"][series] = serialized.map((val) => val.label);
            }
        }

        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("color_series").value;
        if (colorSeries && colorElts.length) {
            // Input is either color mapping or just the series
            plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("_color", "");
                plotConfig["colors"][group] = field.value;
            })
        }

        // Get vlines
        const vlineFields = document.getElementsByClassName("js-vline-field");
        plotConfig["vlines"] = [...vlineFields].map((field) => {
            const vlinePos = field.querySelector(":scope .js-vline-pos").value;
            const vlineStyle = field.querySelector(":scope .js-vline-style-select select").value;
            // Return either objects or nothing (which will be filtered out)
            return vlinePos ?  {"vl_pos":vlinePos, "vl_style":vlineStyle} : null;
        }).filter(x => x !== null);

        // Get data and set up the image area
        const data = await fetchPlotlyData(plotConfig, datasetId, plotType, analysisObj, userId, colorblindMode);
        const {plot_json: plotJson} = data;

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const divElt = generateElements('<div class="container is-max-desktop" id="plotly_preview"></div>');
        plotContainer.append(divElt);
        Plotly.purge("plotly_preview"); // clear old Plotly plots

        // TODO: Throw error if "plotJson" is null
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const curatorDisplayConf = postPlotlyConfig.curator;
        const custonConfig = getPlotlyDisplayUpdates(curatorDisplayConf, plotType, "config");
        Plotly.newPlot("plotly_preview", plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(curatorDisplayConf, plotType, "layout")
        Plotly.relayout("plotly_preview", custonLayout)

    } else if (scanpyPlots.includes(plotType)) {
        for (const elt in scanpyElt2Prop) {
            plotConfig[scanpyElt2Prop[elt]] = document.getElementById(elt).value || undefined;
            if (document.getElementById(elt)?.type == "checkbox") {
                plotConfig[scanpyElt2Prop[elt]] = document.getElementById(elt).checked;
            }
        }
        // Get order
        plotConfig["order"] = {};
        for (const param of ["plot_by_series"]) {
            if (document.getElementById(`${param}_order_list`)) {
                const series = document.getElementById(`${param}_series`).value;
                const serialized = sortable(`#${param}_order_list`, 'serialize')[0].items;
                // Sort by "sortable" index position
                plotConfig["order"][series] = serialized.map((val) => val.label);
            }
        }
        // Get colors
        const colorElts = document.getElementsByClassName("js-plot-color");
        const colorSeries = document.getElementById("colorize_legend_by").textContent;
        if (colorSeries && colorElts.length) {
            plotConfig["colors"] = {};
            [...colorElts].map((field) => {
                const group = field.id.replace("_color", "");
                plotConfig["colors"][group] = field.value;
            })
        }


        // If user did not want to have a colorized annotation, ensure it does not get passed to the scanpy code
        if (!(document.getElementById("show_colorized_legend").checked)) {
            plotConfig["colorize_legend_by"] = null;
            plotConfig["plot_by_group"] = null;
            plotConfig["max_columns"] = null;
            plotConfig["skip_gene_plot"] = false;
            plotConfig["horizontal_legend"] = false;
        }

        const data = await fetchTsneImageData(plotConfig, datasetId, plotType, analysisObj, userId, colorblindMode);
        const {image} = data;
        const imgElt = generateElements('<img class="image" id="tsne_preview"></img>');
        plotContainer.append(imgElt);

        if (image) {
            document.getElementById("tsne_preview").setAttribute("src", `data:image/png;base64,${image}`);
        }
    } else if (plotType === "svg") {
        const data = await fetchSvgData(geneSymbol, datasetId)
        plotConfig["low_color"] = document.getElementById("low_color").value;
        plotConfig["mid_color"] = document.getElementById("mid_color").value;
        plotConfig["high_color"] = document.getElementById("high_color").value;

        // If user did not choose a mid-color, set it as null instead of to black
        if (!(document.getElementById("enable_mid_color").checked)) {
            plotConfig["mid_color"] = null;
        }
        colorSVG(data, plotConfig)
    } else {
        console.warn(`Plot type ${plotType} selected for plotting is not a valid type.`)
        return;
    }

    // Stop loader
    plotBtn.classList.remove("is-loading");

    // Populate params in panel
    for (const param in plotConfig) {
        if (!plotConfig[param]) {
            delete plotConfig[param];
            continue
        };
        if (["object", "array"].includes(typeof(plotConfig[param]))) {
            if (typeof(plotConfig[param]) === "object" && !Object.keys(plotConfig[param]).length) {
                delete plotConfig[param];
                continue
            }
            if (typeof(plotConfig[param]) === "array" && !(plotConfig[param].length)) {
                delete plotConfig[param];
                continue
            };
            const prettyPrint = JSON.stringify(plotConfig[param], null, "\t");
            const paramElt = generateElements(`<p class="columns">
                <span class="column has-text-weight-medium">${capitalized(param)}</span>
                <span class="column has-text-weight-bold">${prettyPrint}</span>
                </p>`);
            paramsContainer.append(paramElt);
        } else {
            const prettyPrint = plotConfig[param];
            const paramElt = generateElements(`<p class=" is-flex is-justify-content-space-between">
                <span class=" has-text-weight-medium">${capitalized(param)}</span>
                <span class=" has-text-weight-bold">${prettyPrint}</span>
                </p>`);
            paramsContainer.append(paramElt);
        }
    }

    // Hide this view
    document.getElementById("content_c").style.display = "none";
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").style.display = "";

    // TODO: Add alert for non-success w/ message

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

const fetchPlotlyData = async (plotConfig, datasetId, plot_type, analysis, analysis_owner_id, colorblind_mode)  => {
    // NOTE: gene_symbol already passed to plotConfig
    const payload = { ...plotConfig, plot_type, analysis, analysis_owner_id, colorblind_mode };
    const { data } = await axios.post(`/api/plot/${datasetId}`, payload);
    return data
}

const fetchSvgData = async (geneSymbol, datasetId) => {
    const { data } = await axios.get(`/api/plot/${datasetId}/svg?gene=${geneSymbol}`);
    return data
};

const fetchTsneImageData = async (plotConfig, datasetId, plot_type, analysis, analysis_owner_id, colorblind_mode) => {
    // NOTE: gene_symbol already passed to plotConfig
    const payload = { ...plotConfig, plot_type, analysis, analysis_owner_id, colorblind_mode };
    const { data } = await axios.post(`/api/plot/${datasetId}/tsne`, payload);
    return data
}

const getPlotlyDisplayUpdates = (plotConfObj, plotType, category) => {
    // Get updates and additions to plot from the plot_display_config JS object
    let updates = {};
    for (const idx in plotConfObj) {
        const conf = plotConfObj[idx];
        // Get config (data and/or layout info) for the plot type chosen, if it exists
        if (conf.plot_type == "all" || conf.plot_type == plotType) {
            const update = category in conf ? conf[category] : {};
            updates = {...updates, ...update};    // Merge updates
        }
    }
    return updates;
  }

const getSelect2Value = (select) => {
    // Get value from select2 element
    return select.selectedOptions[0].data.value;
}

const includePlotParamOptions = async () => {
    const plotType = getSelect2Value(plotSelect);

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

        setupSVGOptions();
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

const renderColorPicker = (allCats, seriesName) => {
    const colorsContainer = document.getElementById("colors_container");

    const seriesNameHtml = generateElements(`<p class="has-text-weight-bold is-underlined">${seriesName}</p>`);
    colorsContainer.append(seriesNameHtml);

    // If there are colors from the observations use them
    if (Object.keys(allCats).includes(`${seriesName}_colors`)) {
        const seriesColorsHtml = generateElements(`<p>Using the color mapping found for this data series within the dataset.</p>`);
        colorsContainer.append(seriesColorsHtml);
        return;
    };


    // Otherwise d3 category10 colors
    const swatchColors = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"];

    let counter = 0;
    for (const group of allCats[seriesName]) {
        const darkerLevel = Math.floor(counter / 10);
        const baseColor = swatchColors[counter%10];
        const groupColor = darkerLevel > 0
            ? d3.color(baseColor).darker(darkerLevel).formatHex()
            : baseColor;    // Cycle through swatch but make darker if exceeding 10 groups
        counter++;
        const groupHtml = generateElements(`
            <p class="is-flex is-justify-content-space-between">
                <span class="has-text-weight-medium">${group}</span>
                <input class="js-plot-color" id="${group}_color" type="color" value="${groupColor}" aria-label="Select a color" />
            </p>
        `);
        colorsContainer.append(groupHtml)
    }
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
                            <footer class="card-footer buttons">
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
                            <footer class="card-footer buttons">
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

const renderOrderSortable = (allCats, param, series) => {
    const orderContainer = document.getElementById("order_container");

    // If param already has HTML defined, remove it
    const paramOrder = document.getElementById(`${param}_order`);
    if (paramOrder) {
        paramOrder.remove();
    }

    // If continouous series, cannot sort.
    if (!(Object.keys(allCats).includes(series))) return;

    // If series is used in another param, also return
    if (document.getElementById(`${series}_order`)) return;

    // Create parent template
    const parentList = `<ul id="${param}_order_list" class="content"></ul>`;
    const template = `
        <div id="${param}_order" class="column is-one-quarter">
            <p class="has-text-weight-bold is-underlined">${param} (<span id="${series}_order">${series}</span>)</p>
            ${parentList}
        </div>
    `;

    const htmlCollection = generateElements(template);
    orderContainer.append(htmlCollection);

    // Add in list elements
    for (const group of allCats[series]) {
        const listElt = `<li class="has-background-grey-lighter has-text-dark">${group}</li>`;
        const listCollection = generateElements(listElt);
        document.getElementById(`${param}_order_list`).append(listCollection);
    }

    // Create sortable for this series
    sortable(`#${param}_order_list`, {
        hoverClass: "has-text-weight-bold"
        , itemSerializer(item, container) {
            item.label = item.node.textContent
            return item
        },
    });

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
            ...plotConfig,  // depending on display type, this object will have different properties
        }),
    };

    if (!displayId) delete payload.id;  // Prevent passing in "null" as a string.

    const data = await axios.post("/cgi/save_dataset_display.cgi", convertToFormData(payload));
    const {display_id} = data;
    return display_id;
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
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
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

        // If x-axis is continuous show vline stuff, otherwise hide
        document.getElementById("x_axis_series").addEventListener("change", (event) => {
            const vLinesField = document.getElementById("vlines_container")
            if (!(catColumns.includes(event.target.value))) {
                vLinesField.style.display = "none";
                // Clear all existing vlines
                const vLineTemplate = document.querySelector(".js-vline-field");
                document.querySelectorAll(".js-vline-field").remove()
                // Add the first vline back
                const vLinesBody = document.getElementById("vlines_body");
                vLinesBody.append(vLineTemplate.cloneNode(true));
                document.querySelector(".js-vline-pos").value = "";
                return;
            }
            vLinesField.style.display = "";

        });

        // Vertical line add and remove events
        const vLinesBody = document.getElementById("vlines_body");
        const vLineTemplate = document.querySelector(".js-vline-field");
        document.getElementById("vline_add_btn").addEventListener("click", (event) => {
            vLinesBody.append(vLineTemplate.cloneNode(true));
            // NOTE: Currently if original is set before cloning, values are copied to clone
            document.getElementById("vline_remove_btn").disabled = false;
        })
        document.getElementById("vline_remove_btn").addEventListener("click", (event) => {
            const lastVLine = document.querySelector(".js-vline-field:nth-last-of-type(1)");   // last sibling div
            lastVLine.remove();
            if (vLineTemplate.length < 2) document.getElementById("vline_remove_btn").disabled = true;
        })
    }

    const validationElts = document.getElementsByClassName("js-plot-req");
    for (const elt of validationElts ) {
        elt.addEventListener("change", () => {
            const xVal = document.getElementById("x_axis_series").value;
            const yVal = document.getElementById("y_axis_series").value;
            document.getElementById("plot_btn").disabled = (xVal && yVal) ? false : true;
        })
    }

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    trigger(document.getElementById("x_axis_series"), "change");

    // If color series is selected, let user choose colors.
    const colorSeries = document.getElementById("color_series");
    colorSeries.addEventListener("change", (event) => {
        document.getElementById("colors_section").style.display = "none";
        document.getElementById("colors_container").replaceChildren();
        if (!event.target.value) {
            return;
        }
        if (!(catColumns.includes(event.target.value))) {
            // ? Continuous series colorbar picker
            return;
        }
        document.getElementById("colors_section").style.display = "";
        renderColorPicker(levels, event.target.value);
        return;
    })

    // Certain elements trigger plot order
    const plotOrderElts = document.getElementsByClassName("js-plot-order");
    for (const elt of plotOrderElts) {
        elt.addEventListener("change", (event) => {
            const paramId = event.target.id;
            const param = paramId.replace("_series", "");
            // NOTE: continuous series will be handled in the function
            renderOrderSortable(levels, param, event.target.value);
            document.getElementById("order_section").style.display = "";
        });
    }

    // TODO: Set up validation checkers
}

/* Set up the scanpy-based plot options, such as "select" elements, events, etc. */
const setupScanpyOptions = async () => {
    const analysisValue = analysisSelect.selectedOptions.length ? getSelect2Value(analysisSelect) : undefined;
    const analysisId = (analysisValue && analysisValue > 0) ? analysisValue : null;
    const plotType = getSelect2Value(plotSelect);
    const {obs_columns: allColumns, obs_levels: levels} = await fetchH5adInfo(datasetId, analysisId);
    const catColumns = Object.keys(levels);

    let xDefaultOption = null;
    let yDefaultOption = null;

    // If these exist, make the default option
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

    updateSeriesOptions("x_axis_series", allColumns, true, xDefaultOption);
    updateSeriesOptions("y_axis_series", allColumns, true, yDefaultOption);
    updateSeriesOptions("colorize_legend_by", allColumns, false);
    updateSeriesOptions("plot_by_series", catColumns, false);

    const validationElts = document.getElementsByClassName("js-plot-req");
    for (const elt of validationElts ) {
        elt.addEventListener("change", () => {
            const xVal = document.getElementById("x_axis_series").value;
            const yVal = document.getElementById("y_axis_series").value;
            document.getElementById("plot_btn").disabled = (xVal && yVal) ? false : true;
        })
    }


    const colorizeLegendBy = document.getElementById("colorize_legend_by");
    const plotBySeries = document.getElementById("plot_by_series");
    const maxColumns = document.getElementById('max_columns');
    const skipGenePlot = document.getElementById("skip_gene_plot");
    const horizontalLegend = document.getElementById("horizontal_legend");

    // Disable colorize_legend options if we are not using it
    document.getElementById("show_colorized_legend").addEventListener("change", (event) => {
        colorizeLegendBy.disabled = true;
        plotBySeries.disabled = true;
        maxColumns.disabled = true;
        skipGenePlot.disabled = true;
        horizontalLegend.disabled = true;
        if (!event.target.checked) {
            return;
        }
        colorizeLegendBy.disabled = false;
        plotBySeries.disabled = false;
        maxColumns.disabled = false;
        skipGenePlot.disabled = false;
        horizontalLegend.disabled = false;
    });

    // Do certain things if the chosen annotation series is categorical or continuous
    colorizeLegendBy.addEventListener("change", (event) => {
        plotBySeries.disabled = false;
        maxColumns.disabled = false;
        horizontalLegend.disabled = false;

        // If colorized legend is continuous, we cannot plot by group
        if ((catColumns.includes(event.target.value))) {
            return;
        }
        plotBySeries.disabled = true;
        maxColumns.disabled = true;
        horizontalLegend.disabled = true;
    });

    // Plotting by group plots gene expression, so cannot skip gene plots.
    plotBySeries.addEventListener("change", (event) => {
        skipGenePlot.disabled = false;
        maxColumns.disabled = true;
        document.getElementById("order_section").style.display = "none";
        document.getElementById("order_container").replaceChildren();
        if (!event.target.value) {
            return;
        }
        skipGenePlot.checked = false;
        skipGenePlot.disabled = true;
        maxColumns.disabled = false;
        document.getElementById("order_section").style.display = "";
        renderOrderSortable(levels, "plot_by_group", event.target.value);

    });

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    trigger(document.getElementById("x_axis_series"), "change");

}

const setupSVGOptions = () => {

    const validationElts = document.getElementsByClassName("js-plot-req");
    for (const elt of validationElts ) {
        elt.addEventListener("change", () => {
            const highVal = document.getElementById("high_color").value;
            const lowVal = document.getElementById("low_color").value;
            document.getElementById("plot_btn").disabled = (highVal && lowVal) ? false : true;
        })
    }

    const enableMidColor = document.getElementById("enable_mid_color");
    const midColorField = document.getElementById("mid_color_field");
    enableMidColor.addEventListener("change", (event) => {
        midColorField.style.display = "none";
        if (event.target.checked) {
            midColorField.style.display = "";
        }
    });

    // Trigger event to enable plot button (in case we switched between plot types, since the HTML vals are saved)
    trigger(document.getElementById("high_color"), "change");

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
    // NOTE: For now, we can just pass analysis id only to tSNE and be fine
    for (const analysis of privateAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.value = analysis.id;
        //option.dataset.type = analysis.type;
        //option.dataset.owner_id = analysis.user_id;
        privateAnalysesElt.append(option);
    }
    for (const analysis of publicAnalyses) {
        const option = document.createElement("option");
        option.textContent = analysis.label;
        option.value = analysis.id;
        //option.dataset.type = analysis.type;
        //option.dataset.owner_id = analysis.user_id;
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
    firstOption.textContent = "Please select a gene";
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
        elt.append(expression);
        expression.textContent = "expression";
        expression.value = "raw_value";
        if ("raw_value" === defaultOption) {
            expression.selected = true;
        }
    }

    // Add categories
    for (const group of seriesArray.sort()) {
        // Skip columns listed as "_colors" as they just provide colors for another series
        if (group.includes("_colors")) continue;

        const option = document.createElement("option");
        option.textContent = group;
        // Change X_pca/X_tsne/X_umap text content to be more user_friendly
        if (group.includes("X_") && (
            group.includes("pca")
            || group.includes("tsne")
            || group.includes("umap")
        )) {
            option.textContent = `${group} (from selected analysis)`;
        }
        option.value = group;
        elt.append(option);
        // NOTE: It is possible for a default option to not be in the list of groups.
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
document.getElementById("analysis_select").addEventListener("change", chooseAnalysis);
document.getElementById("plot_type_select").addEventListener("change", choosePlotType);
document.getElementById("gene_select").addEventListener("change", chooseGene);
document.getElementById("plot_btn").addEventListener("click", createPlot);
document.getElementById("save_json_config").addEventListener("click", () => {
    // Config plot configuration to JSON for sharing (or passing to API by hand)
    const plotType = getSelect2Value(plotSelect);
    const blob = new Blob([JSON.stringify({...plotConfig, plot_type:plotType})]);
    const link = document.createElement("a");
    link.download = "gEAR_plot_configuration.json";
    link.href = window.URL.createObjectURL(blob);
    link.click()
    // Give confirmation
    document.getElementById("save_json_config").classList.add("is-success");
    setTimeout(() => {
        document.getElementById("save_json_config").classList.remove("is-success");
    }, 1000);
});
document.getElementById("save_display_btn").addEventListener("click", async () => {
    // Save new plot display.
    const label = document.getElementById("new_display_label").value;
    const plotType = getSelect2Value(plotSelect);
    const displayId = await saveDatasetDisplay(null, datasetId, userId, label, plotType, plotConfig);
    saveDefaultDisplay(displayId);
    // Give confirmation
    document.getElementById("save_display_btn").classList.add("is-success");
    setTimeout(() => {
        document.getElementById("save_display_btn").classList.remove("is-success");
    }, 1000);
    // TODO: Reload new displays if user goes back to parameter editing
});
document.getElementById("edit_params").addEventListener("click", () => {
    // Hide this view
    document.getElementById("content_c").style.display = "";
    // Generate and display "post-plotting" view/container
    document.getElementById("post_plot_content_c").style.display = "none";
})
