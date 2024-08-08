"use strict";

/*
    Classes representing overall analysis (pipeline) elements and their child classes.
*/

// requires common.js

let analysisLabels = new Set();

// TODO; Standardize runAnalysis() steps and similar functions in terms of UI manipulation (i.e. disable/hide buttons)

class Analysis {
    constructor ({
        id = uuid(),
        datasetObj = null,
        datasetIsRaw = true,
        label = `Unlabeled ${commonDateTime()}`,
        type,
        vetting,
        userSessionId = CURRENT_USER.session_id,
        genesOfInterest = [],
        groupLabels = []
    } = {}) {
        this.id = id;
        this.userSessionId = userSessionId;
        this.dataset = datasetObj;    // The dataset object
        this.type = type;
        this.vetting = vetting;
        this.label = label;
        this.datasetIsRaw = datasetIsRaw;

        // TODO: All this is overwritten if getStoredAnalysis is called, so probably need to adjust the constructor
        // All new analysis start as "primary" analyses, but true primary analyses share the same ID as the dataset
        if (this.id === this.dataset?.id) {
            this.labeledTsne = new AnalysisStepLabeledTsne(this);   // Only for "primary" analyses
            // marker genes next
            // compare genes next
        } else {
            this.primaryFilter = new AnalysisStepPrimaryFilter(this);
            this.qcByMito = new AnalysisStepQCByMito(this);
            this.selectVariableGenes = new AnalysisStepSelectVariableGenes(this);
            this.pca = new AnalysisStepPCA(this);
            this.tsne = new AnalysisSteptSNE(this);
            this.clustering = new AnalysisStepClustering(this); // The old "louvain" step, which is now done with "leiden"
            // marker genes next
            this.clusteringEdit = new AnalysisStepClustering(this, "edit");
            // compare genes next
        }
        this.markerGenes = new AnalysisStepMarkerGenes(this);
        // step-based properties
        this.markerGenes.genesOfInterest = Array.isArray(genesOfInterest) ? new Set(genesOfInterest) : new Set();

        this.compareGenes = new AnalysisStepCompareGenes(this);

        this.groupLabels = groupLabels;

    }

    /**
     * Checks dependencies and runs the specified callback function.
     *
     * @param {Function} callback - The callback function to be executed.
     * @param {Object} opts - The options object to be passed to the callback function.
     */
    async checkDependenciesAndRun(callback, opts) {
        if (this.type === 'public') {
            this.copyToUserUnsaved(callback, opts);
        } else {
            await callback(opts);
        }
    }

    /**
     * Converts a camelCase object to a snake_case object.
     * @param {Object} camelCaseObj - The camelCase object to be converted.
     * @returns {Object} - The snake_case object.
     */
    static convertToJson(camelCaseObj) {
        // If the object is not an object or is null, return it as is
        if (camelCaseObj === null || typeof camelCaseObj !== 'object') {
            return camelCaseObj;
        }

        const toSnakeCase = (str) => {
            return str.replace(/[A-Z]/g, letter => `_${letter.toLowerCase()}`);
        }

        return Object.keys(camelCaseObj).reduce((result, key) => {
            result[toSnakeCase(key)] = camelCaseObj[key];
            return result;
        }, {});
    }

    /**
     * Copies the dataset analysis to a new analysis of the specified type.
     * @param {string} destType - The type of the destination analysis.
     * @returns {Promise<any>} - A promise that resolves with the copied analysis data.
     */
    async copyDatasetAnalysis(destType) {
        const params = {
            session_id: this.userSessionId,
            dataset_id: this.dataset.id,
            source_analysis_id: this.id,
            dest_analysis_id: this.id,
            source_analysis_type: this.type,
            dest_analysis_type: destType
        }

        const {data} = await axios.post("./cgi/copy_dataset_analysis.cgi", convertToFormData(params));
        return data;
    }

    /**
     * Makes a copy of the analysis to the current user's unsaved area.
     * The most common use case is to copy the public analysis of another user so changes can be made.
     *
     * @param {Function} callback - The callback function to be executed after the copy is made.
     * @param {Object} opts - Additional options for the copy operation.
     */
    async copyToUserUnsaved(callback, opts) {
        const newAnalysisId = uuid();

        try {
            const data = await this.copyDatasetAnalysis('user_unsaved');

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.type = 'user_unsaved';
            this.id = newAnalysisId;
            this.userSessionId = CURRENT_USER.session_id;

            document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
            document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");

            await this.getSavedAnalysesList(this.dataset.id, newAnalysisId);

            if (callback) {
                await callback(opts);
            }

        } catch (error) {
            createToast(`Error copying analysis: ${error.message}`);
        }
    }

    /**
     * Deletes the analysis associated with the current instance.
     * @throws {Error} If the deletion is unsuccessful or an error occurs.
     */
    async delete() {
        const {data} = await axios.post("./cgi/delete_dataset_analysis.cgi", convertToFormData({
            session_id: this.userSessionId,
            dataset_id: this.dataset.id,
            analysis_id: this.id,
            analysis_type: this.type
        }));

        if (!data.success || data.success < 1) {
            const error = data.error || "Unknown error. Please contact gEAR support.";
            throw new Error(error);
        }

        // Trigger the selection of a 'New' analysis
        document.querySelector(UI.newAnalysisOptionElt).setAttribute("selected", "selected");
        document.querySelector(UI.analysisSelect).dispatchEvent(new Event("change"));
        await this.getSavedAnalysesList(this.dataset.id, 0);
        resetStepperWithHrefs(UI.primaryFilterSection);
    }

    async download() {

        // create URL parameters
        const params = {
            session_id: this.userSessionId,
            dataset_id: this.dataset.id,
            analysis_id: this.id,
            type: "h5ad",
        }

        // download the h5ad
        const url = `./cgi/download_source_file.cgi?${new URLSearchParams(params).toString()}`;
        try {
            const {data} = await axios.get(url, {responseType: 'blob'});
            const blob = new Blob([data], {type: 'application/octet-stream'});
            const downloadUrl = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = downloadUrl;
            a.download = `${this.dataset.id}.${this.id}.h5ad`;
            a.click();
        } catch (error) {
            console.error("Error downloading analysis h5ad", this);
            logErrorInConsole(error);
            createToast(`Error downloading analysis h5ad`);
        }

    }

    /**
     * Retrieves the stored analysis data from the server.
     * @returns {Promise<void>} A promise that resolves when the analysis data is retrieved.
     */
    async getStoredAnalysis() {

        // Some dataset info (like organism ID) may be lost when loading an analysis from JSON
        const datasetObj = this.dataset;

        try {
            const {data} = await axios.post("./cgi/get_stored_analysis.cgi", convertToFormData({
                analysis_id: this.id,
                analysis_type: this.type,
                session_id: this.userSessionId,
                dataset_id: this.dataset.id
            }));

            // Load the analysis data and assign it to the current instance
            const ana = Analysis.loadFromJson(data);
            Object.assign(this, ana);

            // If tSNE was calculate, show the labeled tSNE section
            // Mainly for primary analyses
            // ? verify claim
            document.querySelector(UI.labeledTsneSection).classList.add("is-hidden");
            if (this.type === "primary" && data['tsne']['tsne_calculated']) {
                document.querySelector(UI.labeledTsneSection).classList.remove("is-hidden");

                // Initialize the labeled tSNE step
                this.labeledTsne = new AnalysisStepLabeledTsne(this);
            }

        } catch (error) {
            logErrorInConsole(`Failed ID was: ${datasetId} because msg: ${error}`);
            createToast(`Error getting stored analysis`);
        }

        // Restore the dataset object
        this.dataset = datasetObj;

    }

    /**
     * Retrieves the list of saved analyses for a given dataset and populates the analysis options in the UI.
     *
     * @param {string} datasetId - The ID of the dataset.
     * @param {string} selectedAnalysisId - The ID of the selected analysis.
     * @param {string} forPage - The page for which the analysis list is being retrieved.
     * @returns {Promise<void>} - A promise that resolves when the analysis list is retrieved and populated in the UI.
     * @throws {Error} - If there is an error retrieving the saved analyses.
     */
    async getSavedAnalysesList(datasetId, selectedAnalysisId, forPage) {
        try {
            const {data} = await axios.post("./cgi/get_stored_analysis_list.cgi", convertToFormData({
                dataset_id: datasetId,
                session_id: this.userSessionId
            }));

            // Create an empty option for the analysis select element
            const emptyAnalysisOption = document.createElement("option");
            emptyAnalysisOption.value = "";
            emptyAnalysisOption.disabled = true;
            emptyAnalysisOption.textContent = "None found";

            const thisAnalysisLabels = new Set();

            /*
                ? 'primary' analysis is too vague.  That could include clustering and/or
                dimensionality reduction.  The workbench needs clustering to be present,
                while the curator can do fine with just UMAP.  Whether 'primary analysis'
                is added to the menu needs to take this into effect.
            */

            const appendAnalysisOption = (parentSelector, analysis) => {
                const option = document.createElement("option");
                option.dataset.analysisId = analysis.id;
                option.dataset.analysisType = analysis.type;
                option.dataset.datasetId = analysis.dataset_id;
                option.textContent = analysis.label || "Unlabeled"
                // ? Using standard HTML, cannot add icons to options, so making icons by vetting status is not possible

                // Add to analysis optgroup
                document.querySelector(parentSelector).appendChild(option);
            }

            // Clear the analysis optgroups
            document.querySelector(UI.analysisPrimaryElt).replaceChildren();
            document.querySelector(UI.analysisUnsavedElt).replaceChildren();
            document.querySelector(UI.analysisSavedElt).replaceChildren();
            document.querySelector(UI.analysisPublicElt).replaceChildren();

            // primary
            if (data.primary.length) {
                if (forPage === 'sc_workbench') {
                    if (data.primary[0].louvain?.calculated || data.primary[0].clustering?.calculated) {
                        for (const analysis of data.primary) {
                            thisAnalysisLabels.add(analysis.label);
                            appendAnalysisOption(UI.analysisPrimaryElt, analysis);

                        }
                    } else {
                        document.querySelector(UI.analysisPrimaryElt).appendChild(emptyAnalysisOption.cloneNode(true));
                    }
                } else {

                    for (const analysis of data.primary) {
                        thisAnalysisLabels.add(analysis.label);
                        appendAnalysisOption(UI.analysisPrimaryElt, analysis);
                    }
                }
            } else {
                document.querySelector(UI.analysisPrimaryElt).appendChild(emptyAnalysisOption.cloneNode(true));
            }

            // unsaved
            if (data.user_unsaved.length) {
                // sort by label
                data.user_unsaved.sort((a, b) => a.label.localeCompare(b.label));

                for (const analysis of data.user_unsaved) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(UI.analysisUnsavedElt, analysis);
                }
            } else {
                document.querySelector(UI.analysisUnsavedElt).appendChild(emptyAnalysisOption.cloneNode(true));
            }

            // saved
            if (data.user_saved.length) {
                // sort by label
                data.user_saved.sort((a, b) => a.label.localeCompare(b.label));

                for (const analysis of data.user_saved) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(UI.analysisSavedElt, analysis);
                }
            } else {
                document.querySelector(UI.analysisSavedElt).appendChild(emptyAnalysisOption.cloneNode(true));
            }

            // public
            if (data.public.length) {
                // sort by label
                data.public.sort((a, b) => a.label.localeCompare(b.label));

                for (const analysis of data.public) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(UI.analysisPublicElt, analysis);
                }
            } else {
                document.querySelector(UI.analysisPublicElt).appendChild(emptyAnalysisOption.cloneNode(true));
            }

            // preselect any analysis ID
            document.querySelector(`#analysis-select option[data-analysis-id="${selectedAnalysisId}"]`).setAttribute("selected", "selected");

            analysisLabels = thisAnalysisLabels;


        } catch (error) {
            createToast(`Error getting saved analyses: ${error}`);
            logErrorInConsole(`Failed ID was: ${datasetId} because msg: ${error}`);
        }

    }

    /**
     * Loads an Analysis object from JSON data.
     *
     * @param {Object} data - The JSON data representing the Analysis object.
     * @returns {Analysis} The loaded Analysis object.
     */
    static async loadFromJson(data) {
        const analysis = new Analysis({
            id: data.id,
            datasetObj: data.dataset,
            datasetIsRaw: data.dataset_is_raw,
            label: data.label,
            type: data.type,
            userSessionId: data.user_session_id || CURRENT_USER.session_id,
            groupLabels: data.group_labels,
            genesOfInterest: data.genesOfInterest
        });

        document.querySelector(UI.btnProgressGuideElt).classList.remove("is-hidden");

        if (analysis.id === analysis.dataset?.id) {
            // If showing a primary display we only want to show marker genes and gene comparison
            //  tools
            document.querySelector(UI.primaryFilterSection).classList.add("is-hidden");
            document.querySelector(UI.qcByMitoSection).classList.add("is-hidden");
            document.querySelector(UI.selectVariableGenesSection).classList.add("is-hidden");
            document.querySelector(UI.pcaSection).classList.add("is-hidden");
            document.querySelector(UI.tsneSection).classList.add("is-hidden");
            document.querySelector(UI.clusteringSection).classList.add("is-hidden");
            document.querySelector(UI.clusteringEditSection).classList.add("is-hidden");

            // labeled tSNE does not have a data object section. Only needs dataset ID.

            // Show the primary analysis stepper
            document.querySelector(UI.primaryStepsElt).classList.remove("is-hidden");
            resetStepperWithHrefs(UI.markerGenesSection);
            // This is initially unclickable due to this step being initially unclickable for de novo analyses
            document.querySelector(UI.markerGenesSection).classList.remove("is-pointer-events-none");

            return analysis
        }

        // Show the de novo analysis stepper
        // Since the stepper relies on finding the "not hidden" stepper, we need to do this first.
        document.querySelector(UI.deNovoStepsElt).classList.remove("is-hidden");
        resetStepperWithHrefs(UI.primaryFilterSection);

        analysis.primaryFilter = AnalysisStepPrimaryFilter.loadFromJson(data.primary_filter, analysis);

        analysis.qcByMito = AnalysisStepQCByMito.loadFromJson(data.qc_by_mito, analysis);

        analysis.selectVariableGenes = AnalysisStepSelectVariableGenes.loadFromJson(data.select_variable_genes, analysis);

        analysis.pca = AnalysisStepPCA.loadFromJson(data.pca, analysis);

        analysis.tsne = AnalysisSteptSNE.loadFromJson(data.tsne, analysis);

        // Generalize the clustering step instead of being specific to louvain
        // Not doing anything with data.clustering yet but would like to
        if (data.louvain) {
            data.clustering = data.louvain;
            if (!data.clustering?.mode) {
                data.clustering.mode = "initial";
            }
        }

        analysis.clustering = AnalysisStepClustering.loadFromJson(data.clustering, analysis);

        analysis.markerGenes = await AnalysisStepMarkerGenes.loadFromJson(data.marker_genes, analysis);

        // Support legacy data.
        const clusteringEditData = data.clustering_edit || data.clustering

        analysis.clusteringEdit = AnalysisStepClustering.loadFromJson(clusteringEditData, analysis);

        analysis.compareGenes = AnalysisStepCompareGenes.loadFromJson(data.compare_genes, analysis);

        document.querySelector(UI.analysisWorkflowElt).classList.remove("is-hidden");

        return analysis;

    }

    /**
     * Loads preliminary figures for the analysis.
     * @async
     * @function loadPreliminaryFigures
     * @memberof Analysis
     * @instance
     * @throws {Error} If failed to access the dataset.
     * @returns {Promise<void>} A promise that resolves when the preliminary figures are loaded.
     */
    async loadPreliminaryFigures() {

        document.querySelector(UI.primaryInitialPlotContainer).classList.remove("is-hidden");
        try {
            const {data} = await axios.post("./cgi/h5ad_preview_primary_filter.cgi", convertToFormData({
                dataset_id: this.dataset.id,
                analysis_id: this.id,
                analysis_type: this.type,
                session_id: this.userSessionId
            }));

            document.querySelector(UI.primaryInitialLoadingPlotElt).classList.add("is-hidden");

            if (!data.success || data.success < 1) {
                document.querySelector(UI.primaryInitialViolinContainer).textContent = "Preliminary figures not yet generated. Continue your analysis.";
                createToast("Preliminary figures not found. You can still continue the analysis though.", "is-warning");
                return;
            }
            document.querySelector(UI.primaryInitialViolinContainer).innerHTML = `<a target="_blank" href="./datasets/${this.dataset.id}.prelim_violin.png"><img src="./datasets/${this.dataset.id}.prelim_violin.png" class="img-fluid img-zoomed" /></a>`;
            document.querySelector(UI.primaryInitialScatterContainer).innerHTML = `<a target="_blank" href="./datasets/${this.dataset.id}.prelim_n_genes.png"><img src="./datasets/${this.dataset.id}.prelim_n_genes.png" class="img-fluid img-zoomed" /></a>`;
            createToast("Preliminary plots displayed", "is-success");

        } catch (error) {
            createToast("Failed to access dataset");
            logErrorInConsole(`Failed ID was: ${datasetId} because msg: ${error.message}`);
        }
    }

    /**
     * Makes a public copy of the analysis.
     *
     * @async
     * @function makePublicCopy
     * @memberof analysis
     * @instance
     *
     * @throws {Error} If there is an error making the analysis public.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is successfully made public.
     */
    async makePublicCopy() {

        try {
            const data = await this.copyDatasetAnalysis('public');

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.type = 'user_saved';
            document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
            createToast("Changes made to this public analysis will create a local copy within your profile.", "is-info", true);
            document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
            document.querySelector(UI.btnMakePublicCopyElt).classList.add("is-hidden");
            document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.add("is-hidden");
            document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.add("is-hidden");

            await this.getSavedAnalysesList(this.dataset.id, this.id);
        } catch (error) {
            createToast(`Error making analysis public: ${error.message}`);
            logErrorInConsole(error);
        }

    }

    /**
     * Places an analysis image in the specified target element.
     *
     * @param {Object} options - The options for placing the analysis image.
     * @param {Object} options.params - The parameters to be sent with the request.
     * @param {string} options.title - The title of the image.
     * @param {string} options.target - The CSS selector for the target element.
     * @returns {Promise<void>} - A promise that resolves when the image is placed successfully.
     */
    async placeAnalysisImage({params, title, target} = {}) {
        const url = "./cgi/get_analysis_image.cgi";
        const response = await axios.get(url, { params });

        if (response.status === 200) {
            const imgSrc = response.request.responseURL;
            const html = `<a target="_blank" href="${imgSrc}"><img src="${imgSrc}" class="image" alt="${title}" /></a>`;
            document.querySelector(target).innerHTML = html;
            return;
        }

        console.error(`Error: ${response.status}`);
        createToast(`Error getting analysis image`);
    }

    /**
     * Resets the Analysis instance to its initial state.
     * This method starts over as a completely new Analysis instance, with a new ID,
     * and resets all existing components.
     */
    reset() {

        // Clear plot images and such
        for (const el of document.getElementsByClassName("js-resetable")) {
            el.replaceChildren();
        }

        // Hide success and failed state icons
        for (const el of document.getElementsByClassName("js-step-success")) {
            el.classList.add("is-hidden");
        }

        for (const el of document.getElementsByClassName("js-step-failure")) {
            el.classList.add("is-hidden");
        }

        // For each step, if it exists, call its reset method
        // This allows us to worry about
        for (const step of Object.values(this)) {
            if (step?.reset) {
                step.reset();
            }
        }

        // Set non-initial steps as unclickable
        document.querySelector(UI.qcByMitoSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.selectVariableGenesSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.pcaSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.tsneSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.clusteringSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.markerGenesSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.clusteringEditSection).classList.add("is-pointer-events-none");
        document.querySelector(UI.compareGenesSection).classList.add("is-pointer-events-none");


        this.datasetIsRaw = true;
        this.id = uuid();
        this.label = null;

        // Hide tSNE section
        document.querySelector(UI.labeledTsneSection).classList.add("is-hidden");

    }

    async save() {

        if (!this.userSessionId) {
            console.warn("Cannot save analysis without a user session ID");
            return;
        }

        // clone this object in such a way to not
        // include the "analysis" property for each step due to circular reference
        const clone = JSON.parse(JSON.stringify(this, (key, value) => {
            if (key === "analysis") {
                return undefined;
            }
            return value;
        }));

        // delete the "analysis" key from each step
        for (const step of Object.values(clone)) {
            if (step?.analysis) {
                delete step.analysis;
            }
        }

        // Save the current analysis parameters to disk
        // they must be in snake_case, the inverse of the loadFromJson method
        // Also do this for nested objects

        const state = Analysis.convertToJson(clone);
        for (const key in state) {
            if (typeof state[key] === 'object') {
                state[key] = Analysis.convertToJson(state[key]);
            }
        }

        // Some legacy things to change around
        state.dataset_id = state.dataset.id;
        if (state.qc_by_mito) {
            state.qc_by_mito.filter_mito_perc = state.qc_by_mito.filter_mito_percent;
            delete state.qc_by_mito.filter_mito_percent;
        }

        try {
            const {data} = await axios.post("./cgi/save_dataset_analysis.cgi", convertToFormData({
                session_id: this.userSessionId,
                dataset_id: this.dataset.id,
                analysis_id: this.id,
                analysis_type: this.type,
                analysis_vetting: this.vetting,
                label: this.label,
                state
            }));
            if ((!data.success) || (data.success < 1)) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            createToast("Analysis saved", "is-success");

            await this.getSavedAnalysesList(this.dataset.id, this.id);

            // Update the current analysis label
            document.querySelector(UI.currentAnalysisElt).textContent = this.label;

        } catch (error) {
            createToast(`Error saving analysis: ${error.message}`);
            logErrorInConsole(error);
        }
    }

    /**
     * Saves the analysis to the user's area.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is saved successfully.
     * @throws {Error} If there is an error saving the analysis.
     */
    async saveToUserArea() {

        try {
            const data = await this.copyDatasetAnalysis('user_saved');

            if ((!data.success) || (data.success < 1)) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.type = 'user_saved';
            document.querySelector(UI.btnSaveAnalysisElt).textContent = "Saved";
            document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
            createToast("This analysis is stored in your profile.", "is-info", true);
            document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
            document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.remove("is-hidden");
            document.querySelector(UI.btnMakePublicCopyElt).classList.remove("is-hidden");
            document.querySelector(UI.newAnalysisLabelContainer).classList.add("is-hidden");

            await this.getSavedAnalysesList(this.dataset.id, this.id);
        } catch (error) {
            createToast(`Error saving analysis: ${error.message}`);
        }
    }

    /**
     * Toggles the visibility of analysis buttons based on the analysis type.
     */
    showHideAnalysisButtons() {
        // Unhide delete buttons, which were hidden for new analyses
        document.querySelector(UI.btnDeleteSavedAnalysisElt).classList.remove("is-hidden");
        document.querySelector(UI.btnDeleteUnsavedAnalysisElt).classList.remove("is-hidden");

        if (this.type === 'user_unsaved') {
            document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
            document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");
        } else if (this.type === 'user_saved') {
            document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
            document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
        }
    }

}

/* Putting these in order of the workbench steps */

class AnalysisStepLabeledTsne {
    /* Special case... this only happens in primary analyses */
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    reset() {
        this.resetUI();
    }

    resetUI() {
        document.querySelector(UI.labeledTsneSection).classList.remove("is-hidden");
    }

    async runAnalysis() {

        // Get the tSNE config
        const dataset = this.analysis.dataset;
        const data = await getEmbeddedTsneDisplay(dataset.id);
        const config = data.plotly_config;

        const img = document.createElement('img');
        img.className = 'image'

        // Generate the tSNE plot
        try {
            const image = await getTsneImageData(document.querySelector(UI.labeledTsneGeneSymbolElt).value, config);
            if (typeof image === 'object' || typeof image === "undefined") {
                throw new Error("No image data returned");
            } else {
                img.src = `data:image/png;base64,${image}`;
                document.querySelector(UI.labeledTsnePlotContainer).appendChild(img);
            }

            createToast("Labeled tSNE plot generated", "is-success");

        } catch (error) {
            createToast(`Error generating tSNE plot: ${error.message}`);
            logErrorInConsole(error);
            failStepWithHref(UI.labeledTsneSection);
        }
    }

    updateUIWithResults() {
        // pass
    }
}


class AnalysisStepPrimaryFilter {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    async applyPrimaryFilter() {

        document.querySelector(UI.primaryFilterSectionFailedElt).classList.add("is-hidden");
        document.querySelector(UI.primaryFilterSectionSuccessElt).classList.add("is-hidden");

        createToast("Applying dataset filters", "is-info");
        for (const elt of document.querySelectorAll(UI.datasetInfoResetableElts)) {
            elt.replaceChildren();
        }

        // If a user redoes the primary filters, it's assumed they want to do this on the primary
        //  datasource again.  This allows them to lessen the stringency of a filter.
        const originalAnalysisType = this.analysis.type;
        this.analysis.type = 'primary';

        this.filterCellsLtNGenesSelected = document.querySelector(UI.filterCellsLtNGenesSelectedElt).checked;
        this.filterCellsLtNGenes = null;
        if (this.filterCellsLtNGenesSelected) {
            this.filterCellsLtNGenes = document.querySelector(UI.filterCellsLtNGenesElt).value;
        }

        this.filterCellsGtNGenesSelected = document.querySelector(UI.filterCellsGtNGenesSelectedElt).checked;
        this.filterCellsGtNGenes = null;
        if (this.filterCellsGtNGenesSelected) {
            this.filterCellsGtNGenes = document.querySelector(UI.filterCellsGtNGenesElt).value;
        }

        this.filterGenesLtNCellsSelected = document.querySelector(UI.filterGenesLtNCellsSelectedElt).checked;
        this.filterGenesLtNCells = null;
        if (this.filterGenesLtNCellsSelected) {
            this.filterGenesLtNCells = document.querySelector(UI.filterGenesLtNCellsElt).value;
        }

        this.filterGenesGtNCellsSelected = document.querySelector(UI.filterGenesGtNCellsSelectedElt).checked;
        this.filterGenesGtNCells = null;
        if (this.filterGenesGtNCellsSelected) {
            this.filterGenesGtNCells = document.querySelector(UI.filterGenesGtNCellsElt).value;
        }

        try {
            const {data} = await axios.post("./cgi/h5ad_apply_primary_filter.cgi", convertToFormData({
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                dataset_id: this.analysis.dataset.id,
                session_id: this.analysis.userSessionId,
                filter_cells_lt_n_genes: this.filterCellsLtNGenes || "",
                filter_cells_gt_n_genes: this.filterCellsGtNGenes || "",
                filter_genes_lt_n_cells: this.filterGenesLtNCells || "",
                filter_genes_gt_n_cells: this.filterGenesGtNCells || ""
            }));

            if (!data.success || data.success < 1) {
                let error = data.error || "Unknown error. Please contact gEAR support.";
                if (data.n_genes === 0) {
                    error = "Filter reduced genes to 0. Try less stringent cutoffs";
                } else if (data.n_obs === 0) {
                    error = "Filter reduced cells to 0. Try less stringent cutoffs";
                }
                throw new Error(error);
            }

            this.filteredGeneCount = data.n_genes;
            this.filteredCellCount = data.n_obs;

            // After this step, we are now working with a curated dataset
            this.analysis.type = originalAnalysisType;
            if(originalAnalysisType === 'primary') {
                this.analysis.type = 'user_unsaved';
            }

            this.calculated = true;
            this.updateUIWithResults();
            createToast("Data filters applied", "is-success");

        } catch (error) {
            createToast(`Error applying primary filters: ${error.message}`);
            logErrorInConsole(error);
            failStepWithHref("#primary-filter-s")
            document.querySelector(UI.primaryFilterSectionFailedElt).classList.remove("is-hidden");
        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }

    }

    /**
     * Loads data from a JSON object into an AnalysisStepPrimaryFilter instance.
     * @param {Object} data - The JSON object containing the data to be loaded.
     * @param {AnalysisStepPrimaryFilter} analysis - The AnalysisStepPrimaryFilter instance to load the data into.
     * @returns {AnalysisStepPrimaryFilter} - The loaded AnalysisStepPrimaryFilter instance.
     */
    static loadFromJson(data, analysis) {
        const step = new AnalysisStepPrimaryFilter(analysis);
        if (!data) return step;
        step.calculated = data['calculated'];

        step.filterCellsGtNGenes = data['filter_cells_gt_n_genes'];
        step.filterCellsGtNGenesSelected = data['filter_cells_gt_n_genes_selected'];
        step.filterCellsLtNGenes = data['filter_cells_lt_n_genes'];
        step.filterCellsLtNGenesSelected = data['filter_cells_lt_n_genes_selected'];

        step.filterGenesGtNCells = data['filter_genes_gt_n_cells'];
        step.filterGenesGtNCellsSelected = data['filter_genes_gt_n_cells_selected'];
        step.filterGenesLtNCells = data['filter_genes_lt_n_cells'];
        step.filterGenesLtNCellsSelected = data['filter_genes_lt_n_cells_selected'];

        step.filteredGeneCount = data['filtered_gene_count'];
        step.filteredCellCount = data['filtered_cell_count'];

        if (step.calculated ) {
            step.updateUIWithResults();
        }

        return step;
    }

    /**
     * Resets the state of the analysis object.
     * Sets all properties to their initial values and calls the `resetUI` method.
     */
    reset() {
        this.calculated = false;
        this.filterCellsGtNGenes = null;
        this.filterCellsGtNGenesSelected = false;
        this.filterCellsLtNGenes = 300;
        this.filterCellsLtNGenesSelected = false;

        this.filterGenesGtNCells = null;
        this.filterGenesGtNCellsSelected = false;
        this.filterGenesLtNCells = 3;
        this.filterGenesLtNCellsSelected = false;

        this.filteredGeneCount = null;
        this.filteredCellCount = null;
        this.resetUI();
    }

    /**
     * Resets the UI by setting the values of various filter inputs and checkboxes,
     * and modifying the visibility of certain elements.
     */
    resetUI() {
        document.querySelector(UI.primaryFilterSection).classList.remove("is-hidden");

        document.querySelector(UI.filterCellsGtNGenesElt).value = this.filterCellsGtNGenes;
        document.querySelector(UI.filterCellsLtNGenesElt).value = this.filterCellsLtNGenes
        document.querySelector(UI.filterGenesLtNCellsElt).value = this.filterGenesLtNCells;
        document.querySelector(UI.filterGenesGtNCellsElt).value = this.filterGenesGtNCells;

        document.querySelector(UI.filterCellsGtNGenesSelectedElt).checked = this.filterCellsGtNGenesSelected;
        document.querySelector(UI.filterCellsLtNGenesSelectedElt).checked = this.filterCellsLtNGenesSelected;
        document.querySelector(UI.filterGenesGtNCellsSelectedElt).checked = this.filterGenesGtNCellsSelected;
        document.querySelector(UI.filterGenesLtNCellsSelectedElt).checked = this.filterGenesLtNCellsSelected;

        for (const elt of document.querySelectorAll(UI.primaryInitialPlotElts)) {
            elt.classList.remove("is-hidden");
        }
        document.querySelector(UI.primaryTopGenesPlotContainer).classList.add("is-hidden");

        // hide previous shown elements
        document.querySelector(UI.selectedDatasetShapeFilteredContainer).classList.add("is-hidden");

        // show the instructions
        document.querySelector(UI.primaryFilterInstructionsElt).classList.remove("is-hidden");
    }

    updateUIWithResults(ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        if (this.calculated) {
            document.querySelector(UI.filterCellsLtNGenesSelectedElt).checked = false;
            if (this.filterCellsLtNGenesSelected) {
                document.querySelector(UI.filterCellsLtNGenesSelectedElt).checked = true;
            }

            document.querySelector(UI.filterCellsGtNGenesSelectedElt).checked = false;
            if (this.filterCellsGtNGenesSelected) {
                document.querySelector(UI.filterCellsGtNGenesSelectedElt).checked = true;
            }

            document.querySelector(UI.filterGenesLtNCellsSelectedElt).checked = false;
            if (this.filterGenesLtNCellsSelected) {
                document.querySelector(UI.filterGenesLtNCellsSelectedElt).checked = true;
            }

            document.querySelector(UI.filterGenesGtNCellsSelectedElt).checked = false;
            if (this.filterGenesGtNCellsSelected) {
                document.querySelector(UI.filterGenesGtNCellsSelectedElt).checked = true;
            }

            document.querySelector(UI.filterCellsLtNGenesElt).value = this.filterCellsLtNGenes || 300;
            document.querySelector(UI.filterCellsGtNGenesElt).value = this.filterCellsGtNGenes;
            document.querySelector(UI.filterGenesLtNCellsElt).value = this.filterGenesLtNCells || 3;
            document.querySelector(UI.filterGenesGtNCellsElt).value = this.filterGenesGtNCells;

            document.querySelector(UI.selectedDatasetShapeFilteredElt).textContent = `${this.filteredGeneCount} genes x ${this.filteredCellCount} obs`;
            document.querySelector(UI.selectedDatasetShapeFilteredContainer).classList.remove("is-hidden");

        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'highest_expr_genes',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            datetime: (new Date()).getTime()
        }

        ana.placeAnalysisImage(
            {'params': params, 'title': 'Highest expressed genes', 'target': UI.primaryTopGenesContainer});

        document.querySelector(UI.primaryTopGenesPlotContainer).classList.remove("is-hidden");

        // Now we can potentially save the analysis if it is a user one
        ana.showHideAnalysisButtons();


        passStepWithHref(UI.primaryFilterSection);
        openNextAnalysisStep([UI.qcByMitoSection], null, true);

        document.querySelector(UI.primaryFilterSectionSuccessElt).classList.remove("is-hidden");
    }
}

class AnalysisStepQCByMito {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    /**
     * Loads data from a JSON object and creates an instance of AnalysisStepQCByMito.
     *
     * @param {Object} data - The JSON object containing the data.
     * @param {Analysis} analysis - The analysis object.
     * @returns {AnalysisStepQCByMito} - The created instance of AnalysisStepQCByMito.
     */
    static loadFromJson(data, analysis) {
        const step = new AnalysisStepQCByMito(analysis);
        if (!data) return step;

        step.calculated = data['calculated'];
        step.genePrefix = data['gene_prefix'];
        step.filterMitoPercent = data['filter_mito_perc'];
        step.filterMitoCount = data['filter_mito_count'];
        step.nGenes = data['n_genes'];
        step.nObs = data['n_obs'];

        if (step.calculated === true) {
            document.querySelector(UI.qbmGenePrefixElt).value = step.genePrefix;
            document.querySelector(UI.qbmFilterMitoPercElt).value = step.filterMitoPercent;
            document.querySelector(UI.qbmFilterMitoCountElt).value = step.filterMitoCount;
            step.updateUIWithResults(true);
        }

        return step;
    }

    /**
     * Resets the analysis state and user interface.
     */
    reset() {
        this.calculated = false;
        this.genePrefix = "mt-";
        this.filterMitoPercent = null;
        this.filterMitoCount = null;
        this.nGenes = null;
        this.nObs = null;
        this.resetUI();
    }

    /**
     * Resets the user interface for analysis.
     */
    resetUI() {
        document.querySelector(UI.qcByMitoSection).classList.remove("is-hidden");

        document.querySelector(UI.qbmGenePrefixElt).value = this.genePrefix;
        disableAndHideElement(document.querySelector(UI.btnQbmSaveElt));
        disableAndHideElement(document.querySelector(UI.qbmSaveWarningElt));
        document.querySelector(UI.btnQbmSaveElt).textContent = 'Save these genes';

        document.querySelector(UI.btnDoAnalysisQcByMitoElt).classList.remove("is-hidden");

        // hide previous shown elements
        document.querySelector(UI.qbmPostShapeContainer).classList.add("is-hidden");
        // show the instructions
        document.querySelector(UI.qbmInstructionsElt).classList.remove("is-hidden");
    }

    /**
     * Runs the analysis.
     *
     * @param {boolean} saveDataset - Indicates whether to save the dataset.
     * @returns {Promise<void>} - A promise that resolves when the analysis is complete.
     * @throws {Error} - If there is an error during the analysis.
     */
    async runAnalysis(saveDataset) {
        document.querySelector(UI.btnQbmSaveElt).disabled = true;
        document.querySelector(UI.qbmSaveWarningElt).classList.add("is-hidden");

        // reset success and failure icons
        document.querySelector(UI.qcByMitoSectionSuccessElt).classList.add("is-hidden");
        document.querySelector(UI.qcByMitoSectionFailedElt).classList.add("is-hidden");

        if (Boolean(saveDataset)) {
            createToast("Applying mitochondrial filter", "is-info");
        } else {
            createToast("Analyzing mitochondrial genes", "is-info");
        }

        // Clear the images
        document.querySelector(UI.qbmViolinContainer).replaceChildren();
        document.querySelector(UI.qbmScatterPercentMitoContainer).replaceChildren();
        document.querySelector(UI.qbmScatterNGenesContainer).replaceChildren();

        this.genePrefix = document.querySelector(UI.qbmGenePrefixElt).value;
        this.filterMitoPercent = document.querySelector(UI.qbmFilterMitoPercElt).value;
        this.filterMitoCount = document.querySelector(UI.qbmFilterMitoCountElt).value;

        try {
            const {data} = await axios.post("./cgi/h5ad_qc_by_mito.cgi", convertToFormData({
                dataset_id: this.analysis.dataset.id,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                genes_prefix: this.genePrefix,
                filter_mito_perc: this.filterMitoPercent,
                filter_mito_count: this.filterMitoCount,
                save_dataset: saveDataset
            }));

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.nGenes = saveDataset ? data.n_genes.toString() : null;
            this.nObs = saveDataset ? data.n_obs.toString() : null;

            this.calculated = Boolean(saveDataset);
            document.querySelector(UI.btnQbmSaveElt).disabled = false;

            this.updateUIWithResults(this.calculated);
            createToast("Mitochondrial plot displayed", "is-success");

        } catch (error) {
            createToast(`Error doing QC analysis: ${error.message}`);
            logErrorInConsole(error);
            failStepWithHref(UI.qcByMitoSection);
            document.querySelector(UI.qcByMitoSectionFailedElt).classList.remove("is-hidden");
        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    /**
     * Updates the UI with the results of the analysis.
     * @param {boolean} [resultsSaved=false] - Indicates whether the results have been saved.
     * @param {Object} [ana=null] - The analysis object.
     */
    updateUIWithResults(resultsSaved=false, ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        enableAndShowElement(document.querySelector(UI.btnQbmSaveElt));
        document.querySelector(UI.qbmSaveWarningElt).classList.remove("is-hidden");

        if (resultsSaved) {
            document.querySelector(UI.btnDoAnalysisQcByMitoElt).classList.add("is-hidden");
            document.querySelector(UI.btnQbmSaveElt).disabled = true;
            document.querySelector(UI.btnQbmSaveElt).textContent = 'Saved';
            document.querySelector(UI.qbmPostShapeContainer).classList.remove("is-hidden");
        }

        document.querySelector(UI.qbmInstructionsElt).classList.add("is-hidden");
        document.querySelector(UI.qbmPostShapeElt).textContent = `${this.nGenes} genes x ${this.nObs} obs`;

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'violin_qc_by_mito',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.placeAnalysisImage(
            {'params': params, 'title': 'QC by mito - Violin',
             'target': UI.qbmViolinContainer});

        params['analysis_name'] = 'scatter_percent_mito'
        ana.placeAnalysisImage(
            {'params': params, 'title': 'QC by mito - Scatter percent mito',
             'target': UI.qbmScatterPercentMitoContainer});

        params['analysis_name'] = 'scatter_n_genes'
        ana.placeAnalysisImage(
            {'params': params, 'title': 'QC by mito - Scatter N genes',
             'target': UI.qbmScatterNGenesContainer});

        // Only pass the step if the results have been saved
        if (resultsSaved) {
            document.querySelector(UI.qcByMitoSectionSuccessElt).classList.remove("is-hidden");
            blockAnalysisStep(UI.primaryFilterSection);
            passStepWithHref(UI.qcByMitoSection);
            blockAnalysisStep(UI.qcByMitoSection);
            openNextAnalysisStep([UI.selectVariableGenesSection], null, true);
        }
    }
}


class AnalysisStepSelectVariableGenes {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    /**
     * Loads data from a JSON object and creates an instance of AnalysisStepSelectVariableGenes.
     *
     * @param {Object} data - The JSON object containing the data.
     * @param {Analysis} analysis - The analysis object.
     * @returns {AnalysisStepSelectVariableGenes} - The created instance of AnalysisStepSelectVariableGenes.
     */
    static loadFromJson(data, analysis) {
        const step = new AnalysisStepSelectVariableGenes(analysis);
        if (!data) return step;

        step.calculated = data['calculated'];
        step.normCountsPerCell = data['norm_counts_per_cell'];
        step.flavor = data['flavor'];
        step.nTopGenes = data['n_top_genes']
        step.minMean = data['min_mean'];
        step.maxMean = data['max_mean'];
        step.minDispersion = data['min_dispersion'];
        step.regressOut = data['regress_out'];
        step.scaleUnitVariance = data['scale_unit_variance'];

        if (step.calculated) {
            step.updateUIWithResults(true);
        }

        return step;
    }

    /**
     * Resets the analysis object to its initial state.
     */
    reset() {
        this.calculated = false;
        this.normCountsPerCell = "1e4";
        this.flavor = 'seurat';
        this.nTopGenes = null;
        this.nGenes = null;
        this.nObs = null;
        this.minMean = 0.0125;
        this.maxMean = 3;
        this.minDispersion = 0.5;
        this.regressOut = false; // no options to change
        this.scaleUnitVariance = false; // no options to change
        this.resetUI();
    }

    /**
     * Resets the UI elements to their default values.
     */
    resetUI() {
        document.querySelector(UI.selectVariableGenesSection).classList.remove("is-hidden");

        document.querySelector(UI.asvgNormCountsPerCellElt).value = this.normCountsPerCell;
        document.querySelector(UI.asvgFlavorElt).value = this.flavor;
        document.querySelector(UI.asvgNTopGenesElt).value = this.nTopGenes;
        document.querySelector(UI.asvgTopGenesListElt).replaceChildren();
        document.querySelector(UI.asvgMinMeanElt).value = this.minMean;
        document.querySelector(UI.asvgMaxMeanElt).value = this.maxMean;
        document.querySelector(UI.asvgMinDispersionElt).value = this.minDispersion;

        disableAndHideElement(document.querySelector(UI.btnAsvgSaveElt));
        document.querySelector(UI.btnAsvgSaveElt).textContent = 'Save these genes';
        document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).classList.remove("is-hidden");

        // Hide preious shown elements
        document.querySelector(UI.asvgPostShapeContainer).classList.add("is-hidden");
        document.querySelector(UI.asvgTopGenesContainer).classList.add("is-hidden");
        // show the instructions
        document.querySelector(UI.asvgInstructionsElt).classList.remove("is-hidden");
    }

    /**
     * Runs the analysis for identifying variable genes.
     *
     * @param {boolean} saveDataset - Indicates whether to save the dataset.
     * @returns {Promise<void>} - A promise that resolves when the analysis is complete.
     */
    async runAnalysis(saveDataset) {
        document.querySelector(UI.btnAsvgSaveElt).disabled = true;
        document.querySelector(UI.asvgSaveWarningElt).classList.add("is-hidden");

        if (Boolean(saveDataset)) {
            createToast("Saving variable genes", "is-info");
        } else {
            createToast("Analyzing variable genes", "is-info");
        }

        document.querySelector(UI.asvgPlotContainer).replaceChildren();

        this.normCountsPerCell = document.querySelector(UI.asvgNormCountsPerCellElt).value;
        this.flavor = document.querySelector(UI.asvgFlavorElt).value;
        this.nTopGenes = document.querySelector(UI.asvgNTopGenesElt).value;
        this.minMean = document.querySelector(UI.asvgMinMeanElt).value;
        this.maxMean = document.querySelector(UI.asvgMaxMeanElt).value;
        this.minDispersion = document.querySelector(UI.asvgMinDispersionElt).value;

        const params = {
            'dataset_id': this.analysis.dataset.id,
            'analysis_id': this.analysis.id,
            'analysis_type': this.analysis.type,
            'session_id': this.analysis.userSessionId,
            'norm_counts_per_cell': this.normCountsPerCell,
            'flavor': this.flavor,
            'n_top_genes': this.nTopGenes,
            'min_mean': this.minMean,
            'max_mean': this.maxMean,
            'min_dispersion': this.minDispersion,
            'regress_out': this.regressOut,
            'scale_unit_variance': this.scaleUnitVariance,
            'save_dataset': saveDataset
        }

        try {
            const {data} = await axios.post("./cgi/h5ad_identify_variable_genes.cgi", convertToFormData(params));

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.nGenes = saveDataset ? data.n_genes.toString() : null;
            this.nObs = saveDataset ? data.n_obs.toString() : null;

            this.calculated = Boolean(saveDataset);

            document.querySelector(UI.btnAsvgSaveElt).disabled = false;
            if (this.calculated) {
                document.querySelector(UI.btnAsvgSaveElt).disabled = true;
            }

            document.querySelector(UI.btnAsvgSaveElt).disabled = false;
            this.updateUIWithResults(this.calculated);
            createToast("Variable genes plot created", "is-success");

            document.querySelector(UI.asvgTopGenesListElt).textContent = `${data['top_genes']}`;
            // pre-populate these genes as PCA genes to plot
            document.querySelector(UI.pcaGenesToColorElt).value = data['top_genes'];
            document.querySelector(UI.asvgTopGenesContainer).classList.remove("is-hidden");

        } catch (error) {
            createToast(`Error identifying variable genes: ${error.message}`);
            logErrorInConsole(error);
            document.getElementById(UI.selectVariableGenesSectionFailedElt).classList.remove("is-hidden");
            failStepWithHref(UI.selectVariableGenesSection);
        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    /**
     * Updates the UI with the results of the analysis.
     *
     * @param {boolean} [resultsSaved=false] - Indicates whether the results have been saved.
     * @param {object} [ana=null] - The analysis object.
     */
    updateUIWithResults(resultsSaved=false, ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }


        document.querySelector(UI.asvgNormCountsPerCellElt).value = this.normCountsPerCell;
        document.querySelector(UI.asvgFlavorElt).value = this.flavor;
        document.querySelector(UI.asvgNTopGenesElt).value = this.nTopGenes;
        document.querySelector(UI.asvgMinMeanElt).value = this.minMean;
        document.querySelector(UI.asvgMaxMeanElt).value = this.maxMean;
        document.querySelector(UI.asvgMinDispersionElt).value = this.minDispersion;

        enableAndShowElement(document.querySelector(UI.btnAsvgSaveElt));
        document.querySelector(UI.asvgSaveWarningElt).classList.remove("is-hidden");

        if (resultsSaved) {
            document.querySelector(UI.btnAsvgSaveElt).disabled = true;
            document.querySelector(UI.btnAsvgSaveElt).textContent = 'Saved';
            document.querySelector(UI.asvgPostShapeContainer).classList.remove("is-hidden");
            document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).classList.add("is-hidden");
        }

        document.querySelector(UI.asvgPostShapeElt).textContent = `${this.nGenes} genes x ${this.nObs} obs`;


        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'filter_genes_dispersion',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,

            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.placeAnalysisImage(
            {'params': params, 'title': 'Variable genes', 'target': UI.asvgPlotContainer});

        document.querySelector(UI.asvgInstructionsElt).classList.add("is-hidden");

        // Only pass the step if the results have been saved
        if (resultsSaved) {
            document.querySelector(UI.selectVariableGenesSectionSuccessElt).classList.remove("is-hidden");
            passStepWithHref(UI.selectVariableGenesSection);
            blockAnalysisStep(UI.selectVariableGenesSection);
            openNextAnalysisStep([UI.pcaSection], null, true);
        }

    }
}
class AnalysisStepPCA {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    /**
     * Loads data from a JSON object and creates an instance of AnalysisStepPCA.
     *
     * @param {Object} data - The JSON object containing the data to load.
     * @param {Analysis} analysis - The analysis object to associate with the created instance.
     * @returns {AnalysisStepPCA} - The created instance of AnalysisStepPCA.
     */
    static loadFromJson(data, analysis) {
        const step = new AnalysisStepPCA(analysis);
        if (!data) return step;

        step.calculated = data['calculated'];
        step.genesToColor = data['genes_to_color'];

        if (step.calculated ) {
            step.updateUIWithResults();
        }

        return step;
    }

    /**
     * Resets the analysis state and user interface.
     */
    reset() {
        this.calculated = false;
        this.genesToColor = null;
        this.resetUI();
    }

    /**
     * Resets the user interface by clearing the values of the pcaGenesToColorElt and topPcaGenesElt elements.
     */
    resetUI() {
        document.querySelector(UI.pcaSection).classList.remove("is-hidden");

        document.querySelector(UI.pcaGenesToColorElt).value = this.genesToColor;
        document.querySelector(UI.topPcaGenesElt).value = null;

        // hide previously shown elements
        document.querySelector(UI.pcaMissingGeneContainer).classList.add("is-hidden");
        document.querySelector(UI.pcaPcToTopGenesContainer).classList.add("is-hidden");
        document.querySelector(UI.pcaGeneListContainer).classList.add("is-hidden");
        // show the instructions
        document.querySelector(UI.pcaInstructionsElt).classList.remove("is-hidden");

        // disable the "save as gene list" button
        document.querySelector(UI.btnSavePcaGeneListElt).disabled = true;
    }

    /**
     * Runs the analysis by computing principal components and variance plot.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is completed.
     * @throws {Error} If there is an error running the PCA.
     */
    async runAnalysis() {
        createToast("Computing principal components and variance plot", "is-info");
        document.querySelector(UI.pcaMissingGeneElt).classList.add("is-hidden");

        document.querySelector(UI.pcaSectionSuccessElt).classList.add("is-hidden");
        document.querySelector(UI.pcaSectionFailedElt).classList.add("is-hidden");

        document.querySelector(UI.pcaMissingGeneContainer).classList.add("is-hidden");

        for (const container of document.querySelectorAll(UI.pcaResetableElts)) {
            container.replaceChildren();
        }

        const computePCA = this.calculated ? false : true;

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_pca.cgi", convertToFormData({
                dataset_id: this.analysis.dataset.id,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                genes_to_color: document.querySelector(UI.pcaGenesToColorElt).value,
                compute_pca: computePCA
            }));
            if (!data.success || data.success < 1) {
                document.querySelector(UI.pcaMissingGeneElt).textContent = data['missing_gene'] || "";

                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = true;
            this.genesToColor = document.querySelector(UI.pcaGenesToColorElt).value;
            this.updateUIWithResults();
            createToast("PCA and variance computed", "is-success");

            // Reveals more UI form elements
            document.querySelector(UI.pcaPcToTopGenesContainer).classList.remove("is-hidden");
            document.querySelector(UI.pcaGeneListContainer).classList.remove("is-hidden");

        } catch (error) {
            createToast(`Error running PCA: ${error.message}`);
            logErrorInConsole(error);
            document.querySelector(UI.pcaMissingGeneContainer).classList.remove("is-hidden");
            failStepWithHref(UI.pcaSection);
            document.querySelector(UI.pcaSectionFailedElt).classList.remove("is-hidden");

        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    /**
     * Runs the analysis to compute the top genes for principal components.
     *
     * @async
     * @function runAnalysisTopGenes
     * @memberof analysis
     * @returns {Promise<void>} A Promise that resolves when the analysis is complete.
     * @throws {Error} If there is an error computing the top PCA genes.
     */
    async runAnalysisTopGenes() {
        createToast("Computing top genes for principal components", "is-info");
        document.querySelector(UI.pcaTopGenesPlotContainer).replaceChildren();

        try {

            const {data} = await axios.post("./api/analysis/plotTopGenesPCA", convertToFormData({
                dataset_id: this.analysis.dataset.id,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                pcs: document.querySelector(UI.topPcaGenesElt).value
            }));

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            const params = {
                'analysis_id': this.analysis.id,
                'analysis_name': 'pca_loadings',
                'analysis_type': this.analysis.type,
                'dataset_id': this.analysis.dataset.id,
                'session_id': this.analysis.userSessionId,
                datetime: (new Date()).getTime()
            }

            this.analysis.placeAnalysisImage({
                'params': params,
                'title': 'Top PCA Genes',
                'target': UI.pcaTopGenesPlotContainer
            });

            document.querySelector(UI.pcaTopGenesPlotContainer).classList.remove("is-hidden");
            createToast("Top PCA genes computed", "is-success");
        } catch (error) {
            createToast(`Error computing top PCA genes: ${error.message}`);
            logErrorInConsole(error);
        }
    }

    /**
     * Updates the UI with analysis images based on the given analysis object.
     * If no analysis object is provided, it uses the analysis object stored in the class.
     *
     * @param {Object} [ana] - The analysis object to update the UI with.
     */
    updateUIWithResults(ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        document.querySelector(UI.pcaInstructionsElt).classList.add("is-hidden");

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'pca',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.placeAnalysisImage(
            {'params': params, 'title': 'PCA scatter', 'target': UI.pcaScatterContainer});

        params['analysis_name'] = 'pca_variance_ratio';
        ana.placeAnalysisImage(
            {'params': params, 'title': 'PCA variance', 'target': UI.pcaVarianceContainer});

        // PCA calculation enables tSNE/UMAP
        document.querySelector(UI.pcaGenesToColorElt).value = this.genesToColor;
        document.querySelector(UI.tsneGenesToColorElt).value = this.genesToColor;   // obvious 1st choice
        document.querySelector(UI.pcaPcToTopGenesContainer).classList.remove("is-hidden");

        document.querySelector(UI.pcaSectionSuccessElt).classList.remove("is-hidden");
        passStepWithHref(UI.pcaSection);
        openNextAnalysisStep([UI.tsneSection], null, true);
    }
}

class AnalysisSteptSNE {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    static loadFromJson(data, analysis) {
        const step = new AnalysisSteptSNE(analysis);
        if (!data) return step;

        step.calculated = data['calculated'];
        step.neighborsCalculated = data['neighbors_calculated'];
        step.tsneCalculated = data['tsne_calculated'];
        step.umapCalculated = data['umap_calculated'];
        step.genesToColor = data['genes_to_color'];
        step.nPcs = data['n_pcs'];
        step.nNeighbors = data['n_neighbors'];
        step.randomState = data['random_state'];
        step.plotTsne = data['plot_tsne'];
        step.plotUmap = data['plot_umap'];

        if (step.calculated) {
            step.updateUIWithResults();
        }

        return step;
    }

    /**
     * Resets the state of the analysis object.
     */
    reset() {
        this.calculated = false;
        this.neighborsCalculated = false;
        this.tsneCalculated = false;
        this.umapCalculated = false;
        this.genesToColor = false;
        this.nPcs = 2;
        this.nNeighbors = 5;
        this.randomState = 2;   // no option to change
        this.useScaled = false; // no option to change
        this.plotTsne = 0;
        this.plotUmap = 0;
        this.resetUI();
    }

    /**
     * Resets the UI by clearing the values of various input elements.
     */
    resetUI() {
        document.querySelector(UI.tsneSection).classList.remove("is-hidden");

        document.querySelector(UI.tsneGenesToColorElt).value = '';
        document.querySelector(UI.dimReductionNNeighborsElt).value = this.nNeighbors;
        document.querySelector(UI.tsneNPcsElt).value = this.nPcs;
        document.querySelector(UI.dimReductionMethodTsneElt).checked = false;
        document.querySelector(UI.dimReductionMethodUmapElt).checked = true;

        // hide elements that were revealed by previous steps
        document.querySelector(UI.tsneMissingGeneContainer).classList.add("is-hidden");
        // show the instructions
        document.querySelector(UI.tsneInstructionsElt).classList.remove("is-hidden");
    }

    async runAnalysis() {
        createToast("Computing tSNE/UMAP and generating plot", "is-info");

        // Reset success and failure icons
        document.querySelector(UI.tsneSectionSuccessElt).classList.add("is-hidden");
        document.querySelector(UI.tsneSectionFailedElt).classList.add("is-hidden");

        document.querySelector(UI.tsneMissingGeneContainer).classList.add("is-hidden");

        for (const container of document.querySelectorAll(UI.tsneResetableElts)) {
            container.replaceChildren();
        }

        // Anytime we run this there are three things which might need to be computed depending on what
        //  has happened.  neighborhood, tSNE and UMAP need to be done if it's the first time or if
        //  the settings have changed.
        // Also, both tSNE and UMAP can be replotted (with different gene coloring) and not recomputed

        let computeNeighbors = 1;
        let computeTsne = 0;
        let computeUmap = 0;
        let plotTsne = 0;
        let plotUmap = 0;

        const tsneChecked = document.querySelector(UI.dimReductionMethodTsneElt).checked;
        const umapChecked = document.querySelector(UI.dimReductionMethodUmapElt).checked;

        if (tsneChecked) {
            computeTsne = 1;
            plotTsne = 1;
        }

        if (umapChecked) {
            computeUmap = 1;
            plotUmap = 1;
        }

        // We don't have to recompute if none of the plotting parameters have changed.  This allows
        //  us to just do something like recolor.
        if (this.neighborsCalculated) {
            if (this.nPcs === document.querySelector(UI.tsneNPcsElt).value &&
                this.nNeighbors === document.querySelector(UI.dimReductionNNeighborsElt).value &&
                this.randomState === document.querySelector(UI.tsneRandomStateElt).value) {
                computeNeighbors = 0;
            }
        }

        // don't recompute tSNE if we already have and nothing has changed
        if (this.tsneCalculated) {
            if (this.nPcs === document.querySelector(UI.tsneNPcsElt).value &&
                this.randomState === document.querySelector(UI.tsneRandomStateElt).value) {
                computeTsne = 0;
            }
        }

        // Don't recompute UMAP if we already have
        if (this.umapCalculated) {
            computeUmap = 0;
        }

        // If plotting of either is requested and neighbors are being recomputed, we need to
        //  also recompute tSNE/UMAP
        if (computeNeighbors === 1) {
            if (plotUmap === 1) {
                computeUmap = 1;
            }

            if (plotTsne === 1) {
                computeTsne = 1;
            }
        }

        const params = {
            'dataset_id': this.analysis.dataset.id,
            'analysis_id': this.analysis.id,
            'analysis_type': this.analysis.type,
            'session_id': this.analysis.userSessionId,
            'genes_to_color': document.querySelector(UI.tsneGenesToColorElt).value,
            'n_pcs': document.querySelector(UI.tsneNPcsElt).value,
            'n_neighbors': document.querySelector(UI.dimReductionNNeighborsElt).value,
            'random_state': this.randomState,
            'use_scaled': this.useScaled,
            'compute_neighbors': computeNeighbors,
            'compute_tsne': computeTsne,
            'compute_umap': computeUmap,
            'plot_tsne': plotTsne,
            'plot_umap': plotUmap
        }

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_tsne.cgi", convertToFormData(params));

            if (!data.success || data.success < 1) {
                document.querySelector(UI.tsneMissingGeneContainer).textContent = data['missing_gene'] || "";
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = true;
            this.neighborsCalculated = true;
            this.tsneCalculated = Boolean(computeTsne);
            this.umapCalculated = Boolean(computeUmap);
            this.genesToColor = document.querySelector(UI.tsneGenesToColorElt).value;
            this.plotTsne = plotTsne;
            this.plotUmap = plotUmap;
            this.updateUIWithResults();

            if (computeTsne === 1) {
                this.tsneCalculated = 1;
            }
            if (computeUmap === 1) {
                this.umapCalculated = 1;
            }

            createToast("tSNE/UMAP computed and displayed", "is-success");

        } catch (error) {
            createToast(`Error generating tSNE: ${error.message}`);
            logErrorInConsole(error);
            failStepWithHref(UI.tsneSection);
            document.querySelector(UI.tsneSectionFailedElt).classList.remove("is-hidden");

            document.querySelector(UI.tsneMissingGeneContainer).classList.remove("is-hidden");

        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    updateUIWithResults(ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        document.querySelector(UI.tsneInstructionsElt).classList.add("is-hidden");

        if (!(this.neighborsCalculated || this.tsneCalculated || this.umapCalculated)) {
            console.info("No tSNE or UMAP calculated yet so cannot update UI.")
            return;
        }
        document.querySelector(UI.tsneGenesToColorElt).value = this.genesToColor;
        document.querySelector(UI.dimReductionNNeighborsElt).value = this.nNeighbors;
        document.querySelector(UI.tsneNPcsElt).value = this.nPcs;

        document.querySelector(UI.dimReductionMethodTsneElt).checked = false;
        if (Boolean(this.plotTsne)) {
            document.querySelector(UI.dimReductionMethodTsneElt).checked = true;
        }

        document.querySelector(UI.dimReductionMethodUmapElt).checked = false;
        if (Boolean(this.plotUmap)) {
            document.querySelector(UI.dimReductionMethodUmapElt).checked = true;
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'tsne',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        if (Boolean(this.plotTsne)) {
            ana.placeAnalysisImage(
                {'params': params, 'title': 'tSNE plot', 'target': UI.tsnePlotContainer});
        }

        if (Boolean(this.plotUmap)) {
            params['analysis_name'] = 'umap'
            ana.placeAnalysisImage(
                {'params': params, 'title': 'UMAP plot', 'target': UI.umapPlotContainer});
        }

        // dimensionality reduction enables clustering
        document.querySelector(UI.tsneSectionSuccessElt).classList.remove("is-hidden");
        passStepWithHref(UI.tsneSection);
        openNextAnalysisStep([UI.clusteringSection], null, true);

    }
}
class AnalysisStepClustering {
    constructor(analysis, mode="initial") {
        this.reset();   // TODO: split into two potentially for each step

        this.analysis = analysis;
        this.mode = mode;  // initial, edit
        if (!(["initial", "edit"].includes(mode))) {
            logErrorInConsole("Invalid mode for AnalysisStepClustering. Defaulting to 'initial'.");
            this.mode = "initial";
        }
        if (this.mode === "edit") {
            // Don't run clustering again if we are just editing the labels
            this.calculated = true;
        }

    }

    /**
     * Loads data from a JSON object and creates a new instance of AnalysisStepClustering.
     * @param {Object} data - The JSON object containing the data to load.
     * @param {Analysis} analysis - The analysis object to associate with the new instance.
     * @returns {AnalysisStepClustering} - The newly created instance of AnalysisStepClustering.
     */
    static loadFromJson(data, analysis) {
        const mode = data?.mode || "initial";
        const step = new AnalysisStepClustering(analysis, mode);
        if (!data) return step;

        step.calculated = data['calculated'];
        step.resolution = data['resolution'];
        step.plotTsne = data['plot_tsne'];
        step.plotUmap = data['plot_umap'];

        if (step.calculated ) {
            step.updateUIWithResults();
        }

        return step;
    }

    /**
     * Resets the state of the analysis object.
     */
    reset() {
        this.calculated = false;
        this.resolution = 0.5;
        this.plotUmap = 0;
        this.plotTsne = 0;
        this.resetUI();
    }

    /**
     * Resets the user interface for analysis.
     *
     * This method resets the resolution input back to the default value of 1.3 for all elements
     * specified in the `UI.resolutionElts` selector. It also adds a click listener to sync the
     * resolution inputs.
     */
    resetUI() {

        document.querySelector(UI.clusteringSection).classList.remove("is-hidden");
        document.querySelector(UI.clusteringEditSection).classList.remove("is-hidden");

        // Reset resolution input back to default
        document.querySelector(UI.resolutionElt).value = this.resolution;

        // Hide elements that were revealed by previous steps
        this.moda === "edit" && document.querySelector(UI.groupLabelsContainer).classList.add("is-hidden");

        // Show the instructions
        this.mode === "initial" && document.querySelector(UI.clusteringInstructionsElt).classList.remove("is-hidden");
        this.mode === "edit" && document.querySelector(UI.clusteringEditInstructionsElt).classList.remove("is-hidden");
    }

    /**
     * Runs the analysis by computing Louvain clusters and updating the UI with the results.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is completed.
     * @throws {Error} If there is an error generating clusters.
     */
    async runAnalysis() {
        // Reset success and failed icons
        if (this.mode === "initial") {
            document.querySelector(UI.clusteringSectionSuccessElt).classList.add("is-hidden");
            document.querySelector(UI.clusteringSectionFailedElt).classList.add("is-hidden");
        } else if (this.mode === "edit") {
            document.querySelector(UI.clusteringEditSectionSuccessElt).classList.add("is-hidden");
            document.querySelector(UI.clusteringEditSectionFailedElt).classList.add("is-hidden");
        }

        createToast("Computing clusters using Leiden algorithm", "is-info");

        const uiToReset = (this.mode === "initial") ? UI.clusteringResetableElts : UI.clusteringEditResetableElts;

        for (const container of document.querySelectorAll(uiToReset)) {
            container.replaceChildren();
        }

        const hasMatchingClusteringParams = this.resolution === document.querySelector(UI.resolutionElt).value;

        let computeClustering = true;
        let resolution = document.querySelector(UI.resolutionElt).value;

        const oldLabels = [...this.analysis.groupLabels];  // shallow-copy
        const newLabels = [];
        const keptLabels = [];

        // If the clustering params were changed, we need to allow for the marker genes to be recalculated
        if (hasMatchingClusteringParams && this.mode === "initial") {
            this.analysis.markerGenes.calculated = false;
        }

        const clusterInfo = [];

        // It is not safe to reuse group labels if the clustering params were changed
        if (this.mode === "edit") {
            computeClustering = false;  // we are just updating the labels
            resolution = this.resolution;   // we are not changing the resolution

            // Get the new labels and whether they are kept (checked)
            for (const glElt of document.querySelector(UI.clusterGroupLabelsTableBodyElt).querySelectorAll(".group-user-label input")) {
                newLabels.push(glElt.value);
                const thisRow = glElt.closest("tr");
                const thisCheck = thisRow.querySelector(".group-keep-chk input");
                keptLabels.push(thisCheck.checked);
            }

            // Create the cluster info object
            this.analysis.groupLabels.forEach((v, i) => {
                clusterInfo.push({
                    "old_label": oldLabels[i]
                    , "new_label": newLabels[i]
                    , "keep": keptLabels[i]
                })
            });

        }

        const plotTsne = (document.querySelector(UI.dimReductionMethodTsneElt).checked ? 1 : 0);
        const plotUmap = (document.querySelector(UI.dimReductionMethodUmapElt).checked ? 1 : 0);

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_clusters.cgi", convertToFormData({
                dataset_id: this.analysis.dataset.id,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                resolution,
                compute_clusters: computeClustering,
                plot_tsne: plotTsne,
                plot_umap: plotUmap,
                cluster_info: JSON.stringify(clusterInfo)
            }));

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = true;
            this.resolution = document.querySelector(UI.resolutionElt).value;
            this.plotTsne = plotTsne;
            this.plotUmap = plotUmap;
            this.updateUIWithResults();

            if (data["group_labels"].length) {
                // Update the group labels for the analysis, marker genes, and gene comparison
                this.analysis.groupLabels = []
                if (this.mode === "edit") {
                    this.analysis.markerGenes.populateClusterEditLabels(data.group_labels)
                }

                // Now update the labels so they work with gene comparison
                this.analysis.groupLabels = data['group_labels'].map(x => x.genes);
                this.analysis.compareGenes.populateGroupSelectors(this.analysis.groupLabels);
            }

            createToast("Louvain clusters computed", "is-success");
            // TODO:  $('#clustering_options_c .js-next-step').show();  // auto-select the next collapsable dropdown


        } catch (error) {
            createToast(`Error generating clusters: ${error.message}`);
            logErrorInConsole(error);
            if (this.mode === "initial") {
                failStepWithHref(UI.clusteringSection);
                document.querySelector(UI.clusteringSectionFailedElt).classList.remove("is-hidden");
            } else if (this.mode === "edit") {
                failStepWithHref(UI.clusteringEditSection);
                document.querySelector(UI.clusteringEditSectionFailedElt).classList.remove("is-hidden");
            }

        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    /**
     * Updates the UI with the results of the analysis.
     *
     * @param {Object} [ana=null] - The analysis object. If not provided, the function uses the analysis object of the class.
     * @returns {void}
     */
    updateUIWithResults(ana=null) {

        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        const toolInstructions = (this.mode === "initial") ? UI.clusteringInstructionsElt : UI.clusteringEditInstructionsElt;
        document.querySelector(toolInstructions).classList.add("is-hidden");

        if (this.calculated) {
            document.querySelector(UI.resolutionElt).value = this.resolution;
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'tsne_clustering',
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        // Need to ensure targets are different for the two different clustering steps
        let tsneTarget = UI.clusteringTsnePlotElt;
        let umapTarget = UI.clusteringUmapPlotElt;
        if (this.mode === "edit") {
            tsneTarget = UI.clusteringTsnePlotEditElt;
            umapTarget = UI.clusteringUmapPlotEditElt;
        }

        if (Boolean(this.plotTsne)) {
            ana.placeAnalysisImage(
                {'params': params, 'title': 'Cluster groups', 'target': tsneTarget});
        }

        if (Boolean(this.plotUmap)) {
            params['analysis_name'] = 'umap_clustering'
            ana.placeAnalysisImage(
                {'params': params, 'title': 'Cluster groups', 'target': umapTarget});
        }

        if (this.mode === "initial") {
            document.querySelector(UI.clusteringSectionSuccessElt).classList.remove("is-hidden");
            passStepWithHref(UI.clusteringSection);
            openNextAnalysisStep([UI.markerGenesSection], null, true);
        } else if (this.mode === "edit") {
            document.querySelector(UI.clusteringEditSectionSuccessElt).classList.remove("is-hidden");
            passStepWithHref(UI.clusteringEditSection);
            // block all previous steps
            blockAnalysisStep(UI.pcaSection);
            blockAnalysisStep(UI.tsneSection);
            blockAnalysisStep(UI.clusteringSection);
            blockAnalysisStep(UI.markerGenesSection);
            openNextAnalysisStep([UI.compareGenesSection], null, true);
        }
    }
}
class AnalysisStepMarkerGenes {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    /**
     * Adds or removes gene symbols to/from the genes of interest.
     * @param {Array<string>} geneSymbols - An array of gene symbols to add or remove.
     * @param {boolean} doAdd - A boolean indicating whether to add or remove the gene symbols.
     */
    addRemoveGenesOfInterest(geneSymbols, doAdd) {
        for (const gene of geneSymbols) {
            if (doAdd) {
                this.genesOfInterest.add(gene.trim());
            } else {
                this.genesOfInterest.delete(gene.trim());
            }
        }

        document.querySelector(UI.btnVisualizeMarkerGenesElt).disabled = !this.genesOfInterest.size;

    }

    /**
     * Downloads marker genes as an Excel file.
     */
    downloadMarkerGenes() {

        // Do the header row
        let row = [];

        for (const elt in document.querySelectorAll(UI.markerGenesTableHeadCellElts)) {
            row.push(elt.textContent);
        }
        const fileContentsheaders = row.join("\t") + "\n";

        // Now all the other rows
        let fileContents = fileContentsheaders;
        for (const elt in document.querySelectorAll(UI.markerGenesTableBodyElts)) {
            row = [];
            for (const cell in elt.querySelectorAll("td")) {
                row.push(cell.textContent);
            }
            fileContents += row.join("\t") + "\n";
        }

        // Now download the file
        const element = document.createElement('a');
        element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(fileContents));
        element.setAttribute('download', 'marker_genes.xls');
        element.classList.add('is-hidden');
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
    }

    /**
     * Fetches marker genes for the analysis.
     * @returns {Promise<Object>} The response data containing the marker genes.
     * @throws {Error} If the request fails or returns an error.
     */
    async fetchMarkerGenes() {
        const {data} = await axios.post("./cgi/h5ad_find_marker_genes.cgi", convertToFormData({
            'dataset_id': this.analysis.dataset.id,
            'analysis_id': this.analysis.id,
            'analysis_type': this.analysis.type,
            'session_id': this.analysis.userSessionId,
            'n_genes': this.nGenes,
            'compute_marker_genes': this.computeMarkerGenes
        }));

        if ((!data.success) || (data.success < 1)) {
            const error = data.error || "Unknown error. Please contact gEAR support.";
            throw new Error(error);
        }
        return data;
    }

    static async loadFromJson(data, analysis) {
        const step = new AnalysisStepMarkerGenes(analysis);
        if (!data) return step;

        step.calculated = data['calculated']
        step.nGenes = data['n_genes']
        step.groupLabels = data['group_labels']
        step.clusterLabel = data['cluster_label'] || "louvain";

        if (step.calculated ) {
            await step.updateUIWithResults(data);
        }

        return step;
    }

    /**
     * Performs marker gene visualization by making a POST request to the server and updating the UI with the results.
     * @async
     * @function performMarkerGeneVisualization
     * @memberof analysis
     * @throws {Error} If there is an error during the visualization process.
     */
    async performMarkerGeneVisualization() {

        document.querySelector(UI.markerGenesDotplotContainer).replaceChildren();
        document.querySelector(UI.markerGenesViolinContainer).replaceChildren();

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_marker_gene_visualization.cgi", convertToFormData({
                'dataset_id': this.analysis.dataset.id,
                'analysis_id': this.analysis.id,
                'analysis_type': this.analysis.type,
                'session_id': this.analysis.userSessionId,
                'marker_genes': JSON.stringify([...this.genesOfInterest])
            }));

            if ((!data.success) || (data.success < 1)) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            const params = {
                'analysis_id': this.analysis.id,
                'analysis_name': 'dotplot_goi',
                'analysis_type': this.analysis.type,
                'dataset_id': this.analysis.dataset.id,
                'session_id': this.analysis.userSessionId,
                // this saves the user from getting a cached image each time
                datetime: (new Date()).getTime()
            }

            this.analysis.placeAnalysisImage(
                {'params': params, 'title': 'Marker genes (dotplot)', 'target': UI.markerGenesDotplotContainer}
            );

            params['analysis_name'] = 'stacked_violin_goi';
            this.analysis.placeAnalysisImage(
                {'params': params, 'title': 'Marker genes (stacked violin)', 'target': UI.markerGenesViolinContainer}
            );

        } catch (error) {
            createToast("Error visualizing marker genes");
            logErrorInConsole(error);
        }
    }

    /**
     * Populates the marker genes labels in the analysis.
     *
     * @param {Array} groupLabels - The array of group labels.
     */
    populateClusterEditLabels(groupLabels) {
        this.groupLabels = [];

        // If the user has saved labels before, put them in the table here.  Else leave it
        //  as the top gene.
        let i = 0;
        if (this.analysis.groupLabels.length > 0) {
            for (i=0; i < this.analysis.groupLabels.length; i++) {
                groupLabels[i]['new_group_label'] = this.analysis.groupLabels[i];
                this.groupLabels.push(this.analysis.groupLabels[i]);
            }
        } else {
            for (i=0; i < groupLabels.length; i++) {
                groupLabels[i].new_group_label = groupLabels[i].genes;
                // For the overall labels, do the group number rather than gene since that's what's
                //  displayed by scanpy in the images
                this.groupLabels.push(i);
            }
        }

        document.querySelector(UI.clusterGroupLabelsTableBodyElt).replaceChildren();

        // show the abbreviated table in the louvain analysis block
        for (const group of groupLabels) {
            const clusterGroupLabelsHtml = document.querySelector(UI.clusterGroupLabelsTmpl).content.cloneNode(true);
            clusterGroupLabelsHtml.querySelector(".group-orig-label").textContent = group.group_label;
            clusterGroupLabelsHtml.querySelector(".group-num-cells").textContent = group.num_cells;
            clusterGroupLabelsHtml.querySelector(".group-marker").textContent = group.genes;
            clusterGroupLabelsHtml.querySelector(".group-user-label input").value = group.new_group_label;

            document.querySelector(UI.clusterGroupLabelsTableBodyElt).appendChild(clusterGroupLabelsHtml);
        }
        document.querySelector(UI.groupLabelsContainer).classList.remove("is-hidden");
    }

    /**
     * Populates the marker genes table with data retrieved from the server.
     *
     * @returns {Promise<void>} - A promise that resolves when the marker genes table is populated.
     * @throws {Error} - If there is an error retrieving the marker genes data.
     */
    populateMarkerGenesTable(table) {
        try {
            // Add the first th element (empty)
            const markerGenesHeaderHtml = document.querySelector(UI.markerGenesTableHeadTmpl).content.cloneNode(true);
            markerGenesHeaderHtml.querySelector("th").textContent = "";
            document.querySelector(UI.markerGenesTableHeadRowElt).appendChild(markerGenesHeaderHtml);

            // add the table header
            for (const column of table.columns) {
                const markerGenesHeaderHtml = document.querySelector(UI.markerGenesTableHeadTmpl).content.cloneNode(true);
                markerGenesHeaderHtml.querySelector("th").textContent = column.label;
                document.querySelector(UI.markerGenesTableHeadRowElt).appendChild(markerGenesHeaderHtml);
            }

            // add the table rows
            for (const row of table.rows) {
                // Create the row label
                const markerGenesRowHtml = document.querySelector(UI.markerGenesTableRowTmpl).content.cloneNode(true);
                markerGenesRowHtml.querySelector(".js-row-idx").textContent = row.rowid
                const currRow = markerGenesRowHtml.querySelector("tr");

                // create new td objects for each column

                for (const column of row.columns) {
                    const cellElt = document.createElement("td");
                    cellElt.textContent = column.label
                    currRow.appendChild(cellElt)
                }
                document.querySelector(UI.markerGenesTableBodyElt).appendChild(markerGenesRowHtml);


            }

            // Show the table and associated stuff
            document.querySelector(UI.markerGenesTableContainer).classList.remove("is-hidden");
        } catch (error) {
            createToast(`Error getting marker genes: ${error.message}`);
            logErrorInConsole(error);
        }

    }

    /**
     * Resets the analysis state.
     */
    reset() {
        this.calculated = false;
        this.clusterLabel = "louvain";
        this.computeMarkerGenes = !this.calculated;
        this.genesOfInterest = new Set();
        this.groupLabels = [];
        this.nGenes = 5;
        this.resetUI();
    }

    /**
     * Resets the user interface.
     */
    resetUI() {
        document.querySelector(UI.markerGenesSection).classList.remove("is-hidden");


        document.querySelector(UI.markerGenesNGenesElt).value = 5;

        // Affects the clustering block
        document.querySelector(UI.groupLabelsContainer).classList.add("is-hidden");

        // hide elements that were revealed by previous steps
        document.querySelector(UI.markerGenesVisualizationContainer).classList.add("is-hidden");
        document.querySelector(UI.markerGenesListContainer).classList.add("is-hidden");

        // show the instructions
        document.querySelector(UI.markerGenesInstructionsElt).classList.remove("is-hidden");

        // disable the "save as gene list" button and "visualize" button
        document.querySelector(UI.btnSaveMarkerGeneListElt).disabled = true;
        document.querySelector(UI.btnVisualizeMarkerGenesElt).disabled = true;
    }

    /**
     * Runs the analysis to compute marker genes.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is complete.
     * @throws {Error} If there is an error computing marker genes.
     */
    async runAnalysis() {
        createToast("Computing marker genes", "is-info");
        document.querySelector(UI.markerGenesPlotContainer).replaceChildren();

        document.querySelector(UI.markerGenesSectionSuccessElt).classList.add("is-hidden");
        document.querySelector(UI.markerGenesSectionFailedElt).classList.add("is-hidden");

        // Clear the table
        document.querySelector(UI.markerGenesTableHeadRowElt).replaceChildren();
        document.querySelector(UI.markerGenesTableBodyElt).replaceChildren();

        document.querySelector(UI.markerGenesManuallyEnteredElt).value = '';
        document.querySelector(UI.markerGenesSelectedCountElt).textContent = 0;
        document.querySelector(UI.markerGenesEnteredCountElt).textContent = 0;
        document.querySelector(UI.markerGenesUniqueCountElt).textContent = 0;

        this.genesOfInterest = new Set();
        this.clickedMarkerGenes = new Set();
        this.enteredMarkerGenes = new Set();

        // If marker genes are run, need to ensure the cluster label is populated with the "marker gene" label.
        this.analysis.groupLabels = [];

        this.computeMarkerGenes = true;
        if (this.nGenes === document.querySelector(UI.markerGenesNGenesElt).value) {
            this.computeMarkerGenes = false;
        }

        this.nGenes = document.querySelector(UI.markerGenesNGenesElt).value;

        try {
            const data = await this.fetchMarkerGenes();

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.clusterLabel = data.cluster_label;

            this.calculated = true;
            this.updateUIWithResults(data);
            this.analysis.groupLabels = data.group_labels.map(x => x.group_label);

        } catch (error) {
            createToast(`Error computing marker genes: ${error.message}`);
            logErrorInConsole(error);
            document.querySelector(UI.markerGenesSectionFailedElt).classList.remove("is-hidden");

            // mark in stepper
            failStepWithHref(UI.markerGenesSection)

        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    /**
     * Updates the UI with the results of the analysis.
     *
     * @param {Object} data - The data containing the analysis results.
     * @param {Object} [ana=null] - The analysis object. If not provided, the method uses the analysis object of the class instance.
     * @returns {void}
     */
    async updateUIWithResults(data, ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': `rank_genes_groups_${data.cluster_label}`,
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        document.querySelector(UI.markerGenesNGenesElt).value = this.nGenes;

        ana.placeAnalysisImage(
            {'params': params, 'title': 'Marker genes', 'target': UI.markerGenesPlotContainer});

        // May need to fetch marker gene data to build table if loading from JSON
        this.computeMarkerGenes = true;
        data = data?.table ? data : await this.fetchMarkerGenes();
        this.computeMarkerGenes = false;

        this.populateMarkerGenesTable(data.table);

        const groupLabels = data.group_labels.map(x => x.group_label);

        document.querySelector(UI.markerGenesVisualizationContainer).classList.remove("is-hidden");
        document.querySelector(UI.markerGenesListContainer).classList.remove("is-hidden");

        // marker gene calculation enables cluster comparison
        ana.compareGenes.populateGroupSelectors(groupLabels);

        // mark success
        document.querySelector(UI.markerGenesSectionSuccessElt).classList.remove("is-hidden");


        const nextSteps = [UI.compareGenesSection]

        // This can only be done with a non-primary analysis
        if (ana.type != "primary") {
            this.populateClusterEditLabels(data.group_labels);
            nextSteps.push(UI.clusteringEditSection)
        }

        // move stepper to next step
        passStepWithHref(UI.markerGenesSection)
        openNextAnalysisStep(nextSteps, UI.compareGenesSection)

    }
}

class AnalysisStepCompareGenes {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
    }

    /**
     * Loads an instance of AnalysisStepCompareGenes from JSON data.
     *
     * @param {Object} data - The JSON data to load from.
     * @returns {AnalysisStepCompareGenes} - The loaded AnalysisStepCompareGenes instance.
     */
    static loadFromJson(data, analysis) {
        const step = new AnalysisStepCompareGenes(analysis);
        if (!data) {
            data = {
                "calculated": false,
                "n_genes": 0,
                "query_cluster": null,
                "reference_cluster": null,
                "cluster_label": "louvain",
                "method": 't-test_overestim_var',
                "corr_method": 'benjamini-hochberg',
                "table_json_f": null,
                "table_json_r": null
            };
        }

        step.calculated = data['calculated'];
        step.nGenes = data['n_genes'];
        step.queryCluster = data['query_cluster'];
        step.referenceCluster = data['reference_cluster'];
        step.clusterLabel = data['cluster_label'];
        step.method = 't-test_overestim_var';
        step.corrMethod = 'benjamini-hochberg';

        if (data.hasOwnProperty('cluster_label')) {
            step.clusterLabel = data['metcluster_labelhod'];
        }

        if (data.hasOwnProperty('method')) {
            step.method = data['method'];
        }

        if (data.hasOwnProperty('corr_method')) {
            step.corrMethod = data['corr_method'];
        }

        if (data.hasOwnProperty('table_json')) {
            step.tableJson = data['table_json'];
        }

        if (step.calculated) {
            step.updateUIWithResults(data);
        }

        return step;

    }

    /**
     * Populates a comparison table with data.
     *
     * @param {string} tableId - The ID of the table element.
     * @param {Array<Array<string>>} tableJson - The JSON data representing the table rows and cells.
     */
    populateComparisonTable(tableId, tableJson) {
        const tableElt = document.querySelector(`${tableId} tbody`);
        tableElt.replaceChildren();

        for (const row of tableJson) {
            const rowElt = document.createElement("tr");
            // Append the key as the first cell
            const keyElt = document.createElement("th");
            keyElt.textContent = row[0];

            for (const cell of row) {
                const cellElt = document.createElement("td");
                cellElt.textContent = cell;
                rowElt.appendChild(cellElt);
            }
            tableElt.appendChild(rowElt);
        }
    }

    /**
     * Populates the group selectors with the given group labels.
     *
     * @param {Array<string>} groupLabels - The labels of the groups.
     */
    populateGroupSelectors(groupLabels) {

        // Clear the options
        document.querySelector(UI.queryClusterOptionsElt).replaceChildren();
        document.querySelector(UI.referenceClusterOptionsElt).replaceChildren();

        for (const label of groupLabels) {

            const option = document.createElement("option");
            option.textContent = label;
            option.value = label;

            document.querySelector(UI.queryClusterOptionsElt).appendChild(option);
            document.querySelector(UI.referenceClusterOptionsElt).appendChild(option.cloneNode(true));
        }

        document.querySelector(UI.queryClusterSelectElt).value = this.queryCluster;
        document.querySelector(UI.referenceClusterSelectElt).value = this.referenceCluster;
    }

    /**
     * Resets the analysis object to its initial state.
     */
    reset() {
        this.calculated = false;
        this.nGenes = 0;
        this.queryCluster = null;
        this.referenceCluster = null;
        this.clusterLabel = null;
        this.method = 't-test_overestim_var';
        this.corrMethod = 'benjamini-hochberg';
        this.tableJsonF = null;
        this.tableJsonR = null;
        this.resetUI(); // last to ensure the cluster select values are reset correctly

    }

    /**
     * Resets the UI elements to their default values.
     */
    resetUI() {
        document.querySelector(UI.compareGenesSection).classList.remove("is-hidden");

        document.querySelector(UI.queryClusterSelectElt).value = "";
        document.querySelector(UI.referenceClusterSelectElt).value = "all-reference-clusters";

        document.querySelector(UI.compareGenesNGenesElt).value = 7;
        document.querySelector(UI.compareGenesMethodSelectElt).value = "t-test_overestim_var";
        document.querySelector(UI.comapreGenesCorrMethodSelectElt).value = "benjamini-hochberg";

        // Hide elements that were revealed by previous steps
        document.querySelector(UI.compareGenesRankedContainer).classList.add("is-hidden");

        // show the instructions
        document.querySelector(UI.compareGenesInstructionsElt).classList.remove("is-hidden");
    }

    /**
     * Runs the analysis by sending a POST request to the server and updating the UI with the results.
     *
     * @async
     * @returns {Promise<void>} A promise that resolves when the analysis is completed.
     * @throws {Error} If there is an error computing the comparison.
     */
    async runAnalysis() {
        createToast("Computing comparison", "is-info");

        // Hide things
        document.querySelector(UI.compareGenesSectionSuccessElt).classList.add("is-hidden");
        document.querySelector(UI.compareGenesSectionFailedElt).classList.add("is-hidden");
        document.querySelector(UI.compareGenesResultsContainer).classList.add("is-hidden");

        // Reset the plots
        for (const container of document.querySelectorAll(UI.compareGenesResetableElts)) {
            container.replaceChildren();
        }

        try {
            const {data} = await axios.post("./cgi/h5ad_compare_genes.cgi", convertToFormData({
                'dataset_id': this.analysis.dataset.id,
                'analysis_id': this.analysis.id,
                'analysis_type': this.analysis.type,
                'session_id': this.analysis.userSessionId,
                'n_genes': document.querySelector(UI.compareGenesNGenesElt).value,
                'group_labels': JSON.stringify(this.analysis.groupLabels),
                'query_cluster': document.querySelector(UI.queryClusterSelectElt).value,
                'reference_cluster': document.querySelector(UI.referenceClusterSelectElt).value,
                'method': document.querySelector(UI.compareGenesMethodSelectElt).value,
                'corr_method': document.querySelector(UI.comapreGenesCorrMethodSelectElt).value
            }));

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            document.querySelector(UI.compareGenesInstructionsElt).classList.add("is-hidden");
            document.querySelector(UI.compareGenesResultsContainer).classList.remove("is-hidden");

            this.calculated = true;
            this.nGenes = document.querySelector(UI.compareGenesNGenesElt).value;
            this.queryCluster = document.querySelector(UI.queryClusterSelectElt).value;
            this.referenceCluster = document.querySelector(UI.referenceClusterSelectElt).value;
            this.method = document.querySelector(UI.compareGenesMethodSelectElt).value;
            this.corrMethod = document.querySelector(UI.comapreGenesCorrMethodSelectElt).value;

            this.clusterLabel = data.cluster_label;
            this.updateUIWithResults(data);
            createToast("Comparison computed", "is-success");
        } catch (error) {
            createToast(`Error computing comparison: ${error.message}`);
            logErrorInConsole(error);
            document.querySelector(UI.compareGenesSectionFailedElt).classList.remove("is-hidden");
        } finally {
            if (this.analysis.type !== 'primary') {
                this.analysis.save();
            }
        }
    }

    /**
     * Updates the UI with the analysis results.
     *
     * @param {Object} data - The analysis data.
     * @param {Object} [ana=null] - The analysis object.
     */
    updateUIWithResults(data, ana=null) {

        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        document.querySelector(UI.compareGenesInstructionsElt).classList.add("is-hidden");
        document.querySelector(UI.compareGenesRankedContainer).classList.remove("is-hidden");

        if (this.calculated ) {

            document.querySelector(UI.compareGenesNGenesElt).value = this.nGenes;
            document.querySelector(UI.compareGenesMethodSelectElt).value = this.method;
            document.querySelector(UI.comapreGenesCorrMethodSelectElt).value = this.corrMethod;

            // the query and reference cluster values have to be loaded after the option lists
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': `rank_genes_groups_${data.cluster_label}_comp_ranked`,
            'analysis_type': ana.type,
            'dataset_id': ana.dataset.id,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesRankedContainer});

        params['analysis_name'] = `rank_genes_groups_${data.cluster_label}_${this.queryCluster}_comp_violin`

        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesViolinContainer});

        if (data.hasOwnProperty('table_json_f')) {
            this.tableJsonF = data.table_json_f;
            data.table_json_f = JSON.parse(data.table_json_f);
            this.populateComparisonTable(UI.compareGenesTableFElt, data.table_json_f.data);
        }

        if (this.referenceCluster === 'all-reference-clusters') {
            return;
        }

        // Now do the reverse comparison
        params['analysis_name'] =`rank_genes_groups_${data.cluster_label}_comp_ranked_rev`
        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesRankedRevContainer});

        params["analysis_name"] = `rank_genes_groups_${data.cluster_label}_${this.referenceCluster}_comp_violin_rev`
        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesViolinRevContainer});

        if (data.hasOwnProperty('table_json_r')) {
            this.tableJsonR = data.table_json_r;
            data.table_json_r = JSON.parse(data.table_json_r);
            this.populateComparisonTable(UI.compareGenesTableRElt, data.table_json_r.data);
        }

        // mark success
        document.querySelector(UI.compareGenesSectionSuccessElt).classList.remove("is-hidden");

        passStepWithHref(UI.compareGenesSection);

    }
}
