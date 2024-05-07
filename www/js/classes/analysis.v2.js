"use strict";

/*
  Classes representing overall analysis (pipeline) elements and their child classes.
*/

// requires common.js

let analysisLabels;

class Analysis {
    constructor ({
        id = uuid(),
        datasetId,
        datasetIsRaw = true,
        label = `Unlabeled ${common_datetime()}`,
        type,
        vetting,
        userSessionId = CURRENT_USER.session_id,
        genesOfInterest = [],
        groupLabels = []
    } = {}) {
        this.id = id;
        this.userSessionId = userSessionId;
        this.datasetId = datasetId;
        this.type = type;
        this.vetting = vetting;
        this.label = label;
        this.datasetIsRaw = datasetIsRaw;

        this.primaryFilter = new AnalysisStepPrimaryFilter();
        this.qcByMito = new AnalysisStepQCByMito();
        this.selectVariableGenes = new AnalysisStepSelectVariableGenes();
        this.pca = new AnalysisStepPCA();
        this.tsne = new AnalysisSteptSNE();
        this.clustering = new AnalysisStepClustering(); // The old "louvain" step, which is now done with "leiden"
        this.markerGenes = new AnalysisStepMarkerGenes();
        this.geneComparison = new AnalysisStepCompareGenes();

        this.groupLabels = groupLabels;
        this.genesOfInterest = Array.isArray(genesOfInterest) ? new Set(genesOfInterest) : new Set();

        this.emptyAnalysisOptionTmpl = document.getElementById("analyses-list-empty-tmpl")
        this.analysisOptionTmpl = document.getElementById("analyses-list-tmpl");
        this.analysisSelectElt = document.getElementById("analysis-id");
        this.analysisPrimaryElt = document.getElementById("analyses-primary");
        this.analysisUnsavedElt = document.getElementById("analyses-unsaved");
        this.analysisSavedElt = document.getElementById("analyses-saved");
        this.analysisPublicElt = document.getElementById("analyses-public");
    }

    /**
     * Adds a gene of interest to the set of genes.
     *
     * @param {string} geneSymbol - The symbol of the gene to add.
     */
    addGeneOfInterest(geneSymbol) {
        geneSymbol = geneSymbol.trim();
        this.genesOfInterest.add(geneSymbol);

        if (this.genesOfInterest.length) {
            this.markerGenes.visualizeMarkerGenesBtnElt.classList.remove("is-hidden");
        }
    }

    /**
     * Makes a copy of the analysis to the current user's unsaved area.
     * The most common use case is to copy the public analysis of another user so changes can be made.
     *
     * @param {Function} callback - The callback function to be executed after the copy is made.
     * @param {Object} opts - Additional options for the copy operation.
     */
    async copyToUserUnsaved(callback, opts) {

        const thisAnalysis = this;  // ? Not sure why we are doing this as opposed to just using "this" directly.
        const newAnalysisId = uuid();

        try {
            const {data} = await axios.post("./cgi/copy_dataset_analysis.cgi", {
                session_id: thisAnalysis.userSessionId,
                dataset_id: thisAnalysis.datasetId,
                source_analysis_id: thisAnalysis.id,
                dest_analysis_id: newAnalysisId,
                dest_analysis_type: 'user_unsaved',
                source_analysis_type: thisAnalysis.type
            });

            if (data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.type = 'user_unsaved';
            thisAnalysis.id = newAnalysisId;
            thisAnalysis.userSessionId = CURRENT_USER.session_id;

            // TODO$("#analysis_action_c").show();
            // $("#analysis_status_info_c").hide();

            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, newAnalysisId);

            if (callback) {
                callback(opts);
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
        const thisAnalysis = this;

        try {
            const {data} = await axios.post("./cgi/delete_dataset_analysis.cgi", {
                session_id: thisAnalysis.userSessionId,
                dataset_id: thisAnalysis.datasetId,
                analysis_id: thisAnalysis.id,
                analysis_type: thisAnalysis.type
            });

            if (data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            // Trigger the selection of a 'New' analysis
            // TODO $('#analysis_id option[data-analysis-id="0"]').prop("selected", true).change();
            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, 0);

        } catch (error) {
            createToast(`Error deleting analysis: ${error.message}`);
        }
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
            const {data} = await axios.post("./cgi/get_stored_analysis_list.cgi", {
                dataset_id: datasetId,
                session_id: this.userSessionId
            });

            const emptyAnalysisListHtml = this.emptyAnalysisOptionTmpl.content.cloneNode(true);
            const thisAnalysisLabels = new Set();

            /*
                'primary' analysis is too vague.  That could include clustering and/or
                dimensionality reduction.  The workbench needs clustering to be present,
                while the curator can do fine with just UMAP.  Whether 'primary analysis'
                is added to the menu needs to take this into effect.
            */

            const appendAnalysisOption = (template, parentElt, analysis) => {
                const analysisOptionHtml = template.content.cloneNode(true);
                const option = analysisOptionHtml.querySelector("option");
                option.dataset.analysisId = analysis.id;
                option.dataset.analysisType = analysis.type;
                option.dataset.datasetId = analysis.dataset_id;
                option.textContent = analysis.label;

                // Add to analysis optgroup
                parentElt.appendChild(analysisOptionHtml);
            }

            // primary
            if (data.primary.length) {
                if (forPage == 'sc_workbench') {
                    if (data.primary[0].louvain.calculated) {
                        for (const analysis of data.primary) {
                            thisAnalysisLabels.add(analysis.label);
                            appendAnalysisOption(this.analysisOptionTmpl, this.analysisPrimaryElt, analysis);

                        }
                    } else {
                        this.analysisPrimaryElt.appendChild(emptyAnalysisListHtml);
                    }
                } else {

                    for (const analysis of data.primary) {
                        thisAnalysisLabels.add(analysis.label);
                        appendAnalysisOption(this.analysisOptionTmpl, this.analysisPrimaryElt, analysis);
                    }
                }
            } else {
                this.analysisPrimaryElt.appendChild(emptyAnalysisListHtml);
            }

            // unsaved
            if (data.user_unsaved.length) {
                for (const analysis of data.user_unsaved) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(this.analysisOptionTmpl, this.analysisUnsavedElt, analysis);
                }
            } else {
                this.analysisUnsavedElt.appendChild(emptyAnalysisListHtml);
            }

            // saved
            if (data.user_saved.length) {
                for (const analysis of data.user_saved) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(this.analysisOptionTmpl, this.analysisSavedElt, analysis);
                }
            } else {
                this.analysisSavedElt.appendChild(emptyAnalysisListHtml);
            }

            // public
            if (data.public.length) {
                for (const analysis of data.public) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(this.analysisOptionTmpl, this.analysisPublicElt, analysis);
                }
            } else {
                this.analysisPublicElt.appendChild(emptyAnalysisListHtml);
            }

            // preselect any analysis ID
            document.querySelector(`#analysis-id option[data-analysis-id="${selectedAnalysisId}"]`).setAttribute("selected", "selected");
            // TODO $("#stored_analyses_c").show(10);

            analysisLabels = thisAnalysisLabels;


        } catch (error) {
            createToast(`Error getting saved analyses: ${error.message}`);
            logErrorInConsole(`Failed ID was: ${datasetId} because msg: ${error.message}`);
        }

    }

    /**
     * Loads an Analysis object from JSON data.
     *
     * @param {Object} data - The JSON data representing the Analysis object.
     * @returns {Analysis} The loaded Analysis object.
     */
    static loadFromJson(data) {
        const analysis = new Analysis({
            id: data.id,
            datasetId: data.dataset_id,
            datasetIsRaw: data.dataset_is_raw,
            label: data.label,
            type: data.type,
            userSessionId: data.user_session_id,
            groupLabels: data.group_labels,
            genesOfInterest: data.genesOfInterest
        });

        if (analysis.type == 'primary') {
            // If showing a primary display we only want to show marker genes and gene comparison
            //  tools
            /*
            $("#toggle_marker_genes").bootstrapToggle('enable');
            $("#toggle_marker_genes").bootstrapToggle('on');
            $("#dataset_info").hide();
            */
            return analysis
        }

        analysis.primaryFilter = AnalysisStepPrimaryFilter.loadFromJson(analysis, data.primaryFilter);
        analysis.primaryFilter.updateUI(analysis);

        analysis.qcByMito = AnalysisStepQCByMito.loadFromJson(analysis, data.qcByMito);
        analysis.selectVariableGenes = AnalysisStepSelectVariableGenes.loadFromJson(analysis, data.selectVariableGenes);
        analysis.pca = AnalysisStepPCA.loadFromJson(analysis, data.pca);

        analysis.tsne = AnalysisSteptSNE.loadFromJson(analysis, data.tsne);
        analysis.tsne.updateUI(analysis);

        // Generalize the clustering step instead of being specific to louvain
        // Not doing anything with data.clustering yet but would like to
        if (data.louvain) {
            data.clustering = data.louvain;
        }

        analysis.clustering = AnalysisStepClustering.loadFromJson(analysis, data.louvain);
        analysis.clustering.updateUI(analysis);

        analysis.markerGenes = AnalysisStepMarkerGenes.loadFromJson(analysis, data.markerGenes);
        analysis.geneComparison = AnalysisStepCompareGenes.loadFromJson(analysis, data.geneComparison);

        return analysis;

    }

    /**
     * Makes a public copy of the analysis.
     *
     * @async
     * @function makePublicCopy
     * @memberof analysis.v2
     * @instance
     *
     * @throws {Error} If there is an error making the analysis public.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is successfully made public.
     */
    async makePublicCopy() {
        const thisAnalysis = this;

        try {
            const {data} = await axios.post("./cgi/copy_dataset_analysis.cgi", {
                session_id: thisAnalysis.userSessionId,
                dataset_id: thisAnalysis.datasetId,
                source_analysis_id: thisAnalysis.id,
                dest_analysis_id: thisAnalysis.id,
                source_analysis_type: thisAnalysis.type,
                dest_analysis_type: 'public'
            });

            if (data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.type = 'user_saved';
            // TODO $("#analysis_action_c").hide();
            // $("#analysis_status_info").text("Changes made to this public analysis will spawn a local copy within your profile.");
            // $("#analysis_status_info_c").show();
            // $("#btn_make_public_copy").hide();
            // $("#btn_delete_saved_analysis").hide();
            // $("#btn_delete_unsaved_analysis").hide();

            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, thisAnalysis.id);
        } catch (error) {
            createToast(`Error making analysis public: ${error.message}`);
        }

    }

    /**
     * Grabs an image from the server, usually generated temporarily by a module like scanpy,
     * and places it into the target location as a binary stream.
     * This prevents us from having to keep a lot of temporary images on the server
     * and also ensures the user never gets a cached image.
     *
     * @async
     * @param {Object} options - The options for placing the analysis image.
     * @param {Object} options.params - The parameters for the image request.
     * @param {string} options.title - The title of the image.
     * @param {HTMLElement} [options.target] - The target element where the image will be placed.
     * @returns {Promise<void>} A promise that resolves when the image is placed.
     */
    async placeAnalysisImage ({params, title, target = []} = {}) {
        const imgSrc = await axios.get("./cgi/get_analysis_image.cgi", {
            params: params
        });
        const html = `<a target="_blank" href="${imgSrc}"><img src="${imgSrc}" class="image" alt="${title}" /></a>`;
        target.appendChild(html);
    }

    /**
     * Removes a gene of interest from the set of genes.
     *
     * @param {string} geneSymbol - The symbol of the gene to be removed.
     */
    removeGeneOfInterest(geneSymbol) {
        this.genesOfInterest.delete(geneSymbol);

        if (!this.genesOfInterest.length) {
            this.markerGenes.visualizeMarkerGenesBtnElt.classList.add("is-hidden");
        }
    }


    /**
     * Resets the Analysis instance to its initial state.
     * This method starts over as a completely new Analysis instance, with a new ID,
     * and resets all existing components.
     */
    reset() {
        this.clustering.reset();
        this.markerGenes.reset();
        this.pca.reset();
        this.tsne.reset();
        this.primaryFilter.reset();
        this.qcByMito.reset();
        this.selectVariableGenes.reset();
        this.geneComparison.reset();
        this.datasetIsRaw = true;
        this.id = uuid();
        this.label = null;
    }

    async save() {
        /*
          Saves the current analysis parameters to disk.
         */
        const state = JSON.stringify(this);
        const thisAnalysis = this;

        try {
            const {data} = await axios.post("./cgi/save_dataset_analysis.cgi", {
                session_id: thisAnalysis.userSessionId,
                dataset_id: thisAnalysis.datasetId,
                analysis_id: thisAnalysis.id,
                analysis_type: thisAnalysis.type,
                analysis_vetting: thisAnalysis.vetting,
                label: thisAnalysis.label,
                state: state
            });
            if (data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, thisAnalysis.id);
        } catch (error) {
            createToast(`Error saving analysis: ${error.message}`);
        }
    }

    async saveToUserArea() {

        const thisAnalysis = this;

        try {
            const {data} = await axios.post("./cgi/copy_dataset_analysis.cgi", {
                session_id: thisAnalysis.userSessionId,
                dataset_id: thisAnalysis.datasetId,
                source_analysis_id: thisAnalysis.id,
                dest_analysis_id: thisAnalysis.id,
                source_analysis_type: thisAnalysis.type,
                dest_analysis_type: 'user_saved'
            });
            if (data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.type = 'user_saved';
            // TODO $("#btn_save_analysis").text("Saved");
            // $("#analysis_action_c").hide();
            // $("#analysis_status_info").text("This analysis is stored in your profile.");
            // $("#analysis_status_info_c").show(500);
            // $("#btn_delete_saved_analysis").show();
            // $("#btn_make_public_copy").show();
            // $("#new_analysis_label_c").hide();

            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, thisAnalysis.id);

        } catch (error) {
            createToast(`Error saving analysis: ${error.message}`);
        }
    }

}

/* Putting these in order of the workbench steps */

class AnalysisStepPrimaryFilter {
    constructor({}) {
        this.reset();

        this.primaryFilterSectionElt = document.getElementById("primary-filter-s");
        this.primaryFilterBtnElt = document.getElementById("btn-primary-filter");
    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisStepQCByMito {
    constructor({}) {
        this.reset();

        this.qcByMitoSectionElt = document.getElementById("qc-by-mito-s");
        this.qcByMitoBtnElt = document.getElementById("btn-qc-by-mito");
    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisStepSelectVariableGenes {
    constructor({}) {
        this.reset();

        this.selectVariableGenesSectionElt = document.getElementById("select-variable-genes-s");
        this.selectVariableGenesBtnElt = document.getElementById("btn-select-variable-genes");
    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisStepPCA {
    constructor({}) {
        this.reset();

        this.pcaSectionElt = document.getElementById("pca-s");
        this.pcaBtnElt = document.getElementById("btn-pca");
    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisSteptSNE {
    constructor({}) {
        this.reset();

        this.tsneSectionElt = document.getElementById("tsne-s");
        this.tsneBtnElt = document.getElementById("btn-tsne");
    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisStepClustering {
    constructor({}) {
        this.reset();   // TODO: split into two potentially for each step

        this.clusteringSectionElt = document.getElementById("clustering-s");
        this.clusteringBtnElt = document.getElementById("btn-clustering");

        // Rename, Merge, Delete clusters from first clustering step
        this.editClusteringSectionElt = document.getElementById("edit-clustering-s");
        this.editClusteringBtnElt = document.getElementById("btn-edit-clustering");

    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisStepMarkerGenes {
    constructor({}) {
        this.reset();

        this.visualMarkerGenesSectionElt = document.getElementById("visualize-marker-genes-s");
        this.visualizeMarkerGenesBtnElt = document.getElementById("btn-visualize-marker-genes");
    }

    static loadFromJason(ana, data) {

    }
}
class AnalysisStepCompareGenes {
    constructor({}) {
        this.reset();

        this.compareGenesSectionElt = document.getElementById("gene-comparison-s");
        this.compareGenesBtnElt = document.getElementById("btn-gene-comparison");
        this.clusterOptsTmpl = document.getElementById("cluster-list-tmpl");
        this.queryClusterOptionsElt = document.getElementById("query-cluster-options");
        this.referenceClusterOptionsElt = document.getElementById("reference-cluster-options");
        this.queryClusterSelectElt = document.getElementById("query-cluster");
        this.referenceClusterSelectElt = document.getElementById("reference-cluster");

        this.compareGenesNGenesElt = document.getElementById("compare-genes-n-genes");
        this.compareGenesMethodElt = document.getElementById("compare-genes-method");
        this.comapreGenesCorrMethodElt = document.getElementById("compare-genes-corr-method");

    }

    /**
     * Loads an instance of AnalysisStepCompareGenes from JSON data.
     *
     * @param {Analysis} ana - The Analysis instance.
     * @param {Object} data - The JSON data to load from.
     * @returns {AnalysisStepCompareGenes} - The loaded AnalysisStepCompareGenes instance.
     */
    static loadFromJason(ana, data) {
        const step = new AnalysisStepCompareGenes();
        if (data === undefined) {
            data = {
                "calculated": false,
                "n_genes": 0,
                "query_cluster": null,
                "reference_cluster": null,
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
        step.method = 't-test_overestim_var';
        step.corrMethod = 'benjamini-hochberg';

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
            step.updateUiWithResults(ana, data);
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
        const tableElt = document.querySelector(`#${tableId} tbody`);

        for (const row of tableJson) {
            const rowElt = document.createElement("tr");
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

        let groupNum = 0;

        for (const label of groupLabels) {

            const clusterOptionHtml = this.clusterOptsTmpl.content.cloneNode(true);
            const option = clusterOptionHtml.querySelector("option");
            option.textContent = label;
            option.value = groupNum;

            this.queryClusterOptionsElt.appendChild(clusterOptionHtml);
            this.referenceClusterOptionsElt.appendChild(clusterOptionHtml);

            groupNum++;
        }

        this.queryClusterSelectElt.value = this.queryCluster;
        this.referenceClusterSelectElt.value = this.referenceCluster;
    }

    /**
     * Resets the analysis object to its initial state.
     */
    reset() {
        this.calculated = false;
        this.resetUi();
        this.nGenes = 0;
        this.queryCluster = null;
        this.referenceCluster = null;
        this.method = 't-test_overestim_var';
        this.corrMethod = 'benjamini-hochberg';
        this.tableJsonF = null;
        this.tableJsonR = null;
    }

    /**
     * Resets the UI elements to their default values.
     */
    resetUi() {

        this.queryClusterSelectElt.value = null;
        this.referenceClusterSelectElt.value = "all-reference-clusters";

        this.compareGenesNGenesElt.value = null;
        this.compareGenesMethodElt.value = "t-test_overestim_var";
        this.comapreGenesCorrMethodElt.value = "benjamini-hochberg";

        // TODO $('#analysis_compare_genes .tool_instructions').show();
        // $('#compare_genes_results_c').hide();
    }


    /**
     * Updates the UI with the analysis results.
     *
     * @param {Object} ana - The analysis object.
     * @param {Object} data - The data object containing the analysis results.
     */
    updateUiWithResults(ana, data) {

        // TODO $('#analysis_compare_genes .tool_instructions').hide();
        // $('#compare_genes_results_c').show(500);

        if (!this.calculated ) {
            return;
        }

        // TODO $('#toggle_compare_genes').bootstrapToggle('on');

        this.compareGenesNGenesElt.value = this.nGenes;
        this.compareGenesMethodElt.value = this.method;
        this.comapreGenesCorrMethodElt.value = this.corrMethod

        // the query and reference cluster values have to be loaded after the option lists

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'rank_genes_groups_' + data['cluster_label'] + '_comp_ranked',
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            datetime: (new Date()).getTime()
        }

        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_ranked_c'});

        params['analysis_name'] = `rank_genes_groups_${data['cluster_label']}_${this.queryCluster}_comp_violin`

        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_violin_c'});

        if (data.hasOwnProperty('table_json_f')) {
            this.tableJsonF = data['table_json_f'];
            data['table_json_f'] = JSON.parse(data['table_json_f']);
            this.populateComparisonTable(ana, 'compare_genes_table_f', data['table_json_f']['data']);
        }

        if (data.hasOwnProperty('table_json_r')) {
            this.tableJsonR = data['table_json_r'];
            data['table_json_r'] = JSON.parse(data['table_json_r']);
            this.populateComparisonTable(ana, 'compare_genes_table_r', data['table_json_r']['data']);
        }

        if (this.referenceCluster != 'all-reference-clusters') {
            params['analysis_name'] ='rank_genes_groups_' + data['cluster_label'] + '_comp_ranked_rev'
            ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_ranked_rev_c'});

            params["analysis_name"] = `rank_genes_groups_${data['cluster_label']}_${this.referenceCluster}_comp_violin_rev`
            ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': '#compare_genes_violin_rev_c'});
        }


    }
}
