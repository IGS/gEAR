"use strict";

/*
    Classes representing overall analysis (pipeline) elements and their child classes.
*/

// requires common.js

let analysisLabels;

class UI {
    // This class is a singleton that manages the UI elements of the analysis pipeline.

    // general analysis
    emptyAnalysisOptionTmpl = "#analyses-list-empty-tmpl"
    analysisOptionTmpl = "#analyses-list-tmpl"
    analysisSelectElt = "#analysis-id"
    newAnalysisOptionElt = `${analysisSelectElt} option[data-analysis-id='0']`
    newAnalysisLabelContainer = "#new-analysis-label-c"
    analysisPrimaryElt = "#analyses-primary"
    analysisUnsavedElt = "#analyses-unsaved"
    analysisSavedElt = "#analyses-saved"
    analysisPublicElt = "#analyses-public"
    analysisActionContainer = "#analysis-action-c"
    analysisStatusInfoContainer = "#analysis-status-info-c"
    analysisStatusInfoElt = "#analysis-status-info"
    storedAnalysesContainer = "#stored-analyses-c"
    datasetInfoElt = "#dataset-info"
    btnSaveAnalysisElt = "#btn-save-analysis"
    btnDeleteSavedAnalysis = "#btn-delete-saved-analysis"
    btnMakePublicCopy = "#btn-make-public-copy"

    // primary filter
    primaryFilterSectionElt = "#primary-filter-s"
    primaryFilterBtnElt = "#btn-primary-filter"
    filterCellsGtNGenesElt = "#filter-cells-gt-n-genes"
    filterCellsLtNGenesElt = "#filter-cells-lt-n-genes"
    filterGenesGtNCellsElt = "#filter-genes-gt-n-cells"
    filterGenesLtNCellsElt = "#filter-genes-lt-n-cells"
    filterCellsGtNGenesSelectedElt = "#filter-cells-gt-n-genes-selected"
    filterCellsLtNGenesSelectedElt = "#filter-cells-lt-n-genes-selected"
    filterGenesGtNCellsSelectedElt = "#filter-genes-gt-n-cells-selected"
    filterGenesLtNCellsSelectedElt = "#filter-genes-lt-n-cells-selected"
    selectedDatasetShapeFilteredElt = "#selected-dataset-shape-filtered"
    selectedDatasetShapeFilteredContainerElt = "#selected-dataset-shape-filtered-c"
    primaryInitialPlotElt = ".primary-initial-plot"
    primaryInitialPlotContainerElt = "#primary-initial-plot-c"
    primaryTopGenesContainer = "#primary-top-genes-c"
    primaryTopGenesPlotContainer = "#primary-top-genes-plot-c"
    primaryInitialPlotElts = ".primary-initial-plot"

    // QC by mito
    qcByMitoToggleElt = "#toggle-qc-by-mito"    // Temporary
    qcByMitoSectionElt = "#qc-by-mito-s"
    qcByMitoBtnElt = "#btn-qc-by-mito"
    qbmGenePrefixElt = "#qbm-gene-prefix"
    qbmFilterMitoPercElt = "#qbm-filter-mito-perc"
    qbmFilterMitoCountElt = "#qbm-filter-mito-count"
    btnQbmSaveElt = "#btn-qbm-save"
    btnDoAnalysisQcByMitoElt = "#btn-do-analysis-qc-by-mito"
    qbmPostShapeElt = "#qbm-post-shape"
    qbmGeneCountElt = "#qbm-gene-count"
    qbmObsCountElt = "#qbm-obs-count"
    qbmInstructionsElt = "#analysis-qc-by-mito .tool-instructions"
    qbmViolinContainer = "#qbm-violin-c"
    qbmScatterPercentMitoContainer = "#qbm-scatter-percent-mito-c"
    qbmScatterNGenesContainer = "#qbm-scatter-n-genes-c"

    // select variable genes
    selectVariableGenesToggleElt = "#toggle-select-variable-genes"  // Temporary
    selectVariableGenesSectionElt = "#select-variable-genes-s"
    selectVariableGenesBtnElt = "#btn-select-variable-genes"
    asvgNormCountsPerCellElt = "#asvg-norm-counts-per-cell"
    asvgFlavorElt = "#asvg-flavor"
    asvgNTopGenesElt = "#asvg-n-top-genes"
    asvgMinMeanElt = "#asvg-min-mean"
    asvgMaxMeanElt = "#asvg-max-mean"
    asvgMinDispersionElt = "#asvg-min-dispersion"
    asvgRegressOutElt = "#asvg-regress-out"
    asvgScaleUnitVarianceElt = "#asvg-scale-unit-variance"
    btnAsvgSaveElt = "#btn-asvg-save"
    btnDoAnalysisSelectVariableGenesElt = "#btn-do-analysis-select-variable-genes"
    asvgSaveOptionsElts = ".asvg-save-options"
    asvgPlotContainer = "#asvg-plot-c"
    asvgPlotNormContainer = "#asvg-plot-norm-c"
    asvgInstructionsElt = "#analysis-select-variable-genes .tool-instructions"

    // PCA
    pcaToggleElt = "#toggle-pca"   // Temporary
    pcaSectionElt = "#pca-s"
    btnPcaRunElt = "#btn-pca-run"
    genesToColorElt = "#pca-genes-to-color"
    pcaOptionsGroupElt = "#pca-options-group"
    topPcaGenesElt = "#pca-top-genes"
    pcaInstructionsElt = "#analysis-pca .tool-instructions"
    pcaScatterContainer = "#pca-scatter-c"
    pcaVarianceContainer = "#pca-variance-c"
    pcaMissingGeneContainer = "#pca-missing-gene-c"
    pcaMissingGeneElt = "#pca-missing-gene"
    pcaImageResultContainers = "#analysis-pca div.image-result-c"
    pcaOptionsGroupElt = "#pca-options-g"
    weightedGeneCartGroupElt = "#weighted-gene-cart-g"


    // tSNE
    tsneToggleElt = "#toggle-tsne"   // Temporary
    tsneSectionElt = "#tsne-s"
    btnTsneRunElt = "#btn-tsne-run"
    tsneGenesToColorElt = "#tsne-genes-to-color"
    dimReductionNNeighborsElt = "#dredux-n-neighbors"
    tsneNPcsElt = "#tsne-n-pcs"
    tsneRandomStateElt = "#tsne-random-state"
    dimensionalityReductionMethodTsneElt = "#dimensionality-reduction-method-tsne"
    dimensionalityReductionMethodUmapElt = "#dimensionality-reduction-method-umap"
    tsneInstructionsElt = "#analysis-tsne .tool-instructions"
    tsnePlotContainer = "#tsne-plot-c"
    umapPlotContainer = "#umap-plot-c"
    tsneMissingGeneContainer = "#tsne-missing-gene-c"
    tsneUseScaledElt = "#tsne-use-scaled"

    // Clustering
    clusteringToggleElt = "#toggle-clustering"   // Temporary
    clusteringSectionElt = "#clustering-s"
    btnClusteringRunElt = "#btn-clustering-run"
    btnClusteringRerunWithGroupsElt = "#btn-clustering-rerun-with-groups"
    resolutionElt = "#clustering-resolution"
    clusteringNNeighborsElt = "#clustering-n-neighbors"
    clusterTsnePlotElt = "#cluster-tsne-plot-c"
    clusterUmapPlotElt = "#cluster-umap-plot-c"
    clusteringInstructionsElt = "#analysis-clustering .tool-instructions"
    clusteringImageResultContainers = "#analysis-clustering div.image-result-c"

    // Marker Genes
    markerGenesToggleElt = "#toggle_marker_genes"   // Temporary
    btnMarkerGenesRunElt = "#btn-marker-genes-run"
    visualMarkerGenesSectionElt = "#visualize-marker-genes-s"
    visualizeMarkerGenesBtnElt = "#btn-visualize-marker-genes"
    downloadMarkerGenesBtnElt = "#btn-download-marker-genes"
    markerGenesNGenesElt = "#marker-genes-n-genes"
    markerGenesTableElt = "#marker-genes-table"
    markerGenesTableHeadTmpl = "#marker-genes-table-head-tmpl"
    markerGenesTableHeadRowElt = `${markerGenesTableElt} thead tr`
    markerGenesTableBodyTmpl = "#marker-genes-table-body-tmpl"
    markerGenesTableBodyElt = `${markerGenesTableElt} tbody`
    markerGenesTableHeader = "#marker-genes-table-header"
    markerGenesPlotContainer = "#marker-genes-plot-c"
    markerGenesVisualizationContainer = "#marker-genes-visualization-c"
    markerGenesDotplotContainer = "#marker-genes-dotplot-c"
    markerGenesViolinContainer = "#marker-genes-violin-c"
    markerGenesManuallyEnteredElt = "#marker-genes-manually-entered"
    markerGenesSelectedCountElt = "#marker-genes-selected-count"
    markerGenesEnteredCountElt = "#marker-genes-entered-count"
    markerGenesUniqueCountElt = "#marker-genes-unique-count"

    groupLabelsContainer = "#group-labels-c"
    clusterGroupLabelsTmpl = "#cluster-group-labels-tmpl"
    clusterGroupLabelsTableBodyElt = "#cluster-group-labels tbody"
    clusterGroupLabelsInputElts = '#cluster-group-labels td.group-user-label input'

    // Clustering (edit mode)
    clusteringToggleElt = "#toggle-clustering-edit"   // Temporary
    editClusteringSectionElt = "#edit-clustering-s"
    editClusteringBtnElt = "#btn-edit-clustering"
    // -- using resolutionElts values from the non-edit clustering
    // -- using clusterNNeighborsElt values from the non-edit clustering
    clusterTsnePlotEditElt = "#cluster-tsne-plot-edit-c"
    clusterUmapPlotEditElt = "#cluster-umap-plot-edit-c"
    clusteringEditInstructionsElt = "#analysis-clustering-edit .tool-instructions"
    // Compare Genes
    compareGenesToggleElt = "#toggle-compare-genes"   // Temporary
    geneComparisonSectionElt = "#gene-comparison-s"
    geneComparisonBtnElt = "#btn-gene-comparison"
    clusterOptsTmpl = "#cluster-list-tmpl"
    queryClusterOptionsElt = "#query-cluster-options"
    referenceClusterOptionsElt = "#reference-cluster-options"
    queryClusterSelectElt = "#query-cluster"
    referenceClusterSelectElt = "#reference-cluster"
    compareGenesNGenesElt = "#compare-genes-n-genes"
    compareGenesMethodElt = "#compare-genes-method"
    comapreGenesCorrMethodElt = "#compare-genes-corr-method"
    compareGenesRankedContainer = "#compare-genes-ranked-c"
    compareGenesViolinContainer = "#compare-genes-violin-c"
    compareGenesRankedRevContainer = "#compare-genes-ranked-rev-c"
    compareGenesViolinRevContainer = "#compare-genes-violin-rev-c"
    compareGenesInstructionsElt = "#analysis-compare-genes .tool-instructions"
    compareGenesResultsContainer = "#compare-genes-results-c"
    compareGenesImagesResultsContainers = `${compareGenesResultsContainer} div.image-result-c`

}

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

        this.primaryFilter = new AnalysisStepPrimaryFilter(this);
        this.qcByMito = new AnalysisStepQCByMito(this);
        this.selectVariableGenes = new AnalysisStepSelectVariableGenes(this);
        this.pca = new AnalysisStepPCA(this);
        this.tsne = new AnalysisSteptSNE(this);
        this.clustering = new AnalysisStepClustering(this); // The old "louvain" step, which is now done with "leiden"
        this.markerGenes = new AnalysisStepMarkerGenes(this);
        this.clusteringEdit = new AnalysisStepClusteringEdit(this, "edit");
        this.compareGenes = new AnalysisStepCompareGenes(this);

        this.groupLabels = groupLabels;
        this.genesOfInterest = Array.isArray(genesOfInterest) ? new Set(genesOfInterest) : new Set();
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

            if (data?.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.type = 'user_unsaved';
            thisAnalysis.id = newAnalysisId;
            thisAnalysis.userSessionId = CURRENT_USER.session_id;

            document.querySelector(UI.analysisActionContainer).classList.remove("is-hidden");
            document.querySelector(UI.analysisStatusInfoContainer).classList.add("is-hidden");

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

            if (data?.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            // Trigger the selection of a 'New' analysis
            document.querySelector(UI.newAnalysisOptionElt).setAttribute("selected", "selected");
            document.querySelector(UI.analysisSelectElt).dispatchEvent(new Event("change"));
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

            const emptyAnalysisListHtml = document.querySelector(UI.emptyAnalysisOptionTmpl).content.cloneNode(true);
            const thisAnalysisLabels = new Set();

            /*
                'primary' analysis is too vague.  That could include clustering and/or
                dimensionality reduction.  The workbench needs clustering to be present,
                while the curator can do fine with just UMAP.  Whether 'primary analysis'
                is added to the menu needs to take this into effect.
            */

            const appendAnalysisOption = (template, parentSelector, analysis) => {
                const analysisOptionHtml = document.querySelector(template).content.cloneNode(true);
                const option = analysisOptionHtml.querySelector("option");
                option.dataset.analysisId = analysis.id;
                option.dataset.analysisType = analysis.type;
                option.dataset.datasetId = analysis.dataset_id;
                option.textContent = analysis.label;

                // Add to analysis optgroup
                document.querySelector(parentSelector).appendChild(analysisOptionHtml);
            }

            // primary
            if (data.primary.length) {
                if (forPage == 'sc_workbench') {
                    if (data.primary[0].louvain.calculated) {
                        for (const analysis of data.primary) {
                            thisAnalysisLabels.add(analysis.label);
                            appendAnalysisOption(UI.analysisOptionTmpl, UI.analysisPrimaryElt, analysis);

                        }
                    } else {
                        document.querySelector(UI.analysisPrimaryElt).appendChild(emptyAnalysisListHtml);
                    }
                } else {

                    for (const analysis of data.primary) {
                        thisAnalysisLabels.add(analysis.label);
                        appendAnalysisOption(UI.analysisOptionTmpl, UI.analysisPrimaryElt, analysis);
                    }
                }
            } else {
                document.querySelector(UI.analysisPrimaryElt).appendChild(emptyAnalysisListHtml);
            }

            // unsaved
            if (data.user_unsaved.length) {
                for (const analysis of data.user_unsaved) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(UI.analysisOptionTmpl, UI.analysisUnsavedElt, analysis);
                }
            } else {
                document.querySelector(UI.analysisUnsavedElt).appendChild(emptyAnalysisListHtml);
            }

            // saved
            if (data.user_saved.length) {
                for (const analysis of data.user_saved) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(UI.analysisOptionTmpl, UI.analysisSavedElt, analysis);
                }
            } else {
                document.querySelector(UI.analysisSavedElt).appendChild(emptyAnalysisListHtml);
            }

            // public
            if (data.public.length) {
                for (const analysis of data.public) {
                    thisAnalysisLabels.add(analysis.label);
                    appendAnalysisOption(UI.analysisOptionTmpl, UI.analysisPublicElt, analysis);
                }
            } else {
                document.querySelector(UI.analysisPublicElt).appendChild(emptyAnalysisListHtml);
            }

            // preselect any analysis ID
            document.querySelector(`#analysis-id option[data-analysis-id="${selectedAnalysisId}"]`).setAttribute("selected", "selected");
            document.querySelector(UI.storedAnalysesContainer).classList.remove("is-hidden");

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
            document.querySelector(UI.markerGenesToggleElt).classList.remove("is-hidden");
            document.querySelector(UI.datasetInfoElt).classList.add("is-hidden");
            return analysis
        }

        analysis.primaryFilter = AnalysisStepPrimaryFilter.loadFromJson(data.primaryFilter, analysis);
        analysis.primaryFilter.updateUIWithResults(analysis);

        analysis.qcByMito = AnalysisStepQCByMito.loadFromJson(data.qcByMito, analysis);
        analysis.selectVariableGenes = AnalysisStepSelectVariableGenes.loadFromJson(data.selectVariableGenes, analysis);
        analysis.pca = AnalysisStepPCA.loadFromJson(data.pca, analysis);

        analysis.tsne = AnalysisSteptSNE.loadFromJson(data.tsne, analysis);
        analysis.tsne.updateUIWithResults(analysis);

        // Generalize the clustering step instead of being specific to louvain
        // Not doing anything with data.clustering yet but would like to
        if (data.louvain) {
            data.clustering = data.louvain;
        }

        analysis.clustering = AnalysisStepClustering.loadFromJson(data.clustering, analysis);
        analysis.clustering.updateUIWithResults(analysis);

        analysis.markerGenes = AnalysisStepMarkerGenes.loadFromJson(data.markerGenes, analysis);

        analysis.clusteringEdit = AnalysisStepClusteringEdit.loadFromJson(data.clustering, analysis);
        analysis.clusteringEdit.mode = "edit";

        analysis.compareGenes = AnalysisStepCompareGenes.loadFromJson(data.compareGenes, analysis);

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

            if (data?.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.type = 'user_saved';
            document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
            document.querySelector(UI.analysisStatusInfoElt).textContent = "Changes made to this public analysis will spawn a local copy within your profile.";
            document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
            document.querySelector(UI.btnMakePublicCopy).classList.add("is-hidden");
            document.querySelector(UI.btnDeleteSavedAnalysis).classList.add("is-hidden");
            document.querySelector(UI.btnDeleteUnsavedAnalysis).classList.add("is-hidden");

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
        this.compareGenes.reset();
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
                state
            });
            if ((!data.success) || (data.success < 1)) {
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
            if ((!data.success) || (data.success < 1)) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            thisAnalysis.type = 'user_saved';
            document.querySelector(UI.btnSaveAnalysisElt).textContent = "Saved";
            document.querySelector(UI.analysisActionContainer).classList.add("is-hidden");
            document.querySelector(UI.analysisStatusInfoElt).textContent = "This analysis is stored in your profile.";
            document.querySelector(UI.analysisStatusInfoContainer).classList.remove("is-hidden");
            document.querySelector(UI.btnDeleteSavedAnalysis).classList.remove("is-hidden");
            document.querySelector(UI.btnMakePublicCopy).classList.remove("is-hidden");
            document.querySelector(UI.newAnalysisLabelContainer).classList.add("is-hidden");

            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, thisAnalysis.id);

        } catch (error) {
            createToast(`Error saving analysis: ${error.message}`);
        }
    }

}

/* Putting these in order of the workbench steps */

class AnalysisStepPrimaryFilter {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
        this.irreversible = true;
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

        step.filteredGeneCound = data['filtered_gene_count'];
        step.filteredCellCount = data['filtered_cell_count'];

        return step;
    }

    /**
     * Resets the state of the analysis object.
     * Sets all properties to their initial values and calls the `resetUI` method.
     */
    reset() {
        this.calculated = false;
        this.filterCellsGtNGenes = null;
        this.filterCellsGtNGenesSelected = null;
        this.filterCellsLtNGenes = null;
        this.filterCellsLtNGenesSelected = null;

        this.filterGenesGtNCells = null;
        this.filterGenesGtNCellsSelected = null;
        this.filterGenesLtNCells = null;
        this.filterGenesLtNCellsSelected = null;

        this.filteredGeneCound = null;
        this.filteredCellCount = null;
        this.resetUI();
    }

    /**
     * Resets the UI by setting the values of various filter inputs and checkboxes,
     * and modifying the visibility of certain elements.
     */
    resetUI() {
        document.querySelector(UI.filterCellsGtNGenesElt).value = 300;
        document.querySelector(UI.filterCellsLtNGenesElt).value = '';
        document.querySelector(UI.filterGenesLtNCellsElt).value = 3;
        document.querySelector(UI.filterGenesGtNCellsElt).value = '';

        document.querySelector(UI.filterCellsLtNGenesSelectedElt).checked = false;
        document.querySelector(UI.filterCellsGtNGenesSelectedElt).checked = false;
        document.querySelector(UI.filterGenesLtNCellsSelectedElt).checked = false;
        document.querySelector(UI.filterGenesGtNCellsSelectedElt).checked = false;

        document.querySelector(UI.primaryInitialPlotElt).classList.remove("is-hidden");
        document.querySelector(UI.primaryInitialPlotContainerElt).classList.add("is-hidden");
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

            document.querySelector(UI.filterCellsLtNGenesElt).value = this.filterCellsLtNGenes;
            document.querySelector(UI.filterCellsGtNGenesElt).value = this.filterCellsGtNGenes;
            document.querySelector(UI.filterGenesLtNCellsElt).value = this.filterGenesLtNCells;
            document.querySelector(UI.filterGenesGtNCellsElt).value = this.filterGenesGtNCells;

            document.querySelector(UI.selectedDatasetShapeFilteredElt).textContent = `${this.filteredGeneCound} genes x ${this.filteredCellCount} obs`;
            document.querySelector(UI.selectedDatasetShapeFilteredContainerElt).classList.remove("is-hidden");

        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'highest_expr_genes',
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            datetime: (new Date()).getTime()
        }

        ana.placeAnalysisImage(
            {'params': params, 'title': 'Highest expressed genes', 'target': UI.primaryTopGenesContainer});

        for (const elt of document.querySelectorAll(UI.primaryInitialPlotElts)) {
            elt.classList.add("is-hidden");
        }
        document.querySelector(UI.primaryInitialPlotContainerElt).classList.remove("is-hidden");

        document.querySelector(UI.qcByMitoToggleElt).classList.remove("is-hidden");
        document.querySelector(UI.selectVariableGenesToggleElt).classList.remove("is-hidden");

    }
}
class AnalysisStepQCByMito {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
        this.irreversible = false;   // If "saved", then this is true as filtering happens in the backend
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

        if (step.calculated == true) {
            document.querySelector(UI.qcByMitoToggleElt).checked = true;
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
        this.genePrefix = null;
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
        document.querySelector(UI.qbmGenePrefixElt).value = 'mt-';
        disableAndHideElement(document.querySelector(UI.btnQbmSaveElt));
        document.querySelector(UI.btnQbmSaveElt).textContent = 'Save these genes';

        document.querySelector(UI.btnDoAnalysisQcByMitoElt).classList.remove("is-hidden");
    }

    /**
     * Runs the analysis.
     *
     * @param {boolean} saveDataset - Indicates whether to save the dataset.
     * @returns {Promise<void>} - A promise that resolves when the analysis is complete.
     * @throws {Error} - If there is an error during the analysis.
     */
    async runAnalysis(saveDataset) {
        disableAndHideElement(document.querySelector(UI.btnQbmSaveElt));
        if (Boolean(saveDataset)) {
            createToast("Applying mitochondrial filter", "is-info");
            document.querySelector(UI.btnQbmSaveElt).textContent = 'Saving';
        } else {
            document.querySelector(UI.btnQbmSaveElt).classlist.remove("is-hidden");
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
            const {data} = await axios.post("./cgi/h5ad_qc_by_mito.cgi", {
                dataset_id: this.analysis.datasetId,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                genes_prefix: this.genePrefix,
                filter_mito_perc: this.filterMitoPercent,
                filter_mito_count: this.filterMitoCount,
                save_dataset: saveDataset
            });

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.nGenes = data.n_genes;
            this.nObs = data.n_obs;

            document.querySelector(UI.btnQbmSaveElt).disabled = false;
            this.calculated = false;
            if (Boolean(saveDataset)) {
                document.querySelector(UI.btnQbmSaveElt).disabled = true;
                this.calculated = true;
            }

            this.updateUIWithResults(this.calculated, data);
            createToast("Mitochondrial plot displayed", "is-success");
            // TODO - Show next step collapsable

        } catch (error) {
            createToast(`Error doing QC analysis: ${error.message}`);
        } finally {
            document.querySelector(UI.btnDoAnalysisQcByMitoElt).disabled = false;
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

        enableAndShowElement(UI.btnQbmSaveElt);

        if (resultsSaved) {
            document.querySelector(UI.btnDoAnalysisQcByMitoElt).classList.add("is-hidden");
            document.querySelector(UI.btnQbmSaveElt).disabled = true;
            document.querySelector(UI.btnQbmSaveElt).textContent = 'Saved';
            document.querySelector(UI.qbmPostShapeElt).classList.remove("is-hidden");
        }

        document.querySelector(UI.qbmInstructionsElt).classList.add("is-hidden");
        document.querySelector(UI.qbmGeneCountElt).textContent = this.nGenes;
        document.querySelector(UI.qbmObsCountElt).textContent = this.nObs;

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'violin_qc_by_mito',
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
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

    }
}
class AnalysisStepSelectVariableGenes {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
        this.irreversible = false;  // TODO: Add "adata.layers["counts"] = adata.X.copy()" to the backend
        // At this step, we theoretically could save original counts to a layer before normalizing
        // Then if we want to back up to this step, we could have a script to restore the original counts
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
        this.normCountsPerCell = null;
        this.flavor = 'seurat';
        this.nTopGenes = null;
        this.minMean = null;
        this.maxMean = null;
        this.minDispersion = null;
        this.regressOut = true;
        this.scaleUnitVariance = true;
        this.resetUI();
    }

    /**
     * Resets the UI elements to their default values.
     */
    resetUI() {
        document.querySelector(UI.asvgNormCountsPerCellElt).value = '1e4';
        document.querySelector(UI.asvgFlavorElt).value = 'seurat';
        document.querySelector(UI.asvgNTopGenesElt).value = '';
        document.querySelector(UI.asvgMinMeanElt).value = 0.0125;
        document.querySelector(UI.asvgMaxMeanElt).value = 3;
        document.querySelector(UI.asvgMinDispersionElt).value = 0.5;
        document.querySelector(UI.asvgRegressOutElt).checked = true;
        document.querySelector(UI.asvgScaleUnitVarianceElt).checked = true;

        disableAndHideElement(document.querySelector(UI.btnAsvgSaveElt));
        document.querySelector(UI.btnAsvgSaveElt).textContent = 'Save these genes';
        document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).classList.remove("is-hidden");
        for (const elt of document.querySelectorAll(UI.asvgSaveOptionsElts)) {
            elt.classList.add("is-hidden");
        }
    }

    /**
     * Runs the analysis for identifying variable genes.
     *
     * @param {boolean} saveDataset - Indicates whether to save the dataset.
     * @returns {Promise<void>} - A promise that resolves when the analysis is complete.
     */
    async runAnalysis(saveDataset) {
        if (Boolean(saveDataset)) {
            createToast("Saving variable genes", "is-info");
        } else {
            createToast("Analyzing variable genes", "is-info");
        }

        document.querySelector(UI.asvgPlotContainerElt).replaceChildren();
        document.querySelector(UI.asvgPlotNormContainerElt).replaceChildren();
        disableAndHideElement(document.querySelector(UI.btnAsvgSaveElt));

        this.normCountsPerCell = document.querySelector(UI.asvgNormCountsPerCellElt).value;
        this.flavor = document.querySelector(UI.asvgFlavorElt).value;
        this.nTopGenes = document.querySelector(UI.asvgNTopGenesElt).value;
        this.minMean = document.querySelector(UI.asvgMinMeanElt).value;
        this.maxMean = document.querySelector(UI.asvgMaxMeanElt).value;
        this.minDispersion = document.querySelector(UI.asvgMinDispersionElt).value;

        this.regressOut = document.querySelector(UI.asvgRegressOutElt).checked;
        this.scaleUnitVariance = document.querySelector(UI.asvgScaleUnitVarianceElt).checked;

        const params = {
            'dataset_id': this.analysis.datasetId,
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
            const {data} = await axios.post("./cgi/h5ad_identify_variable_genes.cgi", params);

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = Boolean(saveDataset);

            this.updateUIWithResults(data, this.calculated);
            createToast("Variable genes image created", "is-success");

            document.querySelector(UI.asvgResultCountElt).textContent = `(${data['n_genes']})`;
            enableAndShowElement(document.querySelector(UI.btnAsvgSaveElt));
            document.querySelector(UI.topGenesElt).textContent = `Suggested highly-variable genes:\n${data['top_genes']}`;
            document.querySelector(UI.topGenesElt).classList.remove("is-hidden");
            // TODO:  $('#asvg_options_c .js-next-step').show();  // auto-select the next collapsable dropdown

            document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).classList.add("is-hidden");
        } catch (error) {
            createToast(`Error identifying variable genes: ${error.message}`);
        } finally {
            document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).disabled = false;
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

        if (this.calculated) {
            document.querySelector(UI.selectVariableGenesToggleElt).checked = true;
            document.querySelector(UI.asvgNormCountsPerCellElt).value = this.normCountsPerCell;
            document.querySelector(UI.asvgFlavorElt).value = this.flavor;
            document.querySelector(UI.asvgNTopGenesElt).value = this.nTopGenes;
            document.querySelector(UI.asvgMinMeanElt).value = this.minMean;
            document.querySelector(UI.asvgMaxMeanElt).value = this.maxMean;
            document.querySelector(UI.asvgMinDispersionElt).value = this.minDispersion;

            document.querySelector(UI.asvgRegressOutElt).checked = false;
            if (this.regressOut) {
                document.querySelector(UI.asvgRegressOutElt).checked = true;
            }

            document.querySelector(UI.asvgScaleUnitVarianceElt).checked = false;
            if (this.scaleUnitVariance) {
                document.querySelector(UI.asvgScaleUnitVarianceElt).checked = true;
            }

            enableAndShowElement(document.querySelector(UI.btnAsvgSaveElt));
            if (resultsSaved) {
                document.querySelector(UI.btnAsvgSaveElt).disabled = true;
                document.querySelector(UI.btnAsvgSaveElt).textContent = 'Saved';

            }
            document.querySelector(UI.btnDoAnalysisSelectVariableGenesElt).classList.add("is-hidden");


            // Variable gene calculation enables PCA
            document.querySelector(UI.pcaToggleElt).classList.remove("is-hidden");

            for (const elt of document.querySelectorAll(UI.asvgSaveOptionsElts)) {
                elt.classList.remove("is-hidden");
            }
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'filter_genes_dispersion',
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
            'session_id': ana.userSessionId,

            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.placeAnalysisImage(
            {'params': params, 'title': 'Variable genes', 'target': UI.asvgPlotContainer});

        UI.asvgInstructionsElt.classList.add("is-hidden");


    }
}
class AnalysisStepPCA {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
        this.irreversible = false;
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
        this.genesToColor = false;
        this.resetUI();
    }

    /**
     * Resets the user interface by clearing the values of the genesToColorElt and topPcaGenesElt elements.
     */
    resetUI() {
        document.querySelector(UI.genesToColorElt).value = null;
        document.querySelector(UI.topPcaGenesElt).value = null;
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

        for (const container of document.querySelector(UI.pcaImageResultContainers)) {
            container.replaceChildren();
        }

        const computePCA = this.calculated ? false : true;

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_pca.cgi", {
                dataset_id: this.analysis.datasetId,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                genes_to_color: document.querySelector(UI.genesToColorElt).value,
                compute_pca: computePCA
            });
            if (!data.success || data.success < 1) {
                document.querySelector(UI.pcaMissingGeneElt).textContent = data['missing_gene'] || "";

                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = true;
            this.genesToColor = document.querySelector(UI.genesToColorElt).value;
            this.updateUIWithResults();
            createToast("PCA and variance computed", "is-success");

            document.querySelector(UI.pcaOptionsElt).classList.remove("is-hidden");
            document.querySelector(UI.weightedGeneCartElt).classList.remove("is-hidden");
            // TODO:  $('#pca_options_c .js-next-step').show();  // auto-select the next collapsable dropdown

        } catch (error) {
            createToast(`Error running PCA: ${error.message}`);
            document.querySelector(UI.pcaMissingGeneContainer).classList.remove("is-hidden");

        } finally {
            document.querySelector(UI.btnPcaRunElt).disabled = false;
        }
    }

    /**
     * Runs the analysis to compute the top genes for principal components.
     *
     * @async
     * @function runAnalysisTopGenes
     * @memberof analysis.v2
     * @returns {Promise<void>} A Promise that resolves when the analysis is complete.
     * @throws {Error} If there is an error computing the top PCA genes.
     */
    async runAnalysisTopGenes() {
        createToast("Computing top genes for principal components", "is-info");
        document.querySelector(UI.pcaTopGenesContainer).replaceChildren();

        try {

            const {data} = await axios.post("./cgi/h5ad_top_pca_genes.cgi", {
                dataset_id: this.analysis.datasetId,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                pcs: document.querySelector(UI.topPcaGenesElt).value
            });

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.updateUIWithResultsTopGenes(data);
            createToast("Top PCA genes computed", "is-success");

        } catch (error) {
            createToast(`Error computing top PCA genes: ${error.message}`);
        } finally {
            document.querySelector(UI.btnPcaTopGenesElt).disabled = false;
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

        document.querySelector(UI.instructionsElt).classList.add("is-hidden");

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'pca',
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
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
        document.querySelector(UI.tsneToggleElt).classList.remove("is-hidden");
        document.querySelector(UI.pcaToggleElt).checked = true;
        document.querySelector(UI.genesToColorElt).value = step.genesToColor;
        document.querySelector(UI.pcaOptionsDivElt).classList.remove("is-hidden");
    }
}

class AnalysisSteptSNE {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
        this.irreversible = false;
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
        this.nPcs = null;
        this.nNeighbors = null;
        this.randomState = null;
        this.plotTsne = 0;
        this.plotUmap = 0;
        this.resetUI();
    }

    /**
     * Resets the UI by clearing the values of various input elements.
     */
    resetUI() {
        document.querySelector(UI.tsneGenesToColorElt).value = '';
        document.querySelector(UI.dimReductionNNeighborsElt).value = '';
        document.querySelector(UI.tsneNPcsElt).value = '';
        document.querySelector(UI.tsneRandomStateElt).value = 2;
        document.querySelector(UI.dimReductionMethodTsneElt).checked = false;
        document.querySelector(UI.dimReductionMethodUmapElt).checked = true;
    }

    async runAnalysis() {
        createToast("Computing tSNE/UMAP and generating plot", "is-info");

        document.querySelector(UI.tsneMissingGeneContainer).classList.add("is-hidden");
        document.querySelector(UI.tsnePlotContainerElt).replaceChildren();
        document.querySelector(UI.umapPlotContainerElt).replaceChildren();

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

        const useScaled = document.querySelector(UI.tsneUseScaledElt).checked;

        const params = {
            'dataset_id': this.analysis.datasetId,
            'analysis_id': this.analysis.id,
            'analysis_type': this.analysis.type,
            'session_id': this.analysis.userSessionId,
            'genes_to_color': document.querySelector(UI.tsneGenesToColorElt).value,
            'n_pcs': document.querySelector(UI.tsneNPcsElt).value,
            'n_neighbors': document.querySelector(UI.dimReductionNNeighborsElt).value,
            'random_state': document.querySelector(UI.tsneRandomStateElt).value,
            'use_scaled': useScaled,
            'compute_neighbors': computeNeighbors,
            'compute_tsne': computeTsne,
            'compute_umap': computeUmap,
            'plot_tsne': plotTsne,
            'plot_umap': plotUmap
        }

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_tsne.cgi", params);

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
            // TODO:  $('#tsne_options_c .js-next-step').show();  // Reveal next collapsable dropdown


        } catch (error) {
            createToast(`Error generating tSNE: ${error.message}`);

            document.querySelector(UI.tsneMissingGeneContainer).classList.remove("is-hidden");

        } finally {
            document.querySelector(UI.btnTsneRunElt).disabled = false;
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
        document.querySelector(UI.tsneToggleElt).classList.remove("is-hidden");
        document.querySelector(UI.tsneGenesToColorElt).value = this.genesToColor;
        document.querySelector(UI.dimReductionNNeighborsElt).value = this.nNeighbors;
        document.querySelector(UI.tsneNPcsElt).value = this.nPcs;
        document.querySelector(UI.tsneRandomStateElt).value = this.randomState;

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
            'dataset_id': ana.datasetId,
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
        document.querySelector(UI.clusterToggleElt).classList.remove("is-hidden");

    }
}
class AnalysisStepClustering {
    constructor(analysis, mode="initial") {
        this.reset();   // TODO: split into two potentially for each step

        this.analysis = analysis;
        this.irreversible = false;
        this.mode = mode;  // initial, edit
        if (!(["initial", "edit"].includes(mode))) {
            logErrorInConsole("Invalid mode for AnalysisStepClustering. Defaulting to 'initial'.");
            this.mode = "initial";
        }
        if (mode == "edit") {
            this.calculated = true;
            this.irreversible = true;
        }

    }

    /**
     * Loads data from a JSON object and creates a new instance of AnalysisStepClustering.
     * @param {Object} data - The JSON object containing the data to load.
     * @param {Analysis} analysis - The analysis object to associate with the new instance.
     * @returns {AnalysisStepClustering} - The newly created instance of AnalysisStepClustering.
     */
    static loadFromJson(data, analysis) {
        const step = new AnalysisStepClustering(analysis);
        if (!step) return step;

        step.calculated = data['calculated'];
        step.nNeighbors = data['n_neighbors'];
        step.resolution = data['resolution'];
        step.plotTsne = data['plot_tsne'];
        step.plotUmap = data['plot_umap'];

        return step;
    }

    /**
     * Resets the state of the analysis object.
     */
    reset() {
        this.calculated = false;
        this.nNeighbors = null;
        this.resolution = null;
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

        // Reset resolution input back to default
        document.querySelector(UI.resolutionElt).value = 1.3;

        // TODO: On page js file, add click listener to sync the resolution inputs

    }

    /**
     * Runs the analysis by computing Louvain clusters and updating the UI with the results.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is completed.
     * @throws {Error} If there is an error generating clusters.
     */
    async runAnalysis() {
        // TODO: Fine-tune this based on if type is initial or edit

        createToast("Computing Louvain clusters", "is-info");
        for (const container of document.querySelector(UI.clusteringImageResultContainers)) {
            container.replaceChildren();
        }

        document.querySelector(UI.btnClusteringRunElt).disabled = true;
        if (this.type == "edit") {
            document.querySelector(UI.btnClusteringRerunWithGroupsElt).disabled = true;
        }

        const isSameLouvainParams = (this.nNeighbors == document.querySelector(UI.clusteringNNeighborsElt).value
            && this.resolution == document.querySelector(UI.resolutionElt).value);

        let computeClustering = true;

        const oldLabels = [...this.groupLabels];  // shallow-copy
        const newLabels = [];
        const keptLabels = [];

        // It is not safe to reuse group labels if the clustering params were changed
        if (isSameLouvainParams) {
            for (const glElt of clusterGroupLabelsHtml.querySelector(".group-user-label input")) {
                newLabels.push(glElt.value);
                const thisRow = glElt.closest("tr");
                const thisCheck = thisRow.querySelector(".group-keep-chk input");
                keptLabels.push(thisCheck.checked);
            }
        }

        const clusterInfo = [];

        if (isSameLouvainParams && this.calculated && this.type == "edit") {
            computeClustering = false;
            this.groupLabels.forEach((v, i) => {
                clusterInfo.push({
                    "old_label": oldLabels[i]
                    , "new_label": newLabels[i]
                    , "keep": keptLabels[i]
                })
            });
        }

        if (computeClustering) {
            document.querySelector(UI.groupLabelsContainer).classList.add("is-hidden");
        }

        const plotTsne = (document.querySelector(UI.dimReductionMethodTsneElt).checked ? 1 : 0);
        const plotUmap = (document.querySelector(UI.dimReductionMethodUmapElt).checked ? 1 : 0);

        try {
            const {data} = await axios.post("./cgi/h5ad_generate_clusters.cgi", {
                dataset_id: this.analysis.datasetId,
                analysis_id: this.analysis.id,
                analysis_type: this.analysis.type,
                session_id: this.analysis.userSessionId,
                resolution: document.querySelector(UI.resolutionElt).value,
                compute_clusters: computeClustering,
                plot_tsne: plotTsne,
                plot_umap: plotUmap,
                cluster_info: JSON.stringify(clusterInfo)
            });

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = true;
            this.nNeighbors = document.querySelector(UI.clusteringNNeighborsElt).value;
            this.resolution = document.querySelector(UI.resolutionElt).value;
            this.plotTsne = plotTsne;
            this.plotUmap = plotUmap;
            this.updateUIWithResults();

            if (data["group_labels"].length) {
                // Update the group labels for the analysis, marker genes, and gene comparison
                this.analysis.groupLabels = []
                this.analysis.markerGenes.populateMarkerGenesLabels(data)

                // Now update the labels so they work with gene comparison
                this.analysis.groupLabels = data['group_labels'].map(x => x.genes);
                this.analysis.geneComparison.populateGroupSelectors(this.analysis.groupLabels);
            }

            createToast("Louvain clusters computed", "is-success");
            // TODO:  $('#clustering_options_c .js-next-step').show();  // auto-select the next collapsable dropdown


        } catch (error) {
            createToast(`Error generating clusters: ${error.message}`);

        } finally {
            document.querySelector(UI.btnClusteringRunElt).disabled = false;
            if (this.type == "edit") {
                document.querySelector(UI.btnClusteringRerunWithGroupsElt).disabled = false;
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

        document.querySelector(UI.clusterInstructionsElt).classList.add("is-hidden");

        if (this.calculated) {
            document.querySelector(UI.clusterToggleElt).classList.checked = true

            document.querySelector(UI.clusteringNNeighborsElt).value = this.nNeighbors;
            document.querySelector(UI.resolutionElt).value = this.resolution;

        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': 'tsne_clustering',
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        // Need to ensure targets are different for the two different clustering steps
        let tsneTarget = UI.clusterTsnePlotElt;
        let umapTarget = UI.clusterUmapPlotElt;
        if (this.mode == "edit") {
            tsneTarget += UI.clusterTsnePlotEditElt;
            umapTarget += UI.clusterUmapPlotEditElt;
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

        // clustering enables marker gene identification
        document.querySelector(UI.markerGenesToggleElt).classList.remove("is-hidden");
    }
}
class AnalysisStepMarkerGenes {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;

        this.irreversible = false;
    }

    /**
     * Counts and highlights duplicate values in the clustering group labels step.
     * @returns {number} The number of duplicate values found.
     */
    countAndHighlightDuplicates() {
        const allValues = [];
        const dupValues = [];

        // ? Should this actually be in the cluster class

        // first remove any duplicate-labeled ones
        for (const elt of document.querySelectorAll(UI.clusterGroupLabelsInputElts)) {
            elt.classList.remove('duplicate');

            const clusterLabel = elt.value.trim();

            // this means it WAS found
            if (allValues.includes(clusterLabel)) {
                dupValues.push(clusterLabel);
                elt.classList.add('duplicate');
            } else {
                allValues.push(clusterLabel);
            }
        }

        return dupValues.length;
    }

    /**
     * Fetches marker genes for the analysis.
     * @returns {Promise<Object>} The response data containing the marker genes.
     * @throws {Error} If the request fails or returns an error.
     */
    async fetchMarkerGenes() {
        const {data} = await axios.post("./cgi/h5ad_find_marker_genes.cgi", {
            'dataset_id': this.analysis.datasetId,
            'analysis_id': this.analysis.id,
            'analysis_type': this.analysis.type,
            'session_id': this.analysis.userSessionId,
            'n_genes': this.nGenes,
            'compute_marker_genes': true
        });

        if ((!data.success) || (data.success < 1)) {
            const error = data.error || "Unknown error. Please contact gEAR support.";
            throw new Error(error);
        }
        return data;
    }

    static loadFromJson(data, analysis) {
        const step = new AnalysisStepMarkerGenes(analysis);
        if (!data) return step;

        step.calculated = data['calculated']
        step.nGenes = data['n_genes']
        step.groupNames = data['group_names']

        if (step.calculated ) {
            document.querySelector(UI.markerGenesToggleElt).checked = true;
            step.updateUIWithResults(data);
        }

        return step;
    }

    /**
     * Performs marker gene visualization by making a POST request to the server and updating the UI with the results.
     * @async
     * @function performMarkerGeneVisualization
     * @memberof analysis.v2
     * @throws {Error} If there is an error during the visualization process.
     */
    async performMarkerGeneVisualization() {

        document.querySelector(UI.markerGenesDotplotContainer).replaceChildren();
        document.querySelector(UI.markerGenesViolinContainer).replaceChildren();

        try {
            const data = await axios.post("./cgi/h5ad_generate_marker_gene_visualization.cgi", {
                'dataset_id': this.analysis.datasetId,
                'analysis_id': this.analysis.id,
                'analysis_type': this.analysis.type,
                'session_id': this.analysis.userSessionId,
                'marker_genes': JSON.stringify([...this.genesOfInterest])
            });

            if ((!data.success) || (data.success < 1)) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            const params = {
                'analysis_id': this.analysis.id,
                'analysis_name': 'dotplot_goi',
                'analysis_type': this.analysis.type,
                'dataset_id': this.analysis.datasetId,
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

        } catch {
            createToast("Error visualizing marker genes");
        }
    }

    /**
     * Populates the marker genes labels in the analysis.
     *
     * @param {Object} data - The data containing the group labels.
     */
    populateMarkerGenesLabels(data) {
        this.groupLabels = [];

        // If the user has saved labels before, put them in the table here.  Else leave it
        //  as the top gene.
        let i = 0;
        if (this.analysis.groupLabels.length > 0) {
            for (i=0; i < this.analysis.groupLabels.length; i++) {
                data['group_labels'][i]['new_group_label'] = this.analysis.groupLabels[i];
                this.groupLabels.push(this.analysis.groupLabels[i]);
            }
        } else {
            for (i=0; i < data['group_labels'].length; i++) {
                data['group_labels'][i]['new_group_label'] = data['group_labels'][i]['genes'];
                // For the overall labels, do the group number rather than gene since that's what's
                //  displayed by scanpy in the images
                this.groupLabels.push(i);
            }
        }

        // show the abbreviated table in the louvain analysis block
        for (const group of data['group_labels']) {
            const clusterGroupLabelsHtml = document.querySelector(UI.clusterGroupLabelsTmpl).cloneNode(true);
            clusterGroupLabelsHtml.querySelector(".group-orig-label").textContent = group['group_label'];
            clusterGroupLabelsHtml.querySelector(".group-num-cells").value = group['num_cells'];
            clusterGroupLabelsHtml.querySelector(".group-marker").value = group['genes'];
            clusterGroupLabelsHtml.querySelector(".group-user-label input").value = group['new_group_label'];

            document.querySelector(UI.clusterGroupLabelsTableBodyElt).appendChild(clusterGroupLabelsHtml);
        }
        this.groupLabelsContainer.classList.remove("is-hidden");
    }

    /**
     * Populates the marker genes table with data retrieved from the server.
     *
     * @returns {Promise<void>} - A promise that resolves when the marker genes table is populated.
     * @throws {Error} - If there is an error retrieving the marker genes data.
     */
    async populateMarkerGenesTable(table) {
        try {

            if (!table) {
                const data = await this.fetchMarkerGenes();
                table = data['table'];
            }

            // Add the first th element (empty)
            const markerGenesHeaderHtml = document.querySelector(UI.markerGenesTableHeadTmpl).content.cloneNode(true);
            markerGenesHeaderHtml.querySelector("th").textContent = "";
            document.querySelector(UI.markerGenesTableHeadRowElt).appendChild(markerGenesHeaderHtml);

            // add the table header
            for (const column of data['table']['columns']) {
                const markerGenesHeaderHtml = document.querySelector(UI.markerGenesTableHeadTmpl).content.cloneNode(true);
                markerGenesHeaderHtml.querySelector("th").textContent = column;
                document.querySelector(UI.markerGenesTableHeadRowElt).appendChild(markerGenesHeaderHtml);
            }

            // add the table rows
            for (const row of data['table']['rows']) {
                // Create the row label
                const markerGenesBodyHtml = document.querySelector(UI.markerGenesTableBodyTmpl).content.cloneNode(true);
                markerGenesBodyHtml.querySelector(".js-row-idx").textContent = row["rowid"]

                // create new td objects for each column
                for (const column of row["columns"]) {
                    const cellElt = document.createElement("td");
                    cellElt.textContent = column["label"]
                    markerGenesBodyHtml.appendChild(cellElt)
                }

                document.querySelector(UI.markerGenesTableBodyElt).appendChild(markerGenesBodyHtml);
            }

            this.populateMarkerGenesLabels(data);
            const groupLabels = data['group_labels'].map(x => x.group_label);
            this.analysis.compareGenes.populateGroupSelectors(groupLabels);

            document.querySelector(UI.downloadMarkerGenesBtnElt).classList.remove("is-hidden");
        } catch (error) {
            createToast(`Error getting marker genes: ${error.message}`);
        }

    }

    /**
     * Resets the analysis state.
     */
    reset() {
        this.calculated = false;
        this.genesOfInterest = new Set();
        this.nGenes = false;
        this.resetUI();
    }

    /**
     * Resets the user interface.
     */
    resetUI() {
        document.querySelector(UI.markerGenesNGenesElt).value = 5;
    }

    /**
     * Runs the analysis to compute marker genes.
     *
     * @returns {Promise<void>} A promise that resolves when the analysis is complete.
     * @throws {Error} If there is an error computing marker genes.
     */
    async runAnalysis() {
        createToast("Computing marker genes", "is-info");
        document.querySelector(UI.btnMarkerGenesRunElt).disabled = true;
        document.querySelector(UI.markerGenesPlotContainer).replaceChildren();
        document.querySelector(UI.markerGenesTableElt).replaceChildren();

        document.querySelector(UI.markerGenesManuallyEnteredElt).value = '';
        document.querySelector(UI.markerGenesSelectedCountElt).textContent = 0;
        document.querySelector(UI.markerGenesEnteredCountElt).textContent = 0;
        document.querySelector(UI.markerGenesUniqueCountElt).textContent = 0;

        this.genesOfInterest = new Set();
        this.clickedMarkerGenes = new Set();
        this.enteredMarkerGenes = new Set();

        const computeMarkerGenes = true
        if (this.calculated && this.nGenes == document.querySelector(UI.markerGenesNGenesElt).value) {
            computeMarkerGenes = false;
        }

        try {
            const {data} = await axios.post("./cgi/h5ad_find_marker_genes.cgi", {
                'dataset_id': this.analysis.datasetId,
                'analysis_id': this.analysis.id,
                'analysis_type': this.analysis.type,
                'session_id': this.analysis.userSessionId,
                'n_genes': document.querySelector(UI.markerGenesNGenesElt).value,
                'compute_marker_genes': computeMarkerGenes
            });

            if (!data.success || data.success < 1) {
                const error = data.error || "Unknown error. Please contact gEAR support.";
                throw new Error(error);
            }

            this.calculated = true;
            this.nGenes = document.querySelector(UI.markerGenesNGenesElt).value;
            this.updateUIWithResults(data);
            this.groupLabels = data['group_labels'].map(x => x.group_label);
        } catch (error) {
            createToast(`Error computing marker genes: ${error.message}`);
        } finally {
            document.querySelector(UI.btnMarkerGenesRunElt).disabled = false;
        }
    }

    /**
     * Updates the UI with the results of the analysis.
     *
     * @param {Object} data - The data containing the analysis results.
     * @param {Object} [ana=null] - The analysis object. If not provided, the method uses the analysis object of the class instance.
     * @returns {void}
     */
    updateUIWithResults(data, ana=null) {
        if (!ana) {
            ana = this.analysis;
        }

        if (!ana) {
            createToast("No analysis object found. Cannot update UI.")
            return;
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': `rank_genes_groups_${data['cluster_label']}`,
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        document.querySelector(UI.markerGenesNGenesElt).value = this.nGenes;
        document.querySelector(UI.markerGenesTableHeader).classList.remove("is-hidden");

        ana.placeAnalysisImage(
            {'params': params, 'title': 'Marker genes', 'target': UI.markerGenesPlotContainer});

        if (data['table']) {
            this.populateMarkerGenesTable(data['table']);

            this.populateMarkerGenesLabels(data);
            const groupLabels = data['group_labels'].map(x => x.group_label);
            ana.gene_comparison.populateGroupSelectors(groupLabels);
        } else {
            this.populateMarkerGenesTable();    // This will fetch the data
        }

        document.querySelector(UI.markerGenesVisualizationContainer).classList.remove("is-hidden");

        // marker gene calculation enables cluster comparison
        document.querySelector(UI.compareGenesToggleElt).classList.remove("is-hidden");

    }
}

class AnalysisStepCompareGenes {
    constructor(analysis) {
        this.reset();
        this.analysis = analysis;
        this.irreversible = false;
    }

    /**
     * Loads an instance of AnalysisStepCompareGenes from JSON data.
     *
     * @param {Object} data - The JSON data to load from.
     * @returns {AnalysisStepCompareGenes} - The loaded AnalysisStepCompareGenes instance.
     */
    static loadFromJson(data) {
        const step = new AnalysisStepCompareGenes(analysis);
        if (!data) {
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

            const clusterOptionHtml = document.querySelector(UI.clusterOptsTmpl).content.cloneNode(true);
            const option = clusterOptionHtml.querySelector("option");
            option.textContent = label;
            option.value = groupNum;

            document.querySelector(UI.queryClusterOptionsElt).appendChild(clusterOptionHtml);
            document.querySelector(UI.referenceClusterOptionsElt).appendChild(clusterOptionHtml);

            groupNum++;
        }

        document.querySelector(UI.queryClusterSelectElt).value = this.queryCluster;
        document.querySelector(UI.referenceClusterSelectElt).value = this.referenceCluster;
    }

    /**
     * Resets the analysis object to its initial state.
     */
    reset() {
        this.calculated = false;
        this.resetUI();
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
    resetUI() {

        document.querySelector(UI.queryClusterSelectElt).value = null;
        document.querySelector(UI.referenceClusterSelectElt).value = "all-reference-clusters";

        document.querySelector(UI.compareGenesNGenesElt).value = null;
        document.querySelector(UI.compareGenesMethodElt).value = "t-test_overestim_var";
        document.querySelector(UI.comapreGenesCorrMethodElt).value = "benjamini-hochberg";

        document.querySelector(UI.compareGenesInstructionsElt).classList.remove("is-hidden");
        document.querySelector(UI.compareGenesRankedContainer).classList.add("is-hidden");
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
        for (const container of document.querySelectorAll(UI.compareGenesResultsContainers)) {
            container.replaceChildren();
        }

        const computeGeneComparison = this.calculated ? 1 : 0;

        try {
            const {data} = await axios.post("./cgi/h5ad_compare_genes.cgi", {
                'dataset_id': this.analysis.datasetId,
                'analysis_id': this.analysis.id,
                'analysis_type': this.analysis.type,
                'session_id': this.analysis.userSessionId,
                'n_genes': document.querySelector(UI.compareGenesNGenesElt).value,
                'compute_gene_comparison': computeGeneComparison,
                'group_labels': JSON.stringify(this.analysis.groupLabels),
                'query_cluster': document.querySelector(UI.queryClusterSelectElt).value,
                'reference_cluster': document.querySelector(UI.referenceClusterSelectElt).value,
                'method': document.querySelector(UI.compareGenesMethodElt).value,
                'corr_method': document.querySelector(UI.comapreGenesCorrMethodElt).value
            });

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
            this.method = document.querySelector(UI.compareGenesMethodElt).value;
            this.corrMethod = document.querySelector(UI.comapreGenesCorrMethodElt).value;
            this.updateUIWithResults(data);
            createToast("Comparison computed", "is-success");
        } catch (error) {
            createToast(`Error computing comparison: ${error.message}`);
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
            document.querySelector(UI.compareGenesToggleElt).checked = true;

            document.querySelector(UI.compareGenesNGenesElt).value = this.nGenes;
            document.querySelector(UI.compareGenesMethodElt).value = this.method;
            document.querySelector(UI.comapreGenesCorrMethodElt).value = this.corrMethod;

            // the query and reference cluster values have to be loaded after the option lists
        }

        const params = {
            'analysis_id': ana.id,
            'analysis_name': `rank_genes_groups_${data['cluster_label']}_comp_ranked`,
            'analysis_type': ana.type,
            'dataset_id': ana.datasetId,
            'session_id': ana.userSessionId,
            // this saves the user from getting a cached image each time
            'datetime': (new Date()).getTime()
        }

        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesRankedContainer});

        params['analysis_name'] = `rank_genes_groups_${data['cluster_label']}_${this.queryCluster}_comp_violin`

        ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesViolinContainer});

        if (data.hasOwnProperty('table_json_f')) {
            this.tableJsonF = data['table_json_f'];
            data['table_json_f'] = JSON.parse(data['table_json_f']);
            this.populateComparisonTable('compare_genes_table_f', data['table_json_f']['data']);
        }

        if (data.hasOwnProperty('table_json_r')) {
            this.tableJsonR = data['table_json_r'];
            data['table_json_r'] = JSON.parse(data['table_json_r']);
            this.populateComparisonTable('compare_genes_table_r', data['table_json_r']['data']);
        }

        if (this.referenceCluster != 'all-reference-clusters') {
            params['analysis_name'] =`rank_genes_groups_${data['cluster_label']}_comp_ranked_rev`
            ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesRankedRevContainer});

            params["analysis_name"] = `rank_genes_groups_${data['cluster_label']}_${this.referenceCluster}_comp_violin_rev`
            ana.placeAnalysisImage({'params': params, 'title': 'Comparison with a cluster', 'target': UI.compareGenesViolinRevContainer});
        }


    }
}
