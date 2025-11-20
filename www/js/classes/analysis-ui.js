'use strict';

import { blockStepWithHref, openNextStepWithHrefs } from "../stepper-fxns.js?v=4c8340e";

class AnalysisUI {
    // This class is a singleton that manages the UI elements of the analysis pipeline.

    // TODO: Reorganize these elements into logical groups beyond the steps of the analysis pipeline.
    // This would be something like section, container, button, element, etc.

    // dataset elemnents
    datasetSection = "#dataset-s"
    datasetContainer = "#dataset-c"
    datasetSelectionContainer = "#dataset-selection-c"
    datasetSectionSuccessElt = "#dataset-s-success"
    datasetSectionFailedElt = "#dataset-s-failed"
    currentDatasetContainer = "#current-dataset-c"
    currentDatasetElts = ".js-current-dataset"
    datasetQueryElt = "#dataset-query"
    datasetTreeElt = "#dataset-tree"
    btnToggleDatasetTreeElt = "#btn-toggle-dataset-tree"

    // general analysis
    analysisSection = "#analysis-s"
    analysisContainer = "#analysis-c"
    currentAnalysisElts = ".js-current-analysis"
    analysisSelect = "#analysis-select"
    analysisSelectSuccessElt = "#analysis-type-select-c-success"
    analysisSelectFailedElt = "#analysis-type-select-c-failed"
    newAnalysisOptionElt = `${this.analysisSelect} option[data-analysis-id='0']`
    analysisPrimaryElt = "#analyses-primary"
    analysisPrimaryNotificationElt = "#analysis-primary-notification"
    analysisUnsavedElt = "#analyses-unsaved"
    analysisSavedElt = "#analyses-saved"
    analysisPublicElt = "#analyses-public"
    analysisRenameElts = ".js-show-rename-input"
    analysisDeleteElts = ".js-delete-analysis"
    analysisDownloadElts = ".js-download-analysis"
    analysisActionContainer = "#analysis-action-c"
    btnSaveAnalysisElt = "#btn-save-analysis"
    btnDeleteUnsavedAnalysisElt = "#btn-delete-unsaved-analysis"
    analysisStatusInfoContainer = "#analysis-status-info-c"
    analysisStatusInfoElt = "#analysis-status-info"
    btnMakePublicCopyElt = "#btn-make-public-copy"
    btnDeleteSavedAnalysisElt = "#btn-delete-saved-analysis"
    analysisWorkflowElt = "#analysis-workflow"

    // summary section
    summarySection = "#summary-s"
    primaryStepsElt = "#primary-steps"
    deNovoStepsElt = "#de-novo-steps"
    stepsElts = "#steps-container .steps"
    stepSegmentElts = ".steps-segment"
    btnProgressGuideElt = "#btn-progress-guide"
    progressGuideElt = "#progress-guide"

    // Initial info
    primaryInitialInfoSection = "#initial-info-s"
    primaryInitialPlotContainer = "#initial-plot-c"
    primaryInitialLoadingElt = "#initial-loading-c"
    primaryInitialScatterContainer = "#initial-scatter-c"
    primaryInitialViolinContainer = "#initial-violin-c"
    selectedDatasetShapeInitialElt = "#selected-dataset-shape-initial"

    // labeled tSNE
    labeledTsneSection = "#labeled-tsne-s"
    labeledTsneContainer = "#labeled-tsne-c"
    btnLabeledTsneRunElt = "#btn-labeled-tsne-run"
    labeledTsnePlotContainer = "#labeled-tsne-plot-c"
    labeledTsneGeneSymbolElt = "#labeled-tsne-gene-symbol"

    // primary filter
    primaryFilterSection = "#primary-filter-s"
    primaryFilterSectionSuccessElt = "#primary-filter-s-success"
    primaryFilterSectionFailedElt = "#primary-filter-s-failed"
    btnApplyPrimaryFilterElt = "#btn-apply-primary-filter"
    filterCellsGtNGenesElt = "#filter-cells-gt-n-genes"
    filterCellsLtNGenesElt = "#filter-cells-lt-n-genes"
    filterGenesGtNCellsElt = "#filter-genes-gt-n-cells"
    filterGenesLtNCellsElt = "#filter-genes-lt-n-cells"
    filterCellsGtNGenesSelectedElt = "#filter-cells-gt-n-genes-selected"
    filterCellsLtNGenesSelectedElt = "#filter-cells-lt-n-genes-selected"
    filterGenesGtNCellsSelectedElt = "#filter-genes-gt-n-cells-selected"
    filterGenesLtNCellsSelectedElt = "#filter-genes-lt-n-cells-selected"
    selectedDatasetShapeFilteredElt = "#selected-dataset-shape-filtered"
    selectedDatasetShapeFilteredContainer = "#selected-dataset-shape-filtered-c"

    primaryTopGenesContainer = "#primary-top-genes-c"
    primaryTopGenesPlotContainer = "#primary-top-genes-plot-c"
    primaryFilterInstructionsElt = `${this.primaryFilterSection} .tool-instructions`
    datasetInfoResetableElts = `${this.primaryFilterSection} .js-resetable`

    // QC by mito
    qcByMitoSection = "#qc-by-mito-s"
    qcByMitoSectionSuccessElt = "#qc-by-mito-s-success"
    qcByMitoSectionFailedElt = "#qc-by-mito-s-failed"
    qbmGenePrefixElt = "#qbm-gene-prefix"
    qbmFilterMitoPercElt = "#qbm-filter-mito-perc"
    qbmFilterMitoCountElt = "#qbm-filter-mito-count"
    btnQbmSaveElt = "#btn-qbm-save"
    qbmSaveWarningElt = "#qbm-save-warning"
    btnDoAnalysisQcByMitoElt = "#btn-do-analysis-qc-by-mito"
    qbmPostShapeContainer = "#qbm-post-shape-c"
    qbmPostShapeElt = "#qbm-post-shape"
    qbmInstructionsElt = `${this.qcByMitoSection} .tool-instructions`
    qbmViolinContainer = "#qbm-violin-c"
    qbmScatterPercentMitoContainer = "#qbm-scatter-percent-mito-c"
    qbmScatterNGenesContainer = "#qbm-scatter-n-genes-c"
    qbmResetableElts = `${this.qcByMitoSection} .js-resetable`

    // select variable genes
    selectVariableGenesSection = "#select-variable-genes-s"
    selectVariableGenesSectionSuccessElt = "#select-variable-genes-s-success"
    selectVariableGenesSectionFailedElt = "#select-variable-genes-s-failed"
    asvgNormCountsPerCellElt = "#asvg-norm-counts-per-cell"
    asvgFlavorElt = "#asvg-flavor"
    asvgNTopGenesElt = "#asvg-n-top-genes"
    asvgMinMeanElt = "#asvg-min-mean"
    asvgMaxMeanElt = "#asvg-max-mean"
    asvgMinDispersionElt = "#asvg-min-dispersion"
    btnAsvgSaveElt = "#btn-asvg-save"
    btnDoAnalysisSelectVariableGenesElt = "#btn-do-analysis-select-variable-genes"
    asvgSaveWarningElt = "#asvg-save-warning"
    asvgPlotContainer = "#asvg-plot-c"
    asvgPostShapeContainer = "#asvg-post-shape-c"
    asvgPostShapeElt = "#asvg-post-shape"
    asvgInstructionsElt = `${this.selectVariableGenesSection} .tool-instructions`
    asvgTopGenesContainer = "#asvg-top-genes-c"
    asvgTopGenesListElt = "#asvg-top-genes-list"

    // PCA
    pcaSection = "#pca-s"
    pcaSectionSuccessElt = "#pca-s-success"
    pcaSectionFailedElt = "#pca-s-failed"
    btnPcaRunElt = "#btn-pca-run"
    btnPcaTopGenesElt = "#btn-pca-top-genes"
    btnSavePcaGeneListElt = "#btn-save-pca-gene-list"
    pcaGeneListNameElt = "#pca-gene-list-name"
    pcaGenesToColorElt = "#pca-genes-to-color"
    topPcaGenesElt = "#pca-top-genes"
    pcaInstructionsElt = `${this.pcaSection} .tool-instructions`
    pcaScatterContainer = "#pca-scatter-c"
    pcaVarianceContainer = "#pca-variance-c"
    pcaMissingGeneContainer = "#pca-missing-gene-c"
    pcaMissingGeneElt = "#pca-missing-gene"
    pcaResetableElts  = `${this.pcaSection} .js-resetable`
    pcaPcToTopGenesContainer = "#pca-pc-to-top-genes-c"
    pcaTopGenesPlotContainer = "#pca-top-genes-plot-c"
    pcaGeneListContainer = "#pca-gene-list-c"

    // tSNE
    tsneSection = "#tsne-s"
    tsneSectionSuccessElt = "#tsne-s-success"
    tsneSectionFailedElt = "#tsne-s-failed"
    btnTsneRunElt = "#btn-tsne-run"
    tsneGenesToColorElt = "#tsne-genes-to-color"
    dimReductionNNeighborsElt = "#dredux-n-neighbors"
    tsneNPcsElt = "#tsne-n-pcs"
    dimReductionMethodTsneElt = "#dimensionality-reduction-method-tsne"
    dimReductionMethodUmapElt = "#dimensionality-reduction-method-umap"
    tsneInstructionsElt = `${this.tsneSection} .tool-instructions`
    tsneResetableElts  = `${this.tsneSection} .js-resetable`
    tsnePlotContainer = "#tsne-plot-c"
    umapPlotContainer = "#umap-plot-c"
    tsneMissingGeneContainer = "#tsne-missing-gene-c"
    tsneMissingGeneElt = "#tsne-missing-gene"

    // Clustering
    clusteringSection = "#clustering-s"
    clusteringSectionSuccessElt = "#clustering-s-success"
    clusteringSectionFailedElt = "#clustering-s-failed"
    btnClusteringRunElt = "#btn-clustering-run"
    btnClusteringRerunWithGroupsElt = "#btn-clustering-rerun-with-groups"
    resolutionElt = "#clustering-resolution"
    clusteringTsnePlotElt = "#clustering-tsne-plot-c"
    clusteringUmapPlotElt = "#clustering-umap-plot-c"
    clusteringInstructionsElt = `${this.clusteringSection} .tool-instructions`
    clusteringResetableElts = `${this.clusteringSection} .js-resetable`

    // Marker Genes
    markerGenesSection = "#marker-genes-s"
    markerGenesSectionSuccessElt = "#marker-genes-s-success"
    markerGenesSectionFailedElt = "#marker-genes-s-failed"
    btnMarkerGenesRunElt = "#btn-marker-genes-run"
    btnVisualizeMarkerGenesElt = "#btn-visualize-marker-genes"
    btnDownloadMarkerGenesElt = "#btn-download-marker-genes"
    btnSaveMarkerGeneListElt = "#btn-save-marker-gene-list"
    markerGenesNGenesElt = "#marker-genes-n-genes"
    markerGenesTableContainer = "#marker-genes-table-c"
    markerGenesTableElt = "#marker-genes-table"
    markerGenesTableHeadTmpl = "#marker-genes-table-head-tmpl"
    markerGenesTableHeadRowElt = `${this.markerGenesTableElt} thead tr`
    markerGenesTableHeadCellElts = `${this.markerGenesTableElt} thead th`
    markerGenesTableRowTmpl = "#marker-genes-table-row-tmpl"
    markerGenesTableBodyElt = `${this.markerGenesTableElt} tbody`
    markerGenesTableHighlightedElts = `${this.markerGenesTableElt} td.js-highlighted`
    markerGenesPlotContainer = "#marker-genes-plot-c"
    markerGenesVisualizationContainer = "#marker-genes-visualization-c"
    markerGenesDotplotContainer = "#marker-genes-dotplot-c"
    markerGenesViolinContainer = "#marker-genes-violin-c"
    markerGenesManuallyEnteredElt = "#marker-genes-manually-entered"
    markerGenesSelectedCountElt = "#marker-genes-selected-count"
    markerGenesEnteredCountElt = "#marker-genes-entered-count"
    markerGenesUniqueCountElt = "#marker-genes-unique-count"
    markerGenesListContainer = "#marker-gene-list-c"
    markerGenesListNameElt = "#marker-gene-list-name"
    markerGenesInstructionsElt = `${this.markerGenesSection} .tool-instructions`
    markerGenesResetableElts = `${this.markerGenesSection} .js-resetable`

    // Clustering (edit mode)
    clusteringEditSection = "#clustering-edit-s"
    btnClusteringEditRunElt = "#btn-clustering-edit"
    clusteringEditSectionSuccessElt = "#clustering-edit-s-success"
    clusteringEditSectionFailedElt = "#clustering-edit-s-failed"
    // -- using resolutionElts values from the non-edit clustering
    // -- using clusterNNeighborsElt values from the non-edit clustering
    clusteringTsnePlotEditElt = "#clustering-tsne-plot-edit-c"
    clusteringUmapPlotEditElt = "#clustering-umap-plot-edit-c"
    groupLabelsContainer = "#group-labels-c"
    clusterGroupLabelsTmpl = "#cluster-group-labels-tmpl"
    clusterGroupLabelsTableBodyElt = "#cluster-group-labels tbody"
    clusterGroupLabelsInputElts = '#cluster-group-labels td.group-user-label input'
    clusteringMergeClustersElt = "#clustering-merge-clusters"
    clusteringEditInstructionsElt = `${this.clusteringEditSection} .tool-instructions`
    clusteringEditResetableElts = `${this.clusteringEditSection} .js-resetable`

    // Compare Genes
    compareGenesSection = "#compare-genes-s"
    compareGenesSectionSuccessElt = "#compare-genes-s-success"
    compareGenesSectionFailedElt = "#compare-genes-s-failed"
    btnCompareGenesRunElt = "#btn-compare-genes-run"
    btnCompareGenesShowTableFElt = "#btn-compare-genes-show-table-f"
    btnCompareGenesShowTableRElt = "#btn-compare-genes-show-table-r"
    btnCompareGenesDownloadTableFElt = "#btn-compare-genes-download-table-f"
    btnCompareGenesDownloadTableRElt = "#btn-compare-genes-download-table-r"
    queryClusterOptionsElt = "#query-cluster-options"
    referenceClusterOptionsElt = "#reference-cluster-options"
    queryClusterSelectElt = "#query-cluster"
    referenceClusterSelectElt = "#reference-cluster"
    compareGenesNGenesElt = "#compare-genes-n-genes"
    compareGenesMethodSelectElt = "#compare-genes-method"
    comapreGenesCorrMethodSelectElt = "#compare-genes-corr-method"
    compareGenesRankedContainer = "#compare-genes-ranked-c"
    compareGenesViolinContainer = "#compare-genes-violin-c"
    compareGenesRankedRevContainer = "#compare-genes-ranked-rev-c"
    compareGenesViolinRevContainer = "#compare-genes-violin-rev-c"
    compareGenesTableContainerF = "#compare-genes-table-container-f"
    compareGenesTableContainerR = "#compare-genes-table-container-r"
    compareGenesTableFElt = "#compare-genes-table-f"
    compareGenesTableRElt = "#compare-genes-table-r"
    compareGenesInstructionsElt = `${this.compareGenesSection} .tool-instructions`
    compareGenesResultsContainer = "#compare-genes-visualization-c"
    compareGenesResetableElts = `${this.compareGenesSection} .js-resetable`

}

// Singleton instance of the AnalysisUI class
export const UI = new AnalysisUI();

/**
 * Blocks a step from being opened, such as an irreversible step.
 * @param {string} selector - The CSS selector for the step element.
 */
export const blockAnalysisStep = (selector) => {
    blockStepWithHref(selector);
    // disable all buttons in the blocked step
    document.querySelectorAll(`${selector} button`).forEach((button) => {
        button.disabled = true;
    });
}

/**
 * Opens the next step in the UI.
 *
 * @param {Array} selectors - An array of CSS selectors for the UI elements.
 * @param {string|null} activeSelectorHref - The href of the active selector. Defaults to null.
 * @param {boolean} clickActive - Specifies whether to click the active selector. Defaults to false.
 */
export const openNextAnalysisStep = (selectors, activeSelectorHref=null, clickActive=false) => {
    for (const selector of selectors) {
        document.querySelector(selector).classList.remove("is-pointer-events-none");
    }
    openNextStepWithHrefs(selectors, activeSelectorHref, clickActive);
}
