'use strict';

class AnalysisUI {
    // This class is a singleton that manages the UI elements of the analysis pipeline.

    // dataset elemnents
    datasetSection = "#dataset-s"
    datasetContainer = "#dataset-c"
    datasetSectionSuccessElt = "#dataset-s-success"
    datasetSectionFailedElt = "#dataset-s-failed"
    currentDatasetContainer = "#current-dataset-c"
    currentDatasetElt = "#current-dataset"
    datasetQueryElt = "#dataset-query"
    datasetTreeElt = "#dataset-tree"

    // general analysis
    analysisSection = "#analysis-s"
    analysisContainer = "#analysis-c"
    emptyAnalysisOptionTmpl = "#analyses-list-empty-tmpl"
    analysisOptionTmpl = "#analyses-list-tmpl"
    analysisSelect = "#analysis-select"
    newAnalysisOptionElt = `${this.analysisSelect} option[data-analysis-id='0']`
    newAnalysisLabelContainer = "#new-analysis-label-c"
    newAnalysisLabelElt = "#new-analysis-label"
    btnNewAnalysisLabelSaveElt = "#btn-new-analysis-label-save"
    btnNewAnalysisLabelCancelElt = "#btn-new-analysis-label-cancel"
    duplicateLabelWarningElt = "#duplicate-label-warning"
    analysisPrimaryElt = "#analyses-primary"
    analysisPrimaryNotificationElt = "#analyses-primary-notification"
    analysisUnsavedElt = "#analyses-unsaved"
    analysisSavedElt = "#analyses-saved"
    analysisPublicElt = "#analyses-public"
    analysisRenameElts = ".js-show-rename-input"
    analysisActionContainer = "#analysis-action-c"
    btnSaveAnalysisElt = "#btn-save-analysis"
    btnDeleteUnsavedAnalysisElt = "#btn-delete-unsaved-analysis"
    analysisStatusInfoContainer = "#analysis-status-info-c"
    analysisStatusInfoElt = "#analysis-status-info"
    btnMakePublicCopy = "#btn-make-public-copy"
    btnDeleteSavedAnalysisElt = "#btn-delete-saved-analysis"
    storedAnalysesContainer = "#stored-analyses-c"
    initialInstructionsElt = ".initial-instructions"
    analysisWorkflowElt = "#analysis-workflow"

    // primary filter
    primaryFilterSection = "#primary-filter-s"
    primaryFilterSectionSuccessElt = "#primary-filter-s-success"
    primaryFilterSectionFailedElt = "#primary-filter-s-failed"
    primaryFilterCollapsableElt = "#primary-filter-collapsable"
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
    primaryInitialPlotContainer = "#primary-initial-plot-c"
    primaryInitialScatterContainer = "#primary-initial-scatter-c"
    primaryInitialViolinContainer = "#primary-initial-violin-c"
    primaryTopGenesContainer = "#primary-top-genes-c"
    primaryTopGenesPlotContainer = "#primary-top-genes-plot-c"
    primaryInitialPlotElts = ".primary-initial-plot"
    datasetInfoResetableElts = `${this.primaryFilterCollapsableElt} .js-resetable`

    // QC by mito
    qcByMitoSection = "#qc-by-mito-s"
    qcByMitoSectionSuccessElt = "#qc-by-mito-s-success"
    qcByMitoSectionFailedElt = "#qc-by-mito-s-failed"
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
    qbmResetableElts = `${this.qcByMitoSection} .js-resetable`

    // select variable genes
    //selectVariableGenesToggleElt = "#toggle-select-variable-genes"  // Temporary
    selectVariableGenesSection = "#select-variable-genes-s"
    selectVariableGenesSectionSuccessElt = "#select-variable-genes-s-success"
    selectVariableGenesSectionFailedElt = "#select-variable-genes-s-failed"
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
    asvgTopGenesListElt = "#asvg-top-genes-list"

    // PCA
    //pcaToggleElt = "#toggle-pca"   // Temporary
    pcaSection = "#pca-s"
    pcaSectionSuccessElt = "#pca-s-success"
    pcaSectionFailedElt = "#pca-s-failed"
    btnPcaRunElt = "#btn-pca-run"
    btnSavePcaGeneListElt = "#btn-save-pca-gene-list"
    pcaGeneListNameElt = "#pca-gene-list-name"
    genesToColorElt = "#pca-genes-to-color"
    pcaOptionsGroupElt = "#pca-options-group"
    topPcaGenesElt = "#pca-top-genes"
    pcaInstructionsElt = "#analysis-pca .tool-instructions"
    pcaScatterContainer = "#pca-scatter-c"
    pcaVarianceContainer = "#pca-variance-c"
    pcaMissingGeneContainer = "#pca-missing-gene-c"
    pcaMissingGeneElt = "#pca-missing-gene"
    pcaResetableElts  = `1#analysis-pca .js-resetable`
    pcaOptionsGroupElt = "#pca-options-g"
    weightedGeneListGroupElt = "#weighted-gene-list-g"


    // tSNE
    //tsneToggleElt = "#toggle-tsne"   // Temporary
    tsneSection = "#tsne-s"
    tsneSectionSuccessElt = "#tsne-s-success"
    tsneSectionFailedElt = "#tsne-s-failed"
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
    //clusteringToggleElt = "#toggle-clustering"   // Temporary
    clusteringSection = "#clustering-s"
    clusteringSectionSuccessElt = "#clustering-s-success"
    clusteringSectionFailedElt = "#clustering-s-failed"
    btnClusteringRunElt = "#btn-clustering-run"
    btnClusteringRerunWithGroupsElt = "#btn-clustering-rerun-with-groups"
    resolutionElt = "#clustering-resolution"
    clusteringNNeighborsElt = "#clustering-n-neighbors"
    clusterTsnePlotElt = "#cluster-tsne-plot-c"
    clusterUmapPlotElt = "#cluster-umap-plot-c"
    clusteringInstructionsElt = "#analysis-clustering .tool-instructions"
    clusteringResetableElts = `#analysis-clustering .js-resetable`

    // Marker Genes
    //markerGenesToggleElt = "#toggle_marker_genes"   // Temporary
    markerGenesSection = "#marker-genes-s"
    markerGenesSectionSuccessElt = "#marker-genes-s-success"
    markerGenesSectionFailedElt = "#marker-genes-s-failed"
    btnMarkerGenesRunElt = "#btn-marker-genes-run"
    visualMarkerGenesSection = "#visualize-marker-genes-s"
    btnVisualizeMarkerGenesElt = "#btn-visualize-marker-genes"
    btnDownloadMarkerGenesElt = "#btn-download-marker-genes"
    btnSaveMarkerGeneListElt = "#btn-save-marker-gene-list"
    markerGenesListNameElt = "#marker-gene-list-name"
    markerGenesNGenesElt = "#marker-genes-n-genes"
    markerGenesTableElt = "#marker-genes-table"
    markerGenesTableHeadTmpl = "#marker-genes-table-head-tmpl"
    markerGenesTableHeadRowElt = `${this.markerGenesTableElt} thead tr`
    markerGenesTableHeaderElts = `${this.markerGenesTableElt} thead th`
    markerGenesTableBodyTmpl = "#marker-genes-table-body-tmpl"
    markerGenesTableBodyElt = `${this.markerGenesTableElt} tbody`
    markerGenesTableBodyRowElts = `${this.markerGenesTableBodyElt} tr`
    markerGenesTableHighlightedElts = `${this.markerGenesTableElt} td.highlighted`
    markerGenesTableHeader = "#marker-genes-table-header"
    markerGenesPlotContainer = "#marker-genes-plot-c"
    markerGenesVisualizationContainer = "#marker-genes-visualization-c"
    markerGenesDotplotContainer = "#marker-genes-dotplot-c"
    markerGenesViolinContainer = "#marker-genes-violin-c"
    markerGenesManuallyEnteredElt = "#marker-genes-manually-entered"
    markerGenesSelectedCountElt = "#marker-genes-selected-count"
    markerGenesEnteredCountElt = "#marker-genes-entered-count"
    markerGenesUniqueCountElt = "#marker-genes-unique-count"

    // Clustering (edit mode)
    //clusteringToggleElt = "#toggle-clustering-edit"   // Temporary
    clusteringEditSection = "#edit-clustering-s"
    btnClusteringEditRunElt = "#btn-edit-clustering"
    // -- using resolutionElts values from the non-edit clustering
    // -- using clusterNNeighborsElt values from the non-edit clustering
    clusterTsnePlotEditElt = "#cluster-tsne-plot-edit-c"
    clusterUmapPlotEditElt = "#cluster-umap-plot-edit-c"
    clusteringEditInstructionsElt = "#analysis-clustering-edit .tool-instructions"
    groupLabelsContainer = "#group-labels-c"
    clusterGroupLabelsTmpl = "#cluster-group-labels-tmpl"
    clusterGroupLabelsTableBodyElt = "#cluster-group-labels tbody"
    clusterGroupLabelsInputElts = '#cluster-group-labels td.group-user-label input'
    clusteringMergeClustersElt = "#clustering-merge-clusters"


    // Compare Genes
    //compareGenesToggleElt = "#toggle-compare-genes"   // Temporary
    compareGenesSection = "#compare-genes-s"
    btnCompareGenesRunElt = "#btn-compare-genes-run"
    clusterOptsTmpl = "#cluster-list-tmpl"
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
    compareGenesInstructionsElt = "#analysis-compare-genes .tool-instructions"
    compareGenesResultsContainer = "#compare-genes-results-c"
    compareGenesResetableElts = `${this.compareGenesResultsContainer} .js-resetable`
    btnCompareGenesDownloadTableFElt = "#btn-compare-genes-download-table-f"
    btnCompareGenesDownloadTableRElt = "#btn-compare-genes-download-table-r"
    btnCompareGenesShowTableFElt = "#btn-compare-genes-show-table-f"
    btnCompareGenesShowTableRElt = "#btn-compare-genes-show-table-r"

    // labeled tSNE
    labeledTsneSection = "#labeled-tsne-s"
    labeledTsneContainer = "#labeled-tsne-c"
    //labeledTsneElt = "#analysis-sbs-tsne"
    btnLabeledTsneRunElt = "#btn-labeled-tsne-run"
    labeledTsnePlotContainer = "#labeled-tsne-plot-c"
    labeledTsneGeneNotFoundElt = "#labeled-tsne-gene-not-found"
    labeledTsneGeneSymbolElt = "#labeled-tsne-gene-symbol"
}

// Singleton instance of the AnalysisUI class
const UI = new AnalysisUI();