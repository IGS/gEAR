'use strict';

class UI {
    // This class is a singleton that manages the UI elements of the analysis pipeline.

    // general analysis
    emptyAnalysisOptionTmpl = "#analyses-list-empty-tmpl"
    analysisOptionTmpl = "#analyses-list-tmpl"
    analysisSelectElt = "#analysis-id"
    newAnalysisOptionElt = `${analysisSelectElt} option[data-analysis-id='0']`
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
    analysisActionContainer = "#analysis-action-c"
    analysisStatusInfoContainer = "#analysis-status-info-c"
    analysisStatusInfoElt = "#analysis-status-info"
    storedAnalysesContainer = "#stored-analyses-c"
    datasetInfoElt = "#dataset-info"
    btnSaveAnalysisElt = "#btn-save-analysis"
    btnDeleteSavedAnalysis = "#btn-delete-saved-analysis"
    btnDeleteUnsavedAnalysis = "#btn-delete-unsaved-analysis"
    btnMakePublicCopy = "#btn-make-public-copy"
    initialInstructionsElt = ".initial-instructions"

    // primary filter
    primaryFilterSectionElt = "#primary-filter-s"
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
    primaryInitialPlotContainerElt = "#primary-initial-plot-c"
    primaryInitialScatterContainer = "#primary-initial-scatter-c"
    primaryInitialViolinContainer = "#primary-initial-violin-c"
    primaryTopGenesContainer = "#primary-top-genes-c"
    primaryTopGenesPlotContainer = "#primary-top-genes-plot-c"
    primaryInitialPlotElts = ".primary-initial-plot"
    datasetInfoElt = "#dataset-info"
    datasetInfoResetableElts = `${datasetInfoElt} .js-resetable`

    // QC by mito
    qcByMitoToggleElt = "#toggle-qc-by-mito"    // Temporary
    qcByMitoSectionElt = "#qc-by-mito-s"
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
    qbmResetableElts = `${qcByMitoSectionElt} .js-resetable`

    // select variable genes
    selectVariableGenesToggleElt = "#toggle-select-variable-genes"  // Temporary
    selectVariableGenesSectionElt = "#select-variable-genes-s"
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
    pcaToggleElt = "#toggle-pca"   // Temporary
    pcaSectionElt = "#pca-s"
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
    clusteringResetableElts = `#analysis-clustering .js-resetable`

    // Marker Genes
    markerGenesToggleElt = "#toggle_marker_genes"   // Temporary
    btnMarkerGenesRunElt = "#btn-marker-genes-run"
    visualMarkerGenesSectionElt = "#visualize-marker-genes-s"
    btnVisualizeMarkerGenesElt = "#btn-visualize-marker-genes"
    btnDownloadMarkerGenesElt = "#btn-download-marker-genes"
    btnSaveMarkerGeneListElt = "#btn-save-marker-gene-list"
    markerGenesListNameElt = "#marker-gene-list-name"
    markerGenesNGenesElt = "#marker-genes-n-genes"
    markerGenesTableElt = "#marker-genes-table"
    markerGenesTableHeadTmpl = "#marker-genes-table-head-tmpl"
    markerGenesTableHeadRowElt = `${markerGenesTableElt} thead tr`
    markerGenesTableHeaderElts = `${markerGenesTableElt} thead th`
    markerGenesTableBodyTmpl = "#marker-genes-table-body-tmpl"
    markerGenesTableBodyElt = `${markerGenesTableElt} tbody`
    markerGenesTableBodyRowElts = `${markerGenesTableBodyElt} tr`
    markerGenesTableHighlightedElts = `${markerGenesTableElt} td.highlighted`
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
    clusteringToggleElt = "#toggle-clustering-edit"   // Temporary
    clusteringEditSectionElt = "#edit-clustering-s"
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
    compareGenesToggleElt = "#toggle-compare-genes"   // Temporary
    compareGenesSectionElt = "#compare-genes-s"
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
    compareGenesResetableElts = `${compareGenesResultsContainer} .js-resetable`
    btnCompareGenesDownloadTableFElt = "#btn-compare-genes-download-table-f"
    btnCompareGenesDownloadTableRElt = "#btn-compare-genes-download-table-r"
    btnCompareGenesShowTableFElt = "#btn-compare-genes-show-table-f"
    btnCompareGenesShowTableRElt = "#btn-compare-genes-show-table-r"

    // labeled tSNE
    labeledTsneElt = "#analysis-sbs-tsne"
    btnLabeledTsneRunElt = "#btn-labeled-tsne-run"
    labeledTsnePlotContainer = "#labeled-tsne-plot-c"
    labeledTsneGeneNotFoundElt = "#labeled-tsne-gene-not-found"
    labeledTsneGeneSymbolElt = "#labeled-tsne-gene-symbol"
}
