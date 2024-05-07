"use strict";

/*
  Classes representing overall analysis (pipeline) elements and their child classes.
*/

// requires common.js

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
    }

    /**
     * Adds a gene of interest to the set of genes.
     *
     * @param {string} geneSymbol - The symbol of the gene to add.
     */
    addGeneOfInterest(geneSymbol) {
        geneSymbol = geneSymbol.trim();
        this.genesOfInterest.add(geneSymbol);

        if (this.genesOfInterest.length == 1) {
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
                throw new Error("Error copying analysis: " + error);
            }

            thisAnalysis.type = 'user_unsaved';
            thisAnalysis.id = newAnalysisId;
            thisAnalysis.userSessionId = CURRENT_USER.session_id;

            // $("#analysis_action_c").show();
            // $("#analysis_status_info_c").hide();

            thisAnalysis.getSavedAnalysesList(thisAnalysis.datasetId, newAnalysisId);

            if (callback) {
                callback(opts);
            }

        } catch (error) {
            createToast(`Error creating sandbox copy of the analysis: ${error.message}`);
        }
    }

}

class AnalysisStepPrimaryFilter {
    constructor({}) {}

}
class AnalysisStepQCByMito {
    constructor({}) {}

}
class AnalysisStepSelectVariableGenes {
    constructor({}) {}

}
class AnalysisStepPCA {
    constructor({}) {}

}
class AnalysisSteptSNE {
    constructor({}) {}

}
class AnalysisStepClustering {
    constructor({}) {}
}
class AnalysisStepMarkerGenes {
    constructor({}) {
        this.visualizeMarkerGenesBtnElt = document.getElementById("btn-visualize-marker-genes");
    }

}
class AnalysisStepCompareGenes {
    constructor({}) {}

}
