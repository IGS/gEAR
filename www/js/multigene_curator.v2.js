const isMultigene = 1;

let geneSelect;

class DashHandler extends PlotHandler {
    constructor(plotType) {
        super();
        this.plotType = plotType;
    }

    classElt2Prop = {
        //pass
    };
    configProp2ClassElt = Object.fromEntries(Object.entries(this.classElt2Prop).map(([key, value]) => [value, key]));

    plotConfig = {};  // Plot config that is passed to API

    cloneDisplay() {
        //pass
    }

    async createPlot() {
        //pass
    }

    async loadPlotHtml() {
        //pass
    }

    populatePlotConfig() {
        //pass
    }

    setupPlotSpecificEvents() {
        //pass
    }

}

const geneCartTree = new GeneCartTree({
    element: document.getElementById("genecart_tree")
    , searchElement: document.getElementById("genecart_query")
    , selectCallback: (async (e) => {
        if (e.node.type !== "genecart") {
            return;
        }

        // Get gene symbols from gene cart
        const geneCartId = e.node.data.orig_id;
        const geneCartMembers = await fetchGeneCartMwembers(sessionId, geneCartId);
        const geneCartSymbols = geneCartMembers.map((item) => item.label);

        // Normalize gene symbols to lowercase
        const geneSelectSymbols = geneSelect.data.map((opt) => opt.value);
        const geneCartSymbolsLowerCase = geneCartSymbols.map((x) => x.toLowerCase());

        const geneSelectedOptions = geneSelect.selectedOptions.map((opt) => opt.data.value);

        // Get genes from gene cart that are present in dataset's genes.  Preserve casing of dataset's genes.
        const geneCartIntersection = geneSelectSymbols.filter((x) => geneCartSymbolsLowerCase.includes(x.toLowerCase()));
        // Add in already selected genes (union)
        const geneSelectIntersection = [...new Set(geneCartIntersection.concat(geneSelectedOptions))];

        // change all options to be unselected
        const origSelect = document.getElementById("gene_select");
        for (const opt of origSelect.options) {
            opt.removeAttribute("selected");
        }

        // Assign intersection genes to geneSelect "selected" options
        for (const gene of geneSelectIntersection) {
            const opt = origSelect.querySelector(`option[value="${gene}"]`);
            try {
                opt.setAttribute("selected", "selected");
            } catch (error) {
                // sanity check
                const msg = `Could not add gene ${gene} to gene select.`;
                console.warn(msg);
            }
        }

        geneSelect.update();
        trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to load tags
    })
});

const appendGeneTagButton = (geneTagElt) => {
    // Add delete button
    const deleteBtnElt = document.createElement("button");
    deleteBtnElt.classList.add("delete", "is-small");
    geneTagElt.appendChild(deleteBtnElt);
    deleteBtnElt.addEventListener("click", (event) => {
        // Remove gene from geneSelect
        const gene = event.target.parentNode.textContent;
        const geneSelectElt = document.getElementById("gene_select");
        geneSelectElt.querySelector(`option[value="${gene}"]`).removeAttribute("selected");

        geneSelect.update();
        trigger(document.getElementById("gene_select"), "change"); // triggers chooseGene() to load tags
    });
}

const clearGenes = (event) => {
    document.getElementById("clear_genes_btn").classList.add("is-loading");
    geneSelect.clear();
    document.getElementById("clear_genes_btn").classList.remove("is-loading");
}

const curatorSpecifcChooseGene = (event) => {
    // Triggered when a gene is selected

    // Delete existing tags
    const geneTagsElt = document.getElementById("gene_tags");
    geneTagsElt.replaceChildren();

    if (!geneSelect.selectedOptions.length) return;   // Do not trigger after initial population

    // Update list of gene tags
    const sortedGenes = geneSelect.selectedOptions.map((opt) => opt.data.value).sort();
    for (const opt in sortedGenes) {
        const geneTagElt = document.createElement("span");
        geneTagElt.classList.add("tag", "is-primary");
        geneTagElt.textContent = sortedGenes[opt];
        appendGeneTagButton(geneTagElt);
        geneTagsElt.appendChild(geneTagElt);
    }

    document.getElementById("gene_tags_c").style.display = "";
    if (!geneSelect.selectedOptions.length) {
        document.getElementById("gene_tags_c").style.display = "none";
    }

    // Cannot plot if 2+ genes are not selected
    if (geneSelect.selectedOptions.length < 2) {
        document.getElementById("gene_s_failed").style.display = "";
        document.getElementById("gene_s_success").style.display = "none";
        for (const plotBtn of document.getElementsByClassName("js-plot-btn")) {
            plotBtn.disabled = true;
        }
        document.getElementById("continue_to_plot_options").style.display = "none";
        return;
    }

    document.getElementById("gene_s_failed").style.display = "none";
    document.getElementById("gene_s_success").style.display = "";

    // Force validation check to see if plot button should be enabled
    //trigger(document.querySelector(".js-plot-req"), "change");

    document.getElementById("continue_to_plot_options").style.display = "";

}

const curatorSpecifcCreatePlot = async (plotType) => {

}

const curatorSpecifcDatasetTreeCallback = () => {
    // Creates gene select instance that allows for multiple selection
    geneSelect = createGeneSelectInstance(geneSelect, true)
}

const curatorSpecificOnLoad = async () => {
    // Load gene carts
    await loadGeneCarts();
}

const curatorSpecificUpdateGeneOptions = async (geneSymbols) => {
    //pass
}

const fetchAvailablePlotTypes = async (user_id, session_id, dataset_id, analysis_id) => {
    // Plot types will depend on the number of comparabie categorical conditions
    // Volcano plots must have at least two conditions
    // Quadrant plots must have at least three conditions

    const payload = {user_id, session_id, dataset_id, analysis_id};
    try {
        const {data} = await axios.post(`/api/h5ad/${dataset_id}/mg_availableDisplayTypes`, payload);
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch compatible plot types for this dataset. Please contact the gEAR team.";
        createToast(msg);
        throw new Error(msg);
    }
}

/* Fetch gene collections */
const fetchGeneCarts = async (session_id) => {
    const payload = {session_id};
    try {
        const {data} = await axios.post(`/cgi/get_user_gene_carts.cgi`, convertToFormData(payload));
        return data;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene collections. You can still enter genes manually.";
        createToast(msg);
        throw new Error(msg);
    }
}

/* Fetch gene collection members */
const fetchGeneCartMwembers = async (session_id, geneCartId) => {
    const payload = { session_id, gene_cart_id: geneCartId };
    try {
        const {data} = await axios.post(`/cgi/get_gene_cart_members.cgi`, convertToFormData(payload));
        const {gene_symbols, success} = data;
        if (!success) {
            throw new Error("Could not fetch gene collection members. You can still enter genes manually.");
        }
        return gene_symbols;
    } catch (error) {
        logErrorInConsole(error);
        const msg = "Could not fetch gene collection members. You can still enter genes manually.";
        createToast(msg);
        throw new Error(msg);
    }
}

const includePlotParamOptions = async () => {

}

/* Transform and load gene collection data into a "tree" format */
const loadGeneCarts = async () => {
    try {
        const geneCartData = await fetchGeneCarts(sessionId);
        const carts = {};
        const cartTypes = ['domain', 'user', 'group', 'shared', 'public'];
        let cartsFound = false;

        // Loop through the different types of gene collections and add them to the carts object
        for (const ctype of cartTypes) {
            carts[ctype] = [];

            if (geneCartData[`${ctype}_carts`].length > 0) {
                cartsFound = true;

                for (const item of geneCartData[`${ctype}_carts`]) {
                    carts[ctype].push({value: item.id, text: item.label });
                };
            }
        }

        geneCartTree.domainGeneCarts = carts.domain;
        geneCartTree.userGeneCarts = carts.user;
        geneCartTree.groupGeneCarts = carts.group;
        geneCartTree.sharedGeneCarts = carts.shared;
        geneCartTree.publicGeneCarts = carts.public;
        geneCartTree.generateTree();
        /*if (!cartsFound ) {
            // ? Put some warning if carts not found
            $('#gene_cart_container').show();
        }*/

    } catch (error) {
        document.getElementById("gene_s_failed").style.display = "";
    }

}

document.getElementById("clear_genes_btn").addEventListener("click", clearGenes);