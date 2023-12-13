"use strict";

let firstSearch = true;
let isAddFormOpen = false;
const resultsPerPage = 20;

// TODO - Add transformation code for quick gene collection transformations

// Floating UI function alias. See https://floating-ui.com/docs/getting-started#umd
// Subject to change if we ever need these common names for other things.
const computePosition = window.FloatingUIDOM.computePosition;
const flip = window.FloatingUIDOM.flip;
const shift = window.FloatingUIDOM.shift;
const offset = window.FloatingUIDOM.offset;
const arrow = window.FloatingUIDOM.arrow;

/**
 * Adds event listeners for various actions related to gene collections.
 * @function
 * @returns {void}
 */
const addGeneCollectionEventListeners = () => {

    // Unweighted gene lists
    setupGeneListToggle("js-gc-unweighted-gene-list-toggle", './cgi/get_unweighted_gene_cart_preview.cgi', createUnweightedGeneCollectionTable);

    // Weighted gene lists
    setupGeneListToggle("js-gc-weighted-gene-list-toggle", './cgi/get_weighted_gene_cart_preview.cgi', createWeightedGeneCollectionPreview);

    // Labeled gene lists
    // NOT IMPLEMENTED YET

    // Expand and collapse gene collection view
    for (const classElt of document.getElementsByClassName("js-expand-box")) {
        classElt.addEventListener("click", (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            const selector = `#result_gc_id_${gcId} .js-expandable-view`;
            const expandableViewElts = document.querySelectorAll(selector);
            for (const classElt of expandableViewElts) {
                classElt.classList.toggle("is-hidden");
            }

            // Toggle the icon
            if (e.currentTarget.innerHTML === '<i class="mdi mdi-arrow-expand"></i>') {
                e.currentTarget.innerHTML = '<i class="mdi mdi-arrow-collapse"></i>';
                return;
            }
            e.currentTarget.innerHTML = '<i class="mdi mdi-arrow-expand"></i>';

        });
    }

    // Reformats the <ul> containing the gene symbols into a text file with one gene per row
    for (const classElt of document.getElementsByClassName("js-download-gc")) {
        classElt.addEventListener("click", async (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            const gctype = e.currentTarget.dataset.gcType;
            const shareId = e.currentTarget.dataset.gcShareId;

            if (gctype == "unweighted-list") {
                const gene_symbols = await fetchGeneCartMembers(gcId);
                const fileContents = gene_symbols.map(gene => gene.label).join("\n");

                const element = document.createElement("a");
                element.setAttribute(
                    "href",
                    `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(fileContents)}`
                );

                element.setAttribute("download", `gene_cart.${shareId}.tsv`);
                element.style.display = "none";
                document.body.appendChild(element);
                element.click();
                document.body.removeChild(element);
                return;
            }

            if (gctype == "weighted-list") {
                // Download the source file (returns a content-disposition attachment header)

                const {data} = await axios.post("./cgi/download_weighted_gene_cart.cgi", convertToFormData({
                    'share_id': shareId
                }));

                const element = document.createElement("a");
                element.setAttribute(
                    "href",
                    `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(data)}`
                );

                element.setAttribute("download", `gene_cart.${shareId}.tsv`);
                element.style.display = "none";
                document.body.appendChild(element);
                element.click();
                document.body.removeChild(element);
                return;
            }
            if (gctype == "labeled-list") {
                // Retrieve the cart from the server and put it in a text file
                throw Error("Not implemented yet");

            } else {
                throw Error(`Invalid gene collection type: ${gctype}`);
            }

        });
    }

    for (const classElt of document.getElementsByClassName("js-share-gc")) {
        classElt.addEventListener("click", (e) => {

            const shareId = e.currentTarget.value;
            const currentUrl = window.location.href;
            const currentPage = currentUrl.lastIndexOf("gene_collection_manager.html");
            const shareUrl = `${currentUrl.substring(0, currentPage)}p?c=${shareId}`;
            const gcId = e.currentTarget.dataset.gcId;

            showGcActionNote(gcId, shareUrl);

        });
    }

    // Cancel button for editing a gene collection
    for (const classElt of document.getElementsByClassName("js-edit-gc-cancel")) {
        classElt.addEventListener("click", (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            const selectorBase = `#result_gc_id_${gcId}`;

            // Show editable versions where there are some and hide the display versions
            for (const classElt of document.querySelectorAll(`${selectorBase} .js-editable-version`)) {
                classElt.classList.add("is-hidden");
            };
            for (const classElt of document.querySelectorAll(`${selectorBase} .js-display-version`)) {
                classElt.classList.remove("is-hidden");
            };

            // Reset any unsaved/edited values
            const visibility = document.querySelector(`${selectorBase}_editable_visibility`).dataset.originalVal;
            document.querySelector(`${selectorBase}_editable_visibility`).value = visibility;

            const title = document.querySelector(`${selectorBase}_editable_title`).dataset.originalVal;
            document.querySelector(`${selectorBase}_editable_title`).value = title;

            const orgId = document.querySelector(`${selectorBase}_editable_organism_id`).dataset.originalVal;
            document.querySelector(`${selectorBase}_editable_organism_id`).value = orgId;

            document.querySelector(`${selectorBase} .js-action-links`).classList.remove("is-hidden");

        });
    }

    // Save button for editing a gene collection
    for (const classElt of document.getElementsByClassName("js-edit-gc-save")) {
        classElt.addEventListener("click", async (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            const selectorBase = `#result_gc_id_${gcId}`;
            const newVisibility = document.querySelector(`${selectorBase}_editable_visibility`).value;
            const newTitle = document.querySelector(`${selectorBase}_editable_title`).value;
            const newOrgId = document.querySelector(`${selectorBase}_editable_organism_id`).value;
            const newLdesc = document.querySelector(`${selectorBase}_editable_ldesc`).value;

            try {
                const {data} = await axios.post('./cgi/save_genecart_changes.cgi', convertToFormData({
                    'session_id': CURRENT_USER.session_id,
                    'gc_id': gcId,
                    'visibility': newVisibility,
                    'title': newTitle,
                    'organism_id': newOrgId,
                    'ldesc': newLdesc || ""
                }));

                // Update the UI for the new values
                document.querySelector(`${selectorBase}_editable_visibility`).dataset.isPublic = newVisibility;
                if (newVisibility) {
                    document.querySelector(`${selectorBase}_display_visibility`).textContent = "Public gene collection";
                    document.querySelector(`${selectorBase}_display_visibility`).classList.remove("is-danger");
                    document.querySelector(`${selectorBase}_display_visibility`).classList.add("is-light");
                } else {
                    document.querySelector(`${selectorBase}_display_visibility`).textContent = "Private gene collection";
                    document.querySelector(`${selectorBase}_display_visibility`).classList.remove("is-light");
                    document.querySelector(`${selectorBase}_display_visibility`).classList.add("is-danger");
                }

                document.querySelector(`${selectorBase}_editable_title`).dataset.originalVal = newTitle;
                document.querySelector(`${selectorBase}_display_title`).textContent = newTitle;

                document.querySelector(`${selectorBase}_editable_ldesc`).dataset.originalVal = newLdesc;
                document.querySelector(`${selectorBase}_display_ldesc`).textContent = newLdesc || "No description entered";;

                document.querySelector(`${selectorBase}_display_organism`).textContent =
                    document.querySelector(`${selectorBase}_editable_organism_id > option[value='${newOrgId}']`).textContent;
                document.querySelector(`${selectorBase}_editable_organism_id`).dataset.originalVal = newOrgId;

                // Put interface back to view mode.
                toggleEditableMode(true, selectorBase);

                createToast("Gene collection changes saved", "is-success");

            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to save gene collection changes");
            } finally {
                document.querySelector(`${selectorBase} .js-action-links`).classList.remove("is-hidden");
            }
        });
    }

    // Toggle editable mode when edit button is clicked for a gene collection
    for (const classElt of document.getElementsByClassName("js-edit-gc")) {
        classElt.addEventListener("click", async (e) => {

            const gcId = e.currentTarget.dataset.gcId;
            const selectorBase = `#result_gc_id_${gcId}`;

            // copy the organism selection list for this row
            const editableOrganismIdElt = document.querySelector(`${selectorBase}_editable_organism_id`);
            editableOrganismIdElt.innerHTML = document.getElementById("new_collection_organism_id").innerHTML;

            // Remove the "select an organism" option
            editableOrganismIdElt.removeChild(editableOrganismIdElt.firstChild);

            // set the current value as selected
            editableOrganismIdElt.value = editableOrganismIdElt.dataset.originalVal;

            const editableVisibilityElt = document.querySelector(`${selectorBase}_editable_visibility`);

            const isPublic = parseBool(editableVisibilityElt.dataset.isPublic);

            editableVisibilityElt.checked = isPublic;
            editableVisibilityElt.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";

            // Show editable versions where there are some and hide the display versions
            toggleEditableMode(false, selectorBase);

            // Make sure the view is expanded
            const expandableViewElt = document.querySelector(`${selectorBase} .js-expandable-view`);
            if (expandableViewElt.classList.contains('is-hidden')) {
                document.querySelector(`${selectorBase} span.js-expand-box`).click();
            }

            document.querySelector(`${selectorBase} .js-action-links`).classList.add("is-hidden");

        });
    }

    // Redirect to gene expression search
    for (const classElt of document.getElementsByClassName("js-view-gc")) {
        classElt.addEventListener("click", (e) => {
            window.location = `./p?c=${e.currentTarget.value}`;
        });
    }
}

/**
 * Adds gene collection display elements to DOM.
 * @param {string} geneCollectionId - The ID of the gene collection.
 * @param {string} gctype - The type of gene collection.
 * @param {string} shareId - The ID of the share.
 * @param {number} geneCount - The number of genes.
 */
const addPreviewGenesToGeneCollection = (geneCollectionId, gctype, shareId, geneCount) => {
    const geneCollectionPreviewGenesContainer = document.getElementById(`${geneCollectionId}_preview_genes_container`);
    const button = document.createElement('button');
    button.className = 'button is-small is-dark is-outlined';
    button.title = 'Show gene collection';
    button.dataset.gcId = geneCollectionId;
    button.dataset.gcShareId = shareId;
    button.innerHTML = `<span class="icon is-small"><i class="mdi mdi-format-list-bulleted"></i></span>`;

    if (gctype === "weighted-list") {
        button.classList.add("js-gc-weighted-gene-list-toggle");
        const elt = document.createElement("span");
        elt.id = `btn_gc_${geneCollectionId}_text`;
        elt.textContent = `Info`;
        elt.dataset.offState = elt.textContent
        button.append(elt);
    } else if (gctype === "unweighted-list") {
        button.classList.add("js-gc-unweighted-gene-list-toggle");
        const elt = document.createElement("span");
        elt.id = `btn_gc_${geneCollectionId}_text`;
        elt.textContent = `${geneCount} genes`;
        elt.dataset.offState = elt.textContent
        button.append(elt);
    } else if (gctype === "labeled-list") {
        // Not implemented yet
        button.classList.add("js-gc-labeled-gene-list-toggle");
    } else {
        throw Error(`Invalid gene collection type: ${gctype}`);
    }

    geneCollectionPreviewGenesContainer.append(button);

}

/**
 * Adds public/private visibility information to a gene collection display container in the DOM.
 * @param {string} geneCollectionId - The ID of the gene collection display container.
 * @param {boolean} isPublic - A boolean indicating whether the gene collection is public or private.
 * @returns {void}
 */
const addVisibilityInfoToGeneCollection = (geneCollectionId, isPublic) => {

    // add gene collection public/private info to DOM
    const geneCollectionDisplayContainer = document.getElementById(`${geneCollectionId}_display_visibility`);
    const geneCollectionDisplaySpan = document.createElement("span");
    geneCollectionDisplaySpan.classList.add("tag");
    geneCollectionDisplaySpan.id = `result_gc_id_${geneCollectionId}_display_visibility`;

    if (isPublic) {
        geneCollectionDisplaySpan.classList.add("is-primary", "is-light");
        geneCollectionDisplaySpan.textContent = "Public gene collection";
    } else {
        geneCollectionDisplaySpan.classList.add("is-danger");
        geneCollectionDisplaySpan.textContent = "Private gene collection";
    }
    geneCollectionDisplayContainer.appendChild(geneCollectionDisplaySpan);

    // Toggle switch (public is checked, private is unchecked)
    const visibilitySwitch = document.getElementById(`result_gc_id_${geneCollectionId}_editable_visibility`);

    visibilitySwitch.addEventListener("change", (e) => {
        const isPublic = e.currentTarget.checked;
        e.currentTarget.dataset.isPublic = isPublic;
        e.currentTarget.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";
    });

}

/**
 * Creates a tooltip for a reference element.
 *
 * @param {HTMLElement} referenceElement - The reference element to attach the tooltip to.
 * @param {HTMLElement} tooltip - The tooltip element.
 */
const applyTooltip = (referenceElement, tooltip) => {

    const hideTooltip = (event) => {
        tooltip.classList.add("is-hidden");
    }

    const showTooltip = (event) => {
        // Compute position
        computePosition(event.currentTarget, tooltip, {
            placement: 'top', // Change this to your preferred placement
            middleware: [
                flip(), // flip to bottom if there is not enough space on top
                shift(), // shift the popover to the right if there is not enough space on the left
                offset(5), // offset relative to the button
            ]
        }).then(({ x, y }) => {
            // Position the popover
            Object.assign(tooltip.style, {
                left: `${x}px`,
                top: `${y}px`,
            });
        });

        tooltip.classList.remove("is-hidden");
    }

    [
        ['mouseenter', showTooltip],
        ['mouseleave', hideTooltip],
        ['focus', showTooltip],
        ['blur', hideTooltip],
    ].forEach(([event, listener]) => {

        referenceElement.addEventListener(event, listener);
    });

}

/**
 * Builds a comma-separated string of selected database values for a given group name.
 * @param {string} groupName - The ID of the group to retrieve selected values from.
 * @returns {string} A comma-separated string of selected database values.
 */
const buildFilterString = (groupName) => {
    const selected = document.querySelectorAll(`#${groupName} li.js-selected:not(.js-all-selector)`);
    const dbvals = [];

    for (const li of selected) {
        dbvals.push(li.dataset.dbval);
    }

    return dbvals.join(",");
}

/**
 * Creates a tooltip element and appends it to the body.
 * @param {HTMLElement} referenceElement - The reference element to which the tooltip is associated.
 * @returns {HTMLElement} The created tooltip element.
 */
const createActionTooltips = (referenceElement) => {
    // Create tooltip element
    const tooltip = document.createElement('div');
    tooltip.className = 'tooltip';
    tooltip.textContent = referenceElement.dataset.tooltipContent;
    tooltip.classList.add("has-background-dark", "has-text-white", "is-hidden");

    // Append tooltip to body
    document.body.appendChild(tooltip);
    return tooltip;
}

/**
 * Creates a confirmation popover for deleting a gene collection.
 */
const createDeleteConfirmationPopover = () => {
    const deleteButtons = document.getElementsByClassName("js-delete-gc");
    for (const button of deleteButtons) {
        button.addEventListener('click', (e) => {
            // remove existing popovers
            const existingPopover = document.getElementById('delete_gc_popover');
            if (existingPopover) {
                existingPopover.remove();
            }

            // Create popover content
            const popoverContent = document.createElement('article');
            popoverContent.id = 'delete_gc_popover';
            popoverContent.classList.add("message", "is-danger");
            popoverContent.setAttribute("role", "tooltip");
            popoverContent.innerHTML = `
                <div class='message-header'>
                    <p>Remove collection</p>
                </div>
                <div class='message-body'>
                    <p>Are you sure you want to delete this gene collection?</p>
                    <div class='field is-grouped' style='width:250px'>
                        <p class="control">
                            <button id='confirm_gc_delete' class='button is-danger'>Delete</button>
                        </p>
                        <p class="control">
                            <button id='cancel_gc_delete' class='button' value='cancel_delete'>Cancel</button>
                        </p>
                    </div>
                </div>
                <div id="arrow"></div>
            `;

            // append element to DOM to get its dimensions
            document.body.appendChild(popoverContent);

            const arrowElement = document.getElementById('arrow');

            // Create popover (help from https://floating-ui.com/docs/tutorial)
            computePosition(button, popoverContent, {
                placement: 'top',
                middleware: [
                    flip(), // flip to bottom if there is not enough space on top
                    shift(), // shift the popover to the right if there is not enough space on the left
                    offset(5), // offset relative to the button
                    arrow({ element: arrowElement }) // add an arrow pointing to the button
                ],
            }).then(({ x, y, placement, middlewareData }) => {
                // Position the popover
                Object.assign(popoverContent.style, {
                    left: `${x}px`,
                    top: `${y}px`,
                });
                // Accessing the data
                const { x: arrowX, y: arrowY } = middlewareData.arrow;

                // Position the arrow relative to the popover
                const staticSide = {
                    top: 'bottom',
                    right: 'left',
                    bottom: 'top',
                    left: 'right',
                }[placement.split('-')[0]];

                // Set the arrow position
                Object.assign(arrowElement.style, {
                    left: arrowX != null ? `${arrowX}px` : '',
                    top: arrowY != null ? `${arrowY}px` : '',
                    right: '',
                    bottom: '',
                    [staticSide]: '-4px',
                });
            });

            // Show popover
            document.body.appendChild(popoverContent);

            // Store the gene collection ID to delete
            const gcIdToDelete = e.currentTarget.value;

            // Add event listener to cancel button
            document.getElementById('cancel_gc_delete').addEventListener('click', () => {
                popoverContent.remove();
            });

            // Add event listener to confirm button
            document.getElementById('confirm_gc_delete').addEventListener('click', async () => {

                try {
                    const {data} = await axios.post('./cgi/remove_gene_cart.cgi', convertToFormData({
                        'session_id': CURRENT_USER.session_id,
                        'gene_cart_id': gcIdToDelete
                    }));

                    if (data['success'] == 1) {
                        const resultElement = document.getElementById(`result_gc_id_${gcIdToDelete}`);
                        resultElement.style.transition = 'opacity 1s';
                        resultElement.style.opacity = 0;
                        resultElement.remove();

                        createToast("Gene collection deleted", "is-success");

                        // This can affect page counts, so we need to re-run the search
                        submitSearch();

                    } else {
                        throw new Error(data['error']);
                    }
                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to delete gene collection");
                } finally {
                    popoverContent.remove();
                }
            });
        });
    }
}

/**
 * Creates a pagination button element.
 *
 * @param {number} page - The page number to display on the button.
 * @param {string|null} icon - The icon class name to display on the button. If null, the page number will be displayed instead.
 * @param {Function} clickHandler - The click event handler function for the button.
 * @returns {HTMLLIElement} - The created pagination button element.
 */
const createPaginationButton = (page, icon = null, clickHandler) => {
    const li = document.createElement("li");
    const button = document.createElement("button");
    button.className = "button is-small is-outlined is-dark pagination-link";
    if (icon) {
        button.innerHTML = `<i class="mdi mdi-chevron-${icon}"></i>`;
    } else {
        button.textContent = page;
    }
    button.addEventListener("click", clickHandler);
    li.appendChild(button);
    return li;
}

/**
 * Creates a pagination ellipsis element.
 * @returns {HTMLLIElement} The created list item element containing the pagination ellipsis.
 */
const createPaginationEllipsis = () => {
    const li = document.createElement("li");
    const span = document.createElement("span");
    span.className = "pagination-ellipsis";
    span.textContent = "…";
    li.appendChild(span);
    return li;
}

/**
 * Creates an unweighted gene collection table.
 *
 * @param {HTMLElement} infoContainer - The container element where the table will be appended.
 * @param {Object} data - The gene information data.
 * @param {Object} data.gene_info - The gene information object. The keys are Ensembl IDs and the values are gene symbols and descriptions.
 */
const createUnweightedGeneCollectionTable = (infoContainer, data) => {
    const table = document.createElement("table");
    table.classList.add("table", "is-narrow", "is-hoverable", "is-fullwidth");
    infoContainer.appendChild(table);

    const thead = document.createElement("thead");
    table.appendChild(thead);

    const theadTr = document.createElement("tr");
    thead.appendChild(theadTr);

    const thGeneSymbol = document.createElement("th");
    thGeneSymbol.textContent = "Gene symbol";
    theadTr.appendChild(thGeneSymbol);

    const thEnsemblId = document.createElement("th");
    thEnsemblId.textContent = "Ensembl ID";
    theadTr.appendChild(thEnsemblId);

    const thGeneDescription = document.createElement("th");
    thGeneDescription.textContent = "Gene description";
    theadTr.appendChild(thGeneDescription);

    const tbody = document.createElement("tbody");
    table.appendChild(tbody);

    // sort the gene info by gene symbol
    const sortedGeneInfo = Object.entries(data.gene_info).sort((a, b) => {
        if (a[1].gene_symbol < b[1].gene_symbol) {
            return -1;
        }
        if (a[1].gene_symbol > b[1].gene_symbol) {
            return 1;
        }
        return 0;
    });

    for (const [ensemblId, gene] of sortedGeneInfo) {

        const tr = document.createElement("tr");
        tbody.appendChild(tr);

        const tdGeneSymbol = document.createElement("td");
        tdGeneSymbol.textContent = gene.gene_symbol;
        tr.appendChild(tdGeneSymbol);

        const tdEnsemblId = document.createElement("td");
        tdEnsemblId.textContent = ensemblId;
        tr.appendChild(tdEnsemblId);

        const tdGeneDescription = document.createElement("td");
        tdGeneDescription.textContent = gene.product;
        tr.appendChild(tdGeneDescription);
    }
}

/**
 * Creates a preview of a weighted gene collection.
 *
 * @param {HTMLElement} infoContainer - The container element where the gene collection preview will be appended.
 * @param {Object} data - The data object containing information about the gene collection.
 * @param {number} data.num_genes - The number of genes in the collection.
 * @param {string[]} data.weights - The names of the weights in the gene collection.
 */
const createWeightedGeneCollectionPreview = (infoContainer, data) => {
    // List the number of genes as well as the names of weights in the gene collection
    const geneCount = data.num_genes;
    const weights = data.weights;

    const geneCountElt = document.createElement("p");
    geneCountElt.innerHTML = `<span class="has-text-weight-semibold">Genes:</span> ${geneCount}`;
    infoContainer.appendChild(geneCountElt);
    const numWeightsElt = document.createElement("p");
    numWeightsElt.innerHTML = `<span class="has-text-weight-semibold">Weights:</span> ${weights.length}`;
    infoContainer.appendChild(numWeightsElt);

    const weightsElt = document.createElement("p");
    weightsElt.innerHTML = `<span class="has-text-weight-semibold">Weight names:</span> ${weights.join(", ")}`;
    infoContainer.appendChild(weightsElt);
}

/**
 * Fetches the members of a gene cart.
 * @param {string} geneCartId - The ID of the gene cart.
 * @returns {Promise<Array<string>>} - A promise that resolves to an array of gene symbols.
 * @throws {Error} - If the gene collection members cannot be fetched.
 */
const fetchGeneCartMembers = async (geneCartId) => {
    try {
        const {gene_symbols, success} = await apiCallsMixin.fetchGeneCartMembers(geneCartId);
        if (!success) {
            throw new Error("Could not fetch gene collection members.");
        }
        return gene_symbols;
    } catch (error) {
        logErrorInConsole(error);
        createToast(error.message);
        throw error
    }
}

// Callbacks after attempting to save a gene collection
const geneCollectionFailure = (gc, message) => {
    logErrorInConsole(message);
    createToast("Failed to save gene collection");
}

const geneCollectionSaved = (gc) => {
    document.getElementById("create_new_gene_collection").click(); // resets form also
    createToast("Gene collection saved", "is-success");
    submitSearch();
}

/**
 * Loads the list of organisms from the server and populates the organism choices and new cart organism ID select elements.
 * @function
 * @returns {void}
 */
const loadOrganismList = async () => {
    try {
        const {data} = await axios.get('./cgi/get_organism_list.cgi');
        const organismChoices = document.getElementById("organism_choices");    // <ul> element
        for (const organism of data.organisms) {
            const li = document.createElement("li");
            li.dataset.dbval = organism.id;
            li.textContent = organism.label;
            organismChoices.appendChild(li);
        }
        const newCollectionOrganismSelect = document.getElementById("new_collection_organism_id");    // <select> element
        newCollectionOrganismSelect.innerHTML = "";

        for (const organism of data.organisms) {
            const option = document.createElement("option");
            option.value = organism.id;
            option.textContent = organism.label;
            newCollectionOrganismSelect.appendChild(option);
        }
        // Add default "select an organism" option
        const defaultOption = document.createElement("option");
        defaultOption.value = "";
        defaultOption.textContent = "Select an organism";
        defaultOption.selected = true;
        newCollectionOrganismSelect.prepend(defaultOption);

    } catch (error) {
        logErrorInConsole(error);
        createToast("Failed to load organism list");
    }
}


/**
 * Parses a string representation of a boolean value and returns the corresponding boolean value.
 * @param {string} boolStr - The string representation of the boolean value.
 * @returns {boolean} - The parsed boolean value.
 */
const parseBool = (boolStr) => {
    if ( boolStr === 'true' ) {
        return true;
    }
    return false;
}


/**
 * Processes search results and updates the DOM with the results view.
 * @param {Object} data - The search results data.
 * @param {string} resultLabel - The label to display for the search results.
 */
const processSearchResults = (data) => {

    const resultsContainer = document.getElementById("results_container");

    // If there are no results, display a message
    if (data.gene_carts.length === 0) {
        const noResultsMessage = document.createElement("p");
        noResultsMessage.className = "has-text-centered";
        noResultsMessage.textContent = "No results found.";
        resultsContainer.appendChild(noResultsMessage);
        return;
    }

    const template = document.getElementById("results_view");

    // data.gene_carts is a list of JSON strings
    for (const gcString of data.gene_carts) {
        const gc = JSON.parse(gcString);
        const geneCollectionId = gc.id;
        const gctype = gc.gctype;
        const label = gc.label;
        const longDesc = gc.ldesc || "";
        const shareId = gc.share_id;
        const isPublic = Boolean(gc.is_public);
        const dateAdded = new Date(gc.date_added).toDateString();

        const organismId = gc.organism_id;
        const geneCount = gc.gene_count;
        const userName = gc.user_name;
        const organism = gc.organism;
        const isOwner = gc.is_owner;

        // Clone the template
        const resultsView = template.content.cloneNode(true)

        // Set properties for multiple elements
        setElementProperties(resultsView, ".js-gc-list-element", { id: `result_gc_id_${geneCollectionId}`, dataset: { gcId: geneCollectionId } });
        // title section
        setElementProperties(resultsView, ".js-display-title p", { id: `result_gc_id_${geneCollectionId}_display_title`, textContent: label });
        setElementProperties(resultsView, ".js-editable-title input", { id: `result_gc_id_${geneCollectionId}_editable_title`, dataset: { originalVal: label }, value: label });
        setElementProperties(resultsView, ".js-expand-box", { dataset: { gcId: geneCollectionId } });
        // visibility/other metadata section
        setElementProperties(resultsView, ".js-display-visibility", { id: `${geneCollectionId}_display_visibility` });
        setElementProperties(resultsView, ".js-editable-visibility input", { id: `result_gc_id_${geneCollectionId}_editable_visibility`, checked: isPublic, dataset: { isPublic } });
        setElementProperties(resultsView, ".js-editable-visibility label", { htmlFor: `result_gc_id_${geneCollectionId}_editable_visibility`, textContent: isPublic ? "Public" : "Private" });
        // organism section
        setElementProperties(resultsView, ".js-display-organism span:last-of-type", { id: `result_gc_id_${geneCollectionId}_display_organism`, textContent: organism });
        setElementProperties(resultsView, ".js-editable-organism select", { id: `result_gc_id_${geneCollectionId}_editable_organism_id`, dataset: { originalVal: organismId }, value: organismId });
        setElementProperties(resultsView, ".js-editable-organism label", { htmlFor: `result_gc_id_${geneCollectionId}_editable_organism_id` });
        // owner section
        setElementProperties(resultsView, ".js-display-owner span:last-of-type", { textContent: userName });
        setElementProperties(resultsView, ".js-editable-owner input", { value: userName });
        // date added section
        setElementProperties(resultsView, ".js-display-date-added span:last-of-type", { textContent: dateAdded });
        setElementProperties(resultsView, ".js-editable-date-added input", { value: dateAdded });
        // action buttons section
        setElementProperties(resultsView, ".js-view-gc", { value: shareId });
        setElementProperties(resultsView, ".js-delete-gc", { value: geneCollectionId, dataset: { isOwner } });
        setElementProperties(resultsView, ".js-download-gc", { dataset: { gcShareId: shareId, gcId: geneCollectionId, gcType: gctype } });
        setElementProperties(resultsView, ".js-share-gc", { value: shareId, dataset: { gcId: geneCollectionId } });
        setElementProperties(resultsView, ".js-edit-gc", { value: geneCollectionId, dataset: { gcId: geneCollectionId } });
        setElementProperties(resultsView, ".js-edit-gc-save", { value: geneCollectionId, dataset: { gcId: geneCollectionId } });
        setElementProperties(resultsView, ".js-edit-gc-cancel", { value: geneCollectionId, dataset: { gcId: geneCollectionId } });
        // gene collection type section
        setElementProperties(resultsView, ".js-display-gctype span:last-of-type", { textContent: gctype });
        setElementProperties(resultsView, ".js-editable-gctype input", { value: gctype });
        // long description section
        setElementProperties(resultsView, ".js-display-ldesc", { id: `result_gc_id_${geneCollectionId}_display_ldesc_container` });
        setElementProperties(resultsView, ".js-editable-ldesc textarea", { id: `result_gc_id_${geneCollectionId}_editable_ldesc`, dataset: { originalVal: longDesc }, value: longDesc });
        // preview genes stuff
        setElementProperties(resultsView, ".js-preview-genes-button-container", { id: `${geneCollectionId}_preview_genes_container` });
        setElementProperties(resultsView, ".js-preview-genes-container", { id: `${geneCollectionId}_gene_container` });

        // Append the cloned template to the results container
        resultsContainer.appendChild(resultsView);

        addVisibilityInfoToGeneCollection(geneCollectionId, isPublic);

        addPreviewGenesToGeneCollection(geneCollectionId, gctype, shareId, geneCount);

        // Add weighted or unweighted gene collection to DOM
        const geneCollectionIdContainer = document.getElementById(`${geneCollectionId}_gene_container`);
        const geneInfoContainer = document.createElement("div");
        geneInfoContainer.id = `gc_${geneCollectionId}_gene_info`;
        geneCollectionIdContainer.appendChild(geneInfoContainer);

        // Add ldesc if it exists
        const ldescContainer = document.getElementById(`result_gc_id_${geneCollectionId}_display_ldesc_container`);
        const ldescElt = document.createElement("p");
        ldescElt.id = `result_gc_id_${geneCollectionId}_display_ldesc`;
        ldescElt.textContent = longDesc || "No description entered";
        ldescContainer.appendChild(ldescElt);

    }

    // Hide some buttons if user is not owner
    updateGeneCollectionListButtons();

    for (const tooltipElt of document.getElementsByClassName("tooltip")) {
        // Remove any existing tooltips
        tooltipElt.remove();
    }

    // Create tooltips for all elements with the data-tooltip-content attribute
    // Only creating one set so that they can be reused
    const actionGroupElt = document.querySelector(".js-action-links");
    const tooltips = []
    for (const classElt of actionGroupElt.querySelectorAll("[data-tooltip-content]")) {
        tooltips.push(createActionTooltips(classElt))
    }

    // Then apply each tooltip to the appropriate element for all elements with the data-tooltip-content attribute

    for (const actionElt of document.querySelectorAll(".js-action-links")) {
        const loopTooltips = [...tooltips];
        for (const classElt of actionElt.querySelectorAll("[data-tooltip-content]")) {
            applyTooltip(classElt, loopTooltips.shift());
        }
    }


    // Initiialize delete gene collection popover for each delete button
    createDeleteConfirmationPopover();

    // All event listeners for all gene collection elements
    addGeneCollectionEventListeners();
}

/**
 * Resets the add form by removing the "is-hidden" class from the save button,
 * adding it to the saving button, and resetting all form fields to their default values.
 *
 * Also removes the "has-background-primary" class from the unweighted and weighted headers,
 * sets their color to black, and enables the paste and upload buttons.
 *
 * Finally, adds the "is-hidden" class to the form and pasted genes container, and sets isAddFormOpen to false.
 */
const resetAddForm = () => {
    document.getElementById("btn_new_collection_save").classList.remove("is-loading");

    document.getElementById("new_collection_label").value = "";
    document.getElementById("new_collection_ldesc").value = "";
    document.getElementById("new_collection_pasted_genes").value = "";
    document.getElementById("new_collection_file").value = "";
    document.getElementById("new_collection_file_name").value = "";

    document.getElementById("new_collection_visibility").checked = false;
    document.getElementById("new_collection_visibility").closest(".field").querySelector("label").textContent = "Private";

    for (const classElt of document.getElementsByClassName("js-new-collection-header")) {
        classElt.classList.remove('has-background-primary', 'has-text-white');
        classElt.classList.add("has-text-dark");
    }

    for (const classElt of document.getElementsByClassName("js-upload-gc-btn")) {
        classElt.removeAttribute('disabled');
    }

    // Currently this is not implemented yet
    document.getElementById("btn_gc_upload_labeled_list").setAttribute('disabled', 'disabled');

    document.getElementById("new_collection_upload_c").classList.remove("is-hidden");

    document.getElementById("new_collection_form_c").classList.add("is-hidden");
    document.getElementById("new_collection_pasted_genes_c").classList.add("is-hidden");
    document.getElementById("new_collection_file_name").classList.add("is-hidden");

    for (const classElt of document.querySelectorAll("#new_collection_form_c .js-validation-help")) {
        console.log(classElt);
        classElt.remove();
    }

    // Remove all "is-danger" classes from form fields
    for (const classElt of document.querySelectorAll("#new_collection_form_c .is-danger")) {
        console.log(classElt);
        classElt.classList.remove("is-danger");
    }

    isAddFormOpen = false;
}

/**
 * Sets the dataset of an element based on the provided dataset object.
 * @param {HTMLElement} parentNode - The parent node containing the element.
 * @param {string} selector - The CSS selector to select the element.
 * @param {Object} dataset - The dataset object containing key-value pairs.
 */
const setElementDataset = (parentNode, selector, dataset) => {
    const element = parentNode.querySelector(selector);
    Object.keys(dataset).forEach((key) => {
        element.dataset[key] = dataset[key];
    });
}

/**
 * Sets the properties of an element selected by a given selector within a parent node.
 * @param {HTMLElement} parentNode - The parent node containing the element.
 * @param {string} selector - The CSS selector used to select the element.
 * @param {Object} properties - An object containing the properties to be set on the element.
 */
const setElementProperties = (parentNode, selector, properties) => {
    const element = parentNode.querySelector(selector);
    Object.keys(properties).forEach((property) => {
        if (property === "dataset") {
            setElementDataset(parentNode, selector, properties[property]);
            return;
        }
        element[property] = properties[property];
    });
}

/**
 * Sets up the gene list toggle functionality.
 *
 * @param {string} className - The class name of the elements that trigger the toggle.
 * @param {string} ajaxUrl - The URL for the AJAX request.
 * @param {Function} handleData - The callback function to handle the retrieved data.
 */
const setupGeneListToggle = (className, ajaxUrl, handleData) => {
    // Show genes when button is clicked & adjust button styling
    for (const classElt of document.getElementsByClassName(className)) {
        classElt.addEventListener("click", async (e) => {
            const button = e.currentTarget;
            const gcId = e.currentTarget.dataset.gcId;
            const shareId = e.currentTarget.dataset.gcShareId;
            const geneCollectionIdContainer = document.getElementById(`${gcId}_gene_container`);
            const previewText = document.getElementById(`btn_gc_${gcId}_text`)

            // Toggle gene container visibility
            if (geneCollectionIdContainer.classList.contains("is-hidden")) {
                button.classList.remove("is-outlined");

                // If the preview table already exists, just show it
                if (document.querySelector(`#gc_${gcId}_gene_info > .js-info-container`)) {
                    geneCollectionIdContainer.classList.remove("is-hidden");    // TODO: add animate CSS with fade in

                    button.querySelector("i").classList.add("mdi-eye-off");
                    button.querySelector("i").classList.remove("mdi-format-list-bulleted");
                    previewText.textContent = "Hide";
                    return;
                }

                button.classList.add("is-loading");

                // Create the preview table of the first five genes
                try {
                    const {data} = await axios.post(ajaxUrl, convertToFormData({
                        'share_id': shareId
                    }));

                    if (data.success < 1) {
                        throw Error(data.message);
                    }

                    const infoContainer = document.createElement("div");
                    infoContainer.classList.add("js-info-container");
                    document.getElementById(`gc_${gcId}_gene_info`).appendChild(infoContainer);

                    handleData(infoContainer, data);

                    geneCollectionIdContainer.classList.remove("is-hidden");    // TODO: add animate CSS with fade in

                    button.querySelector("i").classList.add("mdi-eye-off");
                    button.querySelector("i").classList.remove("mdi-format-list-bulleted");
                    previewText.textContent = "Hide";

                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to load gene collection preview")
                } finally {
                    button.classList.remove("is-loading");
                }

                return;
            }
            geneCollectionIdContainer.classList.add("is-hidden");
            button.classList.add("is-outlined");
            button.blur();
            button.querySelector("i").classList.remove("mdi-eye-off");
            button.querySelector("i").classList.add("mdi-format-list-bulleted");
            previewText.textContent = previewText.dataset.offState;

        });
    }
}

const setupPagination = (pagination) => {

        // Update result count and label
        document.getElementById("result_count").textContent = pagination.total_results;
        document.getElementById("result_label").textContent = pagination.total_results == 1 ? " result" : " results";
        document.getElementById("gc_count_label_c").classList.remove("is-hidden");

        const firstResult = pagination.total_results > 0 ? (pagination.current_page - 1) * resultsPerPage + 1 : 0;
        const lastResult = Math.min(pagination.current_page * resultsPerPage, pagination.total_results);
        document.getElementById("result_range").textContent = `${firstResult} - ${lastResult}`;

        // Update pagination buttons
        for (const classElt of document.getElementsByClassName("pagination")) {
            classElt.replaceChildren();
            classElt.classList.remove("is-hidden");

            const paginationList = document.createElement("ul");
            paginationList.className = "pagination-list";
            classElt.appendChild(paginationList);


            // Add previous button
            if (pagination.current_page > 1) {
                paginationList.appendChild(createPaginationButton(null, 'left', () => {
                    submitSearch(pagination.current_page - 1);
                }));
            }

            // Add page buttons but only show 3 at a time (1 ... 3 4 5 ... 10)
            const startPage = Math.max(1, pagination.current_page - 1);
            const endPage = Math.min(pagination.total_pages, pagination.current_page + 1);

            if (startPage > 1) {
                paginationList.appendChild(createPaginationButton(1, null, () => {
                    submitSearch(1);
                }));
            }

            if (startPage > 2) {
                paginationList.appendChild(createPaginationEllipsis());
            }

            for (let i = startPage; i <= endPage; i++) {
                const li = paginationList.appendChild(createPaginationButton(i, null, () => {
                    submitSearch(i);
                }));
                if (i == pagination.current_page) {
                    li.firstChild.classList.add("is-current");
                    li.firstChild.classList.remove("is-outlined");
                }
            }

            if (endPage < pagination.total_pages - 1) {
                paginationList.appendChild(createPaginationEllipsis());
            }

            if (endPage < pagination.total_pages) {
                paginationList.appendChild(createPaginationButton(pagination.total_pages, null, () => {
                    submitSearch(pagination.total_pages);
                }));
            }

            // Add next button
            if (pagination.current_page < pagination.total_pages) {
                paginationList.appendChild(createPaginationButton(null, 'right', () => {
                    submitSearch(pagination.current_page + 1);
                }));
            }
        }

}


/**
 * Displays an action note for a gene collection.
 *
 * @param {string} gcId - The ID of the gene collection.
 * @param {string} shareUrl - The URL to be copied to the clipboard.
 * @returns {void}
 */
const showGcActionNote = (gcId, shareUrl) => {
    // sanitize shareUrl
    shareUrl = shareUrl.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

    const noteSelector = `#result_gc_id_${gcId} span.js-gc-action-note`;
    const noteSelecterElt = document.querySelector(noteSelector);

    noteSelecterElt.innerHTML = copyToClipboard(shareUrl) ? "URL copied to clipboard" : `Failed to copy to clipboard. URL: ${shareUrl}`;
    noteSelecterElt.classList.remove("is-hidden");

    setTimeout(() => {
        noteSelecterElt.classList.add("is-hidden"); // TODO: add animate CSS with fade out
        noteSelecterElt.innerHTML = "";
    }, 5000);
}

/**
 * Submits a search for gene collections based on the user's search terms and filter options.
 * @function
 * @returns {void}
 */
const submitSearch = async (page) => {

    // hide pagination (but keep results were they are)
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.add("is-invisible");
    }

    // Clear any existing results
    const resultsContainer = document.getElementById("results_container");
    for (const elt of resultsContainer.querySelectorAll(":not(#results_view)")) {
        elt.remove()
    }

    const searchTerms = document.getElementById("search_terms").value;

    // If this is the first time searching with terms, set the sort by to relevance
    if (searchTerms && firstSearch) {
        document.getElementById("sort_by").value = 'relevance';
        firstSearch = false;
    }

    const searchCriteria = {
        'session_id': CURRENT_USER.session_id,
        'search_terms': searchTerms,
        'sort_by': document.getElementById("sort_by").value
    };

    // collect the filter options the user defined
    searchCriteria.organism_ids = buildFilterString('controls_organism');
    searchCriteria.date_added = buildFilterString('controls_date_added');
    searchCriteria.ownership = buildFilterString('controls_ownership');
    searchCriteria.limit = resultsPerPage;
    searchCriteria.page = page || 1;

    try {
        const {data} = await axios.post('./cgi/search_gene_carts.cgi', convertToFormData(searchCriteria));
        processSearchResults(data);

        setupPagination(data.pagination);
    } catch (error) {
        logErrorInConsole(error);
        createToast("Failed to search gene collections");
    } finally {
        // show pagination
        for (const classElt of document.getElementsByClassName("pagination")) {
            classElt.classList.remove("is-invisible");
        }
    }
}

/**
 * Toggles the editable mode of elements with the given selector base.
 * @param {boolean} hideEditable - Whether to hide the editable elements or not.
 * @param {string} [selectorBase=""] - The base selector to use for finding the editable and non-editable elements.
 */
const toggleEditableMode = (hideEditable, selectorBase="") => {
    const editableElements = document.querySelectorAll(`${selectorBase} .js-editable-version`);
    const nonEditableElements = document.querySelectorAll(`${selectorBase} .js-display-version`);

    editableElements.forEach(el => el.classList.toggle("is-hidden", hideEditable));
    nonEditableElements.forEach(el => el.classList.toggle("is-hidden", !hideEditable));
}

/**
 * Updates the visibility of edit and delete buttons for each gene collection in the result list based on the current user's ID.
 * @function
 * @returns {void}
 */
const updateGeneCollectionListButtons = () => {
    const gcListElements = document.getElementsByClassName("js-gc-list-element");

    for (const classElt of gcListElements) {

        // The ability to edit and delete and dataset are currently paired
        const deleteButton = classElt.querySelector("button.js-delete-gc");
        const editButton = classElt.querySelector("button.js-edit-gc");

        if (deleteButton.dataset.isOwner === "false") {
            deleteButton.parentElement.remove();    // remove .control element to prevent heavy line where button was
            editButton.parentElement.remove()

            const geneCollectionId = classElt.dataset.gcId;
            const selectorBase = `#result_gc_id_${geneCollectionId}`;

            // Delete all editable elements to prevent editing in the DOM
            for (const editableElt of classElt.querySelectorAll(`${selectorBase} .js-editable-version`)) {
                editableElt.remove();
            }
        }
    };
}

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// User settings has no "active" state for the sidebar
	document.querySelector("#header_bar .navbar-item").textContent = "Gene Collection Manager";
	for (const elt of document.querySelectorAll("#primary_nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}

	document.querySelector("a[tool='manage_genes'").classList.add("is-active");

    const sessionId = CURRENT_USER.session_id;

	if (! sessionId ) {
        document.getElementById("not_logged_in_msg").classList.remove("is-hidden");
        document.getElementById("create_new_gene_collection").classList.add("is-hidden");
        document.getElementById("create_new_gene_collection").setAttribute("disabled", "disabled");
        // only show public gene collections option
        for (const elt of document.querySelectorAll("#controls_ownership li:not([data-dbval='public'])")) {
            elt.remove();
        }
        document.querySelector("#controls_ownership li[data-dbval='public']").classList.add("js-selected");
    }

    await loadOrganismList();

    submitSearch();

    // Settings for selected facets
    for (const elt of document.querySelectorAll("ul.js-expandable-target li")) {
        elt.addEventListener("click", (e) => {
            if (e.target.classList.contains("js-all-selector")) {
                // if the one clicked is the all_selector then highlight it and unclick the rest
                for (const elt of e.target.parentElement.children) {
                    elt.classList.remove("js-selected");
                }

                e.target.classList.add("js-selected");

            } else if (e.target.classList.contains("js-selected")) {
                // If turning off, make sure at least one other option is selected, else set "all" option
                e.target.classList.remove("js-selected");

                if (!e.target.parentElement.querySelectorAll("li.js-selected")) {
                    e.target.parentElement.querySelector("li.js-all-selector").classList.add("js-selected");
                }
            } else {
                // If turning on, make sure all_selector is off
                if (e.target.parentElement.querySelector("li.js-all-selector")) {
                    // In case not logged in and "All" is not an option
                    e.target.parentElement.querySelector("li.js-all-selector").classList.remove("js-selected");
                }

                // If this selection group has the 'only_one' option deselect the rest
                if (e.target.parentElement.classList.contains("js-choose-only-one")) {
                    for (const elt of e.target.parentElement.children) {
                        elt.classList.remove("js-selected");
                    }
                }

                e.target.classList.add("js-selected");
            }

            submitSearch();
        });
    }

};

// validate that #new_collection_label input has a value
document.getElementById("new_collection_label").addEventListener("blur", (e) => {
    e.target.classList.remove("is-danger");
    // Remove small helper text under input
    const helperText = e.target.parentElement.querySelector("p.help");
    if (helperText) {
        helperText.remove();
    }

    if (e.target.value) {
        return;
    }
    e.target.classList.add("is-danger");
    // Add small helper text under input
    const newHelperElt = document.createElement("p");
    newHelperElt.classList.add("help", "has-text-danger-dark", "js-validation-help");
    newHelperElt.textContent = "Please enter a value";
    e.target.parentElement.appendChild(newHelperElt);
});

document.getElementById("search_clear").addEventListener("click", () => {
    document.getElementById("search_terms").value = "";
    submitSearch();
});

// Search for gene collections using the supplied search terms
const searchTermsElt = document.getElementById("search_terms");
searchTermsElt.addEventListener("keyup", (event) => {
    const searchTerms = searchTermsElt.value;
    const searchClearElt = document.getElementById("search_clear");
    searchClearElt.classList.add("is-hidden");
    if (searchTerms) {
        searchClearElt.classList.remove("is-hidden");
    }
    if (event.key === "Enter") {
        submitSearch();
    }
});

// Changing sort by criteria should update the search results
document.getElementById("sort_by").addEventListener("change", () => {
    submitSearch();
});

const btnCreateCartToggle = document.getElementById("create_new_gene_collection");
btnCreateCartToggle.addEventListener("click", () => {
    const createCollectionContainer = document.getElementById("create_collection_container");
    const gcViewport = document.getElementById("results_container");
    const viewControls = document.getElementById("view_controls");
    const newCartVisibility = document.getElementById("new_collection_visibility");
    const belowPagination = document.getElementById("below_pagination_container");

    if (createCollectionContainer.classList.contains("is-hidden")) {
        createCollectionContainer.classList.remove("is-hidden"); // TODO: Add animation with fade in/out
        gcViewport.classList.add("is-hidden");
        viewControls.classList.add("is-hidden");
        belowPagination.classList.add("is-hidden");
        btnCreateCartToggle.innerHTML =`<span class="icon"><i class="mdi mdi-undo-variant"></i></span> <span>Cancel new collection</span>`;
        newCartVisibility.checked = false; // bootstrap toggle off
        newCartVisibility.closest(".field").querySelector("label").textContent = "Private";
        return;
    }
    createCollectionContainer.classList.add("is-hidden");
    gcViewport.classList.remove("is-hidden");   // TODO: Add animation with fade in/out
    viewControls.classList.remove("is-hidden"); // TODO: Add animation with fade in/out
    belowPagination.classList.remove("is-hidden");
    btnCreateCartToggle.innerHTML =`<span class="icon"><i class="mdi mdi-plus"></i></span> <span>Create new gene collection</span>`;
    resetAddForm();
});

// If an upload method is clicked, show the appropriate form and hide the others
document.getElementById("btn_gc_paste_unweighted_list").addEventListener("click", () => {
    document.getElementById("new_collection_unweighted_header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new_collection_unweighted_header").classList.remove("has-text-dark");

    document.getElementById("btn_gc_upload_unweighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_weighted_list").setAttribute('disabled', 'disabled');

    document.getElementById("new_collection_form_c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out
    document.getElementById("new_collection_pasted_genes_c").classList.remove("is-hidden");

    document.getElementById("new_collection_upload_c").classList.add("is-hidden");

    document.getElementById("new_collection_upload_type").value = "pasted_genes";
    isAddFormOpen = true;
});

document.getElementById("btn_gc_upload_unweighted_list").addEventListener("click", () => {
    if (isAddFormOpen) return;
    document.getElementById("new_collection_unweighted_header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new_collection_unweighted_header").classList.remove("has-text-dark");

    document.getElementById("btn_gc_paste_unweighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_weighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_labeled_list").setAttribute('disabled', 'disabled');

    document.getElementById("new_collection_form_c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out

    document.getElementById("new_collection_upload_type").value = "uploaded-unweighted";
    document.getElementById("new_collection_upload_c").classList.remove("is-hidden");
    isAddFormOpen = true;
});

document.getElementById("btn_gc_upload_weighted_list").addEventListener("click", () => {
    if (isAddFormOpen) return;
    document.getElementById("new_collection_weighted_header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new_collection_weighted_header").classList.remove("has-text-dark");

    document.getElementById("btn_gc_paste_unweighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_unweighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_labeled_list").setAttribute('disabled', 'disabled');

    document.getElementById("new_collection_form_c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out

    document.getElementById("new_collection_upload_type").value = "uploaded-weighted";
    document.getElementById("new_collection_upload_c").classList.remove("is-hidden");
    isAddFormOpen = true;
});

document.getElementById("btn_gc_upload_labeled_list").addEventListener("click", () => {
    if (isAddFormOpen) return;
    document.getElementById("new_collection_labeled_header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new_collection_labeled_header").classList.remove("has-text-dark");

    document.getElementById("btn_gc_paste_unweighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_unweighted_list").setAttribute('disabled', 'disabled');
    document.getElementById("btn_gc_upload_weighted_list").setAttribute('disabled', 'disabled');

    document.getElementById("new_collection_form_c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out

    document.getElementById("new_collection_upload_type").value = "uploaded-labeled";
    document.getElementById("new_collection_upload_c").classList.remove("is-hidden");
    isAddFormOpen = true;
});

// If the cancel button is clicked, hide the form and show the upload buttons
document.getElementById("btn_new_collection_cancel").addEventListener("click", () => {
    document.getElementById("create_new_gene_collection").click();
    document.getElementById("new_collection_pasted_genes_c").classList.add("is-hidden");
});


const btnNewCartSave = document.getElementById("btn_new_collection_save");
btnNewCartSave.addEventListener("click", (e) => {
    // disable button and show indicator that it's loading
    btnNewCartSave.classList.add("is-loading");

    let validationFailed = false;

    // check required fields
    const newCartLabel = document.getElementById("new_collection_label");
    if (! newCartLabel.value) {
        newCartLabel.classList.add("is-danger");
        // Add small helper text under input
        const newHelperElt = document.createElement("p");
        newHelperElt.classList.add("help", "has-text-danger-dark", "js-validation-help");
        newHelperElt.textContent = "Please enter a value";
        newCartLabel.parentElement.appendChild(newHelperElt);

        btnNewCartSave.classList.remove("is-loading");
        validationFailed = true;
    }

    const newCartOrganism = document.getElementById("new_collection_organism_id");
    if (! newCartOrganism.value) {
        newCartOrganism.parentElement.classList.add("is-danger");
        // Add small helper text under input
        const newHelperElt = document.createElement("p");
        newHelperElt.classList.add("help", "has-text-danger-dark", "js-validation-help");
        newHelperElt.textContent = "Please select an organism";
        newCartOrganism.parentElement.appendChild(newHelperElt);

        btnNewCartSave.classList.remove("is-loading");
        validationFailed = true;
    }

    // Was file uploaded or genes pasted?
    const uploadType = document.getElementById("new_collection_upload_type").value;
    if (uploadType === "pasted_genes") {
        const newCartPastedGenes = document.getElementById("new_collection_pasted_genes");
        if (! newCartPastedGenes.value) {
            newCartPastedGenes.classList.add("is-danger");
            // Add small helper text under input
            const newHelperElt = document.createElement("p");
            newHelperElt.classList.add("help", "has-text-danger-dark", "js-validation-help");
            newHelperElt.textContent = "Please enter a value";
            newCartPastedGenes.parentElement.appendChild(newHelperElt);
            validationFailed = true;
        }
    } else if (["uploaded-unweighted", "uploaded-weighted", "uploaded-labeled"].includes(uploadType)) {
        const newCartFile = document.getElementById("new_collection_file");
        if (! newCartFile.value) {
            newCartFile.parentElement.parentElement.classList.add("is-danger");
            // Add small helper text under input
            const newHelperElt = document.createElement("p");
            newHelperElt.classList.add("help", "has-text-danger-dark", "js-validation-help", "ml-2");
            newHelperElt.textContent = "Please select a file";
            newCartFile.parentElement.appendChild(newHelperElt);
            validationFailed = true;
        }
    }

    if (validationFailed) {
        return;
    }

    // Remove all "is-danger" classes from form fields
    for (const classElt of document.querySelectorAll("#new_collection_form_c .js-validation-help")) {
        classElt.remove();
    }

    // Remove all "is-danger" classes from form fields
    for (const classElt of document.querySelectorAll("#new_collection_form_c .is-danger")) {
        classElt.classList.remove("is-danger");
    }

    // Passed to CGI script as 1 or 0
    const isPublic = document.getElementById("new_collection_visibility").checked ? 1 : 0;

    const payload = {
        'new_cart_label': newCartLabel.value
        , 'new_cart_organism_id': newCartOrganism.value
        , 'new_cart_ldesc': document.getElementById("new_collection_ldesc").value
        , 'is_public': isPublic
        , 'session_id': CURRENT_USER.session_id
        , 'new_cart_upload_type': document.getElementById("new_collection_upload_type").value
        , "new_cart_pasted_genes": document.getElementById("new_collection_pasted_genes").value
    }

    const gc = new GeneCart();
    gc.addCartToDbFromForm(payload, geneCollectionSaved, geneCollectionFailure);
    btnNewCartSave.classList.remove("is-loading");

});

document.getElementById("btn_list_view_compact").addEventListener("click", () => {
    document.getElementById("btn_list_view_compact").classList.add('is-gear-bg-secondary');
    document.getElementById("btn_list_view_compact").classList.remove('is-dark');
    document.getElementById("btn_list_view_expanded").classList.remove('is-gear-bg-secondary');
    document.getElementById("btn_list_view_expanded").classList.add('is-dark');

    // find all elements with class 'js-expandable-view' and make sure they also have 'expanded-view-hidden'
    for (const elt of document.querySelectorAll(".js-expandable-view")){
        elt.classList.add("is-hidden");
    };

    // Toggle js-expand-box to icon collapse
    for (const elt of document.querySelectorAll(".js-expand-box i")){
        elt.classList.remove("mdi-arrow-collapse");
        elt.classList.add("mdi-arrow-expand");
    }

});

document.getElementById("btn_list_view_expanded").addEventListener("click", () => {
    document.getElementById("btn_list_view_compact").classList.remove('is-gear-bg-secondary');
    document.getElementById("btn_list_view_compact").classList.add('is-dark');
    document.getElementById("btn_list_view_expanded").classList.add('is-gear-bg-secondary');
    document.getElementById("btn_list_view_expanded").classList.remove('is-dark');

    // find all elements with class 'js-expandable-view' and make sure they also have 'expanded-view-hidden'
    for (const elt of document.querySelectorAll(".js-expandable-view")){
        elt.classList.remove("is-hidden");
    };

    // Toggle js-expand-box to icon expand
    for (const elt of document.querySelectorAll(".js-expand-box i")){
        elt.classList.add("mdi-arrow-collapse");
        elt.classList.remove("mdi-arrow-expand");
    }

});


// Collpasible view toggle for left facet menu
for (const elt of document.querySelectorAll(".js-expandable-control")) {
    elt.addEventListener("click", (e) => {
        const exblock = e.currentTarget.parentElement.querySelector(".js-expandable-target");
        const iconElt = e.currentTarget.querySelector(".mdi");

        if (exblock.classList.contains("is-hidden")) {
            iconElt.classList.remove("mdi-plus");
            iconElt.classList.add("mdi-minus");
            exblock.classList.remove("is-hidden");  // TODO: Add animation with fade in/out
            return;
        }
        iconElt.classList.remove("mdi-minus");
        iconElt.classList.add("mdi-plus");
        exblock.classList.add("is-hidden"); // TODO: Add animation with fade in/out
    }
)};


// Adjust the visibility label of the new cart form when the toggle is clicked
document.getElementById("new_collection_visibility").addEventListener("change", (e) => {
    const isPublic = e.currentTarget.checked;
    e.currentTarget.dataset.isPublic = isPublic;
    e.currentTarget.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";
});

// When user uploads file, update the file name in the form
document.getElementById("new_collection_file").addEventListener("change", (e) => {
    const file = e.currentTarget.files[0];
    console.log(file);
    document.getElementById("new_collection_file_name").textContent = file.name;
    console.log(document.getElementById("new_collection_file_name"));
    document.getElementById("new_collection_file_name").classList.remove("is-hidden");
});

// Show modal of upload requirements if clicked
document.getElementById("collection_upload_reqs").addEventListener("click", (e) => {
    e.preventDefault();
    const modalId = e.target.dataset.target;
    const modalElt = document.getElementById(modalId);
    openModal(modalElt);

    const listType = document.getElementById("new_collection_upload_type").value;

    if (listType === "uploaded-unweighted") {
        document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
            document.getElementById("unweighted_example").innerHTML;
    } else if (listType === "uploaded-weighted") {
        document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
            document.getElementById("weighted_example").innerHTML;
    } else if (listType === "uploaded-labeled") {
        document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
            document.getElementById("labeled_example").innerHTML;

    } else {
        throw new Error("Unsupported file-upload list type");
    }

});
