"use strict";

let firstSearch = true;
const animationTime = 200;  // in ms
let isAddFormOpen = false;
const resultsPerPage = 20;

// TODO - Add transformation code for quick gene collection transformations
// TODO - Add "labeled" gene collection type (i.e. marker genes) and integrate into other code snippets

// TODO - Add public/private toggle button code for new genecart
// TODO - Fix "edit genecart" button code
// TODO - Add "new genecart" code back from old version
// TODO - Get filter properties working
// TODO - Add Organism select box

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

    // Show genes when button is clicked & adjust button styling
    for (const classElt of document.getElementsByClassName("js-gc-unweighted-gene-list-toggle")) {
        classElt.addEventListener("click", (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            const geneList = document.getElementById(`${gcId}_gene_list`);

            // see if the gene_list is visible and toggle
            if (geneList.classList.contains("is-hidden")) {
                geneList.classList.remove("is-hidden");
                e.currentTarget.classList.remove("is-outlined");
            } else {
                geneList.classList.add("is-hidden");
                e.currentTarget.classList.add("is-outlined");
            }

        });
    }

    for (const classElt of document.getElementsByClassName("js-gc-weighted-gene-list-toggle")) {
        classElt.addEventListener("click", async (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            const shareId = e.currentTarget.dataset.gcShareId;

            document.getElementById(`btn_gc_${gcId}_preview`).classList.add("is-hidden");
            document.getElementById(`btn_gc_${gcId}_loading`).classList.remove("is-hidden");

            try {
                const {data} = await axios.post('./cgi/get_weighted_gene_cart_preview.cgi', convertToFormData({
                    'share_id': shareId
                }));
                processWeightedGcList(gcId, data['preview_json']);
            } catch (error) {
                logErrorInConsole(error);
                // TODO: Display error notification
            }
        });
    }

    // Hide gene collection view
    for (const classElt of document.getElementsByClassName("js-gc-weighted-gene-list-hider")) {
        classElt.addEventListener("click", (e) => {
            const gcId = e.currentTarget.dataset.gcId;
            document.getElementById(`${gcId}_gene_table`).classList.add("is-hidden");   // TODO: add animate CSS with fade out
            document.getElementById(`btn_gc_${gcId}_hider`).classList.add("is-hidden");
            document.getElementById(`btn_gc_${gcId}_preview`).classList.remove("is-hidden");
        });
    }

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
            } else {
                e.currentTarget.innerHTML = '<i class="mdi mdi-arrow-expand"></i>';
            }


        });
    }

    // Reformats the <ul> containing the gene symbols into a text file with one gene per row
    for (const classElt of document.getElementsByClassName("js-download-gc")) {
        classElt.addEventListener("click", (e) => {
            const gcId = e.currentTarget.dataset.gcId;

            // Take the list elements, separate them, and put them in a text file
            const geneListElt = document.getElementById(`${gcId}_gene_list`);
            const fileContents = Array.from(geneListElt.children).map(li => li.innerText).join("\n");

            const element = document.createElement("a");
            element.setAttribute(
                "href",
                `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(fileContents)}`
            );

            element.setAttribute("download", `gene_cart.${e.currentTarget.dataset.gcShareId}.tsv`);
            element.style.display = "none";
            document.body.appendChild(element);
            element.click();
            document.body.removeChild(element);
        });
    }

    for (const classElt of document.getElementsByClassName("js-share-gc")) {
        classElt.addEventListener("click", (e) => {

            const shareId = e.currentTarget.value;
            const currentUrl = window.location.href;
            const currentPage = currentUrl.lastIndexOf("gene_collection_manager.html");
            const shareUrl = `${currentUrl.substring(0, currentPage)}p?c=${shareId}`;
            const gcId = e.currentTarget.dataset.gcId;

            if (copyToClipboard(shareUrl)) {
                showGcActionNote(gcId, "URL copied to clipboard");
            } else {
                showGcActionNote(gcId, `Failed to copy to clipboard. URL: ${shareUrl}`);
            }

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
            for (const classElt of document.querySelectorAll(`${selectorBase} .js-readonly-version`)) {
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
                    'ldesc': newLdesc
                }));

                // Update the UI for the new values
                document.querySelector(`${selectorBase}_editable_visibility`).dataset.isPublic = newVisibility;
                if (newVisibility) {
                    document.querySelector(`${selectorBase}_display_visibility`).innerHTML = "Public gene collection";
                    document.querySelector(`${selectorBase}_display_visibility`).classList.remove("is-danger");
                    document.querySelector(`${selectorBase}_display_visibility`).classList.add("is-light");
                } else {
                    document.querySelector(`${selectorBase}_display_visibility`).innerHTML = "Private gene collection";
                    document.querySelector(`${selectorBase}_display_visibility`).classList.remove("is-light");
                    document.querySelector(`${selectorBase}_display_visibility`).classList.add("is-danger");
                }

                document.querySelector(`${selectorBase}_editable_title`).dataset.originalVal = newTitle;
                document.querySelector(`${selectorBase}_display_title`).innerHTML = newTitle;

                document.querySelector(`${selectorBase}_editable_ldesc`).dataset.originalVal = newLdesc;
                document.querySelector(`${selectorBase}_display_ldesc`).innerHTML = newLdesc;

                document.querySelector(`${selectorBase}_display_organism`).innerHTML =
                    document.querySelector(`${selectorBase}_editable_organism_id > option[value='${newOrgId}']`);
                document.querySelector(`${selectorBase}_editable_organism_id`).dataset.originalVal = newOrgId;

                // Put interface back to view mode.
                toggleEditableMode(true, selectorBase);

            } catch (error) {
                logErrorInConsole(error);
                // TODO: Display error notification
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

            // set the current value as selected
            editableOrganismIdElt.value = editableOrganismIdElt.dataset.originalVal;

            const editableVisibilityElt = document.querySelector(`${selectorBase}_editable_visibility`);

            const isPublic = parseBool(editableVisibilityElt.dataset.isPublic);

            editableVisibilityElt.checked = isPublic;
            editableVisibilityElt.closest(".field").querySelector("label").innerText = isPublic ? "Public" : "Private";

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
 * Adds a gene collection to the gene collection display container in the DOM.
 * @param {string} geneCollectionId - The ID of the gene collection container.
 * @param {string} gctype - The type of gene collection to add.
 * @param {string[]} genes - An array of genes to add to the gene collection.
 * @returns {void}
 */
const addGeneListToGeneCollection = (geneCollectionId, gctype, genes) => {
    // Add weighted or unweighted gene collection to DOM
    const geneCollectionIdContainer = document.getElementById(`${geneCollectionId}_gene_container`);
    if (gctype == "weighted-list") {
        const weightedGeneListContainer = document.createElement("div");
        weightedGeneListContainer.id = `${geneCollectionId}_gene_table`;
        geneCollectionIdContainer.appendChild(weightedGeneListContainer);
        return;
    }
    const geneListContainer = document.createElement("ul");
    geneListContainer.classList.add("is-hidden");
    geneListContainer.id = `${geneCollectionId}_gene_list`;
    geneCollectionIdContainer.appendChild(geneListContainer);
    // append genes to gene collection
    const geneCollectionIdGeneUlElt = document.getElementById(geneListContainer.id);
    for (const gene of genes) {
        const li = document.createElement("li");
        li.innerText = gene;
        li.classList.add("mt-1");
        geneCollectionIdGeneUlElt.appendChild(li);
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
    if (gctype === "weighted-list") {
        const previewButton = document.createElement('button');
        previewButton.id = `btn_gc_${geneCollectionId}_preview`;
        previewButton.className = 'button is-small is-dark is-outlined js-gc-weighted-gene-list-toggle';
        previewButton.title = 'Show gene collection';
        previewButton.innerHTML = '<i class="mdi mdi-format-list-bulleted"></i> Preview';
        previewButton.dataset.gcId = geneCollectionId;
        previewButton.dataset.gcShareId = shareId;

        const loadingButton = document.createElement('button');
        loadingButton.id = `btn_gc_${geneCollectionId}_loading`;
        loadingButton.className = 'button is-small is-outlined is-hidden is-loading';
        loadingButton.title = 'Loading';
        loadingButton.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Loading';
        loadingButton.dataset.gcId = geneCollectionId;

        const hideButton = document.createElement('button');
        hideButton.id = `btn_gc_${geneCollectionId}_hider`;
        hideButton.className = 'button is-small is-outlined js-gc-weighted-gene-list-hider is-hidden';
        hideButton.title = 'Hide gene collection';
        hideButton.innerHTML = '<i class="mdi mdi-eye-off"></i> Hide';
        hideButton.dataset.gcId = geneCollectionId;

        geneCollectionPreviewGenesContainer.append(previewButton, loadingButton, hideButton);
        return;
    } else if (gctype === "unweighted-list") {
        const geneListButton = document.createElement('button');
        geneListButton.className = 'button is-small is-dark is-outlined js-gc-unweighted-gene-list-toggle';
        geneListButton.title = 'Show gene collection';
        geneListButton.innerHTML = `<i class="mdi mdi-format-list-bulleted"></i> ${geneCount} genes`;
        geneListButton.dataset.gcId = geneCollectionId;
        geneCollectionPreviewGenesContainer.append(geneListButton);

    } else if (gctype === "labeled-list") {
        // Not implemented yet
    }
}

/**
 * Adds public/private visibility information to a gene collection display container in the DOM.
 * @param {string} geneCollectionId - The ID of the gene collection display container.
 * @param {boolean} isPublic - A boolean indicating whether the gene collection is public or private.
 * @returns {void}
 */
const addVisibilityInfoToGeneCollection = (geneCollectionId, isPublic) => {

    // add gene collection public/private info to DOM
    const geneCollectionDisplayContainer = document.getElementById(`${geneCollectionId}_display_container`);
    const geneCollectionDisplaySpan = document.createElement("span");
    geneCollectionDisplaySpan.classList.add("tag");
    geneCollectionDisplaySpan.id = `result_gc_id_${geneCollectionId}_display_visibility`;

    if (isPublic) {
        geneCollectionDisplaySpan.classList.add("is-primary", "is-light");
        geneCollectionDisplaySpan.innerText = "Public gene collection";
    } else {
        geneCollectionDisplaySpan.classList.add("is-danger");
        geneCollectionDisplaySpan.innerText = "Private gene collection";
    }
    geneCollectionDisplayContainer.appendChild(geneCollectionDisplaySpan);

    // Toggle switch (public is checked, private is unchecked)
    const visibilitySwitch = document.getElementById(`result_gc_id_${geneCollectionId}_editable_visibility`);

    visibilitySwitch.addEventListener("change", (e) => {
        const isPublic = e.currentTarget.checked;
        e.currentTarget.dataset.isPublic = isPublic;
        e.currentTarget.closest(".field").querySelector("label").innerText = isPublic ? "Public" : "Private";
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
    const selected = document.querySelectorAll(`#${groupName} ul li.js-selected :not(.js-all-selector)`);
    const dbvals = [];

    for (const li of selected) {
        dbvals.push(li.dataset.dbval);
    }

    return dbvals.join(",");
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

                        // This can affect page counts, so we need to re-run the search
                        submitSearch();

                    } else {
                        throw new Error(data['error']);
                    }
                } catch (error) {
                    logErrorInConsole(error);
                    // TODO: Display error notification
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
 * Creates a tooltip element and appends it to the body.
 * @param {HTMLElement} referenceElement - The reference element to which the tooltip is associated.
 * @returns {HTMLElement} The created tooltip element.
 */
const createActionTooltips = (referenceElement) => {
    // Create tooltip element
    const tooltip = document.createElement('div');
    tooltip.className = 'tooltip';
    tooltip.innerText = referenceElement.dataset.tooltipContent;
    tooltip.classList.add("has-background-dark", "has-text-white", "is-hidden");

    // Append tooltip to body
    document.body.appendChild(tooltip);
    return tooltip;
}

// Callbacks after attempting to save a gene collection
const geneCollectionFailure = (gc, message) => {
    logErrorInConsole(message);
}

const geneCollectionSaved = (gc) => {
    document.getElementById("create_new_gene_collection").click();
    submitSearch();
    resetAddForm();
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
        organismChoices.innerHTML = "";
        for (const organism of data.organisms) {
            const li = document.createElement("li");
            li.dataset.dbval = organism.id;
            li.innerText = organism.label;
            organismChoices.appendChild(li);
        }
        const newCollectionOrganismSelect = document.getElementById("new_collection_organism_id");    // <select> element
        newCollectionOrganismSelect.innerHTML = "";
        for (const organism of data.organisms) {
            const option = document.createElement("option");
            option.value = organism.id;
            option.innerText = organism.label;
            newCollectionOrganismSelect.appendChild(option);
        }
    } catch (error) {
        logErrorInConsole(error);
        // TODO: Display error notification
    }
}

/**
 * Loads all the parts of the page which need initial calls from the server, such as
 * database-driven select boxes.
 * @returns {void}
 */
const loadPreliminaryData = () => {
    loadOrganismList();
    document.getElementById("your_gene_collection_filter").click();
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
    resultsContainer.replaceChildren();

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
        const genes = gc.genes;
        const geneCount = gc.gene_count;
        const userName = gc.user_name;
        const organism = gc.organism;
        const isOwner = gc.is_owner;


        // Build results view and add to DOM
        // TODO: Edit "columns" class to be compatible with mobile
        const resultsViewTmpl = `
            <div class="js-gc-list-element columns" id="result_gc_id_${geneCollectionId}" data-gc-id="${geneCollectionId}">
                <div class="column is-full">
                    <!-- title section -->
                    <div class="columns">
                        <div class="column is-11">
                            <div class="js-readonly-version">
                                <p class="has-text-weight-bold" id="result_gc_id_${geneCollectionId}_display_title">${label}</p>
                            </div>
                            <div id="editable_title_c" class="js-editable-version is-hidden">
                                <div class="field">
                                    <span class="has-text-weight-semibold">Title</span>
                                    <input type="text" class="input js-editable-title" id="result_gc_id_${geneCollectionId}_editable_title" data-original-val="${label}" value="${label}"/>
                                </div>
                            </div>
                        </div>

                        <div class="column is-1">
                            <span class="is-clickable is-pulled-right js-expand-box icon" data-gc-id="${geneCollectionId}"><i class="mdi mdi-arrow-expand"></i></span>
                        </div>
                    </div>
                    <!-- visibility/other metadata section -->
                    <div class="columns is-size-7">
                        <div class="column is-3">
                            <div class="js-readonly-version" id="${geneCollectionId}_display_container"></div>
                            <div class="js-editable-version is-hidden">
                                <div class="field">
                                    <label class="label">Visibility</label>
                                    <!-- toggle switch --->
                                    <input type="checkbox" class="switch is-primary" name="result_gc_id_${geneCollectionId}_editable_visibility" id="result_gc_id_${geneCollectionId}_editable_visibility"
                                        data-is-public="${isPublic}"/>
                                    <label for="result_gc_id_${geneCollectionId}_editable_visibility"></label>
                                </div>
                            </div>
                        </div>

                        <!-- organism section -->
                        <div class="column is-3">
                            <div class="js-readonly-version"><span class="has-text-weight-semibold">Organism</span> <span id="result_gc_id_${geneCollectionId}_display_organism">${organism}</span></div>
                            <div class="js-editable-version is-hidden">
                                <div class="field">
                                    <label class="label" for="result_gc_id_${geneCollectionId}_editable_organism_id">Organism</label>
                                    <div class="control">
                                        <div class="select">
                                            <select id="result_gc_id_${geneCollectionId}_editable_organism_id" data-original-val="${organismId}"></select>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <div class="column is-3">
                            <div class="js-readonly-version"><span class="has-text-weight-semibold">Owner</span> ${userName}</div>
                            <div class="js-editable-version is-hidden">
                                <div class="field">
                                    <label class="label">Owner</label>
                                    <div class="control">
                                        <input class="input is-static" value="${userName}" readonly>
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div class="column is-3">
                            <div class="js-readonly-version"><span class="has-text-weight-semibold">Added</span> ${dateAdded}</div>
                            <div class="js-editable-version is-hidden">
                                <div class="field">
                                    <label class="label">Added</label>
                                    <div class="control">
                                        <input class="input is-static" value="${dateAdded}" readonly>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                    <!-- preview + action buttons section -->
                    <div class="columns is-size-7">
                        <div class="column is-3">
                            <div id="${geneCollectionId}_preview_genes_container"></div>
                        </div>
                        <div class="column is-9">
                            <div class="columns">
                                <div class="column is-two-thirds">
                                    <div class="field has-addons js-action-links" role="group">
                                        <p class="control">
                                            <button class="button is-small is-outlined is-dark js-view-gc" value="${shareId}" data-tooltip-content="View on front page">
                                                <span class="icon is-small"><i class="mdi mdi-eye"></i></span>
                                            </button>
                                        </p>
                                        <p class="control">
                                            <button class="button is-small is-outlined is-danger js-delete-gc" value="${geneCollectionId}" data-is-owner="${isOwner}" data-tooltip-content="Delete collection">
                                                <span class="icon is-small"><i class="mdi mdi-delete"></i></span>
                                            </button>
                                        </p>
                                        <p class="control">
                                            <button class="button is-small is-outlined is-dark js-download-gc" data-gc-share-id="${shareId}" data-gc-id="${geneCollectionId}" data-tooltip-content="Download collection as tsv">
                                                <span class="icon is-small"><i class="mdi mdi-download"></i></span>
                                            </button>
                                        </p>
                                        <p class="control">
                                            <button class="button is-small is-outlined is-dark js-share-gc" value="${shareId}" data-gc-id="${geneCollectionId}" data-tooltip-content="Get shareable link">
                                                <span class="icon is-small"><i class="mdi mdi-share-variant"></i></span>
                                            </button>
                                        </p>
                                        <p class="control">
                                            <button class="button is-small is-outlined is-dark js-edit-gc js-readonly-version" value="${geneCollectionId}" data-gc-id="${geneCollectionId}" data-tooltip-content="Edit collection metadata">
                                                <span class="icon is-small"><i class="mdi mdi-pencil"></i></span>
                                            </button>
                                        </p>
                                    </div>
                                    <div class="field has-addons" role="group">
                                        <p class="control">
                                            <button class="button is-primary js-edit-gc-save js-editable-version is-hidden" value="${geneCollectionId}" data-gc-id="${geneCollectionId}">
                                                <span class="icon"><i class="mdi mdi-content-save-edit-outline"></i></span><span>Save edits</span>
                                            </button>
                                        </p>
                                        <p class="control">
                                            <button class="button is-outlined is-primary js-edit-gc-cancel js-editable-version is-hidden" value="${geneCollectionId}" data-gc-id="${geneCollectionId}">
                                                <span class="icon"><i class="mdi mdi-undo-variant"></i></span><span>Cancel edits</span>
                                            </button>
                                        </p>
                                    </div>
                                </div>
                                <div class="column is-one-third">
                                    <div class="js-readonly-version"><span class="has-text-weight-semibold">Type</span> ${gctype}</div>
                                    <div class="js-editable-version is-hidden">
                                        <div class="field">
                                            <label class="label">Type</label>
                                            <div class="control">
                                                <input class="input is-static" value="${gctype}" readonly>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div> <!-- end .columns -->
                            <!-- Info on actions (e.g. "URL copied to clipboard") -->
                            <span class="js-gc-action-note is-hidden"></span>
                        </div>

                    </div> <!-- end .columns -->

                    <!-- Gene preview section -->
                    <div class="pl-4 is-hidden" id="${geneCollectionId}_gene_container"></div>
                    <hr class="js-expandable-view is-hidden" />
                    <!-- Long description section -->
                    <div class="columns js-expandable-view is-hidden">
                        <div class="column is-3"></div>
                        <div class="column is-9 js-readonly-version" id="${geneCollectionId}_ldesc_container">
                            <p class="has-text-weight-semibold">Long description</p>
                        </div>
                        <div class="column is-9 js-editable-version is-hidden">
                            <div class="field">
                                <label class="label">Long description</label>
                                <div class="control">
                                    <textarea class="textarea" id="result_gc_id_${geneCollectionId}_editable_ldesc" rows="4" data-original-val="${longDesc}">${longDesc}</textarea>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div> <!-- end .js-gc-list-element -->
            <hr class="gc-list-element-divider" />
        `

        // Append to results_container
        resultsContainer.insertAdjacentHTML("beforeend", resultsViewTmpl);

        addVisibilityInfoToGeneCollection(geneCollectionId, isPublic);

        addPreviewGenesToGeneCollection(geneCollectionId, gctype, shareId, geneCount);

        addGeneListToGeneCollection(geneCollectionId, gctype, genes);

        // Add ldesc if it exists
        const ldescContainer = document.getElementById(`${geneCollectionId}_ldesc_container`);
        const ldescElt = document.createElement("p");
        ldescElt.id = `result_gc_id_${geneCollectionId}_display_ldesc`;
        ldescElt.innerText = longDesc || "No description entered";
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
 * Toggles the editable mode of elements with the given selector base.
 * @param {boolean} hideEditable - Whether to hide the editable elements or not.
 * @param {string} [selectorBase=""] - The base selector to use for finding the editable and non-editable elements.
 */
const toggleEditableMode = (hideEditable, selectorBase="") => {
    const editableElements = document.querySelectorAll(`${selectorBase} .js-editable-version`);
    const nonEditableElements = document.querySelectorAll(`${selectorBase} .js-readonly-version`);

    editableElements.forEach(el => el.classList.toggle("is-hidden", hideEditable));
    nonEditableElements.forEach(el => el.classList.toggle("is-hidden", !hideEditable));
}

/**
 * Processes the weighted GC list.
 * @param {string} gcId - The ID of the GC.
 * @param {Object} jdata - The preview JSON weighted data to show.
 * @returns {void}
 */
const processWeightedGcList = (gcId, jdata) => {
    document.getElementById(`btn_gc_${gcId}_loading`).classList.add("is-hidden");
    document.getElementById(`btn_gc_${gcId}_hider`).classList.remove("is-hidden");

    // This creates a table with classes dataframe and weighted-list
    document.getElementById(`${gcId}_gene_table`).innerHTML = jdata;
    document.getElementById(`${gcId}_gene_table`).classList.add("is-hidden");   // TODO: add animate CSS with fade in
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
    document.getElementById("new_collection_visibility").closest(".field").querySelector("label").innerText = "Private";

    for (const classElt of document.getElementsByClassName("js-new-collection-header")) {
        classElt.classList.remove('has-background-primary', 'has-text-white');
        classElt.classList.add("has-text-dark");
    }

    for (const classElt of document.getElementsByClassName("js-upload-gc-btn")) {
        classElt.removeAttribute('disabled');
    }

    document.getElementById("new_collection_form_c").classList.add("is-hidden");
    document.getElementById("new_collection_pasted_genes_c").classList.add("is-hidden");
    isAddFormOpen = false;
}

/**
 * Updates the action note for a given gene collection.
 * @param {string} gcId - The ID of the gene collection.
 * @param {string} msg - The message to display in the action note.
 */
const showGcActionNote = (gcId, msg) => {
    const noteSelector = `#result_gc_id_${gcId} span.js-gc-action-note`;
    const noteSelecterElt = document.querySelector(noteSelector);
    noteSelecterElt.innerHTML = msg;
    noteSelecterElt.classList.remove("is-hidden");

    setTimeout(() => {
        noteSelecterElt.classList.add("is-hidden"); // TODO: add animate CSS with fade out
        noteSelecterElt.innerHTML = "";
    }, 5000);
}

const setupPagination = (pagination) => {

        // Update result count and label
        document.getElementById("result_count").textContent = pagination.total_results;
        document.getElementById("result_label").textContent = pagination.total_results == 1 ? " result" : " results";
        document.getElementById("gc_count_label_c").classList.remove("is-hidden");

        const firstResult = (pagination.current_page - 1) * resultsPerPage + 1;
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
 * Submits a search for gene collections based on the user's search terms and filter options.
 * @function
 * @returns {void}
 */
const submitSearch = async (page) => {
    document.getElementById("results_container").replaceChildren();
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
        // TODO: Display error notification
    }
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
        return;
    }

    // Initialize tooltips and popovers
    //createPopover(document.getElementById("collection_upload_reqs"));

    loadPreliminaryData();
};

// validate that #new_collection_label input has a value
document.getElementById("new_collection_label").addEventListener("blur", (e) => {
    e.target.classList.remove("is-danger-dark");
    // Remove small helper text under input
    const helperText = e.target.parentElement.querySelector("p.help");
    if (helperText) {
        helperText.remove();
    }

    if (e.target.value) {
        return;
    }
    e.target.classList.add("is-danger-dark");
    // Add small helper text under input
    const newHelperElt = document.createElement("p");
    newHelperElt.classList.add("help", "is-danger-dark");
    newHelperElt.innerText = "Please enter a value";
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
        btnCreateCartToggle.textContent = "Cancel cart creation";
        newCartVisibility.checked = false; // bootstrap toggle off
        newCartVisibility.closest(".field").querySelector("label").innerText = "Private";
        return;
    }
    createCollectionContainer.classList.add("is-hidden");
    gcViewport.classList.remove("is-hidden");   // TODO: Add animation with fade in/out
    viewControls.classList.remove("is-hidden"); // TODO: Add animation with fade in/out
    belowPagination.classList.remove("is-hidden");
    btnCreateCartToggle.textContent = "Create new cart";
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

    document.getElementById("new_collection_upload_type").value = "pasted_genes";
    document.getElementById("file_upload_c").classList.add("is-hidden");
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
    document.getElementById("new_collection_pasted_genes_c").classList.remove("is-hidden");

    document.getElementById("new_collection_upload_type").value = "uploaded-unweighted";
    document.getElementById("file_upload_c").classList.remove("is-hidden");
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
    document.getElementById("new_collection_pasted_genes_c").classList.remove("is-hidden");

    document.getElementById("new_collection_upload_type").value = "uploaded-weighted";
    document.getElementById("file_upload_c").classList.remove("is-hidden");
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
    document.getElementById("new_collection_pasted_genes_c").classList.remove("is-hidden");

    document.getElementById("new_collection_upload_type").value = "uploaded-weighted";
    document.getElementById("file_upload_c").classList.remove("is-hidden");
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

    // check required fields
    const newCartLabel = document.getElementById("new_collection_label");
    if (! newCartLabel.value) {
        newCartLabel.classList.add("is-danger-dark");
        // Add small helper text under input
        const newHelperElt = document.createElement("p");
        newHelperElt.classList.add("help", "is-danger-dark");
        newHelperElt.innerText = "Please enter a value";
        newCartLabel.parentElement.appendChild(newHelperElt);

        btnNewCartSave.classList.remove("is-loading");
        return;
    }

    newCartLabel.classList.remove("is-danger-dark");
    // Remove small helper text under input
    const helperText = newCartLabel.parentElement.querySelector("p.help");
    if (helperText) {
        helperText.remove();
    }

    // Passed to CGI script as 1 or 0
    const isPublic = document.getElementById("new_collection_visibility").checked ? 1 : 0;

    const formData = new FormData($(this)[0]);
    formData.append('is_public', isPublic);
    formData.append('session_id', CURRENT_USER.session_id);

    const gc = new GeneCart();
    gc.addCartToDbFromForm(formData, geneCollectionSaved, geneCollectionFailure);
    btnNewCartSave.classList.remove("is-loading");

});

document.getElementById("btn_list_view_compact").addEventListener("click", () => {
    document.getElementById("btn_arrangement_view").classList.remove('active');
    document.getElementById("btn_list_view_compact").classList.add('active');
    document.getElementById("btn_list_view_expanded").classList.remove('active');

    document.getElementById("gc_list_c").classList.remove("is-hidden");

    // find all elements with class 'js-expandable-view' and make sure they also have 'expanded-view-hidden'
    for (const elt of document.querySelectorAll(".js-expandable-view")){
        elt.classList.add("expanded-view-hidden");
    };
});

document.getElementById("btn_list_view_expanded").addEventListener("click", () => {
    document.getElementById("btn_arrangement_view").classList.remove('active');
    document.getElementById("btn_list_view_compact").classList.remove('active');
    document.getElementById("btn_list_view_expanded").classList.add('active');

    document.getElementById("gc_list_c").classList.remove("is-hidden");

    // find all elements with class 'js-expandable-view' and make sure they also have 'expanded-view-hidden'
    for (const elt of document.querySelectorAll(".js-expandable-view")){
        elt.classList.remove("expanded-view-hidden");
    };
});


// Generic function to handle all collapsable menus
// h.expandable_control is clicked and looks for plus/minus icons as siblings
// and an .expandable_target as a direct child

for (const elt of document.querySelectorAll("h4.expandable_control")) {
    elt.addEventListener("click", (e) => {
        const exblock = e.target.nextElementSibling;
        if (exblock.classList.contains("is-hidden")) {
            e.target.querySelector(".mdi-plus").classList.add("is-hidden");
            e.target.querySelector(".mdi-minus").classList.remove("is-hidden");
            exblock.classList.remove("is-hidden");  // TODO: Add animation with fade in/out

            if (e.target.classList.contains("profile_control")) {
                for (const elt of document.getElementsByClassName(".profile_control")) {
                    elt.classList.remove("is-hidden");
                    document.getElementById("btn_arrangement_view").classList.remove("is-hidden");
                }
            }
            return;
        }
        e.target.querySelector(".mdi-plus").classList.remove("is-hidden");
        e.target.querySelector(".mdi-minus").classList.add("is-hidden");
        exblock.classList.add("is-hidden"); // TODO: Add animation with fade in/out

        if (e.target.classList.contains("profile_control")) {
            for (const elt of document.getElementsByClassName(".profile_control")) {
                elt.classList.add("is-hidden");
                document.getElementById("btn_arrangement_view").classList.add("is-hidden");
            }
        }
    }
)};

// Generic function to handle the facet selector choices
//  For any ul.controls_filter_options the list elements can have a class="selected"
//  The groups of <li> also have one/top li with class="all_selector" which
//  toggles the rest of them off since no filter is applied.

for (const elt of document.querySelectorAll("ul.controls_filter_options li")) {
    elt.addEventListener("click", (e) => {
        // if the one clicked is the all_selector then highlight it and unclick the rest
        if (e.target.classList.contains("all_selector")) {
            if (! e.target.classList.contains("selected")) {
                e.target.classList.add("selected");
            }

            for (const elt of e.target.parentElement.children) {
                elt.classList.remove("selected");
            }
        } else if (e.target.classList.contains("selected")) {
            // If turning off, make sure at least one other option is selected, else set
            //  set all_selector on
            e.target.classList.remove("selected");

            if (e.target.parentElement.querySelector("li.selected") == null) {
                e.target.parentElement.querySelector("li.all_selector").classList.add("selected");
            }
        } else {
            // If turning on, make sure all_selector is off
            e.target.parentElement.querySelector("li.all_selector").classList.remove("selected");

            // If this selection group has the 'only_one' option deselect the rest
            if (e.target.parentElement.classList.contains("only_one")) {
                for (const elt of e.target.parentElement.children) {
                    elt.classList.remove("selected");
                }
            }

            e.target.classList.add("selected");
        }

        submitSearch();
    });
}

// Adjust the visibility label of the new cart form when the toggle is clicked
document.getElementById("new_collection_visibility").addEventListener("change", (e) => {
    const isPublic = e.currentTarget.checked;
    e.currentTarget.dataset.isPublic = isPublic;
    e.currentTarget.closest(".field").querySelector("label").innerText = isPublic ? "Public" : "Private";
});

// When user uploads file, update the file name in the form
document.getElementById("new_collection_file").addEventListener("change", (e) => {
    const file = e.currentTarget.files[0];
    document.getElementById("new_collection_file_name").value = file.name;
});
