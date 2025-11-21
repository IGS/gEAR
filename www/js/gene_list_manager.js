"use strict";

import { apiCallsMixin, convertToFormData, copyToClipboard, createToast, getCurrentUser, getRootUrl, initCommonUI, logErrorInConsole, openModal, registerPageSpecificLoginUIUpdates } from "./common.v2.js?v=670b2ed";
import { GeneCart } from "./classes/genecart.v2.js?v=670b2ed";

let firstSearch = true;
let isAddFormOpen = false;
const resultsPerPage = 20;
let listView = "table";

let resultItems = [];   // array of ResultItem objects (with different views of a dataset)

// TODO - Add transformation code for quick gene list transformations

// Floating UI function alias. See https://floating-ui.com/docs/getting-started#umd
// Subject to change if we ever need these common names for other things.
const computePosition = window.FloatingUIDOM.computePosition;
const flip = window.FloatingUIDOM.flip;
const shift = window.FloatingUIDOM.shift;
const offset = window.FloatingUIDOM.offset;
const arrow = window.FloatingUIDOM.arrow;

class ResultItem {
    constructor(data) {

        // UI selectors
        this.tableResultsBody = document.querySelector("#results-table tbody");
        this.tableTemplate = document.getElementById("results-table-view");
        this.resultsListDiv = document.getElementById("results-list-div");
        this.listTemplate = document.getElementById("results-list-view");

        this.geneListId = data.id;
        this.gctype = data.gctype;
        this.gctypeLabel = null
        switch (data.gctype) {
            case "unweighted-list":
                this.gctypeLabel = "Unweighted";
                break;
            case "weighted-list":
                this.gctypeLabel = "Weighted";
                break;
            case "labeled-list":
                this.gctypeLabel = "Labeled";
                break;
            default:
                throw Error(`Invalid gene list type: ${data.gctype}`);
        }
        this.label = data.label;
        this.longDesc = data.ldesc || "";
        this.shareId = data.share_id;
        this.isPublic = Boolean(data.is_public);
        this.dateAdded = new Date(data.date_added).toDateString();
        // as YYYY/MM/DD
        this.shortDateAdded = new Date(data.date_added).toISOString().slice(0, 10);

        this.organismId = data.organism_id;
        this.geneCount = data.num_genes;
        this.userName = data.user_name;
        this.organism = data.organism;
        this.isOwner = data.is_owner;

    }

    createListItem() {
        const makeRandomString = (length) => {
            let result = '';
            const characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
            const charactersLength = characters.length;
            for (let i = 0; i < length; i++) {
                result += characters.charAt(Math.floor(Math.random() * charactersLength));
            }
            return result;
        }

        const geneListId = this.geneListId;

        // Clone the template
        const listItemView = this.listTemplate.content.cloneNode(true)

        // Set properties for multiple elements
        setElementProperties(listItemView, ".js-gc-list-element", { dataset: { geneListId } });
        // title section
        setElementProperties(listItemView, ".js-display-title p", { textContent: this.label });
        setElementProperties(listItemView, ".js-editable-title input", { value: this.label });
        // visibility/other metadata section
        const visibilityId = makeRandomString(10);
        setElementProperties(listItemView, ".js-editable-visibility input", { id: visibilityId, checked: this.isPublic });
        setElementProperties(listItemView, ".js-editable-visibility label", { htmlFor: visibilityId, textContent: this.isPublic ? "Public" : "Private" });
        // organism section
        setElementProperties(listItemView, ".js-display-organism span:last-of-type", { textContent: this.organism });
        setElementProperties(listItemView, ".js-editable-organism select", { value: this.organismId });
        // owner section
        setElementProperties(listItemView, ".js-display-owner span:last-of-type", { textContent: this.userName });
        setElementProperties(listItemView, ".js-editable-owner input", { value: this.userName });
        // date added section
        setElementProperties(listItemView, ".js-display-date-added span:last-of-type", { textContent: this.dateAdded });
        setElementProperties(listItemView, ".js-editable-date-added input", { value: this.dateAdded });
        // action buttons section
        setElementProperties(listItemView, ".js-view-gc", { value: this.shareId });
        setElementProperties(listItemView, ".js-delete-gc", { value: geneListId });
        setElementProperties(listItemView, ".js-edit-gc", { value: geneListId });
        setElementProperties(listItemView, ".js-edit-gc-save", { value: geneListId });
        setElementProperties(listItemView, ".js-edit-gc-cancel", { value: geneListId });
        // gene list type section
        setElementProperties(listItemView, ".js-display-gctype span:last-of-type", { textContent: this.gctypeLabel });
        setElementProperties(listItemView, ".js-editable-gctype input", { value: this.gctypeLabel });
        // long description section
        setElementProperties(listItemView, ".js-editable-ldesc textarea", { value: this.longDesc });
        return listItemView;
    }

    createListViewItem() {
        const listItemView = this.createListItem();

        // Append the cloned template to the results container
        this.resultsListDiv.appendChild(listItemView);
        this.resultListItem = this.resultsListDiv.querySelector(`.js-gc-list-element[data-gene-list-id="${this.geneListId}"]`);

        this.fixPreviewGenesButton(this.resultListItem);

        this.addListVisibilityInfo(this.resultListItem);
        this.addDescriptionInfo(this.resultListItem);

        this.createDeleteConfirmationPopover(this.resultListItem);
        this.createRenamePermalinkPopover(this.resultListItem);
    }

    createTableExpandRow() {
        const listItemView = this.createListItem();

        this.tableExpandedRow = this.rowItem.parentElement.querySelector(`[data-gene-list-id="${this.geneListId}"] + .js-table-row-expanded`);
        this.tableExpandedCell = this.tableExpandedRow.querySelector("td");

        // Append the cloned template to the expanded cell
        this.tableExpandedCell.appendChild(listItemView);
        this.expandedRowItem = this.tableExpandedCell.querySelector(`.js-gc-list-element[data-gene-list-id="${this.geneListId}"]`);

        this.fixPreviewGenesButton(this.expandedRowItem);

        this.addListVisibilityInfo(this.expandedRowItem);
        this.addDescriptionInfo(this.expandedRowItem);

        this.createDeleteConfirmationPopover(this.expandedRowItem);
        this.createRenamePermalinkPopover(this.expandedRowItem);
    }

    createTableRow() {
        const geneListId = this.geneListId;

        // Clone the template
        const rowItem = this.tableTemplate.content.cloneNode(true)

        // Adding dataset attrubute to be able to key in doing a querySelector action
        setElementProperties(rowItem, ".js-table-row", { dataset: { geneListId } });

        // Set properties for multiple elements
        setElementProperties(rowItem, ".js-display-title", { textContent: this.label });
        setElementProperties(rowItem, ".js-display-organism", { textContent: this.organism });
        setElementProperties(rowItem, ".js-display-owner", { textContent: this.userName });
        setElementProperties(rowItem, ".js-display-gctype", { textContent: this.gctypeLabel });
        setElementProperties(rowItem, ".js-display-num-genes", { textContent: this.geneCount });

        // Append the cloned template to the results container
        this.tableResultsBody.appendChild(rowItem);
        this.rowItem = document.querySelector(`.js-table-row[data-gene-list-id="${geneListId}"]`);
        this.addTableVisibilityInfo();

        // Add the expandable row
        this.createTableExpandRow();

        // Add event listneres
        this.addTableItemEventListeners();
    }

    addTableVisibilityInfo() {
        const tableVisibility = this.rowItem.querySelector(`.js-display-visibility`);
        if (this.isPublic) {
            tableVisibility.classList.add("has-background-primary-light");
            tableVisibility.textContent = "Public";
        } else {
            tableVisibility.classList.add("has-background-danger");
            tableVisibility.textContent = "Private";
        }
    }

    addTableItemEventListeners() {
        // click .js-td-expand in table row to expand/collapse the row
        this.rowItem.querySelector(".js-expand-row").addEventListener("click", (e) => {
            this.tableExpandedRow.classList.toggle("is-hidden");
            // Check if the row is already expanded
            if (this.tableExpandedRow.classList.contains("is-hidden")) {
                // State is expanded
                e.currentTarget.querySelector(".icon").innerHTML = '<i class="mdi mdi-24px mdi-chevron-down"></i>';
                return;
            }
            // Original state is collapsed
            e.currentTarget.querySelector(".icon").innerHTML = '<i class="mdi mdi-24px mdi-chevron-up"></i>';
        });
    }

    addListVisibilityInfo(parentElt) {
        // add dataset public/private info to DOM
        const geneListDisplayContainer = parentElt.querySelector(`.js-display-visibility`);
        const geneListDisplaySpan = document.createElement("span");
        geneListDisplaySpan.classList.add("tag");

        if (this.isPublic) {
            geneListDisplaySpan.classList.add("is-primary", "is-light");
            geneListDisplaySpan.textContent = "Public gene list";
        } else {
            geneListDisplaySpan.classList.add("is-danger");
            geneListDisplaySpan.textContent = "Private gene list";
        }
        geneListDisplayContainer.appendChild(geneListDisplaySpan);

        // Toggle switch (public is checked, private is unchecked)
        const visibilitySwitch = parentElt.querySelector(`.js-editable-visibility input`);
        visibilitySwitch.addEventListener("change", (e) => {
            const isPublic = e.currentTarget.checked;
            e.currentTarget.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";
        });
    }

    addDescriptionInfo(parentElt) {
        // Add ldesc if it exists
        const ldescText = parentElt.querySelector(".js-display-ldesc-text");
        ldescText.textContent = this.longDesc || "No description entered";
    }

    addListItemEventListeners(parentElt) {

        // Unweighted gene lists
        this.setupGeneListToggle(parentElt, ".js-gc-unweighted-gene-list-toggle", './cgi/get_unweighted_gene_cart_preview.cgi', createUnweightedGeneListTable);

        // Weighted gene lists
        this.setupGeneListToggle(parentElt, ".js-gc-weighted-gene-list-toggle", './cgi/get_weighted_gene_cart_preview.cgi', createWeightedGeneListPreview);

        // Labeled gene lists
        // NOT IMPLEMENTED YET

        // Expand and collapse gene list view
        parentElt.querySelector(".js-expand-box").addEventListener("click", (e) => {
            const expandableViewElts = parentElt.querySelectorAll(".js-expandable-view");
            for (const elt of expandableViewElts) {
                elt.classList.toggle("is-hidden");
            }

            // Toggle the icon
            if (e.currentTarget.innerHTML === '<i class="mdi mdi-arrow-expand"></i>') {
                e.currentTarget.innerHTML = '<i class="mdi mdi-arrow-collapse"></i>';
                return;
            }
            e.currentTarget.innerHTML = '<i class="mdi mdi-arrow-expand"></i>';
        });

        // Reformats the <ul> containing the gene symbols into a text file with one gene per row
        parentElt.querySelector(".js-download-gc").addEventListener("click", async (e) => {
            if (this.gctype == "unweighted-list") {
                const geneSymbols = await fetchGeneCartMembers(this.shareId);
                const fileContents = geneSymbols.map(gene => gene.label).join("\n");

                const element = document.createElement("a");
                element.setAttribute(
                    "href",
                    `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(fileContents)}`
                );

                element.setAttribute("download", `gene_cart.${this.shareId}.tsv`);
                element.style.display = "none";
                document.body.appendChild(element);
                element.click();
                document.body.removeChild(element);
                return;
            }

            if (this.gctype == "weighted-list") {
                // Download the source file (returns a content-disposition attachment header)

                const {data} = await axios.post("./cgi/download_weighted_gene_cart.cgi", convertToFormData({
                    'share_id': this.shareId
                }));

                const element = document.createElement("a");
                element.setAttribute(
                    "href",
                    `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(data)}`
                );

                element.setAttribute("download", `gene_cart.${this.shareId}.tsv`);
                element.style.display = "none";
                document.body.appendChild(element);
                element.click();
                document.body.removeChild(element);
                return;
            }
            if (this.gctype == "labeled-list") {
                // Retrieve the cart from the server and put it in a text file
                throw Error("Not implemented yet");

            } else {
                throw Error(`Invalid gene list type: ${this.gctype}`);
            }

        });

        parentElt.querySelector(".js-share-gc").addEventListener("click", (e) => {
            let currentPage = new URL(`${getRootUrl()}/p`);
            const params = new URLSearchParams(currentPage.search);

            params.set('c', this.shareId);

            // if genecart is weighted, can only link to projection.html
            if (this.gctypeLabel === "Weighted") {
                params.set('p', 'p');
            }
            currentPage.search = params.toString();
            const shareUrl = currentPage.toString();
            copyPermalink(shareUrl);
        });

        // Cancel button for editing a gene list
        parentElt.querySelector(".js-edit-gc-cancel").addEventListener("click", (e) => {
            // Show editable versions where there are some and hide the display versions
            for (const classElt of document.querySelectorAll(`.js-editable-version`)) {
                classElt.classList.add("is-hidden");
            };
            for (const classElt of document.querySelectorAll(`.js-display-version`)) {
                classElt.classList.remove("is-hidden");
            };

            // Reset any unsaved/edited values
            parentElt.querySelector(`.js-editable-visibility input`).value = this.isPublic ? "Public" : "Private";
            parentElt.querySelector(`.js-editable-title input`).value = this.title;
            parentElt.querySelector(`.js-editable-ldesc textarea`).value = this.longDesc;
            parentElt.querySelector(`.js-editable-organism select`).value = this.organismId;
            parentElt.querySelector(`.js-action-links`).classList.remove("is-hidden");
        });

        // Save button for editing a gene list
        parentElt.querySelector(".js-edit-gc-save").addEventListener("click", async (e) => {
            const newVisibility = parentElt.querySelector(`.js-editable-visibility input`).checked;
            // convert "true/false" visibility to 1/0
            const intNewVisibility = newVisibility ? 1 : 0;

            const newTitle = parentElt.querySelector(`.js-editable-title input`).value;
            const newLdesc = parentElt.querySelector(`.js-editable-ldesc textarea`).value;
            const newOrgId = parentElt.querySelector(`.js-editable-organism select`).value;
            const newOrgText = parentElt.querySelector(`.js-editable-organism select option[value='${newOrgId}']`).textContent;


            try {
                const data = await apiCallsMixin.saveGeneListInfoChanges(this.geneListId, intNewVisibility, newTitle, newOrgId, newLdesc);
                createToast("Gene list changes saved", "is-success");

            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to save gene list changes");
                return;
            } finally {
                parentElt.querySelector(`.js-action-links`).classList.remove("is-hidden");
            }

            this.isPublic = newVisibility;
            this.title = newTitle;
            this.longDesc = newLdesc;
            this.organismId = newOrgId;

            // Update the UI for the new values for the expanded row and the list item
            for (const selector of [this.expandedRowItem, this.resultListItem]) {

                if (newVisibility) {
                    selector.querySelector(`.js-display-visibility span`).textContent = "Public gene list";
                    selector.querySelector(`.js-display-visibility span`).classList.remove("is-danger");
                    selector.querySelector(`.js-display-visibility span`).classList.add("is-light", "is-primary");
                } else {
                    selector.querySelector(`.js-display-visibility span`).textContent = "Private gene list";
                    selector.querySelector(`.js-display-visibility span`).classList.remove("is-light", "is-primary");
                    selector.querySelector(`.js-display-visibility span`).classList.add("is-danger");
                }

                selector.querySelector(`.js-display-title p`).textContent = newTitle;
                selector.querySelector(`.js-display-ldesc-text`).textContent = newLdesc || "No description entered";

                selector.querySelector(`.js-display-organism span:last-of-type`).textContent = newOrgText;
            }


            // Update the UI for the table row values
            if (newVisibility) {
                this.rowItem.querySelector(`.js-display-visibility`).textContent = "Public";
                this.rowItem.querySelector(`.js-display-visibility`).classList.replace("has-background-danger", "has-background-primary-light");

            } else {
                this.rowItem.querySelector(`.js-display-visibility`).textContent = "Private";
                this.rowItem.querySelector(`.js-display-visibility`).classList.replace("has-background-primary-light", "has-background-danger");
            }

            this.rowItem.querySelector(`.js-display-title`).textContent = newTitle;
            this.rowItem.querySelector(`.js-display-organism`).textContent = this.expandedRowItem.querySelector(`.js-display-organism span:last-of-type`).textContent;

            // Put interface back to view mode.
            toggleEditableMode(true, parentElt);

        });

        // Toggle editable mode when edit button is clicked for a gene list
        const editSelector = parentElt.querySelector(".js-edit-gc")
        if (editSelector) {
            editSelector.addEventListener("click", async (e) => {
                const editableVisibilityElt = parentElt.querySelector(`.js-editable-visibility input`);

                editableVisibilityElt.checked = this.isPublic;
                editableVisibilityElt.closest(".field").querySelector("label").textContent = this.isPublic ? "Public" : "Private";

                // copy the organism selection list (from "create new gene list" form) for this row
                const editableOrganismIdElt = parentElt.querySelector(".js-editable-organism select");
                editableOrganismIdElt.innerHTML = document.getElementById("new-list-organism-id").innerHTML;

                // Remove the "select an organism" option
                editableOrganismIdElt.removeChild(editableOrganismIdElt.firstChild);

                // set the current value as selected
                editableOrganismIdElt.value = this.organismId;

                // Show editable versions where there are some and hide the display versions
                toggleEditableMode(false, parentElt);

                // Make sure the view is expanded
                const expandableViewElt = parentElt.querySelector(`.js-expandable-view`);
                if (expandableViewElt.classList.contains('is-hidden')) {
                    parentElt.querySelector(`span.js-expand-box`).click();
                }
                parentElt.querySelector(`.js-action-links`).classList.add("is-hidden");
            });
        }

        // Redirect to gene expression search
        parentElt.querySelector(".js-view-gc").addEventListener("click", (e) => {
            let currentPage = `${getRootUrl()}/p?`;

            // if genecart is weighted, can only link to projection.html
            if (this.gctypeLabel === "Weighted") {
                currentPage += "p=p&";
            }

            window.open(`${currentPage}c=${this.shareId}`, '_blank');
        });
    }

    // Adds gene list display elements to DOM.
    fixPreviewGenesButton(parentElt){
        const buttonElt = parentElt.querySelector(`.js-preview-genes-button-container button`);
        const buttonLabel = parentElt.querySelector(".js-gc-button-label");

        if (this.gctype === "weighted-list") {
            buttonElt.classList.add("js-gc-weighted-gene-list-toggle");

            buttonLabel.textContent = `Info`;
        } else if (this.gctype === "unweighted-list") {
            buttonElt.classList.add("js-gc-unweighted-gene-list-toggle");

            buttonLabel.textContent = `${this.geneCount} genes`;
        } else if (this.gctype === "labeled-list") {
            // Not implemented yet
            buttonElt.classList.add("js-gc-labeled-gene-list-toggle");
        } else {
            throw Error(`Invalid gene list type: ${this.gctype}`);
        }
        buttonLabel.dataset.offState = buttonLabel.textContent

    }


    //Sets up the gene list toggle functionality.
    setupGeneListToggle(parentElt, className, ajaxUrl, handleData) {

        // Gene List item will either be weighted or unweighted, not both.
        if (!parentElt.querySelector(className)) {
            return;
        }

        // Show genes when button is clicked & adjust button styling
        parentElt.querySelector(className).addEventListener("click", async (e) => {
            const button = e.currentTarget;
            const previewGenesContainer = parentElt.querySelector(".js-preview-genes-container");
            const buttonLabel = parentElt.querySelector(".js-gc-button-label");

            // Toggle gene container visibility
            if (previewGenesContainer.classList.contains("is-hidden")) {
                button.classList.remove("is-outlined");

                // If the preview table already exists, just show it
                if (parentElt.querySelector(`.js-info-container`)) {
                    previewGenesContainer.classList.remove("is-hidden");    // TODO: add animate CSS with fade in

                    button.querySelector("i").classList.add("mdi-eye-off");
                    button.querySelector("i").classList.remove("mdi-format-list-bulleted");
                    buttonLabel.textContent = "Hide";
                    return;
                }

                button.classList.add("is-loading");

                // Create the preview table of the first five genes
                try {
                    const {data} = await axios.post(ajaxUrl, convertToFormData({
                        'share_id': this.shareId
                    }));

                    if (data.success < 1) {
                        throw Error(data.message);
                    }

                    const infoContainer = document.createElement("div");
                    infoContainer.classList.add("js-info-container");
                    previewGenesContainer.appendChild(infoContainer);

                    // ? Should we re-add the preview table for weighted gene lists?
                    handleData(infoContainer, data);

                    previewGenesContainer.classList.remove("is-hidden");    // TODO: add animate CSS with fade in

                    button.querySelector("i").classList.add("mdi-eye-off");
                    button.querySelector("i").classList.remove("mdi-format-list-bulleted");
                    buttonLabel.textContent = "Hide";

                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to load gene list preview")
                    button.classList.add("is-outlined");
                } finally {
                    button.classList.remove("is-loading");
                }

                return;
            }
            previewGenesContainer.classList.add("is-hidden");
            button.classList.add("is-outlined");
            button.blur();
            button.querySelector("i").classList.remove("mdi-eye-off");
            button.querySelector("i").classList.add("mdi-format-list-bulleted");
            buttonLabel.textContent = buttonLabel.dataset.offState;

        });
    }

    updateGeneListButtons(parentElt) {
        //unhide all buttons
        for (const actionLinks of parentElt.querySelectorAll(".js-action-links .control")) {
            actionLinks.classList.remove("is-hidden");
        }

        // The ability to edit and delete and dataset are currently paired
        const deleteButton = parentElt.querySelector("button.js-delete-gc");
        const editButton = parentElt.querySelector("button.js-edit-gc");
        const editPermalinkButton = parentElt.querySelector("button.js-edit-gc-permalink");

        if (this.isOwner) {
            return;
        }

        // If user is not the owner of the dataset, remove the delete and edit buttons so user cannot manipulate
        // These will be regenerated when a search triggers processSearchResults
        deleteButton.parentElement.remove();    // remove .control element to prevent heavy line where button was
        editButton.parentElement.remove()
        editPermalinkButton.parentElement.remove()

        // Remove all editable elements to prevent editing in the DOM
        for (const editableElt of parentElt.querySelectorAll(`.js-editable-version`)) {
            editableElt.classList.remove()
        }
    }

    /**
     * Creates a confirmation popover for deleting a gene list.
     */
    createDeleteConfirmationPopover(parentElt) {
        parentElt.querySelector(".js-delete-gc").addEventListener('click', (e) => {
            const button = e.currentTarget;

            // remove existing popovers
            const existingPopover = document.getElementById('delete-gc-popover');
            if (existingPopover) {
                existingPopover.remove();
            }

            // Create popover content
            const popoverContent = document.createElement('article');
            popoverContent.id = 'delete-gc-popover';
            popoverContent.classList.add("message", "is-danger");
            popoverContent.setAttribute("role", "tooltip");
            popoverContent.style.width = "500px";
            popoverContent.innerHTML = `
                <div class='message-header'>
                    <p>Remove list</p>
                </div>
                <div class='message-body'>
                    <p>Are you sure you want to delete this gene list?</p>
                    <div class='field is-grouped' style='width:250px'>
                        <p class="control">
                            <button id='confirm-gc-delete' class='button is-danger'>Delete</button>
                        </p>
                        <p class="control">
                            <button id='cancel-gc-delete' class='button' value='cancel_delete'>Cancel</button>
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
                placement: 'bottom',
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

            // Add event listener to cancel button
            document.getElementById('cancel-gc-delete').addEventListener('click', () => {
                popoverContent.remove();
            });

            // Add event listener to confirm button
            document.getElementById('confirm-gc-delete').addEventListener('click', async (event) => {
                event.target.classList.add("is-loading");
                try {
                    const data = await apiCallsMixin.deleteGeneList(this.shareId);
                    if (data['success'] == 1) {

                        // ? Is this necessary if we are blowing away the object anyways
                        for (const selector of [this.expandedRowItem, this.resultListItem]) {
                            // Remove the dataset from the DOM
                            selector.remove();
                        }

                        createToast("Gene list deleted", "is-success");

                        // This can affect page counts, so we need to re-run the search
                        await submitSearch();

                    } else {
                        throw new Error(data['error']);
                    }
                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to delete gene list");
                } finally {
                    event.target.classList.remove("is-loading");
                    popoverContent.remove();
                }
            });
        });
    }

    /**
     * Creates a popover for renaming gene list permalink.
     */
    createRenamePermalinkPopover(parentElt){
        parentElt.querySelector(".js-edit-gc-permalink").addEventListener('click', (e) => {
            const button = e.currentTarget;

            // remove existing popovers
            const existingPopover = document.getElementById('rename-gc-link-popover');
            if (existingPopover) {
                existingPopover.remove();
            }

            let currentPage = `${getRootUrl()}/p?`;

            // if genecart is weighted, can only link to projection.html
            if (this.gctypeLabel === "Weighted") {
                currentPage += "p=p&";
            }


            // Create popover content
            const popoverContent = document.createElement('article');
            popoverContent.id = 'rename-gc-link-popover';
            popoverContent.classList.add("message", "is-primary");
            popoverContent.setAttribute("role", "tooltip");
            popoverContent.innerHTML = `
                <div class='message-header'>
                    <p>Rename gene list permalink</p>
                </div>
                <div class='message-body'>
                    <p>Please provide a new name for the gene list short-hand permalink.</p>
                    <div class='field has-addons'>
                        <div class='control'>
                            <a class="button is-static">
                                ${currentPage}c=
                            </a>
                        </div>
                        <div class='control'>
                            <input id='gc-link-name' class='input' type='text' placeholder='permalink' value=${this.shareId}>
                        </div>
                    </div>
                    <div class='field is-grouped' style='width:250px'>
                        <p class="control">
                            <button id='confirm-gc-link-rename' class='button is-primary' disabled>Update</button>
                        </p>
                        <p class="control">
                            <button id='cancel-gc-link-rename' class='button' value='cancel_rename'>Cancel</button>
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
                placement: 'bottom',
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

            document.getElementById("gc-link-name").addEventListener("keyup", () => {
                const newLinkName = document.getElementById("gc-link-name");
                const confirmRenameLink = document.getElementById("confirm-gc-link-rename");

                if (newLinkName.value.length === 0 || newLinkName.value === this.shareId) {
                    confirmRenameLink.disabled = true;
                    return;
                }
                confirmRenameLink.disabled = false;
            });

            // Add event listener to cancel button
            document.getElementById('cancel-gc-link-rename').addEventListener('click', () => {
                popoverContent.remove();
            });

            // Add event listener to confirm button
            document.getElementById('confirm-gc-link-rename').addEventListener('click', async (event) => {
                event.target.classList.add("is-loading");
                const newShareId = document.getElementById("gc-link-name").value;

                try {
                    const data = await apiCallsMixin.updateShareId(this.shareId, newShareId, "genecart");

                    if ((!data.success) || (data.success < 1)) {
                        const error = data.error || "Unknown error. Please contact gEAR support.";
                        throw new Error(error);
                    }

                    createToast("Gene list permalink renamed", "is-success");

                    this.shareId = newShareId;

                    popoverContent.remove();

                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to rename gene list permalink: " + error);
                } finally {
                    event.target.classList.remove("is-loading");
                }
            });
        });
    }


}

/**
 * Applies a tooltip to a reference element.
 *
 * @param {HTMLElement} referenceElement - The element to which the tooltip is applied.
 * @param {HTMLElement} tooltip - The tooltip element.
 * @param {string} [position="top"] - The preferred position of the tooltip relative to the reference element.
 */
const applyTooltip = (referenceElement, tooltip, position="top") => {

    const hideTooltip = (event) => {
        tooltip.classList.add("is-hidden");
    }

    const showTooltip = (event) => {
        // Compute position
        computePosition(event.currentTarget, tooltip, {
            placement: position, // Change this to your preferred placement
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
 * Clears the results views by hiding pagination, removing existing results, and removing the "no results" message if it exists.
 */
const clearResultsViews = () => {
    // hide pagination (but keep results were they are)
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.add("is-invisible");
    }

    // Clear any existing results
    const resultsListDiv = document.getElementById("results-list-div");
    for (const elt of resultsListDiv.querySelectorAll(":not(#results-list-view)")) {
        elt.remove()
    }

    const resultsTableBody = document.querySelector("#results-table tbody");
    for (const elt of resultsTableBody.querySelectorAll(":not(#results-table-view)")) {
        elt.remove()
    }

    // remove "no results" message if it exists
    if (document.getElementById("no-results-message")) {
        document.getElementById("no-results-message").remove();
    }
}

/**
 * Copies a permalink to the clipboard.
 *
 * @param {string} shareUrl - The sanitized URL to be copied to the clipboard.
 * @returns {void}
 */
const copyPermalink = (shareUrl) => {

    if(copyToClipboard(shareUrl)) {
        createToast("URL copied to clipboard", "is-info");
    } else {
        createToast(`Failed to copy to clipboard. URL: ${shareUrl}`);
    }

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
 * Creates an unweighted gene list table.
 *
 * @param {HTMLElement} infoContainer - The container element where the table will be appended.
 * @param {Object} data - The gene information data.
 * @param {Object} data.gene_info - The gene information object. The keys are Ensembl IDs and the values are gene symbols and descriptions.
 */
const createUnweightedGeneListTable = (infoContainer, data) => {
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
 * Creates a preview of a weighted gene list.
 *
 * @param {HTMLElement} infoContainer - The container element where the gene list preview will be appended.
 * @param {Object} data - The data object containing information about the gene list.
 * @param {number} data.num_genes - The number of genes in the list.
 * @param {string[]} data.weights - The names of the weights in the gene list.
 */
const createWeightedGeneListPreview = (infoContainer, data) => {
    // List the number of genes as well as the names of weights in the gene list
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
 * @param {string} geneCartShareId - The share ID of the gene cart.
 * @returns {Promise<Array<string>>} - A promise that resolves to an array of gene symbols.
 * @throws {Error} - If the gene list members cannot be fetched.
 */
const fetchGeneCartMembers = async (geneCartShareId) => {
    try {
        const {gene_symbols, success} = await apiCallsMixin.fetchGeneCartMembers(geneCartShareId);
        if (!success) {
            throw new Error("Could not fetch gene list members.");
        }
        return gene_symbols;
    } catch (error) {
        logErrorInConsole(error);
        createToast(error.message);
        throw error
    }
}

// Callbacks after attempting to save a gene list
const geneListFailure = (gc, message) => {
    logErrorInConsole(message);
    createToast(message);
}

const geneListSaved = async (gc) => {
    document.getElementById("create-new-gene-list").click(); // resets form also
    createToast("Gene list saved", "is-success");
    await submitSearch();
}

/**
 * Loads the organism list from the server and populates the organism choices and new list organism select elements.
 * @returns {Promise<void>} A promise that resolves when the organism list is loaded and elements are populated.
 */
const loadOrganismList = async () => {
    try {
        const data = await apiCallsMixin.fetchOrganismList();
        const organismChoices = document.getElementById("organism-choices");    // <ul> element
        for (const organism of data.organisms) {
            const li = document.createElement("li");
            li.dataset.dbval = organism.id;
            li.textContent = organism.label;
            organismChoices.appendChild(li);
        }
        const newListOrganismSelect = document.getElementById("new-list-organism-id");    // <select> element in "create new gene list" form
        newListOrganismSelect.innerHTML = "";

        for (const organism of data.organisms) {
            const option = document.createElement("option");
            option.value = organism.id;
            option.textContent = organism.label;
            newListOrganismSelect.appendChild(option);
        }
        // Add default "select an organism" option
        const defaultOption = document.createElement("option");
        defaultOption.value = "";
        defaultOption.textContent = "Select an organism";
        defaultOption.selected = true;
        newListOrganismSelect.prepend(defaultOption);

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
 */
const processSearchResults = (data) => {

    const resultsContainer = document.getElementById("results-container");

    // If there are no results, display a message
    if (data.gene_carts.length === 0) {
        const noResultsMessage = document.createElement("p");
        noResultsMessage.id = "no-results-message";

        noResultsMessage.className = "has-text-centered";
        noResultsMessage.textContent = "No results found.";
        resultsContainer.appendChild(noResultsMessage);
        return;
    }

    // data.gene_carts is a list of JSON strings (different from data.datasets in dataset explorer)
    for (const gcString of data.gene_carts) {
        const gc = JSON.parse(gcString);
        const resultItem = new ResultItem(gc);

        // TABLE VIEW AND EXPANDABLE VIEW
        resultItem.createTableRow();
        // LIST VIEW
        resultItem.createListViewItem();

        resultItems.push(resultItem);
    }

    setupGeneListItemActionTooltips();

    // Now that tooltips have been populated we can remove buttons and add event listeners
    for (const resultItem of resultItems) {
        for (const selector of [resultItem.expandedRowItem, resultItem.resultListItem]) {
            // Add event listeners for all gene list elements
            resultItem.addListItemEventListeners(selector);

            // Show/hide gene list buttons based on user ownership
            resultItem.updateGeneListButtons(selector);

        }
    }
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
    document.getElementById("btn-new-list-save").classList.remove("is-loading");

    document.getElementById("new-list-label").value = "";
    document.getElementById("new-list-ldesc").value = "";
    document.getElementById("new-list-pasted-genes").value = "";
    document.getElementById("new-list-file").value = "";
    document.getElementById("new-list-file-name").value = "";

    document.getElementById("new-list-visibility").checked = false;
    document.getElementById("new-list-visibility").closest(".field").querySelector("label").textContent = "Private";

    for (const classElt of document.getElementsByClassName("js-new-list-header")) {
        classElt.classList.remove('has-background-primary', 'has-text-white');
        classElt.classList.add("has-text-dark");
    }

    for (const classElt of document.getElementsByClassName("js-upload-gc-btn")) {
        classElt.removeAttribute('disabled');
    }

    // Currently this is not implemented yet
    document.getElementById("btn-gc-upload-labeled-list").setAttribute('disabled', 'disabled');

    document.getElementById("new-list-upload-c").classList.remove("is-hidden");

    document.getElementById("new-list-form-c").classList.add("is-hidden");
    document.getElementById("new-list-pasted-genes-c").classList.add("is-hidden");
    document.getElementById("new-list-file-name").classList.add("is-hidden");

    for (const classElt of document.querySelectorAll("#new-list-form-c .js-validation-help")) {
        classElt.remove();
    }

    // Remove all "is-danger" classes from form fields
    for (const classElt of document.querySelectorAll("#new-list-form-c .is-danger")) {
        classElt.classList.remove("is-danger");
    }

    isAddFormOpen = false;
}

/**
 * Sets properties on a target element selected by a selector within a given element.
 *
 * @param {HTMLElement} element - The parent element to query within.
 * @param {string} selector - The CSS selector to find the target element.
 * @param {Object} properties - An object containing the properties to set on the target element.
 * @param {Object} [properties.dataset] - An optional object containing data attributes to set on the target element.
 */
const setElementProperties = (element, selector, properties) => {
    const targetElement = element.querySelector(selector);
    if (targetElement) {
        for (const [key, value] of Object.entries(properties)) {
            if (key === 'dataset') {
                for (const [dataKey, dataValue] of Object.entries(value)) {
                    targetElement.dataset[dataKey] = dataValue;
                }
            } else {
                targetElement[key] = value;
            }
        }
    }
}

/**
 * Sets up tooltips for various elements on the dataset explorer page.
 *
 * This function performs the following steps:
 * 1. Removes existing tooltips from elements, except those with the class "js-collection-tooltip".
 * 2. Creates new tooltips for elements with the "data-tooltip-content" attribute within the action links group.
 * 3. Creates a tooltip for the table display group element with the "data-tooltip-content" attribute.
 * 4. Applies the created tooltips to the appropriate elements in the action links and table display groups.
 * 5. Creates and applies tooltips for dataset view buttons.
 */
const setupGeneListItemActionTooltips = () => {
    // Configure tooltips before manipulating action link buttons
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
    const tableDisplayGroupElt = document.querySelector(".js-td-expand");
    const tableDisplayGroupTooltip = createActionTooltips(tableDisplayGroupElt.querySelector("[data-tooltip-content]"))

    // Then apply each tooltip to the appropriate element for all elements with the data-tooltip-content attribute
    // NOTE: It's important to do this after creating all the tooltips so that the tooltips are created in the correct order

    for (const actionElt of document.getElementsByClassName("js-action-links")) {
        const loopTooltips = [...tooltips];
        for (const classElt of actionElt.querySelectorAll("[data-tooltip-content]")) {
            applyTooltip(classElt, loopTooltips.shift());
        }
    }
    for (const actionElt of document.getElementsByClassName("js-td-expand")) {
        for (const classElt of actionElt.querySelectorAll("[data-tooltip-content]")) {
            applyTooltip(classElt, tableDisplayGroupTooltip);
        }
    }

}

/**
 * Sets up the pagination UI based on the provided pagination data.
 *
 * @param {Object} pagination - The pagination data object.
 * @param {number} pagination.total_results - The total number of results.
 * @param {number} pagination.current_page - The current page number.
 * @param {number} pagination.total_pages - The total number of pages.
 * @param {number} resultsPerPage - The number of results per page.
 */
const setupPagination = (pagination) => {

        // Update result count and label
        document.getElementById("result-count").textContent = pagination.total_results;
        document.getElementById("result-label").textContent = pagination.total_results == 1 ? " result" : " results";
        document.getElementById("count-label-c").classList.remove("is-hidden");

        const firstResult = pagination.total_results > 0 ? (pagination.current_page - 1) * resultsPerPage + 1 : 0;
        const lastResult = Math.min(pagination.current_page * resultsPerPage, pagination.total_results);
        document.getElementById("result-range").textContent = `${firstResult} - ${lastResult}`;

        // Update pagination buttons
        for (const classElt of document.getElementsByClassName("pagination")) {
            classElt.replaceChildren();
            classElt.classList.remove("is-hidden");

            const paginationList = document.createElement("ul");
            paginationList.className = "pagination-list";
            classElt.appendChild(paginationList);


            // Add previous button
            if (pagination.current_page > 1) {
                paginationList.appendChild(createPaginationButton(null, 'left', async () => {
                    await submitSearch(pagination.current_page - 1);
                }));
            }

            // Add page buttons but only show 3 at a time (1 ... 3 4 5 ... 10)
            const startPage = Math.max(1, pagination.current_page - 1);
            const endPage = Math.min(pagination.total_pages, pagination.current_page + 1);

            if (startPage > 1) {
                paginationList.appendChild(createPaginationButton(1, null, async () => {
                    await submitSearch(1);
                }));
            }

            if (startPage > 2) {
                paginationList.appendChild(createPaginationEllipsis());
            }

            for (let i = startPage; i <= endPage; i++) {
                const li = paginationList.appendChild(createPaginationButton(i, null, async () => {
                    await submitSearch(i);
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
                paginationList.appendChild(createPaginationButton(pagination.total_pages, null, async () => {
                    await submitSearch(pagination.total_pages);
                }));
            }

            // Add next button
            if (pagination.current_page < pagination.total_pages) {
                paginationList.appendChild(createPaginationButton(null, 'right', async () => {
                    await submitSearch(pagination.current_page + 1);
                }));
            }
        }

}

/**
 * Submits a search request with the specified search terms and filters.
 *
 * @param {number} page - The page number of the search results.
 * @returns {Promise<void>} - A promise that resolves when the search results are processed.
 */
const submitSearch = async (page) => {

    // destroy current ResultItem objects
    resultItems = [];

    const searchTerms = document.getElementById("search-terms").value;

    // If this is the first time searching with terms, set the sort by to relevance
    if (searchTerms && firstSearch) {
        document.getElementById("sort-by").value = 'relevance';
        firstSearch = false;
    }

    const searchCriteria = {
        'session_id': getCurrentUser()?.session_id,
        'search_terms': searchTerms,
        'sort_by': document.getElementById("sort-by").value
    };

    // collect the filter options the user defined
    searchCriteria.organism_ids = buildFilterString('controls-organism');
    searchCriteria.date_added = buildFilterString('controls-date-added');
    searchCriteria.ownership = buildFilterString('controls-ownership');
    searchCriteria.limit = resultsPerPage;
    searchCriteria.page = page || 1;

    try {
        const data = await apiCallsMixin.fetchGeneLists(searchCriteria)

        // This is added here to prevent duplicate elements in the results generation if the user hits enter too quickly
        clearResultsViews();

        processSearchResults(data);
        setupPagination(data.pagination);
    } catch (error) {
        logErrorInConsole(error);
        createToast("Failed to search gene lists");
    } finally {
        // show pagination
        for (const classElt of document.getElementsByClassName("pagination")) {
            classElt.classList.remove("is-invisible");
        }
    }

    Cookies.set("default_gene_list_ownership_view", searchCriteria.ownership);
    Cookies.set("default_gene_list_organism_view", searchCriteria.organism_ids);
    Cookies.set("default_gene_list_date_added_view", searchCriteria.date_added);

    // restore previous list view
    if (listView === "table") {
        document.getElementById("btn-table-view").click();
    } else if (listView === "list-compact") {
        document.getElementById("btn-list-view-compact").click();
    } else if (listView === "list-expanded") {
        document.getElementById("btn-list-view-expanded").click();
    }
}

/**
 * Toggles the editable mode of elements with the given selector base.
 * @param {boolean} hideEditable - Whether to hide the editable elements or not.
 * @param {string} [target=""] - The base selector to use for finding the editable and non-editable elements.
 */
const toggleEditableMode = (hideEditable, target="") => {
    const editableElements = target.getElementsByClassName(`js-editable-version`);
    const nonEditableElements = target.getElementsByClassName(`js-display-version`);

    [...editableElements].forEach(el => el.classList.toggle("is-hidden", hideEditable));
    [...nonEditableElements].forEach(el => el.classList.toggle("is-hidden", !hideEditable));
}

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {
	// User settings has no "active" state for the sidebar
	document.getElementById("page-header-label").textContent = "Gene List Manager";

    const sessionId = getCurrentUser()?.session_id;

	if (! sessionId ) {
        document.getElementById("not-logged-in-msg").classList.remove("is-hidden");
        document.getElementById("create-new-gene-list").classList.add("is-hidden");
        document.getElementById("create-new-gene-list").setAttribute("disabled", "disabled");
        // only show public gene lists option
        for (const elt of document.querySelectorAll("#controls-ownership li:not([data-dbval='public'])")) {
            elt.remove();
        }
        document.querySelector("#controls-ownership li[data-dbval='public']").classList.add("js-selected");
    }

    await loadOrganismList();

    // Select the user's last remembered filter options
    const defaultOwnershipView = Cookies.get("default_gene_list_ownership_view");
    const defaultOrganismView = Cookies.get("default_gene_list_organism_view");
    const defaultDateAddedView = Cookies.get("default_gene_list_date_added_view");

    if (defaultOwnershipView && getCurrentUser()?.session_id) {
        // deselect All
        document.querySelector("#controls-ownership li.js-all-selector").classList.remove("js-selected");

        for (const ownership of defaultOwnershipView.split(",")) {
            document.querySelector(`#controls-ownership li[data-dbval='${ownership}']`).classList.add("js-selected");
        }
    }
    if (defaultOrganismView && getCurrentUser()?.session_id) {
        // deselect All
        document.querySelector("#controls-organism li.js-all-selector").classList.remove("js-selected");

        for (const organism of defaultOrganismView.split(",")) {
            document.querySelector(`#controls-organism li[data-dbval='${organism}']`).classList.add("js-selected");
        }
    }
    if (defaultDateAddedView && getCurrentUser()?.session_id) {
        document.querySelector(`#controls-date-added li[data-dbval='${defaultDateAddedView}']`).classList.add("js-selected");
    }

    await submitSearch();

    // Settings for selected facets
    for (const elt of document.querySelectorAll("ul.js-expandable-target li")) {
        elt.addEventListener("click", async (e) => {
            if (e.currentTarget.classList.contains("js-all-selector")) {
                // if the one clicked is the all_selector then highlight it and unclick the rest
                for (const elt of e.currentTarget.parentElement.children) {
                    elt.classList.remove("js-selected");
                }

                e.currentTarget.classList.add("js-selected");

            } else if (e.currentTarget.classList.contains("js-selected")) {
                // If turning off, make sure at least one other option is selected, else set "all" option
                e.currentTarget.classList.remove("js-selected");

                if (!e.currentTarget.parentElement.querySelectorAll("li.js-selected")) {
                    e.currentTarget.parentElement.querySelector("li.js-all-selector").classList.add("js-selected");
                }
            } else {
                // If turning on, make sure all_selector is off
                if (e.currentTarget.parentElement.querySelector("li.js-all-selector")) {
                    // In case not logged in and "All" is not an option
                    e.currentTarget.parentElement.querySelector("li.js-all-selector").classList.remove("js-selected");
                }

                // If this selection group has the 'only_one' option deselect the rest
                if (e.currentTarget.parentElement.classList.contains("js-choose-only-one")) {
                    for (const elt of e.currentTarget.parentElement.children) {
                        elt.classList.remove("js-selected");
                    }
                }

                e.currentTarget.classList.add("js-selected");
            }

            await submitSearch();
        });
    }
};
registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

// Pre-initialize some stuff
await initCommonUI();

// validate that #new-list-label input has a value
document.getElementById("new-list-label").addEventListener("blur", (e) => {
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

document.getElementById("search-clear").addEventListener("click", async () => {
    document.getElementById("search-terms").value = "";
    await submitSearch();
});

// Search for gene lists using the supplied search terms
const searchTermsElt = document.getElementById("search-terms");
searchTermsElt.addEventListener("keyup", async (event) => {
    const searchTerms = searchTermsElt.value;
    const searchClearElt = document.getElementById("search-clear");
    searchClearElt.classList.add("is-hidden");
    if (searchTerms) {
        searchClearElt.classList.remove("is-hidden");
    }
    if (event.key === "Enter") {
        await submitSearch();
    }
});

// Changing sort by criteria should update the search results
document.getElementById("sort-by").addEventListener("change", async () => {
    await submitSearch();
});

const btnCreateCartToggle = document.getElementById("create-new-gene-list");
btnCreateCartToggle.addEventListener("click", () => {
    const createListContainer = document.getElementById("create-list-container");
    const gcViewport = document.getElementById("results-container");
    const viewControls = document.getElementById("view-controls");
    const newCartVisibility = document.getElementById("new-list-visibility");
    const belowPagination = document.getElementById("below-pagination-container");

    if (createListContainer.classList.contains("is-hidden")) {
        createListContainer.classList.remove("is-hidden"); // TODO: Add animation with fade in/out
        gcViewport.classList.add("is-hidden");
        viewControls.classList.add("is-hidden");
        belowPagination.classList.add("is-hidden");
        btnCreateCartToggle.innerHTML =`<span class="icon"><i class="mdi mdi-undo-variant"></i></span> <span>Cancel new list</span>`;
        newCartVisibility.checked = false; // bootstrap toggle off
        newCartVisibility.closest(".field").querySelector("label").textContent = "Private";
        return;
    }
    createListContainer.classList.add("is-hidden");
    gcViewport.classList.remove("is-hidden");   // TODO: Add animation with fade in/out
    viewControls.classList.remove("is-hidden"); // TODO: Add animation with fade in/out
    belowPagination.classList.remove("is-hidden");
    btnCreateCartToggle.innerHTML =`<span class="icon"><i class="mdi mdi-plus"></i></span> <span>Create new gene list</span>`;
    resetAddForm();
});

// If an upload method is clicked, show the appropriate form and hide the others
document.getElementById("btn-gc-paste-unweighted-list").addEventListener("click", () => {
    document.getElementById("new-list-unweighted-header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new-list-unweighted-header").classList.remove("has-text-dark");

    document.getElementById("btn-gc-upload-unweighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-weighted-list").setAttribute('disabled', 'disabled');

    document.getElementById("new-list-form-c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out
    document.getElementById("new-list-pasted-genes-c").classList.remove("is-hidden");

    document.getElementById("new-list-upload-c").classList.add("is-hidden");

    document.getElementById("new-list-upload-type").value = "pasted_genes";
    isAddFormOpen = true;
});

document.getElementById("btn-gc-upload-unweighted-list").addEventListener("click", () => {
    if (isAddFormOpen) return;
    document.getElementById("new-list-unweighted-header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new-list-unweighted-header").classList.remove("has-text-dark");

    document.getElementById("btn-gc-paste-unweighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-weighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-labeled-list").setAttribute('disabled', 'disabled');

    document.getElementById("new-list-form-c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out

    document.getElementById("new-list-upload-type").value = "uploaded-unweighted";
    document.getElementById("new-list-upload-c").classList.remove("is-hidden");
    isAddFormOpen = true;
});

document.getElementById("btn-gc-upload-weighted-list").addEventListener("click", () => {
    if (isAddFormOpen) return;
    document.getElementById("new-list-weighted-header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new-list-weighted-header").classList.remove("has-text-dark");

    document.getElementById("btn-gc-paste-unweighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-unweighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-labeled-list").setAttribute('disabled', 'disabled');

    document.getElementById("new-list-form-c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out

    document.getElementById("new-list-upload-type").value = "uploaded-weighted";
    document.getElementById("new-list-upload-c").classList.remove("is-hidden");
    isAddFormOpen = true;
});

document.getElementById("btn-gc-upload-labeled-list").addEventListener("click", () => {
    if (isAddFormOpen) return;
    document.getElementById("new-list-labeled-header").classList.add('has-background-primary', "has-text-white");
    document.getElementById("new-list-labeled-header").classList.remove("has-text-dark");

    document.getElementById("btn-gc-paste-unweighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-unweighted-list").setAttribute('disabled', 'disabled');
    document.getElementById("btn-gc-upload-weighted-list").setAttribute('disabled', 'disabled');

    document.getElementById("new-list-form-c").classList.remove("is-hidden");   // TODO: Add animation with fade in/out

    document.getElementById("new-list-upload-type").value = "uploaded-labeled";
    document.getElementById("new-list-upload-c").classList.remove("is-hidden");
    isAddFormOpen = true;
});

// If the cancel button is clicked, hide the form and show the upload buttons
document.getElementById("btn-new-list-cancel").addEventListener("click", () => {
    document.getElementById("create-new-gene-list").click();
    document.getElementById("new-list-pasted-genes-c").classList.add("is-hidden");
});


const btnNewCartSave = document.getElementById("btn-new-list-save");
btnNewCartSave.addEventListener("click", (e) => {
    // disable button and show indicator that it's loading
    btnNewCartSave.classList.add("is-loading");

    let validationFailed = false;

    // check required fields
    const newCartLabel = document.getElementById("new-list-label");
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

    const newCartOrganism = document.getElementById("new-list-organism-id"); // <select> element in "create new gene list" form
    if (! newCartOrganism.value) {
        newCartOrganism.parentElement.classList.add("is-danger");
        // Add small helper text under input
        const newHelperElt = document.createElement("p");
        newHelperElt.classList.add("help", "has-text-danger-dark", "js-validation-help");
        newHelperElt.textContent = "Please select an organism";
        newCartOrganism.parentElement.appendChild(newHelperElt);

        validationFailed = true;
    }

    // Was file uploaded or genes pasted?
    const uploadType = document.getElementById("new-list-upload-type").value;
    if (uploadType === "pasted_genes") {
        const newCartPastedGenes = document.getElementById("new-list-pasted-genes");
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
        const newCartFile = document.getElementById("new-list-file");
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
        btnNewCartSave.classList.remove("is-loading");
        return;
    }

    // sleep for a second to show the loading indicator
    setTimeout(() => {
        // Remove all "is-danger" classes from form fields
        for (const classElt of document.querySelectorAll("#new-list-form-c .js-validation-help")) {
            classElt.remove();
        }

        // Remove all "is-danger" classes from form fields
        for (const classElt of document.querySelectorAll("#new-list-form-c .is-danger")) {
            classElt.classList.remove("is-danger");
        }

        // Passed to CGI script as 1 or 0
        const isPublic = document.getElementById("new-list-visibility").checked ? 1 : 0;

        const payload = {
            'new_cart_label': newCartLabel.value
            , 'new_cart_organism_id': newCartOrganism.value
            , 'new_cart_ldesc': document.getElementById("new-list-ldesc").value
            , 'is_public': isPublic
            , 'session_id': getCurrentUser()?.session_id
            , 'new_cart_upload_type': document.getElementById("new-list-upload-type").value
            , "new_cart_pasted_genes": document.getElementById("new-list-pasted-genes").value
            , "new_cart_file": document.getElementById("new-list-file").files[0]
        }

        const gc = new GeneCart();
        gc.addCartToDbFromForm(payload, geneListSaved, geneListFailure);
        btnNewCartSave.classList.remove("is-loading");
    }, 1000);
});

document.getElementById("btn-table-view").addEventListener("click", () => {
    listView = "table";
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-table-view").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-table-view").classList.remove('is-dark');

    document.getElementById("results-table").classList.remove("is-hidden");
    document.getElementById("results-list-div").classList.add("is-hidden");

})

document.getElementById("btn-list-view-compact").addEventListener("click", () => {
    listView = "list-compact";
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-list-view-compact").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-list-view-compact").classList.remove('is-dark');

    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.remove("is-hidden");

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

document.getElementById("btn-list-view-expanded").addEventListener("click", () => {
    listView = "list-expanded";
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-list-view-expanded").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-list-view-expanded").classList.remove('is-dark');

    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.remove("is-hidden");

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
document.getElementById("new-list-visibility").addEventListener("change", (e) => {
    const isPublic = e.currentTarget.checked;
    e.currentTarget.dataset.isPublic = isPublic;
    e.currentTarget.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";
});

// When user uploads file, update the file name in the form
document.getElementById("new-list-file").addEventListener("change", (e) => {
    const file = e.currentTarget.files[0];
    document.getElementById("new-list-file-name").textContent = file.name;
    document.getElementById("new-list-file-name").classList.remove("is-hidden");
});

// Show modal of upload requirements if clicked
document.getElementById("list-upload-reqs").addEventListener("click", (e) => {
    e.preventDefault();
    const modalId = e.target.dataset.target;
    const modalElt = document.getElementById(modalId);
    openModal(modalElt);

    const listType = document.getElementById("new-list-upload-type").value;

    if (listType === "uploaded-unweighted") {
        document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
            document.getElementById("unweighted-example").innerHTML;
    } else if (listType === "uploaded-weighted") {
        document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
            document.getElementById("weighted-example").innerHTML;
    } else if (listType === "uploaded-labeled") {
        document.querySelector(`#${modalId} .modal-card-body .content`).innerHTML =
            document.getElementById("labeled-example").innerHTML;

    } else {
        throw new Error("Unsupported file-upload list type");
    }

});
