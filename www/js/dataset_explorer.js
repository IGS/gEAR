"use strict";

import { apiCallsMixin, closeModal, createToast, getCurrentUser, getRootUrl, disableAndHideElement, enableAndShowElement, getUrlParameter, initCommonUI, logErrorInConsole, openModal, registerPageSpecificLoginUIUpdates } from "./common.v2.js?v=90ca195";
import { datasetCollectionState, fetchDatasetCollections, selectDatasetCollection } from "../include/dataset-collection-selector/dataset-collection-selector.js?v=90ca195";

// Pre-initialize some stuff
initCommonUI();

/* Imported variables
let datasetCollectionState.data; // from dataset-collection-selector
let datasetCollectionState.selectedShareId; // from dataset-collection-selector
*/

let firstSearch = true;
let searchByCollection = false;
let includePublicMembership = false;
const resultsPerPage = 20;
let listView = "table";

let flatDatasetCollectionData = {};   // flattened version of all dataset collections availabe to user

let resultItems = [];   // array of ResultItem objects (with different views of a dataset)

// TODO - Add transformation code for quick dataset transformations
// TODO - Add superuser edit abilities for "is_curator" users

// Floating UI function alias. See https://floating-ui.com/docs/getting-started#umd
// Subject to change if we ever need these common names for other things.
const computePosition = window.FloatingUIDOM.computePosition;
const flip = window.FloatingUIDOM.flip;
const shift = window.FloatingUIDOM.shift;
const offset = window.FloatingUIDOM.offset;
const arrow = window.FloatingUIDOM.arrow;

let singleArrangement;
let multiArrangement;

class ResultItem {
    constructor(data) {
        this.layouts = data.layouts;

        // UI selectors
        this.tableResultsBody = document.querySelector("#results-table tbody");
        this.tableTemplate = document.getElementById("results-table-view");
        this.resultsListDiv = document.getElementById("results-list-div");
        this.listTemplate = document.getElementById("results-list-view");

        this.datasetId = data.id;
        this.datasetType = data.dtype;
        this.label = data.title;

        this.longDesc = data.ldesc || "";
        this.shareId = data.share_id;
        this.isPublic = Boolean(data.is_public);
        this.isDownloadable = Boolean(data.is_downloadable);
        this.hasH5ad = Boolean(data.has_h5ad);
        this.dateAdded = new Date(data.date_added).toDateString();
        // as YYYY/MM/DD
        this.shortDateAdded = new Date(data.date_added).toISOString().slice(0, 10);

        this.userName = data.user_name;
        this.organism = data.organism;
        this.isOwner = data.is_owner

        this.annotSource = data.annotation_source || "Not given";
        this.annotVersion = data.annotation_release || "Not given";

        this.pubmedId = data.pubmed_id || null;
        this.geoId = data.geo_id || null;

        this.previewImageUrl = data.preview_image_url || "/img/dataset_previews/missing.png";

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

        const datasetId = this.datasetId;

        // Clone the template
        const listItemView = this.listTemplate.content.cloneNode(true)

        // Adding dataset attrubute to be able to key in doing a querySelector action
        setElementProperties(listItemView, ".js-dataset-list-element", { dataset: { datasetId } });

        // title section
        setElementProperties(listItemView, ".js-display-title p", { textContent: this.label });
        setElementProperties(listItemView, ".js-editable-title input", { value: this.label });
        // visibility/other metadata section
        const visibilityId = makeRandomString(10);
        setElementProperties(listItemView, ".js-editable-visibility input", { id: visibilityId, checked: this.isPublic });
        setElementProperties(listItemView, ".js-editable-visibility label", { htmlFor: visibilityId, textContent: this.isPublic ? "Public" : "Private" });

        // downloadable section
        const downloadableId = makeRandomString(10);
        setElementProperties(listItemView, ".js-editable-downloadable input", { id: downloadableId, checked: this.isDownloadable });
        setElementProperties(listItemView, ".js-editable-downloadable label", { htmlFor: downloadableId, textContent: this.isDownloadable ? "Yes" : "No" });
        // if h5ad is not present, then disable input and set to "No"
        if (!this.hasH5ad) {
            setElementProperties(listItemView, ".js-editable-downloadable input", { disabled: true });
            setElementProperties(listItemView, ".js-editable-downloadable label", { textContent: "No" });
        }

        // organism section
        setElementProperties(listItemView, ".js-display-organism span:last-of-type", { textContent: this.organism });
        setElementProperties(listItemView, ".js-editable-organism input", {value: this.organism });
        // owner section
        setElementProperties(listItemView, ".js-display-owner span:last-of-type", { textContent: this.userName });
        setElementProperties(listItemView, ".js-editable-owner input", { value: this.userName });
        // date added section
        setElementProperties(listItemView, ".js-display-date-added span:last-of-type", { textContent: this.dateAdded });
        setElementProperties(listItemView, ".js-editable-date-added input", { value: this.dateAdded });

        // annotation source section
        setElementProperties(listItemView, ".js-display-annot-source span:last-of-type", { textContent: this.annotSource });
        setElementProperties(listItemView, ".js-editable-annot-source input", { value: this.annotSource });

        // annotation version section
        setElementProperties(listItemView, ".js-display-annot-version span:last-of-type", { textContent: this.annotVersion });
        setElementProperties(listItemView, ".js-editable-annot-version input", { value: this.annotVersion });

        // pubmed id section
        const pubmedProp = this.pubmedId ? {
            // Not adding "icon-text" and "icon" classes because it aligns text and icon under the label
            innerHTML: `<a href="https://pubmed.ncbi.nlm.nih.gov/${this.pubmedId}" target="_blank">
                    <span>${this.pubmedId}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>`
        } : { textContent: "Not available" };

        setElementProperties(listItemView, ".js-display-pubmed-id span:last-of-type", pubmedProp);
        setElementProperties(listItemView, ".js-editable-pubmed-id input", { value: this.pubmedId });

        // geo id section
        const geoProp = this.geoId ? {
            innerHTML: `<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${this.geoId}" target="_blank">
                    <span>${this.geoId}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>`
        } : { textContent: "Not available" };

        setElementProperties(listItemView, ".js-display-geo-id span:last-of-type", geoProp);
        setElementProperties(listItemView, ".js-editable-geo-id input", { value: this.geoId });

        // action buttons section
        setElementProperties(listItemView, ".js-view-dataset", { value: this.shareId });
        setElementProperties(listItemView, ".js-view-projection-dataset", { value: this.shareId });
        setElementProperties(listItemView, ".js-delete-dataset", { value: datasetId });
        setElementProperties(listItemView, ".js-edit-dataset-permalink", { value: datasetId });

        setElementProperties(listItemView, ".js-share-dataset", { value: this.shareId });
        setElementProperties(listItemView, ".js-edit-dataset", { value: datasetId });
        setElementProperties(listItemView, ".js-edit-dataset-save", { value: datasetId });
        setElementProperties(listItemView, ".js-edit-dataset-cancel", { value: datasetId });

        setElementProperties(listItemView, ".js-dataset-curator", { href: `./dataset_curator.html?share_id=${this.shareId}`});
        setElementProperties(listItemView, ".js-multigene-viewer", { href: `./multigene_curator.html?share_id=${this.shareId}`});
        setElementProperties(listItemView, ".js-compare-tool", { href: `./compare_datasets.html?share_id=${this.shareId}`});
        setElementProperties(listItemView, ".js-sc-workbench", { href: `./sc_workbench.html?share_id=${this.shareId}`});


        // dataset type section
        setElementProperties(listItemView, ".js-display-dataset-type span:last-of-type", { textContent: this.datasetType });
        setElementProperties(listItemView, ".js-editable-dataset-type input", { value: this.datasetType });
        // long description section
        setElementProperties(listItemView, ".js-editable-ldesc textarea", { value: this.longDesc });

        return listItemView

    }

    createListViewItem() {
        const listItemView = this.createListItem();

        // Append the cloned template to the results container
        this.resultsListDiv.appendChild(listItemView);
        this.resultListItem = this.resultsListDiv.querySelector(`.js-dataset-list-element[data-dataset-id="${this.datasetId}"]`);

        this.createDeleteDatasetConfirmationPopover(this.resultListItem);
        this.createRenameDatasetPermalinkPopover(this.resultListItem);
    }

    createTableExpandRow() {
        const listItemView = this.createListItem();

        this.tableExpandedRow = this.rowItem.parentElement.querySelector(`[data-dataset-id="${this.datasetId}"] + .js-table-row-expanded`);
        this.tableExpandedCell = this.tableExpandedRow.querySelector("td");

        // Append the cloned template to the expanded cell
        this.tableExpandedCell.appendChild(listItemView);
        this.expandedRowItem = this.tableExpandedCell.querySelector(`.js-dataset-list-element[data-dataset-id="${this.datasetId}"]`);

        this.createDeleteDatasetConfirmationPopover(this.expandedRowItem);
        this.createRenameDatasetPermalinkPopover(this.expandedRowItem);
    }


    createTableRow() {
        const datasetId = this.datasetId;

        // Clone the template
        const rowItem = this.tableTemplate.content.cloneNode(true);


        // Adding dataset attrubute to be able to key in doing a querySelector action
        setElementProperties(rowItem, ".js-table-row", { dataset: { datasetId } });

        // TODO: Truncate titles with tooltip for full title
        // Set properties for multiple elements
        setElementProperties(rowItem, ".js-display-title", {textContent: this.label });
        setElementProperties(rowItem, ".js-display-organism", {textContent: this.organism });
        setElementProperties(rowItem, ".js-display-owner", { textContent: this.userName });
        setElementProperties(rowItem, ".js-display-date-added", { textContent: this.shortDateAdded });
        setElementProperties(rowItem, ".js-display-dataset-type", { textContent: this.datasetType });

        // Append the cloned template to the results container
        this.tableResultsBody.appendChild(rowItem);
        this.rowItem = document.querySelector(`.js-table-row[data-dataset-id="${datasetId}"]`);
        this.addTableVisibilityInfo();

        // Add the expandable row
        this.createTableExpandRow();

        // Add event listneres
        this.addTableItemEventListeners();

    }

    /**
     * Asynchronously renders and opens a modal for displaying dataset information.
     *
     * This function calls `renderDisplaysModal` with the dataset's ID, title, and
     * public status, then retrieves the modal element by its ID and opens it.
     *
     * @async
     * @function displaysModalCallback
     * @returns {Promise<void>} A promise that resolves when the modal has been rendered and opened.
     */
    async displaysModalCallback() {
        await renderDisplaysModal(this.datasetId, this.title, this.isPublic);
        const modalElt = document.getElementById(`displays-modal-${this.datasetId}`);
        openModal(modalElt);
    }

    // *** Table row stuff ***
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

        this.rowItem.querySelector(".js-view-displays").addEventListener("click", async (e) => {
            await this.displaysModalCallback();
        });

    }

    // *** List item stuff ***
    addCollectionMembershipInfo(parentElt) {
        parentElt.querySelector(`.js-found-in-collection-text`).innerHTML = "Found in these owned or highlighted collections";
        if (includePublicMembership) {
            parentElt.querySelector(`.js-found-in-collection-text`).innerHTML = "Found in these owned or public collections";
        }


        const collections = this.layouts;
        for (const collectionString of collections) {
            const collection = JSON.parse(collectionString);
            const shareId = collection.share_id;
            const label = collection.label;

            const collectionListElt = parentElt.querySelector(`.js-found-in-collections-list`);
            // Create list element that can link out to the collection
            const collectionListItem = document.createElement("li");
            collectionListItem.innerHTML = `
            <span>
                <a href="./p?l=${shareId}" target="_blank">
                    <span>${label}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>
            </span>
            `;
            collectionListElt.appendChild(collectionListItem);
        }
    }

    addDownloadableInfoToDataset(parentElt){
        const datasetDisplayContainer = parentElt.querySelector(`.js-display-downloadable`);
        const datasetDisplaySpan = document.createElement("span");
        datasetDisplaySpan.classList.add("tag");

        if (this.hasH5ad && this.isDownloadable) {
            datasetDisplaySpan.classList.add("is-success");
            datasetDisplaySpan.textContent = "Downloadable";
        } else {
            datasetDisplaySpan.classList.add("is-dark");
            datasetDisplaySpan.textContent = "Not downloadable";
        }
        datasetDisplayContainer.appendChild(datasetDisplaySpan);
    }

    addListVisibilityInfo(parentElt){
        // add dataset public/private info to DOM
        const datasetDisplayContainer = parentElt.querySelector(`.js-display-visibility`);
        const datasetDisplaySpan = document.createElement("span");
        datasetDisplaySpan.classList.add("tag");

        if (this.isPublic) {
            datasetDisplaySpan.classList.add("is-primary", "is-light");
            datasetDisplaySpan.textContent = "Public dataset";
        } else {
            datasetDisplaySpan.classList.add("is-danger");
            datasetDisplaySpan.textContent = "Private dataset";
        }
        datasetDisplayContainer.appendChild(datasetDisplaySpan);

        // Toggle switch (public is checked, private is unchecked)
        const visibilitySwitch = parentElt.querySelector(`.js-editable-visibility input`);
        visibilitySwitch.addEventListener("change", (e) => {
            const isPublic = e.currentTarget.checked;
            e.currentTarget.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";
        });

        const downloadableSwitch = parentElt.querySelector(`.js-editable-downloadable input`);

        downloadableSwitch.addEventListener("change", (e) => {
            const isDownloadable = e.currentTarget.checked;
            e.currentTarget.closest(".field").querySelector("label").textContent = isDownloadable ? "Yes" : "No";
        });
    }

    addImagePreview(parentElt) {
        const datasetImageContainer = parentElt.querySelector(`figure.is-square img`);
        datasetImageContainer.src = this.previewImageUrl;
    }

    addDescriptionInfo(parentElt) {
        // Add ldesc if it exists
        const ldescText = parentElt.querySelector(".js-display-ldesc-text");
        ldescText.textContent = this.longDesc || "No description entered";
    }

    addListItemEventListeners(parentElt) {

        // Add event listener to analysis dropdown trigger
        parentElt.querySelector(".js-analysis-dropdown .dropdown-trigger").addEventListener("click", (e) => {
            const item = e.currentTarget;
            item.closest(".dropdown").classList.toggle('is-active');
        });

        // Expand and collapse dataset view
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

        // Download button just needs href set. File will be <dataset_id>.h5ad
        const downloadSelector = parentElt.querySelector(".js-download-dataset");
        if (downloadSelector) {

            downloadSelector.addEventListener("click", async (e) => {
                e.currentTarget.classList.add("is-loading");
                try {
                    // download the h5ad
                    const datasetId = this.datasetId;
                    const url = `./cgi/download_source_file.cgi?type=h5ad&share_id=${this.shareId}`;
                    const a = document.createElement('a');
                    a.href = url;
                    a.click();
                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to download dataset");
                } finally {
                    e.currentTarget.classList.remove("is-loading");
                }
            });

        }

        parentElt.querySelector(".js-share-dataset").addEventListener("click", (e) => {

            let currentPage = new URL(`${getRootUrl()}/p`);
            const params = new URLSearchParams(currentPage.search);

            params.set('s', this.shareId);

            currentPage.search = params.toString();
            const shareUrl = currentPage.toString();
            copyPermalink(shareUrl);
        });

        parentElt.querySelector(".js-view-displays").addEventListener("click", async (e) => {
            await this.displaysModalCallback();
        });

        // Cancel button for editing a dataset
        parentElt.querySelector(".js-edit-dataset-cancel").addEventListener("click", (e) => {

            // Show editable versions where there are some and hide the display versions
            for (const classElt of parentElt.querySelectorAll(`.js-editable-version`)) {
                classElt.classList.add("is-hidden");
            };
            for (const classElt of parentElt.querySelectorAll(`.js-display-version`)) {
                classElt.classList.remove("is-hidden");
            };

            // Reset any unsaved/edited values
            parentElt.querySelector(`.js-editable-visibility input`).value = this.isPublic ? "Public" : "Private";
            parentElt.querySelector(`.js-editable-downloadable input`).value = this.isDownloadable ? "Yes" : "No";
            parentElt.querySelector(`.js-editable-title input`).value = this.title;
            parentElt.querySelector(`.js-editable-ldesc textarea`).value = this.longDesc;
            parentElt.querySelector(`.js-editable-pubmed-id input`).value = this.pubmedId;
            parentElt.querySelector(`.js-editable-geo-id input`).value = this.geoId;
            parentElt.querySelector(`.js-action-links`).classList.remove("is-hidden");

        });

        // Save button for editing a dataset
        parentElt.querySelector(".js-edit-dataset-save").addEventListener("click", async (e) => {
            const newVisibility = parentElt.querySelector(`.js-editable-visibility input`).checked;
            // convert "true/false" visibility to 1/0
            const intNewVisibility = newVisibility ? 1 : 0;

            const isDownloadable = parentElt.querySelector(`.js-editable-downloadable input`).checked;
            // convert "true/false" visibility to 1/0
            const intIsDownloadable = isDownloadable ? 1 : 0;

            const newTitle = parentElt.querySelector(`.js-editable-title input`).value;
            const newLdesc = parentElt.querySelector(`.js-editable-ldesc textarea`).value;
            const newPubmedId = parentElt.querySelector(`.js-editable-pubmed-id input`).value;
            const newGeoId = parentElt.querySelector(`.js-editable-geo-id input`).value;

            try {
                const data = await apiCallsMixin.saveDatasetInfoChanges(this.datasetId, intNewVisibility, intIsDownloadable, newTitle, newPubmedId, newGeoId, newLdesc);
                createToast("Dataset changes saved", "is-success");

            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to save dataset changes");
                return;
            } finally {
                parentElt.querySelector(`.js-action-links`).classList.remove("is-hidden");
            }

            this.isPublic = newVisibility;
            this.isDownloadable = isDownloadable;
            this.title = newTitle;
            this.longDesc = newLdesc;
            this.pubmedId = newPubmedId;
            this.geoId = newGeoId;

            // Update the UI for the new values for the expanded row and the list item
            for (const selector of [this.expandedRowItem, this.resultListItem]) {

                if (newVisibility) {
                    selector.querySelector(`.js-display-visibility span`).textContent = "Public dataset";
                    selector.querySelector(`.js-display-visibility span`).classList.remove("is-danger");
                    selector.querySelector(`.js-display-visibility span`).classList.add("is-light", "is-primary");
                } else {
                    selector.querySelector(`.js-display-visibility span`).textContent = "Private dataset";
                    selector.querySelector(`.js-display-visibility span`).classList.remove("is-light", "is-primary");
                    selector.querySelector(`.js-display-visibility span`).classList.add("is-danger");
                }


                if (isDownloadable) {
                    selector.querySelector(`.js-display-downloadable span`).textContent = "Downloadable";
                    selector.querySelector(`.js-display-downloadable span`).classList.remove("is-dark");
                    selector.querySelector(`.js-display-downloadable span`).classList.add("is-success");
                } else {
                    selector.querySelector(`.js-display-downloadable span`).textContent = "Not downloadable";
                    selector.querySelector(`.js-display-downloadable span`).classList.remove("is-success");
                    selector.querySelector(`.js-display-downloadable span`).classList.add("is-dark");
                }

                const downloadButton = selector.querySelector(`.js-download-dataset`);
                if (downloadButton) {
                    // If button exists (has h5ad), update the button visibility
                    disableAndHideElement(downloadButton, true);
                    if (isDownloadable) {
                        enableAndShowElement(downloadButton, true);
                    }
                }

                selector.querySelector(`.js-display-title p`).textContent = newTitle;

                selector.querySelector(`.js-display-ldesc-text`).textContent = newLdesc || "No description entered";

                // pubmed and geo display are links if they exist
                selector.querySelector(`.js-editable-pubmed-id input`).value = newPubmedId;
                if (newPubmedId) {
                    selector.querySelector(`.js-display-pubmed-id span:last-of-type`).innerHTML = `<a href="https://pubmed.ncbi.nlm.nih.gov/${newPubmedId}" target="_blank">
                        <span>${newPubmedId}</span>
                        <i class="mdi mdi-open-in-new"></i>
                    </a>`;
                } else {
                    selector.querySelector(`.js-display-pubmed-id span:last-of-type`).textContent = "Not available";
                }

                selector.querySelector(`.js-editable-geo-id input`).value = newGeoId;
                if (newGeoId) {
                    selector.querySelector(`.js-display-geo-id span:last-of-type`).innerHTML = `<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${newGeoId}" target="_blank">
                        <span>${newGeoId}</span>
                        <i class="mdi mdi-open-in-new"></i>
                    </a>`;
                } else {
                    selector.querySelector(`.js-display-geo-id span:last-of-type`).textContent = "Not available";
                }
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

            // Put interface back to view mode for current list item.
            toggleEditableMode(true, parentElt);

        });

        const editSelector = parentElt.querySelector(".js-edit-dataset");
        if (editSelector) {
            // Toggle editable mode when edit button is clicked for a dataset
            editSelector.addEventListener("click", async (e) => {
                const editableVisibilityElt = parentElt.querySelector(`.js-editable-visibility input`);

                editableVisibilityElt.checked = this.isPublic;
                editableVisibilityElt.closest(".field").querySelector("label").textContent = this.isPublic ? "Public" : "Private";

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
        parentElt.querySelector(".js-view-dataset").addEventListener("click", (e) => {
            window.open(`./p?s=${this.shareId}`, '_blank');
        });

        // Redirect to gene expression search
        parentElt.querySelector(".js-view-projection-dataset").addEventListener("click", (e) => {
            window.open(`./p?p=p&s=${this.shareId}`, '_blank');
        });
    }

    /**
     * Updates the dataset list buttons based on the dataset properties.
     */
    updateDatasetListButtons(parentElt) {
        //unhide all buttons
        for (const actionLinks of parentElt.querySelectorAll(".js-action-links .control")) {
            actionLinks.classList.remove("is-hidden");
        }

        // If the dataset has no h5ad file, remove the download button since it cannot be downloaded regardless
        const downloadButton = parentElt.querySelector("button.js-download-dataset");
        if (downloadButton && !this.hasH5ad) {
            downloadButton.parentElement.remove();
        }

        // If button still exists, update its visibility if the dataset is downloadable
        if (downloadButton) {
            disableAndHideElement(downloadButton, true);
            if (this.isDownloadable) {
                enableAndShowElement(downloadButton, true);
            }
        }

        // The ability to edit and delete and dataset are currently paired
        const deleteButton = parentElt.querySelector("button.js-delete-dataset");
        const editButton = parentElt.querySelector("button.js-edit-dataset");
        const editPermalinkButton = parentElt.querySelector("button.js-edit-dataset-permalink");

        if (this.isOwner) {
            return;
        }

        // If user is not the owner of the dataset, remove the delete and edit buttons so user cannot manipulate
        // These will be regenerated when a search triggers processSearchResults
        deleteButton.parentElement.remove(); // remove .control element to prevent heavy line where button was
        editButton.parentElement.remove()
        editPermalinkButton.parentElement.remove();

        // Remove all editable elements to prevent editing in the DOM
        for (const editableElt of parentElt.querySelectorAll(`.js-editable-version`)) {
            editableElt.classList.remove()
        }
    }

    /**
     * Creates a confirmation popover for deleting a dataset.
     */
    createDeleteDatasetConfirmationPopover(parentElt) {
        parentElt.querySelector(".js-delete-dataset").addEventListener('click', (e) => {
            const button = e.currentTarget;

            // remove existing popovers
            const existingPopover = document.getElementById('delete-dataset-popover');
            if (existingPopover) {
                existingPopover.remove();
            }

            // Create popover content
            const popoverContent = document.createElement('article');
            popoverContent.id = 'delete-dataset-popover';
            popoverContent.classList.add("message", "is-danger");
            popoverContent.setAttribute("role", "tooltip");
            popoverContent.style.width = "500px";
            popoverContent.innerHTML = `
                <div class='message-header'>
                    <p>Delete dataset</p>
                </div>
                <div class='message-body'>
                    <p>Are you sure you want to delete this dataset? This will affect any saved displays and dataset collections for yourself and for other users.</p>
                    <div class='field is-grouped' style='width:250px'>
                        <p class="control">
                            <button id='confirm-dataset-delete' class='button is-danger'>Delete</button>
                        </p>
                        <p class="control">
                            <button id='cancel-dataset-delete' class='button' value='cancel_delete'>Cancel</button>
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
            document.getElementById('cancel-dataset-delete').addEventListener('click', () => {
                popoverContent.remove();
            });

            // Add event listener to confirm button
            document.getElementById('confirm-dataset-delete').addEventListener('click', async (event) => {
                event.target.classList.add("is-loading");

                try {
                    const data = await apiCallsMixin.deleteDataset(this.datasetId);

                    if (data['success'] === 1) {

                        // ? Is this necessary if we are blowing away the object anyways
                        for (const selector of [this.expandedRowItem, this.resultListItem]) {
                            // Remove the dataset from the DOM
                            selector.remove();
                        }

                        createToast("Dataset deleted", "is-success");

                        // This can affect page counts, so we need to re-run the search
                        await submitSearch();

                    } else {
                        throw new Error(data['error']);
                    }
                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to delete dataset");
                } finally {
                    event.target.classList.remove("is-loading");
                    popoverContent.remove();
                }
            });
        });
    }

    /**
     * Creates a popover for renaming dataset permalink.
     */
    createRenameDatasetPermalinkPopover(parentElt) {
        parentElt.querySelector(".js-edit-dataset-permalink").addEventListener('click', (e) => {
            const button = e.currentTarget;

            // remove existing popovers
            const existingPopover = document.getElementById('rename-dataset-link-popover');
            if (existingPopover) {
                existingPopover.remove();
            }

            // Create popover content
            const popoverContent = document.createElement('article');
            popoverContent.id = 'rename-dataset-link-popover';
            popoverContent.classList.add("message", "is-primary");
            popoverContent.setAttribute("role", "tooltip");
            popoverContent.innerHTML = `
                <div class='message-header'>
                    <p>Rename dataset permalink</p>
                </div>
                <div class='message-body'>
                    <p>Please provide a new name for the dataset short-hand permalink.</p>
                    <div class='field has-addons'>
                        <div class='control'>
                            <a class="button is-static">
                                ${getRootUrl()}/p?s=
                            </a>
                        </div>
                        <div class='control'>
                            <input id='dataset-link-name' class='input' type='text' placeholder='permalink' value=${this.shareId}>
                        </div>
                    </div>
                    <div class='field is-grouped' style='width:250px'>
                        <p class="control">
                            <button id='confirm-dataset-link-rename' class='button is-primary' disabled>Update</button>
                        </p>
                        <p class="control">
                            <button id='cancel-dataset-link-rename' class='button' value='cancel_rename'>Cancel</button>
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

            document.getElementById("dataset-link-name").addEventListener("keyup", () => {
                const newLinkName = document.getElementById("dataset-link-name");
                const confirmRenameLink = document.getElementById("confirm-dataset-link-rename");

                if (newLinkName.value.length === 0 || newLinkName.value === this.shareId) {
                    confirmRenameLink.disabled = true;
                    return;
                }
                confirmRenameLink.disabled = false;
            });

            // Add event listener to cancel button
            document.getElementById('cancel-dataset-link-rename').addEventListener('click', () => {
                popoverContent.remove();
            });

            // Add event listener to confirm button
            document.getElementById('confirm-dataset-link-rename').addEventListener('click', async (event) => {
                event.target.classList.add("is-loading");
                const newShareId = document.getElementById("dataset-link-name").value;

                try {
                    const data = await apiCallsMixin.updateShareId(this.shareId, newShareId, "dataset");

                    if ((!data.success) || (data.success < 1)) {
                        const error = data.error || "Unknown error. Please contact gEAR support.";
                        throw new Error(error);
                    }

                    createToast("Dataset permalink renamed", "is-success");

                    // Update the share_id in the object, since the previous share_id is now invalid
                    this.shareId = newShareId;

                    popoverContent.remove();

                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to rename dataset permalink: " + error);
                } finally {
                    event.target.classList.remove("is-loading");
                }
            });
        });
    }
}

class LayoutArrangement {

    constructor(isMulti=false) {
        this.type = isMulti ? "multi" : "single";
        this.arrangement = [];
        this.arrangementDiv = document.getElementById(`dataset-arrangement-${this.type}`);

        this.arrangementWidth = 960; // NOTE: Originally used singleGeneArrangementDiv.offsetWidth, but it is 0 unless the element is visible
        this.rowWidth = this.arrangementWidth / 12; // Split width into 12 columns
        this.colHeight = this.rowWidth * 4; // A unit of height for us is 4 units of width (to make a square)
    }

    addMember(member) {
        this.arrangement.push(member);
    }

    /**
     * Sets up the adjustable arrangement grid for displaying tiles.
     *
     * - Calculates the maximum number of rows required based on the arrangement.
     * - Updates the CSS grid template rows for the arrangement container.
     * - Sorts the arrangement tiles by their starting row and column.
     * - Clears the arrangement container and appends each tile in order.
     * - Initializes interactable behavior for each tile.
     *
     * @returns {void}
     */
    setupArrangementAdjustable() {

        this.maxRow = Math.max(...this.arrangement.map((tile) => tile.startRow + tile.gridHeight)) -1;
        this.arrangementDiv.style.gridTemplateRows = `repeat(${this.maxRow}, ${this.colHeight}px)`;

        // Sort by row and column
        this.arrangement.sort((a, b) => {
            if (a.startRow === b.startRow) {
                return a.startCol - b.startCol;
            }
            return a.startRow - b.startRow;
        });

        // Clear the arrangement div
        this.arrangementDiv.innerHTML = "";

        // Add each member to the arrangement div
        for (const member of this.arrangement) {
            this.arrangementDiv.appendChild(member.createArrangementTile());
            member.createInteractable();
        }

    }
}

class LayoutArrangementMember {

    constructor(arrangementObj, displayId, gridPosition, startCol, startRow, gridWidth, gridHeight) {
        this.parentArrangement = arrangementObj;
        this.displayId = displayId;
        this.gridPosition = gridPosition;
        this.startCol = startCol;
        this.startRow = startRow;
        this.gridWidth = gridWidth;
        this.gridHeight = gridHeight;
        this.datasetTitle = "";
        this.image = "";
        this.tileTemplate = document.getElementById("dataset-arrangement-tile-template");
        this.selector = `.js-sortable-tile[data-display-id="${this.displayId}"]`;
        this.snapWidth = this.parentArrangement.rowWidth * 2; // Snap to 1/6 increments for width
        this.snapHeight = this.parentArrangement.colHeight; // Snap to 1/4 increments for height;
    }

    /**
     * Creates a new arrangement tile element by cloning the tile template,
     * setting its properties, grid area, title, and preview image.
     *
     * @returns {DocumentFragment} The newly created arrangement tile element.
     */
    createArrangementTile() {
        const tile = this.tileTemplate.content.cloneNode(true);

        setElementProperties(tile, ".js-sortable-tile", { dataset: {displayId: this.displayId} });
        // add gridAra to tile
        tile.querySelector(".js-sortable-tile").style.gridArea = `${this.startRow} / ${this.startCol} / span ${this.gridHeight} / span ${this.gridWidth}`;
        // add dataset title
        setElementProperties(tile, ".js-sortable-tile-title", { textContent: this.datasetTitle });
        // add preview image
        const tileImg = tile.querySelector(".js-sortable-tile-img");
        tileImg.src = this.image;
        return tile;
    }

    /**
     * Initializes interact.js draggable and resizable functionality on the element specified by `this.selector`.
     *
     * - Enables dragging within the parent container, snapping to a grid defined by `snapWidth` and `snapHeight`.
     * - Updates the CSS Grid layout properties (`gridArea`, `gridRowStart`, `gridColumnStart`, etc.) as the element is moved or resized.
     * - Dynamically adds new rows to the parent grid if the element is moved or resized beyond the current grid bounds.
     * - Calls `this.determineGridOverlap` at the end of drag or resize to handle overlap logic.
     *
     * @returns {void}
     */
    createInteractable() {
        const that = this;  // Preserve this context for event listeners

        this.interactable = interact(this.selector);

        this.interactable.draggable({
            // keep the element within the area of it's parent
            modifiers: [
                interact.modifiers.restrictRect({
                    restriction: 'parent',  // keep the drag within the parent. If dragged outside, it will snap back
                    endOnly: true
                }),
                interact.modifiers.snap({
                    // snap to the 12 possible grid columns, and n+1 grid rows.
                    targets: [
                        interact.createSnapGrid({ x: that.snapWidth, y: that.snapHeight }),
                    ],
                    range: Infinity,
                    relativePoints: [{ x: 1, y: 1 }],   // snap to the bottom right corner of the element (allows for a 1-to-1 with the grid start row/column)
                    offset: 'parent'
                })
            ],
            listeners: {
                move (event) {
                    const target = event.target

                    // get current style gridArea values (preserve span values)
                    const rowStart = parseInt(target.style.gridRowStart);
                    const colStart = parseInt(target.style.gridColumnStart);

                    // These are reported as "span N" in the gridArea style
                    const rowSpan = parseInt(target.style.gridRowEnd.split(" ")[1]);
                    const colSpan = parseInt(target.style.gridColumnEnd.split(" ")[1]);

                    const dStartCol = (Math.round(event.delta.x / that.snapWidth) * 2);  // x2 converts to 12 column grid
                    const dStartRow = Math.round(event.delta.y / that.snapHeight);

                    const newStartCol = colStart + dStartCol;
                    const newStartRow = rowStart + dStartRow;

                    target.style.gridArea = `${newStartRow} / ${newStartCol} / span ${rowSpan} / span ${colSpan}`;

                    // When dragged to the n+1 row, add a new row to the grid.
                    const arrangementDiv = target.parentElement;    // AKA this.parentArrangement.arrangementDiv
                    // get current number of rows in the grid
                    // it is stored as a string like "repeat(3, 100px)"
                    const arrangementRows = parseInt(arrangementDiv.style.gridTemplateRows.split(",")[0].split("(")[1]);

                    const lastTileRow = newStartRow + rowSpan -1;
                    // If the dragged tile's bottom edge is below the last row, add another row to the grid template
                    if (lastTileRow > arrangementRows) {
                        arrangementDiv.style.gridTemplateRows = `repeat(${lastTileRow}, ${that.parentArrangement.colHeight}px)`;
                    }
                },
                end: (event) => {
                    this.determineGridOverlap(event);
                }
            }

        }).resizable({
            // resize from only the right and bottom edges (adding top and left adds complexity)
            edges: { left: true, right: true, bottom: true, top: true },
            listeners: {
                move (event) {
                    const target = event.target;

                    // get current style gridArea values (at least the ones to preserve)
                    let rowStart = parseInt(target.style.gridRowStart);
                    let colStart = parseInt(target.style.gridColumnStart);
                    let rowSpan = parseInt(target.style.gridRowEnd.split(" ")[1]);
                    let colSpan = parseInt(target.style.gridColumnEnd.split(" ")[1]);

                    // Snap to grid in 1/6 increments for width
                    const newSpanWidth = (Math.round(event.rect.width / that.snapWidth) * 2);  // x2 multiplier converts to 12 column grid
                    const newSpanHeight = Math.round(event.rect.height / that.snapHeight);

                    // Calculate deltas for left/top resize
                    let deltaCol = 0;
                    let deltaRow = 0;
                    if (event.edges.left) {
                        deltaCol = colSpan - newSpanWidth;
                        colStart += deltaCol;
                    }
                    colSpan = newSpanWidth;
                    if (event.edges.top) {
                        deltaRow = rowSpan - newSpanHeight;
                        rowStart += deltaRow;
                    }
                    rowSpan = newSpanHeight;

                    // Clamp the start positions to ensure they are within the grid bounds
                    colStart = Math.max(1, colStart);
                    rowStart = Math.max(1, rowStart);

                    // Set new grid area for the tile
                    event.target.style.gridArea = `${rowStart} / ${colStart} / span ${newSpanHeight} / span ${newSpanWidth}`;

                    // When resized to the n+1 row, add a new row to the grid.
                    const arrangementDiv = target.parentElement;    // AKA this.parentArrangement.arrangementDiv
                    // get current number of rows in the grid
                    // it is stored as a string like "repeat(3, 100px)"
                    const arrangementRows = parseInt(arrangementDiv.style.gridTemplateRows.split(",")[0].split("(")[1]);

                    const lastTileRow = rowStart + newSpanHeight -1;
                    // If the dragged tile's bottom edge is below the last row, add another row to the grid template
                    if (lastTileRow > arrangementRows) {
                        arrangementDiv.style.gridTemplateRows = `repeat(${lastTileRow}, ${that.parentArrangement.colHeight}px)`;
                    }
                },
                end: (event) => {
                    this.determineGridOverlap(event);
                }
            }
        });
    }

    addOverlapState(tile) {
        // Add the overlap state for a specific tile
        tile.classList.add("js-is-overlapping");
        tile.style.opacity = 0.5;
        tile.querySelector(".js-sortable-tile-title").classList.add("has-text-danger-light", "has-background-danger-dark");
        tile.querySelector(".js-sortable-tile-title").classList.remove("has-background-primary-light");
    }

    removeOverlapState(tile) {
        // Clear the overlap state for a specific tile
        tile.classList.remove("js-is-overlapping");
        tile.style.opacity = 1;
        tile.querySelector(".js-sortable-tile-title").classList.remove("has-text-danger-light", "has-background-danger-dark");
        tile.querySelector(".js-sortable-tile-title").classList.add("has-background-primary-light");
    }

    /**
     * Determines if the currently moved/resized grid tile overlaps with any other tiles in the parent grid.
     * Provides a visual cue for overlapping tiles and disables the save button if any overlap is detected.
     *
     * @param {Event} event - The event object from the tile movement or resize action. Expects `event.target` to be a grid tile element with CSS grid properties.
     *
     * @returns {void}
     */
    determineGridOverlap(event) {

        // Get all tiles in the arrangement
        const arrangementDiv = event.target.parentElement;
        const arrangementTiles = arrangementDiv.querySelectorAll(".js-sortable-tile");

        // 1. Initialize num_overlaps for all tiles
        for (const tile of arrangementTiles) {
            tile.num_overlaps = 0;
        }

        // 2. Check all pairs for overlap
        // This is a brute-force O(n^2) check, but for a small number of tiles, it should be fine.
        for (let i = 0; i < arrangementTiles.length; i++) {
            for (let j = i + 1; j < arrangementTiles.length; j++) {
                const tileA = arrangementTiles[i];
                const tileB = arrangementTiles[j];

                // ...calculate grid positions for tileA and tileB...
                const tileARowStart = parseInt(tileA.style.gridRowStart);
                const tileARowEnd = tileARowStart + parseInt(tileA.style.gridRowEnd.split(" ")[1]);
                const tileAColStart = parseInt(tileA.style.gridColumnStart);
                const tileAColEnd = tileAColStart + parseInt(tileA.style.gridColumnEnd.split(" ")[1]);

                const tileBRowStart = parseInt(tileB.style.gridRowStart);
                const tileBRowEnd = tileBRowStart + parseInt(tileB.style.gridRowEnd.split(" ")[1]);
                const tileBColStart = parseInt(tileB.style.gridColumnStart);
                const tileBColEnd = tileBColStart + parseInt(tileB.style.gridColumnEnd.split(" ")[1]);

                // Check for overlap
                // This reminds me of the video game detection algorithms
                if (tileARowStart < tileBRowEnd
                    && tileARowEnd > tileBRowStart
                    && tileAColStart < tileBColEnd
                    && tileAColEnd > tileBColStart) {
                    tileA.num_overlaps++;
                    tileB.num_overlaps++;
                }
            }
        }

        // 3. Update overlap state and save button
        let anyOverlap = false;
        for (const tile of arrangementTiles) {
            if (tile.num_overlaps > 0) {
                this.addOverlapState(tile);
                anyOverlap = true;
            } else {
                this.removeOverlapState(tile);
            }
        }
        document.getElementById("btn-save-arrangement").disabled = anyOverlap;
    }
}

/**
 * Adds the dataset title to the modal.
 *
 * @param {HTMLElement} modalHTML - The HTML element representing the modal.
 */
const addDatasetTitleToModal = (modalHTML, title) => {
    const modalContent = modalHTML.querySelector('.modal-content');
    const datasetTitle = modalContent.querySelector("h5");
    datasetTitle.replaceChildren();
    datasetTitle.textContent = title;
}

const addModalEventListeners = () => {
    for (const classElt of document.getElementsByClassName("js-collection-add-display")) {
        classElt.addEventListener("click", async (e) => {
            const thisElt = e.currentTarget;

            // get nearest parent .js-modal-display
            const displayElement = thisElt.closest(".js-modal-display");
            const displayId = parseInt(displayElement.dataset.displayId);

            try {
                const data = await apiCallsMixin.addDisplayToCollection(datasetCollectionState.selectedShareId, displayId);
                if (!data.success) {
                    throw new Error(data.error);
                }

                // Update count of how many times the display is in the collection (.js-collection-display-count)
                const displayCount = displayElement.querySelector('.js-collection-display-count');
                const currentCount = parseInt(displayCount.textContent);
                displayCount.textContent = currentCount + 1;
                // Enable ".js-collection-remove-display" button
                const removeButton = displayElement.querySelector(".js-collection-remove-display");
                enableAndShowElement(removeButton);

                createToast("Display added to collection", "is-success");

                // Update the layout arrangement views
                await updateDatasetCollections();

            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to add dataset to collection");
            }
        });
    }

    for (const classElt of document.getElementsByClassName("js-collection-remove-display")) {
        classElt.addEventListener("click", async (e) => {
            const thisElt = e.currentTarget;
            // get nearest parent .js-modal-display
            const displayElement = thisElt.closest(".js-modal-display");
            const displayId = parseInt(displayElement.dataset.displayId);

            try {
                const data = await apiCallsMixin.deleteDisplayFromCollection(datasetCollectionState.selectedShareId, displayId);
                if (!data.success) {
                    throw new Error(data.error);
                }

                // Update count of how many times the display is in the collection (.js-collection-display-count)
                const displayCount = displayElement.querySelector('.js-collection-display-count');
                const currentCount = parseInt(displayCount.textContent);
                displayCount.textContent = currentCount - 1;
                // Disable ".js-collection-remove-display" button if this was the last instance of the display
                if (currentCount === 1) {
                    disableAndHideElement(thisElt);
                }

                createToast("Display removed from collection", "is-success");

                // Update the layout arrangement views
                await updateDatasetCollections();
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to remove dataset from collection");
            }
        });
    }
}

/**
 * Adds a section title to the specified element in a modal display.
 *
 * @param {HTMLElement} element - The element to which the section title will be added.
 * @param {string} titleText - The text content of the section title.
 */
const addModalDisplaySectionTitle = (element, titleText) => {
    if (!element.children.length) {
        return;
    }
    const title = document.createElement("p");
    title.classList.add("has-text-weight-bold", "is-underlined", "column", "is-full");
    title.textContent = titleText;
    element.prepend(title);
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
 * Clears the results views and hides pagination.
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
 * Creates a confirmation popover for deleting a dataset collection.
 * @returns {void}
 */
const createDeleteCollectionConfirmationPopover = () => {
    const button = document.getElementById("btn-delete-collection");
    button.addEventListener('click', (e) => {
        // remove existing popovers
        const existingPopover = document.getElementById('delete-collection-popover');
        if (existingPopover) {
            existingPopover.remove();
        }

        // Create popover content
        const popoverContent = document.createElement('article');
        popoverContent.id = 'delete-collection-popover';
        popoverContent.classList.add("message", "is-danger");
        popoverContent.setAttribute("role", "tooltip");
        popoverContent.innerHTML = `
            <div class='message-header'>
                <p>Delete collection</p>
            </div>
            <div class='message-body'>
                <p>Are you sure you want to delete this dataset collection?</p>
                <div class='field is-grouped' style='width:250px'>
                    <p class="control">
                        <button id='confirm-collection-delete' class='button is-danger'>Delete</button>
                    </p>
                    <p class="control">
                        <button id='cancel-collection-delete' class='button' value='cancel_delete'>Cancel</button>
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
        document.getElementById('cancel-collection-delete').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-collection-delete').addEventListener('click', async (event) => {
            event.target.classList.add("is-loading");
            try {
                const data = await apiCallsMixin.deleteDatasetCollection(datasetCollectionState.selectedShareId);

                if (data['success'] === 1) {
                    datasetCollectionState.selectedShareId = getCurrentUser().layout_share_id;

                    // This will trigger
                    // a) selectDatasetCollection
                    // b) datasetCollectionSelectorCallback
                    await updateDatasetCollections();
                    setActiveDCCategory("user");    // Update the user category to reflect the datasets without the deleted collection
                    createToast("Dataset collection deleted", "is-success");

                } else {
                    const error = data['error'] || "Failed to delete collection";
                    throw new Error(error);
                }
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to delete collection");
            } finally {
                event.target.classList.remove("is-loading");
                popoverContent.remove();

            }
        });
    });
}

/**
 * Creates a new collection popover.
 * This function attaches an event listener to the "Add Collection" button,
 * which, when clicked, creates a popover with a form to add a new dataset collection.
 * The popover is positioned relative to the button and includes input fields for collection name,
 * as well as buttons to confirm or cancel the addition of the collection.
 * Upon confirming the addition, the function makes an API call to create the dataset collection,
 * updates the UI, and displays a success message.
 * If an error occurs during the process, an error message is displayed.
 */
const createNewCollectionPopover = () => {
    const button = document.getElementById("btn-add-collection");
    button.addEventListener('click', (e) => {
        // remove existing popovers
        const existingPopover = document.getElementById('add-collection-popover');
        if (existingPopover) {
            existingPopover.remove();
        }

        // Create popover content
        const popoverContent = document.createElement('article');
        popoverContent.id = 'add-collection-popover';
        popoverContent.classList.add("message", "is-primary");
        popoverContent.setAttribute("role", "tooltip");
        popoverContent.innerHTML = `
            <div class='message-header'>
                <p>New collection</p>
            </div>
            <div class='message-body'>
                <p>Please provide a name for the new dataset collection</p>
                <div class='field'>
                    <div class='control'>
                        <input id='collection-name' class='input' type='text' placeholder='Collection name'>
                    </div>
                </div>
                <div class='field is-grouped' style='width:250px'>
                    <p class="control">
                        <button id='confirm-collection-add' class='button is-primary'>Add</button>
                    </p>
                    <p class="control">
                        <button id='cancel-collection-add' class='button' value='cancel_add'>Cancel</button>
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

        document.getElementById("collection-name").addEventListener("keyup", () => {
            const newCollectionName = document.getElementById("collection-name");
            const confirmAddCollection = document.getElementById("confirm-collection-add");

            if (newCollectionName.value.length === 0) {
                confirmAddCollection.disabled = true;
                return;
            }
            confirmAddCollection.disabled = false;
        });

        // Add event listener to cancel button
        document.getElementById('cancel-collection-add').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-collection-add').addEventListener('click', async (event) => {
            event.target.classList.add("is-loading");
            const newName = document.getElementById("collection-name").value;

            try {
                const data = await apiCallsMixin.createDatasetCollection(newName);

                if (data['layout_share_id']) {
                    datasetCollectionState.selectedShareId = data['layout_share_id'];
                    // This will trigger
                    // a) selectDatasetCollection
                    // b) datasetCollectionSelectorCallback
                    await updateDatasetCollections();
                    setActiveDCCategory("user");    // Update the user category to reflect the new collection item
                    createToast("Dataset collection created", "is-success");

                } else {
                    const error = data['error'] || "Failed to create new collection";
                    throw new Error(error);
                }
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to create new collection");
            } finally {
                event.target.classList.remove("is-loading");
                popoverContent.remove();
            }
        });
    });
}

/**
 * Creates a popover for renaming a dataset collection.
 * @returns {void}
 */
const createRenameCollectionPopover = () => {
    const button = document.getElementById("btn-rename-collection");
    button.addEventListener('click', (e) => {
        // remove existing popovers
        const existingPopover = document.getElementById('rename-collection-popover');
        if (existingPopover) {
            existingPopover.remove();
        }

        // Create popover content
        const popoverContent = document.createElement('article');
        popoverContent.id = 'rename-collection-popover';
        popoverContent.classList.add("message", "is-primary");
        popoverContent.setAttribute("role", "tooltip");
        popoverContent.innerHTML = `
            <div class='message-header'>
                <p>Rename collection</p>
            </div>
            <div class='message-body'>
                <p>Please provide a new name for the dataset collection</p>
                <div class='field'>
                    <div class='control'>
                        <input id='collection-name' class='input' type='text' placeholder='Collection name'>
                    </div>
                </div>
                <div class='field is-grouped' style='width:250px'>
                    <p class="control">
                        <button id='confirm-collection-rename' class='button is-primary' disabled>Update</button>
                    </p>
                    <p class="control">
                        <button id='cancel-collection-rename' class='button' value='cancel_rename'>Cancel</button>
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

        document.getElementById("collection-name").addEventListener("keyup", () => {
            const newCollectionName = document.getElementById("collection-name");
            const confirmRenameCollection = document.getElementById("confirm-collection-rename");

            if (newCollectionName.value.length === 0 || newCollectionName.value === datasetCollectionState.selectedLabel) {
                confirmRenameCollection.disabled = true;
                return;
            }
            confirmRenameCollection.disabled = false;
        });

        // Add event listener to cancel button
        document.getElementById('cancel-collection-rename').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-collection-rename').addEventListener('click', async (event) => {
            event.target.classList.add("is-loading");
            const newName = document.getElementById("collection-name").value;

            try {
                const data = await apiCallsMixin.renameDatasetCollection(datasetCollectionState.selectedShareId, newName);

                if (data['layout_label']) {
                    datasetCollectionState.selectedShareId = data['layout_share_id'];
                    // This will trigger
                    // a) selectDatasetCollection
                    // b) datasetCollectionSelectorCallback
                    await updateDatasetCollections();
                    setActiveDCCategory("user");    // Update the user category to reflect the new collection name
                    createToast("Dataset collection renamed", "is-success");
                } else {
                    const error = data['error'] || "Failed to rename collection";
                    throw new Error(error);
                }
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to rename collection");
            } finally {
                event.target.classList.remove("is-loading");
                popoverContent.remove();
            }
        });
    });
}

/**
 * Creates a popover for renaming the collection permalink.
 * @function createRenameCollectionPermalinkPopover
 * @returns {void}
 */
const createRenameCollectionPermalinkPopover = () => {
    const button = document.getElementById("btn-rename-collection-permalink");
    button.addEventListener('click', (e) => {
        // remove existing popovers
        const existingPopover = document.getElementById('rename-collection-link-popover');
        if (existingPopover) {
            existingPopover.remove();
        }

        // Create popover content
        const popoverContent = document.createElement('article');
        popoverContent.id = 'rename-collection-link-popover';
        popoverContent.classList.add("message", "is-primary");
        popoverContent.setAttribute("role", "tooltip");
        popoverContent.innerHTML = `
            <div class='message-header'>
                <p>Rename collection permalink</p>
            </div>
            <div class='message-body'>
                <p>Please provide a new name for the dataset collection short-hand permalink.</p>
                <div class='field has-addons'>
                    <div class='control'>
                        <a class="button is-static">
                            ${getRootUrl()}/p?l=
                        </a>
                    </div>
                    <div class='control'>
                        <input id='collection-link-name' class='input' type='text' placeholder='permalink' value=${datasetCollectionState.selectedShareId}>
                    </div>
                </div>
                <div class='field is-grouped' style='width:250px'>
                    <p class="control">
                        <button id='confirm-collection-link-rename' class='button is-primary' disabled>Update</button>
                    </p>
                    <p class="control">
                        <button id='cancel-collection-link-rename' class='button' value='cancel_rename'>Cancel</button>
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

        document.getElementById("collection-link-name").addEventListener("keyup", () => {
            const newLinkName = document.getElementById("collection-link-name");
            const confirmRenameLink = document.getElementById("confirm-collection-link-rename");

            if (newLinkName.value.length === 0 || newLinkName.value === datasetCollectionState.selectedShareId) {
                confirmRenameLink.disabled = true;
                return;
            }
            confirmRenameLink.disabled = false;
        });

        // Add event listener to cancel button
        document.getElementById('cancel-collection-link-rename').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-collection-link-rename').addEventListener('click', async (event) => {
            event.target.classList.add("is-loading");
            const newShareId = document.getElementById("collection-link-name").value;

            try {
                const data = await apiCallsMixin.updateShareId(datasetCollectionState.selectedShareId, newShareId, "layout");

                if ((!data.success) || (data.success < 1)) {
                    const error = data.error || "Unknown error. Please contact gEAR support.";
                    throw new Error(error);
                }

                await updateDatasetCollections();

                createToast("Dataset collection permalink renamed", "is-success");

                selectDatasetCollection(newShareId); // performs DatasetCollectionSelectorCallback when label is set
                popoverContent.remove();

            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to rename collection permalink: " + error);
            } finally {
                event.target.classList.remove("is-loading");
            }
        });
    });
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
 * Callback function for dataset collection selection.
 * This function updates the dataset list, shows action buttons, resets arrangement views,
 * finds the dataset collection data for the selected share_id, creates popovers for collection actions,
 * fetches collection members, renders layout arranger if the user owns the collection,
 * updates action buttons for the dataset collection or datasets, and hides/shows js-view-displays buttons.
 * @returns {Promise<void>} A promise that resolves when the function completes.
 */
const datasetCollectionSelectionCallback = async () => {
    // Update dataset list (in case user has "only in collection" toggle set)
    if (searchByCollection) {
        await submitSearch();
    }

    // Reset the arrangement views
    const arrangementViewSingle = document.getElementById("dataset-arrangement-single");
    const arrangementViewMulti = document.getElementById("dataset-arrangement-multi");
    arrangementViewSingle.innerHTML = "";
    arrangementViewMulti.innerHTML = "";

    // Get collection with displays
    const data = await apiCallsMixin.fetchDatasetCollectionMembers(datasetCollectionState.selectedShareId);
    document.getElementById("btn-arrangement-view").classList.add("is-hidden");
    // If user owns collection, show layout arranger
    if (data.is_owner) {
        await renderLayoutArranger(data);
        document.getElementById("btn-arrangement-view").classList.remove("is-hidden");
    } else {
        // restore previous list view since user should not be in arrangement view
        if (!searchByCollection) {
            if (listView === "table") {
                document.getElementById("btn-table-view").click();
            } else if (listView === "list-compact") {
                document.getElementById("btn-list-view-compact").click();
            } else if (listView === "list-expanded") {
                document.getElementById("btn-list-view-expanded").click();
            }
        }
    }

    // Update action buttons for the dataset collection or datasets
    updateDatasetCollectionButtons(data);

    // Hide js-view-displays buttons if user is not the owner of the collection
    // Also hide the "Displays" table column under the same condition
    updateViewDisplayAccess(data);

    // If the selected dataset collection is the current collection, make it look like the primary collection
    // ! Currently the selector will auto-make that collection the primary collection
    document.getElementById("btn-set-primary-collection").classList.add("is-outlined");
    if (datasetCollectionState.selectedShareId === getCurrentUser().layout_share_id) {
        document.getElementById("btn-set-primary-collection").classList.remove("is-outlined");
    }
}

/**
 * Retrieves all dataset displays for a given dataset ID.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @returns {Promise<{userDisplays: Array, ownerDisplays: Array}>} - An object containing arrays of user displays and owner displays.
 */
const getAllDatasetDisplays = async (datasetId) => {
    // Worth noting, both single and multigene displays are returned. Can add each to their respective layout arrangements
    const {user: userDisplays, owner: ownerDisplays} = await apiCallsMixin.fetchDatasetDisplays(datasetId);
    return { userDisplays, ownerDisplays };
}

/**
 * Initializes the dataset collection selection.
 *
 * @returns {void}
 */
const initializeDatasetCollectionSelection = () => {
    // Add mutation observer to watch if #dropdown-dc-selector-label changes
    const observer = new MutationObserver(async (mutationsList, observer) => {
        for (const mutation of mutationsList) {
            if (mutation.type === 'childList') {
                await datasetCollectionSelectionCallback();
            }
        }
    });

    observer.observe(document.getElementById("dropdown-dc-selector-label"), { childList: true });

    // Trigger the default dataset collection to be selected at the start
    if (getCurrentUser().layout_share_id) {
        selectDatasetCollection(getCurrentUser().layout_share_id);
    }

    // Show action buttons
    document.getElementById("collection-actions-c").classList.remove("is-hidden");

    // Create popvers for collection actions (reused when collection is changed)
    createDeleteCollectionConfirmationPopover();
    createNewCollectionPopover();
    createRenameCollectionPopover();
    createRenameCollectionPermalinkPopover();

}


/**
 * Loads the list of organisms from the server and populates the organism choices and new cart organism ID select elements.
 * @function
 * @returns {void}
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
    if (data.datasets.length === 0) {
        const noResultsMessage = document.createElement("p");
        noResultsMessage.id = "no-results-message";

        noResultsMessage.className = "has-text-centered";
        noResultsMessage.textContent = "No results found.";
        resultsContainer.appendChild(noResultsMessage);
        return;
    }

    // data.datasets is a list of JSON strings
    for (const dataset of data.datasets) {

        const resultItem = new ResultItem(dataset);

        // TABLE VIEW + EXPANDED ROW
        resultItem.createTableRow();

        // LIST VIEW
        resultItem.createListViewItem();

        for (const selector of [resultItem.expandedRowItem, resultItem.resultListItem]) {
            // Add collection membership info
            resultItem.addCollectionMembershipInfo(selector);

            // Add downloadable info
            resultItem.addDownloadableInfoToDataset(selector);

            // Add list visibility info
            resultItem.addListVisibilityInfo(selector);

            // Add image preview
            resultItem.addImagePreview(selector);

            // Add description info
            resultItem.addDescriptionInfo(selector);

        }

        resultItems.push(resultItem);

    }

    setupDatasetItemActionTooltips();

    // Normally this is done in the datasetCollectionCallback but we also need it when searching
    // to ensure the table-view button is shown/hid when filters are applied
    let collection = null;
    try {
        collection = flatDatasetCollectionData.find((collection) => collection.share_id === datasetCollectionState.selectedShareId);
    } catch (error) {
        // pass
    }

    updateViewDisplayAccess(collection);

    // Now that tooltips have been populated we can remove buttons and add event listeners
    for (const resultItem of resultItems) {
        for (const selector of [resultItem.expandedRowItem, resultItem.resultListItem]) {
            // Add event listeners for all dataset elements
            resultItem.addListItemEventListeners(selector);

            // Show/hide dataset list buttons based on user ownership
            resultItem.updateDatasetListButtons(selector);
        }
    }

}

/**
 * Removes existing modals from the document.
 */
const removeExistingModals = () => {
    const existingModals = document.querySelectorAll('.js-displays-modal');
    for (const modal of existingModals) {
        modal.remove();
    }

    // clear modal tooltips
    for (const tooltip of document.querySelectorAll(".js-modal-tooltip")) {
        tooltip.remove();
    }
}

/**
 * Renders the displays modal for a dataset.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} title - The title of the dataset.
 * @param {boolean} isPublic - Indicates whether the dataset is public or not.
 */
const renderDisplaysModal = async (datasetId, title, isPublic) => {
    removeExistingModals();

    // Clone the #tmpl-displays-modal template
    const modalTemplate = document.getElementById("tmpl-displays-modal");
    const modalHTML = modalTemplate.content.cloneNode(true);

    const modalDiv = modalHTML.querySelector('.modal');
    modalDiv.id = `displays-modal-${datasetId}`;
    modalDiv.dataset.datasetId = datasetId;

    addDatasetTitleToModal(modalHTML, title);

    const { userDisplays, ownerDisplays } = await getAllDatasetDisplays(datasetId);

    const modalContent = modalHTML.querySelector('.modal-content');
    // Add user and owner displays
    const userDisplaysElt = modalContent.querySelector(".js-modal-user-displays");
    userDisplaysElt.replaceChildren();
    const ownerDisplaysElt = modalContent.querySelector(".js-modal-owner-displays");
    ownerDisplaysElt.replaceChildren();

    const collection = flatDatasetCollectionData.find((collection) => collection.share_id === datasetCollectionState.selectedShareId);
    if (collection) {
        const layoutMemberData = await apiCallsMixin.fetchDatasetCollectionMembers(datasetCollectionState.selectedShareId);
        collection.members = layoutMemberData.layout_members.single.concat(layoutMemberData.layout_members.multi);
    }

    // Did it this way so we didn't have to pass async/await up the chain
    await Promise.allSettled([
        renderDisplaysModalDisplays(userDisplays, collection, userDisplaysElt, datasetId),
        renderDisplaysModalDisplays(ownerDisplays, collection, ownerDisplaysElt, datasetId)
    ]);

    // Add after displays are rendered
    addModalDisplaySectionTitle(userDisplaysElt, "Your Displays");
    addModalDisplaySectionTitle(ownerDisplaysElt, "Displays by Dataset Owner");

    // Close button event listener
    const closeButton = modalDiv.querySelector(".modal-close");
    closeButton.addEventListener("click", (event) => {
        closeModal(modalDiv);
    });
    const modalBackground = modalDiv.querySelector(".modal-background");
    modalBackground.addEventListener("click", (event) => {
        closeModal(modalDiv);
    });

    // If there are no displays, display a message
    if (userDisplays.length === 0 && ownerDisplays.length === 0) {
        const noDisplaysMessage = document.createElement("p");
        noDisplaysMessage.className = "has-text-centered";
        noDisplaysMessage.textContent = "No displays found for this dataset.";
        const modalContentBox = modalContent.querySelector(".box");
        modalContentBox.append(noDisplaysMessage);
        document.body.append(modalHTML);
        return;
    }

    // Add modal to DOM
    document.body.append(modalHTML);

    // Set state of initial display add/remove buttons based on the initial dataset collection.
    updateDisplayAddRemoveToCollectionButtons(modalDiv.id, collection);

    addModalEventListeners();

    // Create tooltips for all elements with the data-tooltip-content attribute
    // Only creating one set so that they can be reused
    const actionGroupElt = document.querySelector(".js-collection-add-remove-group");
    const tooltips = []

    if (actionGroupElt) {
        // Push into master list, the first instance of tooltips
        for (const classElt of actionGroupElt.querySelectorAll("[data-tooltip-content]")) {
            tooltips.push(createActionTooltips(classElt))
        }

        for (const tooltip of tooltips) {
            tooltip.classList.add("js-modal-tooltip"); // prevent removal of tooltip
        }

        // Then apply each tooltip to the appropriate element for all elements with the data-tooltip-content attribute
        for (const actionElt of document.querySelectorAll(".js-collection-add-remove-group")) {
            const loopTooltips = [...tooltips];
            for (const classElt of actionElt.querySelectorAll("[data-tooltip-content]")) {
                applyTooltip(classElt, loopTooltips.shift());
            }
        }
    }


    // Render warning about public collection and private dataset and disable "add to collection" button
    if (!(collection.is_public && !isPublic)) {
        return;
    }
    const warningElt = document.createElement("p");
    warningElt.className.add("has-background-warning-light", "has-text-warning-dark", "has-text-centered", "py-2", "px-3", "mb-3");
    warningElt.textContent = "This collection is public, but the dataset is private. You cannot add displays from private datasets to public collections.";
    modalContent.prepend(warningElt);
    disableAndHideElement(modalContent.querySelector(".js-add-to-collection"));
    modalContent.querySelector(".js-collection-add-remove-group").classList.add("is-hidden");
}

/**
 * Renders the displays in a modal window.
 *
 * @param {Array} displays - The array of displays to render.
 * @param {Object} collection - The collection object.
 * @param {HTMLElement} displayElt - The element where the displays will be appended.
 * @param {string} datasetId - The ID of the dataset.
 * @returns {Promise<void>} - A promise that resolves when the displays are rendered.
 */
const renderDisplaysModalDisplays = async (displays, collection, displayElt, datasetId) => {
    const displayTemplate = document.getElementById('tmpl-displays-modal-display');

    // Add user displays
    for (const display of displays) {
        const displayHTML = displayTemplate.content.cloneNode(true);

        const displayId = display.id;

        const displayElement = displayHTML.querySelector('.js-modal-display');
        displayElement.dataset.displayId = displayId;
        displayElement.dataset.datasetId = datasetId;

        // Add display image
        let displayUrl = "";
        try {
            // NOTE: SVGs are colorless
            displayUrl = await apiCallsMixin.fetchDatasetDisplayImage(datasetId, displayId);
        } catch (error) {
            logErrorInConsole(error);
            // Realistically we should try to plot, but I assume most saved displays will have an image present.
            displayUrl = "/img/dataset_previews/missing.png";
            if (display.plot_type === "epiviz") {
                displayUrl = "/img/epiviz_mini_screenshot.jpg"; // TODO: Replace with real logo
            }
        }

        const displayImage = displayElement.querySelector('figure > img');
        displayImage.src = displayUrl;

        // Add tag indicating plot type
        const displayType = displayElement.querySelector('.js-modal-display-type');

        const multiGeneDisplay = ["heatmap", "dotplot", "mg_violin", "volcano", "quadrant"];

        // Add color tags to displayType depending on plot type
        if (multiGeneDisplay.includes(display.plot_type)) {
            displayType.classList.add("is-danger");
            displayType.textContent = `${display.plot_type} [multigene]`;

        } else {
            displayType.classList.add("is-info");
            displayType.textContent = display.plot_type;
        }

        // Determine number of times display is in current layout
        const displayCount = displayElement.querySelector('.js-collection-display-count');

        if (collection?.members) {
            const displayCountValue = collection.members.filter((member) => JSON.parse(member).display_id === displayId).length
            displayCount.textContent = displayCountValue;
        } else {
            // maybe a new collection
            displayCount.textContent = 0;
        }

        // Append display to modal
        displayElt.append(displayHTML);
    }
}

/**
 * Renders the layout arranger for a given collection.
 *
 * @param {Object} collection - The collection object containing layout members.
 * @returns {Promise<void>} - A promise that resolves when the layout arranger is rendered.
 */
const renderLayoutArranger = async (collection) => {
    // Loading indication
    document.getElementById("dataset-arrangement-loading-notification").classList.remove("is-hidden");

    // The share_id should be updated in the component when a new dataset collection is selected
    const datasetData = await apiCallsMixin.fetchDatasets({layout_share_id: datasetCollectionState.selectedShareId, sort_by: "date_added"})

    // Get the titles of the datasets
    const titles = {};
    for (const dataset of datasetData.datasets) {
        titles[dataset.dataset_id] = dataset.title;
    }

    // JSON parse every layout member
    const layoutMembers = collection.layout_members;
    const singleLayoutMembers = layoutMembers.single || [];
    const multiLayoutMembers = layoutMembers.multi || [];

    singleArrangement = new LayoutArrangement();
    multiArrangement = new LayoutArrangement(true);

    // Single-gene displays

    document.getElementById("dataset-arrangement-single").classList.remove("is-hidden");
    document.getElementById("dataset-arrangement-single-no-displays-notification").classList.add("is-hidden");
    if (!singleLayoutMembers.length) {
        document.getElementById("dataset-arrangement-single").classList.add("is-hidden");
        document.getElementById("dataset-arrangement-single-no-displays-notification").classList.remove("is-hidden");
    }

    for (const display of singleLayoutMembers) {
        const member = JSON.parse(display);
        const displayId = member.display_id;
        const datasetId = member.dataset_id;

        const singleMember = new LayoutArrangementMember(singleArrangement, displayId, member.grid_position, member.start_col, member.start_row, member.grid_width, member.grid_height);

        singleMember.image = await apiCallsMixin.fetchDatasetDisplayImage(datasetId, displayId)

        singleMember.datasetTitle = titles[datasetId];
        singleArrangement.addMember(singleMember);
    }

    // Multi-gene displays

    document.getElementById("dataset-arrangement-multi").classList.remove("is-hidden");
    document.getElementById("dataset-arrangement-multi-no-displays-notification").classList.add("is-hidden");
    if (!multiLayoutMembers.length) {
        document.getElementById("dataset-arrangement-multi").classList.add("is-hidden");
        document.getElementById("dataset-arrangement-multi-no-displays-notification").classList.remove("is-hidden");
    }

    for (const display of multiLayoutMembers) {
        const member = JSON.parse(display);
        const displayId = member.display_id;
        const datasetId = member.dataset_id;

        const multiMember = new LayoutArrangementMember(multiArrangement, displayId, member.grid_position, member.start_col, member.start_row, member.grid_width, member.grid_height);

        multiMember.image = await apiCallsMixin.fetchDatasetDisplayImage(datasetId, displayId)

        multiMember.datasetTitle = titles[datasetId];
        multiArrangement.addMember(multiMember);

    }

    singleArrangement.setupArrangementAdjustable();
    multiArrangement.setupArrangementAdjustable();

    // Hide loading indication
    document.getElementById("dataset-arrangement-loading-notification").classList.add("is-hidden");
}


/**
 * Sets properties on a target element found within a given element using a selector.
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
const setupDatasetItemActionTooltips = () => {
    // Configure tooltips before manipulating action link buttons
    for (const tooltipElt of document.getElementsByClassName("tooltip")) {
        // Do not remove collection tooltips
        if (tooltipElt.classList.contains("js-collection-tooltip")) {
            continue;
        }
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

    // Create tooltips for dataset view buttons
    const viewBtns = document.getElementsByClassName("js-view-btn");
    for (const classElt of viewBtns) {
        applyTooltip(classElt, createActionTooltips(classElt), "bottom");
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
        document.getElementById("result-label").textContent = pagination.total_results === 1 ? " result" : " results";
        document.getElementById("count-label-c").classList.remove("is-hidden");

        const firstResult = pagination.total_results > 0 ? (pagination.current_page - 1) * resultsPerPage + 1 : 0;
        const lastResult = Math.min(pagination.current_page * resultsPerPage, pagination.total_results);
        document.getElementById("result-range").textContent = `${firstResult} - ${lastResult}`;

        // Update pagination buttons
        for (const classElt of document.getElementsByClassName("pagination")) {
            classElt.replaceChildren();
            classElt.classList.remove("is-invisible");

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
                if (i === pagination.current_page) {
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
 * Submits a search request with the specified search terms and pagination options.
 *
 * @param {number} [page=1] - The page number to retrieve.
 * @returns {Promise<void>} - A promise that resolves when the search request is completed.
 */
const submitSearch = async (page=1) => {

    // destroy current ResultItem objects
    resultItems = [];

    const searchTerms = document.getElementById("search-terms").value;

    // If this is the first time searching with terms, set the sort by to relevance
    if (searchTerms && firstSearch) {
        document.getElementById("sort-by").value = 'relevance';
        firstSearch = false;
    }

    const searchCriteria = {
        'session_id': getCurrentUser().session_id,
        'search_terms': searchTerms,
        'sort_by': document.getElementById("sort-by").value
    };

    if (searchByCollection) {
        searchCriteria.layout_share_id = datasetCollectionState.selectedShareId;
    }

    if (includePublicMembership) {
        searchCriteria.include_public_collection_membership = true;
    }

    // collect the filter options the user defined
    searchCriteria.organism_ids = buildFilterString('controls-organism');
    searchCriteria.date_added = buildFilterString('controls-date-added');
    searchCriteria.ownership = buildFilterString('controls-ownership');
    searchCriteria.dtypes = buildFilterString('controls-dataset-type');
    searchCriteria.limit = resultsPerPage;
    searchCriteria.page = page;

    try {
        const data = await apiCallsMixin.fetchDatasets(searchCriteria)
        // This is added here to prevent duplicate elements in the results generation if the user hits enter too quickly
        clearResultsViews();

        processSearchResults(data);
        setupPagination(data.pagination);
    } catch (error) {
        logErrorInConsole(error);
        createToast("Failed to search datasets");
    } finally {
        // show pagination
        for (const classElt of document.getElementsByClassName("pagination")) {
            classElt.classList.remove("is-invisible");
        }
    }

    // Update current user profile with last "controls" settings
    Cookies.set("default_collection_ownership_view", searchCriteria.ownership);
    Cookies.set("default_collection_organism_view", searchCriteria.organism_ids);
    Cookies.set("default_collection_date_added_view", searchCriteria.date_added);
    Cookies.set("default_collection_dataset_type_view", searchCriteria.dtypes);

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

/**
 * Updates the dataset collection buttons based on the ownership status of the collection.
 * If the user is not the owner of the collection, the delete and rename buttons are removed.
 *
 * @param {Object} collection - The dataset collection object.
 */
const updateDatasetCollectionButtons = (collection=null) => {

    // If user is not the owner of the collection, remove the delete and rename buttons
    const isOwner = (collection && Boolean(collection.is_owner)) || false;
    const isPublic = (collection && Boolean(collection.is_public)) || false;

    const collectionRenameButton = document.getElementById("btn-rename-collection");
    const collectionRenamePermalinkButton = document.getElementById("btn-rename-collection-permalink");
    const collectionDeleteButton = document.getElementById("btn-delete-collection");
    const collectionVisibilityContainer = document.getElementById("collection-visibility-c");

    enableAndShowElement(collectionRenameButton);
    enableAndShowElement(collectionRenamePermalinkButton);
    enableAndShowElement(collectionDeleteButton);
    enableAndShowElement(collectionVisibilityContainer);
    if (!isOwner) {
        disableAndHideElement(collectionRenameButton);
        disableAndHideElement(collectionRenamePermalinkButton);
        disableAndHideElement(collectionDeleteButton);
        disableAndHideElement(collectionVisibilityContainer);
        return;
    }

    // Set the visibility of the collection
    const collectionVisibilityInput = document.getElementById("collection-visibility");
    collectionVisibilityInput.closest(".field").querySelector("label").textContent = "Private collection";
    if (isPublic) {
        // If the checkbox is checked, change label accordingly
        collectionVisibilityInput.closest(".field").querySelector("label").textContent = "Public collection";
    }
    collectionVisibilityInput.checked = isPublic;

    // Add event to update the collection visibility on the server
    collectionVisibilityInput.addEventListener("change", async (event) => {
        const visibility = event.target.checked;
        try {
            await apiCallsMixin.updateDatasetCollectionVisibility(datasetCollectionState.selectedShareId, visibility);
            createToast("Collection visibility updated", "is-success");
        } catch (error) {
            logErrorInConsole(error);
            createToast("Failed to update collection visibility");
        }
        // update label
        event.target.closest(".field").querySelector("label").textContent = visibility ? "Public collection" : "Private collection";
    });

}

/**
 * Updates the dataset collections by fetching the dataset collections and updating the dataset collection selector.
 *
 * @returns {Promise<void>} A promise that resolves when the dataset collections are updated.
 */
const updateDatasetCollections = async () => {

    // Fetch the dataset collections, which will update the dataset collection selector
    await fetchDatasetCollections()

    // Uses dataset-collection-selector.js variable
    const datasetCollectionData = datasetCollectionState.data;
    // merge all dataset collection data from domain_layouts, group_layouts, public_layouts, shared_layouts, and user_layouts into one array
    flatDatasetCollectionData = [...datasetCollectionData.domain_layouts, ...datasetCollectionData.group_layouts, ...datasetCollectionData.public_layouts, ...datasetCollectionData.shared_layouts, ...datasetCollectionData.user_layouts];

}

/**
 * Updates the display of "Add to Collection" and "Remove from Collection" buttons in a modal.
 * @param {string} modalDivId - The ID of the modal's div element.
 * @param {object|null} collection - The collection object. If null, the buttons will be disabled.
 */
const updateDisplayAddRemoveToCollectionButtons = (modalDivId, collection=null) => {
    const modalDisplayBox = document.getElementById(modalDivId);

    const modalDisplayElts = modalDisplayBox.getElementsByClassName("js-modal-display");

    for (const modalDisplay of modalDisplayElts) {
        const displayId = parseInt(modalDisplay.dataset.displayId);

        const collectionAddRemoveGroup = modalDisplay.querySelector(".js-collection-add-remove-group");
        const addToCollectionButton = modalDisplay.querySelector("button.js-collection-add-display");
        const removeFromCollectionButton = modalDisplay.querySelector("button.js-collection-remove-display");

        collectionAddRemoveGroup.classList.remove("is-hidden");
        enableAndShowElement(addToCollectionButton);
        enableAndShowElement(removeFromCollectionButton);

        // if dataset collection, remove the "add to collection" and "remove from collection" buttons
        if (!collection) {
            collectionAddRemoveGroup.classList.add("is-hidden");
            disableAndHideElement(addToCollectionButton);
            disableAndHideElement(removeFromCollectionButton);
            continue;
        }

        // if display is not in the currently selected collection, hide the "remove from collection" button
        const displayInCollection = collection.members.some(member => JSON.parse(member).display_id === displayId);
        if (!displayInCollection) {
            disableAndHideElement(removeFromCollectionButton);
        }
    }
}

/**
 * Updates the display access for view display buttons and table columns based on the provided data.
 *
 * @param {Object} data - The data object containing user information.
 * @param {boolean} data.is_owner - Indicates if the current user is the owner.
 */
const updateViewDisplayAccess = (data) => {
    const viewDisplayButtons = document.getElementsByClassName("js-view-displays");
    for (const classElt of viewDisplayButtons) {
        disableAndHideElement(classElt);
        if (data?.is_owner) {
            enableAndShowElement(classElt);
        }
    }
    const viewDisplaysTableHeader = document.getElementById("view-displays-header");
    viewDisplaysTableHeader.classList.add("is-hidden");
    if (data?.is_owner) {
        viewDisplaysTableHeader.classList.remove("is-hidden");
    }
    const columnIndex = viewDisplaysTableHeader.cellIndex;
    // get all rows in the table
    const rows = document.querySelectorAll("#results-table tbody tr.js-table-row");
    for (const row of rows) {
        const cell = row.cells[columnIndex];
        cell.classList.add("is-hidden");
        if (data?.is_owner) {
            cell.classList.remove("is-hidden");
        }
    }
}

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// User settings has no "active" state for the sidebar
	document.getElementById("page-header-label").textContent = "Dataset Explorer";

    const sessionId = getCurrentUser().session_id;
	if (! sessionId ) {
        // ? Technically we can show profiles, but I would need to build in "logged out controls".
        document.getElementById("collection-management").classList.add("is-hidden");
        document.getElementById("filter-collection-container").classList.add("is-hidden");
        // Ensure this toggle remains hidden when user switches views
        document.getElementById("filter-collection-container").classList.remove("js-trigger-dataset-search");
        // only show public datasets option
        for (const elt of document.querySelectorAll("#controls-ownership li:not([data-dbval='public'])")) {
            elt.remove();
        }
        document.querySelector("#controls-ownership li[data-dbval='public']").classList.add("js-selected");

    }

    for (const actionElt of document.querySelectorAll(".js-collection-action-links")) {
        for (const classElt of actionElt.querySelectorAll("[data-tooltip-content]")) {
            const tooltip = createActionTooltips(classElt);
            tooltip.classList.add("js-collection-tooltip"); // prevent removal of tooltip
            applyTooltip(classElt, tooltip, "bottom");
        }
    }

    // Prep filters
    await loadOrganismList();

    // Select the user's last remembered filter options
    const defaultOwnershipView = Cookies.get("default_collection_ownership_view");
    const defaultOrganismView = Cookies.get("default_collection_organism_view");
    const defaultDateAddedView = Cookies.get("default_collection_date_added_view");
    const defaultDatasetTypeView = Cookies.get("default_collection_dataset_type_view");

    if (defaultOwnershipView && getCurrentUser().session_id) {
        // deselect All
        document.querySelector("#controls-ownership li.js-all-selector").classList.remove("js-selected");

        for (const ownership of defaultOwnershipView.split(",")) {
            document.querySelector(`#controls-ownership li[data-dbval='${ownership}']`).classList.add("js-selected");
        }
    }
    if (defaultOrganismView && getCurrentUser().session_id) {
        // deselect All
        document.querySelector("#controls-organism li.js-all-selector").classList.remove("js-selected");
        for (const organism of defaultOrganismView.split(",")) {
            document.querySelector(`#controls-organism li[data-dbval='${organism}']`).classList.add("js-selected");
        }
    }
    if (defaultDateAddedView && getCurrentUser().session_id) {
        document.querySelector(`#controls-date-added li[data-dbval='${getCurrentUser().default_date_added_view}']`).classList.add("js-selected");
    }
    if (defaultDatasetTypeView && getCurrentUser().session_id) {
        // deselect All
        document.querySelector("#controls-dataset-type li.js-all-selector").classList.remove("js-selected");
        for (const dtype of defaultDatasetTypeView.split(",")) {
            document.querySelector(`#controls-dataset-type li[data-dbval='${dtype}']`).classList.add("js-selected");
        }
    }

    // If they passed search_string URL parameter, set that
    let search_string = getUrlParameter("search_string");
    if (search_string) {
        document.getElementById("search-terms").value = search_string;
    }

    await submitSearch();

    // Load dataset collections
    await updateDatasetCollections();
    initializeDatasetCollectionSelection();
    document.getElementById("dropdown-dc").classList.remove("is-right");    // Cannot see the dropdown if it is right aligned

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

document.getElementById("search-clear").addEventListener("click", async () => {
    document.getElementById("search-terms").value = "";
    await submitSearch();
});

// Search for datasets using the supplied search terms
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

document.getElementById("btn-table-view").addEventListener("click", () => {
    listView = "table";
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-table-view").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-table-view").classList.remove('is-dark');

    document.getElementById("title-filter-controls").classList.remove("is-hidden");

    for (const classElt of document.getElementsByClassName("js-trigger-dataset-search")) {
        classElt.classList.remove("is-hidden");
    }
    document.getElementById("results-table").classList.remove("is-hidden");
    document.getElementById("results-list-div").classList.add("is-hidden");
    document.getElementById("dataset-arrangement-c").classList.add("is-hidden");

    document.getElementById("include-public-membership-c").classList.add("is-hidden");

    // Show pagination in case arrangement view hid the pagination
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.remove("is-invisible");
    }
    document.getElementById("count-label-c").classList.remove("is-hidden");

})

document.getElementById("btn-list-view-compact").addEventListener("click", () => {
    listView = "list-compact";
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-list-view-compact").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-list-view-compact").classList.remove('is-dark');

    document.getElementById("title-filter-controls").classList.remove("is-hidden");

    for (const classElt of document.getElementsByClassName("js-trigger-dataset-search")) {
        classElt.classList.remove("is-hidden");
    }
    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.remove("is-hidden");
    document.getElementById("dataset-arrangement-c").classList.add("is-hidden");

    document.getElementById("include-public-membership-c").classList.add("is-hidden");

    // find all elements with class 'js-expandable-view' and make sure they also have 'expanded-view-hidden'
    for (const elt of document.querySelectorAll(".js-expandable-view")){
        elt.classList.add("is-hidden");
    };

    // Toggle js-expand-box to icon collapse
    for (const elt of document.querySelectorAll(".js-expand-box i")){
        elt.classList.remove("mdi-arrow-collapse");
        elt.classList.add("mdi-arrow-expand");
    }

    // Show pagination in case arrangement view hid the pagination
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.remove("is-invisible");
    }
    document.getElementById("count-label-c").classList.remove("is-hidden");

});

document.getElementById("btn-list-view-expanded").addEventListener("click", () => {
    listView = "list-expanded";
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-list-view-expanded").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-list-view-expanded").classList.remove('is-dark');

    document.getElementById("title-filter-controls").classList.remove("is-hidden");

    for (const classElt of document.getElementsByClassName("js-trigger-dataset-search")) {
        classElt.classList.remove("is-hidden");
    }
    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.remove("is-hidden");
    document.getElementById("dataset-arrangement-c").classList.add("is-hidden");

    document.getElementById("include-public-membership-c").classList.remove("is-hidden");

    // find all elements with class 'js-expandable-view' and make sure they also have 'expanded-view-hidden'
    for (const elt of document.querySelectorAll(".js-expandable-view")){
        elt.classList.remove("is-hidden");
    };

    // Toggle js-expand-box to icon expand
    for (const elt of document.querySelectorAll(".js-expand-box i")){
        elt.classList.add("mdi-arrow-collapse");
        elt.classList.remove("mdi-arrow-expand");
    }

    // Show pagination in case arrangement view hid the pagination
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.remove("is-invisible");
    }
    document.getElementById("count-label-c").classList.remove("is-hidden");

});

document.getElementById("btn-arrangement-view").addEventListener("click", () => {
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-arrangement-view").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-arrangement-view").classList.remove('is-dark');

    document.getElementById("include-public-membership-c").classList.add("is-hidden");

    document.getElementById("title-filter-controls").classList.add("is-hidden");

    // Elements that would trigger submitSearch() are hidden so that pagination and count label won't appear
    for (const classElt of document.getElementsByClassName("js-trigger-dataset-search")) {
        classElt.classList.add("is-hidden");
    }

    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.add("is-hidden");

    document.getElementById("dataset-arrangement-c").classList.remove("is-hidden");

    // hide .pagination and #count-label-c
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.add("is-invisible");
    }
    document.getElementById("count-label-c").classList.add("is-hidden");

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

// If checkbox is changed, set the searchByCollection
document.getElementById("filter-only-in-collection").addEventListener("change", (e) => {
    searchByCollection = e.currentTarget.checked;
    // If the checkbox is checked, change label accordingly
    e.currentTarget.closest(".field").querySelector("label").textContent = searchByCollection ? "Yes" : "No"

    // trigger search
    submitSearch();
});

// If checkbox is changed, set the searchByCollection
document.getElementById("include-public-membership").addEventListener("change", (e) => {
    includePublicMembership = e.currentTarget.checked;
    // If the checkbox is checked, change label accordingly
    e.currentTarget.closest(".field").querySelector("label").textContent = includePublicMembership ? "Yes" : "No"

    // trigger search
    submitSearch();
});

document.getElementById("btn-save-arrangement").addEventListener("click", async () => {

    const singlearrangementTiles = document.querySelectorAll("#dataset-arrangement-single .js-sortable-tile");
    const multiarrangementTiles = document.querySelectorAll("#dataset-arrangement-multi .js-sortable-tile");

    // Get the current grid states for all tiles and put in an object

    const layoutArrangement = {"single": [], "multi": []};  // dataset ids as keys

    for (const tile of singlearrangementTiles) {
        const displayId = tile.dataset.displayId;
        layoutArrangement.single.push( {
            "display_id": displayId,
            "start_row": parseInt(tile.style.gridRowStart),
            "start_col": parseInt(tile.style.gridColumnStart),
            "grid_height": parseInt(tile.style.gridRowEnd.split(" ")[1]),
            "grid_width": parseInt(tile.style.gridColumnEnd.split(" ")[1])
        })
    }

    for (const tile of multiarrangementTiles) {
        const displayId = tile.dataset.displayId;
        layoutArrangement.multi.push( {
            "display_id": displayId,
            "start_row": parseInt(tile.style.gridRowStart),
            "start_col": parseInt(tile.style.gridColumnStart),
            "grid_height": parseInt(tile.style.gridRowEnd.split(" ")[1]),
            "grid_width": parseInt(tile.style.gridColumnEnd.split(" ")[1])
        })
    }


    const data = await apiCallsMixin.saveDatasetCollectionArrangement(datasetCollectionState.selectedShareId, layoutArrangement)
    if (data.success) {
        createToast("Layout arrangement saved successfully", "is-success");
    } else {
        createToast("Failed to save layout arrangement");
    }
});

document.getElementById("btn-set-primary-collection").addEventListener("click", async () => {
    try {
        const data = await apiCallsMixin.setUserPrimaryDatasetCollection(datasetCollectionState.selectedShareId)
        if (data.success) {
            createToast("Primary collection set successfully", "is-success");

            Cookies.set('gear_default_domain', datasetCollectionState.selectedShareId);

            // Make button outlined to look "official"
            document.getElementById("btn-set-primary-collection").classList.remove("is-outlined");

        } else {
            throw new Error("Failed to set primary collection");
        }
    } catch (e) {
        createToast(e.message);
    }
});

document.getElementById("btn-share-collection").addEventListener("click", (e) => {
    let currentPage = new URL(`${getRootUrl()}/p`);
    let params = new URLSearchParams(currentPage.search);
    params.set("l", datasetCollectionState.selectedShareId);
    currentPage.search = params.toString();
    const shareUrl = currentPage.toString();
    copyPermalink(shareUrl);
});

