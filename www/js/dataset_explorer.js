"use strict";

/* Imported variables
let dataset_collection_data; // from dataset-collection-selector
let selected_dc_share_id; // from dataset-collection-selector

*/

let firstSearch = true;
let searchByCollection = false;
const resultsPerPage = 20;

let flatDatasetCollectionData = {};   // flattened version of all dataset collections availabe to user

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

class LayoutArrangement {

    constructor(isMulti=false) {
        this.type = isMulti ? "multi" : "single";
        this.arrangement = [];
        this.arrangementDiv = document.getElementById(`dataset-arrangement-${this.type}`);

        this.arrangementWidth = 1080; // NOTE: Originally used singleGeneArrangementDiv.offsetWidth, but it is 0 unless the element is visible
        this.rowWidth = this.arrangementWidth / 12; // Split width into 12 columns
        this.colHeight = this.rowWidth * 4; // A unit of height for us is 4 units of width (to make a square)
    }

    addMember(member) {
        this.arrangement.push(member);
    }

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

    constructor(arrangementObj, datasetId, gridPosition, startCol, startRow, gridWidth, gridHeight) {
        this.parentArrangement = arrangementObj;
        this.datasetId = datasetId;
        this.gridPosition = gridPosition;
        this.startCol = startCol;
        this.startRow = startRow;
        this.gridWidth = gridWidth;
        this.gridHeight = gridHeight;
        this.datasetTitle = "";
        this.image = "";
        this.tileTemplate = document.getElementById("dataset-arrangement-tile-template");
        this.selector = `.js-sortable-tile[data-dataset-id="${this.datasetId}"]`;
        this.snapWidth = this.parentArrangement.rowWidth * 2; // Snap to 1/6 increments for width
        this.snapHeight = this.parentArrangement.colHeight; // Snap to 1/4 increments for height;
    }

    createArrangementTile() {
        const tile = this.tileTemplate.content.cloneNode(true);

        setElementProperties(tile, ".js-sortable-tile", { dataset: {datasetId: this.datasetId} });
        // add gridAra to tile
        tile.querySelector(".js-sortable-tile").style.gridArea = `${this.startRow} / ${this.startCol} / span ${this.gridHeight} / span ${this.gridWidth}`;
        // add dataset title
        setElementProperties(tile, ".js-sortable-tile-title", { textContent: this.datasetTitle });
        // add preview image
        const tileImg = tile.querySelector(".js-sortable-tile-img");
        tileImg.src = this.image;
        return tile;
    }

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
                end: this.determineGridOverlap
            }

        }).resizable({
            // resize from only the right and bottom edges (adding top and left adds complexity)
            edges: { left: false, right: true, bottom: true, top: false },
            listeners: {
                move (event) {

                    // get current style gridArea values (at least the ones to preserve)
                    const rowStart = parseInt(event.target.style.gridRowStart);
                    const colStart = parseInt(event.target.style.gridColumnStart);

                    // Snap to grid in 1/6 increments for width
                    const newSpanWidth = (Math.round(event.rect.width / that.snapWidth) * 2);  // x2 converts to 12 column grid
                    const newSpanHeight = Math.round(event.rect.height / that.snapHeight);

                    event.target.style.gridArea = `${rowStart} / ${colStart} / span ${newSpanHeight} / span ${newSpanWidth}`;

                    // When resized to the n+1 row, add a new row to the grid.
                    const arrangementDiv = event.target.parentElement;    // AKA this.parentArrangement.arrangementDiv
                    // get current number of rows in the grid
                    // it is stored as a string like "repeat(3, 100px)"
                    const arrangementRows = parseInt(arrangementDiv.style.gridTemplateRows.split(",")[0].split("(")[1]);

                    const lastTileRow = rowStart + newSpanHeight -1;
                    // If the dragged tile's bottom edge is below the last row, add another row to the grid template
                    if (lastTileRow > arrangementRows) {
                        arrangementDiv.style.gridTemplateRows = `repeat(${lastTileRow}, ${that.parentArrangement.colHeight}px)`;
                    }
                },
                end: this.determineGridOverlap
            }
        });

    }

    determineGridOverlap(event) {
        // Determine if any tiles overlap with one another on the parent grid and provide visual cue

        const eventRowStart = parseInt(event.target.style.gridRowStart);
        const eventRowEnd = eventRowStart + parseInt(event.target.style.gridRowEnd.split(" ")[1]);
        const eventColStart = parseInt(event.target.style.gridColumnStart);
        const eventColEnd = eventColStart + parseInt(event.target.style.gridColumnEnd.split(" ")[1]);

        // Get all tiles in the arrangement
        const arrangementDiv = event.target.parentElement;
        const arrangementTiles = arrangementDiv.querySelectorAll(".js-sortable-tile");

        // Check for overlap with other tiles
        for (const tile of arrangementTiles) {
            if (tile === event.target) {
                continue;
            }

            // Get the grid area of the current tile and event tile
            const tileRowStart = parseInt(tile.style.gridRowStart);
            const tileRowEnd = tileRowStart + parseInt(tile.style.gridRowEnd.split(" ")[1]);
            const tileColStart = parseInt(tile.style.gridColumnStart);
            const tileColEnd = tileColStart + parseInt(tile.style.gridColumnEnd.split(" ")[1]);

            // Check for overlap
            // This reminds me of the video game detection algorithms
            event.target.classList.remove("js-is-overlapping")
            event.target.style.opacity = 1;
            event.target.querySelector(".js-sortable-tile-title").classList.remove("has-text-danger-light", "has-background-danger-dark");
            event.target.querySelector(".js-sortable-tile-title").classList.add("has-background-primary-light");
            document.getElementById("btn-save-arrangement").disabled = false;

            if (tileRowStart < eventRowEnd
                && tileRowEnd > eventRowStart
                && tileColStart < eventColEnd
                && tileColEnd > eventColStart) {
                // Overlap detected
                event.target.classList.add("js-is-overlapping");
                event.target.style.opacity = 0.5;
                event.target.querySelector(".js-sortable-tile-title").classList.add("has-text-danger-light", "has-background-danger-dark");
                event.target.querySelector(".js-sortable-tile-title").classList.remove("has-background-primary-light");
                break;
            }

        }

        // only enable when no tiles have "js-is-overlapping" class
        if( [...document.querySelectorAll(".js-sortable-tile")].some(tile => tile.classList.contains("js-is-overlapping")) ) {
            document.getElementById("btn-save-arrangement").disabled = true;
        } else {
            // All tiles are valid, clear the "overlap" appearances
            for (const tile of document.querySelectorAll(".js-sortable-tile")) {
                tile.classList.remove("js-is-overlapping");
                tile.style.opacity = 1;
                tile.querySelector(".js-sortable-tile-title").classList.remove("has-text-danger-light", "has-background-danger-dark");
                tile.querySelector(".js-sortable-tile-title").classList.add("has-background-primary-light");
            }
        }

    }

}

/**
 * Adds event listeners for various actions related to datasets.
 * @function
 * @returns {void}
 */
const addDatasetListEventListeners = () => {

    // Add event listener to analysis dropdown trigger
    for (const classElt of document.querySelectorAll(".js-analysis-dropdown .dropdown-trigger")) {
        classElt.addEventListener("click", (event) => {
            const item = event.currentTarget;
            item.closest(".dropdown").classList.toggle('is-active');
        });
    };

    // Expand and collapse dataset view
    for (const classElt of document.getElementsByClassName("js-expand-box")) {
        classElt.addEventListener("click", (e) => {
            const datasetId = e.currentTarget.dataset.datasetId;
            const selector = `#result-dataset-id-${datasetId} .js-expandable-view`;
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

    for (const classElt of document.getElementsByClassName("js-download-dataset")) {
        classElt.addEventListener("click", async (e) => {
            // download the h5ad
            const datasetId = e.currentTarget.dataset.datasetId;
            const url = `./cgi/download_source_file.cgi?type=h5ad&dataset_id=${datasetId}`;
            const {data} = await axios.get(url, {responseType: 'blob'});
            const blob = new Blob([data], {type: 'application/octet-stream'});
            const downloadUrl = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = downloadUrl;
            a.download = `${datasetId}.h5ad`;
            a.click();
        });
    }

    for (const classElt of document.getElementsByClassName("js-share-dataset")) {
        classElt.addEventListener("click", (e) => {

            const shareId = e.currentTarget.value;
            const currentUrl = window.location.href;
            const currentPage = currentUrl.lastIndexOf("dataset_explorer.v2.html");
            const shareUrl = `${currentUrl.substring(0, currentPage)}p?s=${shareId}`;
            const datasetId = e.currentTarget.dataset.datasetId;

            showDatasetActionNote(datasetId, shareUrl);

        });
    }

    // TODO: Feature - add dataset to collection multiple times, with indicator to show how many times it is in the collection
    for (const classElt of document.getElementsByClassName("js-collection-add-dataset")) {
        classElt.addEventListener("click", async (e) => {
            let collection = flatDatasetCollectionData.find((collection) => collection.share_id === selected_dc_share_id);

            // Currently only collectiosn with datasets members will be fetched.
            // If this isn't one (i.e. brand new one), we can fudge some properties to make the UI work
            if (!collection) {
                collection = {
                    dataset_count: 0,
                    folder_id: null,
                    folder_label: null,
                    folder_parent_id: null,
                    is_current: 0,
                    is_domain: 0,
                    is_owner: 1,
                    is_public: 0,
                    label: selected_dc_label,
                    members: [],
                    share_id: selected_dc_share_id,
                }
            }

            const datasetId = e.currentTarget.dataset.datasetId;
            const datasetIsPublic = e.currentTarget.dataset.isPublic;

            // if dataset is private and the collection is public, then we need to show a popover to confirm to the user the dataset will switch to public.
            if (!parseBool(datasetIsPublic) && Boolean(collection.is_public)) {
                createSwitchDatasetToPublicPopover(classElt, datasetId);
            } else {
                try {
                    const data = await apiCallsMixin.addDatasetToCollection(selected_dc_share_id, datasetId);
                    if (!data.success) {
                        throw new Error(data.error);
                    }

                    changeDatasetCollectionCallback();  // Updates the buttons and UI for the selected collection
                    createToast("Dataset added to collection", "is-success");


                } catch (error) {
                    logErrorInConsole(error);
                    createToast("Failed to add dataset to collection");
                }
            }
        });
    }

    // TODO: remove last instance of dataset from collection
    for (const classElt of document.getElementsByClassName("js-collection-remove-dataset")) {
        classElt.addEventListener("click", async (e) => {
            const datasetId = e.currentTarget.dataset.datasetId;

            try {
                const data = await apiCallsMixin.deleteDatasetFromCollection(selected_dc_share_id, datasetId);
                if (!data.success) {
                    throw new Error(data.error);
                }

                changeDatasetCollectionCallback();  // Updates the buttons and UI for the selected collection

                createToast("Dataset removed from collection", "is-success");


            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to remove dataset from collection");
            }
        });
    }

    // Cancel button for editing a dataset
    for (const classElt of document.getElementsByClassName("js-edit-dataset-cancel")) {
        classElt.addEventListener("click", (e) => {
            const datasetId = e.currentTarget.dataset.datasetId;
            const selectorBase = `#result-dataset-id-${datasetId}`;

            // Show editable versions where there are some and hide the display versions
            for (const classElt of document.querySelectorAll(`${selectorBase} .js-editable-version`)) {
                classElt.classList.add("is-hidden");
            };
            for (const classElt of document.querySelectorAll(`${selectorBase} .js-display-version`)) {
                classElt.classList.remove("is-hidden");
            };

            // Reset any unsaved/edited values
            const visibility = document.querySelector(`${selectorBase}-editable-visibility`).dataset.originalVal;
            document.querySelector(`${selectorBase}-editable-visibility`).value = visibility;

            const title = document.querySelector(`${selectorBase}-editable-title`).dataset.originalVal;
            document.querySelector(`${selectorBase}-editable-title`).value = title;

            document.querySelector(`${selectorBase} .js-action-links`).classList.remove("is-hidden");

        });
    }

    // Save button for editing a dataset
    for (const classElt of document.getElementsByClassName("js-edit-dataset-save")) {
        classElt.addEventListener("click", async (e) => {
            const datasetId = e.currentTarget.dataset.datasetId;
            const selectorBase = `#result-dataset-id-${datasetId}`;
            //
            const newVisibility = document.querySelector(`${selectorBase}-editable-visibility`).checked;
            // convert "true/false" visibility to 1/0
            const intNewVisibility = newVisibility ? 1 : 0;

            const isDownloadable = document.querySelector(`${selectorBase}-editable-downloadable`).checked;
            // convert "true/false" visibility to 1/0
            const intIsDownloadable = isDownloadable ? 1 : 0;

            const newTitle = document.querySelector(`${selectorBase}-editable-title`).value;
            const newPubmedId = document.querySelector(`${selectorBase}-editable-pubmed-id`).value;
            const newGeoId = document.querySelector(`${selectorBase}-editable-geo-id`).value;
            const newLdesc = document.querySelector(`${selectorBase}-editable-ldesc`).value;

            try {
                const data = await apiCallsMixin.saveDatasetInfoChanges(datasetId, intNewVisibility, intIsDownloadable, newTitle, newPubmedId, newGeoId, newLdesc);
                createToast("Dataset changes saved", "is-success");

            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to save dataset changes");
                return;
            } finally {
                document.querySelector(`${selectorBase} .js-action-links`).classList.remove("is-hidden");
            }

            // Update the UI for the new values
            document.querySelector(`${selectorBase}-editable-visibility`).dataset.isPublic = newVisibility;
            if (newVisibility) {
                document.querySelector(`${selectorBase}-display-visibility`).textContent = "Public dataset";
                document.querySelector(`${selectorBase}-table-visibility`).textContent = "Public";
                document.querySelector(`${selectorBase}-display-visibility`).classList.remove("is-danger");
                document.querySelector(`${selectorBase}-display-visibility`).classList.add("is-light", "is-primary");
                document.querySelector(`${selectorBase}-table-visibility`).classList.remove("has-background-danger");
                document.querySelector(`${selectorBase}-table-visibility`).classList.add("has-background-primary-light");

            } else {
                document.querySelector(`${selectorBase}-display-visibility`).textContent = "Private dataset";
                document.querySelector(`${selectorBase}-table-visibility`).textContent = "Private";
                document.querySelector(`${selectorBase}-display-visibility`).classList.remove("is-light", "is-primary");
                document.querySelector(`${selectorBase}-display-visibility`).classList.add("is-danger");
                document.querySelector(`${selectorBase}-table-visibility`).classList.remove("has-background-primary-light");
                document.querySelector(`${selectorBase}-table-visibility`).classList.add("has-background-danger");
            }

            document.querySelector(`${selectorBase}-editable-downloadable`).dataset.isDownloadable = isDownloadable;

            if (isDownloadable) {
                document.querySelector(`${selectorBase}-display-downloadable`).textContent = "Downloadable";
                document.querySelector(`${selectorBase}-display-downloadable`).classList.remove("is-dark");
                document.querySelector(`${selectorBase}-display-downloadable`).classList.add("is-success");
            } else {
                document.querySelector(`${selectorBase}-display-downloadable`).textContent = "Not downloadable";
                document.querySelector(`${selectorBase}-display-downloadable`).classList.remove("is-success");
                document.querySelector(`${selectorBase}-display-downloadable`).classList.add("is-dark");
            }
            const downloadButton = document.querySelector(`${selectorBase}-download-dataset`);
            if (downloadButton) {
                // If button exists (has h5ad), update the button visibility
                downloadButton.dataset.isDownloadable = isDownloadable;
                downloadButton.disabled = !isDownloadable;
                downloadButton.parentElement.classList.toggle("is-hidden", !isDownloadable);
            }

            document.querySelector(`${selectorBase}-editable-title`).dataset.originalVal = newTitle;
            document.querySelector(`${selectorBase}-display-title`).textContent = newTitle;
            document.querySelector(`${selectorBase}-table-title`).textContent = newTitle;

            document.querySelector(`${selectorBase}-editable-ldesc`).dataset.originalVal = newLdesc;
            document.querySelector(`${selectorBase}-display-ldesc`).textContent = newLdesc || "No description entered";

            // pubmed and geo display are links if they exist
            document.querySelector(`${selectorBase}-editable-pubmed-id`).value = newPubmedId;
            if (newPubmedId) {
                document.querySelector(`${selectorBase}-display-pubmed-id`).innerHTML = `<a href="https://pubmed.ncbi.nlm.nih.gov/${newPubmedId}" target="_blank">
                    <span>${newPubmedId}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>`;
            } else {
                document.querySelector(`${selectorBase}-display-pubmed-id`).textContent = "Not available";
            }

            document.querySelector(`${selectorBase}-editable-geo-id`).value = newGeoId;
            if (newGeoId) {
                document.querySelector(`${selectorBase}-display-geo-id`).innerHTML = `<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${newGeoId}" target="_blank">
                    <span>${newGeoId}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>`;
            } else {
                document.querySelector(`${selectorBase}-display-geo-id`).textContent = "Not available";
            }
            document.querySelector(`${selectorBase}-display-geo-id`).value = newGeoId;

            // Put interface back to view mode.
            toggleEditableMode(true, selectorBase);

        });
    }

    // Toggle editable mode when edit button is clicked for a dataset
    for (const classElt of document.getElementsByClassName("js-edit-dataset")) {
        classElt.addEventListener("click", async (e) => {

            const datasetId = e.currentTarget.dataset.datasetId;
            const selectorBase = `#result-dataset-id-${datasetId}`;

            const editableVisibilityElt = document.querySelector(`${selectorBase}-editable-visibility`);

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
    for (const classElt of document.getElementsByClassName("js-view-dataset")) {
        classElt.addEventListener("click", (e) => {
            // ! Currently redirects to blank page
            window.open(`./p?s=${e.currentTarget.value}`, '_blank');
        });
    }
}

/**
 * Adds downloadable information to a dataset.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @param {boolean} isDownloadable - Indicates whether the dataset is downloadable.
 * @param {boolean} hasH5ad - Indicates whether the dataset has an H5ad file.
 */
const addDownloadableInfoToDataset = (datasetId, isDownloadable, hasH5ad) => {
    const datasetDisplayContainer = document.getElementById(`${datasetId}-display-downloadable`);
    const datasetDisplaySpan = document.createElement("span");
    datasetDisplaySpan.classList.add("tag");
    datasetDisplaySpan.id = `result-dataset-id-${datasetId}-display-downloadable`;

    if (hasH5ad && isDownloadable) {
        datasetDisplaySpan.classList.add("is-success");
        datasetDisplaySpan.textContent = "Downloadable";
    } else {
        datasetDisplaySpan.classList.add("is-dark");
        datasetDisplaySpan.textContent = "Not downloadable";
    }
    datasetDisplayContainer.appendChild(datasetDisplaySpan);
}

/**
 * Adds public/private visibility information to a dataset display container in the DOM.
 * @param {string} datasetId - The ID of the datasetdisplay container.
 * @param {boolean} isPublic - A boolean indicating whether the dataset is public or private.
 * @returns {void}
 */
const addVisibilityInfoToDataset = (datasetId, isPublic) => {

    // add dataset public/private info to DOM
    const datasetDisplayContainer = document.getElementById(`${datasetId}-display-visibility`);
    const datasetDisplaySpan = document.createElement("span");
    datasetDisplaySpan.classList.add("tag");
    datasetDisplaySpan.id = `result-dataset-id-${datasetId}-display-visibility`;
    const datasetTableVisibility = document.getElementById(`result-dataset-id-${datasetId}-table-visibility`);

    if (isPublic) {
        datasetDisplaySpan.classList.add("is-primary", "is-light");
        datasetDisplaySpan.textContent = "Public dataset";
        datasetTableVisibility.classList.add("has-background-primary-light");
        datasetTableVisibility.textContent = "Public";
    } else {
        datasetDisplaySpan.classList.add("is-danger");
        datasetDisplaySpan.textContent = "Private dataset";
        datasetTableVisibility.classList.add("has-background-danger");
        datasetTableVisibility.textContent = "Private";
    }
    datasetDisplayContainer.appendChild(datasetDisplaySpan);

    // Toggle switch (public is checked, private is unchecked)
    const visibilitySwitch = document.getElementById(`result-dataset-id-${datasetId}-editable-visibility`);

    visibilitySwitch.addEventListener("change", (e) => {
        const isPublic = e.currentTarget.checked;
        e.currentTarget.dataset.isPublic = isPublic;
        e.currentTarget.closest(".field").querySelector("label").textContent = isPublic ? "Public" : "Private";
    });

    const downloadableSwitch = document.getElementById(`result-dataset-id-${datasetId}-editable-downloadable`);

    downloadableSwitch.addEventListener("change", (e) => {
        const isDownloadable = e.currentTarget.checked;
        e.currentTarget.dataset.isDownloadable = isDownloadable;
        e.currentTarget.closest(".field").querySelector("label").textContent = isDownloadable ? "Yes" : "No";
    });

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
 * Callback function that is triggered when a new dataset collection is selected.
 * It fetches the dataset collections and collection members from the API,
 * updates the UI based on the selected collection, and handles button actions.
 * @returns {Promise<void>} A promise that resolves when the function completes.
 */
const changeDatasetCollectionCallback = async (datasetCollectionData=null) => {
    // Show action buttons
    document.getElementById("collection-actions-c").classList.remove("is-hidden");

    // Reset the arrangement views
    const arrangementViewSingle = document.getElementById("dataset-arrangement-single");
    const arrangementViewMulti = document.getElementById("dataset-arrangement-multi");
    arrangementViewSingle.innerHTML = "";
    arrangementViewMulti.innerHTML = "";

    // The share_id should be updated in the component when a new dataset collection is selected
    const datasetData = await apiCallsMixin.fetchDatasets({layout_share_id: selected_dc_share_id, sort_by: "date_added"})

    // If dataset collection data is not passed in, fetch it from the API
    if (!datasetCollectionData) {
        datasetCollectionData = await apiCallsMixin.fetchDatasetCollections();
    }

    // merge all dataset collection data from domain_layouts, group_layouts, public_layouts, shared_layouts, and user_layouts into one array
    flatDatasetCollectionData = [...datasetCollectionData.domain_layouts, ...datasetCollectionData.group_layouts, ...datasetCollectionData.public_layouts, ...datasetCollectionData.shared_layouts, ...datasetCollectionData.user_layouts]

    // Find the dataset collection data for the selected share_id
    let collection = flatDatasetCollectionData.find((collection) => collection.share_id === selected_dc_share_id);

    // Currently only collectiosn with datasets members will be fetched.
    // If this isn't one (i.e. brand new one), we can fudge some properties to make the UI work
    if (!collection) {
        collection = {
            dataset_count: 0,
            folder_id: null,
            folder_label: null,
            folder_parent_id: null,
            is_current: 0,
            is_domain: 0,
            is_owner: 1,
            is_public: 0,
            label: selected_dc_label,
            members: [],
            share_id: selected_dc_share_id,
        }
    }

    const imageUrls = {"single":{}, "multi":{}};
    const titles = {};
    for (const dataset of datasetData.datasets) {
        const previewImageUrl = dataset.preview_image_url || "/img/dataset_previews/missing.png";
        const mgPreviewImageUrl = dataset.mg_preview_image_url || "/img/dataset_previews/missing.png";
        imageUrls["single"][dataset.dataset_id] = previewImageUrl;
        imageUrls["multi"][dataset.dataset_id] = mgPreviewImageUrl;
        titles[dataset.dataset_id] = dataset.title;
    }

    // Create popvers for collection actions
    createDeleteCollectionConfirmationPopover();
    createNewCollectionPopover();
    createRenameCollectionPopover();

    // Update action buttons for the dataset collection or datasets
    updateDatasetCollectionButtons(collection);
    updateDatasetAddRemoveToCollectionButtons(collection);

    const isDomain = Boolean(collection.is_domain);

    // If selected dataset collection is a "domain" collection, hide the "arrangement" view button.
    // if arrangement view is active (class "gear-bg-secondary") switch to table view
    document.getElementById("btn-arrangement-view").classList.remove("is-hidden");
    if (isDomain) {
        document.getElementById("btn-table-view").click();
        document.getElementById("btn-arrangement-view").classList.add("is-hidden");
    }

    // If the selected dataset collection is the current collection, make it look like the primary collection
    document.getElementById("btn-set-primary-collection").classList.add("is-outlined");
    if (collection.share_id === Cookies.get('gear_default_domain')) {
        document.getElementById("btn-set-primary-collection").classList.remove("is-outlined");
    }

    // JSON parse every layout member
    const layoutMembers = collection.members

    // Legacy mode - if all tiles have startCol = 1, then we are in legacy mode
    // These layouts were generated only with a "width" property
    const legacyMode = layoutMembers.every((dataset) => dataset.start_col === 1);

    if (legacyMode) {
        console.info("Legacy grid mode detected");
    }

    let currentCol = 1;
    let currentRow = 1;
    let currentMgCol = 1;
    let currentMgRow = 1;
    const maxEndCol = 13;   // 12 grid slots. CSS grid is 1-indexed, so 13 is the max

    singleArrangement = new LayoutArrangement();
    multiArrangement = new LayoutArrangement(true);

    for (const dataset of layoutMembers) {

        const datasetId = dataset.dataset_id;

        const singleMember = new LayoutArrangementMember(singleArrangement, datasetId, dataset.grid_position, dataset.start_col, dataset.start_row, dataset.grid_width, dataset.grid_height);
        const multiMember = new LayoutArrangementMember(multiArrangement, datasetId, dataset.mg_grid_position, dataset.mg_start_col, dataset.mg_start_row, dataset.mg_grid_width, dataset.mg_grid_height);

        // If in legacy mode, then we need to calculate the startCol and endCol and startRow and endRow
        // so the arrangement view can be displayed correctly
        if (legacyMode) {
            const width = dataset.grid_width;

            // If endCol is greater than 13, then this tile is in the next row
            if (currentCol + width > maxEndCol) {
                currentCol = 1;
                currentRow++;
            }

            singleMember.startCol = currentCol;
            singleMember.startRow = currentRow;

            currentCol += width;

            // Do this for multi gene view as well
            const mgWidth = dataset.mg_grid_width;

            if (currentMgCol + mgWidth > maxEndCol) {
                currentMgCol = 1;
                currentMgRow++;
            }

            multiMember.startCol = currentMgCol;
            multiMember.startRow = currentMgRow;

            currentMgCol += mgWidth;

        }

        singleMember.datasetTitle = titles[datasetId];
        multiMember.datasetTitle = titles[datasetId];

        singleMember.image = imageUrls.single[datasetId];
        multiMember.image = imageUrls.multi[datasetId];

        singleArrangement.addMember(singleMember);
        multiArrangement.addMember(multiMember);
    }

    singleArrangement.setupArrangementAdjustable();
    multiArrangement.setupArrangementAdjustable();


}

/**
 * Creates a confirmation popover for deleting a dataset.
 */
const createDeleteDatasetConfirmationPopover = () => {
    const deleteButtons = document.getElementsByClassName("js-delete-dataset");
    for (const button of deleteButtons) {
        button.addEventListener('click', (e) => {
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

            // Store the dataset ID to delete
            const datasetIdToDelete = e.currentTarget.value;

            // Add event listener to cancel button
            document.getElementById('cancel-dataset-delete').addEventListener('click', () => {
                popoverContent.remove();
            });

            // Add event listener to confirm button
            document.getElementById('confirm-dataset-delete').addEventListener('click', async () => {

                try {
                    const data = await apiCallsMixin.deleteDataset(datasetIdToDelete);

                    if (data['success'] == 1) {
                        const resultElement = document.getElementById(`result-dataset-id-${datasetIdToDelete}`);
                        resultElement.style.transition = 'opacity 1s';
                        resultElement.style.opacity = 0;
                        resultElement.remove();

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
                    popoverContent.remove();
                }
            });
        });
    }
}

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

        // Add event listener to cancel button
        document.getElementById('cancel-collection-delete').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-collection-delete').addEventListener('click', async () => {

            try {
                const data = await apiCallsMixin.deleteDatasetCollection(selected_dc_share_id);

                if (data['success'] == 1) {
                    // Re-fetch the dataset collections, which will update in the UI via click events
                    await fetchDatasetCollections()

                    createToast("Dataset collection deleted", "is-success");

                    selectDatasetCollection(null);  // performs DatasetCollectionSelectorCallback when label is reset

                    // Override the "show" from the callback. If the collection is deleted, the collection management should be hidden
                    document.getElementById("collection-actions-c").classList.add("is-hidden");

                    // Remove the collection from the dropdown so that it cannot be selected again
                    document.querySelector(`.dropdown-dc-item[data-share-id="${selected_dc_share_id}"]`).remove();

                } else {
                    const error = data['error'] || "Failed to delete collection";
                    throw new Error(error);
                }
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to delete collection");
            } finally {
                popoverContent.remove();

            }
        });
    });
}

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

        document.getElementById("collection-name").addEventListener("keyup", () => {
            const newCollectionName = document.getElementById("collection-name");
            const confirmAddCollection = document.getElementById("confirm-collection-add");

            if (newCollectionName.value.length == 0) {
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
        document.getElementById('confirm-collection-add').addEventListener('click', async () => {
            const newName = document.getElementById("collection-name").value;

            try {
                const data = await apiCallsMixin.createDatasetCollection(newName);

                if (data['layout_share_id']) {
                    // Re-fetch the dataset collections, which will update the UI via click events
                    await fetchDatasetCollections()

                    createToast("Dataset collection created", "is-success");

                    selectDatasetCollection(data['layout_share_id']); // performs DatasetCollectionSelectorCallback when label is set

                } else {
                    const error = data['error'] || "Failed to create new collection";
                    throw new Error(error);
                }
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to create new collection");
            } finally {
                popoverContent.remove();
            }
        });
    });
}

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
                        <button id='confirm-collection-rename' class='button is-primary'>Add</button>
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

        // Show popover
        document.body.appendChild(popoverContent);

        document.getElementById("collection-name").addEventListener("keyup", () => {
            const newCollectionName = document.getElementById("collection-name");
            const confirmAddCollection = document.getElementById("confirm-collection-rename");

            if (newCollectionName.value.length == 0) {
                confirmAddCollection.disabled = true;
                return;
            }
            confirmAddCollection.disabled = false;
        });

        // Add event listener to cancel button
        document.getElementById('cancel-collection-rename').addEventListener('click', () => {
            popoverContent.remove();
        });

        // Add event listener to confirm button
        document.getElementById('confirm-collection-rename').addEventListener('click', async () => {
            const newName = document.getElementById("collection-name").value;

            try {
                const data = await apiCallsMixin.renameCollection(selected_dc_share_id, newName);

                if (data['layout_label']) {
                    // Re-fetch the dataset collections, which will update the UI via click events
                    await fetchDatasetCollections()

                    createToast("Dataset collection renamed", "is-success");

                    selectDatasetCollection(data['layout_share_id']); // performs DatasetCollectionSelectorCallback when label is set

                } else {
                    const error = data['error'] || "Failed to rename collection";
                    throw new Error(error);
                }
            } catch (error) {
                logErrorInConsole(error);
                createToast("Failed to rename collection");
            } finally {
                popoverContent.remove();
            }
        });
    });
}

const createSwitchDatasetToPublicPopover = (addButton, datasetId) => {
    // remove existing popovers
    const existingPopover = document.getElementById('toggle-dataset-to-public-popover');
    if (existingPopover) {
        existingPopover.remove();
    }

    // Create popover content
    const popoverContent = document.createElement('article');
    popoverContent.id = 'toggle-dataset-to-public-popover';
    popoverContent.classList.add("message", "is-warning");
    popoverContent.setAttribute("role", "tooltip");
    popoverContent.style.width = "500px";
    popoverContent.innerHTML = `
        <div class='message-header'>
            <p>Dataset privacy change</p>
        </div>
        <div class='message-body'>
            <p>The selected dataset collection is public, but this dataset is private.
            Adding this dataset to the collection will make the dataset public. Do you want to proceed?</p>
            <div class='field is-grouped' style='width:250px'>
                <p class="control">
                    <button id='confirm-dataset-public-toggle' class='button'>Yes</button>
                </p>
                <p class="control">
                    <button id='cancel-dataset-public-toggle' class='button' value='cancel_add'>No</button>
                </p>
            </div>
        </div>
        <div id="arrow"></div>
    `;

    // append element to DOM to get its dimensions
    document.body.appendChild(popoverContent);

    const arrowElement = document.getElementById('arrow');

    // Create popover (help from https://floating-ui.com/docs/tutorial)
    computePosition(addButton, popoverContent, {
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

    // Add event listener to cancel button
    document.getElementById('cancel-dataset-public-toggle').addEventListener('click', () => {
        popoverContent.remove();
    });

    // Add event listener to confirm button
    document.getElementById('confirm-dataset-public-toggle').addEventListener('click', async () => {
        try {
            // Add to collection (+ make dataset public = true)
            const data = await apiCallsMixin.addDatasetToCollection(selected_dc_share_id, datasetId, true);
            if (!data.success) {
                const error = data['error'] || "Failed to add dataset to collection";
                throw new Error(error);
            }
            await submitSearch();

            changeDatasetCollectionCallback();  // Update the buttons and UI for the selected collection

        } catch (error) {
            logErrorInConsole(error);
            createToast("Failed to add dataset to collection");
        } finally {
            popoverContent.remove();
        }
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

const datasetCollectionSelectCallback = (datasetCollectionData) => {
    // Add mutation observer to watch if #dropdown-dc-selector-label changes
    const observer = new MutationObserver((mutationsList, observer) => {
        for (const mutation of mutationsList) {
            if (mutation.type === 'childList') {
                changeDatasetCollectionCallback(datasetCollectionData);
            }
        }
    });

    observer.observe(document.getElementById("dropdown-dc-selector-label"), { childList: true });

    // Trigger the default dataset collection to be selected in the
    if (CURRENT_USER.default_profile_share_id) {
        selectDatasetCollection(CURRENT_USER.default_profile_share_id);
    }
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

    const tableResultsBody = document.querySelector("#results-table tbody");
    const tableTamplate = document.getElementById("results-table-view")

    const resultsListDiv = document.getElementById("results-list-div");
    const listTemplate = document.getElementById("results-list-view");

    // data.datasets is a list of JSON strings
    for (const dataset of data.datasets) {
        const datasetId = dataset.id;
        const datasetType = dataset.dtype;
        const label = dataset.title;

        const longDesc = dataset.ldesc || "";
        const shareId = dataset.share_id;
        const isPublic = Boolean(dataset.is_public);
        const isDownloadable = Boolean(dataset.is_downloadable);
        const hasH5ad = Boolean(dataset.has_h5ad);
        const dateAdded = new Date(dataset.date_added).toDateString();
        // as YYYY/MM/DD
        const shortDateAdded = new Date(dataset.date_added).toISOString().slice(0, 10);

        const userName = dataset.user_name;
        const organism = dataset.organism;
        const isOwner = dataset.is_owner

        const annotSource = dataset.annotation_source || "Not given";
        const annotVersion = dataset.annotation_release || "Not given";

        const pubmedId = dataset.pubmed_id || null;
        const geoId = dataset.geo_id || null;

        const previewImageUrl = dataset.preview_image_url || "/img/dataset_previews/missing.png";

        const resultDatasetId = `result-dataset-id-${datasetId}`

        // TABLE VIEW

        // Clone the template
        const tableResultsView = tableTamplate.content.cloneNode(true)

        // Set properties for multiple elements
        // TODO: Truncate titles with tooltip for full title
        setElementProperties(tableResultsView, ".js-display-title", { id: `${resultDatasetId}-table-title`, textContent: label });
        setElementProperties(tableResultsView, ".js-display-visibility", { id: `${resultDatasetId}-table-visibility` });
        setElementProperties(tableResultsView, ".js-display-organism", { id: `${resultDatasetId}-table-organism`, textContent: organism });
        setElementProperties(tableResultsView, ".js-display-owner", { textContent: userName });
        setElementProperties(tableResultsView, ".js-display-date-added", { textContent: shortDateAdded });
        setElementProperties(tableResultsView, ".js-display-dataset-type", { textContent: datasetType });

        // Append the cloned template to the results container
        tableResultsBody.appendChild(tableResultsView);

        // LIST VIEW

        // Clone the template
        const listResultsView = listTemplate.content.cloneNode(true)

        // Set properties for multiple elements
        setElementProperties(listResultsView, ".js-dataset-list-element", { id: resultDatasetId, dataset: { datasetId } });

        // Figure section
        setElementProperties(listResultsView, ".js-dataset-list-element figure", { id: `${resultDatasetId}-figure` });

        // title section
        setElementProperties(listResultsView, ".js-display-title p", { id: `${resultDatasetId}-display-title`, textContent: label });
        setElementProperties(listResultsView, ".js-editable-title input", { id: `${resultDatasetId}-editable-title`, dataset: { originalVal: label }, value: label });
        setElementProperties(listResultsView, ".js-expand-box", { dataset: { datasetId } });
        // visibility/other metadata section
        setElementProperties(listResultsView, ".js-display-visibility", { id: `${datasetId}-display-visibility` });
        setElementProperties(listResultsView, ".js-editable-visibility input", { id: `${resultDatasetId}-editable-visibility`, checked: isPublic, dataset: { isPublic } });
        setElementProperties(listResultsView, ".js-editable-visibility label", { htmlFor: `${resultDatasetId}-editable-visibility`, textContent: isPublic ? "Public" : "Private" });
        // downloadable section
        setElementProperties(listResultsView, ".js-display-downloadable", { id: `${datasetId}-display-downloadable`});
        setElementProperties(listResultsView, ".js-editable-downloadable input", { id: `${resultDatasetId}-editable-downloadable`, checked: isDownloadable, dataset: { downloadable: dataset.is_downloadable } });
        setElementProperties(listResultsView, ".js-editable-downloadable label", { htmlFor: `${resultDatasetId}-editable-downloadable`, textContent: isDownloadable ? "Yes" : "No" });
        // if h5ad is not present, then disable input and set to "No"
        if (!hasH5ad) {
            setElementProperties(listResultsView, ".js-editable-downloadable input", { disabled: true });
            setElementProperties(listResultsView, ".js-editable-downloadable label", { textContent: "No" });
        }

        // organism section
        setElementProperties(listResultsView, ".js-display-organism span:last-of-type", { id: `${resultDatasetId}-display-organism`, textContent: organism });
        setElementProperties(listResultsView, ".js-editable-organism input", {value: organism });
        // owner section
        setElementProperties(listResultsView, ".js-display-owner span:last-of-type", { textContent: userName });
        setElementProperties(listResultsView, ".js-editable-owner input", { value: userName });
        // date added section
        setElementProperties(listResultsView, ".js-display-date-added span:last-of-type", { textContent: dateAdded });
        setElementProperties(listResultsView, ".js-editable-date-added input", { value: dateAdded });

        // annotation source section
        setElementProperties(listResultsView, ".js-display-annot-source span:last-of-type", { textContent: annotSource });
        setElementProperties(listResultsView, ".js-editable-annot-source input", { value: annotSource });

        // annotation version section
        setElementProperties(listResultsView, ".js-display-annot-version span:last-of-type", { textContent: annotVersion });
        setElementProperties(listResultsView, ".js-editable-annot-version input", { value: annotVersion });

        // pubmed id section
        const pubmedProp = pubmedId ? {
            // Not adding "icon-text" and "icon" classes because it aligns text and icon under the label
            innerHTML: `<a href="https://pubmed.ncbi.nlm.nih.gov/${pubmedId}" target="_blank">
                    <span>${pubmedId}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>`
        } : { textContent: "Not available" };
        pubmedProp.id = `${resultDatasetId}-display-pubmed-id`;

        setElementProperties(listResultsView, ".js-display-pubmed-id span:last-of-type", pubmedProp);
        setElementProperties(listResultsView, ".js-editable-pubmed-id input", { id: `${resultDatasetId}-editable-pubmed-id`, dataset: { originalVal: pubmedId }, value: pubmedId });

        // geo id section
        const geoProp = geoId ? {
            innerHTML: `<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${geoId}" target="_blank">
                    <span>${geoId}</span>
                    <i class="mdi mdi-open-in-new"></i>
                </a>`
        } : { textContent: "Not available" };
        geoProp.id = `${resultDatasetId}-display-geo-id`;

        setElementProperties(listResultsView, ".js-display-geo-id span:last-of-type", geoProp);
        setElementProperties(listResultsView, ".js-editable-geo-id input", { id: `${resultDatasetId}-editable-geo-id`, dataset: { originalVal: geoId }, value: geoId });

        // action buttons section
        setElementProperties(listResultsView, ".js-view-dataset", { value: shareId });
        setElementProperties(listResultsView, ".js-delete-dataset", { value: datasetId, dataset: { isOwner } });
        setElementProperties(listResultsView, ".js-collection-add-dataset", { dataset: { datasetId, isPublic } });
        setElementProperties(listResultsView, ".js-collection-remove-dataset", { dataset: { datasetId } });
        setElementProperties(listResultsView, ".js-download-dataset", { id: `${resultDatasetId}-download-dataset`, dataset: { datasetId, isDownloadable, hasH5ad } });
        setElementProperties(listResultsView, ".js-share-dataset", { value: shareId, dataset: { datasetId } });
        setElementProperties(listResultsView, ".js-edit-dataset", { value: datasetId, dataset: { datasetId } });
        setElementProperties(listResultsView, ".js-edit-dataset-save", { value: datasetId, dataset: { datasetId } });
        setElementProperties(listResultsView, ".js-edit-dataset-cancel", { value: datasetId, dataset: { datasetId } });

        setElementProperties(listResultsView, ".js-dataset-curator", { href: `./dataset_curator.html?dataset_id=${datasetId}`});
        setElementProperties(listResultsView, ".js-multigene-viewer", { href: `./multigene_curator.html?dataset_id=${datasetId}`});
        setElementProperties(listResultsView, ".js-compare-tool", { href: `./compare_datasets.html?dataset_id=${datasetId}`});
        setElementProperties(listResultsView, ".js-sc-workbench", { href: `./analyze_dataset.html?dataset_id=${datasetId}`});


        // dataset type section
        setElementProperties(listResultsView, ".js-display-dataset-type span:last-of-type", { textContent: datasetType });
        setElementProperties(listResultsView, ".js-editable-dataset-type input", { value: datasetType });
        // long description section
        setElementProperties(listResultsView, ".js-display-ldesc", { id: `${resultDatasetId}-display-ldesc-container` });
        setElementProperties(listResultsView, ".js-editable-ldesc textarea", { id: `${resultDatasetId}-editable-ldesc`, dataset: { originalVal: longDesc }, value: longDesc });

        // Append the cloned template to the results container
        resultsListDiv.appendChild(listResultsView);

        // add collection membership info to dataset
        const collections = dataset.layouts;
        for (const collectionString of collections) {
            const collection = JSON.parse(collectionString);
            const shareId = collection.share_id;
            const label = collection.label;

            const collectionListElt = document.querySelector(`#${resultDatasetId} .js-found-in-collections-list`);
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

        // EXTRA STUFF TO BOTH VIEWS
        addDownloadableInfoToDataset(datasetId, isDownloadable, hasH5ad);
        addVisibilityInfoToDataset(datasetId, isPublic);

        const datasetImageContainer = document.querySelector(`#${resultDatasetId}-figure img`);
        datasetImageContainer.src = previewImageUrl;

        // Add ldesc if it exists
        const ldescContainer = document.getElementById(`${resultDatasetId}-display-ldesc-container`);
        const ldescElt = document.createElement("p");
        ldescElt.id = `${resultDatasetId}-display-ldesc`;
        ldescElt.textContent = longDesc || "No description entered";
        ldescContainer.appendChild(ldescElt);


    }

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

    // Then apply each tooltip to the appropriate element for all elements with the data-tooltip-content attribute

    for (const actionElt of document.querySelectorAll(".js-action-links")) {
        const loopTooltips = [...tooltips];
        for (const classElt of actionElt.querySelectorAll("[data-tooltip-content]")) {
            applyTooltip(classElt, loopTooltips.shift());
        }
    }

    // Create tooltips for dataset view buttons
    const viewBtns = document.getElementsByClassName("js-view-btn");
    for (const classElt of viewBtns) {
        applyTooltip(classElt, createActionTooltips(classElt), "bottom");
    }

    // Hide/Remove some buttons if user is not owner
    updateDatasetListButtons();

    // Initiialize delete dataset popover for each delete button
    createDeleteDatasetConfirmationPopover();

    // All event listeners for all dataset elements
    addDatasetListEventListeners();
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
 * Displays an action note for a dataset.
 *
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} shareUrl - The URL to be copied to the clipboard.
 * @returns {void}
 */
const showDatasetActionNote = (datasetId, shareUrl) => {
    // sanitize shareUrl
    shareUrl = shareUrl.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

    const noteSelector = `#result-dataset-id-${datasetId} span.js-dataset-action-note`;
    const noteSelecterElt = document.querySelector(noteSelector);

    noteSelecterElt.innerHTML = copyToClipboard(shareUrl) ? "URL copied to clipboard" : `Failed to copy to clipboard. URL: ${shareUrl}`;
    noteSelecterElt.classList.remove("is-hidden");

    setTimeout(() => {
        noteSelecterElt.classList.add("is-hidden"); // TODO: add animate CSS with fade out
        noteSelecterElt.innerHTML = "";
    }, 5000);
}

/**
 * Submits a search request with the specified search terms and pagination options.
 *
 * @param {number} [page=1] - The page number to retrieve.
 * @returns {Promise<void>} - A promise that resolves when the search request is completed.
 */
const submitSearch = async (page=1) => {

    const searchTerms = document.getElementById("search-terms").value;

    // If this is the first time searching with terms, set the sort by to relevance
    if (searchTerms && firstSearch) {
        document.getElementById("sort-by").value = 'relevance';
        firstSearch = false;
    }

    const searchCriteria = {
        'session_id': CURRENT_USER.session_id,
        'search_terms': searchTerms,
        'sort_by': document.getElementById("sort-by").value
    };

    if (searchByCollection) {
        searchCriteria.layout_share_id = selected_dc_share_id;
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

const updateDatasetCollectionButtons = (collection=null) => {

    // If user is not the owner of the collection, remove the delete and rename buttons
    const isOwner = (collection && Boolean(collection.is_owner)) || false;

    const collectionRenameButton = document.getElementById("btn-rename-collection");
    const collectionDeleteButton = document.getElementById("btn-delete-collection");

    collectionRenameButton.classList.remove("is-hidden");
    collectionDeleteButton.classList.remove("is-hidden");
    collectionRenameButton.disabled = false;
    collectionDeleteButton.disabled = false;
    if (!isOwner) {
        collectionRenameButton.classList.add("is-hidden");
        collectionDeleteButton.classList.add("is-hidden");
        // disable the buttons
        collectionRenameButton.disabled = true;
        collectionDeleteButton.disabled = true;
    }
}

/**
 * Updates the dataset list buttons based on the dataset properties.
 */
const updateDatasetListButtons = () => {

    const datasetListElements = document.getElementsByClassName("js-dataset-list-element");
    for (const classElt of datasetListElements) {
        const datasetId = classElt.dataset.datasetId;

        //unhide all buttons
        for (const actionLinks of classElt.querySelectorAll(".js-action-links .control")) {
            actionLinks.classList.remove("is-hidden");
        }

        // If the dataset has no h5ad file, remove the download button since it cannot be downloaded regardless
        const downloadButton = classElt.querySelector("button.js-download-dataset");
        if (downloadButton) {
            const hasH5ad = parseBool(downloadButton.dataset.hasH5ad);
            if (!hasH5ad) {
                downloadButton.parentElement.remove();
            }
        }

        // If button still exists, update its visibility if the dataset is downloadable
        if (downloadButton) {
            const isDownloadable = parseBool(downloadButton.dataset.isDownloadable);
            downloadButton.parentElement.classList.toggle("is-hidden", !isDownloadable);
            downloadButton.disabled = !isDownloadable;
        }

        // The ability to edit and delete and dataset are currently paired
        const deleteButton = classElt.querySelector("button.js-delete-dataset");
        const editButton = classElt.querySelector("button.js-edit-dataset");

        const selectorBase = `#result-dataset-id-${datasetId}`;

        if (deleteButton) {
            // If user is not the owner of the dataset, remove the delete and edit buttons so user cannot manipulate
            // These will be regenerated when a search triggers processSearchResults
            if (deleteButton.dataset.isOwner === "false") {
                deleteButton.parentElement.remove();
                editButton.parentElement.remove()

                // Remove all editable elements to prevent editing in the DOM
                for (const editableElt of classElt.querySelectorAll(`${selectorBase} .js-editable-version`)) {
                    editableElt.classList.remove()
                }
            };
        }
    }
}


/**
 * Updates the "Add to Collection" and "Remove from Collection" buttons for each dataset element in the dataset list.
 * If a collection is provided, the buttons will be enabled or disabled based on the dataset's membership in the collection.
 * If the collection is a domain or no collection is provided, the buttons will be hidden and disabled.
 *
 * @param {Object} collection - The collection object. If null, the buttons will be hidden and disabled for all datasets.
 */
const updateDatasetAddRemoveToCollectionButtons = (collection=null) => {
    const datasetListElements = document.getElementsByClassName("js-dataset-list-element");
    for (const classElt of datasetListElements) {
        const datasetId = classElt.dataset.datasetId;

        const addToCollectionButton = classElt.querySelector("button.js-collection-add-dataset");
        const removeFromCollectionButton = classElt.querySelector("button.js-collection-remove-dataset");

        // is-hidden removed from previous loop

        addToCollectionButton.parentElement.classList.remove("is-hidden")
        removeFromCollectionButton.parentElement.classList.remove("is-hidden")
        addToCollectionButton.disabled = false;
        removeFromCollectionButton.disabled = false;

        // if dataset collection is domain, remove the "add to collection" and "remove from collection" buttons
        if (!collection || Boolean(collection.is_domain)) {
            addToCollectionButton.parentElement.classList.add("is-hidden")
            removeFromCollectionButton.parentElement.classList.add("is-hidden")
            // disable the buttons
            addToCollectionButton.disabled = true;
            removeFromCollectionButton.disabled = true;
            continue;
        }

        // if dataset is not in the currently selected collection, remove the "remove from collection" button (and vice versa)
        const datasetInCollection = collection.members.some(member => member.dataset_id === datasetId);
        if (!datasetInCollection) {
            removeFromCollectionButton.parentElement.classList.add("is-hidden")
            removeFromCollectionButton.disabled = true;
        }
    };
}

/* --- Entry point --- */
const handlePageSpecificLoginUIUpdates = async (event) => {

	// User settings has no "active" state for the sidebar
	document.getElementById("page-header-label").textContent = "Dataset Explorer";
	for (const elt of document.querySelectorAll("#primary-nav .menu-list a.is-active")) {
		elt.classList.remove("is-active");
	}

	document.querySelector("a[tool='search_datasets'").classList.add("is-active");

    const sessionId = CURRENT_USER.session_id;

	if (! sessionId ) {
        document.getElementById("not-logged-in-msg").classList.remove("is-hidden");
        document.getElementById("collection-management").classList.add("is-hidden");
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

    await Promise.all([
        loadOrganismList(),
        fetchDatasetCollections(datasetCollectionSelectCallback)
    ]);

    document.getElementById("dropdown-dc").classList.remove("is-right");    // Cannot see the dropdown if it is right aligned

    await submitSearch();

    // Set state of initial dataset buttons based on the initial dataset collection.
    const collection = flatDatasetCollectionData.find((collection) => collection.share_id === selected_dc_share_id);
    updateDatasetAddRemoveToCollectionButtons(collection);

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
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-table-view").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-table-view").classList.remove('is-dark');

    document.getElementById("results-table").classList.remove("is-hidden");
    document.getElementById("results-list-div").classList.add("is-hidden");
    document.getElementById("dataset-arrangement-c").classList.add("is-hidden");

    // Show pagination in case arrangement view hid the pagination
    for (const classElt of document.getElementsByClassName("pagination")) {
        classElt.classList.remove("is-invisible");
    }
    document.getElementById("count-label-c").classList.remove("is-hidden");

})

document.getElementById("btn-list-view-compact").addEventListener("click", () => {
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-list-view-compact").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-list-view-compact").classList.remove('is-dark');

    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.remove("is-hidden");
    document.getElementById("dataset-arrangement-c").classList.add("is-hidden");


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
    for (const classElt of document.getElementsByClassName("js-view-btn")) {
        classElt.classList.remove('is-gear-bg-secondary');
        classElt.classList.add('is-dark');
    }

    document.getElementById("btn-list-view-expanded").classList.add('is-gear-bg-secondary');
    document.getElementById("btn-list-view-expanded").classList.remove('is-dark');

    document.getElementById("results-table").classList.add("is-hidden");
    document.getElementById("results-list-div").classList.remove("is-hidden");
    document.getElementById("dataset-arrangement-c").classList.add("is-hidden");


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

// If checkbox is changed, set the searchByCollection flag
document.getElementById("filter-only-in-collection").addEventListener("change", (e) => {
    searchByCollection = e.currentTarget.checked;
    // If the checkbox is checked, change label accordingly
    e.currentTarget.closest(".field").querySelector("label").textContent = searchByCollection ? "Yes" : "No"
});

document.getElementById("btn-save-arrangement").addEventListener("click", async () => {

    const singlearrangementTiles = document.querySelectorAll("#dataset-arrangement-single .js-sortable-tile");
    const multiarrangementTiles = document.querySelectorAll("#dataset-arrangement-multi .js-sortable-tile");

    // Get the current grid states for all tiles and put in an object

    const layoutArrangement = {"single": {}, "multi": {}};  // dataset ids as keys

    for (const tile of singlearrangementTiles) {
        const datasetId = tile.dataset.datasetId;
        layoutArrangement.single[datasetId] = {
            "start_row": parseInt(tile.style.gridRowStart),
            "start_col": parseInt(tile.style.gridColumnStart),
            "grid_height": parseInt(tile.style.gridRowEnd.split(" ")[1]),
            "grid_width": parseInt(tile.style.gridColumnEnd.split(" ")[1])
        }
    }

    for (const tile of multiarrangementTiles) {
        const datasetId = tile.dataset.datasetId;
        layoutArrangement.multi[datasetId] = {
            "start_row": parseInt(tile.style.gridRowStart),
            "start_col": parseInt(tile.style.gridColumnStart),
            "grid_height": parseInt(tile.style.gridRowEnd.split(" ")[1]),
            "grid_width": parseInt(tile.style.gridColumnEnd.split(" ")[1])
        }
    }


    const data = await apiCallsMixin.saveDatasetCollectionArrangement(selected_dc_share_id, layoutArrangement)
    if (data.success) {
        createToast("Layout arrangement saved successfully", "is-success");
    } else {
        createToast("Failed to save layout arrangement");
    }
});

document.getElementById("btn-set-primary-collection").addEventListener("click", async () => {
    try {
        const data = await apiCallsMixin.setUserPrimaryDatasetCollection(selected_dc_share_id)
        if (data.success) {
            createToast("Primary collection set successfully", "is-success");

            Cookies.set('gear_default_domain', selected_dc_share_id);

            // Make button outlined to look "official"
            document.getElementById("btn-set-primary-collection").classList.remove("is-outlined");

        } else {
            throw new Error("Failed to set primary collection");
        }
    } catch (e) {
        createToast(e.message);
    }


});
