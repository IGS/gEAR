'use strict';

/* Given a passed-in layout_id, genereate a 2-dimensional tile-based grid object.
This uses Bulma CSS for stylings (https://bulma.io/documentation/layout/tiles/)
For the given layout, a single-gene grid and a multi-gene grid are generated.
*/

const plotlyPlots = ["bar", "line", "scatter", "tsne/umap_dynamic", "violin"];  // "tsne_dynamic" is a legacy option
const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];   // "tsne" is a legacy option

// Epiviz overrides the <script> d3 version when it loads so we save as a new variable to preserve it
const new_d3 = d3;

class TileGrid {

    constructor(shareId, type="layout", selector ) {
        this.shareId = shareId;
        this.type = type;
        this.datasets = []; // this.getDatasets();
        this.layout = {};   // this.getLayout();

        this.maxCols = 12 // highest number of columns in a row
        this.arrangementWidth = 1080;
        this.rowWidth = this.arrangementWidth / this.maxCols; // Split width into 12 columns
        this.colHeight = this.rowWidth * 4; // A unit of height for us is 4 units of width (to make a square)

        this.tiles = [];
        this.selector = selector;

    }

    /**
     * Adds all displays to the tile grid.
     * @returns {Promise<void>} A promise that resolves when all displays are added.
     */
    async addAllDisplays() {

        for (const dataset of this.datasets) {
            dataset.userDisplays = [];
            dataset.ownerDisplays = [];

            const {user: userDisplays, owner: ownerDisplays} = await apiCallsMixin.fetchDatasetDisplays(dataset.id);
            dataset.userDisplays = userDisplays;
            dataset.ownerDisplays = ownerDisplays;
        }
    }

    /**
     * Applies the tile grid layout to the specified selector element.
     *
     * @param {boolean} [isMulti=false] - Indicates whether the grid layout supports multiple tiles per grid slot.
     */
    applyTileGrid(isMulti=false) {
        const selector = this.selector;

        // Clear selector element
        const selectorElt = document.querySelector(selector);
        selectorElt.replaceChildren();

        if (this.type === "dataset") {
            const datasetTile = new DatasetTile(this, null, this.datasets[0], isMulti); // default display will be currentDisplayId
            this.applySingleTileGrid(datasetTile, selectorElt);
            this.tiles = [datasetTile];
            return;
        }

        const tiles = [];

        const layout = isMulti ? this.layout.multi : this.layout.single;

        // clear any no-displays message
        const noDisplaysElt = selectorElt.querySelector(".no-displays");
        if (noDisplaysElt) {
            noDisplaysElt.remove();
        }
        if (!layout || layout.length === 0) {
            // create element to say this dataset collection has no saved displays.
            const noDisplays = document.createElement("div");
            noDisplays.classList.add("no-displays", "notification", "is-warning", "has-text-centered");
            noDisplays.style.gridColumn = "1 / -1";
            const keyword = isMulti ? "multi" : "single";
            noDisplays.textContent = `No ${keyword}-gene displays were saved for this dataset collection. Add some in the Dataset Explorer.`;
            selectorElt.append(noDisplays);
            return;
        }

        for (const memberString of layout) {
            const member = JSON.parse(memberString);

            const dataset = this.datasets.find(d => d.id === member.dataset_id);
            if (!dataset) {
                console.warn(`Dataset with ID ${member.dataset_id} not found.`);
                continue;
            }
            const datasetTile = new DatasetTile(this, member, dataset, isMulti);
            tiles.push(datasetTile);
        }
        this.tiles = tiles;



        // sort by grid position
        this.tiles.sort((a, b) => a.tile.gridPosition - b.tile.gridPosition);

        // Legacy mode - if all tiles have startCol = 1, then we are in legacy mode
        // These layouts were generated only with a "width" property
        const legacyMode = this.tiles.every(tile => tile.tile.startCol === 1);

        // If in legacy mode, then we need to calculate the startCol and endCol and startRow and endRow
        if (legacyMode) {
            let currentCol = 1;
            let currentRow = 1;
            const maxEndCol = 13;   // 12 grid slots. CSS grid is 1-indexed, so 13 is the max
            for (const tile of this.tiles) {
                const tileWidth = Number(tile.tile.width);
                const tileHeight = Number(tile.tile.height);

                // If endCol is greater than 13, then this tile is in the next row
                if (currentCol + tileWidth > maxEndCol) {
                    currentCol = 1;
                    currentRow++;
                }

                // TODO: Eventually change "end" values to a "span" value
                tile.tile.startCol = currentCol;
                tile.tile.endCol = currentCol + tileWidth;
                tile.tile.startRow = currentRow;
                tile.tile.endRow = currentRow + tileHeight;

                currentCol += tileWidth;
            }
        }

        // Add grid template information to the selector element.
        // Max col is going to be 12 (since width is stored as 12-col format in db),
        // and max row is max(endRow)

        const maxRows = Math.max(...this.tiles.map(tile => tile.tile.endRow)) - 1; // Subtract 1 because endRow is first row tile is not in
        selectorElt.style.gridTemplateRows = `repeat(${maxRows}, ${this.colHeight}px)`;
        // if going from single dataset to layout, need to unset gridTemplateColumns. Otherwise the grid will be messed up
        selectorElt.style.gridTemplateColumns = "";

        // Build the CSS grid using the startRow, startCol, endRow, and endCol properties of each tile
        for (const datasetTile of this.tiles) {
            // Get HTML for the tile
            const tile = datasetTile.tile;
            const tileChildHTML = tile.html;

            // Add the tile to the selector element
            selectorElt.append(tileChildHTML);

            // Set the grid-area property of the tile. Must be added after the tile is appended to the DOM
            const tileElement = document.getElementById(`tile-${tile.tileId}`);
            tileElement.style.gridArea = `${tile.startRow} / ${tile.startCol} / ${tile.endRow} / ${tile.endCol}`;
        }
    }

    /**
     * Applies a single tile grid layout to the specified selector element.
     *
     * @param {Object} datasetTile - The DatasetTile object to apply the grid layout to.
     * @param {HTMLElement} selectorElt - The selector element to apply the grid layout to.
     * @param {boolean} [isZoomed=false] - Indicates whether the grid layout is a zoomed dataset.
     * @returns {void}
     */
    applySingleTileGrid(datasetTile, selectorElt, isZoomed=false) {


        selectorElt.style.gridTemplateColumns = `repeat(1, 1fr)`;
        selectorElt.style.gridTemplateRows = `repeat(1, min-content)`; // this prevents extra space at the bottom of each grid row

        // if zoomed, create new DatasetTile object with zoomed flag
        // SAdkins - I tried to use a "clone node and change ids" approach, but it was not working.
        const createZoomedTile = () => {
            const zoomedDatasetTile = new DatasetTile(this, datasetTile.display, datasetTile.dataset, datasetTile.type === "multi", true);
            // add gene symbol to zoomed tile
            zoomedDatasetTile.geneSymbol = datasetTile.geneSymbol;
            zoomedDatasetTile.currentDisplayId = datasetTile.currentDisplayId;
            zoomedDatasetTile.svgScoringMethod = datasetTile.svgScoringMethod;
            return zoomedDatasetTile;
        }

        const zoomedDatasetTile = createZoomedTile();

        const tile = isZoomed ? zoomedDatasetTile.tile : datasetTile.tile;

        const tileChildHTML = tile.html;

        // Add the tile to the selector element
        selectorElt.append(tileChildHTML);

        // Set the grid-area property of the tile. Must be added after the tile is appended to the DOM
        const tileElement = document.getElementById(`tile-${tile.tileId}`);
        tileElement.style.gridArea = "auto";

        if (isZoomed) {
            zoomedDatasetTile.renderDisplay(zoomedDatasetTile.geneSymbol, zoomedDatasetTile.currentDisplayId, zoomedDatasetTile.svgScoringMethod);
        }

    }

    /**
     * Retrieves the datasets associated with the tile grid.
     * @returns {Promise<Array>} A promise that resolves to an array of datasets.
     * @throws {Error} If there is an error fetching the dataset list information.
     */
    async getDatasets() {
        const args = this.type === "dataset" ? {permalink_share_id: this.shareId} : {layout_share_id: this.shareId};

        try {
            const data = await apiCallsMixin.fetchDatasetListInfo(args);
            return data['datasets'];
        } catch (error) {
            throw error;
        }
    }

    /**
     * Retrieves the layout of the tile grid.
     * @returns {Promise<Array>} A promise that resolves to an array representing the layout of the tile grid.
     * @throws {Error} If there is an error while fetching the layout.
     */
    async getLayout() {
        try {
            const data = await apiCallsMixin.fetchDatasetCollectionMembers(this.shareId);
            return data["layout_members"];
        } catch (error) {
            throw error;
        }
    }

    /**
     * Renders the displays for the given gene symbols in the tile grid.
     *
     * @param {string|string[]} geneSymbols - The gene symbol or an array of gene symbols to render displays for.
     * @param {boolean} [isMultigene=false] - Indicates whether the gene symbols represent multiple genes.
     * @param {string} svgScoringMethod - The SVG scoring method to use for rendering the displays.
     * @param {object} [projectionOpts={}] - The options for performing projection.
     * @param {string} [projectionOpts.patternSource] - The pattern source for projection.
     * @param {string} [projectionOpts.algorithm] - The algorithm for projection.
     * @param {string} [projectionOpts.gctype] - The gene correlation type for projection.
     * @returns {Promise<void>} - A promise that resolves when all displays have been rendered.
     * @throws {Error} - If geneSymbols is not provided or if an error occurs during rendering.
     */
    async renderDisplays(geneSymbols, isMultigene = false, svgScoringMethod, projectionOpts={}) {
        if (!geneSymbols) {
            throw new Error("Gene symbol or symbols are required to render displays.");
        }

        // Adapt geneSymbols to an array if it is not already
        let geneSymbolInput = geneSymbols;
        if (!isMultigene) {
            geneSymbolInput = Array.isArray(geneSymbols) ? geneSymbols : [geneSymbols];
        }

        // sort tiles by height, ascending.  This should help cases where taller plots render as same height as shorter plots
        this.tiles.sort((a, b) => a.tile.height - b.tile.height);

        // Sometimes fails to render due to OOM errors, so we want to try each tile individually
        // Orthology mapping also seems to fail due to file locking as well.
        this.tiles.map(async tile => {
            await tile.processTileForRenderingDisplay(projectionOpts, geneSymbolInput, svgScoringMethod);
        });

    }
};

class DatasetTile {
    constructor(thisTileGrid, display, dataset, isMulti=true, isZoomed=false) {
        this.display = display; // has the layout member info
        this.dataset = dataset; // Has dataset metadata info

        this.type = isMulti ? 'multi' : 'single';
        this.typeInt = isMulti ? 1 : 0;

        this.tile = this.generateTile();

        this.isZoomed = isZoomed;
        if (this.isZoomed) {
            this.tile.tileId = `zoomed-${this.tile.tileId}`;
        }

        this.parentTileGrid = thisTileGrid;

        // Set the end row and end col based on the startow and startCol
        this.tile.endRow = this.tile.startRow + this.tile.height,
        this.tile.endCol = this.tile.startCol + this.tile.width,

        this.tile.html = this.generateTileHTML();

        this.controller = new AbortController(); // Create new controller for new set of frames

        this.orthologs = null;  // Mapping of all orthologs for all gene symbol inputs for this dataset
        this.orthologsToPlot = null;    // A flattened list of all orthologs to plot

        this.geneSymbol = null;
        this.currentDisplayId = this.display?.display_id || this.addDefaultDisplay().then((displayId) => this.currentDisplayId = displayId);

        // Projection information
        // modeEnabled: boolean - Indicates whether projection mode is enabled for this tile
        // projectionId: string - The ID of the projection returned
        // projectionInfo: string - The information about the projection (like common genes, etc.)
        // performingProjection: boolean - Indicates whether a projection is currently being performed.  We do not want to waste resources by performing the same projection multiple times.
        // success: boolean - Indicates whether the projection was successful
        this.projectR = {modeEnabled: false, projectionId: null, projectionInfo: null, performingProjection: false, success: false};
    }

    /**
     * Generates a tile object based on the current dataset and tile type.
     * @returns {Object} The generated tile object.
     */
    generateTile() {
        if (this.display) {
            const { display_id, grid_position, grid_width, grid_height, start_row, start_col } = this.display;
            return {
                gridPosition: grid_position,
                width: grid_width,
                // Heights are not bound by the 12-spaced grid, so just use 1, 2, 3, etc.
                height: grid_height,
                startRow: start_row,
                startCol: start_col,
                tileId: `${display_id}-${grid_position}-${this.type}`,  // Alternatively could use "id" which is layout_displays.id
            };
        }

        // This was a single dataset, so use dataset ID instead tile ID
        return {
            gridPosition: 1,
            width: 12,
            height: 12,
            startRow: 1,
            startCol: 1,
            tileId: `${this.dataset.id}-1-${this.type}`,
        };

    }

    /**
     * Generates the HTML representation of a tile.
     *
     * @returns {DocumentFragment} The generated HTML fragment.
     */
    generateTileHTML() {
        const tile = this.tile;

        const template = document.getElementById('tmpl-tile-grid-tile');
        const tileHTML = template.content.cloneNode(true);

        // Set tile id & title
        const tileElement = tileHTML.querySelector('.js-tile');
        tileElement.id = `tile-${tile.tileId}`;

        const tileTitle = tileHTML.querySelector('.card-header-title');
        tileTitle.textContent = this.dataset.title;

        this.addDropdownInformation(tileElement, tileElement.id, this.dataset);

        return tileHTML;
    }

    /**
     * Enables the Project R mode.
     */
    enableProjectR() {
        this.projectR.modeEnabled = true;
    }

    /**
     * Adds a default display to the tile grid.
     * @returns {Promise<void>} A promise that resolves when the default display is added.
     */
    async addDefaultDisplay() {

        const dataset = this.dataset;
        try {
            const {default_display_id: defaultDisplayId} = await apiCallsMixin.fetchDefaultDisplay(dataset.id, this.typeInt);
            return defaultDisplayId;
        } catch (error) {
            //pass
            console.warn("Could not fetch default display for dataset", dataset.id);
        }
    }

    /**
     * Processes a tile for rendering display.
     * @param {Object} projectionOpts - The projection options.
     * @param {string} geneSymbolInput - The gene symbol input.
     * @param {string} svgScoringMethod - The SVG scoring method.
     * @returns {Promise<void>} - A promise that resolves when the rendering is complete.
     */
    async processTileForRenderingDisplay(projectionOpts, geneSymbolInput, svgScoringMethod) {
        const tileId = this.tile.tileId;
        const tileElement = document.getElementById(`tile-${tileId}`);

        // If projection mode is enabled, then perform projection
        if (this.projectR.modeEnabled) {
            if (!this.projectR.projectionId) {
                const {patternSource, algorithm, gctype} = projectionOpts;
                await this.getProjection(patternSource, algorithm, gctype);
            }

            if (!this.projectR.success) {
                return;
            }
        }

        // Clear the tile content so ortholog dropdown does not duplicate
        const optionContent = tileElement.querySelector('.content');
        optionContent.replaceChildren();

        if (this.projectR.modeEnabled) {
            if (this.projectR.projectionInfo) {
                this.renderProjectionInfo();
            }

            this.resizeCardImage();

            // Plot and render the display
            await this.renderDisplay(geneSymbolInput, null, svgScoringMethod);
            return;
        }

        // Resize the card image (plots) to accomodate the title header
        this.resizeCardImage();

        // If the dataset type is epiviz, then give warning that it hasn't been implemented yet
        if (this.dataset.dtype === "epiviz") {
            createCardMessage(tileId, "warning", "Epiviz datasets are not yet supported.");
            return;
        }

        // fail fast if no h5ad file
        if (!this.dataset.has_h5ad) {
            createCardMessage(tileId, "danger", "No h5ad file found for this dataset. Please contact the gEAR team.");
            return;
        }

        // Not projection mode, so get orthologs
        await this.getOrthologs(geneSymbolInput)

        // If no genes were found, then raise an error
        // This should never happen as geneSymbolInput should be a key in the orthologs object
        if (!this.orthologs || Object.keys(this.orthologs).length === 0) {

            createCardMessage(tileId, "danger", "No orthologs were mapped for this dataset. This should not have happened.");
            throw new Error("Should never happen. Please contact the gEAR team.");
        }

        // Make a flattened list of all orthologs to plot
        this.orthologsToPlot = Object.keys(this.orthologs).map(g => this.orthologs[g].sort()).flat();

        if (this.orthologsToPlot.length === 0) {
            const message = "Neither the given gene symbol(s) nor corresponding orthologs were found in this dataset.";
            // Fill in card-image with error message
            createCardMessage(tileId, "danger", message);
            // skip to the next tile
            return;
        }

        // Get the first ortholog for each gene (needed for multigene plots, and initial single-gene plots).
        const orthologs = Object.keys(this.orthologs).map(g => this.orthologs[g].sort()[0]).flat();

        // Render ortholog dropdown if in single-gene view and there is more than one ortholog for the gene
        if (this.type === "single" && this.orthologsToPlot.length > 1) {
            this.renderOrthologDropdown();
        }

        // Plot and render the display
        await this.renderDisplay(orthologs, null, svgScoringMethod);
    }

    /**
     * Retrieves orthologs for the given gene symbols.
     * @param {Array<string>} geneSymbols - The gene symbols to retrieve orthologs for.
     * @returns {Promise<void>} - A promise that resolves when the orthologs are fetched.
     */
    async getOrthologs(geneSymbols) {
        if (this.orthologs) {
            this.orthologs = null;
        }

        const geneOrganismId = CURRENT_USER.default_org_id || null;

        try {
            const data = await apiCallsMixin.fetchOrthologs(this.dataset.id, geneSymbols, geneOrganismId);

            this.orthologs = data.mapping;
        } catch (error) {
            logErrorInConsole(error);
        }
    }

    /**
     * Retrieves the projection for the given pattern source, algorithm, and gctype.
     * @param {string} patternSource - The pattern source.
     * @param {string} algorithm - The algorithm to use for projection.
     * @param {string} gctype - The gctype.
     * @returns {Promise<void>} A promise that resolves when the projection is retrieved.
     */
    async getProjection(patternSource, algorithm, gctype) {
        this.resetAbortController();
        const otherOpts = {}
        if (this.controller) {
            otherOpts.signal = this.controller.signal;
        }

        try {
            const data = await apiCallsMixin.checkForProjection(this.dataset.id, patternSource, algorithm);
            // If file was not found, put some loading text in the plot
            if (! data.projection_id) {
                createCardMessage(this.tile.tileId, "info", "Creating projection. This may take a few minutes.");
            }
            this.projection_id = data.projection_id || null;
        } catch (error) {
            logErrorInConsole(error);
            return;
        }

        this.projectR.performingProjection = true;

        try {
            const data = await apiCallsMixin.fetchProjection(this.dataset.id, this.projection_id, patternSource, algorithm, gctype, otherOpts);
            const message = data.message || null;
            if (data.success < 1) {
                // throw error with message
                throw new Error(message);
            }
            this.projectR.projectionId = data.projection_id;
            this.projectR.projectionInfo = data.message;
            this.projectR.success = true;

        } catch (error) {
            if (error.name == "CanceledError") {
                console.info("display draw canceled for previous request");
                return;
            }

            createCardMessage(this.tile.tileId, "danger", error);

            logErrorInConsole(error);
        } finally {
            this.projectR.performingProjection = false;
        }
    }

    /**
     * Adds the dataset title to the modal.
     *
     * @param {HTMLElement} modalHTML - The HTML element representing the modal.
     */
    addDatasetTitleToModal(modalHTML) {
        const modalContent = modalHTML.querySelector('.modal-content');
        const datasetTitle = modalContent.querySelector("h5");
        datasetTitle.replaceChildren();
        datasetTitle.textContent = this.dataset.title;
    }

    /**
     * Renders the projection information for the tile.
     */
    renderProjectionInfo() {
        const template = document.getElementById('tmpl-tile-grid-projection-info');
        const projectionInfoHTML = template.content.cloneNode(true);
        projectionInfoHTML.querySelector("p").textContent = this.projectR.projectionInfo;
        const tileElement = document.getElementById(`tile-${this.tile.tileId}`);
        const cardContent = tileElement.querySelector('.js-card-extras');

        cardContent.append(projectionInfoHTML);
    }

    /**
     * Renders the ortholog dropdown for the tile grid.
     */
    renderOrthologDropdown() {
        const orthoTemplate = document.getElementById('tmpl-tile-grid-ortholog-dropdown');
        const orthoHTML = orthoTemplate.content.cloneNode(true);

        const orthoElement = orthoHTML.querySelector('.js-ortholog-dropdown');
        orthoElement.dataset.datasetId = this.dataset.id;

        // Populate the dropdown-items with the orthologs
        const dropdownContent = orthoHTML.querySelector('.dropdown-content');
        dropdownContent.replaceChildren();

        for (const ortholog of this.orthologsToPlot) {
            const orthologItem = document.createElement("a");
            orthologItem.classList.add("dropdown-item");
            orthologItem.textContent = ortholog;

            // If ortholog is clicked, render display using that ortholog
            orthologItem.addEventListener("click", async (event) => {
                // Set clicked item as active
                const allItems = dropdownContent.querySelectorAll('.dropdown-item');
                for (const item of allItems) {
                    item.classList.remove("is-active");
                }
                orthologItem.classList.add("is-active");

                const ortholog = event.currentTarget.textContent;

                // close dropdown
                orthoElement.classList.remove('is-active');

                await this.renderDisplay([ortholog], this.currentDisplayId, this.svgScoringMethod);

            });

            dropdownContent.append(orthologItem);
        }

        // Make first dropdown item active
        const firstItem = dropdownContent.querySelector('.dropdown-item');
        firstItem.classList.add("is-active");

        // Add event listener to dropdown trigger
        orthoElement.querySelector(".dropdown-trigger button").addEventListener("click", (event) => {
            const item = event.currentTarget;
            const triggerElt = item.parentElement;

            // close other dropdowns
            for (const dropdown of document.querySelectorAll('#result-panel-grid .dropdown')) {
                if (triggerElt !== dropdown.querySelector('.dropdown-trigger')) {
                    dropdown.classList.remove('is-active');
                }
            };

            const arrow = item.querySelector('i');

            // Close dropdown if already active
            if (orthoElement.classList.contains('is-active')) {
                orthoElement.classList.remove('is-active');
                arrow.classList.replace('mdi-chevron-up', 'mdi-chevron-down');
                return;
            }

            orthoElement.classList.add('is-active');
            arrow.classList.replace('mdi-chevron-down', 'mdi-chevron-up');

        });

        // Click off dropdown to close
        document.addEventListener('click', (event) => {
            if (!event.target.closest('.dropdown')) {
                for (const dropdown of document.querySelectorAll('#result-panel-grid .dropdown')) {
                    dropdown.classList.remove('is-active');
                }
            }
        });

        // Add dropdown to tile
        const tileElement = document.getElementById(`tile-${this.tile.tileId}`);
        const cardContent = tileElement.querySelector('.js-card-extras');
        cardContent.append(orthoHTML);

        // Adjust height of card image to account for the dropdown
        document.querySelector(`#tile-${this.tile.tileId} .card-image`).style.height = "calc(100% - 80px)";
    }

    /**
     * Adds dropdown information to a tile element.
     *
     * @param {HTMLElement} tileElement - The tile element to add dropdown information to.
     * @param {string} tileId - The ID of the tile.
     * @param {object} dataset - The dataset object containing information for the dropdown items.
     */
    addDropdownInformation(tileElement, tileId, dataset) {

        const dropdownMenu = tileElement.querySelector(`#${tileId} .dropdown-menu`);
        dropdownMenu.id = `dropdown-menu-${tileId}`;
        const dropdownContent = dropdownMenu.querySelector(".dropdown-content");
        const dropdownItems = dropdownContent.querySelectorAll('.dropdown-item');

        const datasetId = dataset.id;
        const pubmedId = dataset.pubmed_id;
        const geoId = dataset.geo_id;
        const hasTarball = dataset.has_tarball;
        const hasH5ad = dataset.has_h5ad;
        const isDownloadable = dataset.is_downloadable;
        const links = dataset.links;

        for (const item of dropdownItems) {
            switch (item.dataset.tool) {
                case "display":
                    // Add event listener to dropdown item
                    item.addEventListener("click", (event) => {
                        // Create and open modal for all user and owner displays
                        this.renderChooseDisplayModal();    // TODO: Need to add after displays are retrieved
                        const modalElt = document.getElementById(`choose-display-modal-${this.tile.tileId}`);
                        openModal(modalElt);
                    });
                    break;
                case "info":
                    // Modal for dataset information
                    item.addEventListener("click", (event) => {
                        const existingModals = document.querySelectorAll('.js-infobox');
                        for (const modal of existingModals) {
                            modal.remove();
                        }
                        const modalHTML = this.createModalInfobox();
                        // Add modal to DOM
                        document.body.append(modalHTML);
                        const modalElt = document.getElementById(`infobox-modal-${this.tile.tileId}`);
                        openModal(modalElt);
                    });
                    break;
                case "publication":
                    // Link to publication if it exists
                    if (pubmedId) {
                        item.href = `http://www.ncbi.nlm.nih.gov/pubmed/?term=${pubmedId}`;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "geo":
                    // Link to GEO entry if it exists
                    if (geoId) {
                        item.href = `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${geoId}`;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "notes":
                    // Modal for notes
                    item.addEventListener("click", (event) => {
                        // Create and open modal for all user and owner displays
                        createToast("This feature is not yet implemented.", "is-warning");
                    });
                    break;
                case "single-cell":
                    // Redirect to single-cell analysis workbench
                    if (hasH5ad) {
                        const url = `./sc_workbench.html?dataset_id=${datasetId}`;
                        //const url = `./sc_workbench.html?dataset_id=${datasetId}`;
                        item.href = url;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "compare":
                    // Redirect to comparison tool
                    if (hasH5ad) {
                        const url = `./compare_datasets.html?dataset_id=${datasetId}`;
                        item.href = url;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "single-gene":
                    // Redirect to single-gene curation tool
                    // TODO: It would be cool to pass in the gene symbol to the curation tool to auto-populate the gene symbol field
                    if (hasH5ad) {
                        const url = `./dataset_curator.html?dataset_id=${datasetId}`;
                        item.href = url;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "multi-gene":
                    // Redirect to multi-gene curation tool
                    // TODO: It would be cool to pass in the gene symbols to the curation tool to auto-populate the gene symbol field
                    if (hasH5ad) {
                        const url = `./multigene_curator.html?dataset_id=${datasetId}`;
                        item.href = url;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "download-bundle":
                    // Download dataset bundle
                    if (hasTarball && isDownloadable) {
                        try {
                            const url = `./cgi/download_source_file.cgi?type=tarball&dataset_id=${datasetId}`;
                            item.href = url;
                        } catch (error) {
                            logErrorInConsole(error);
                            createToast("An error occurred while trying to download the dataset bundle.");
                        }

                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "download-h5ad":
                    // Download h5ad file
                    if (hasH5ad && isDownloadable) {
                        try {
                            const url = `./cgi/download_source_file.cgi?type=h5ad&dataset_id=${datasetId}`;
                            item.href = url;
                        } catch (error) {
                            logErrorInConsole(error);
                            createToast("An error occurred while trying to download the h5ad file.");
                        }
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "download-png":
                    // Handle when plot type is known
                    item.classList.add("is-hidden");
                    break;
                default:
                    console.warn(`Unknown dropdown item ${item.dataset.tool} for dataset ${datasetId}.`);
                    break;
            }

            // Remove "is-active" after click
            item.addEventListener("click", (event) => {
                const item = event.currentTarget;
                item.classList.remove("is-active");
            });
        }

        // Add links to the dropdown above the download-bundle item
        if (links.length) {
            const downloadBundle = dropdownContent.querySelector('.dropdown-item[data-tool="download-bundle"]');
            // create new download-item for each link
            for (const link of links) {
                const linkTemplate = document.getElementById('tmpl-tile-grid-dropdown-link');
                const linkHTML = linkTemplate.content.cloneNode(true);

                const linkItem = linkHTML.querySelector('.dropdown-item');
                const linkLabel = linkHTML.querySelector('.link-label');
                const linkA = linkHTML.querySelector('a');

                linkItem.classList.add(`dataset-link-${link.resource}`);    // In case there is custom CSS for this link
                linkA.href = link.url;
                linkLabel.textContent = link.label;

                downloadBundle.insertAdjacentElement("beforebegin", linkItem);
            }

            // add divider between links and download items
            const divider = document.createElement("hr");
            divider.classList.add("dropdown-divider");
            downloadBundle.insertAdjacentElement("beforebegin", divider);
        }

        if (this.isZoomed) {
            tileElement.querySelector('.js-expand-display').classList.add("is-hidden");
        } else {
            tileElement.querySelector('.js-shrink-display').classList.add("is-hidden");
        }
        tileElement.querySelector('.js-expand-display').addEventListener("click", (event) => {
            // Apply a zoomed-in display
            document.getElementById("result-panel-grid").classList.add("is-hidden");
            document.getElementById("zoomed-panel-grid").replaceChildren();
            document.getElementById("zoomed-panel-grid").classList.remove("is-hidden");

            // Apply single tile grid
            this.parentTileGrid.applySingleTileGrid(this, document.getElementById("zoomed-panel-grid"), true);

        });
        tileElement.querySelector('.js-shrink-display').addEventListener("click", (event) => {
            // Revert back to "#result-panel-grid" display
            document.getElementById("result-panel-grid").classList.remove("is-hidden");
            document.getElementById("zoomed-panel-grid").classList.add("is-hidden");
        });

        // Add event listener to dropdown trigger
        tileElement.querySelector("button.dropdown-trigger").addEventListener("click", (event) => {
            const item = event.currentTarget;
            if (!item.classList.contains('dropdown-trigger')) {
                return;
            }
            // close other dropdowns
            for (const dropdown of document.querySelectorAll('#result-panel-grid .dropdown')) {
                if (item !== dropdown.querySelector('.dropdown-trigger')) {
                    dropdown.classList.remove('is-active');
                }
            };

            item.closest(".dropdown").classList.toggle('is-active');
        });

        // Click off dropdown to close
        document.addEventListener('click', (event) => {
            if (!event.target.closest('.dropdown')) {
                for (const dropdown of document.querySelectorAll('#result-panel-grid .dropdown')) {
                    dropdown.classList.remove('is-active');
                }
            }
        });
    }

    /**
     * Adds a modal display section title to the specified element.
     *
     * @param {HTMLElement} element - The element to which the title will be added.
     * @param {string} titleText - The text content of the title.
     */
    addModalDisplaySectionTitle(element, titleText) {
        if (!element.children.length) {
            return;
        }
        const title = document.createElement("p");
        title.classList.add("has-text-weight-bold", "is-underlined", "column", "is-full");
        title.textContent = titleText;
        element.prepend(title);
    }

    /**
     * Creates the HTML for the modal used in the tile grid to choose display options.
     * @returns {DocumentFragment} The HTML fragment representing the modal.
     */
    createModalHTML() {
        const modalTemplate = document.getElementById('tmpl-tile-grid-choose-display-modal');
        const modalHTML = modalTemplate.content.cloneNode(true);

        const modalDiv = modalHTML.querySelector('.modal');
        modalDiv.id = `choose-display-modal-${this.tile.tileId}`;

        return modalHTML;
    }

    /**
     * Creates a modal infobox for the tile.
     * @returns {DocumentFragment} The cloned infobox HTML.
     */
    createModalInfobox() {
        const infoboxTemplate = document.getElementById('tmpl-tile-infobox');
        const infoboxHTML = infoboxTemplate.content.cloneNode(true);

        const modalDiv = infoboxHTML.querySelector('.modal');
        modalDiv.id = `infobox-modal-${this.tile.tileId}`;

        modalDiv.querySelector('.modal-card-body .js-infobox-title').textContent = this.dataset.title;
        modalDiv.querySelector('.modal-card-body .js-infobox-ldesc').innerHTML = this.dataset.ldesc;

        // Some .js-infobox-ldesc may have embedded table data. If so, need to display this properly
        const tableElts = modalDiv.querySelectorAll('.modal-card-body .js-infobox-ldesc table');
        for (const tableElt of tableElts) {
            tableElt.classList.add("table", "is-fullwidth");
        }


        const img = modalDiv.querySelector('.modal-card-body .js-infobox-schematic');
        if (this.dataset.schematic_image) {
            img.src = this.dataset.schematic_image;
        } else {
            // remove image element
            img.remove();
        }

        // Close button event listener
        const closeButton = modalDiv.querySelector(".delete");
        closeButton.addEventListener("click", (event) => {
            closeModal(modalDiv);
        });
        const modalBackground = modalDiv.querySelector(".modal-background");
        modalBackground.addEventListener("click", (event) => {
            closeModal(modalDiv);
        });

        return infoboxHTML;
    }

    /**
     * Retrieves all displays from the dataset based on the type of tile grid.
     *
     * @returns {Object} An object containing user displays and owner displays.
     */
    getAllDisplays() {
        const filterKey = this.type === "single" ? "gene_symbol" : "gene_symbols";
        const userDisplays = this.dataset.userDisplays.filter((d) => d.plotly_config.hasOwnProperty(filterKey));
        const ownerDisplays = this.dataset.ownerDisplays.filter((d) => d.plotly_config.hasOwnProperty(filterKey));

        if (this.type === "single") {
            // Add userEpivizDisplays to userDisplays...
            const userEpivizDisplays = this.dataset.userDisplays.filter((d) => d.plot_type === "epiviz");
            const ownerEpivizDisplays = this.dataset.ownerDisplays.filter((d) => d.plot_type === "epiviz");
            userDisplays.push(...userEpivizDisplays);
            ownerDisplays.push(...ownerEpivizDisplays);
        }

        return { userDisplays, ownerDisplays };
    }

    /**
     * Removes existing modals from the document.
     */
    removeExistingModals() {
        const existingModals = document.querySelectorAll('.js-choose-display-modal');
        for (const modal of existingModals) {
            modal.remove();
        }
    }

    /**
     * Renders a modal to choose a display from the user or owner display lists.
     *
     * @returns {void}
     */
    renderChooseDisplayModal() {

        this.removeExistingModals();
        const modalHTML = this.createModalHTML();
        const modalDiv = modalHTML.querySelector('.modal');

        this.addDatasetTitleToModal(modalHTML);

        const { userDisplays, ownerDisplays } = this.getAllDisplays();

        const modalContent = modalHTML.querySelector('.modal-content');
        // Add user and owner displays
        const userDisplaysElt = modalContent.querySelector(".js-modal-user-displays");
        userDisplaysElt.replaceChildren();
        const ownerDisplaysElt = modalContent.querySelector(".js-modal-owner-displays");
        ownerDisplaysElt.replaceChildren();

        // Did it this way so we didn't have to pass async/await up the chain
        Promise.allSettled([
            this.renderChooseDisplayModalDisplays(userDisplays, userDisplaysElt),
            this.renderChooseDisplayModalDisplays(ownerDisplays, ownerDisplaysElt)
        ]).then(() => {

            this.addModalDisplaySectionTitle(userDisplaysElt, "Your Displays");
            this.addModalDisplaySectionTitle(ownerDisplaysElt, "Displays by Dataset Owner");


            const currentDisplayElt = modalContent.querySelector(`.js-modal-display[data-display-id="${this.currentDisplayId}"]`);

            // remove tag from all other displays
            const allDisplayElts = modalContent.querySelectorAll(".js-modal-display");
            for (const displayElt of allDisplayElts) {
                displayElt.classList.remove("is-selected");
            }

            // add tag to the currently selected display
            if (currentDisplayElt) {
                currentDisplayElt.classList.add("is-selected");
            }

            // Add event listeners to all display elements to render the display
            for (const displayElt of allDisplayElts) {
                displayElt.addEventListener("click", (event) => {
                    const displayElement = event.currentTarget;
                    const displayId = parseInt(displayElement.dataset.displayId);
                    // Render display
                    if (!this.svgScoringMethod) this.svgScoringMethod = "gene";
                    this.renderDisplay(this.geneSymbol, displayId, this.svgScoringMethod);

                    // Close modal
                    closeModal(modalDiv);
                });
            }
        });


        // Close button event listener
        const closeButton = modalDiv.querySelector(".modal-close");
        closeButton.addEventListener("click", (event) => {
            closeModal(modalDiv);
        });
        const modalBackground = modalDiv.querySelector(".modal-background");
        modalBackground.addEventListener("click", (event) => {
            closeModal(modalDiv);
        });

        // Add modal to DOM
        document.body.append(modalHTML);

    }

    /**
     * Renders the choose display modal with the given displays.
     *
     * @param {Array} displays - The array of displays to render.
     * @param {HTMLElement} displayElt - The element to append the rendered displays to.
     * @returns {Promise<void>} - A promise that resolves when the rendering is complete.
     */
    async renderChooseDisplayModalDisplays(displays, displayElt) {
        const displayTemplate = document.getElementById('tmpl-tile-grid-choose-display-modal-display');
        // Add user displays
        for (const display of displays) {
            const displayHTML = displayTemplate.content.cloneNode(true);

            const displayElement = displayHTML.querySelector('.js-modal-display');
            displayElement.dataset.displayId = display.id;
            displayElement.dataset.datasetId = this.dataset.id;

            // Add display image
            let displayUrl = "";
            try {
                // NOTE: SVGs are colorless
                displayUrl = await apiCallsMixin.fetchDatasetDisplayImage(this.dataset.id, display.id);
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
            displayType.textContent = display.plot_type;

            displayElt.append(displayHTML);
        }
    }

    /**
     * Renders the display for a given gene symbol.
     * @param {string} geneSymbolInput - The gene symbol(s) to render the display for.
     * @param {string|null} displayId - The ID of the display to render. If null, the default display ID will be used.
     * @param {string} svgScoringMethod - The SVG scoring method to use.
     * @throws {Error} If geneSymbol is not provided.
     * @returns {Promise<void>} A promise that resolves when the display is rendered.
     */
    async renderDisplay(geneSymbolInput, displayId=null, svgScoringMethod="gene") {
        if (!geneSymbolInput) {
            throw new Error("Gene symbol or symbols are required to render this display.");
        }

        // Store gene symbol for future use (i.e. changing display, etc.)
        // Since this method manipulates the geneSymbolInput, we need to store the original input
        this.geneSymbol = geneSymbolInput;

        createCardMessage(this.tile.tileId, "info", "Loading display...");

        if (displayId === null) {
            displayId = this.currentDisplayId;  // this should be the defaultDisplayId
        };

        // For use for zoomed-in views
        this.currentDisplayId = displayId;

        geneSymbolInput = this.type === "single" ? geneSymbolInput[0] : geneSymbolInput;

        this.svgScoringMethod = svgScoringMethod;

        this.resetAbortController();
        const otherOpts = {}
        if (this.controller) {
            otherOpts.signal = this.controller.signal;
        }

        const filterKey = this.type === "single" ? "gene_symbol" : "gene_symbols";

        // find the display config in the user or owner display lists
        let userDisplay = this.dataset.userDisplays.find((d) => d.id === displayId && d.plotly_config.hasOwnProperty(filterKey));
        let ownerDisplay = this.dataset.ownerDisplays.find((d) => d.id === displayId && d.plotly_config.hasOwnProperty(filterKey));

        // There is a chance this display is a member of a dataset collection
        // but not owned by the current user or dataset owner (i.e. curator account made it)
        let layouts = [];
        if (this.parentTileGrid.layout) {
            layouts = this.type === "single" ? this.parentTileGrid.layout.single : this.parentTileGrid.layout.multi;
        }
        const layoutDisplay = layouts.find((d) => JSON.parse(d).display_id === displayId);

        // Try epiviz display if no plotly display was found
        if (this.type === "single") {
            if (!userDisplay) userDisplay = this.dataset.userDisplays.find((d) => d.id === displayId && d.plot_type === "epiviz");

            if (!ownerDisplay) ownerDisplay = this.dataset.ownerDisplays.find((d) => d.id === displayId && d.plot_type === "epiviz");
        }

        // add console warning if default display id was not found in the user or owner display lists
        if (!userDisplay && !ownerDisplay && !layoutDisplay) {
            // This can happen if the display ID for a layout member is owned by a different user that is not the dataset owner.
            console.warn(`Selected display id ${this.currentDisplayId} for dataset ${this.dataset.title} was not found. Will show first available.`);

            // last chance... if still no display config (i.e default display was not found), then use the first display config
            if (!userDisplay) userDisplay = this.dataset.userDisplays.find((d) => d.plotly_config.hasOwnProperty(filterKey));
            if (!ownerDisplay) ownerDisplay = this.dataset.ownerDisplays.find((d) => d.plotly_config.hasOwnProperty(filterKey));
        }

        // if the display config was not found, then do not render
        if (!userDisplay && !ownerDisplay && !layoutDisplay) {
            console.warn(`Display config for dataset ${this.dataset.title} was not found.`)
            // Let the user know that the display config was not found
            const message = `This dataset has no viewable curations for this view. Create a new curation in the ${this.type === "single" ? "Single-gene" : "Multi-gene"} Curator to view this dataset.`;
            createCardMessage(this.tile.tileId, "warning", message);
            return;
        }

        // if the display config was found, then render
        let display = userDisplay || ownerDisplay;

        // if the display is a layout member, then need to retrieve the actual display
        if (!display) {
            console.warn(`Shown display for dataset ${this.dataset.title} is coming from a dataset collection member that was not owned by the current user or dataset owner.`)
            const layoutDisplayId = JSON.parse(layoutDisplay).display_id;
            display = await apiCallsMixin.fetchDisplay(layoutDisplayId);
        }

        this.currentDisplayId = display.id;

        // handle legacy plot types
        if (display.plot_type === "tsne_dynamic") {
            display.plot_type = "tsne/umap_dynamic";
        }
        if (display.plot_type === "tsne") {
            display.plot_type = "tsne_static";
        }

        // Add gene or genes to plot config
        if (this.type === "multi") {
            display.plotly_config.gene_symbols = geneSymbolInput;
        } else {
            display.plotly_config.gene_symbol = geneSymbolInput;
        }

        // if projection ran, add the projection info to the plotly config
        if (this.projectR.modeEnabled && this.projectR.projectionId) {
            display.plotly_config.projection_id = this.projectR.projectionId;
        }

        try {
            if (plotlyPlots.includes(display.plot_type)) {
                await this.renderPlotlyDisplay(display, otherOpts);
            } else if (scanpyPlots.includes(display.plot_type)) {
                await this.renderScanpyDisplay(display, otherOpts);

                // Determine how "download_png" is handled for scanpy plots
                const downloadPNG = document.querySelector(`#tile-${this.tile.tileId} .dropdown-item[data-tool="download-png"]`);
                if (downloadPNG) {

                    // If I use the existing "download image" button after switching displays, all previous tsne-static displays will
                    // also be downloaded becuase event listeners are not removed. So, I will remove the button and re-add it.
                    // Source -> https://stackoverflow.com/a/9251864

                    const newDownloadPNG = downloadPNG.cloneNode(true);
                    downloadPNG.parentNode.replaceChild(newDownloadPNG, downloadPNG);

                    newDownloadPNG.classList.remove("is-hidden");
                    newDownloadPNG.addEventListener("click", async (event) => {
                        // get the download URL
                        await this.getScanpyPNG(display);
                    });

                }

            } else if (display.plot_type === "svg") {
                await this.renderSVG(display, this.svgScoringMethod, otherOpts);
            } else if (display.plot_type === "epiviz") {
                await this.renderEpivizDisplay(display, otherOpts);
            } else if (this.type === "multi") {
                await this.renderMultiGeneDisplay(display, otherOpts);
            } else {
                throw new Error(`Display config for dataset ${this.dataset.id} has an invalid plot type ${display.plot_type}.`);
            }
        } catch (error) {
            // we want to ensure other plots load even if one fails
            logErrorInConsole(error);
            // Fill in card-image with error message
            createCardMessage(this.tile.tileId, "danger", error.message);
        }

    }

    /**
     * Renders the Epiviz display on the tile grid.
     * @param {Object} display - The display object containing dataset_id and plotly_config.
     * @returns {Promise<void>} - A promise that resolves when the rendering is complete.
     */
    async renderEpivizDisplay(display, otherOpts) {
        const datasetId = display.dataset_id;
        const {gene_symbol: geneSymbol} = display.plotly_config;

        createCardMessage(this.tile.tileId, "warning", "Epiviz displays have not been implemented yet.");
        return;

        let genome = null;
        const genesTrack = display.plotly_config.tracks["EPIVIZ-GENES-TRACK"];
        if (genesTrack.length > 0) {
            const gttrack = genesTrack[0];
            genome = gttrack.measurements ? gttrack.measurements[0].id : gttrack.id[0].id;
        }

        // Get data and set up the image area
        const data = await apiCallsMixin.fetchEpivizDisplay(datasetId, geneSymbol, genome, otherOpts);
        if (data.hasOwnProperty("success") && data.success === -1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }

        const extendRangeRatio = 10;

        // generate the epiviz panel + tracks
        const epiviznav = document.querySelector(`#epiviznav_${this.tile.tileId}`);
        const plotContainer = document.querySelector(`#tile-${this.tile.tileId} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned

        if (!epiviznav) {
            // epiviz container already exists, so only update gneomic position in the browser

            plotContainer.replaceChildren();    // erase plot
            try {
                plotContainer.append(this.renderEpivizTemplate(data, display.plotly_config, extendRangeRatio));
            } catch (error) {
                logErrorInConsole(error);
                throw new Error(`Could not render Epiviz display. Please contact gEAR support`);
            }
            return;
        }
        const nStart = epivizNavStart(data, extendRangeRatio);
        const nEnd = epivizNavEnd(data, extendRangeRatio);

        // epiviz container already exists, so only update gneomic position in the browser
        epiviznav.setAttribute("chr", data.chr);
        epiviznav.setAttribute("start", nStart);
        epiviznav.setAttribute("end", nEnd);
        epiviznav.range = epiviznav.getGenomicRange(data.chr, nstart, nend);    // function is imported from epiviz JS
    }

    /**
     * Renders an Epiviz template with the provided data and configuration.
     * @param {Object} data - The data to be rendered in the template.
     * @param {Object} plotConfig - The configuration for the plot.
     * @param {number} extendRangeRatio - The ratio by which to extend the range.
     * @returns {HTMLElement} - The rendered Epiviz template.
     */
    renderEpivizTemplate(data, plotConfig, extendRangeRatio) {
        const template = document.getElementById('tmpl-epiviz-container');
        const epivizHTML = template.content.cloneNode(true);
        // Add in properties to the epiviz container
        const epivizContainer = epivizHTML.querySelector('.epiviz-container');
        epivizContainer.id = `epiviz_${this.tile.tileId}`;
        const epivizDataSource = epivizHTML.querySelector('epiviz-data-source');
        epivizDataSource.id = `${this.tile.tileId}epivizds`;
        epivizDataSource.setAttribute("provider-url", plotConfig.dataserver);

        const epivizNavigation = epivizHTML.querySelector('epiviz-navigation');
        epivizNavigation.id = `${this.tile.tileId}_epiviznav`;
        // the chr, start and end should come from query - map gene to genomic position.
        epivizNavigation.setAttribute("chr", data.chr);
        epivizNavigation.setAttribute("start", epivizNavStart(data, extendRangeRatio));
        epivizNavigation.setAttribute("end", epivizNavEnd(data, extendRangeRatio));
        epivizNavigation.setAttribute("viewer", `/epiviz.html?dataset_id=${this.dataset.id}&chr=${data.chr}&start=${data.start}&end=${data.end}`);
        epivizNavigation.innerHTML(this.renderEpivizTracks(plotConfig));
        return epivizHTML;
    }

    /**
     * Renders Epiviz tracks based on the provided plot configuration.
     * @param {Object} plotConfig - The plot configuration object.
     * @returns {string} - The HTML template for the Epiviz tracks.
     */
    renderEpivizTracks(plotConfig) {
        //Create the tracks
        let epivizTracksTemplate = "";
        for (const track in plotConfig.tracks) {
            const trackConfig = plotConfig.tracks[track];
            trackConfig.forEach((tc) => {
                let tempTrack = `<${track} slot='charts' `;
                tempTrack += Object.keys(tc).includes("id") ? ` dim-s='${JSON.stringify(tc.id)}' ` : ` measurements='${JSON.stringify(tc.measurements)}' `;

                if (tc.colors != null) {
                    tempTrack += ` chart-colors='${JSON.stringify(tc.colors)}' `;
                }

                if (tc.settings != null) {
                    tempTrack += ` chart-settings='${JSON.stringify(tc.settings)}' `;
                }

                tempTrack += ` style='min-height:200px;'></${track}> `;

                epivizTracksTemplate += tempTrack;
            });
        }

        return epivizTracksTemplate;
    }

    /**
     * Renders the multi-gene display.
     *
     * @param {Object} display - The display object.
     * @param {Object} otherOpts - Other options.
     * @returns {Promise<void>} - A promise that resolves when the rendering is complete.
     * @throws {Error} - If there is an error fetching the data or rendering the plot.
     */
    async renderMultiGeneDisplay(display, otherOpts) {

        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        // Get data and set up the image area
        const data = await apiCallsMixin.fetchDashData(datasetId, analysisObj, plotType, plotConfig, otherOpts);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {plot_json: plotJson} = data;

        const plotContainer = document.querySelector(`#tile-${this.tile.tileId} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const plotlyPreview = document.createElement("div");
        plotlyPreview.classList.add("container");
        plotlyPreview.id = `tile-${this.tile.tileId}-plotly-preview`;
        plotContainer.append(plotlyPreview);
        Plotly.purge(plotlyPreview.id); // clear old Plotly plots

        if (!plotJson) {
            console.warn(`Could not retrieve plot information for dataset display ${display.id}. Cannot make plot.`);
            return;
        }

        if (plotType === 'heatmap') {
            // These modify the plotJson object in place
            // TODO: Adjust these functions
            adjustExpressionColorbar(plotJson.data);
            adjustClusterColorbars(plotJson.data);
        }

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const expressionDisplayConf = postPlotlyConfig.expression;
        const customConfig = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "config");
        Plotly.newPlot(plotlyPreview.id , plotJson.data, plotJson.layout, customConfig);
        // ! Occasionally get a "something went wrong with axis scaling" error. This seems to arise if the colorbar is too small

        const customLayout = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "layout");
        Plotly.relayout(plotlyPreview.id , customLayout);


    }

    /**
     * Renders a Plotly display on the tile grid.
     *
     * @param {Object} display - The display object containing the plotly configuration.
     * @param {Object} otherOpts - Additional options for rendering the plotly display.
     * @returns {Promise<void>} - A promise that resolves when the plotly display is rendered.
     * @throws {Error} - If there is an error fetching the plotly data or rendering the plot.
     */
    async renderPlotlyDisplay(display, otherOpts) {
        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        // Get data and set up the image area
        const data = await apiCallsMixin.fetchPlotlyData(datasetId, analysisObj, plotType, plotConfig, otherOpts);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {plot_json: plotJson} = data;

        const plotContainer = document.querySelector(`#tile-${this.tile.tileId} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const plotlyPreview = document.createElement("div");
        plotlyPreview.classList.add("container");
        plotlyPreview.id = `tile-${this.tile.tileId}-plotly-preview`;
        plotContainer.append(plotlyPreview);
        Plotly.purge(plotlyPreview.id); // clear old Plotly plots

        if (!plotJson) {
            console.warn(`Could not retrieve plot information for dataset display ${display.id}. Cannot make plot.`);
            return;
        }
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const expressionDisplayConf = postPlotlyConfig.expression;
        const customConfig = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "config");
        Plotly.newPlot(plotlyPreview.id, plotJson.data, plotJson.layout, customConfig);
        const customLayout = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "layout");
        Plotly.relayout(plotlyPreview.id, customLayout);
    }

    /**
     * Renders the Scanpy display on the tile grid.
     *
     * @param {Object} display - The display object containing the dataset and plot information.
     * @param {Object} otherOpts - Additional options for rendering the display.
     * @returns {Promise<void>} - A promise that resolves when the display is rendered.
     * @throws {Error} - If there is an error fetching the image data or if the image data is not available.
     */
    async renderScanpyDisplay(display, otherOpts) {

        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        const tileElement = document.getElementById(`tile-${this.tile.tileId}`);
        if (!this.isZoomed) {
            plotConfig.grid_spec = tileElement.style.gridArea   // add grid spec to plot config
            if (plotConfig.grid_spec === "auto") delete plotConfig.grid_spec;   // single dataset grid spec
        }

        const plotContainer = document.querySelector(`#tile-${this.tile.tileId} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        const tsnePreview = document.createElement("img");
        tsnePreview.classList.add("image", "is-fullwidth");
        tsnePreview.id = `tile-${this.tile.tileId}-tsne-preview`;
        plotContainer.append(tsnePreview);

        const data = await apiCallsMixin.fetchTsneImage(datasetId, analysisObj, plotType, plotConfig, otherOpts);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {image} = data;

        if (!image) {
            console.warn(`Could not retrieve plot image for dataset display ${display.id}. Cannot make plot.`);
            return;
        }

        const blob = await fetch(`data:image/webp;base64,${image}`).then(r => r.blob());

        // decode base64 image and set as src
        tsnePreview.src = URL.createObjectURL(blob);

        tsnePreview.onload = () => {
            // Revoke the object URL to free up memory
            // ! This does prevent right-click saving though
            //URL.revokeObjectURL(tsnePreview.src);
        }
        return;
    }


    /**
     * Retrieves a PNG image for a given display and initiates its download.
     *
     * @async
     * @param {Object} display - The display object containing information about the dataset and plot configuration.
     * @returns {Promise<void>} - A promise that resolves when the PNG image is downloaded.
     * @throws {Error} - If the image retrieval is unsuccessful or encounters an unknown error.
     */
    async getScanpyPNG(display) {
        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.analysis_id ? {id: display.analysis_id} : display.analysis || null;
        const plotType = display.plot_type;
        const geneSymbol = display.plotly_config.gene_symbol;

        // deep copy plotly_config to avoid modifying the original
        const plotConfig = JSON.parse(JSON.stringify(display.plotly_config));
        plotConfig.high_dpi = true;

        const tileElement = document.getElementById(`tile-${this.tile.tileId}`);
        if (!this.isZoomed) {
            plotConfig.grid_spec = tileElement.style.gridArea   // add grid spec to plot config
            if (plotConfig.grid_spec === "auto") delete plotConfig.grid_spec;   // single dataset grid spec
        }

        const data = await apiCallsMixin.fetchTsneImage(datasetId, analysisObj, plotType, plotConfig);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {image} = data;
        if (!image) {
            console.warn(`Could not retrieve downloadable image for dataset display ${display.id}.`);
            return;
        }

        const blob = await fetch(`data:image/png;base64,${image}`).then(r => r.blob());
        const download = URL.createObjectURL(blob);

        // create a hidden element that will be clicked to download the PNG
        const hiddenLink = document.createElement('a');
        document.body.appendChild(hiddenLink);
        hiddenLink.classList.add("is-hidden");

        // download URL
        hiddenLink.download = `${this.dataset.id}_${geneSymbol}_${display.plot_type}.png`;
        hiddenLink.href = download;

        hiddenLink.setAttribute('target', '_blank');

        // click the hidden link to download the PNG
        hiddenLink.click();

        // save memory (but breaks download)
        URL.revokeObjectURL(download);
        hiddenLink.remove();
    }


    /**
     * Renders an SVG for the display.
     *
     * @param {Object} display - The display object.
     * @param {string} [svgScoringMethod="gene"] - The SVG scoring method.
     * @param {Object} otherOpts - Other options.
     * @returns {Promise<void>} - A promise that resolves when the SVG is rendered.
     * @throws {Error} - If there is an error rendering the SVG.
     */
    async renderSVG(display, svgScoringMethod="gene", otherOpts) {
        const datasetId = display.dataset_id;
        const plotConfig = display.plotly_config;
        const {gene_symbol: geneSymbol, projection_id: projectionId} = plotConfig;

        const data = await apiCallsMixin.fetchSvgData(datasetId, geneSymbol, projectionId, otherOpts)
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const plotContainer = document.querySelector(`#tile-${this.tile.tileId} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        colorSVG(data, plotConfig.colors, datasetId, this.tile.tileId, svgScoringMethod);

    }

    /**
     * Resets the AbortController and cancels any previous axios requests.
     * Creates a new AbortController for a new set of frames.
     */
    resetAbortController() {
        if (this.controller && !this.projectR.performingProjection) {
            this.controller.abort(); // Cancel any previous axios requests (such as drawing plots for a previous dataset)
        }
        this.controller = new AbortController(); // Create new controller for new set of frames

    }

    /**
     * Resizes the card image based on the sibling card header height.
     */
    resizeCardImage() {
        const cardImage = document.querySelector(`#tile-${this.tile.tileId} .card-image`);
        // resize based on the sibling .card-header height
        const cardHeader = document.querySelector(`#tile-${this.tile.tileId} .card-header`);
        const cardExtras = document.querySelector(`#tile-${this.tile.tileId} .js-card-extras`); // this is used in projections
        let offsetHeight = cardHeader.offsetHeight;
        if (cardExtras) {
            offsetHeight += cardExtras.offsetHeight;
        }
        cardImage.style.height = `calc(100% - ${offsetHeight}px)`;

    }
}

/**
 * Color the SVG based on the chart data and plot configuration.
 * @param {Object} chartData - The data for the chart.
 * @param {Object} plotConfig - The configuration for the plot.
 * @param {string} datasetId - The ID of the dataset.
 * @param {string} tileId - The ID of the tile.
 * @param {string} [svgScoringMethod="gene"] - The scoring method for coloring the SVG.
 * @returns {Promise<void>} - A promise that resolves when the SVG is colored.
 */
const colorSVG = async (chartData, plotConfig, datasetId, tileId, svgScoringMethod="gene") => {
    // I found adding the mid color for the colorblind mode  skews the whole scheme towards the high color
    const colorblindMode = CURRENT_USER.colorblind_mode;
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : (plotConfig?.low_color || '#e7d1d5');
    const midColor = colorblindMode ? null : (plotConfig?.mid_color || null);
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : (plotConfig?.high_color || '#401362');

    // Fill in tissue classes with the expression colors
    const {data: expression} = chartData;
    const tissues = Object.keys(expression);   // dataframe

    const score = chartData.scores[svgScoringMethod];

    // for those fields which have no reading, a specific value is sometimes put in instead
    // These are colored a neutral color
    const NA_FIELD_PLACEHOLDER = -0.012345679104328156;
    const NA_FIELD_COLOR = '#808080';

    // Load SVG file and set up the window
    const cardImage = document.querySelector(`#tile-${tileId} .card-image`);

    // create a legend div
    const legendDiv = document.createElement('div');
    legendDiv.classList.add('legend');
    legendDiv.style.zIndex = 1;
    legendDiv.style.height = "40px";    // match the viewbox height of child
    cardImage.append(legendDiv);

    // create a svg div (CSS for margin-top will now work nicely)
    const svgDiv = document.createElement('div');
    svgDiv.classList.add('svg');
    // higher z-index so we can mouseover the svg
    svgDiv.style.zIndex = 2;
    svgDiv.style.height = "calc(100% - 40px)";
    cardImage.append(svgDiv);

    const snap = Snap(svgDiv);
    const svg_path = `datasets_uploaded/${datasetId}.svg`;

    await Snap.load(svg_path, async (path) => {
        await snap.append(path);
        const svg = snap.select("svg");

        svg.attr({
            width: "100%"
        });

        // TODO: Set viewbar just like the legend.
        // TODO: Set at bottom of card-image

        // Get all paths, circles, rects, and ellipses
        const paths = Snap.selectAll("path, circle, rect, ellipse");

        if (svgScoringMethod === 'gene' || svgScoringMethod === 'dataset') {
            const { min, max } = score;
            let color = null;

            // are we doing a three- or two-color gradient?
            if (midColor) {
                if (min >= 0) {
                    // All values greater than 0, do right side of three-color
                    color = new_d3
                        .scaleLinear()
                        .domain([min, max])
                        .range([midColor, highColor]);
                } else if (max <= 0) {
                    // All values under 0, do left side of three-color
                    color = new_d3
                        .scaleLinear()
                        .domain([min, max])
                        .range([lowColor, midColor]);
                } else {
                    // We have a good value range, do the three-color
                    color = new_d3
                        .scaleLinear()
                        .domain([min, 0, max])
                        .range([lowColor, midColor, highColor]);
                }
            } else {
                color = new_d3
                    .scaleLinear()
                    .domain([min, max])
                    .range([lowColor, highColor]);
            }


            // NOTE: This must use the SnapSVG API Set.forEach function to iterate
            paths.forEach(path => {
                const tissue_classes = path.node.className.baseVal.split(' ');
                tissue_classes.forEach(tissue => {
                    if (!tissues.includes(tissue)) {
                        return;
                    }

                    if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                        path.attr('fill', NA_FIELD_COLOR);
                    } else {
                        path.attr('fill', color(expression[tissue]));
                    }

                    // log-transfom the expression score
                    const math = "raw";
                    let score;
                    // Apply math transformation to expression score
                    if (math == 'log2') {
                        score = new_d3.format('.2f')(Math.log2(expression[tissue]));
                    } else if (math == 'log10') {
                        score = new_d3.format('.2f')(Math.log10(expression[tissue]));
                    } else {
                        //math == 'raw'
                        score = new_d3.format('.2f')(expression[tissue]);
                    }

                    // Place tissue in score in a nice compact tooltip
                    const tooltipText = `${tissue}: ${score}`;

                    // Add mouseover and mouseout events to create and destroy the tooltip
                    path.mouseover(() => {
                        // get position of path node relative to page
                        const yOffset = path.node.getBoundingClientRect().top - svgDiv.getBoundingClientRect().top;

                        const tooltip = document.createElement('div');
                        tooltip.classList.add('tooltip');
                        tooltip.textContent = tooltipText;
                        tooltip.style.position = 'absolute';
                        tooltip.style.top = `${yOffset}px`;
                        tooltip.style.left = `${0}px`;
                        tooltip.style.backgroundColor = 'white';
                        tooltip.style.color = 'black';
                        tooltip.style.padding = '5px';
                        tooltip.style.border = '1px solid black';
                        tooltip.style.zIndex = 3;
                        svgDiv.appendChild(tooltip);

                    });
                    path.mouseout(() => {
                        svgDiv.querySelector('.tooltip').remove();
                    });

                });
            });
            return;
            // Draw the tooltip
        } else if (svgScoringMethod === 'tissue') {
            // tissues scoring
            const tissues = Object.keys(score);

            const color = {};

            if (midColor) {
                tissues.forEach(tissue => {
                    let {
                        min,
                        max
                    } = score[tissue];

                    if (min >= 0) {
                        color[tissue] = new_d3
                            .scaleLinear()
                            .domain([min, max])
                            .range([midColor, highColor]);
                    } else if (max <= 0) {
                        color[tissue] = new_d3
                            .scaleLinear()
                            .domain([min, max])
                            .range([lowColor, midColor]);
                    } else {
                        color[tissue] = new_d3
                            .scaleLinear()
                            .domain([min, 0, max])
                            .range([lowColor, midColor, highColor]);
                    }
                });
            } else {
                tissues.forEach(tissue => {
                    let {
                        min,
                        max
                    } = score[tissue];

                    color[tissue] = new_d3
                        .scaleLinear()
                        .domain([min, max])
                        .range([lowColor, highColor]);
                });
            }

            paths.forEach(path => {
                const tissue_classes = path.node.className.baseVal.split(' ');
                tissue_classes.forEach(tissue => {
                    if (!(tissue && color[tissue])) {
                        return;
                    }
                    const color_scale = color[tissue];
                    path.attr('fill', color_scale(expression[tissue]));

                    // log-transfom the expression score
                    const math = "raw";
                    let score;
                    // Apply math transformation to expression score
                    if (math == 'log2') {
                        score = new_d3.format('.2f')(Math.log2(expression[tissue]));
                    } else if (math == 'log10') {
                        score = new_d3.format('.2f')(Math.log10(expression[tissue]));
                    } else {
                        //math == 'raw'
                        score = new_d3.format('.2f')(expression[tissue]);
                    }

                    // Place tissue in score in a nice compact tooltip
                    const tooltipText = `${tissue}: ${score}`;

                    // Add mouseover and mouseout events to create and destroy the tooltip
                    path.mouseover(() => {
                        const yOffset = svgDiv.getBoundingClientRect().top - path.node.getBoundingClientRect().top;

                        const tooltip = document.createElement('div');
                        tooltip.classList.add('tooltip');
                        tooltip.textContent = tooltipText;
                        tooltip.style.position = 'absolute';
                        tooltip.style.top = `${yOffset}px`;
                        tooltip.style.left = `${0}px`;
                        tooltip.style.backgroundColor = 'white';
                        tooltip.style.color = 'black';
                        tooltip.style.padding = '5px';
                        tooltip.style.border = '1px solid black';
                        tooltip.style.zIndex = 3;
                        svgDiv.appendChild(tooltip);

                    });
                    path.mouseout(() => {
                        svgDiv.querySelector('.tooltip').remove();
                    });

                });
            });
        } else {
            throw new Error(`Invalid svgScoringMethod ${svgScoringMethod}.`);
        }

    });

    drawSVGLegend(plotConfig, tileId, score);

}

/**
 * Draws a legend for a SVG image
 *
 * @param {Object} plotConfig - The configuration for the plot.
 * @param {string} tileId - The ID of the tile.
 * @param {Object} score - The score object containing the minimum and maximum values.
 */
const drawSVGLegend = (plotConfig, tileId, score) => {
    const colorblindMode = CURRENT_USER.colorblind_mode;
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : plotConfig["low_color"];
    const midColor = colorblindMode ? null : plotConfig["mid_color"];
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : plotConfig["high_color"];

    const card = document.querySelector(`#tile-${tileId}.card`);
    const node = document.querySelector(`#tile-${tileId} .legend`);
    // Create our legend svg
    const legend = new_d3.select(node)  // returns document.documentElement
        .append('svg')
        .style('position', 'absolute')
        .style('width', '100%')
        .style("height", "40px")    // Without a fixed heigh, the box is too tall and prevents mouseover of the svg image
        .attr('viewbox', `0 0 ${node.getBoundingClientRect().width} 40`)
        .attr('class', 'svg-gradient-container');
    const defs = legend.append('defs');
    // Define our gradient shape
    const linearGradient = defs
        .append('linearGradient')
        .attr('id', `tile-${tileId}-linear-gradient`)
        .attr('x1', '0%')
        .attr('y1', '0%')
        .attr('x2', '100%')
        .attr('y2', '0%');

    const { min, max } = score;

    // Create the gradient points for either three- or two-color gradients
    if (midColor) {
        // Even if a midpoint is called for, it doesn't make sense if the values are
        //  all less than or all greater than 0
        if (min >= 0) {
            linearGradient
                .append('stop')
                .attr('offset', '0%')
                .attr('stop-color', midColor);
                linearGradient
                .append('stop')
                .attr('offset', '100%')
                .attr('stop-color', highColor);
        } else if (max <= 0) {
            linearGradient
                .append('stop')
                .attr('offset', '0%')
                .attr('stop-color', lowColor);
            linearGradient
                .append('stop')
                .attr('offset', '100%')
                .attr('stop-color', midColor);
        } else {
            // This means we've got a good distribution of min under 0 and max above
            //  it, so we can do a proper three-color range
            // midpoint offset calculation, so the mid color is at 0
            //var mid_offset = (1 - (min / max - min))*100;
            const midOffset = (Math.abs(min) / (max + Math.abs(min))) * 100;

            linearGradient
                .append('stop')
                .attr('offset', '0%')
                .attr('stop-color', lowColor);
            linearGradient
                .append('stop')
                .attr('offset', `${midOffset}%`)
                .attr('stop-color', midColor);
            linearGradient
                .append('stop')
                .attr('offset', '100%')
                .attr('stop-color', highColor);
        }
    } else {
        linearGradient
            .append('stop')
            .attr('offset', '0%')
            .attr('stop-color', lowColor);
        linearGradient
            .append('stop')
            .attr('offset', '100%')
            .attr('stop-color', highColor);
    }

    const width = node.getBoundingClientRect().width;

    // Draw the rectangle using the linear gradient
    legend
        .append('rect')
        .attr('width', "50%")
        .attr('y', 10)
        .attr('x', "25%")
        .attr('height', 10) // quarter of viewport height
        .style(
            'fill',
            `url(#tile-${tileId}-linear-gradient)`
        );

    const xScale = new_d3
        .scaleLinear()
        .domain([min, max])
        .range([0, width / 2]);

    const xAxis = new_d3
        .axisBottom()
        .ticks(3)
        .scale(xScale)

    legend
        .append('g')
        .attr('class', 'axis')
        .attr('transform', `translate(${width / 4}, 20)`)   // start quarter from left, and 10 px below rectangle
        .attr("stroke", "black")
        .call(xAxis);

    // Ensure axis is responsive
    window.addEventListener('resize', () => {
        legend.attr('viewbox', `0 0 ${card.getBoundingClientRect().width} 40`);

        // purge old axis
        legend.select('.axis').remove();

        // redraw axis
        // TODO: this is hacky, but it works
        const xScale = new_d3
            .scaleLinear()
            .domain([min, max])
            .range([0, card.getBoundingClientRect().width / 2]);

        const xAxis = new_d3
            .axisBottom()
            .ticks(3)
            .scale(xScale)

        legend
            .append('g')
            .attr('class', 'axis')
            .attr('transform', `translate(${card.getBoundingClientRect().width / 4}, 20)`)   // start quarter from left, and 10 px below rectangle
            .attr("stroke", "black")
            .call(xAxis);
    });
}

/**
 * Calculates the start position for the epiviz navigation based on the given data and extend range ratio.
 * @param {Object} data - The data object containing the start and end positions.
 * @param {number} extendRangeRatio - The ratio by which to extend the range.
 * @returns {number} - The calculated start position.
 */
const epivizNavStart = (data, extendRangeRatio) => {
    return data.start - Math.round((data.end - data.start) * extendRangeRatio);
}

/**
 * Calculates the end position of a navigation based on the given data and extend range ratio.
 * @param {Object} data - The data object containing the start and end positions.
 * @param {number} extendRangeRatio - The ratio by which to extend the range.
 * @returns {number} - The calculated end position.
 */
const epivizNavEnd = (data, extendRangeRatio) => {
    return data.end + Math.round((data.end - data.start) * extendRangeRatio);
}

/**
 * Retrieves updates and additions to the plot from the plot_display_config JS object.
 *
 * @param {Object[]} plotConfObj - The plot configuration object.
 * @param {string} plotType - The type of plot.
 * @param {string} category - The category of updates to retrieve.
 * @returns {Object} - The updates and additions to the plot.
 */
const getPlotlyDisplayUpdates = (plotConfObj, plotType, category) => {
    let updates = {};
    for (const idx in plotConfObj) {
        const conf = plotConfObj[idx];
        // Get config (data and/or layout info) for the plot type chosen, if it exists
        if (conf.plot_type == "all" || conf.plot_type == plotType) {
            const update = category in conf ? conf[category] : {};
            updates = {...updates, ...update};    // Merge updates
        }
    }
    return updates;
}

/**
 * Creates a card message for a tile.
 *
 * @param {string} tileId - The ID of the dataset tile.
 * @param {string} level - The level of the message to apply to Bulma CSS class (e.g., "info", "warning", "danger").
 * @param {string} message - The message to be displayed in the card.
 */
const createCardMessage = (tileId, level, message) => {
    const cardContent = document.querySelector(`#tile-${tileId} .card-image`);
    cardContent.replaceChildren();
    const messageElt = document.createElement("p");
    const textLevel = `has-text-${level}-dark`;
    const bgLevel = `has-background-${level}-light`;
    messageElt.classList.add(textLevel, bgLevel, "p-2", "m-2", "has-text-weight-bold");
    messageElt.textContent = message;
    cardContent.append(messageElt);
}
