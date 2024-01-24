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

    constructor(layoutShareId, selector ) {
        this.layoutShareId = layoutShareId;
        this.layout = [];   // this.getLayout();

        this.maxCols = 12 // highest number of columns in a row
        this.maxRows = 3 // highest number of rows in a column

        this.tiles = [];
        this.tilegrid = [] // this.generateTileGrid();
        this.selector = selector;

        //this.applyTileGrid();

    }

    /**
     * Adds all displays to the tile grid.
     * @returns {Promise<void>} A promise that resolves when all displays are added.
     */
    async addAllDisplays() {

        for (const dataset of this.layout) {
            dataset.userDisplays = [];
            dataset.ownerDisplays = [];

            const {user: userDisplays, owner: ownerDisplays} = await apiCallsMixin.fetchDatasetDisplays(dataset.id);
            dataset.userDisplays = userDisplays;
            dataset.ownerDisplays = ownerDisplays;
        }
    }

    /**
     * Adds default displays to all tiles in the tile grid.
     * @returns {Promise<void>} A promise that resolves when all default displays have been added.
     */
    async addDefaultDisplays() {
        // Each tile has "single" or "multi" type stored, so we can use that to determine the correct default display
        await Promise.allSettled(this.tiles.map( async tile => await tile.addDefaultDisplay()));
    }

    /**
     * Applies the tile grid to the specified selector element.
     */
    applyTileGrid() {
        const tilegrid = this.tilegrid;
        const selector = this.selector;

        // Clear selector element
        document.querySelector(selector).replaceChildren();

        // All child tiles fit into a single parent in a vertical sense.
        // If "is-vertical" is present, all children tiles will be stacked vertically
        // If "tile" class does not have "is-vertical", then children tiles will be stacked horizontally
        // A fully horizontal row of tiles will have a "is-ancestor" class

        for (const row of tilegrid) {
            const tilegridHTML = document.createElement('div');
            tilegridHTML.classList.add('tile', 'is-ancestor');
            for (const col of row) {
                const tile = col.tile;
                const tileChildHTML = tile.html;

                const tileParentHTML = document.createElement('div');
                tileParentHTML.classList.add('tile', 'is-parent', 'p-1', `is-${tile.width}`);
                tileParentHTML.append(tileChildHTML);

                tilegridHTML.append(tileParentHTML);
            }
            document.querySelector(selector).append(tilegridHTML);
        }

        // Make all card-header titles the same height
        //const cardHeaderTitles = document.querySelectorAll(`${selector} .card-header`);
        //const cardHeaderTitleHeight = Math.max(...Array.from(cardHeaderTitles).map(e => e.offsetHeight));
        //cardHeaderTitles.forEach(e => e.style.height = `${cardHeaderTitleHeight}px`);

    }

    // NOTE: This may change if data is returned previously and can be loaded
    /**
     * Retrieves the layout information from the server.
     * @returns {Promise<Array>} A promise that resolves to an array of datasets.
     * @throws {Error} If an error occurs during the API call.
     */
    async getLayout() {
        try {
            const data = await apiCallsMixin.fetchDatasetListInfo({layout_share_id: this.layoutShareId})
            return data['datasets'];
        } catch (error) {
            throw error;
        }
    }


    /**
     * Generates the tile grid based on the layout and dataset information.
     *
     * @param {boolean} is_multigene - Indicates whether the grid is for a multigene view.
     */
    generateTileGrid(is_multigene = false) {

        const tiles = [];
        const tilegrid = [];

        for (const dataset of this.layout) {
            const datasetTile = new DatasetTile(dataset, is_multigene);
            tiles.push(datasetTile);
        }

        // sort by grid position
        tiles.sort((a, b) => a.dataset.grid_position - b.dataset.grid_position);

        this.tiles = tiles;

        for (const datasetTile of tiles) {
            if (datasetTile.used) {
                continue;
            }

            const width = datasetTile.tile.width;
            const height = datasetTile.tile.height;

            if (width === this.maxCols) {
                // tile spans the entire row
                const tileRow = [];
                tileRow.push(datasetTile);
                tilegrid.push(tileRow);
                datasetTile.used = true;
                continue;
            }

            // tile does not span the entire row
            const tileRow = [];
            datasetTile.used = true;
            const usedTiles = [datasetTile];

            let remainingWidth = this.maxCols - width;

            // find tiles that fit into the remaining width
            while (remainingWidth > 0) {
                const tile = tiles.find((t) => !t.used && t.tile.width <= remainingWidth);
                if (!tile) {
                    break;
                }

                tile.used = true;
                usedTiles.push(tile);
                remainingWidth -= tile.tile.width;
            }

            tileRow.push(...usedTiles);
            tilegrid.push(tileRow);

            // check if all tiles are the same height
            const allSameHeight = usedTiles.every((t) => t.tile.height === height);
            if (allSameHeight) {
                // all tiles are the same height
                //tileRow.push(...usedTiles);
                //tilegrid.push(tileRow);
            }
        }

        this.tilegrid = tilegrid;

        // TODO: Create a subgrid for variable heights
    }

    /**
     * Renders the displays for the given gene symbols.
     *
     * @param {string|string[]} geneSymbols - The gene symbol or an array of gene symbols.
     * @param {boolean} [isMultigene=false] - Indicates whether multiple gene symbols are provided.
     * @param {string} svgScoringMethod - The SVG scoring method.
     * @throws {Error} If geneSymbols is not provided.
     * @returns {Promise<void>} A promise that resolves when all displays are rendered.
     */
    async renderDisplays(geneSymbols, isMultigene = false, svgScoringMethod) {
        if (!geneSymbols) {
            throw new Error("Gene symbol or symbols are required to render displays.");
        }

        let geneSymbolInput = geneSymbols;
        if (!isMultigene) {
            geneSymbolInput = Array.isArray(geneSymbols) ? geneSymbols[0] : geneSymbols;
        }

        // Sometimes fails to render due to OOM errors, so we want to try each tile individually
        for (const tile of this.tiles) {
            await tile.renderDisplay(geneSymbolInput, null, svgScoringMethod);
        }

        // await Promise.allSettled(this.tiles.map( async tile => await tile.renderDisplay(geneSymbolInput, null, svgScoringMethod)));
    }
};

class DatasetTile {
    constructor(dataset, isMulti = true) {
        this.dataset = dataset;
        this.type = isMulti ? 'multi' : 'single';
        this.typeInt = isMulti ? 1 : 0;

        this.performingProjection = false;  // Indicates whether a projection is currently being performed

        this.controller = new AbortController(); // Create new controller for new set of frames


        this.tile = this.generateTile();
        this.tile.html = this.generateTileHTML();
    }

    /**
     * Adds a default display to the tile grid.
     * @returns {Promise<void>} A promise that resolves when the default display is added.
     */
    async addDefaultDisplay() {

        const dataset = this.dataset;
        this.defaultDisplayId = null;
        try {
            const {default_display_id: defaultDisplayId} = await apiCallsMixin.fetchDefaultDisplay(dataset.id, this.typeInt);
            this.defaultDisplayId = defaultDisplayId;
        } catch (error) {
            //pass
        }
    }

    // TODO: Refactor since both of these functions are the same except for the using "grid_width" vs "mg_grid_width"
    // TODO: Also add "grid_height" and "mg_grid_height" to the dataset object
    /**
     * Generates a tile object based on the current dataset and tile type.
     * @returns {Object} The generated tile object.
     */
    generateTile() {
        const tile = {};
        const dataset = this.dataset;
        tile.width = this.type === "single" ? dataset.grid_width : dataset.mg_grid_width;
        tile.height = 1; // dataset.grid_height || dataset.mg_grid_height;  // Heights are not bound by the 12-spaced grid, so just use 1, 2, 3, etc.
        tile.tile_id = `${dataset.id}_${dataset.grid_position}_${this.type}`
        tile.used = false;
        tile.title = dataset.title;
        return tile;
    }

    /**
     * Generates the HTML representation of a tile.
     * @returns {HTMLElement} The HTML element representing the tile.
     *
     * Note: Tile html template comes from /include/tile-grid/tile.html
     */
    generateTileHTML() {
        const tile = this.tile;

        const template = document.getElementById('tmpl-tile-grid-tile');
        const tileHTML = template.content.cloneNode(true);

        // Set tile id & title
        const tileElement = tileHTML.querySelector('.tile');
        tileElement.id = `tile_${tile.tile_id}`;

        const tileTitle = tileHTML.querySelector('.card-header-title');
        tileTitle.textContent = tile.title;

        this.addDropdownInformation(tileElement, tileElement.id, this.dataset);

        return tileHTML;
    }

    /**
     * Renders a modal to choose a display from the user or owner display lists.
     *
     * @returns {void}
     */
    renderChooseDisplayModal() {

        // Remove any existing modals
        const existingModals = document.querySelectorAll('.js-choose-display-modal');
        for (const modal of existingModals) {
            modal.remove();
        }

        const modalTemplate = document.getElementById('tmpl-tile-grid-choose-display-modal');
        const modalHTML = modalTemplate.content.cloneNode(true);

        const modalDiv = modalHTML.querySelector('.modal');
        modalDiv.id = `choose-display-modal_${this.tile.tile_id}`;

        const modalContent = modalHTML.querySelector('.modal-content');

        // Add dataset title
        const datasetTitle = modalContent.querySelector("h5");
        datasetTitle.replaceChildren();
        datasetTitle.textContent = this.dataset.title;

        // Add user and owner displays
        const userDisplaysElt = modalContent.querySelector(".js-modal-user-displays");
        userDisplaysElt.replaceChildren();
        const ownerDisplaysElt = modalContent.querySelector(".js-modal-owner-displays");
        ownerDisplaysElt.replaceChildren();

        // Get all user and owner displays for a single-gene or multi-gene view
        const filterKey = this.type === "single" ? "gene_symbol" : "gene_symbols";

        // Find all the display config in the user or owner display lists
        const userDisplays = this.dataset.userDisplays.filter((d) => d.plotly_config.hasOwnProperty(filterKey));
        const ownerDisplays = this.dataset.ownerDisplays.filter((d) => d.plotly_config.hasOwnProperty(filterKey));

        // Append epiviz displays to user and owner displays
        if (this.type === "single") {
            const userEpivizDisplays = this.dataset.userDisplays.filter((d) => d.plot_type === "epiviz");
            const ownerEpivizDisplays = this.dataset.ownerDisplays.filter((d) => d.plot_type === "epiviz");
            userDisplays.push(...userEpivizDisplays);
            ownerDisplays.push(...ownerEpivizDisplays);
        }

        // Add titles to each section if there are displays
        if (userDisplays.length) {
            const userTitle = document.createElement("p");
            userTitle.classList.add("has-text-weight-bold", "is-underlined", "column", "is-full");
            userTitle.textContent = "Your Displays";
            userDisplaysElt.append(userTitle);

        }

        if (ownerDisplays.length) {
            const ownerTitle = document.createElement("p");
            ownerTitle.classList.add("has-text-weight-bold", "is-underlined", "column", "is-full");
            ownerTitle.textContent = "Displays by Dataset Owner";
            ownerDisplaysElt.append(ownerTitle);
        }

        // Did it this way so we didn't have to pass async/await up the chain
        Promise.allSettled([
            this.renderChooseDisplayModalDisplays(userDisplays, userDisplaysElt),
            this.renderChooseDisplayModalDisplays(ownerDisplays, ownerDisplaysElt)
        ]).then(() => {
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
        // Add user displays
        for (const display of displays) {
            const displayTemplate = document.getElementById('tmpl-tile-grid-choose-display-modal-display');
            const displayHTML = displayTemplate.content.cloneNode(true);

            const displayElement = displayHTML.querySelector('.js-modal-display');
            displayElement.dataset.displayId = display.id;
            displayElement.dataset.datasetId = this.dataset.id;

            // Add display image
            let displayUrl = "";
            try {
                // TODO: SVGs are colorless
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
     * Adds dropdown information to a tile element.
     *
     * @param {HTMLElement} tileElement - The tile element to add dropdown information to.
     * @param {string} tileId - The ID of the tile.
     * @param {object} dataset - The dataset object containing information for the dropdown items.
     */
    addDropdownInformation(tileElement, tileId, dataset) {

        const dropdownMenu = tileElement.querySelector(`#${tileId} .dropdown-menu`);
        dropdownMenu.id = `dropdown-menu_${tileId}`;
        const dropdownContent = dropdownMenu.querySelector(".dropdown-content");
        const dropdownItems = dropdownContent.querySelectorAll('.dropdown-item');

        const datasetId = dataset.id;
        const pubmedId = dataset.pubmed_id;
        const geoId = dataset.geo_id;
        const hasTarball = dataset.has_tarball;
        const hasH5ad = dataset.has_h5ad;
        const links = dataset.links;

        for (const item of dropdownItems) {
            switch (item.dataset.tool) {
                case "display":
                    // Add event listener to dropdown item
                    item.addEventListener("click", (event) => {
                        // Create and open modal for all user and owner displays
                        this.renderChooseDisplayModal();    // TODO: Need to add after displays are retrieved
                        const modalElt = document.getElementById(`choose-display-modal_${this.tile.tile_id}`);
                        openModal(modalElt);
                    });
                    break;
                case "expand":
                    // Zoom panel to take up all of "#result-panel-grid"
                    break;
                case "info":
                    // Modal for dataset information
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
                    break;
                case "single-cell":
                    // Redirect to single-cell analysis workbench
                    if (hasH5ad) {
                        const url = `./analyze_dataset.html?dataset_id=${datasetId}`;
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
                    if (hasTarball) {
                        const url = `./cgi/download_source_file.cgi?type=tarball&dataset_id=${datasetId}`;
                        item.href = url;
                    } else {
                        item.classList.add("is-hidden");
                    }
                    break;
                case "download-h5ad":
                    // Download h5ad file
                    if (hasH5ad) {
                        const url = `./cgi/download_source_file.cgi?type=h5ad&dataset_id=${datasetId}`;
                        item.href = url;
                    } else {
                        item.classList.add("is-hidden");
                    }
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

        // Add event listener to dropdown trigger
        tileElement.querySelector("button.dropdown-trigger").addEventListener("click", (event) => {
            const item = event.currentTarget;
            if (item.classList.contains('dropdown-trigger')) {
                // close other dropdowns
                for (const dropdown of document.querySelectorAll('#result-panel-grid .dropdown')) {
                    if (item !== dropdown.querySelector('.dropdown-trigger')) {
                        dropdown.classList.remove('is-active');
                    }
                };

                item.closest(".dropdown").classList.toggle('is-active');
            }
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
     * Renders the display for a given gene symbol.
     * @param {string} geneSymbol - The gene symbol to render the display for.
     * @param {string|null} displayId - The ID of the display to render. If null, the default display ID will be used.
     * @param {string} svgScoringMethod - The SVG scoring method to use.
     * @throws {Error} If geneSymbol is not provided.
     * @returns {Promise<void>} A promise that resolves when the display is rendered.
     */
    async renderDisplay(geneSymbol, displayId=null, svgScoringMethod="gene") {
        if (!geneSymbol) {
            throw new Error("Gene symbol or symbols are required to render this display.");
        }

        if (displayId === null) {
            displayId = this.defaultDisplayId;
        };

        this.geneSymbol = geneSymbol;
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

        // Try epiviz display if no plotly display was found
        if (this.type === "single") {
            if (!userDisplay) userDisplay = this.dataset.userDisplays.find((d) => d.id === displayId && d.plot_type === "epiviz");

            if (!ownerDisplay) ownerDisplay = this.dataset.ownerDisplays.find((d) => d.id === displayId && d.plot_type === "epiviz");
        }

        // add console warning if default display id was not found in the user or owner display lists
        if (!userDisplay && !ownerDisplay) {
            console.warn(`Selected display config for dataset ${this.dataset.title} was not found. Will show first available.`);

            // last chance... if still no display config (i.e default display was not found), then use the first display config
            if (!userDisplay) userDisplay = this.dataset.userDisplays.find((d) => d.plotly_config.hasOwnProperty(filterKey));
            if (!ownerDisplay) ownerDisplay = this.dataset.ownerDisplays.find((d) => d.plotly_config.hasOwnProperty(filterKey));
        }

        // if the display config was not found, then do not render
        if (!userDisplay && !ownerDisplay) {
            console.warn(`Display config for dataset ${this.dataset.title} was not found.`)
            // Let the user know that the display config was not found
            const cardContent = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
            cardContent.replaceChildren();
            const warningMessage = document.createElement("p");
            warningMessage.classList.add("has-text-warning-dark", "has-background-warning-light", "p-2", "m-2", "has-text-weight-bold");
            // Add 200 px height and center vertically
            warningMessage.style.height = "200px";
            warningMessage.style.display = "flex";
            warningMessage.style.alignItems = "center";
            warningMessage.style.justifyContent = "center";
            warningMessage.textContent = `This dataset has no viewable curations for this view. Create a new curation in the ${this.type === "single" ? "Single-gene" : "Multi-gene"} Curator to view this dataset.`;
            cardContent.append(warningMessage);
            return;
        }

        // if the display config was found, then render
        const display = userDisplay || ownerDisplay;

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
            display.plotly_config.gene_symbols = geneSymbol;
        } else {
            display.plotly_config.gene_symbol = geneSymbol;
        }
        const cardContent = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        //cardContent.classList.add("loader");

        try {
            if (plotlyPlots.includes(display.plot_type)) {
                await this.renderPlotlyDisplay(display, otherOpts);
            } else if (scanpyPlots.includes(display.plot_type)) {
                await this.renderScanpyDisplay(display, otherOpts);
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
            console.error(error);
            // Fill in card-image with error message
            cardContent.replaceChildren();

            const template = document.getElementById('tmpl-tile-grid-error');
            const errorHTML = template.content.cloneNode(true);
            const errorElement = errorHTML.querySelector('p');
            errorElement.textContent = error.message;
            cardContent.append(errorHTML);
        } finally {
           // cardContent.classList.remove("loader");
        }

    }

    epivizNavStart(data, extendRangeRatio) {
        return data.start - Math.round((data.end - data.start) * extendRangeRatio);
    }

    epivizNavEnd(data, extendRangeRatio) {
        return data.end + Math.round((data.end - data.start) * extendRangeRatio);
    }

    /**
     * Renders the Epiviz display on the tile grid.
     * @param {Object} display - The display object containing dataset_id and plotly_config.
     * @returns {Promise<void>} - A promise that resolves when the rendering is complete.
     */
    async renderEpivizDisplay(display, otherOpts) {
        const datasetId = display.dataset_id;
        const {gene_symbol: geneSymbol} = display.plotly_config;

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
        const epiviznav = document.querySelector(`#epiviznav_${this.tile.tile_id}`);
        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
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
        const nStart = this.epivizNavStart(data, extendRangeRatio);
        const nEnd = this.epivizNavEnd(data, extendRangeRatio);

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
        epivizContainer.id = `epiviz_${this.tile.tile_id}`;
        const epivizDataSource = epivizHTML.querySelector('epiviz-data-source');
        epivizDataSource.id = `${this.tile.tile_id}epivizds`;
        epivizDataSource.setAttribute("provider-url", plotConfig.dataserver);

        const epivizNavigation = epivizHTML.querySelector('epiviz-navigation');
        epivizNavigation.id = `${this.tile.tile_id}_epiviznav`;
        // the chr, start and end should come from query - map gene to genomic position.
        epivizNavigation.setAttribute("chr", data.chr);
        epivizNavigation.setAttribute("start", this.epivizNavStart(data, extendRangeRatio));
        epivizNavigation.setAttribute("end", this.epivizNavEnd(data, extendRangeRatio));
        epivizNavigation.setAttribute("viewer", `/epiviz.html?dataset_id=${this.dataset.id}&chr=${data.chr}&start=${data.start}&end=${data.end}`);
        epivizNavigation.innerHTML(this.renderEpivizTracks(plotConfig));
        return epivizHTML;
    }

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

        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const plotlyPreview = document.createElement("div");
        plotlyPreview.classList.add("container");
        plotlyPreview.id = `tile_${this.tile.tile_id}_plotly_preview`;
        plotContainer.append(plotlyPreview);
        Plotly.purge(plotlyPreview.id); // clear old Plotly plots

        if (!plotJson) {
            console.warn(`Could not retrieve plot information for dataset display ${display.id}. Cannot make plot.`);
            return;
        }

        if (plotType === 'heatmap') {
            // These modify the plotJson object in place
            // TODO: Adjust these functions
            //adjustExpressionColorbar(plotJson.data);
            //adjustClusterColorbars(plotJson.data);
        }

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const expressionDisplayConf = postPlotlyConfig.expression;
        const custonConfig = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "config");
        Plotly.newPlot(plotlyPreview.id , plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "layout")
        Plotly.relayout(plotlyPreview.id , custonLayout)

        const legendTitle = document.getElementById("legend_title_container");
        if (legendTitle) {
            legendTitle.classList.remove("is-hidden");
            if (plotType === "dotplot") {
                legendTitle.classList.add("is-hidden");
            }
        }

    }

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

        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        // NOTE: Plot initially is created to a default width but is responsive.
        // Noticed container within our "column" will make full-width go beyond the screen
        const plotlyPreview = document.createElement("div");
        plotlyPreview.classList.add("container");
        plotlyPreview.id = `tile_${this.tile.tile_id}_plotly_preview`;
        plotContainer.append(plotlyPreview);
        Plotly.purge(plotlyPreview.id); // clear old Plotly plots

        if (!plotJson) {
            console.warn(`Could not retrieve plot information for dataset display ${display.id}. Cannot make plot.`);
            return;
        }
        // Update plot with custom plot config stuff stored in plot_display_config.js
        const expressionDisplayConf = postPlotlyConfig.expression;
        const custonConfig = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "config");
        Plotly.newPlot(plotlyPreview.id, plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "layout")
        Plotly.relayout(plotlyPreview.id, custonLayout)
    }

    async renderScanpyDisplay(display, otherOpts) {

        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        const data = await apiCallsMixin.fetchTsneImage(datasetId, analysisObj, plotType, plotConfig, otherOpts);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {image} = data;

        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        const tsnePreview = document.createElement("img");
        tsnePreview.classList.add("image");
        tsnePreview.id = `tile_${this.tile.tile_id}_tsne_preview`;;
        plotContainer.append(tsnePreview);

        if (image) {
            document.getElementById(tsnePreview.id ).setAttribute("src", `data:image/png;base64,${image}`);
        } else {
            console.warn(`Could not retrieve plot image for dataset display ${display.id}. Cannot make plot.`);
            return;
        }
    }

    async renderSVG(display, svgScoringMethod="gene", otherOpts) {
        const datasetId = display.dataset_id;
        const plotConfig = display.plotly_config;
        const {gene_symbol: geneSymbol} = plotConfig;

        const data = await apiCallsMixin.fetchSvgData(datasetId, geneSymbol, otherOpts)
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        if (!plotContainer) return; // tile was removed before data was returned
        plotContainer.replaceChildren();    // erase plot

        colorSVG(data, plotConfig.colors, this, svgScoringMethod);

    }

    resetAbortController() {
        if (this.controller && !this.performingProjection) {
            this.controller.abort(); // Cancel any previous axios requests (such as drawing plots for a previous dataset)
        }
        this.controller = new AbortController(); // Create new controller for new set of frames

    }
}

const colorSVG = async (chartData, plotConfig, datasetTile, svgScoringMethod="gene") => {
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
    const cardImage = document.querySelector(`#tile_${datasetTile.tile.tile_id} .card-image`);

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
    cardImage.append(svgDiv);

    const snap = Snap(svgDiv);
    const svg_path = `datasets_uploaded/${datasetTile.dataset.id}.svg`;

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

    drawLegend(plotConfig, datasetTile, score)

}

const drawLegend = (plotConfig, datasetTile, score) => {
    const colorblindMode = CURRENT_USER.colorblind_mode;
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : plotConfig["low_color"];
    const midColor = colorblindMode ? null : plotConfig["mid_color"];
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : plotConfig["high_color"];

    const card = document.querySelector(`#tile_${datasetTile.tile.tile_id}.card`);
    const node = document.querySelector(`#tile_${datasetTile.tile.tile_id} .legend`);
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
        .attr('id', `tile_${datasetTile.tile.tile_id}-linear-gradient`)
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
            `url(#tile_${datasetTile.tile.tile_id}-linear-gradient)`
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