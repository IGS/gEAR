'use strict';

/* Given a passed-in layout_id, genereate a 2-dimensional tile-based grid object.
This uses Bulma CSS for stylings (https://bulma.io/documentation/layout/tiles/)
For the given layout, a single-gene grid and a multi-gene grid are generated.
*/

class TileGrid {

    constructor(layoutShareId, selector ) {
        this.layoutShareId = layoutShareId;
        this.layout = [];   // this.getLayout();

        this.maxCols = 12 // highest number of columns in a row
        this.maxRows = 12 // highest number of rows in a column

        this.singleGeneTiles = [];  // collection of DatasetTile objects
        this.multiGeneTiles = [];

        this.tilegrid = [] // this.generateTileGrid();
        this.selector = selector;

        //this.applyTileGrid();

    }

    async addDefaultDisplays() {

    }

    /**
     * Applies the tile grid layout to the specified element.
     *
     * @param {boolean} [isMulti=false] - Indicates whether the multi-tile grid should be used.
     */
    applyTileGrid(isMulti = false) {
        const tilegrid = isMulti ? this.tilegrid.multi : this.tilegrid.single;
        const selector = this.selector;

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
            console.error(error);
            throw error;
        }
    }

    /**
     * Generates a tile grid.
     * @returns {Object} The generated tile grid object.
     */
    generateTileGrid() {
        const tilegrid = {};

        tilegrid.single = this.generateSingleTileGrid();
        tilegrid.multi = this.generateMultiTileGrid();
        return tilegrid;
    }

    /**
     * Generates a multi-gene tile grid based on the layout.
     * @returns {Array<Array<DatasetTile>>} The generated tile grid.
     */
    generateMultiTileGrid() {
        const tilegrid = [];
        const tiles = [];

        for (const dataset of this.layout) {
            const datasetTile = new DatasetTile(dataset, true);
            tiles.push(datasetTile);
        }

        // sort by grid position
        tiles.sort((a, b) => a.dataset.grid_position - b.dataset.grid_position);

        this.multiGeneTiles = tiles;

        for (const datasetTile of tiles) {
            if (datasetTile.used) {
                continue;
            }

            const width = datasetTile.tile.width;
            const height = datasetTile.tile.height;

            if (width === 12) {
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

            let remainingWidth = 12 - width;

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

        return tilegrid;
    }

    /**
     * Generates a single-gene tile grid based on the layout.
     * @returns {Array<Array<DatasetTile>>} The generated tile grid.
     */
    generateSingleTileGrid() {
        const tilegrid = [];
        const tiles = [];

        for (const dataset of this.layout) {
            const datasetTile = new DatasetTile(dataset, false);
            tiles.push(datasetTile);
        }

        // sort by grid position
        tiles.sort((a, b) => a.dataset.grid_position - b.dataset.grid_position);

        this.singleGeneTiles = tiles;


        for (const datasetTile of tiles) {
            if (datasetTile.used) {
                continue;
            }

            const width = datasetTile.tile.width;
            const height = datasetTile.tile.height;

            if (width === 12) {
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

            let remainingWidth = 12 - width;

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

        return tilegrid;;

        // TODO: Create a subgrid for variable heights
    }
};

class DatasetTile {
    constructor(dataset, isMulti = true) {
        this.dataset = dataset;
        this.type = isMulti ? 'multi' : 'single';
        this.tile = this.generateTile();
        this.tile.html = this.generateTileHTML();
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
     */
    generateTileHTML() {
        const tile = this.tile;

        const tileHTML = document.createElement('div');
        tileHTML.classList.add('tile', 'is-child', 'card');
        tileHTML.id = `tile_${tile.tile_id}`;

        // card header (title + icon)
        const cardHeader = document.createElement('div');
        cardHeader.classList.add('card-header');
        const cardHeaderTitle = document.createElement('div');
        cardHeaderTitle.classList.add('card-header-title');
        cardHeaderTitle.textContent = tile.title;

        const cardHeaderIcon = document.createElement('div');
        cardHeaderIcon.classList.add('card-header-icon');

        const cardHeaderIconSpan = document.createElement('span');
        cardHeaderIconSpan.classList.add('icon');

        const cardHeaderIconSpanI = document.createElement('i');
        cardHeaderIconSpanI.classList.add('mdi', 'mdi-18px', 'mdi-dots-vertical');

        cardHeaderIconSpan.appendChild(cardHeaderIconSpanI);
        cardHeaderIcon.appendChild(cardHeaderIconSpan);
        cardHeader.appendChild(cardHeaderTitle);
        cardHeader.appendChild(cardHeaderIcon);
        tileHTML.appendChild(cardHeader);

        // card content
        const cardContent = document.createElement('div');
        cardContent.classList.add('card-content');
        cardContent.textContent = 'Card Graphic here';  // TODO: Add graphic here
        tileHTML.appendChild(cardContent);

        return tileHTML;
    }

}