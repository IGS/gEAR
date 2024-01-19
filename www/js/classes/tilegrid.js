'use strict';

/* Given a passed-in layout_id, genereate a 2-dimensional tile-based grid object.
This uses Bulma CSS for stylings (https://bulma.io/documentation/layout/tiles/)
For the given layout, a single-gene grid and a multi-gene grid are generated.
*/

const plotlyPlots = ["bar", "line", "scatter", "tsne/umap_dynamic", "violin"];  // "tsne_dynamic" is a legacy option
const scanpyPlots = ["pca_static", "tsne_static", "umap_static"];   // "tsne" is a legacy option


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

    async addAllDisplays() {

        for (const dataset of this.layout) {
            dataset.userDisplays = [];
            dataset.ownerDisplays = [];

            const {user: userDisplays, owner: ownerDisplays} = await apiCallsMixin.fetchDatasetDisplays(dataset.id);
            dataset.userDisplays = userDisplays;
            dataset.ownerDisplays = ownerDisplays;
        }
    }

    async addDefaultDisplays() {
        await Promise.allSettled(this.singleGeneTiles.map( async tile => await tile.addDefaultDisplay()));
        await Promise.allSettled(this.multiGeneTiles.map( async tile => await tile.addDefaultDisplay()));
    }

    /**
     * Applies the tile grid layout to the specified element.
     *
     * @param {boolean} [isMulti=false] - Indicates whether the multi-tile grid should be used.
     */
    applyTileGrid(isMulti = false) {
        const tilegrid = isMulti ? this.tilegrid.multi : this.tilegrid.single;
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

    async renderDisplays(geneSymbols, isMultigene = false) {
        if (!geneSymbols) {
            throw new Error("Gene symbol or symbols are required to render displays.");
        }

        if (isMultigene) {
            await Promise.allSettled(this.multiGeneTiles.map( async tile => await tile.renderDisplay(geneSymbols)));
        } else {
            const geneSymbol = Array.isArray(geneSymbols) ? geneSymbols[0] : geneSymbols;
            await Promise.allSettled(this.singleGeneTiles.map( async tile => await tile.renderDisplay(geneSymbol)));
        }


    }

};

class DatasetTile {
    constructor(dataset, isMulti = true) {
        this.dataset = dataset;
        this.type = isMulti ? 'multi' : 'single';
        this.typeInt = isMulti ? 1 : 0;

        this.tile = this.generateTile();
        this.tile.html = this.generateTileHTML();
    }

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
     */
    generateTileHTML() {
        const tile = this.tile;

        const tileHTML = document.createElement('div');
        tileHTML.classList.add('tile', 'is-child', 'card', 'has-background-white');
        tileHTML.id = `tile_${tile.tile_id}`;

        // card header (title + icon)
        const cardHeader = document.createElement('div');
        cardHeader.classList.add('card-header', 'has-background-primary');
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
        cardContent.classList.add('card-image');
        tileHTML.appendChild(cardContent);

        return tileHTML;
    }

    async renderDisplay(geneSymbol, displayId=null) {
        if (!geneSymbol) {
            throw new Error("Gene symbol or symbols are required to render this display.");
        }

        if (displayId === null) {
            displayId = this.defaultDisplayId;
        };

        const filterKey = this.type === "single" ? "gene_symbol" : "gene_symbols";

        // find the display config in the user or owner display lists
        const userDisplay = this.dataset.userDisplays.find((d) => d.id === displayId && d.plotly_config.hasOwnProperty(filterKey));
        const ownerDisplay = this.dataset.ownerDisplays.find((d) => d.id === displayId && d.plotly_config.hasOwnProperty(filterKey));

        // if the display config was not found, then do not render
        if (!userDisplay && !ownerDisplay) {
            console.warn(`Display config for dataset ${this.dataset.id} was not found.`)
            return;
        }

        // if the display config was found, then render
        const display = userDisplay || ownerDisplay;

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
        cardContent.classList.add("is-loading");

        try {
            if (plotlyPlots.includes(display.plot_type)) {
                await this.renderPlotlyDisplay(display);
            } else if (scanpyPlots.includes(display.plot_type)) {
                await this.renderScanpyDisplay(display);
            } else if (display.plot_type === "svg") {
                await this.renderSVG(display);
            } else if (this.type === "multi") {
                await this.renderMultiGeneDisplay(display);
            } else {
                throw new Error(`Display config for dataset ${this.dataset.id} has an invalid plot type ${display.plot_type}.`);
            }
        } catch (error) {
            console.error(error);
            // we want to ensure other plots load even if one fails
        } finally {
            cardContent.classList.remove("is-loading");
        }

    }

    async renderMultiGeneDisplay(display) {

        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        // Get data and set up the image area
        const data = await apiCallsMixin.fetchDashData(datasetId, analysisObj, plotType, plotConfig);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {plot_json: plotJson} = data;

        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
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
            setHeatmapHeightBasedOnGenes(plotJson.layout, plotConfig.gene_symbols);
        } else if (plotType === "mg_violin" && plotConfig.stacked_violin){
            adjustStackedViolinHeight(this.plotJson.layout);
        }

        // Update plot with custom plot config stuff stored in plot_display_config.js
        const expressionDisplayConf = postPlotlyConfig.expression;
        const custonConfig = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "config");
        Plotly.newPlot(plotlyPreview.id , plotJson.data, plotJson.layout, custonConfig);
        const custonLayout = getPlotlyDisplayUpdates(expressionDisplayConf, this.plotType, "layout")
        Plotly.relayout(plotlyPreview.id , custonLayout)

        document.getElementById("legend_title_container").classList.remove("is-hidden");
        if (plotType === "dotplot") {
            document.getElementById("legend_title_container").classList.add("is-hidden");
        }

    }

    async renderPlotlyDisplay(display) {
        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        // Get data and set up the image area
        const data = await apiCallsMixin.fetchPlotlyData(datasetId, analysisObj, plotType, plotConfig);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {plot_json: plotJson} = data;

        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
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

    async renderScanpyDisplay(display) {

        const datasetId = display.dataset_id;
        // Create analysis object if it exists.  Also supports legacy "analysis_id" string
        const analysisObj = display.plotly_config.analysis_id ? {id: display.plotly_config.analysis_id} : display.plotly_config.analysis || null;
        const plotType = display.plot_type;
        const plotConfig = display.plotly_config;

        const data = await apiCallsMixin.fetchTsneImage(datasetId, analysisObj, plotType, plotConfig);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const {image} = data;

        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
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

    async renderSVG(display) {
        const datasetId = display.dataset_id;
        const plotConfig = display.plotly_config;
        const {gene_symbol: geneSymbol} = plotConfig;


        const data = await apiCallsMixin.fetchSvgData(datasetId, geneSymbol)
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        plotContainer.replaceChildren();    // erase plot

        colorSVG(data, plotConfig.colors, this);

    }
}

const colorSVG = (chartData, plotConfig, datasetTile) => {
    // I found adding the mid color for the colorblind mode  skews the whole scheme towards the high color
    const colorblindMode = CURRENT_USER.colorblind_mode;
    const lowColor = colorblindMode ? 'rgb(254, 232, 56)' : plotConfig["low_color"];
    const midColor = colorblindMode ? null : plotConfig["mid_color"];
    const highColor = colorblindMode ? 'rgb(0, 34, 78)' : plotConfig["high_color"];

    // for those fields which have no reading, a specific value is sometimes put in instead
    // These are colored a neutral color
    const NA_FIELD_PLACEHOLDER = -0.012345679104328156;
    const NA_FIELD_COLOR = '#808080';

    //const scoreMethod = document.getElementById("scoring_method").value;
    const score = chartData.scores["gene"]
    const { min, max } = score;
    let color = null;
    // are we doing a three- or two-color gradient?
    if (midColor) {
        if (min >= 0) {
            // All values greater than 0, do right side of three-color
            color = d3
                .scaleLinear()
                .domain([min, max])
                .range([midColor, highColor]);
        } else if (max <= 0) {
            // All values under 0, do left side of three-color
            color = d3
                .scaleLinear()
                .domain([min, max])
                .range([lowColor, midColor]);
        } else {
            // We have a good value range, do the three-color
            color = d3
                .scaleLinear()
                .domain([min, 0, max])
                .range([lowColor, midColor, highColor]);
        }
    } else {
        color = d3
            .scaleLinear()
            .domain([min, max])
            .range([lowColor, highColor]);
    }


    // Load SVG file and set up the window
    const svg = document.querySelector(`#tile_${datasetTile.tile.tile_id} .card-image`);
    const snap = Snap(svg);
    const svg_path = `datasets_uploaded/${datasetTile.dataset.id}.svg`;
    Snap.load(svg_path, async (path) => {
        await snap.append(path)

        snap.select("svg").attr({
            width: "100%",
        });

        // Fill in tissue classes with the expression colors
        const {data: expression} = chartData;
        const tissues = Object.keys(chartData.data);   // dataframe
        const paths = Snap.selectAll("path, circle");

        // NOTE: This must use the SnapSVG API Set.forEach function to iterate
        paths.forEach(path => {
            const tissue = path.node.className.baseVal;
            if (tissues.includes(tissue)) {
                if (expression[tissue] == NA_FIELD_PLACEHOLDER) {
                    path.attr('fill', NA_FIELD_COLOR);
                } else {
                    path.attr('fill', color(expression[tissue]));
                }
            }
        });

        // TODO: Potentially replicate some of the features in display.js like log-transforms and tooltips
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