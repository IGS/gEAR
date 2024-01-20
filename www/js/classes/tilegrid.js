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
        cardHeaderTitle.classList.add('card-header-title', "py-0");
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

    async renderDisplay(geneSymbol, displayId=null, svgScoringMethod="gene") {
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
            // TODO: Add Epiviz
            if (plotlyPlots.includes(display.plot_type)) {
                await this.renderPlotlyDisplay(display);
            } else if (scanpyPlots.includes(display.plot_type)) {
                await this.renderScanpyDisplay(display);
            } else if (display.plot_type === "svg") {
                await this.renderSVG(display, svgScoringMethod);
            } else if (this.type === "multi") {
                await this.renderMultiGeneDisplay(display);
            } else {
                throw new Error(`Display config for dataset ${this.dataset.id} has an invalid plot type ${display.plot_type}.`);
            }
        } catch (error) {
            // we want to ensure other plots load even if one fails
            console.error(error);
            // Fill in card-image with error message
            cardContent.replaceChildren();
            const errorMessage = document.createElement("p");
            errorMessage.classList.add("has-text-danger-dark", "has-background-danger-light", "p-2", "m-2", "has-text-weight-bold");
            // Add 200 px height and center vertically
            errorMessage.style.height = "200px";
            errorMessage.style.display = "flex";
            errorMessage.style.alignItems = "center";
            errorMessage.style.justifyContent = "center";
            errorMessage.textContent = error.message;
            cardContent.append(errorMessage);
        } finally {
           // cardContent.classList.remove("loader");
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
            // These modify the plotJson object in place
            adjustExpressionColorbar(plotJson.data);
            adjustClusterColorbars(plotJson.data);
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

    async renderSVG(display, svgScoringMethod="gene") {
        const datasetId = display.dataset_id;
        const plotConfig = display.plotly_config;
        const {gene_symbol: geneSymbol} = plotConfig;

        const data = await apiCallsMixin.fetchSvgData(datasetId, geneSymbol)
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        const plotContainer = document.querySelector(`#tile_${this.tile.tile_id} .card-image`);
        plotContainer.replaceChildren();    // erase plot

        colorSVG(data, plotConfig.colors, this, svgScoringMethod);

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

    // create a legend div about 7% of the card-image height
    const legendDiv = document.createElement('div');
    legendDiv.classList.add('legend');
    //legendDiv.style.height = '10%';
    cardImage.append(legendDiv);

    // create a svg div about 92% of the card-image height
    const svgDiv = document.createElement('div');
    svgDiv.classList.add('svg');
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

                    // Add data-tooltip and "has-tooltip-top" class to each path so that Bulma can render the tooltip
                    path.attr('data-tooltip', tooltipText);
                    path.addClass('has-tooltip-top');
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

                    // Add data-tooltip and "has-tooltip-top" class to each path so that Bulma can render the tooltip
                    // TODO: doesn't work
                    path.attr('data-tooltip', tooltipText);
                    path.addClass('has-tooltip-top');
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
        .style('max-width', '100%')
        .attr('viewbox', `0 0 ${card.getBoundingClientRect().width} 40`)
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

    // TODO: Issues to resolve here.  The 'atf4' gene in this datasets:
    //  The Adult Cochlea Response to PTS-Inducing Noise - Summary View
    // Has a data range of around -0.34 up but here the min is returning as
    // 0. Is this being converted to an into somewhere?
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

    const width = card.getBoundingClientRect().width;

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