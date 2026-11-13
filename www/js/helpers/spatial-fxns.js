"use strict";

/*
Spatial panel rendering for dataset tiles.

The spatial panel is a Panel (Bokeh server) app in services/spatial. The API returns a <script> that loads
the app into the tile. The app talks back to this page through globals and events set in
services/spatial/panel_common.py:
- window.gearSpatialViewState[datasetId]: current view range, so switching genes keeps the view
- window.gearSpatialUrlParams: current app parameters, used by the HTML download
- "save_spatial_display_<datasetId>" window event: the app's "save display" button (zoomed view only)
*/

import { apiCallsMixin, createToast, logErrorInConsole } from "../common.v2.js";

const RENDER_POLL_INTERVAL = 100;   // ms between checks for the rendered app
const RENDER_MAX_ATTEMPTS = 1800;   // 3 minute maximum wait (100 ms * 1800)

// "save display" listeners, by tile ID, so a re-render replaces the tile's previous listener
const saveHandlers = new Map();  // tileId -> {eventName, handler}

/**
 * Renders a spatial panel display into a tile's card image.
 * Errors are shown as a toast instead of being thrown.
 *
 * @param {HTMLElement} cardImage - The tile's .card-image element. Its contents are replaced.
 * @param {Object} plotConfig - The display's plotly_config (gene_symbol, projection_id, vmin, vmax, ...).
 * @param {Object} opts - Render options.
 * @param {string} opts.tileId - The ID of the dataset tile.
 * @param {string} opts.datasetId - The dataset ID.
 * @param {boolean} opts.isZoomed - True for the zoomed (expanded) view, which allows saving the display.
 * @param {boolean} opts.isMultigene - True for multi-gene tiles (saving is not supported there).
 * @param {function(string, string): void} opts.showCardMessage - Shows a message in the tile, given (level, message).
 * @param {Object} [otherOpts={}] - Request options (e.g. an abort signal).
 * @returns {Promise<void>}
 */
export const renderSpatialPanel = async (cardImage, plotConfig, { tileId, datasetId, isZoomed, isMultigene, showCardMessage }, otherOpts = {}) => {
    showCardMessage("info", "Loading spatial display...");

    const prepPlotConfig = {
        gene_symbol: plotConfig.gene_symbol,
        projection_id: plotConfig.projection_id,
        is_zoomed: isZoomed,
        expression_min_clip: plotConfig.expression_min_clip,
        vmin: plotConfig.vmin,  // User-defined expression colorbar bounds
        vmax: plotConfig.vmax,
        disable_save: !apiCallsMixin.sessionId && isZoomed,  // If not logged in, then do not allow saving the display
    }

    // Grab the saved range from the global state if it exists and apply it to the prepPlotConfig
    const savedRange = window.gearSpatialViewState?.[datasetId];
    if (savedRange) {
        Object.assign(prepPlotConfig, savedRange);
    }

    try {
        const data = await apiCallsMixin.prepSpatialPanelData(datasetId, prepPlotConfig, otherOpts);
        if (data?.success < 1) {
            throw new Error (data?.message ? data.message : "Unknown error.")
        }
        if (! data?.script) {
            throw new Error("Could not retrieve script for spatial display.");
        }

        const parser = new DOMParser();
        const doc = parser.parseFromString(data.script, 'text/html');
        const parsedScript = doc.querySelector('script');

        cardImage.replaceChildren();

        // Enforce boundaries so Panel's stretch_both has walls to push against
        cardImage.style.display = "block";
        cardImage.style.width = "100%";
        if (!isZoomed) {
            cardImage.style.minHeight = "360px";
        }

        // inject the script into the card to trigger loading the Panel app (Bokeh server_document)
        if (parsedScript) {
            injectPanelScript(cardImage, parsedScript, showCardMessage);
        }

        registerSaveHandler({ tileId, datasetId, isZoomed, isMultigene, plotConfig });

    } catch (error) {
        createToast(`Error rendering spatial display: ${error.message}`, "is-danger");
        console.trace(error);
    }
}

/**
 * Adds the Panel app script to the card, behind a loading overlay that is removed once the app has rendered.
 *
 * @param {HTMLElement} cardImage - The tile's .card-image element.
 * @param {HTMLScriptElement} parsedScript - The script parsed from the API response.
 * @param {function(string, string): void} showCardMessage - Shows a message in the tile, given (level, message).
 */
const injectPanelScript = (cardImage, parsedScript, showCardMessage) => {
    // A parsed <script> does not run when inserted, so copy it into a new element
    const scriptElement = document.createElement('script');

    // Copy attributes over
    Array.from(parsedScript.attributes).forEach(attr => {
        scriptElement.setAttribute(attr.name, attr.value);
    });

    // Copy the actual inline JavaScript payload
    scriptElement.textContent = parsedScript.textContent;

    // Create a absolute overlay to hide the empty white layout
    const loader = document.createElement('div');
    loader.innerHTML = `
        <div style="width: 3rem; height: 3rem; border: 4px solid #f3f3f3; border-top: 4px solid #544A8E; border-radius: 50%; animation: spin 1s linear infinite;"></div>
        <div style="margin-top: 15px; font-weight: bold; color: #444;">Loading spatial architecture...</div>
        <style>@keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }</style>
    `;
    loader.style.position = "absolute";
    loader.style.top = "0";
    loader.style.left = "0";
    loader.style.width = "100%";
    loader.style.height = "100%";
    loader.style.backgroundColor = "white"; // Solid white perfectly hides the rendering boxes
    loader.style.display = "flex";
    loader.style.flexDirection = "column";
    loader.style.justifyContent = "center";
    loader.style.alignItems = "center";
    loader.style.zIndex = "9999";

    // The parent container MUST be relative so the absolute overlay fits perfectly inside it
    cardImage.style.position = "relative";
    cardImage.appendChild(loader);

    // Poll for the rendered app, with a timeout if it takes too long.
    let attempts = 0;
    const checkRender = setInterval(() => {
        attempts++;

        // Look for ANY physical plot boundary Bokeh creates, not just a canvas.
        const plotElements = cardImage.querySelectorAll('.bk-panel-models-layout-Column');

        if (plotElements.length > 0) {
            // The DOM structure has arrived! Stop polling.
            clearInterval(checkRender);

            loader.style.transition = "opacity 0.3s ease";
            loader.style.opacity = "0";

            // Give Datashader/WebGL an extra 300ms to push the final pixel colors
            // into those containers before lifting the curtain.
            setTimeout(() => loader.remove(), 300);
        }
        else if (attempts >= RENDER_MAX_ATTEMPTS) {
            // FAILSAFE: If the backend crashed silently, do not trap the user.
            clearInterval(checkRender);
            showCardMessage("danger", "Spatial display failed to load. Please check the dataset or refresh the page.");
            setTimeout(() => loader.remove(), 200);
        }
    }, RENDER_POLL_INTERVAL);

    // Append and execute
    cardImage.appendChild(scriptElement);
}

/**
 * Listens for the Panel app's "save display" event. Only the zoomed view can save.
 * Any listener from the tile's previous render is removed first, to prevent double saves.
 *
 * @param {Object} opts - Options.
 * @param {string} opts.tileId - The ID of the dataset tile.
 * @param {string} opts.datasetId - The dataset ID.
 * @param {boolean} opts.isZoomed - True for the zoomed (expanded) view.
 * @param {boolean} opts.isMultigene - True for multi-gene tiles.
 * @param {Object} opts.plotConfig - The display's plotly_config, saved as the new display's config.
 */
const registerSaveHandler = ({ tileId, datasetId, isZoomed, isMultigene, plotConfig }) => {
    const previous = saveHandlers.get(tileId);
    if (previous) {
        window.removeEventListener(previous.eventName, previous.handler);
        saveHandlers.delete(tileId);
    }

    if (!isZoomed) return;

    const eventName = `save_spatial_display_${datasetId}`;
    const handler = async (event) => {
        const { displayName, makeDefault } = event.detail;

        try {
            if (!apiCallsMixin.sessionId) {
                createToast("Must be logged in to save as a display.");
                throw new Error("Must be logged in to save as a display.");
            }
            await saveSpatialParameters({ datasetId, isMultigene, displayName, makeDefault, plotConfig });
            createToast("Spatial display saved successfully!", "is-success");
        } catch (error) {
            console.error(error);
        }
    };
    window.addEventListener(eventName, handler);
    saveHandlers.set(tileId, { eventName, handler });
}

/**
 * Downloads the current spatial panel as an HTML file.
 *
 * @param {string} shareId - The dataset share ID, used in the file name.
 * @param {string} geneSymbol - The gene symbol, used in the file name.
 * @returns {Promise<void>}
 */
export const downloadSpatialHTML = async (shareId, geneSymbol) => {
    // These are grabbed from the state sync within the Panel app
    if (!window.gearSpatialUrlParams) {
        createToast("Cannot download HTML because spatial display parameters are not available.", "is-warning");
        return;
    }

    const urlParams = window.gearSpatialUrlParams;
    urlParams.append("_", Date.now());   // add timestamp to prevent caching issues

    // Must hit regular URL and not websocket (ws) version
    // because the HTTP version is unidirectional and will return response data
    const url = `/panel/spatial_download?${urlParams.toString()}`;
    try {
        const response = await fetch(url, {
            method: "GET",
            headers: {
                "Content-Type": "text/html",
                "X-Session-ID": apiCallsMixin.sessionId || "",
            },
        });

        if (!response.ok) {
            throw new Error(`Error downloading HTML: ${response.statusText}`);
        }

        // Wanted to convert to PDF but the plots are in JS, which is not executed inside the returned HTML
        // when added to the DOM. So HTML will have to do.

        // The response will be a blob representing the HTML file
        const blob = await response.blob();
        const downloadUrl = URL.createObjectURL(blob);

        // Create a hidden link to trigger the download
        const hiddenLink = document.createElement("a");
        hiddenLink.href = downloadUrl;
        hiddenLink.download = `${shareId}_${geneSymbol}_spatial.html`;
        document.body.appendChild(hiddenLink);
        hiddenLink.click();

        // Clean up
        URL.revokeObjectURL(downloadUrl);
        hiddenLink.remove();

    } catch (error) {
        console.error(error);
        createToast(`Error downloading PNG: ${error.message}`, "is-danger");
    }
}

/**
 * Saves the current spatial panel view as a new display, and optionally makes it the user's default.
 *
 * @param {Object} opts - Options.
 * @param {string} opts.datasetId - The dataset ID.
 * @param {boolean} opts.isMultigene - True for multi-gene tiles (not supported).
 * @param {string} opts.displayName - Name for the new display.
 * @param {boolean} opts.makeDefault - True to make the new display the user's default.
 * @param {Object} opts.plotConfig - The display's plotly_config. The current view range is added to it.
 * @returns {Promise<void>}
 * @throws {Error} For multi-gene tiles.
 */
export const saveSpatialParameters = async ({ datasetId, isMultigene, displayName, makeDefault, plotConfig }) => {
    if (isMultigene) {
        throw new Error("Saving spatial display is not supported for multi-gene displays.");
    }
    const plotType = "spatial_panel";

    const savedRange = window.gearSpatialViewState?.[datasetId];
    if (savedRange) {
        Object.assign(plotConfig, savedRange);
    }

    try {
        const {display_id: displayId, success: saveSuccess} = await apiCallsMixin.saveDatasetDisplay(datasetId, null, displayName, plotType, plotConfig);
        if (!saveSuccess) {
            throw new Error("Could not save this new display. Please contact the gEAR team.");
        }
        createToast("Display saved.", "is-success");

        if (!makeDefault) return;

        const {success: defaultSuccess} = await apiCallsMixin.saveDefaultDisplay(datasetId, displayId, 0);  // 0 = single-gene
        if (!defaultSuccess) {
            throw new Error("Could not make this display as your default, but it is saved. Please contact the gEAR team.");
        }

    } catch (error) {
        logErrorInConsole(error);
        createToast(error);
    }
}
