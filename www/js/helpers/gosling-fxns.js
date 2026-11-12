"use strict";

/*
Gosling (epigenome) plot rendering for dataset tiles.

renderGoslingPlot() embeds a Gosling spec built by /api/plot/<dataset_id>/gosling
(www/api/resources/gosling_spec.py) and adds the controls around it:
- pan and zoom buttons for each panel, plus keyboard shortcuts (←/→ pan, +/− zoom)
- in the zoomed (expanded) view, the "View in UCSC Genome Browser" button and the Panel B gene search
*/

import { apiCallsMixin, createToast, logErrorInConsole } from "../common.v2.js";

const GOSLING_MODULE_URL = "https://esm.sh/gosling.js@1.0.6";

// Pan and zoom step sizes match the UCSC Genome Browser defaults.
const GOSLING_PAN_FRACTION = 0.1;   // pan by 10% of the visible region
const GOSLING_ZOOM_FACTOR = 3;      // zoom in or out 3x around the center
const GOSLING_NAV_DURATION = 300;   // animation length (ms)
const GOSLING_MIN_REGION = 20;      // narrowest region (bp) that zooming in will reach
const GOSLING_BASE_PADDING = 1500;  // matches zoom_view_to_domain() in www/api/resources/gosling_spec.py
const GOSLING_HIC_PADDING_MULTIPLIER = 500;

// Gosling's "location" PubSub is shared by every embedded Gosling plot on the page, and its unsubscribe()
// removes every subscriber. So we subscribe once and route events to the tile whose id namespace prefixes the event id.
const goslingLocationHandlers = new Map();  // id namespace -> callback(eventData)
let goslingLocationSubscribed = false;

/**
 * Loads the Gosling library and returns its embed function.
 *
 * @returns {Promise<Function>} Gosling's embed(element, spec, options).
 * @throws {Error} If the library cannot be loaded.
 */
export const loadGoslingEmbed = async () => {
    const gos = await import(GOSLING_MODULE_URL);
    // prefer the named export, then try the default export
    const embedFn = gos.embed ?? gos.default?.embed;
    if (!embedFn) throw new Error("Gosling embed function not found.");
    return embedFn;
}

/**
 * Stops routing Gosling "location" events to a tile, e.g. when the tile shows a different display.
 *
 * @param {string} tileId - The ID of the dataset tile.
 */
export const clearGoslingTile = (tileId) => {
    goslingLocationHandlers.delete(getGoslingIdNamespace(tileId));
}

/**
 * Embeds a Gosling plot in a tile and adds the pan/zoom controls, keyboard shortcuts and,
 * for the zoomed view, the UCSC export button and Panel B gene search.
 *
 * @param {Function} embedFn - Gosling's embed function, from loadGoslingEmbed().
 * @param {HTMLElement} plotContainer - The tile's .card-image element. Its contents are replaced.
 * @param {Object} spec - The Gosling spec from the API. Its ids are prefixed in place.
 * @param {Object} opts - Plot options.
 * @param {string} opts.tileId - The ID of the dataset tile.
 * @param {boolean} opts.zoom - True for the zoomed view (Panel A and Panel B).
 * @param {string} opts.position - Gene position from the API ("assembly.chromosome:start-end"), or "NA".
 * @param {boolean} opts.hicFound - True if the spec has a Hi-C track (which is shown with more padding).
 * @param {string} opts.assembly - Genome assembly (e.g. "mm10").
 * @param {string} [opts.hubUrl] - UCSC track hub URL, for the UCSC export button.
 * @returns {Promise<Object>} The Gosling API for the embedded plot.
 * @throws {Error} If Gosling fails to render the plot.
 */
export const renderGoslingPlot = async (embedFn, plotContainer, spec, { tileId, zoom, position, hicFound, assembly, hubUrl }) => {
    const positionArr = [position, zoom ? position : ""]; // [leftPosition, rightPosition], for the UCSC export

    // Give this tile's track and view ids a unique prefix, so its "location" events can be told apart
    // from those of other Gosling tiles on the page (the event PubSub is shared).
    const idNamespace = getGoslingIdNamespace(tileId);
    namespaceGoslingIds(spec, idNamespace);
    const chromSizes = getGoslingChromSizes(spec.assembly);

    // Current visible region of each panel, kept in sync with mouse zooming and panning through "location" events
    const padding = hicFound ? GOSLING_BASE_PADDING * GOSLING_HIC_PADDING_MULTIPLIER : GOSLING_BASE_PADDING;
    const initialRegion = parseGoslingPosition(position, padding);
    const panels = zoom ? ["left", "right"] : ["left"];
    const regions = { left: initialRegion, right: zoom ? initialRegion : null };

    let goslingApi = null;  // set once the plot is embedded

    // Pan or zoom a panel. The annotation track is zoom-linked to the rest of its panel, so moving it moves them all.
    const navigate = (panel, action) => {
        if (!goslingApi) return;
        const target = computeGoslingNavRegion(regions[panel], action, chromSizes);
        if (!target) return;
        // Update now so quick repeated clicks build on each other; location events overwrite this as the view moves
        regions[panel] = { chromosome: target.chromosome, start: target.start, endChromosome: target.chromosome, end: target.end };
        try {
            goslingApi.zoomTo(`${idNamespace}${panel}-annotation`, `${target.chromosome}:${target.start}-${target.end}`, 0, GOSLING_NAV_DURATION);
        } catch (error) {
            logErrorInConsole(error);
        }
    };

    const navBar = document.createElement("div");
    navBar.classList.add("gosling-nav-bar", "is-flex", "is-flex-wrap-wrap", "is-justify-content-space-between", "is-align-items-center", "px-2");
    const navGroups = {};
    for (const panel of panels) {
        const label = zoom ? (panel === "left" ? "Panel A" : "Panel B") : "";
        navGroups[panel] = createGoslingNavControls(label, (action) => navigate(panel, action));
        navBar.append(navGroups[panel]);
    }
    // A div, not a p: ".js-tile .card-image p" styles the tile's loading/error messages
    const navHint = document.createElement("div");
    navHint.classList.add("help", "has-text-grey-dark", "px-2", "mt-0");

    // Panel that the keyboard shortcuts act on. In the zoomed view, clicking a panel selects it.
    let keyboardPanel = "left";
    const setKeyboardPanel = (panel) => {
        keyboardPanel = panel;
        if (!zoom) {
            navHint.textContent = "Drag a track to pan and scroll to zoom. Or click the plot, then use ←/→ to pan and +/− to zoom.";
            return;
        }
        const label = panel === "left" ? "Panel A" : "Panel B";
        navHint.textContent = `Drag a track to pan and scroll to zoom. Or click Panel A or Panel B, then use ←/→ to pan and +/− to zoom it (now: ${label}).`;
        // Highlight the label of the panel the keyboard acts on
        for (const p of panels) {
            const labelElt = navGroups[p].querySelector(".js-gosling-nav-label");
            labelElt?.classList.toggle("has-text-link", p === panel);
            labelElt?.classList.toggle("has-text-grey-dark", p !== panel);  // Bulma's grey helpers override has-text-link
        }
    };
    setKeyboardPanel("left");

    const setNavEnabled = (panel) => {
        for (const button of navGroups[panel].querySelectorAll("button")) {
            button.disabled = !regions[panel];
        }
    };
    panels.forEach(setNavEnabled);

    const goslingContainer = document.createElement("div");
    goslingContainer.id = `tile-${tileId}-gosling`;
    goslingContainer.classList.add("gosling-nav-target");
    goslingContainer.tabIndex = 0;  // focusable, for keyboard panning and zooming
    goslingContainer.style.width = "100%";
    goslingContainer.style.boxSizing = "border-box";
    goslingContainer.style.marginTop = "5px";
    plotContainer.replaceChildren();    // erase card message
    plotContainer.append(navBar, navHint, goslingContainer);

    // Register before embedding so events from the first render are not missed once subscribed
    goslingLocationHandlers.set(idNamespace, (eventData) => {
        const range = eventData?.genomicRange;
        if (!Array.isArray(range) || range.length !== 2) return;
        const panel = panels.find((p) => eventData.id.startsWith(`${idNamespace}${p}-`));
        if (!panel) return;  // e.g. the static genome-wide and chromosome views
        const hadRegion = Boolean(regions[panel]);
        regions[panel] = {
            chromosome: range[0].chromosome,
            start: range[0].position,
            endChromosome: range[1].chromosome,
            end: range[1].position,
        };
        if (!hadRegion) setNavEnabled(panel);
    });

    // Themes -> https://gosling-lang.org/themes/
    const embedOpts = { "padding": 0, "theme": null };
    // NOTE: re-embedding does work but it causes some stability issues
    try {
        goslingApi = await embedFn(goslingContainer, spec, embedOpts);
    } catch (error) {
        goslingLocationHandlers.delete(idNamespace);
        throw error;
    }

    if (!goslingLocationSubscribed) {
        goslingLocationSubscribed = true;
        goslingApi.subscribe("location", (_msg, eventData) => {
            if (typeof eventData?.id !== "string") return;
            for (const [namespace, handler] of goslingLocationHandlers) {
                if (eventData.id.startsWith(namespace)) handler(eventData);
            }
        });
    }

    /**
     * Returns the panel ("left" or "right") under a pointer event, or null if the pointer is not over
     * the panels (e.g. over the genome-wide or chromosome views above them).
     * Track shapes are in pixels from the top-left of the Gosling component.
     */
    const getPanelAtPointer = (event) => {
        const component = goslingContainer.querySelector(".gosling-component");
        const rightShape = goslingApi.getTrack(`${idNamespace}right-annotation`)?.shape;
        const leftShape = goslingApi.getTrack(`${idNamespace}left-annotation`)?.shape;
        if (!component || !rightShape || !leftShape) return null;
        const rect = component.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;
        if (y < Math.min(leftShape.y, rightShape.y)) return null;
        return x >= rightShape.x ? "right" : "left";
    };

    // HiGlass may cancel mousedown while dragging, which would block the click from focusing the plot, so focus it here.
    goslingContainer.addEventListener("pointerdown", (event) => {
        if (event.target.closest("input, textarea, select, button")) return;
        if (zoom) {
            const panel = getPanelAtPointer(event);
            if (panel && panel !== keyboardPanel) setKeyboardPanel(panel);
        }
        goslingContainer.focus({ preventScroll: true });
    }, { capture: true });
    const navKeys = { ArrowLeft: "left", ArrowRight: "right", "+": "in", "=": "in", "-": "out", "_": "out" };
    goslingContainer.addEventListener("keydown", (event) => {
        const action = navKeys[event.key];
        // Leave browser shortcuts (e.g. Ctrl/Cmd +/- page zoom) and typing in the Panel B search box alone
        if (!action || event.ctrlKey || event.metaKey || event.altKey) return;
        if (event.target.closest("input, textarea, select")) return;
        event.preventDefault();  // keep the arrow keys from scrolling the page
        navigate(keyboardPanel, action);
    });

    if (!zoom) {
        return goslingApi;
    }

    // If the view is a zoomed view extra controls and events are added.
    const exportButton = createExportButton(tileId, goslingContainer.id);
    const searchButton = createPanelBSearchBox(tileId, exportButton.id);

    if (hubUrl) {
        // Add event listener to export button to open UCSC Genome Browser with hub URL and position
        exportButton.addEventListener('click', () => {
            const url = "https://genome.ucsc.edu/cgi-bin/hgTracks";
            // add DB and Hub parameters
            const urlParams = new URLSearchParams({
                db: assembly,
                hubUrl,    // preconfigured hub file
                ignoreCookie: 1, // Ignore cookie to ensure trackHub changes are respected. Unfortunately, adds default tracks.
            })

            // if an element in positionArr is not empty, add it to the URL, separated by a pipe
            const highlightedPositions = positionArr.filter(pos => pos).join('|');
            urlParams.set('highlight', highlightedPositions);

            // open in a new tab
            window.open(`${url}?${urlParams.toString()}`, '_blank');
        });
    }

    // Add event listener to search button to zoom to gene in Panel B
    searchButton.addEventListener('click', async () => {
        const geneInput = document.getElementById(`tile-${tileId}-panel-b-gene-input`);
        const gene = geneInput.value.trim();
        if (!gene) {
            alert("Please enter a gene name.");
            return;
        }

        let panelBData = null;
        try {
            panelBData = await apiCallsMixin.fetchHiglassGeneCoords(gene, assembly);
        } catch (error) {
            console.error("Error searching for gene:", error);
            createToast(`An error occurred while searching for gene: ${gene}`);
        }

        if (!panelBData) {
            alert(`Gene ${gene} not found for assembly ${assembly}.`);
            return;
        }

        const chr = panelBData.chr
        const start = panelBData.txStart
        const end = panelBData.txEnd

        const rightPosition = `${chr}:${start}-${end}`;
        positionArr[1] = `${assembly}.${rightPosition}`;
        await goslingApi.zoomTo(`${idNamespace}right-annotation`, rightPosition, padding); // track name, position, padding, duration (ms)

        // TODO:  if Hi-C data is found, add an annotation to the gene position
    });

    return goslingApi;
}

/**
 * Returns the prefix added to a tile's Gosling track and view ids.
 *
 * @param {string} tileId - The ID of the dataset tile.
 * @returns {string} The id namespace.
 */
export const getGoslingIdNamespace = (tileId) => `tile-${tileId}--`;  // no ":" or "." in case HiGlass uses ids in CSS selectors

/**
 * Prefixes every "id" in a Gosling spec, in place. "linkingId" values are left alone, since they only
 * need to be unique within one embedded plot.
 *
 * @param {*} node - The spec, or a part of it.
 * @param {string} namespace - The prefix to add.
 */
export const namespaceGoslingIds = (node, namespace) => {
    if (Array.isArray(node)) {
        node.forEach((child) => namespaceGoslingIds(child, namespace));
        return;
    }
    if (!node || typeof node !== "object") return;
    for (const [key, value] of Object.entries(node)) {
        if (key === "id" && typeof value === "string") {
            node[key] = `${namespace}${value}`;
        } else {
            namespaceGoslingIds(value, namespace);
        }
    }
}

/**
 * Builds a chromosome -> size map from a Gosling spec's "assembly", which the server sends as
 * an array of [chromosome, size] pairs.
 *
 * @param {Array|Object|string} assembly - The spec's assembly.
 * @returns {Map<string, number>} Chromosome sizes. Empty if the assembly has no sizes.
 */
export const getGoslingChromSizes = (assembly) => {
    const pairs = Array.isArray(assembly) ? assembly : assembly?.chromSizes;
    return new Map(Array.isArray(pairs) ? pairs.map(([chrom, size]) => [chrom, Number(size)]) : []);
}

/**
 * Parses a position string from the Gosling API ("assembly.chromosome:start-end") into a region.
 *
 * @param {string} positionStr - The position string, or "NA" when the gene was not found.
 * @param {number} [padding=0] - Bases to add to each side, to match the server's initial zoom.
 * @returns {Object|null} The region ({chromosome, start, endChromosome, end}), or null if it cannot be parsed.
 */
export const parseGoslingPosition = (positionStr, padding = 0) => {
    if (!positionStr || positionStr === "NA") return null;
    // The assembly comes before the first "."; chromosome names may contain "." themselves
    const location = positionStr.slice(positionStr.indexOf(".") + 1);
    const match = location.match(/^(.+):(\d+)-(\d+)$/);
    if (!match) return null;
    const [, chromosome, start, end] = match;
    return {
        chromosome,
        start: Math.max(1, Number(start) - padding),
        endChromosome: chromosome,
        end: Number(end) + padding,
    };
}

/**
 * Computes the region to show after a pan or zoom step, kept within one chromosome.
 *
 * @param {Object|null} region - The current region ({chromosome, start, endChromosome, end}).
 * @param {string} action - "left", "right", "in", or "out".
 * @param {Map<string, number>} chromSizes - Chromosome sizes from getGoslingChromSizes().
 * @returns {Object|null} The new region ({chromosome, start, end}), or null if there is no current region.
 */
export const computeGoslingNavRegion = (region, action, chromSizes) => {
    if (!region) return null;
    const chromosome = region.chromosome;
    const chromSize = chromSizes.get(chromosome) ?? Infinity;
    // If the view spans chromosomes (only reachable with the mouse), stay on the left-most one
    let start = region.start;
    let end = region.endChromosome === chromosome ? region.end : chromSize;
    if (!Number.isFinite(end)) return null;
    const width = Math.max(1, end - start);

    if (action === "left" || action === "right") {
        const shift = Math.max(1, Math.round(width * GOSLING_PAN_FRACTION));
        start += action === "left" ? -shift : shift;
        end = start + width;
    } else if (action === "in" || action === "out") {
        const center = (start + end) / 2;
        let newWidth = action === "in" ? width / GOSLING_ZOOM_FACTOR : width * GOSLING_ZOOM_FACTOR;
        newWidth = Math.min(Math.max(newWidth, GOSLING_MIN_REGION), chromSize);
        start = Math.round(center - newWidth / 2);
        end = Math.round(start + newWidth);
    } else {
        return null;
    }

    // Slide back inside the chromosome while keeping the width
    if (start < 1) {
        end += 1 - start;
        start = 1;
    }
    if (end > chromSize) {
        start = Math.max(1, start - (end - chromSize));
        end = chromSize;
    }
    return { chromosome, start, end };
}

/**
 * Creates a group of pan and zoom buttons for one Gosling panel.
 *
 * @param {string} label - Text shown before the buttons (e.g. "Panel A"), or "" for none.
 * @param {function(string): void} onNavigate - Called with "left", "right", "in", or "out".
 * @returns {HTMLElement} The button group.
 */
const createGoslingNavControls = (label, onNavigate) => {
    const group = document.createElement("div");
    group.classList.add("is-flex", "is-align-items-center", "my-1");

    if (label) {
        const labelElt = document.createElement("span");
        labelElt.classList.add("js-gosling-nav-label", "has-text-weight-semibold", "has-text-grey-dark", "is-size-7", "mr-2");
        labelElt.textContent = label;
        group.append(labelElt);
    }

    const buttons = document.createElement("div");
    buttons.classList.add("buttons", "has-addons", "mb-0");
    const actions = [
        ["left", "mdi-chevron-left", `Pan ${label || "view"} left`],
        ["right", "mdi-chevron-right", `Pan ${label || "view"} right`],
        ["in", "mdi-magnify-plus-outline", `Zoom ${label || "view"} in`],
        ["out", "mdi-magnify-minus-outline", `Zoom ${label || "view"} out`],
    ];
    for (const [action, iconClass, title] of actions) {
        const button = document.createElement("button");
        button.type = "button";
        button.classList.add("button", "is-small", "mb-0");
        button.title = title;
        button.setAttribute("aria-label", title);

        const icon = document.createElement("span");
        icon.classList.add("icon");
        const iconElt = document.createElement("i");
        iconElt.classList.add("mdi", "mdi-18px", iconClass);
        icon.append(iconElt);
        button.append(icon);

        button.addEventListener("click", () => onNavigate(action));
        buttons.append(button);
    }
    group.append(buttons);
    return group;
}

/**
 * Adds a "View in UCSC Genome Browser" button to the top of an element.
 *
 * @param {string} tileId - The ID of the dataset tile.
 * @param {string} selectorId - ID of the element to add the button to.
 * @returns {HTMLButtonElement} The button.
 */
const createExportButton = (tileId, selectorId) => {
    const exportButton = document.createElement('button');
    exportButton.id = `tile-${tileId}-ucsc-export-button`;
    exportButton.textContent = 'View in UCSC Genome Browser';
    // stylize the button
    exportButton.style.zIndex = '1000';
    exportButton.style.padding = '10px';
    exportButton.style.marginBottom = '10px';
    exportButton.style.marginRight = '10px';
    exportButton.style.backgroundColor = 'purple'; // Purple background
    exportButton.style.color = 'white'; // White text
    exportButton.style.border = 'none';
    exportButton.style.borderRadius = '5px';
    exportButton.style.cursor = 'pointer';
    document.getElementById(selectorId).prepend(exportButton);
    return exportButton
}

/**
 * Adds a gene search box for Panel B, just before an element.
 *
 * @param {string} tileId - The ID of the dataset tile.
 * @param {string} selectorId - ID of the element to put the search box before.
 * @returns {HTMLButtonElement} The search button.
 */
const createPanelBSearchBox = (tileId, selectorId) => {
    const searchBox = document.createElement('div');
    searchBox.id = `tile-${tileId}-panel-b-gene-input-search;`
    searchBox.style.float = "right";

    const searchInput = document.createElement('input');
    searchInput.id = `tile-${tileId}-panel-b-gene-input`;
    searchInput.type = 'text';
    searchInput.placeholder = 'Enter gene name';
    searchInput.style.marginRight = '4px'; // for some reason the label and input are not aligned propely with the panel A search box outside the gosling container

    // Add label for the input
    const searchLabel = document.createElement('label');
    searchLabel.setAttribute('for', searchInput.id);
    searchLabel.textContent = 'Search for a gene in Panel B:';
    searchLabel.style.fontWeight = 'bold';
    searchLabel.style.color = "black";
    searchLabel.style.marginRight = '4px';

    const searchButton = document.createElement('button');
    searchButton.id = `tile-${tileId}-panel-b-search-button`;
    searchButton.textContent = 'Search';

    searchBox.appendChild(searchLabel);
    searchBox.appendChild(searchInput);
    searchBox.appendChild(searchButton);
    document.getElementById(selectorId).before (searchBox);
    return searchButton;
}
