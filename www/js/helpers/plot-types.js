/**
 * plot-types.js - Browser lookups into the shared plot-type registry, /plot_types.json (#1380).
 *
 * The registry lists every display plot type with its label, family (which backend draws it) and
 * legacy aliases. The Python backend reads the same file through lib/gear/plot_types.py, so the
 * lookups here mirror that module.
 *
 * The registry is fetched once, when this module is first imported. Modules that import it wait for
 * the fetch (top-level await), so the lookups below can be used synchronously.
 */

// Families whose plot types show several genes at once
const MULTIGENE_FAMILIES = ["mg_plotly", "mg_scanpy"];

/**
 * Fetches the plot-type registry.
 * "no-cache" asks the browser to check for a newer copy each time, so an edited registry is picked up
 * without bumping cache_version.json (only page entry scripts are cache-busted).
 * @returns {Promise<Object>} Plot type entries keyed by canonical name, or {} if the registry could not be loaded.
 */
const loadPlotTypeRegistry = async () => {
    try {
        const response = await fetch("/plot_types.json", { cache: "no-cache" });
        if (!response.ok) {
            throw new Error(`HTTP ${response.status}`);
        }
        const { plot_types: plotTypes } = await response.json();
        return plotTypes;
    } catch (error) {
        console.error("Could not load the plot-type registry (/plot_types.json).", error);
        return {};
    }
};

const registry = await loadPlotTypeRegistry();

// Every canonical name and alias, mapped to its canonical name
const canonicalNames = new Map();
for (const [name, entry] of Object.entries(registry)) {
    canonicalNames.set(name, name);
    for (const alias of entry.aliases || []) {
        canonicalNames.set(alias, name);
    }
}

/**
 * Returns the canonical name for a plot type or one of its aliases, ignoring case and whitespace.
 * @param {string} plotType - A plot type name or alias (e.g. "tsne", "tsne_dyna").
 * @returns {string|null} The canonical name, or null for an empty or unknown plot type.
 */
export const normalizePlotType = (plotType) => {
    if (!plotType) return null;
    return canonicalNames.get(plotType.trim().toLowerCase()) ?? null;
};

/**
 * Returns the family (plotly, scanpy, mg_plotly, mg_scanpy, svg, gosling or spatial) of a plot type.
 * @param {string} plotType - A plot type name or alias.
 * @returns {string|null} The family, or null if the plot type is unknown.
 */
export const getPlotTypeFamily = (plotType) => registry[normalizePlotType(plotType)]?.family ?? null;

/**
 * Returns true if the plot type shows several genes at once.
 * @param {string} plotType - A plot type name or alias.
 * @returns {boolean}
 */
export const isMultigenePlotType = (plotType) => MULTIGENE_FAMILIES.includes(getPlotTypeFamily(plotType));

/**
 * Returns the human-readable label of a plot type, falling back to the given name if it is unknown.
 * @param {string} plotType - A plot type name or alias.
 * @returns {string}
 */
export const getPlotTypeLabel = (plotType) => registry[normalizePlotType(plotType)]?.label ?? plotType;
