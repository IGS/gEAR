"use strict";

/*
Dataset gene lists: fetching the genes in a dataset (or analysis), finding genes similar to a searched one,
and links to view the whole list or search for one of its genes.

Used by dataset tiles (gene-not-found help, Dataset Information modal) and the curators (gene pickers).
*/

import { apiCallsMixin } from "../common.v2.js";

// Dataset types with no gene expression matrix to list genes from
const NO_GENE_MATRIX_TYPES = ["epiviz", "gosling"];

// Primary-dataset gene lists, by dataset ID. Primary data does not change during a page visit,
// so tiles of the same dataset share one request. Analyses are not cached: they can be re-saved in the workbench.
const primaryGeneCache = new Map();  // datasetId -> Promise<string[]>

// Shortest shared prefix that counts as a similar gene (e.g. "MYO" for MYO7 and MYO6).
// Many symbols are only 3 characters (Ace, Bmf, Cel), so a 3-character search needs only 2 (e.g. "Soz" -> Sox2).
const MIN_SHARED_PREFIX = 3;
const SHORT_QUERY_LENGTH = 3;
const SHORT_QUERY_MIN_SHARED_PREFIX = 2;
// Largest edit distance that counts as a similar gene (typos). Short symbols allow fewer edits,
// since 2 edits can turn a 4-letter symbol into an unrelated one.
const MAX_EDIT_DISTANCE = 2;
const SHORT_SYMBOL_LENGTH = 4;
const SHORT_SYMBOL_MAX_EDIT_DISTANCE = 1;

/**
 * Returns whether a dataset has a gene expression matrix whose genes can be listed.
 *
 * @param {Object} dataset - Dataset metadata (dtype, has_h5ad).
 * @returns {boolean}
 */
export const hasGeneList = (dataset) => {
    if (!dataset || NO_GENE_MATRIX_TYPES.includes(dataset.dtype)) return false;
    return Boolean(dataset.has_h5ad) || dataset.dtype === "spatial";
}

/**
 * Fetches the unique gene symbols in a dataset or one of its analyses, in file order.
 * Missing symbols (sent as "None") are left out.
 *
 * @param {string} datasetId - The dataset ID.
 * @param {string|null} [analysisId=null] - The analysis ID, or null for the primary dataset.
 * @returns {Promise<string[]>} The gene symbols.
 * @throws {Error} If the server could not list the genes.
 */
export const fetchDatasetGeneSymbols = (datasetId, analysisId = null) => {
    if (analysisId) {
        return requestGeneSymbols(datasetId, analysisId);
    }

    if (!primaryGeneCache.has(datasetId)) {
        const request = requestGeneSymbols(datasetId, null);
        // Do not keep a failed request, so a later call can try again
        request.catch(() => primaryGeneCache.delete(datasetId));
        primaryGeneCache.set(datasetId, request);
    }
    return primaryGeneCache.get(datasetId);
}

/**
 * Requests the gene symbols from the API.
 *
 * @param {string} datasetId - The dataset ID.
 * @param {string|null} analysisId - The analysis ID, or null for the primary dataset.
 * @returns {Promise<string[]>} The unique gene symbols, in file order.
 */
const requestGeneSymbols = async (datasetId, analysisId) => {
    const data = await apiCallsMixin.fetchGeneSymbols(datasetId, analysisId);
    if (data?.success !== 1) {
        throw new Error(data?.message || "Could not list the genes in this dataset.");
    }
    // A dataset may repeat a gene symbol
    return [...new Set(data.gene_symbols)].filter((symbol) => symbol !== "None");
}

/**
 * URL of the dataset's gene list as plain text, one gene per line (for opening in a new tab).
 *
 * @param {string} datasetId - The dataset ID.
 * @returns {string}
 */
export const getDatasetGeneListUrl = (datasetId) => {
    return `./api/h5ad/${encodeURIComponent(datasetId)}/genes?format=text`;
}

/**
 * URL of a single-gene expression search for one gene in one dataset.
 *
 * @param {string} geneSymbol - The gene symbol to search for.
 * @param {string} shareId - The dataset share ID.
 * @returns {string}
 */
export const getGeneSearchUrl = (geneSymbol, shareId) => {
    const url = new URL("/expression.html", window.location.origin);
    url.searchParams.append("gene_symbol", geneSymbol);
    url.searchParams.append("gene_symbol_exact_match", "1");
    url.searchParams.append("is_multigene", "0");
    url.searchParams.append("share_id", shareId);
    return url.toString();
}

/**
 * Finds the dataset genes most similar to a searched gene, ignoring case.
 * Ranked by: longest shared prefix (at least 3 characters, or 2 for searches of 3 characters or fewer), then containing the search (or contained in it),
 * then small edit distance (typos; 1 edit for symbols of 4 characters or fewer, else 2). Exact matches are left out, since the search would have found them.
 *
 * @param {string[]} datasetGenes - The dataset's gene symbols.
 * @param {string} query - The searched gene symbol.
 * @param {number} [limit=8] - Most genes to return.
 * @returns {string[]} Similar gene symbols, best first, in their dataset casing.
 */
export const findSimilarGenes = (datasetGenes, query, limit = 8) => {
    const search = (query ?? "").trim().toLowerCase();
    if (!search) return [];

    const maxDistance = search.length <= SHORT_SYMBOL_LENGTH ? SHORT_SYMBOL_MAX_EDIT_DISTANCE : MAX_EDIT_DISTANCE;
    const minPrefix = search.length <= SHORT_QUERY_LENGTH ? SHORT_QUERY_MIN_SHARED_PREFIX : MIN_SHARED_PREFIX;
    const scored = [];
    for (const gene of datasetGenes) {
        const candidate = gene.toLowerCase();
        if (candidate === search) continue;

        const prefix = sharedPrefixLength(candidate, search);
        let rank = null;
        if (prefix >= minPrefix) {
            rank = [0, -prefix];
        } else if (search.length >= MIN_SHARED_PREFIX && (candidate.includes(search) || search.includes(candidate))) {
            rank = [1, 0];
        } else {
            const distance = editDistance(candidate, search, maxDistance);
            if (distance <= maxDistance) {
                rank = [2, distance];
            }
        }
        if (rank) {
            scored.push({ gene, rank, lengthGap: Math.abs(candidate.length - search.length) });
        }
    }

    scored.sort((a, b) =>
        a.rank[0] - b.rank[0]
        || a.rank[1] - b.rank[1]
        || a.lengthGap - b.lengthGap
        || a.gene.localeCompare(b.gene)
    );
    return scored.slice(0, limit).map((item) => item.gene);
}

/**
 * Length of the prefix two strings share.
 *
 * @param {string} a
 * @param {string} b
 * @returns {number}
 */
const sharedPrefixLength = (a, b) => {
    const max = Math.min(a.length, b.length);
    let i = 0;
    while (i < max && a[i] === b[i]) i++;
    return i;
}

/**
 * Levenshtein edit distance between two strings, stopping early once it is over maxDistance.
 *
 * @param {string} a
 * @param {string} b
 * @param {number} maxDistance - Any distance over this is returned as maxDistance + 1.
 * @returns {number}
 */
const editDistance = (a, b, maxDistance) => {
    if (Math.abs(a.length - b.length) > maxDistance) return maxDistance + 1;

    let previous = Array.from({ length: b.length + 1 }, (_, j) => j);
    for (let i = 1; i <= a.length; i++) {
        const current = [i];
        let rowMin = i;
        for (let j = 1; j <= b.length; j++) {
            const cost = a[i - 1] === b[j - 1] ? 0 : 1;
            current[j] = Math.min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + cost);
            rowMin = Math.min(rowMin, current[j]);
        }
        if (rowMin > maxDistance) return maxDistance + 1;
        previous = current;
    }
    return previous[b.length];
}
