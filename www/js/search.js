"use strict";

/*
 * search.js - The unified search page, search.html (#1383, Stage 1).
 *
 * Results come from GET /api/search (www/api/resources/search.py). The search state is kept in
 * the URL (q, types, organism_ids, ownership, date_added, sort, page, sel), so a search can be
 * shared or bookmarked, and the browser's Back and Forward buttons move between searches.
 * Picking a result only changes "sel", without adding a history entry.
 */

import { apiCallsMixin, copyToClipboard, createToast, getCurrentUser, initCommonUI, logErrorInConsole, registerPageSpecificLoginUIUpdates } from "./common.v2.js";

// Result types, in the order of the type filter buttons
const TYPES = ["dataset", "collection", "genelist", "study", "projection", "ccc"];

const TYPE_INFO = {
    dataset: {label: "Dataset", plural: "datasets", icon: "mdi-database-outline"},
    collection: {label: "Collection", plural: "collections", icon: "mdi-layers-outline"},
    genelist: {label: "Gene list", plural: "gene lists", icon: "mdi-format-list-bulleted"},
    study: {label: "Study", plural: "studies", icon: "mdi-file-document-outline"},
    projection: {label: "Projection", plural: "projections", icon: "mdi-hexagon-multiple"},
    ccc: {label: "Cell-cell communication", plural: "cell-cell communication", icon: "mdi-swap-horizontal"},
};

const RESULTS_PER_PAGE = 20;

// URL parameters that make up a search. "sel" (the selected result) is kept apart, because
//  picking a result doesn't run a new search.
const SEARCH_PARAMS = ["q", "types", "organism_ids", "ownership", "date_added", "sort", "page"];

// The current search and its response
let state = {q: "", types: [], organism_ids: "", ownership: "", date_added: "", sort: "relevance", page: 1, sel: null};
let lastResponse = null;

// Ignores the responses of searches overtaken by a newer one (e.g. quick clicks on type buttons)
let searchCounter = 0;

/**
 * Reads the search state from the URL.
 *
 * @returns {Object} - The state, with defaults for anything missing.
 */
const stateFromUrl = () => {
    const params = new URLSearchParams(window.location.search);
    const types = (params.get("types") || "").split(",").filter((type) => TYPES.includes(type));
    const page = Number.parseInt(params.get("page"), 10);
    return {
        q: params.get("q") || "",
        // Every type selected is the same as All
        types: types.length === TYPES.length ? [] : [...new Set(types)],
        organism_ids: params.get("organism_ids") || "",
        ownership: params.get("ownership") || "",
        date_added: params.get("date_added") || "",
        sort: params.get("sort") || "relevance",
        page: page > 0 ? page : 1,
        sel: params.get("sel"),
    };
}

/**
 * Writes the state to the URL, leaving out defaults so links stay short.
 *
 * @param {boolean} addHistoryEntry - true for a new search (Back returns to the previous one),
 *  false for changes within the same search (picking a result).
 */
const writeUrl = (addHistoryEntry) => {
    const params = new URLSearchParams();
    if (state.q) params.set("q", state.q);
    if (state.types.length) params.set("types", state.types.join(","));
    if (state.organism_ids) params.set("organism_ids", state.organism_ids);
    if (state.ownership) params.set("ownership", state.ownership);
    if (state.date_added) params.set("date_added", state.date_added);
    if (state.sort !== "relevance") params.set("sort", state.sort);
    if (state.page > 1) params.set("page", state.page);
    if (state.sel) params.set("sel", state.sel);

    const query = params.toString();
    const url = `${window.location.pathname}${query ? `?${query}` : ""}`;
    if (addHistoryEntry) {
        window.history.pushState({}, "", url);
    } else {
        window.history.replaceState({}, "", url);
    }
}

/**
 * Puts the state into the search box, the type buttons and the filters.
 */
const showStateInControls = () => {
    document.getElementById("search-q").value = state.q;
    document.getElementById("search-sort").value = state.sort;
    document.getElementById("search-date-added").value = state.date_added;
    document.getElementById("search-ownership").value = state.ownership;
    // The organism list may not have this organism (e.g. an edited URL); fall back to All
    const organismSelect = document.getElementById("search-organism");
    organismSelect.value = state.organism_ids;
    if (organismSelect.value !== state.organism_ids) {
        organismSelect.value = "";
    }

    for (const button of document.getElementsByClassName("type-button")) {
        const type = button.dataset.type;
        const selected = type === "all" ? state.types.length === 0 : state.types.includes(type);
        button.setAttribute("aria-pressed", String(selected));
    }
}

/**
 * Runs the search in the state and shows the response.
 */
const runSearch = async () => {
    const thisSearch = ++searchCounter;
    const params = {limit: RESULTS_PER_PAGE};
    for (const name of SEARCH_PARAMS) {
        const value = name === "types" ? state.types.join(",") : state[name];
        if (value) params[name] = value;
    }

    // Show the search is running where the eye is: the results heading gets a spinner and
    //  "Searching…", and the old results fade until the new ones arrive
    const submitButton = document.getElementById("search-submit");
    const loadingIcon = document.getElementById("results-loading");
    submitButton.classList.add("is-loading");
    loadingIcon.classList.remove("is-hidden");
    document.getElementById("results-summary").textContent = "Searching…";
    document.getElementById("results-section").setAttribute("aria-busy", "true");

    try {
        const data = await apiCallsMixin.fetchUnifiedSearch(params);
        if (thisSearch !== searchCounter) return;
        lastResponse = data;
        showResponse();
    } catch (error) {
        if (thisSearch !== searchCounter) return;
        logErrorInConsole(error);
        // The API explains a bad request (400) or a failed search (500) in "problem"
        createToast(error.response?.data?.problem || "The search could not be run. Please try again.");
        document.getElementById("results-summary").textContent = "The search could not be run";
    } finally {
        if (thisSearch === searchCounter) {
            submitButton.classList.remove("is-loading");
            loadingIcon.classList.add("is-hidden");
            document.getElementById("results-section").removeAttribute("aria-busy");
        }
    }
}

/**
 * A short description of the results, e.g. "12 results for “hair cell” in datasets + gene lists".
 */
const summaryText = (total) => {
    let summary = `${total} result${total === 1 ? "" : "s"}`;
    if (state.q) summary += ` for “${state.q}”`;
    if (state.types.length) summary += ` in ${state.types.map((type) => TYPE_INFO[type].plural).join(" + ")}`;
    return summary;
}

/**
 * Shows the counts, results, preview and pagination of the last response.
 */
const showResponse = () => {
    const {counts, results, pagination} = lastResponse;

    for (const countElt of document.querySelectorAll("[data-count-for]")) {
        countElt.textContent = counts[countElt.dataset.countFor] ?? 0;
    }

    document.getElementById("results-summary").textContent = summaryText(pagination.total_results);

    const hasResults = results.length > 0;
    document.getElementById("no-results").classList.toggle("is-hidden", hasResults);
    document.getElementById("results-layout").classList.toggle("is-hidden", !hasResults);
    for (const pager of document.getElementsByClassName("results-pager")) {
        pager.classList.toggle("is-hidden", !hasResults);
    }

    // Keep the selected result if it's on this page, otherwise select the first
    const keys = results.map(resultKey);
    if (!keys.includes(state.sel)) {
        state.sel = keys[0] || null;
        writeUrl(false);
    }

    showResultsList(results);
    showPreview();
    showPagination(pagination);
}

/**
 * The key of a result in the "sel" URL parameter: "<type>:<id>".
 */
const resultKey = (result) => `${result.type}:${result.id}`;

/**
 * Builds the results list.
 */
const showResultsList = (results) => {
    const list = document.getElementById("results-list");
    const template = document.getElementById("result-row-template");
    list.replaceChildren();

    for (const result of results) {
        const row = template.content.cloneNode(true);
        const button = row.querySelector(".result-row");
        button.dataset.key = resultKey(result);
        button.dataset.type = result.type;
        button.setAttribute("aria-pressed", String(resultKey(result) === state.sel));
        row.querySelector(".type-icon i").classList.add(TYPE_INFO[result.type].icon);
        row.querySelector(".result-title").textContent = result.title;
        row.querySelector(".result-subline").textContent = result.subline || "";
        list.appendChild(row);
    }
}

/**
 * Creates an element with the given class names and text.
 */
const makeElement = (tag, classNames = [], text = null) => {
    const elt = document.createElement(tag);
    elt.classList.add(...classNames);
    if (text !== null) elt.textContent = text;
    return elt;
}

/**
 * A titled block for the type-specific part of the preview.
 */
const makeSection = (heading) => {
    const section = makeElement("div", ["preview-section"]);
    section.appendChild(makeElement("h4", ["preview-section-title"], heading));
    return section;
}

/**
 * Makes a link open in a new tab, and says so: an icon for sighted users and hidden text for
 * screen readers (WCAG 3.2.5).
 *
 * @param {HTMLAnchorElement} link
 */
const openInNewTab = (link) => {
    link.target = "_blank";
    link.rel = "noopener";
    const icon = makeElement("i", ["mdi", "mdi-open-in-new", "new-tab-icon"]);
    icon.setAttribute("aria-hidden", "true");
    link.append(icon, makeElement("span", ["is-sr-only"], " (opens in a new tab)"));
}

/**
 * A list of links, e.g. the datasets in a collection.
 *
 * @param {Array<{title: string, url: string}>} items
 * @param {number} total - How many there are in all; any beyond the list are mentioned.
 */
const makeLinkList = (items, total = items.length) => {
    const list = makeElement("ul", ["preview-link-list"]);
    for (const item of items) {
        const li = document.createElement("li");
        const link = makeElement("a", [], item.title);
        link.href = item.url;
        openInNewTab(link);
        li.appendChild(link);
        list.appendChild(li);
    }
    if (total > items.length) {
        list.appendChild(makeElement("li", ["has-text-grey"], `and ${total - items.length} more`));
    }
    return list;
}

/**
 * The type-specific part of the preview.
 */
const makeExtra = (result) => {
    const extra = result.extra || {};
    const parts = [];

    if (result.type === "genelist" && extra.genes?.length) {
        const section = makeSection("Genes");
        const chips = makeElement("div", ["chips"]);
        for (const gene of extra.genes) {
            chips.appendChild(makeElement("span", ["chip", "gene-chip"], gene));
        }
        if (extra.gene_count > extra.genes.length) {
            chips.appendChild(makeElement("span", ["more-label"], `+ ${extra.gene_count - extra.genes.length} more`));
        }
        section.appendChild(chips);
        parts.push(section);
    }
    // Other types' chips are their organisms, which the details above already list

    if (result.type === "collection" && extra.datasets?.length) {
        const section = makeSection("Datasets in this collection");
        section.appendChild(makeLinkList(extra.datasets, extra.dataset_count));
        parts.push(section);
    }

    if (result.type === "study") {
        if (extra.thumbnail_url) {
            const image = makeElement("img", ["preview-thumbnail"]);
            image.src = extra.thumbnail_url;
            image.alt = "";     // decorative: the title is right above it
            image.loading = "lazy";
            parts.push(image);
        }
        if (extra.collections?.length) {
            const section = makeSection("Linked collections");
            section.appendChild(makeLinkList(extra.collections.map((c) => ({title: c.label, url: c.url}))));
            parts.push(section);
        }
    }

    if (result.type === "projection" && extra.runs?.length) {
        const section = makeSection("Saved runs");
        section.appendChild(makeLinkList(extra.runs.map((run) => ({title: run.label, url: run.url}))));
        parts.push(section);
    }

    if (result.type === "dataset" && extra.preview_image_url) {
        const section = makeSection("Default display");
        const image = makeElement("img", ["preview-thumbnail"]);
        image.src = extra.preview_image_url;
        image.alt = `Default display of ${result.title}`;
        image.loading = "lazy";
        section.appendChild(image);
        parts.push(section);
    }

    return parts;
}

/**
 * Shows the selected result in the preview pane.
 */
const showPreview = () => {
    const preview = document.getElementById("result-preview");
    const result = lastResponse?.results.find((r) => resultKey(r) === state.sel);
    preview.replaceChildren();
    if (!result) return;

    const content = document.getElementById("result-preview-template").content.cloneNode(true);
    const info = TYPE_INFO[result.type];
    preview.dataset.type = result.type;

    content.querySelector(".preview-type i").classList.add(info.icon);
    content.querySelector(".preview-type-label").textContent = info.label;
    content.querySelector(".access-badge").textContent = result.access || "";

    content.querySelector(".preview-title").textContent = result.title;
    const desc = content.querySelector(".preview-desc");
    desc.textContent = result.desc || "";
    desc.classList.toggle("is-hidden", !result.desc);

    // Actions: the main one, an optional second one, and the share link
    const {primary, secondary, share_url: shareUrl} = result.actions || {};
    for (const [selector, action] of [[".preview-primary", primary], [".preview-secondary", secondary]]) {
        const link = content.querySelector(selector);
        if (action) {
            link.textContent = action.label;
            link.href = action.url;
            openInNewTab(link);
        } else {
            link.remove();
        }
    }
    const copyButton = content.querySelector(".preview-copy");
    if (shareUrl) {
        copyButton.addEventListener("click", async () => {
            // The result URLs are relative to the site root (see <base href="/">)
            const absoluteUrl = new URL(shareUrl, document.baseURI).href;
            if (await copyToClipboard(absoluteUrl)) {
                createToast("Link copied to the clipboard", "is-success");
            } else {
                createToast(`Could not copy the link. Here it is: ${absoluteUrl}`);
            }
        });
    } else {
        copyButton.remove();
    }

    const meta = content.querySelector(".preview-meta");
    for (const [name, value] of result.meta || []) {
        const row = makeElement("div", ["preview-meta-row"]);
        row.appendChild(makeElement("dt", [], name));
        row.appendChild(makeElement("dd", [], value));
        meta.appendChild(row);
    }
    meta.classList.toggle("is-hidden", !(result.meta || []).length);

    content.querySelector(".preview-extra").append(...makeExtra(result));

    const links = content.querySelector(".preview-links");
    for (const item of result.links || []) {
        const link = makeElement("a", [], item.label);
        link.href = item.url;
        openInNewTab(link);
        links.appendChild(link);
    }
    links.classList.toggle("is-hidden", !(result.links || []).length);

    preview.appendChild(content);
}

/**
 * Fills the pagers above and below the results.
 */
const showPagination = (pagination) => {
    const {current_page: current, total_results: total, limit} = pagination;
    const first = (current - 1) * limit + 1;
    const last = Math.min(current * limit, total);
    for (const range of document.getElementsByClassName("results-range")) {
        range.textContent = total ? `Showing ${first}–${last} of ${total}` : "";
    }
    for (const nav of document.getElementsByClassName("results-pagination")) {
        fillPageNav(nav, pagination);
    }
}

/**
 * Builds one set of page links: previous, a window of page numbers, next.
 *
 * @param {HTMLElement} nav
 */
const fillPageNav = (nav, pagination) => {
    const {current_page: current, total_pages: totalPages} = pagination;
    nav.replaceChildren();
    if (totalPages <= 1) return;

    const pageButton = (label, page, classNames, ariaLabel) => {
        const button = makeElement("button", ["button", "is-small", ...classNames], label);
        button.type = "button";
        button.dataset.page = page;
        if (ariaLabel) button.setAttribute("aria-label", ariaLabel);
        return button;
    }

    const previous = pageButton("Previous", current - 1, ["pagination-previous"], "Previous page");
    previous.disabled = !pagination.prev_page;
    const next = pageButton("Next", current + 1, ["pagination-next"], "Next page");
    next.disabled = !pagination.next_page;

    const list = makeElement("ul", ["pagination-list"]);
    // The first and last pages, and up to two on either side of the current one
    const pages = [...new Set([1, current - 2, current - 1, current, current + 1, current + 2, totalPages])]
        .filter((page) => page >= 1 && page <= totalPages)
        .sort((a, b) => a - b);
    let previousPage = 0;
    for (const page of pages) {
        if (page - previousPage > 1) {
            const ellipsisItem = document.createElement("li");
            ellipsisItem.appendChild(makeElement("span", ["pagination-ellipsis"], "…"));
            list.appendChild(ellipsisItem);
        }
        const item = document.createElement("li");
        const isCurrent = page === current;
        const button = pageButton(String(page), page, ["pagination-link", ...(isCurrent ? ["is-current"] : [])],
            `Page ${page}`);
        if (isCurrent) button.setAttribute("aria-current", "page");
        item.appendChild(button);
        list.appendChild(item);
        previousPage = page;
    }

    nav.append(previous, next, list);
}

/**
 * Starts a new search from the controls: back to page 1, nothing selected.
 */
const searchFromControls = async () => {
    state.page = 1;
    state.sel = null;
    writeUrl(true);
    await runSearch();
}

/**
 * Fills the organism filter from the organism list.
 */
const loadOrganisms = async () => {
    try {
        const data = await apiCallsMixin.fetchOrganismList();
        const select = document.getElementById("search-organism");
        for (const organism of data.organisms) {
            const option = document.createElement("option");
            option.value = String(organism.id);
            option.textContent = organism.label;
            select.appendChild(option);
        }
    } catch (error) {
        // The page still works without the organism filter
        logErrorInConsole(error);
    }
}

const handlePageSpecificLoginUIUpdates = async (event) => {
    document.getElementById("page-header-label").textContent = "Search gEAR";

    // Logged-out searches see only public items, so the ownership filter would do nothing
    const loggedIn = Boolean(getCurrentUser()?.session_id);
    document.getElementById("search-ownership-field").classList.toggle("is-hidden", !loggedIn);

    await loadOrganisms();
    state = stateFromUrl();
    showStateInControls();
    await runSearch();
}

registerPageSpecificLoginUIUpdates(handlePageSpecificLoginUIUpdates);

await initCommonUI();

// ------------------------------------------------------------------ event listeners

document.getElementById("search-form").addEventListener("submit", async (event) => {
    event.preventDefault();
    state.q = document.getElementById("search-q").value.trim();
    await searchFromControls();
});

// Type buttons: All clears the selection; the others add to it. Picking every type is All.
document.getElementById("type-filter").addEventListener("click", async (event) => {
    const button = event.target.closest(".type-button");
    if (!button) return;

    const type = button.dataset.type;
    if (type === "all") {
        state.types = [];
    } else if (state.types.includes(type)) {
        state.types = state.types.filter((t) => t !== type);
    } else {
        state.types = [...state.types, type];
    }
    if (state.types.length === TYPES.length) {
        state.types = [];
    }
    showStateInControls();
    await searchFromControls();
});

for (const [id, name] of [["search-organism", "organism_ids"], ["search-ownership", "ownership"],
                          ["search-date-added", "date_added"], ["search-sort", "sort"]]) {
    document.getElementById(id).addEventListener("change", async (event) => {
        state[name] = event.target.value;
        await searchFromControls();
    });
}

document.getElementById("results-list").addEventListener("click", (event) => {
    const row = event.target.closest(".result-row");
    if (!row) return;

    state.sel = row.dataset.key;
    writeUrl(false);
    for (const other of document.getElementsByClassName("result-row")) {
        other.setAttribute("aria-pressed", String(other === row));
    }
    showPreview();

    // On narrow screens the preview is below the list, out of sight; bring it into view
    const preview = document.getElementById("result-preview");
    if (preview.offsetTop > document.getElementById("results-list").offsetTop) {
        preview.scrollIntoView({behavior: "smooth", block: "start"});
    }
});

for (const nav of document.getElementsByClassName("results-pagination")) {
    nav.addEventListener("click", async (event) => {
        const button = event.target.closest("button[data-page]");
        if (!button || button.disabled) return;

        state.page = Number(button.dataset.page);
        state.sel = null;
        writeUrl(true);
        await runSearch();
        document.getElementById("results-summary").scrollIntoView({behavior: "smooth", block: "start"});
    });
}

document.getElementById("clear-search").addEventListener("click", async () => {
    state.q = "";
    state.types = [];
    showStateInControls();
    await searchFromControls();
});

// Back and Forward move between searches
window.addEventListener("popstate", async () => {
    state = stateFromUrl();
    showStateInControls();
    await runSearch();
});
