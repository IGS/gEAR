'use strict';

/*
 * landing.js - landing.html, the full list of study landing pages for this portal.
 *
 * Pages come from landing_pages in stats_archive/intro_stats.json (built by
 * bin/generate_intro_stats.py from each www/landing/<name>/landing.json), narrowed to this
 * portal the same way as the "Featured studies" scroller on index.html and intro.html.
 * Unlike the scroller, every page is listed, not only the featured ones.
 */

import { initCommonUI, logErrorInConsole } from "./common.v2.js";
import { landingPagesForThisPortal } from "./featured_studies.js";

const STATS_URL = '/stats_archive/intro_stats.json';

const list = document.getElementById('studies-list');
const filterInput = document.getElementById('studies-filter');
const sortSelect = document.getElementById('studies-sort');
const countLabel = document.getElementById('studies-count');
const noMatch = document.getElementById('studies-no-match');

let pages = [];

const fmtDate = (iso) => iso
    ? new Date(`${iso}T00:00:00`).toLocaleDateString(undefined, { year: 'numeric', month: 'short', day: 'numeric' })
    : '';

const el = (tag, className, text) => {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text) node.textContent = text;
    return node;
};
const icon = (name) => el('i', `mdi ${name}`);

const placeholderMedia = (media) => {
    media.classList.add('is-placeholder');
    media.replaceChildren(icon('mdi-book-open-page-variant-outline'));
};

/** One study as a list item: figure on the left, details and links on the right. */
const buildCard = (page) => {
    const card = el('li', 'study-card');

    const media = el('div', 'study-media');
    media.setAttribute('aria-hidden', 'true');
    if (page.thumbnail) {
        const img = el('img');
        img.src = page.thumbnail;
        img.alt = '';
        img.loading = 'lazy';
        img.decoding = 'async';
        img.addEventListener('error', () => placeholderMedia(media));
        media.append(img);
    } else {
        placeholderMedia(media);
    }

    const body = el('div', 'study-body');

    const meta = el('p', 'study-meta');
    if (page.featured) {
        const badge = el('span', 'study-badge', ' Featured');
        badge.prepend(icon('mdi-star-four-points'));
        meta.append(badge);
    }
    if (page.date) {
        const added = el('span', 'study-date', `Added ${fmtDate(page.date)}`);
        meta.append(added);
    }

    const title = el('h3', 'study-title');
    const link = el('a', null, page.title);
    link.href = page.url;
    title.append(link);

    body.append(meta, title);
    if (page.summary) body.append(el('p', 'study-summary', page.summary));

    const actions = el('div', 'study-actions');
    const explore = el('span', 'study-cta', 'Explore this study ');
    explore.append(icon('mdi-arrow-right'));
    explore.setAttribute('aria-hidden', 'true');
    actions.append(explore);
    if (page.pmid) {
        const pubmed = el('a', 'study-pubmed', `PubMed ${page.pmid} `);
        pubmed.href = `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(page.pmid)}/`;
        pubmed.target = '_blank';
        pubmed.rel = 'noopener';
        pubmed.setAttribute('aria-label', `PubMed record ${page.pmid} (opens in a new tab)`);
        pubmed.append(icon('mdi-open-in-new'));
        actions.append(pubmed);
    }
    body.append(actions);

    card.append(media, body);
    // Searchable text, lowercased once
    card.dataset.search = `${page.title} ${page.summary ?? ''}`.toLowerCase();
    return card;
};

const SORTS = {
    newest: (a, b) => (b.date ?? '').localeCompare(a.date ?? '') || a.title.localeCompare(b.title),
    oldest: (a, b) => (a.date ?? '').localeCompare(b.date ?? '') || a.title.localeCompare(b.title),
    title: (a, b) => a.title.localeCompare(b.title, undefined, { sensitivity: 'base' }),
};

const render = () => {
    const term = filterInput.value.trim().toLowerCase();
    const sorted = [...pages].sort(SORTS[sortSelect.value] ?? SORTS.newest);

    let shown = 0;
    list.replaceChildren(...sorted.map((page) => {
        const card = page.card;
        const match = !term || card.dataset.search.includes(term);
        card.hidden = !match;
        if (match) shown++;
        return card;
    }));

    const noun = pages.length === 1 ? 'study' : 'studies';
    countLabel.textContent = term ? `${shown} of ${pages.length} ${noun}` : `${pages.length} ${noun}`;
    noMatch.hidden = shown > 0;
    document.getElementById('studies-no-match-term').textContent = filterInput.value.trim();
};

const showUnavailable = (message) => {
    list.replaceChildren();
    list.hidden = true;
    if (message) document.getElementById('studies-unavailable-text').textContent = message;
    document.getElementById('studies-unavailable').hidden = false;
};

const loadStudies = async () => {
    try {
        const response = await fetch(STATS_URL, { cache: 'no-cache' });
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        pages = await landingPagesForThisPortal((await response.json()).landing_pages ?? []);
    } catch (error) {
        logErrorInConsole(`Study list unavailable: ${error}`);
        showUnavailable("The study list couldn't be loaded right now. Please try again later.");
        return;
    } finally {
        list.removeAttribute('aria-busy');
    }

    if (!pages.length) {
        showUnavailable();
        return;
    }

    for (const page of pages) page.card = buildCard(page);

    filterInput.disabled = false;
    sortSelect.disabled = false;
    filterInput.addEventListener('input', render);
    sortSelect.addEventListener('change', render);
    document.getElementById('studies-clear-filter').addEventListener('click', () => {
        filterInput.value = '';
        render();
        filterInput.focus();
    });

    render();
};

document.getElementById('page-header-label').textContent = 'Studies';

await initCommonUI();
loadStudies();
