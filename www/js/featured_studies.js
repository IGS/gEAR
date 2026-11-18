/*
 * featured_studies.js - "Featured studies" side-scroller, shared by index.html and intro.html.
 *
 * Expects this markup on the page (styles in css/featured_studies.css):
 *
 *   <section id="featured-studies" hidden>
 *     ... .featured-nav with .featured-prev / .featured-next buttons ...
 *     <div class="featured-viewport"><ul id="featured-track" class="featured-track"></ul></div>
 *   </section>
 *
 * Cards come from the landing pages in stats_archive/intro_stats.json (built by
 * bin/generate_intro_stats.py), which lists every portal's pages. Only the featured ones whose
 * "portals" include domain_label from site_domain_prefs.json are shown, so switching portals
 * doesn't require regenerating the stats file.
 */

import { logErrorInConsole } from "./common.v2.js";

const STATS_URL = '/stats_archive/intro_stats.json';
const PREFS_URL = '/site_domain_prefs.json';

/** True if a landing page's "portals" list includes this portal's domain_label, or "*". */
const isForPortal = (page, domainLabel) => {
    const portals = page.portals ?? [];
    return portals.includes('*') || portals.includes(domainLabel);
};

/**
 * Narrows intro_stats.json's landing_pages (every portal's) to this portal's, using domain_label
 * from site_domain_prefs.json. Throws if the prefs file can't be read. Also used by intro.js to
 * count them.
 */
export const landingPagesForThisPortal = async (pages) => {
    const domainLabel = (await (await fetch(PREFS_URL, { cache: 'no-cache' })).json()).domain_label;
    return pages.filter((page) => isForPortal(page, domainLabel));
};

/**
 * Fills the side-scroller and un-hides the section. It stays hidden if the stats file is
 * missing or no page is marked featured.
 *
 * @param {Array|null} landingPages - landing_pages from intro_stats.json if the caller already
 *     fetched it; otherwise the file is fetched here.
 */
export const loadFeaturedStudies = async (landingPages = null) => {
    const section = document.getElementById('featured-studies');
    const track = document.getElementById('featured-track');
    if (!section || !track) return;

    let pages = landingPages;
    if (!pages) {
        try {
            const response = await fetch(STATS_URL, { cache: 'no-cache' });
            if (!response.ok) throw new Error(`HTTP ${response.status}`);
            pages = (await response.json()).landing_pages ?? [];
        } catch (error) {
            logErrorInConsole(`Featured studies unavailable: ${error}`);
            return;
        }
    }
    // intro_stats.json lists every portal's pages; keep this portal's
    try {
        pages = (await landingPagesForThisPortal(pages)).filter((page) => page.featured);
    } catch (error) {
        logErrorInConsole(`Featured studies unavailable, couldn't read domain_label: ${error}`);
        return;
    }
    if (!pages.length) return;

    const el = (tag, className, text) => {
        const node = document.createElement(tag);
        if (className) node.className = className;
        if (text) node.textContent = text;
        return node;
    };
    const icon = (name) => el('i', `mdi ${name}`);

    for (const page of pages) {
        const card = el('li', 'featured-card');

        const media = el('div', 'featured-media');
        if (page.thumbnail) {
            const img = el('img');
            img.src = page.thumbnail;
            img.alt = '';
            img.loading = 'lazy';
            img.decoding = 'async';
            // A broken image falls back to the placeholder
            img.addEventListener('error', () => {
                media.classList.add('is-placeholder');
                media.replaceChildren(icon('mdi-book-open-page-variant-outline'));
            });
            media.append(img);
        } else {
            media.classList.add('is-placeholder');
            media.append(icon('mdi-book-open-page-variant-outline'));
        }
        media.setAttribute('aria-hidden', 'true');

        const body = el('div', 'featured-body');
        const title = el('h4', 'featured-title');
        const link = el('a', null, page.title);
        link.href = page.url;
        title.append(link);
        body.append(title);
        if (page.summary) body.append(el('p', 'featured-summary', page.summary));

        const meta = el('div', 'featured-meta');
        const cta = el('span', 'featured-cta', 'Explore ');
        cta.append(icon('mdi-arrow-right'));
        cta.setAttribute('aria-hidden', 'true');
        meta.append(cta);
        if (page.pmid) {
            const pubmed = el('a', 'featured-pubmed', 'PubMed ');
            pubmed.href = `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(page.pmid)}/`;
            pubmed.target = '_blank';
            pubmed.rel = 'noopener';
            pubmed.setAttribute('aria-label', `PubMed record ${page.pmid} (opens in a new tab)`);
            pubmed.append(icon('mdi-open-in-new'));
            meta.append(pubmed);
        }
        body.append(meta);

        card.append(media, body);
        track.append(card);
    }

    section.hidden = false;

    // Arrows and edge fades only matter when the cards overflow
    const viewport = section.querySelector('.featured-viewport');
    const nav = section.querySelector('.featured-nav');
    const prev = section.querySelector('.featured-prev');
    const next = section.querySelector('.featured-next');

    const updateScrollState = () => {
        const maxScroll = track.scrollWidth - track.clientWidth;
        const overflowing = maxScroll > 1;
        nav.hidden = !overflowing;
        prev.disabled = track.scrollLeft <= 1;
        next.disabled = track.scrollLeft >= maxScroll - 1;
        viewport.classList.toggle('can-scroll-left', overflowing && !prev.disabled);
        viewport.classList.toggle('can-scroll-right', overflowing && !next.disabled);
    };

    // Page by however many whole cards are visible
    const scrollByPage = (direction) => {
        const card = track.querySelector('.featured-card');
        const step = card ? card.getBoundingClientRect().width + parseFloat(getComputedStyle(track).columnGap || 0) : track.clientWidth;
        const perPage = Math.max(1, Math.floor((track.clientWidth + 1) / step));
        track.scrollBy({ left: direction * perPage * step });
    };

    prev.addEventListener('click', () => scrollByPage(-1));
    next.addEventListener('click', () => scrollByPage(1));
    track.addEventListener('scroll', updateScrollState, { passive: true });
    new ResizeObserver(updateScrollState).observe(track);
    updateScrollState();
};
