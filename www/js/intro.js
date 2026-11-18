/*
 * intro.js - Intro page for new gEAR visitors.
 *
 * All statistics and "What's new" content come from one pre-generated file,
 * stats_archive/intro_stats.json (see bin/generate_intro_stats.py), so this page
 * never queries the database. If the file is missing, those sections hide and
 * the rest of the page (search, tools, guides) still works.
 */

import { loadFeaturedStudies } from './featured_studies.js';

const STATS_URL = '/stats_archive/intro_stats.json';

// ---------------------------------------------------------------------------
// Visitor state (mirrors intro.html so index.html's first-visit redirect works)
// ---------------------------------------------------------------------------
document.cookie = 'gear_returning_visitor=1; path=/; max-age=31536000; SameSite=Lax';

if (/(^|;\s*)gear_session_id=/.test(document.cookie)) {
    for (const id of ['header-create-account', 'cta-create-account']) {
        document.getElementById(id)?.classList.add('is-hidden');
    }
    document.querySelector('.manage-callout')?.classList.add('is-hidden');
    // Logged-in visitors get a button to their workspace in place of the Log in link
    const workspaceLink = document.getElementById('header-login');
    workspaceLink.textContent = 'Go to my workspace';
    workspaceLink.classList.add('button', 'is-small', 'is-rounded', 'intro-nav-cta');
}

// ---------------------------------------------------------------------------
// Site menu: disclosure dropdowns (click/tap or Enter to open; also on hover with a mouse)
// ---------------------------------------------------------------------------
const navEl = document.getElementById('intro-nav');
const navToggle = navEl.querySelector('.intro-nav-toggle');
const navTriggers = [...navEl.querySelectorAll('.intro-nav-trigger')];
const canHover = window.matchMedia('(hover: hover) and (pointer: fine)');

const setDropdown = (trigger, open) => {
    trigger.setAttribute('aria-expanded', String(open));
    trigger.closest('.intro-nav-item').classList.toggle('is-open', open);
};
const closeDropdowns = (except = null) => {
    for (const t of navTriggers) if (t !== except) setDropdown(t, false);
};

for (const trigger of navTriggers) {
    const item = trigger.closest('.intro-nav-item');
    trigger.addEventListener('click', () => {
        const open = trigger.getAttribute('aria-expanded') !== 'true';
        closeDropdowns(trigger);
        setDropdown(trigger, open);
    });
    item.addEventListener('mouseenter', () => {
        if (!canHover.matches || navEl.classList.contains('is-collapsed')) return;
        closeDropdowns(trigger);
        setDropdown(trigger, true);
    });
    item.addEventListener('mouseleave', () => {
        if (canHover.matches && !navEl.classList.contains('is-collapsed')) setDropdown(trigger, false);
    });
    // Close once focus leaves the menu entirely (e.g. tabbing past the last link)
    item.addEventListener('focusout', (event) => {
        if (!item.contains(event.relatedTarget)) setDropdown(trigger, false);
    });
}

// The bar is transparent over the hero and turns solid (with a shadow) once the page scrolls
const updateNavSolid = () => navEl.classList.toggle('is-scrolled', window.scrollY > 8);
window.addEventListener('scroll', updateNavSolid, { passive: true });
// Also catch a scroll position the browser restores on reload or back/forward
for (const type of ['load', 'pageshow', 'resize']) window.addEventListener(type, updateNavSolid);
updateNavSolid();

// The narrow-screen "Menu" button shows or hides the whole list
navToggle.addEventListener('click', () => {
    const open = navToggle.getAttribute('aria-expanded') !== 'true';
    navToggle.setAttribute('aria-expanded', String(open));
    navEl.classList.toggle('is-menu-open', open);
});

// Following a link (including the in-page guides link) closes everything
navEl.addEventListener('click', (event) => {
    if (!event.target.closest('a')) return;
    closeDropdowns();
    navToggle.setAttribute('aria-expanded', 'false');
    navEl.classList.remove('is-menu-open');
});

document.addEventListener('click', (event) => {
    if (!navEl.contains(event.target)) closeDropdowns();
});
document.addEventListener('keydown', (event) => {
    if (event.key !== 'Escape') return;
    const open = navTriggers.find((t) => t.getAttribute('aria-expanded') === 'true');
    if (open) {
        setDropdown(open, false);
        open.focus();
    }
});

// Narrow layout: the list collapses behind the Menu button and dropdowns open inline
const narrow = window.matchMedia('(max-width: 768px)');
const syncNavLayout = () => {
    navEl.classList.toggle('is-collapsed', narrow.matches);
    navToggle.setAttribute('aria-expanded', 'false');
    navEl.classList.remove('is-menu-open');
    closeDropdowns();
};
narrow.addEventListener('change', syncNavLayout);
syncNavLayout();

// ---------------------------------------------------------------------------
// Portal branding from site_domain_prefs.json
// (colors need nothing here: intro.css reads the --gear-* variables the theme build sets)
// ---------------------------------------------------------------------------

/**
 * Applies the portal's name, tagline, logo and example genes to the page.
 *
 * Uses domain_label, domain_short_display_label and domain_tagline, plus optionally
 * intro_example_genes: a list of searches for the "Try:" chips, each one gene or a
 * comma-separated set (e.g. ["Sox2", "Myo7a,Pou4f3,Gfi1"]). Without it the chips fall
 * back to default_demo_gene, and are hidden if neither is set.
 *
 * @param {Object} prefs - Result of getDomainPreferences()
 */
export const applyDomainBranding = (prefs) => {
    const label = prefs.domain_short_display_label || 'gEAR';

    document.title = document.title.replace('gEAR', label);
    for (const el of document.querySelectorAll('.js-domain-label')) {
        el.textContent = label;
    }

    // The tagline (e.g. "gene Expression Analysis Resource") is a small line above the hero question
    if (prefs.domain_tagline?.trim()) {
        const tagline = document.getElementById('hero-tagline');
        tagline.textContent = prefs.domain_tagline.trim();
        tagline.hidden = false;
    }

    // Same dark-background logo the portal sidebar uses; show the name as text if a portal has none
    const logo = document.getElementById('intro-logo');
    logo.alt = `${label} home`;
    logo.addEventListener('error', () => {
        logo.hidden = true;
        document.getElementById('intro-brand-text').hidden = false;
    }, { once: true });
    logo.src = `/img/by_domain/${encodeURIComponent(prefs.domain_label)}/logo-main-normal.png`;

    renderExampleGenes(prefs);
};

const renderExampleGenes = (prefs) => {
    let examples = Array.isArray(prefs.intro_example_genes) ? prefs.intro_example_genes : [];
    if (!examples.length && prefs.default_demo_gene) examples = [prefs.default_demo_gene];

    const searches = examples
        .map((example) => String(example).split(/[\s,;]+/).filter(Boolean))
        .filter((genes) => genes.length);
    if (!searches.length) return;

    const container = document.getElementById('hero-example-genes');
    for (const genes of searches) {
        const params = new URLSearchParams({ gene_symbol: genes.join(',') });
        if (genes.length > 1) params.set('is_multigene', '1');
        const link = document.createElement('a');
        link.href = `/expression.html?${params}`;
        link.textContent = genes.join(' + ');
        container.append(' ', link);
    }
    container.hidden = false;
    document.getElementById('hero-examples-sep').hidden = false;

    const firstGenes = searches.find((genes) => genes.length === 1) ?? searches[0];
    searchInput.placeholder = `Enter one or more gene symbols, e.g. ${firstGenes.join(', ')}`;
};

// ---------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------
const searchForm = document.getElementById('hero-search-form');
const searchInput = document.getElementById('hero-gene-symbol');

searchForm.addEventListener('submit', (event) => {
    event.preventDefault();
    // Accept commas, spaces or semicolons between symbols
    const genes = searchInput.value.split(/[\s,;]+/).filter(Boolean);
    if (!genes.length) return;

    // Always single-gene display; with several symbols, expression.html lists them and shows
    // the first, and the user can switch to multi-gene there
    const params = new URLSearchParams({ gene_symbol: genes.join(','), is_multigene: '0' });
    window.location.href = `/expression.html?${params}`;
});

document.getElementById('cta-search').addEventListener('click', (event) => {
    event.preventDefault();
    window.scrollTo({ top: 0, behavior: 'smooth' });
    searchInput.focus({ preventScroll: true });
});

// ---------------------------------------------------------------------------
// What's new tabs (WAI-ARIA tabs pattern, arrow-key navigable)
// ---------------------------------------------------------------------------
const tabs = [...document.querySelectorAll('.new-tab')];

const selectTab = (tab) => {
    for (const t of tabs) {
        const selected = t === tab;
        t.classList.toggle('is-active', selected);
        t.setAttribute('aria-selected', selected);
        t.tabIndex = selected ? 0 : -1;
        document.getElementById(t.getAttribute('aria-controls')).hidden = !selected;
    }
};

tabs.forEach((tab, i) => {
    tab.addEventListener('click', () => selectTab(tab));
    tab.addEventListener('keydown', (event) => {
        const step = { ArrowRight: 1, ArrowLeft: -1 }[event.key];
        if (!step) return;
        const next = tabs[(i + step + tabs.length) % tabs.length];
        selectTab(next);
        next.focus();
    });
});

// ---------------------------------------------------------------------------
// Rendering helpers
// ---------------------------------------------------------------------------
const fmt = new Intl.NumberFormat();
const fmtDate = (iso) => iso
    ? new Date(`${iso}T00:00:00`).toLocaleDateString(undefined, { year: 'numeric', month: 'short', day: 'numeric' })
    : '';

const escapeHtml = (s) => String(s ?? '').replace(/[&<>"']/g,
    (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

// Many titles are entered in ALL CAPS; soften those for display only
const tidyTitle = (s) => {
    const letters = s.replace(/[^A-Za-z]/g, '');
    if (letters.length > 8 && letters === letters.toUpperCase()) {
        return s.toLowerCase().replace(/^\w/, (c) => c.toUpperCase())
            .replace(/\b(snrna|scrna|rna|dna|atac)(-seq)?\b/gi, (m) => m.replace(/rna|dna|atac/i, (x) => x.toUpperCase()));
    }
    return s;
};

const tooltip = document.getElementById('chart-tooltip');

const attachTooltip = (el, html) => {
    const show = (event) => {
        tooltip.innerHTML = html;
        tooltip.hidden = false;
        const rect = el.getBoundingClientRect();
        const x = event?.clientX ?? rect.left + rect.width / 2;
        const tipRect = tooltip.getBoundingClientRect();
        let left = x - tipRect.width / 2;
        left = Math.max(8, Math.min(left, window.innerWidth - tipRect.width - 8));
        tooltip.style.left = `${left + window.scrollX}px`;
        tooltip.style.top = `${rect.top + window.scrollY - tipRect.height - 8}px`;
    };
    const hide = () => { tooltip.hidden = true; };
    el.addEventListener('mouseenter', show);
    el.addEventListener('mousemove', show);
    el.addEventListener('mouseleave', hide);
    el.addEventListener('focus', () => show());
    el.addEventListener('blur', hide);
};

// Horizontal bars: label | bar | value. One series, so no legend; values are direct-labeled.
const renderHBars = (container, rows, total) => {
    const max = Math.max(...rows.map((r) => r.count));
    container.innerHTML = '';
    for (const row of rows) {
        const pct = total ? Math.round((row.count / total) * 100) : 0;
        const item = document.createElement('div');
        item.className = 'hbar-row';
        item.tabIndex = 0;
        item.innerHTML = `
            <span class="hbar-label">${escapeHtml(row.label)}</span>
            <span class="hbar-track"><span class="hbar-fill" style="--w:${(row.count / max) * 100}%"></span></span>
            <span class="hbar-value">${fmt.format(row.count)}</span>`;
        attachTooltip(item, `<strong>${escapeHtml(row.label)}</strong><br>${fmt.format(row.count)} datasets · ${pct}% of public`);
        container.appendChild(item);
    }
};

// Column chart of cumulative public datasets per year
const renderGrowth = (container, rows) => {
    const max = Math.max(...rows.map((r) => r.cumulative));
    container.innerHTML = '';
    const plot = document.createElement('div');
    plot.className = 'col-plot';
    rows.forEach((row, i) => {
        const col = document.createElement('div');
        col.className = 'col-item';
        col.tabIndex = 0;
        const isLast = i === rows.length - 1;
        col.innerHTML = `
            ${isLast ? `<span class="col-value">${fmt.format(row.cumulative)}</span>` : ''}
            <span class="col-fill${isLast ? ' is-current' : ''}" style="--h:${(row.cumulative / max) * 100}%"></span>`;
        attachTooltip(col, `<strong>${row.year}</strong><br>${fmt.format(row.cumulative)} total · +${fmt.format(row.count)} that year`);
        plot.appendChild(col);
    });
    const axis = document.createElement('div');
    axis.className = 'col-axis';
    // Label only the first, middle and last year so labels never collide
    const labelIdx = new Set([0, Math.floor((rows.length - 1) / 2), rows.length - 1]);
    rows.forEach((row, i) => {
        const tick = document.createElement('span');
        tick.textContent = labelIdx.has(i) ? row.year : '';
        axis.appendChild(tick);
    });
    container.append(plot, axis);
};

const renderCards = (panel, items, kind) => {
    if (!items?.length) {
        panel.innerHTML = '<p class="muted">Nothing to show yet. Check back soon.</p>';
        return;
    }
    const icon = { datasets: 'mdi-database-outline', collections: 'mdi-folder-open-outline' }[kind];
    panel.innerHTML = `<div class="new-grid">${items.map((item) => {
        let meta = [];
        if (kind === 'datasets') meta = [item.organism, item.dtype, fmtDate(item.date)];
        if (kind === 'collections') meta = [`${item.dataset_count} datasets`, item.organisms];
        return `
            <a class="new-card" href="${escapeHtml(item.url)}">
                <span class="new-card-kind"><i class="mdi ${icon}" aria-hidden="true"></i></span>
                <span class="new-card-title">${escapeHtml(tidyTitle(item.title))}</span>
                ${item.summary ? `<span class="new-card-summary">${escapeHtml(item.summary)}</span>` : ''}
                <span class="new-card-meta">${meta.filter(Boolean).map(escapeHtml).join(' · ')}</span>
            </a>`;
    }).join('')}</div>`;
};

// ---------------------------------------------------------------------------
// Load stats
// ---------------------------------------------------------------------------
const loadStats = async () => {
    let data;
    try {
        const response = await fetch(STATS_URL, { cache: 'no-cache' });
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        data = await response.json();
    } catch (error) {
        console.warn('Intro stats unavailable; hiding stats and What\'s new sections.', error);
        document.getElementById('by-the-numbers').hidden = true;
        document.getElementById('whats-new').hidden = true;
        document.querySelector('.hero-ticker').hidden = true;
        return;
    }

    for (const el of document.querySelectorAll('.js-stat')) {
        const value = data.totals?.[el.dataset.stat];
        if (value != null) el.textContent = fmt.format(value);
    }

    if (data.generated) {
        document.getElementById('stats-generated').textContent = `Updated ${fmtDate(data.generated.slice(0, 10))}.`;
    }

    const total = data.totals?.public_datasets;
    renderHBars(document.getElementById('chart-organism'), data.by_organism ?? [], total);
    renderHBars(document.getElementById('chart-dtype'), data.by_dtype ?? [], total);
    renderGrowth(document.getElementById('chart-growth'), data.by_year ?? []);

    renderCards(document.getElementById('new-panel-datasets'), data.recent_datasets, 'datasets');
    renderCards(document.getElementById('new-panel-collections'), data.recent_collections, 'collections');

    // Reuse the data already fetched rather than requesting the file a second time
    loadFeaturedStudies(data.landing_pages ?? []);
};

loadStats();
