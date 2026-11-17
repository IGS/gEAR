/*
 * landing_portal_check.js - Shared by every page under www/landing/.
 *
 * All portals (gEAR, NeMO Analytics, cancer gEAR, ...) deploy the same www/landing/ folder, so each
 * landing page exists on every portal even though its datasets only live on some of them. Each page's
 * landing.json lists the portals it belongs to. When the page is opened on another portal, this shows
 * a banner at the top of the content pointing to the page's canonical_url.
 */
(async () => {
    try {
        const [prefs, manifest] = await Promise.all([
            fetch('/site_domain_prefs.json').then((r) => r.json()),
            fetch(new URL('landing.json', window.location.href)).then((r) => r.json()),
        ]);

        const portals = manifest.portals ?? [];
        if (portals.includes('*') || portals.includes(prefs.domain_label) || !manifest.canonical_url) {
            return;
        }

        const target = new URL(manifest.canonical_url);
        // Already on the canonical host (e.g. a misconfigured domain_label); nothing to point to
        if (target.host === window.location.host) return;

        const banner = document.createElement('div');
        banner.className = 'notification is-warning is-light m-3';
        banner.setAttribute('role', 'status');

        const text = document.createElement('span');
        text.textContent = `This study's data is hosted on ${target.host}, so some links on this page may not work here. `;

        const link = document.createElement('a');
        link.href = target.href + window.location.hash;
        link.textContent = `View it on ${target.host} →`;
        link.className = 'has-text-weight-semibold';

        banner.append(text, link);

        const content = document.getElementById('content-c') ?? document.body;
        content.prepend(banner);
    } catch (error) {
        // Missing landing.json or prefs: leave the page as it is
        console.warn('Landing portal check skipped:', error);
    }
})();
