"use strict";

// allows gEAR citation to be generated
// eventually other citation formats may be added here as well

export class Citation {
    static getDOI(shareId) {
        return `https://umgear.org/p?id=d.${shareId}`;
    }

    /**
     * Builds the citation/watermark lines stamped below downloaded plots.
     * Mirrors build_stamp_lines in lib/gear/plot_stamp.py, which stamps the server-rendered downloads.
     *
     * @param {Object} dataset - Dataset being plotted (title, pubmed_id, share_id).
     * @param {Object} [sitePrefs={}] - Site domain preferences (domain_url, domain_short_display_label).
     * @returns {string[]} Lines of plain text, top to bottom.
     */
    static plotStampLines(dataset, sitePrefs={}) {
        const maxTitleLength = 100;
        const siteUrl = (sitePrefs.domain_url || window.location.origin).replace(/\/+$/, "");
        const siteLabel = sitePrefs.domain_short_display_label || "gEAR";

        let title = (dataset.title || "").trim();
        if (title.length > maxTitleLength) {
            title = `${title.slice(0, maxTitleLength - 1).trimEnd()}…`;
        }

        // Point back to the dataset permalink when there is no publication to cite
        const pubmedId = String(dataset.pubmed_id ?? "").trim();
        let siteLine = `Made with ${siteLabel} · ${siteUrl}`;
        if (/^\d+$/.test(pubmedId)) {
            siteLine = `PMID ${pubmedId} · Made with ${siteLabel} · ${siteUrl.replace(/^https?:\/\//, "")}`;
        } else if (dataset.share_id) {
            siteLine = `Made with ${siteLabel} · ${siteUrl}/p?id=d.${dataset.share_id}`;
        }

        return [title, siteLine].filter(Boolean);
    }

    static gEAR(authors, year, title, shareId, accessDate, license) {
        if (authors.length > 2) {
            authors = `${authors[0]} et al.`;
        } else if (authors.length === 2) {
            authors = `${authors[0]} and ${authors[1]}`;
        } else {
            authors = authors[0];
        }

        const accessTimeStamp = accessDate.toLocaleDateString('en-US', { day: 'numeric', month: 'short', year: 'numeric' });

        const licenseToUse = license ? ` Licensed under ${license}.` : "";

        return {
            orig: `${authors} (${year}). ${title}. Available from ${Citation.getDOI(shareId)} (Accessed ${accessTimeStamp}).${licenseToUse}`,
            format: `${authors} (${year}). <i>${title}</i>. Available from ${Citation.getDOI(shareId)} (Accessed ${accessTimeStamp}).${licenseToUse}`
        }
    }

    static APA(authors, year, title, shareId, accessDate, license) {
        const accessTimestamp = accessDate.toLocaleDateString('en-US', { day: 'numeric', month: 'short', year: 'numeric' });

        const licenseToUse = license ? ` Licensed under ${license}.` : "";

        // If there are no authors, we should start with the title
        if (authors === null) {
            return {
                orig: `${title} (${year}). [Data set]. gEAR Portal. Retrieved ${accessTimestamp}, from ${Citation.getDOI(shareId)}.${licenseToUse}`,
                format: `<i>${title}</i> (${year}). [Data set]. gEAR Portal. Retrieved ${accessTimestamp}, from ${Citation.getDOI(shareId)}.${licenseToUse}`
            };
        }

        // Convert authors to "Last, F. M." format
        authors = authors.map(author => {
            const names = author.split(' ').map(s => s.trim());
            const lastName = names.pop();

            const initials = names.map(n => n[0].toUpperCase() + '.').join(' ');
            return `${lastName}, ${initials}`;
        });
        if (authors.length === 1) {
            authors = authors[0];
        } else if (authors.length <= 20) {
            authors = `${authors.slice(0, -1).join(', ')} & ${authors.slice(-1)}`;
        } else {
            authors = `${authors.slice(0, 19).join(', ')}, ... & ${authors.slice(-1)}`;
        }

        return {
            orig: `${authors} (${year}). ${title} [Data set]. gEAR Portal. Retrieved ${accessTimestamp}, from ${Citation.getDOI(shareId)}.${licenseToUse}`,
            format: `${authors} (${year}). <i>${title}</i> [Data set]. gEAR Portal. Retrieved ${accessTimestamp}, from ${Citation.getDOI(shareId)}.${licenseToUse}`
        };
    }
}