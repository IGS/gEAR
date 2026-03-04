"use strict";

// allows gEAR citation to be generated
// eventually other citation formats may be added here as well

export class Citation {
    static getDOI(shareId) {
        return `https://umgear.org/p?id=d.${shareId}`;
    }

    static gEAR(authors, year, title, shareId, accessDate, license) {
        if (authors.length > 2) {
            authors = `${authors[0]} et al.`;
        } else if (authors.length === 2) {
            authors = `${authors[0]} and ${authors[1]}`;
        } else {
            authors = authors[0];
        }

        accessDate = accessDate.toLocaleDateString('en-US', { day: 'numeric', month: 'short', year: 'numeric' });

        if (license)
            license = ` Licensed under ${license}.`;

        return {
            orig: `${authors} (${year}). ${title}. Available from ${Citation.getDOI(shareId)} (Accessed ${accessDate}).${license || ""}`,
            format: `${authors} (${year}). <i>${title}</i>. Available from ${Citation.getDOI(shareId)} (Accessed ${accessDate}).${license || ""}`
        }
    }
}