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

    static APA(authors, year, title, shareId, accessDate, license) {
        authors = authors.map(author => {
            const names = author.split(' ').map(s => s.trim());
            const lastName = names.pop();

            const initials = names.map(n => n[0].toUpperCase() + '.').join(' ');
            return `${lastName}, ${initials}`;
        });
        if (authors.length === 1) {
            authors = authors[0];
        } else if (authors.length <= 20) {
            authors = authors.slice(0, -1).join(', ') + ' & ' + authors.slice(-1);
        } else {
            authors = authors.slice(0, 19).join(', ') + ', ... & ' + authors.slice(-1);
        }

        accessDate = accessDate.toLocaleDateString('en-US', { day: 'numeric', month: 'short', year: 'numeric' });

        if (license)
            license = ` Licensed under ${license}.`;

        return {
            orig: `${authors} (${year}). ${title} [Data set]. gEAR Portal. Retrieved ${accessDate}, from ${Citation.getDOI(shareId)}.${license || ""}`,
            format: `${authors} (${year}). <i>${title}</i> [Data set]. gEAR Portal. Retrieved ${accessDate}, from ${Citation.getDOI(shareId)}.${license || ""}`
        };
    }
}