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