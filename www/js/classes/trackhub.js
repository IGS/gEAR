// Hub class to represent a UCSC Track Hub
export class Hub {
    /**
     * @param {string} identifier - Unique identifier for the hub.
     * @param {string} shortLabel - Short label for the hub.
     * @param {string} longLabel - Long label for the hub.
     * @param {string} email - Contact email for the hub.
     * @param {string} genome - Genome associated with the hub. Currently only support one genome.
     */
    constructor(identifier, shortLabel, longLabel, email, genome) {
        this.identifier = identifier;
        this.shortLabel = shortLabel;
        this.longLabel = longLabel;
        this.email = email;
        this.genome = genome;
    }

    /**
     * Generates the hub.txt content for the UCSC Track Hub.
     * @returns {string} The hub.txt content.
     */
    generateHubTxt() {
        new Error("generateHubTxt() is not implemented yet.");
        return;
        return `hub ${this.identifier}
shortLabel ${this.shortLabel}
longLabel ${this.longLabel}
genomesFile genomes.txt
email ${this.email}`;
    }
}

// Track class to represent a UCSC Track
export class Track {
    /**
     * @param {string} identifier - Unique identifier for the track.
     * @param {string} shortLabel - Short label for the track.
     * @param {string} longLabel - Optional. Long label for the track.
     * @param {string} tracktype - Type of the track (e.g., bigWig, bigBed).
     * @param {string} url - URL to the track data file.
     * @param {string} visibility - Optional. Visibility setting for the track (e.g., "dense", "full").
     * @param {string} color - Optional. RGB color for the track (e.g., "255,0,0").
     * @param {string|null} parent - Identifier of the parent track (if applicable).
     */
    constructor(identifier, shortLabel, longLabel, tracktype, url, visibility, color, parent = null) {
        this.identifier = identifier;
        this.shortLabel = shortLabel;
        this.longLabel = longLabel || "";
        this.tracktype = tracktype;
        this.url = url;
        this.visibility = visibility || "dense";
        this.color = color || "0,0,0";
        this.parent = parent;
    }

    /**
     * Generates the trackDb.txt entry for this track.
     * @returns {string} The trackDb.txt entry.
     */
    generateTrackDbEntry() {
        new Error("generateTrackDbEntry() is not implemented yet.");
        return;
        let entry = `track ${this.identifier}
shortLabel ${this.shortLabel}
longLabel ${this.longLabel}
type ${this.tracktype}
bigDataUrl ${this.url}
visibility ${this.visibility}
color ${this.color}`;

        if (this.parent) {
            entry += `\nparent ${this.parent}`;
        }

        return entry;
    }
}