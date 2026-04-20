import { createToast } from "../common.v2.js";

export const TRACK_STATUS_COLORS = {
    downloading: { color: 'is-info', label: 'Downloading' },
    converting: { color: 'is-info', label: 'Converting' },
    ingesting: { color: 'is-info', label: 'Ingesting' },
    downloaded: { color: 'is-info', label: 'Downloaded' },
    completed: { color: 'is-success', label: 'Completed' },
    failed: { color: 'is-danger', label: 'Failed' },
};

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
        this.genomesFile = null; // URL to the genomes.txt file for this hub
        this.genome = genome;
        this.trackDbUrl = null; // URL to the trackDb.txt file for this hub
        this.extraKeys = {}; // Store any extra key-value pairs from hub.txt that are not part of the standard uploader form. Can be passed to UCSC Genome Browser.
    }

    /**
     * Generates the hub.txt content as a JSON object based on the current properties of the Hub instance.
     * @returns {Object} The hub.txt content as a JSON object.
     */
    generateHubJson(useOneFile=false) {
        // Neither "genome" nor "genomesFile" are used in the API but it is a sanity check.
        let genomeBit;
        if (useOneFile) {
            genomeBit = {
                useOneFile: "on",
                genome: this.genome
            }
        } else {
            genomeBit = {
                useOneFile: "off",
                genomesFile: "genomes.txt"
            }
        }

        return {
            hub: this.identifier,
            shortLabel: this.shortLabel,
            longLabel: this.longLabel,
            email: this.email,
            extraKeys: this.extraKeys,
            ...genomeBit
        };
    }
}

export class HubContainer {

    constructor() {
        this.hub = new Hub();
        this.oneFileMode = false;
        this.addHubContainerEvents()
        this.trackContainerObj;
        this.hubContent = null;
    }

    addHubContainerEvents() {
        document.getElementById("hub-identifier").addEventListener('input', (e) => {
            document.getElementById("hub-identifier").classList.remove("is-danger");
            this.hub.identifier = e.target.value.trim();
        });

        document.getElementById("hub-short-label").addEventListener('input', (e) => {
            document.getElementById("hub-short-label").classList.remove("is-danger");
            this.hub.shortLabel = e.target.value.trim();
        });

        document.getElementById("hub-long-label").addEventListener('input', (e) => {
            document.getElementById("hub-long-label").classList.remove("is-danger");
            this.hub.longLabel = e.target.value.trim();
        });

        document.getElementById("hub-email").addEventListener('input', (e) => {
            document.getElementById("hub-email").classList.remove("is-danger");
            this.hub.email = e.target.value.trim();
        });

        document.getElementById("hub-genome-select").addEventListener('change', (e) => {
            document.getElementById("hub-genome-select").classList.remove("is-danger");
            this.hub.genome = e.target.value.trim();
            document.querySelector(".js-human-assembly-warning").style.display = "none";
            if (["hg19", "hg38"].includes(this.hub.genome)) {
                document.querySelector(".js-human-assembly-warning").style.display = "block";
            }
            // Restrict track types based on assembly selection
            if (this.trackContainerObj) {
                this.trackContainerObj.restrictPrivacyAwareTrackTypes(this.hub.genome);
            }
        });
    }

    /**
     * Parses raw hub.txt content and populates the Hub object with extracted metadata.
     *
     * Handles hub.txt key-value pairs, detects oneFileMode, validates genome entries, and updates form fields.
     * If oneFileMode is enabled and genome mismatches the selected assembly, switches to genome found in hub.txt.
     *
     * @param {string} assembly - Selected genome assembly (validated against hub.txt in oneFileMode).
     * @returns {void}
     * @throws {Error} If oneFileMode is enabled but no genome entry found in hub.txt.
     *
     * @description
     * Expected hub.txt format:
     * ```
     * hub myHub
     * shortLabel My Hub
     * longLabel My Hub Long Label
     * email user@example.com
     * genomesFile genomes.txt
     * ```
     *
     * @see {@link parseHubFile} for file-based parsing
     * @see {@link parseHubUrl} for URL-based parsing
     * @see {@link populateHubData} for updating form fields
     */
    parseHubContent(assembly) {
        const hubJson = {
            hub: null,
            shortLabel: null,
            longLabel: null,
            email: null,
            genomesFile: null,
            genome: assembly,
            trackDbUrl: null,
            extraKeys: {}
        }

        // TODO: Instead of passing "assembly" from previous page, read genomes.txt or the genome entry to determine what genomes are selectable.

        let genomeInHubTxt = "";

        lineLoop: for (const line of this.hubContent.split('\n')) {
            const [key, ...rest] = line.split(' ');
            const value = rest.join(' ').trim();
            switch (key) {
                case 'hub':
                    hubJson.hub = value;
                    break;
                case 'shortLabel':
                    hubJson.shortLabel = value;
                    break;
                case 'longLabel':
                    hubJson.longLabel = value;
                    break;
                case 'genomesFile':
                    hubJson.genomesFile = value;
                    break;
                case 'email':
                    hubJson.email = value;
                    break;
                case "useOneFile":
                    if (value === "on")
                        this.oneFileMode = true;
                    break;
                case "genome":
                    // break out of the loop, so that existing track properties do not overwrite the hub ones with the same name
                    genomeInHubTxt = value;
                    break lineLoop;
                default:
                    // Store other keys
                    if (key) {
                        hubJson.extraKeys[key] = value;
                    }
                    break;
            }
        }

        if (this.oneFileMode) {
            if (! genomeInHubTxt) {
                throw new Error("Hub is in one-file mode but no genome entry found in hub.txt.");
            }
            if (genomeInHubTxt !== assembly) {
                createToast(`Hub is in one-file mode but genome entry in hub.txt (${genomeInHubTxt}) does not match the specified assembly (${assembly}). Switching to genome found in hub.txt.`, 'is-warning');
                hubJson.genome = genomeInHubTxt;
            }
        }

        this.populateHubData(hubJson);
    }

    /**
     * Parses hub.txt file content and populates the Hub object.
     *
     * @async
     * @param {File} file - File object from HTML file input.
     * @param {string} assembly - Genome assembly (e.g., "hg38", "mm10").
     * @returns {Promise<void>}
     * @throws {Error} If file or assembly is missing.
     *
     * @example
     * await hubContainer.parseHubFile(fileInput.files[0], 'hg38');
     *
     * @see {@link parseHubContent} for parsing logic
     * @see {@link parseHubUrl} for URL-based parsing
     */
    async parseHubFile(file, assembly) {
        if (!file) {
            throw new Error("No file provided for parsing.");
        }
        if (!assembly) {
            throw new Error("Assembly is required to parse trackDb information.");
        }
        const fileContent = await file.text();
        this.hubContent = fileContent;
        this.parseHubContent(assembly);
    }

    /**
     * Fetches and parses hub.txt file from a remote URL.
     *
     * @async
     * @param {string} hubUrl - Complete URL to hub.txt (e.g., "https://example.com/hub.txt").
     * @param {string} assembly - Genome assembly (e.g., "hg38", "mm10").
     * @returns {Promise<void>}
     * @throws {Error} If hubUrl/assembly is missing or URL is unreachable.
     *
     * @example
     * await hubContainer.parseHubUrl('https://example.com/hub.txt', 'hg38');
     *
     * @see {@link parseHubContent} for parsing logic
     * @see {@link parseHubFile} for file upload parsing
     * @see {@link retrieveTrackDbPath} for fetching trackDb.txt URL from genomes.txt
     */
    async parseHubUrl(hubUrl, assembly) {
        if (!hubUrl) {
            throw new Error("Hub URL is required to parse trackDb information.");
        }
        if (!assembly) {
            throw new Error("Assembly is required to parse trackDb information.");
        }

        // Test if URL is reachable
        let hubContent;
        try {
            const hubResp = await fetch(hubUrl);
            if (!hubResp.ok) {
                throw new Error(`Hub URL returned status ${hubResp.status}`);
            }
            hubContent = await hubResp.text();
        } catch (error) {
            console.debug(error);
            throw new Error("Provided hub URL is not reachable.");
        }

        this.hubContent = hubContent;
        this.parseHubContent(assembly);
    }

    /**
     * Populates the Hub object and associated form fields with data from a JSON object.
     * This method parses the provided `hubJson` object and updates both the Hub instance
     * and the corresponding HTML form fields with the extracted values.
     *
     * @param {Object} hubJson - JSON object containing hub data.
     *
     * @returns {void}
     */
    populateHubData(hubJson) {
        // Parse hub.txt json data and populate the hub object and form fields
        if (!hubJson) {
            console.warn("No hub data found to populate.");
            return;
        }

        this.hub.identifier = hubJson.hub || "";
        this.hub.shortLabel = hubJson.shortLabel || "";
        this.hub.longLabel = hubJson.longLabel || "";
        this.hub.email = hubJson.email || "";
        this.hub.genomesFile = hubJson.genomesFile || null;
        this.hub.genome = hubJson.genome || "";
        this.hub.trackDbUrl = hubJson.trackDbUrl || null;
        this.hub.extraKeys = hubJson.extraKeys || {};

        // Update form fields with hub data
        document.getElementById("hub-identifier").value = this.hub.identifier;
        document.getElementById("hub-short-label").value = this.hub.shortLabel;
        document.getElementById("hub-long-label").value = this.hub.longLabel;
        document.getElementById("hub-email").value = this.hub.email;
        document.getElementById("hub-genome-select").value = this.hub.genome;
        // Make genome read-only, and check for human assembly.
        document.getElementById("hub-genome-select").setAttribute("disabled", "disabled");
        document.querySelector(".js-human-assembly-warning").style.display = "none";
        if (["hg19", "hg38"].includes(this.hub.genome)) {
            document.querySelector(".js-human-assembly-warning").style.display = "block";
            // Will restrict when creating tracks, since track types are in the track form, not the hub form.
        }
    }

    /**
     * Asynchronously retrieves the trackDb.txt file path from a UCSC Track Hub's genomes.txt file.
     *
     * Constructs the URL to the genomes.txt file, fetches and parses it to locate the entry
     * matching the specified genome assembly, then extracts the trackDb URL for that assembly.
     * The retrieved trackDb URL is stored in the hub's trackDbUrl property.
     *
     * @async
     * @param {string} hubUrl - The URL of the hub.txt file (e.g., "https://example.com/hubs/hub1/hub.txt").
     * @returns {Promise<void>} Resolves when the trackDb URL has been successfully retrieved and stored.
     * @throws {Error} If the genomes.txt file cannot be fetched or parsed.
     * @throws {Error} If the specified assembly is not found in genomes.txt or has no trackDb entry.
     *
     * @description
     * This method performs the following steps:
     * 1. Parses the hubUrl to extract the origin and parent directory path.
     * 2. Constructs the full genomes.txt URL using the origin, parent path, and genomesFile property.
     * 3. Fetches the genomes.txt file from the constructed URL.
     * 4. Parses the genomes.txt content line-by-line to find the entry matching this.hub.genome.
     * 5. Scans the lines following the matching genome entry for a "trackDb" key.
     * 6. Constructs the absolute trackDb URL and stores it in this.hub.trackDbUrl.
     * 7. Stops scanning at the next "genome" entry or end of file.
     *
     * Expected genomes.txt format:
     * ```
     * genome hg38
     * trackDb genome_hg38/trackDb.txt
     *
     * genome mm10
     * trackDb genome_mm10/trackDb.txt
     * ```
     *
     * @see {@link parseHubUrl} for populating the genomesFile property before calling this method
     * @see {@link Hub.trackDbUrl} for the property where the retrieved URL is stored
     */
    async retrieveTrackDbPath(hubUrl) {
        const url = new URL(hubUrl);
        const hubUrlPathname = url.pathname.replace(/\/?$/, '/'); // Ensure pathname ends with a slash
        const hubUrlParent = hubUrlPathname.split('/').slice(0, -2).join('/'); // Get parent directory of the hub URL
        const assembly = this.hub.genome

        // Fetch and parse the genomes.txt file to validate the assembly and get trackDb URL
        const genomesTxtUrl = `${url.origin}/${hubUrlParent}/${this.hub.genomesFile}`;
        try {
            const resp = await fetch(genomesTxtUrl);
            if (!resp.ok) {
                throw new Error(`genomes.txt URL returned status ${resp.status}`);
            }
            const genomesTxtContent = await resp.text();
            const lines = genomesTxtContent.split('\n');
            for (const line of lines) {
                const [key, value] = line.split(' ').map(s => s.trim());
                if (key === 'genome' && value === assembly) {
                    // look for a trackDb "key" line in the following lines until the next "genome" line or end of file
                    for (const subLine of lines.slice(lines.indexOf(line) + 1)) {
                        const [subKey, subValue] = subLine.split(' ').map(s => s.trim());
                        if (subKey === 'genome') {
                            break; // stop searching if we reach the next genome entry
                        }
                        if (subKey === 'trackDb') {
                            const trackDbUrl = `${url.origin}${hubUrlParent}/${subValue}`;
                            this.hub.trackDbUrl = trackDbUrl;
                            return;
                        }
                    }
                }
            }
            throw new Error(`Assembly ${assembly} not found in genomes.txt or trackDb entry missing for assembly.`);
        } catch (error) {
            console.debug(error);
            throw new Error("Could not fetch or parse genomes.txt. Please check the hub URL and assembly.");
        }
    }

    validateHub() {
        // Validate that required hub fields are present and correctly formatted.
        // Return a list of errors or an empty list if valid.
        // In addition, higlight the form fields with errors for the user to fix
        const errors = [];
        let hasLocalFileReferences = false;

        if (this.hub.identifier) {
            // Hub identifier must be a single word with no spaces
            if (/\s/.test(this.hub.identifier)) {
                document.getElementById("hub-identifier").classList.add("is-danger");
                errors.push("Hub identifier must be a single word with no spaces.");
            }
        } else {
            document.getElementById("hub-identifier").classList.add("is-danger");
            errors.push("Hub identifier is required.");
        }

        if (!this.hub.shortLabel) {
            document.getElementById("hub-short-label").classList.add("is-danger");
            errors.push("Hub short label is required.");
        }

        if (!this.hub.longLabel) {
            document.getElementById("hub-long-label").classList.add("is-danger");
            errors.push("Hub long label is required.");
        }

        if (!this.hub.email) {
            document.getElementById("hub-email").classList.add("is-danger");
            errors.push("Hub email is required.");
        }

        if (!this.hub.genome) {
            document.getElementById("hub-genome-select").classList.add("is-danger");
            errors.push("Hub genome assembly is required.");
        }

        // Check for local file references in track stanzas
        if (this.trackContainerObj) {
            for (const trackId in this.trackContainerObj.tracks) {
                const track = this.trackContainerObj.tracks[trackId];
                const bigDataUrl = track.url || "";

                // Check if URL is relative or local (not http/https and doesn't start with /)
                if (bigDataUrl && !bigDataUrl.startsWith("http://") && !bigDataUrl.startsWith("https://") && !bigDataUrl.startsWith("/")) {
                    hasLocalFileReferences = true;
                    break;
                }
            }
        }

        return {
            errors,
            hasLocalFileReferences
        };
    }

    generateHubJson() {
        return this.hub.generateHubJson(this.oneFileMode);
    }

    getAssembly() {
        return this.hub.genome;
    }

    getTrackDbUrl() {
        return this.hub.trackDbUrl;
    }

    setTrackContainer(trackContainer) {
        this.trackContainerObj = trackContainer;
    }

    /**
     * Static method to validate hub.txt file content for local file references.
     * Does not require a HubContainer instance.
     *
     * @static
     * @param {string} hubContent - The raw text content of the hub.txt file.
     * @returns {boolean} True if the hub contains local file references, false otherwise.
     */
    static hasLocalFileReferences(hubContent) {
        const lines = hubContent.split('\n');
        for (const line of lines) {
            const [key, ...rest] = line.split(' ');
            const value = rest.join(' ').trim();

            // If there is a reference to a genomesFile, then we can assume there are local file references,
            // since the genomes.txt file must be local to the hub.txt file and we do not support external URLs for genomes.txt.
            if (key === 'genomesFile') {
                return true
            }

            // The oneFileMode property is probably "on". Look for bigDataUrl entries (track property)
            if (key === 'bigDataUrl') {
                // Check if URL is relative or local (not http/https and doesn't start with /)
                if (value && !value.startsWith('http://') && !value.startsWith('https://') && !value.startsWith('/')) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Static method to check if hub.txt content indicates useOneFile mode.
     * Does not require a HubContainer instance.
     */
    static hasUseOneFileMode(hubContent) {
        const lines = hubContent.split('\n');
        for (const line of lines) {
            const [key, value] = line.split(' ');
            if (key === 'useOneFile' && value === 'on') {
                return true;
            }
        }
        return false;
    }

}

// Track class to represent a UCSC Track
export class Track {

    static trackCount = 0; // Global track count for unique identifiers. Immediately increments when new track added.

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
        this.extraKeys = {}; // Store any extra key-value pairs from trackDb.txt that are not part of the standard uploader form. Can be passed to UCSC Genome Browser.
    }

    /**
     * Generates the trackDb.txt entry for this track as a JSON object based on the current properties of the Track instance.
     * @returns {Object} The trackDb.txt entry as a JSON object.
     */
    generateTrackDbEntry() {
        return {
            track: this.identifier,
            shortLabel: this.shortLabel,
            longLabel: this.longLabel,
            type: this.tracktype,
            bigDataUrl: this.url,
            visibility: this.visibility,
            color: this.color,
            parent: this.parent,
            extraKeys: this.extraKeys
        };
    }

    // Convert passed in color value to UCSC-style RGB requirement ("")
    static convertColorToRGB(color) {
        if (color.startsWith('#')) {
            // Convert hex color to RGB
            const r = parseInt(color.slice(1, 3), 16);
            const g = parseInt(color.slice(3, 5), 16);
            const b = parseInt(color.slice(5, 7), 16);
            return `${r},${g},${b}`;
        }
        throw new Error("Invalid color format. Expected hex (e.g., #ff0000)");
    }

}

export class TrackContainer {

    constructor(containerId, addBtnId) {
        this.containerId = containerId || 'tracks-container';
        this.addBtnId = addBtnId || 'add-track-btn';

        this.tracksContainer = document.getElementById(this.containerId);
        this.addTrackBtn = document.getElementById(this.addBtnId);

        this.trackDbUrl = null;
        this.tracks = {}; // Store track data with track ID as key
        this.trackFilesMap = new Map()

        this.hubContainerObj;

        this.createTrackItem = this.createTrackItem.bind(this); // Bind the method
        this.addTrackContainerEvents()
    }

    addTrackContainerEvents() {
        // Add new track on button click
        this.addTrackBtn.addEventListener('click', this.createTrackItem);

        // Drag-and-drop reordering
        this.tracksContainer.addEventListener('dragstart', (e) => {
            if (e.target.classList.contains('track-item')) {
                e.target.classList.add('dragging');
            }
        });

        this.tracksContainer.addEventListener('dragend', (e) => {
            if (e.target.classList.contains('track-item')) {
                e.target.classList.remove('dragging');
            }
        });

        this.tracksContainer.addEventListener('dragover', (e) => {
            e.preventDefault();
            const draggingItem = document.querySelector('.dragging');
            const afterElement = this.getDragAfterElement(e.clientY);
            if (afterElement == null) {
                this.tracksContainer.appendChild(draggingItem);
            } else {
                this.tracksContainer.insertBefore(draggingItem, afterElement);
            }
        });
    }

    getDragAfterElement = (y) => {
        const draggableElements = [...this.tracksContainer.querySelectorAll('.track-item:not(.dragging)')];
        return draggableElements.reduce((closest, child) => {
            const box = child.getBoundingClientRect();
            const offset = y - box.top - box.height / 2;
            if (offset < 0 && offset > closest.offset) {
                return { offset, element: child };
            } else {
                return closest;
            }
        }, { offset: Number.NEGATIVE_INFINITY }).element;
    };

    // Function to create a new track item with event listeners
    createTrackItem() {
        Track.trackCount++;
        const trackItem = document.createElement('div');
        trackItem.classList.add('track-item');
        trackItem.setAttribute('draggable', 'true');
        const trackId = String(Track.trackCount);
        trackItem.id = `track-${trackId}`;

        this.tracks[trackId] = new Track();

        // Collapsible header
        const collapsibleHeader = document.createElement('div');
        collapsibleHeader.classList.add('collapsible-header', 'my-3');
        collapsibleHeader.innerHTML = `
            <span data-orig-name="Track ${trackId}" class="track-title">Track ${trackId}</span>
            <div class="track-status-container" style="display: flex; align-items: center; gap: 0.5rem;">
                <span class="track-status-danger has-text-danger-dark icon is-hidden">
                    <i class="mdi mdi-alert"></i>
                </span>
                <span class="track-status-arrow icon">
                    <i class="mdi mdi-chevron-down"></i>
                </span>
            </div>
        `;

        // Collapsible content
        const replaceTemplatePlaceholders = (documentFragment, trackId) => {
            const tempDiv = document.createElement('div');
            tempDiv.classList.add('collapsible-content');
            tempDiv.appendChild(documentFragment);
            const html = tempDiv.innerHTML.replace(/{trackId}/g, trackId);
            tempDiv.innerHTML = html;
            return tempDiv;
        };

        // Clone the track template and append to collapsible content
        const template = document.getElementById('tmpl-trackhub-track');
        const templateClone = template.content.cloneNode(true);

        const collapsibleContent = replaceTemplatePlaceholders(templateClone, trackId);

        // Add event listener for removing the track
        collapsibleContent.querySelector('.remove-track-btn').addEventListener('click', () => {
            trackItem.remove();
            // remove from tracks object
            delete this.tracks[trackId];
            this.trackFilesMap.delete(trackId);
        });

        // Add toggle functionality for collapsible content
        collapsibleHeader.addEventListener('click', (event) => {
            const isOpen = collapsibleContent.style.display === 'block';
            collapsibleContent.style.display = isOpen ? 'none' : 'block';
            collapsibleHeader.querySelector('.track-status-arrow .mdi').classList.toggle('mdi-chevron-down', isOpen);
            collapsibleHeader.querySelector('.track-status-arrow .mdi').classList.toggle('mdi-chevron-up', !isOpen);
        });

        // Update the header title when the identifier changes
        collapsibleContent.querySelector('.js-track-identifier').addEventListener('input', (e) => {
            const value = e.target.value.trim();
            const labelElt = collapsibleHeader.querySelector('.track-title');
            labelElt.textContent = value || labelElt.dataset.origName;

            document.querySelector(`#track-${trackId} .js-track-identifier`).classList.remove("is-danger");
            // update the track object
            this.tracks[trackId].identifier = value;
        });

        collapsibleContent.querySelector('.js-track-short-label').addEventListener('input', (e) => {
            document.querySelector(`#track-${trackId} .js-track-short-label`).classList.remove("is-danger");

            this.tracks[trackId].shortLabel = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-long-label').addEventListener('input', (e) => {
            document.querySelector(`#track-${trackId} .js-track-long-label`).classList.remove("is-danger");
            this.tracks[trackId].longLabel = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-type').addEventListener('input', (e) => {
            document.querySelector(`#track-${trackId} .js-track-type`).classList.remove("is-danger");
            this.tracks[trackId].tracktype = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-url').addEventListener('input', (e) => {
            document.querySelector(`#track-${trackId} .js-track-url`).classList.remove("is-danger");
            this.tracks[trackId].url = e.target.value.trim();
            if (this.tracks[trackId].url) {
                // If there is a URL, remove any file that was selected for this track
                this.trackFilesMap.delete(trackId);
                collapsibleContent.querySelector('.js-track-file').value = ""; // Clear the file input

                // remove warning from header
                const dangerElt = document.querySelector(`#track-${trackId} .track-status-danger`);
                dangerElt.classList.add("is-hidden");
            }
        });

        // Capture file selections without uploading
        collapsibleContent.querySelector('.js-track-file').addEventListener('change', (e) => {
            if (e.target.files.length > 0) {
                this.trackFilesMap.set(trackId, e.target.files[0]);
                document.querySelector(`#track-${trackId} .file-name`).textContent = e.target.files[0].name;

                // If a file is selected, clear the URL field for this track
                this.tracks[trackId].url = "";
                collapsibleContent.querySelector('.js-track-url').value = ""; // Clear the URL input

                // remove warning from header
                const dangerElt = document.querySelector(`#track-${trackId} .track-status-danger`);
                dangerElt.classList.add("is-hidden");

            } else {
                this.trackFilesMap.delete(trackId);
            }
        });

        collapsibleContent.querySelector('.js-track-visibility').addEventListener('input', (e) => {
            this.tracks[trackId].visibility = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-color').addEventListener('input', (e) => {
            try {
                this.tracks[trackId].color = Track.convertColorToRGB(e.currentTarget.value.trim());
                e.currentTarget.classList.remove("is-danger");
            } catch (error) {
                e.currentTarget.classList.add("is-danger");
                console.warn("Error converting color... will ignore:", error);
            }
        });

        trackItem.appendChild(collapsibleHeader);
        trackItem.appendChild(collapsibleContent);
        this.tracksContainer.appendChild(trackItem);

        // Restrict track types for the newly added track
        if (this.hubContainerObj) {
            this.restrictPrivacyAwareTrackTypes(this.hubContainerObj.getAssembly(), collapsibleContent);
        }

        // Start with the content expanded for new tracks
        collapsibleContent.style.display = 'block';
        collapsibleHeader.querySelector('.track-status-arrow .mdi').classList.remove('mdi-chevron-down');
        collapsibleHeader.querySelector('.track-status-arrow .mdi').classList.add('mdi-chevron-up');
    };

    createNewTracksFromData(trackData) {
        // Create new track objects and form fields based on parsed trackDb.txt data
        const trackEntries = trackData.split('track ')
        for (const trackEntry of trackEntries) {
            // skip empty entries
            if (!trackEntry.trim()) {
                continue;
            }

            // Extract track type from stanza to check privacy awareness
            const trackTypeMatch = trackEntry.match(/type\s+(\S+)/);
            const trackType = trackTypeMatch ? trackTypeMatch[1] : '';

            // if assembly is human, skip creation of privacy-aware tracks.
            if (this.isAssemblyPrivacyAware(this.hubContainerObj.getAssembly()) && this.checkPrivacyTrackType(trackType)) {
                continue;
            }

            // add "track " back into the stanza since we split it out (it's a field name)
            const stanza = `track ${trackEntry}`;
            this.createTrackItem(); // This will create a new track item and increment the track count
            const trackId = String(Track.trackCount); // Get the current track ID after incrementing
            this.populateTrackData(trackId, stanza); // Populate the new track item with data
        };
    }

    parseHubTracksFromContent(hubContent) {
        // In this situation, the tracks are in the current hub.txt file.
        // File is a concatenation of the following order: hub.txt, genomes.txt, trackDb.txt
        // Extract the trackDb.txt portion of the file by finding the "track " keyword that starts the track entries
        const trackDbIndex = hubContent.indexOf('track ');
        if (trackDbIndex === -1) {
            throw new Error("No track entries found in the hub.txt file.");
        }
        const trackDbContent = hubContent.slice(trackDbIndex);
        this.createNewTracksFromData(trackDbContent);
    }

    async parseHubTracks(hubUrl) {
        // In this situation, the tracks are in the current hub.txt file.
        // File is a concatenation of the following order: hub.txt, genomes.txt, trackDb.txt
        let hubContent;
        try {
            const hubResp = await fetch(hubUrl);
            if (!hubResp.ok) {
                throw new Error(`Hub URL returned status ${hubResp.status}`);
            }
            hubContent = await hubResp.text();
        } catch (error) {
            console.debug(error);
            throw new Error("Provided hub URL is not reachable.");
        }

        this.trackDbUrl = hubUrl; // Nneeded for filepath resolution of bigDataUrl
        this.parseHubTracksFromContent(hubContent);
    }

    async parseTrackDbUrl(trackDbUrl) {
        if (!trackDbUrl) {
            throw new Error("TrackDb URL is required to parse track data.");
        }

        this.trackDbUrl = trackDbUrl;

        let trackDbContent;
        try {
            const trackDbResp = await fetch(trackDbUrl);
            if (!trackDbResp.ok) {
                throw new Error(`TrackDb URL returned status ${trackDbResp.status}`);
            }
            trackDbContent = await trackDbResp.text();
        } catch (error) {
            console.debug(error);
            throw new Error("Constructed trackDb URL is not reachable. Please check the hub URL and assembly.");
        }

        this.createNewTracksFromData(trackDbContent);
    }

    populateTrackData(trackId, stanza) {
        // Parse trackDb.txt entry text data and populate the track object and form fields
        const lines = stanza.split('\n');
        const track = this.tracks[trackId];

        for (const line of lines) {
            // if only whitespace or a comment, skip
            if (!line.trim() || line.trim().startsWith('#')) {
                return;
            }

            // Split to "key value"
            const [key, ...rest] = line.split(' ');
            const value = rest.join(" ").trim();
            switch (key) {
                case 'track':
                    track.identifier = value;
                    document.querySelector(`#track-${trackId} .js-track-identifier`).value = value;
                    const labelElt = document.querySelector(`#track-${trackId} .track-title`);
                    labelElt.textContent = track.identifier;
                    break;
                case 'shortLabel':
                    track.shortLabel = value;
                    document.querySelector(`#track-${trackId} .js-track-short-label`).value = value;
                    break;
                case 'longLabel':
                    track.longLabel = value;
                    document.querySelector(`#track-${trackId} .js-track-long-label`).value = value;
                    break;
                case 'type':
                    track.tracktype = value;
                    document.querySelector(`#track-${trackId} .js-track-type`).value = value;
                    break;
                case 'bigDataUrl':
                    if (value.startsWith('http://') || value.startsWith('https://')) {
                        // If URL is absolute, use as is
                        track.url = value;
                    } else if (this.trackDbUrl) {
                        // If URL is relative, convert to absolute based on trackDb URL
                        const trackDbUrlObj = new URL(this.trackDbUrl);
                        const trackDbBasePath = trackDbUrlObj.pathname.split('/').slice(0, -1).join('/');
                        const absoluteUrl = `${trackDbUrlObj.origin}${trackDbBasePath}/${value}`;
                        track.url = absoluteUrl;
                    } else {
                        // If there is no trackDbUrl, this must be a relative URL from a hub.txt file with local file references.
                        // Add some warning indicator to the track header so that the user looks.
                        track.url = "";
                        const dangerElt = document.querySelector(`#track-${trackId} .track-status-danger`);
                        dangerElt.classList.remove("is-hidden");
                    }
                    document.querySelector(`#track-${trackId} .js-track-url`).value = track.url;
                    break;
                case 'visibility':
                    track.visibility = value;
                    document.querySelector(`#track-${trackId} .js-track-visibility`).value = value;
                    break;
                case 'color':
                    track.color = value;
                    // Convert RGB back to hex for the color input field
                    const [r, g, b] = value.split(',').map(Number);
                    const hexColor = `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1)}`;
                    document.querySelector(`#track-${trackId} .js-track-color`).value = hexColor;
                    break;
                //case 'parent':
                //    track.parent = value
                //    break;
                default:
                    // populate extra keys
                    if (key) {
                        track.extraKeys[key] = value;
                    }
                    break;
            }
        };

        // Collapse the track item after populating data for better UX, especially for large trackhubs
        const collapsibleHeader = document.querySelector(`#track-${trackId} .collapsible-header`);
        collapsibleHeader.click(); // Simulate a click to collapse the item
    }

    validateTracks() {
        // For each track item in the container, validate that required fields are present and correctly formatted.
        // Return a list of errors or an empty list if all tracks are valid.
        // In addition, highlight the form fields with errors for the user to fix
        const errors = [];
        const missingFileOrUrl = {}; // Map trackId → true if both URL and file are missing

        for (const trackId in this.tracks) {
            const track = this.tracks[trackId];
            if (!track.identifier) {
                document.querySelector(`#track-${trackId} .js-track-identifier`).classList.add("is-danger");
                errors.push(`Track ${trackId}: Identifier is required.`);
            }
            if (!track.shortLabel) {
                document.querySelector(`#track-${trackId} .js-track-short-label`).classList.add("is-danger");
                errors.push(`Track ${trackId}: Short label is required.`);
            }
            if (!track.longLabel) {
                document.querySelector(`#track-${trackId} .js-track-long-label`).classList.add("is-danger");
                errors.push(`Track ${trackId}: Long label is required.`);
            }
            if (!track.tracktype) {
                document.querySelector(`#track-${trackId} .js-track-type`).classList.add("is-danger");
                errors.push(`Track ${trackId}: Track type is required.`);
            }
            // Check for either URL or uploaded file
            const hasUrl = track.url?.trim();
            const hasFile = this.trackFilesMap.has(String(trackId));

            if (!hasUrl && !hasFile) {
                document.querySelector(`#track-${trackId} .js-track-url`).classList.add("is-danger");
                errors.push(`Track ${trackId}: Must provide either a URL or upload a file.`);
                missingFileOrUrl[trackId] = true;
            }
        }

        return {
            errors,
            missingFileOrUrl
        };
    }

    /**
     * Generates a list of trackDb.txt entries as JSON objects based on the current track data in the container.
     * Each entry corresponds to a track and includes properties such as identifier, labels, type, URL, visibility, and color.
     * @returns {Object[]} List of trackDb.txt entries as JSON objects.
     */
    generateTrackDbEntries() {
        const entries = [];
        for (const trackId in this.tracks) {
            const track = this.tracks[trackId];
            entries.push(track.generateTrackDbEntry());
        }
        return entries;
    }

    /**
     * Builds FormData containing track file uploads keyed by track ID.
     * @returns {FormData} FormData with entries like `tracks[trackId][file]`
     */
    buildTrackFilesFormData() {
        const formData = new FormData();

        for (const [trackId, file] of this.trackFilesMap.entries()) {
            formData.append(`tracks[${trackId}][file]`, file);
        }

        return formData;
    }

    isAssemblyPrivacyAware(assembly) {
        return ["hg19", "hg38"].includes(assembly);
    }

    /**
     * Determines if a given track type is considered privacy-aware and should be restricted for human assemblies.
     * @returns {boolean} True if the track type is privacy-aware, false otherwise.
     */
    checkPrivacyTrackType(trackType) {
        const privacyAwareTypes = ["vcfTabix", "hic"];
        return privacyAwareTypes.includes(trackType);
    }

    /**
     * Disables selection of privacy-aware track types (e.g., VCF, Hi-C) in the track type dropdowns if the assembly is human.
     * This is to comply with federal standards regarding personally-identifiable data.
     * @param {string} assembly - The genome assembly to check for privacy awareness (e.g., "hg38", "mm10").
     * @param {HTMLElement} parent - Optional. The parent element to search within for track type select elements. Defaults to the entire document.
     * @returns {void}
     */
    restrictPrivacyAwareTrackTypes(assembly, parent=document) {
        const trackTypeSelectElts = parent.getElementsByClassName('js-track-type');
        for (const trackTypeSelect of trackTypeSelectElts) {
            trackTypeSelect.querySelector('option[value="vcfTabix"]').disabled = false;
            trackTypeSelect.querySelector('option[value="hic"]').disabled = false;
            if (this.isAssemblyPrivacyAware(assembly)) {
                trackTypeSelect.querySelector('option[value="vcfTabix"]').disabled = true;
                trackTypeSelect.querySelector('option[value="hic"]').disabled = true;
            }
        }
    }

    setHubContainer(hubContainer) {
        this.hubContainerObj = hubContainer;
    }

    registerTrackDataFileListener(callback) {
        // Call immediately to capture initial state
        callback();

        // Listen for file input changes
        this.tracksContainer.addEventListener('change', (e) => {
            if (e.target.classList.contains('js-track-file')) {
                callback();
            }
        });

        // Listen for URL input changes
        this.tracksContainer.addEventListener('input', (e) => {
            if (e.target.classList.contains('js-track-url')) {
                callback();
            }
        });

        // If a track is removed, we should also check if there are still any tracks with file references
        this.tracksContainer.addEventListener('click', (e) => {
            if (e.target.classList.contains('remove-track-btn')) {
                callback();
            }
        });

    }

    hasAtLeastOneTrackWithDataFile() {
        for (const trackId in this.tracks) {
            const track = this.tracks[trackId];
            const hasUrl = track.url?.trim();
            const hasFile = this.trackFilesMap.has(String(trackId));
            if (hasUrl || hasFile) {
                return true;
            }
        }
        return false;
    }

    /**
     * Static method to validate trackDb.txt file content for local file references.
     * Does not require a TrackContainer instance.
     *
     * @static
     * @param {string} trackDbContent - The raw text content of the trackDb.txt file.
     * @returns {boolean} True if trackDb contains local file references, false otherwise.
     */
    static hasLocalFileReferences(trackDbContent) {
        const lines = trackDbContent.split('\n');
        for (const line of lines) {
            const [key, ...rest] = line.split(' ');
            const value = rest.join(' ').trim();

            // Look for bigDataUrl entries
            if (key === 'bigDataUrl') {
                // Check if URL is relative or local
                if (value && !value.startsWith('http://') && !value.startsWith('https://') && !value.startsWith('/')) {
                    return true;
                }
            }
        }
        return false;
    }
}