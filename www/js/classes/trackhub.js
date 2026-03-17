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
        this.trackDbUrl = null; // URL to the trackDb.txt file for this hub (to be populated after parsing the hub URL)
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

export class HubContainer {

    constructor() {
        this.hub = new Hub();
        this.addHubContainerEvents()
    }

    addHubContainerEvents() {
        document.getElementById("hub-identifier").addEventListener('input', (e) => {
            this.hub.identifier = e.target.value.trim();
        });

        document.getElementById("hub-short-label").addEventListener('input', (e) => {
            this.hub.shortLabel = e.target.value.trim();
        });

        document.getElementById("hub-long-label").addEventListener('input', (e) => {
            this.hub.longLabel = e.target.value.trim();
        });

        document.getElementById("hub-email").addEventListener('input', (e) => {
            this.hub.email = e.target.value.trim();
        });

        document.getElementById("hub-genome-select").addEventListener('change', (e) => {
            this.hub.genome = e.target.value.trim();
        });
    }

    /**
     * Parses the provided UCSC Track Hub URL and assembly to extract hub metadata and trackDb information.
     * Validates the hub URL, fetches its content, and constructs the trackDb URL based on the provided assembly.
     * Populates the hub data and updates the associated form fields.
     *
     * @async
     * @function parseHubUrl
     * @param {string} hubUrl - The URL of the UCSC Track Hub (e.g., "http://example.com/hub.txt").
     * @param {string} assembly - The genome assembly to use (e.g., "hg38").
     * @throws {Error} If the hub URL or assembly is missing, or if the hub URL is unreachable.
     * @returns {Promise<void>} Resolves when the hub data is successfully parsed and populated.
     */
    async parseHubUrl(hubUrl, assembly) {
        if (!hubUrl) {
            throw new Error("Hub URL is required to parse trackDb information.");
        }
        if (!assembly) {
            throw new Error("Assembly is required to parse trackDb information.");
        }

        // Test if URL is reachable
        let hubContent
        try {
            const hubResp = await fetch(hubUrl);
            if (!hubResp.ok) {
                throw new Error(`Hub URL returned status ${hubResp.status}`);
            }
            hubContent = await hubResp.text();
        } catch (error) {
            throw new Error("Provided hub URL is not reachable.");
        }

        const hubJson = {
            hub: null,
            shortLabel: null,
            longLabel: null,
            email: null,
            genome: assembly,
            trackDbUrl: null
        }

        for (const line of hubContent.split('\n')) {
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
                case 'email':
                    hubJson.email = value;
                    break;
                default:
                    // Ignore other keys for now
                    break;
            }
        }

        const url = new URL(hubUrl);
        const hubUrlPathname = url.pathname.replace(/\/?$/, '/'); // Ensure pathname ends with a slash
        const hubUrlParent = hubUrlPathname.split('/').slice(0, -2).join('/'); // Get parent directory of the hub URL
        const trackDbUrl = `${url.origin}${hubUrlParent}/${assembly}/trackDb.txt`;
        try {
            const trackDbResp = await fetch(trackDbUrl);
            if (!trackDbResp.ok) {
                throw new Error(`TrackDb URL returned status ${trackDbResp.status}`);
            }
            hubJson.trackDbUrl = trackDbUrl;
        } catch (error) {
            throw new Error("Constructed trackDb URL is not reachable. Please check the hub URL and assembly.");
        }


        this.populateHubData(hubJson);
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
        this.hub.genome = hubJson.genome || "";
        this.hub.trackDbUrl = hubJson.trackDbUrl || null;

        // Update form fields with hub data
        document.getElementById("hub-identifier").value = this.hub.identifier;
        document.getElementById("hub-short-label").value = this.hub.shortLabel;
        document.getElementById("hub-long-label").value = this.hub.longLabel;
        document.getElementById("hub-email").value = this.hub.email;
        document.getElementById("hub-genome-select").value = this.hub.genome;
    }
}

// Track class to represent a UCSC Track
export class Track {

    static trackCount = 1; // Global track count for unique identifiers

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

    // Convert passed in color value to UCSC-style RGB requirement ("")
    static convertColorToRGB(color) {
        if (color.startsWith('#')) {
            // Convert hex color to RGB
            const r = parseInt(color.slice(1, 3), 16);
            const g = parseInt(color.slice(3, 5), 16);
            const b = parseInt(color.slice(5, 7), 16);
            return `${r},${g},${b}`;
        } else {
            throw new Error("Invalid color format. Use hex (e.g., #ff0000) or RGB (e.g., 255,0,0).");
        }
    }

}

export class TrackContainer {

    constructor(containerId, addBtnId) {
        this.containerId = containerId || 'tracks-container';
        this.addBtnId = addBtnId || 'add-track-btn';

        this.tracksContainer = document.getElementById(this.containerId);
        this.addTrackBtn = document.getElementById(this.addBtnId);
        this.tracks = {}; // Store track data with track ID as key

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
        const trackItem = document.createElement('div');
        trackItem.classList.add('track-item');
        trackItem.setAttribute('draggable', 'true');
        const trackId = String(Track.trackCount);
        this.tracks[trackId] = new Track();
        Track.trackCount++;

        // Collapsible header
        const collapsibleHeader = document.createElement('div');
        collapsibleHeader.classList.add('collapsible-header', 'my-3');
        collapsibleHeader.innerHTML = `
            <span data-orig-name="Track ${trackId}" class="track-title">Track ${trackId}</span>
            <span class="icon">
                <i class="mdi mdi-chevron-down"></i>
            </span>
        `;

        // Collapsible content
        const collapsibleContent = document.createElement('div');
        collapsibleContent.classList.add('collapsible-content');

        // Clone the track template and append to collapsible content
        const template = document.getElementById('tmpl-trackhub-track');
        const templateHTML = template.content.cloneNode(true);
        const templateElement = templateHTML.querySelector('.js-track-item');
        templateElement.id = `track-${trackId}`;

        collapsibleContent.appendChild(templateHTML);

        // Add event listener for removing the track
        collapsibleContent.querySelector('.remove-track-btn').addEventListener('click', () => {
            trackItem.remove();
            // remove from tracks object
            delete this.tracks[trackId];
        });

        // Add toggle functionality for collapsible content
        collapsibleHeader.addEventListener('click', () => {
            const isOpen = collapsibleContent.style.display === 'block';
            collapsibleContent.style.display = isOpen ? 'none' : 'block';
            collapsibleHeader.querySelector('.mdi').classList.toggle('mdi-chevron-down', isOpen);
            collapsibleHeader.querySelector('.mdi').classList.toggle('mdi-chevron-up', !isOpen);
        });

        // Update the header title when the identifier changes
        collapsibleContent.querySelector('.js-track-identifier').addEventListener('input', (e) => {
            const value = e.target.value.trim();
            const labelElt = collapsibleHeader.querySelector('.track-title');
            labelElt.textContent = value || labelElt.dataset.origName;

            // update the track object
            this.tracks[trackId].identifier = value;
        });

        collapsibleContent.querySelector('.js-track-short-label').addEventListener('input', (e) => {
            this.tracks[trackId].shortLabel = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-long-label').addEventListener('input', (e) => {
            this.tracks[trackId].longLabel = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-type').addEventListener('input', (e) => {
            this.tracks[trackId].tracktype = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-url').addEventListener('input', (e) => {
            this.tracks[trackId].url = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-visibility').addEventListener('input', (e) => {
            this.tracks[trackId].visibility = e.target.value.trim();
        });

        collapsibleContent.querySelector('.js-track-color').addEventListener('input', (e) => {
            try {
                this.tracks[trackId].color = this.tracks[trackId].convertColorToRGB(e.target.value.trim());
            } catch (error) {
                console.error("Error converting color... will ignore:", error);
            }
        });

        trackItem.appendChild(collapsibleHeader);
        trackItem.appendChild(collapsibleContent);
        this.tracksContainer.appendChild(trackItem);

        // Start with the content expanded for new tracks
        collapsibleContent.style.display = 'block';
        collapsibleHeader.querySelector('.mdi').classList.remove('mdi-chevron-down');
        collapsibleHeader.querySelector('.mdi').classList.add('mdi-chevron-up');
    };

    createNewTracksFromData(trackData) {
        // Create new track objects and form fields based on parsed trackDb.txt data
        const trackEntries = trackData.split('\ntrack ').slice(1); // Split into individual track entries
        trackEntries.forEach((entry, index) => {
            const trackId = String(Track.trackCount + index);
            this.createTrackItem(); // This will create a new track item and increment the track count
            this.populateTrackData(trackId, entry); // Populate the new track item with data
        });
    }

    parseTrackDbUrl(trackDbUrl) {
        if (!trackDbUrl) {
            throw new Error("TrackDb URL is required to parse track data.");
        }
    }

    populateTrackData(trackId, trackEntry) {
        // Parse trackDb.txt entry text data and populate the track object and form fields
        const lines = trackEntry.split('\n');
        const track = this.tracks[trackId];
        lines.forEach(line => {
            const [key, ...rest] = line.split(' ');
            const value = rest.join(' ').trim();
            switch (key) {
                case 'track':
                    track.identifier = value;
                    document.querySelector(`#track-${trackId} .js-track-identifier`).value = value;
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
                    track.url = value;
                    document.querySelector(`#track-${trackId} .js-track-url`).value = value;
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
                    console.warn(`Unknown key in trackDb.txt: ${key}`);
            }
        });
    }
}