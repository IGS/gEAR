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
}