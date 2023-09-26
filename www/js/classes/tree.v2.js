"use strict";

/*
Component to manage tree views for various data
*/

class Tree {
        /**
     * Initialize tree.
     * @constructor
     * @param {Object} - Tree data.
     */
    constructor({
        element // element node, not ID
        , searchElement // same as above
        , selectCallback
    } = {}) {
        this.element = element; // Element to generate the tree structure on
        this.searchElement = searchElement  // Element that is used to query and filter datasets (i.e. search text box)
        this.selectCallback = selectCallback || ((e) => {})   // Callback to execute if node is "activated"
    }

    // Generate the tree structure for the DOM and return the JSTree object.
    generateTree () {
        this.generateTreeData();

        this.tree = new mar10.Wunderbaum({
            // https://mar10.github.io/wunderbaum/api/interfaces/wb_options.WunderbaumOptions.html
            element: this.element,
            source: this.treeData,
            filter: {
                // https://github.com/mar10/wunderbaum/blob/c393a7a4e01a8c6bb9084002d3a70b9e18bcc974/src/wb_ext_filter.ts
                connectInput: this.searchElement,
                autoExpand: true,
                mode: "hide",
            },
            types: {
                // See https://stackoverflow.com/a/11508530 for how to use variables as object keys
                [this.nodeType]: {icon: this.leafIcon}
            },
            autoCollapse: true,
            debugLevel: 2,  // set to "warn" level
            //event handlers
            init: (e) => {
                // After tree is created and data is loaded
                e.tree.setFocus();
            },
            render: (e) => {
                // e.node was rendered. We may now modify the markup...
            },
            click: (e) => {
                /* Single-click folders to expand them  (previous default dbl-click) */
                // NOTE: Clicking fast will sometimes not cause the expand/collapse
               //if (e.node.type === "folder") e.node.setExpanded(!(e.node.isExpanded())); // toggle open/close
                // BUG: Enabling the above code breaks clickable chevrons so we need to fix those
            },
            activate: (e) => {
                this.selectCallback(e)
            }
        });
    }

    // Create a nested folder node.
    addFolder(treeData, folder, kwargs) {
        // We skip folder ID 0, since it handled as id:domain_node already
        if (folder.id == 0) {
            return false;
        }

        // folder.parent_id == null is acceptable and is the top-level.

        this.addNode(treeData, {}, folder.id, folder.parent_id, folder.label, "folder", {...kwargs})
    }

    // Add a node to the tree. Edits "treeData" inplace.
    addNode(treeData, usedIDs, id, parentID, text, defaultNodeType, kwargs={}) {

        if (defaultNodeType !== "folder") {
            kwargs["orig_id"] = id; // Keep original ID in case value needs to be passed downstream
            id = `${defaultNodeType}__${id}`;

            // Modifies the ID to prevent duplicates.
            while (usedIDs.hasOwnProperty(id)) {
                // Add a dash to the exiting ID
                id += '-';
            }
        }

        // parentID == null is the top-level.

        // [PARENT_ID, [POSITIONAL_ARGS], {KEY_VALUE_ARGS}]
        const nodeData = [parentID, [], {key: id, title: text, type: defaultNodeType, ...kwargs}];

        usedIDs[id] = true; // Used this property
        treeData.push(nodeData)

    }

    // Find first node that matches condition. Implements Wunderbaum's findFirst() method.
    findFirst(match) {
        return this.tree.findFirst(match);
    }

}

/**
 * Class representing a Projection Source selection tree
 * @extends Tree
 */
class ProjectionSourceTree extends Tree {
    /**
     * Initialize ProjectionSourceTree
     * @constructor
     * @param {Object} Data - Tree data
     */
    constructor({
        ...args
    }={}, weightedDomainGeneCarts, weightedGroupGeneCarts, weightedUserGeneCarts, weightedSharedGeneCarts, weightedPublicGeneCarts,
        unweightedDomainGeneCarts, unweightedGroupGeneCarts, unweightedUserGeneCarts, unweightedSharedGeneCarts, unweightedPublicGeneCarts) {
        super(args);

        this["weighted-list"] = {
            domainGeneCarts: weightedDomainGeneCarts || []
            , groupGeneCarts: weightedGroupGeneCarts || []
            , userGeneCarts: weightedUserGeneCarts || []
            , sharedGeneCarts: weightedSharedGeneCarts || []
            , publicGeneCarts: weightedPublicGeneCarts || []
        };
        this["unweighted-list"] = {
            domainGeneCarts: unweightedDomainGeneCarts || []
            , groupGeneCarts: unweightedGroupGeneCarts || []
            , userGeneCarts: unweightedUserGeneCarts || []
            , sharedGeneCarts: unweightedSharedGeneCarts|| []
            , publicGeneCarts: unweightedPublicGeneCarts || []
        };
        // This is needed so we can add folders with labels to the tree
        this["weighted-list"].folders = [];
        this["unweighted-list"].folders = [];

    }

    nodeType = "genecart";
    leafIcon = 'mdi mdi-cart-outline';
    usedIDs = {};


    getTotalWeightedCarts() {
        // get the total number of weighted gene carts
        return Object.keys(this["weighted-list"]).reduce((acc, curr) => acc + this["weighted-list"][curr].length, 0);
    }

    getTotalUnweightedCarts() {
        // get the total number of unweighted gene carts
        return Object.keys(this["unweighted-list"]).reduce((acc, curr) => acc + this["unweighted-list"][curr].length, 0);
    }

    generateTreeData() {
        // Due to loading custom folders from db, we cannot guarantee order of loading.
        // So we need to load a flat-file with parents instead of nested children
        // https://mar10.github.io/wunderbaum/#/tutorial/tutorial_initialize?id=flat-parent-referencing-list
        // [PARENT_ID, [POSITIONAL_ARGS], {KEY_VALUE_ARGS}]
        const treeData = [
            [null, [], {key: "unweighted_genes_node", title: `Unweighted Genes (${this.getTotalUnweightedCarts()})`, type: "folder"}],
            ["unweighted_genes_node", [], {key: "uw_domain_node", title: `Highlighted gene carts (${this["unweighted-list"].domainGeneCarts.length})`, type: "folder"}],
            ["unweighted_genes_node", [], {key: "uw_user_node", title: `Your gene carts (${this["unweighted-list"].userGeneCarts.length})`, type: "folder"}],
            ["unweighted_genes_node", [], {key: "uw_group_node", title: `Group gene carts (${this["unweighted-list"].groupGeneCarts.length})`, type: "folder"}],
            ["unweighted_genes_node", [], {key: "uw_shared_node", title: `Gene carts shared with you (${this["unweighted-list"].sharedGeneCarts.length})`, type: "folder"}],
            ["unweighted_genes_node", [], {key: "uw_public_node", title: `Public carts from other users (${this["unweighted-list"].publicGeneCarts.length})`, type: "folder"}],
            [null, [], {key: "weighted_genes_node", title: `Weighted Genes (${this.getTotalWeightedCarts()})`, type: "folder"}],
            ['weighted_genes_node', [], {key: 'w_domain_node', title: `Highlighted gene carts (${this["weighted-list"].domainGeneCarts.length})`, type: "folder"}],
            ['weighted_genes_node', [], {key: "w_user_node", title: `Your gene carts (${this["weighted-list"].userGeneCarts.length})`, type: "folder"}],
            ['weighted_genes_node', [], {key: "w_group_node", title: `Group gene carts (${this["weighted-list"].groupGeneCarts.length})`, type: "folder"}],
            ['weighted_genes_node', [], {key: "w_shared_node", title: `Gene carts shared with you (${this["weighted-list"].sharedGeneCarts.length})`, type: "folder"}],
            ['weighted_genes_node', [], {key: "w_public_node", title: `Public carts from other users (${this["weighted-list"].publicGeneCarts.length})`, type: "folder"}],
        ];

        // ? Currently item.value is the gene cart share ID,
        // ? but for unweighted gene carts it is the db ID.
        // ? Is it worth refactoring so that the item.value is the db ID for both types?
        // ? If so, then item.value will need to ensure that no duplicates occur.
        // ? Could make use of "data" attribute in the tree node to store the original value and cart share ID.

        // Add all the folders first
        /*this["weighted-list"].folders.forEach( (item) => {
            this.addFolder(treeData, item);
        });
        /*this["unweighted-list"].folders.forEach( (item) => {
            this.addFolder(treeData, item);
        });*/

        // Sort the cart contents alphabetically
        ["domainGeneCarts", "userGeneCarts", "groupGeneCarts", "sharedGeneCarts", "publicGeneCarts"].forEach(e => {
            this["weighted-list"][e].sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
            this["unweighted-list"][e].sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        })

        this["weighted-list"].domainGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'w_domain_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            });
        });

        this["weighted-list"].userGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'w_user_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["weighted-list"].groupGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'w_group_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["weighted-list"].sharedGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'w_shared_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["weighted-list"].publicGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'w_public_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["unweighted-list"].domainGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'uw_domain_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            });
        });

        this["unweighted-list"].userGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'uw_user_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["unweighted-list"].groupGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'uw_group_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["unweighted-list"].sharedGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'uw_shared_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        this["unweighted-list"].publicGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'uw_public_node', item.text, this.nodeType, {
                'gctype': item.gctype,
            })
        });

        // Create JSON tree structure for the data
        const fullTreeData = {
            "_format": "flat",  // Specify that data is a flat, parent-referencing list
            "_positional": [],  // Fails if this is not provided, though we are processing everything as kwargs
            "children": treeData
        };

        this.treeData = fullTreeData;
        return this.treeData;

    }

    generateTree() { super.generateTree(); }

    findFirst(match) { return super.findFirst(match); }

}

/**
 * Class representing a gene cart selection tree
 * @extends Tree
 */
class GeneCartTree extends Tree {
    /**
     * Initialize GeneCartTree
     * @constructor
     * @param {Object} Data - Tree data
     */

    constructor({
        ...args
    }={}, domainGeneCarts, groupGeneCarts, userGeneCarts, sharedGeneCarts, publicGeneCarts) {
        super(args);
        this.domainGeneCarts = domainGeneCarts || [];
        this.userGeneCarts = userGeneCarts || [];
        this.groupGeneCarts = groupGeneCarts || [];
        this.sharedGeneCarts = sharedGeneCarts || [] ;
        this.publicGeneCarts = publicGeneCarts || [] ;
        // This is needed so we can add folders with labels to the tree
        this.folders = [];
    }

    nodeType = "genecart";
    leafIcon = 'mdi mdi-cart-outline';
    usedIDs = {};

    generateTreeData() {
        // Due to loading custom folders from db, we cannot guarantee order of loading.
        // So we need to load a flat-file with parents instead of nested children
        // https://mar10.github.io/wunderbaum/#/tutorial/tutorial_initialize?id=flat-parent-referencing-list
        // [PARENT_ID, [POSITIONAL_ARGS], {KEY_VALUE_ARGS}]

        const treeData = [
            [null, [], {key: "domain_node", title: `Highlighted gene carts (${this.domainGeneCarts.length})`, type: "folder"}],
            [null, [], {key: "user_node", title: `Your gene carts (${this.userGeneCarts.length})`, type: "folder"}],
            [null, [], {key: "group_node", title: `Group gene carts (${this.groupGeneCarts.length})`, type: "folder"}],
            [null, [], {key: "shared_node", title: `Gene carts shared with you (${this.sharedGeneCarts.length})`, type: "folder"}],
            [null, [], {key: "public_node", title: `Public carts from other users (${this.publicGeneCarts.length})`, type: "folder"}],
        ];

        // Add all the folders first
        /*this.folders.forEach( (item) => {
            this.addFolder(treeData, item);
        });*/

        // Sort the cart contents alphabetically
        ["domainGeneCarts", "userGeneCarts", "groupGeneCarts", "sharedGeneCarts", "publicGeneCarts"].forEach(e => {
            this[e].sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        });

        this.domainGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'domain_node', item.text, this.nodeType)
        });

        this.userGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'user_node', item.text, this.nodeType)
        });

        this.groupGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'group_node', item.text, this.nodeType)
        });

        this.sharedGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'shared_node', item.text, this.nodeType)
        });

        this.publicGeneCarts.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'public_node', item.text, this.nodeType)
        });

        // Create JSON tree structure for the data
        const fullTreeData = {
            "_format": "flat",  // Specify that data is a flat, parent-referencing list
            "_positional": [],  // Fails if this is not provided, though we are processing everything as kwargs
            "children": treeData
        };

        this.treeData = fullTreeData;
        return this.treeData;

    }

    generateTree() { super.generateTree(); }

    findFirst(match) { return super.findFirst(match); }

}

/**
 * Class representing a profile selection tree
 * @extends Tree
 */
class ProfileTree extends Tree {
    /**
     * Initialize ProfileTree
     * @constructor
     * @param {Object} Data - Tree data
     */
    constructor({
        ...args
    }={}, domainProfiles, userProfiles, groupProfiles, sharedProfiles, publicProfiles) {
        super(args);
        this.domainProfiles = domainProfiles || [];
        this.userProfiles = userProfiles || [];
        this.groupProfiles = groupProfiles || [];
        this.sharedProfiles = sharedProfiles || [];
        this.publicProfiles = publicProfiles || [];
        // This is needed so we can add folders with labels to the tree
        this.folders = [];
    }

    nodeType = 'profile';
    leafIcon = "mdi mdi-view-grid-outline";
    usedIDs = {};

    generateTreeData() {
        // Create JSON tree structure for the data
        const treeData = [];

        // Add all the folders first
        $.each(this.folders, (_i, item) => {
            this.addFolder(treeData, item);
        });

        // Sort the cart contents alphabetically
        ["domainProfiles", "userProfiles", "groupProfiles", "sharedProfiles", "publicProfiles"].forEach(e => {
            this[e].sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        });

        // Load profiles into the tree data property
        this.domainProfiles.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 101, item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        this.userProfiles.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 102, item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        this.groupProfiles.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 103, item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            });
        });

        this.sharedProfiles.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 104, item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        this.publicProfiles.forEach( (item) => {
            this.addNode(treeData, this.usedIDs, item.value, 105, item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        // Create JSON tree structure for the data
        const fullTreeData = {
            "_format": "flat",  // Specify that data is a flat, parent-referencing list
            "_positional": [],  // Fails if this is not provided, though we are processing everything as kwargs
            "children": treeData
        };

        this.treeData = fullTreeData;
        return this.treeData;
    }

    generateTree() { super.generateTree(); }

    findFirst(match) { return super.findFirst(match); }

}

/**
 * Class representing a dataset selection tree
 * @extends Tree
 */
class DatasetTree extends Tree {
    /**
     * Initialize ProfileTree
     * @constructor
     * @param {Object} Data - Tree data
     */
    constructor({
        ...args
    }={}, domainDatasets, sharedDatasets, userDatasets) {
        super(args);
        this.domainDatasets = domainDatasets || [];
        this.sharedDatasets = sharedDatasets || [];
        this.userDatasets = userDatasets || [];
        // This is needed so we can add folders with labels to the tree
        this.folders = [];
    }

    nodeType = 'dataset';
    treeKeys = {"domain_node": true, "shared_node": true, "user_node": true};
    leafIcon = "mdi mdi-card-account-details-outline";
    usedIDs = {};

    generateTreeData() {
        // Due to loading custom folders from db, we cannot guarantee order of loading.
        // So we need to load a flat-file with parents instead of nested children
        // https://mar10.github.io/wunderbaum/#/tutorial/tutorial_initialize?id=flat-parent-referencing-list
        // [PARENT_ID, [POSITIONAL_ARGS], {KEY_VALUE_ARGS}]

        const treeData = [
            [null, [], {key: "domain_node", title: `Public datasets (${this.domainDatasets.length})`, type: "folder"}],
            [null, [], {key: "shared_node", title: `Shared datasets (${this.sharedDatasets.length})`, type: "folder"}],
            [null, [],{key: "user_node", title: `User datasets (${this.userDatasets.length})`, type: "folder"}]
        ];

        // Load datasets into the tree data property
        // NOTE - Datasets can appear in multiple lists, so dataset IDs cannot be used as the node ID

        // Add all the folders first
        /*this.folders.forEach( (item) => {
            this.addFolder(treeData, item);
        });*/

        // Sort the cart contents alphabetically
        ["domainDatasets", "userDatasets", "sharedDatasets"].forEach(e => {
            this[e].sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        });

        this.domainDatasets.forEach((item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'domain_node', item.text, this.nodeType, {
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });

        this.sharedDatasets.forEach((item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'shared_node', item.text, this.nodeType, {
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });

        this.userDatasets.forEach((item) => {
            this.addNode(treeData, this.usedIDs, item.value, 'user_node', item.text, this.nodeType, {
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });

        // Create JSON tree structure for the data
        const fullTreeData = {
            "_format": "flat",  // Specify that data is a flat, parent-referencing list
            "_positional": [],  // Fails if this is not provided, though we are processing everything as kwargs
            "children": treeData
        };

        this.treeData = fullTreeData;
        return this.treeData;
    }

    generateTree() { super.generateTree(); }

    findFirst(match) { return super.findFirst(match); }
}