"use strict";

/*
Gene Cart and Profile tree stuff

Just about everything here uses the JSTree library - jstree.com
*/

/** Base class representing a tree structure */
class Tree {
    /**
     * Initialize tree.
     * @constructor
     * @param {Object} - Tree data.
     */
    constructor({
        treeDiv
    } = {}) {
        this.treeDiv = treeDiv;

        this.setTree();
    }

    isTree() {
        return Boolean(this.tree)
    }

    setTree() {
        this.tree = (this.treeDiv) ? $.jstree.reference(this.treeDiv) : undefined;
        // $(this.treeDiv).jstree(true) returns the same as $.jstree.reference(this.treeDiv)
        // so they may be interchanged in the codebase
    }

    register_search() {
    // Code from "search" section of https://www.jstree.com/plugins/
    // Sets text input to search as tree search box.
    const self = this;
    let to = false;
    // Requires searchbox to be named #{treeDiv}_q
    $(`${this.treeDiv}_q`).keyup(() => {
        if (to) { clearTimeout(to); }
        to = setTimeout(() => {
        const v = $(`${self.treeDiv}_q`).val();
        self.tree.search(v);
        }, 250);
    });
    }

    // Update the contents of a JSTree object with new data.
    updateTreeData(newData=null) {
        (this.tree).settings.core.data = newData ? newData : this.treeData;
        (this.tree).refresh();
    }

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
        this.domainGeneCarts = (domainGeneCarts) ? domainGeneCarts : [];
        this.userGeneCarts = (userGeneCarts) ? userGeneCarts : [] ;
        this.groupGeneCarts = (groupGeneCarts) ? groupGeneCarts : [];
        this.sharedGeneCarts = (sharedGeneCarts) ? sharedGeneCarts : [] ;
        this.publicGeneCarts = (publicGeneCarts) ? publicGeneCarts : [] ;

    }

    generateTreeData() {
        // Create JSON tree structure for the data
        const treeData = [
            {'id':'domain_node', 'parent':'#', 'text':`Highlighted gene carts (${this.domainGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'user_node', 'parent':'#', 'text':`Your gene carts (${this.userGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'group_node', 'parent':'#', 'text':`Group gene carts (${this.groupGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'shared_node', 'parent':'#', 'text':`Gene carts shared with you (${this.sharedGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'public_node', 'parent':'#', 'text':`Public carts from other users (${this.publicGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
        ];

        $.each(this.domainGeneCarts, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'domain_node',
                'text':item.text,
                'type': 'genecart',
                'a_attr': {
                    'class': "py-0",
                }
            })
        });

        $.each(this.userGeneCarts, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'user_node',
                'text':item.text,
                'type': 'genecart',
                'a_attr': {
                    'class': "py-0",
                }
            })
        });

        $.each(this.groupGeneCarts, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'group_node',
                'text':item.text,
                'type': 'genecart',
                'a_attr': {
                    'class': "py-0",
                }
            })
        });

        $.each(this.sharedGeneCarts, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'shared_node',
                'text':item.text,
                'type': 'genecart',
                'a_attr': {
                    'class': "py-0",
                }
            })
        });

        $.each(this.publicGeneCarts, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'public_node',
                'text':item.text,
                'type': 'genecart',
                'a_attr': {
                    'class': "py-0",
                }
            })
        });

        this.treeData = treeData;
        return this.treeData;
    }

    // Load all saved gene carts for the current user
    // TODO: Change based on gene cart manager page code
    generateTree () {
        this.generateTreeData();

        // Update existing tree or generate new tree if it doesn't exist
        if (this.tree) {
            this.updateTreeData()
        } else {
            $(this.treeDiv).jstree({
                'core':{
                    'data':this.treeData,
                },
                'plugins': ["search", "types", "wholerow"],
                /* Plugins
                    search - search for matching items and expand tree if found
                    types - Allows you to define node types with nesting rules and icons
                    wholerow - makes each node block-level for easier selection
                */
                'search': {
                    "show_only_matches": true
                },
                'types': {
                    'default': {
                        'icon': 'fa fa-folder-o'
                    },
                    'genecart': {
                        'icon': 'fa fa-shopping-cart',
                        'valid_children':[]
                    }
                }
            })
            this.setTree();
        }

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first()
        this.register_events();

    }

    register_events() {
        const self = this;
        this.register_search();

        // Get genes from the selected gene cart
        $(this.treeDiv).on('select_node.jstree', (_e, data) => {
            // Though you can select multiple nodes in the tree, let's only select the first
            const geneCartId = data.selected[0];  // Returns node 'id' property
            if (data.node.type === "default") {
                // Do not toggle if user is navigating a branch node
                // NOTE: If tree is inside a <form>, which cannot be nested inside another <form>, this could toggle closed anyways due to the conflict.
                return;
            }
            const selectedNode = data.instance.get_node(geneCartId);
            $(self.dropdownToggleElt).text(selectedNode.text);
            $(self.dropdownToggleElt).val(geneCartId);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(self.dropdownToggleElt).change(); // Force the change event to fire, triggering downstream things like getting cart members
        }).jstree(true);
    }

    loadFromDB() {
        //pass
    }

    saveToDB() {
        //pass
    }
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
    }={}, domainProfiles, userProfiles, groupProfiles, sharedProfiles) {
        super(args);
        this.domainProfiles = (domainProfiles) ? domainProfiles : [];
        this.userProfiles = (userProfiles) ? userProfiles : [];
        this.groupProfiles = (groupProfiles) ? groupProfiles : [];
        this.sharedProfiles = (sharedProfiles) ? sharedProfiles : [];
    }

    generateTreeData() {
        // Create JSON tree structure for the data
        const treeData = [
            {'id':'domain_node', 'parent':'#', 'text':`Highlighted profiles (${this.domainProfiles.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'user_node', 'parent':'#', 'text':`Your profiles (${this.userProfiles.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'group_node', 'parent':'#', 'text':`Group profiles (${this.groupProfiles.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'shared_node', 'parent':'#', 'text':`Profiles shared with you (${this.sharedProfiles.length})`, 'a_attr': {'class':'jstree-ocl'}},
        ];

        // user_profiles/domain_profiles properties - value, text, share_id

        // Load profiles into the tree data property
        $.each(this.domainProfiles, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'domain_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': "py-0",
                },
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        $.each(this.userProfiles, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'user_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': "py-0",
                },
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        $.each(this.groupProfiles, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'group_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': "py-0",
                },
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        $.each(this.sharedProfiles, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'shared_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': "py-0",
                },
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        this.treeData = treeData;
        return this.treeData;
    }

    generateTree() {

        this.generateTreeData();

        // Update existing tree or generate new tree if it doesn't exist
        if (this.tree) {
            this.updateTreeData()
        } else {
            // Instantiate the tree
            $(this.treeDiv).jstree({
                'core':{
                    'data':this.treeData,
                },
                'plugins': ["search", "types", "wholerow"],
                'search': {
                    "show_only_matches": true
                },
                'types': {
                    'default': {
                        'icon': 'fa fa-folder-o'
                    },
                    'profile': {
                        'icon': 'fa fa-address-card-o',
                        'valid_children':[]
                    }
                }
            })
            this.setTree();
        }

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first()

        this.register_events();
    }

    // Register various ProfileTree events as object properties are updated.
    register_events() {
        const self = this;
        this.register_search();

        // Get layout from the selected node and close dropdown
        $(this.treeDiv).on('select_node.jstree', (_e, data) => {

            // Though you can select multiple nodes in the tree, let's only select the first
            const layoutId = data.selected[0];  // Returns node 'id' property
            if (data.node.type === "default") {
                // Do not toggle if user is navigating a branch node
                // NOTE: If tree is inside a <form>, which cannot be nested inside another <form>, this could toggle closed anyways due to the conflict.
                return;
            }
            // The dropdown toggle text/val change already happens in DatasetCollectionPanel->set_layouts() for the index page,
            // but this should be set to assist with other pages.
            const selectedNode = data.instance.get_node(layoutId);
            $(self.dropdownToggleElt).text(selectedNode.text);
            $(self.dropdownToggleElt).val(layoutId);
            $(self.dropdownToggleElt).data("profile-id", selectedNode.original.profile_id);
            $(self.dropdownToggleElt).data("profile-label", selectedNode.original.profile_label);
            $(self.dropdownToggleElt).data("profile-share-id", selectedNode.original.profile_share_id);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(self.dropdownToggleElt).trigger('change');   // Force the change event to fire, triggering downstream things

        }).jstree(true);
    }

    loadFromDB() {
        //pass
    }

    saveToDB() {
        //pass
    }

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
        this.domainDatasets = (domainDatasets) ? domainDatasets : [];
        this.sharedDatasets = (sharedDatasets) ? sharedDatasets : [];
        this.userDatasets = (userDatasets) ? userDatasets : [];
    }

    generateTreeData() {
        // Create JSON tree structure for the data
        const treeData = [
            {'id':'domain_node', 'parent':'#', 'text':`Public Datasets (${this.domainDatasets.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'shared_node', 'parent':'#', 'text':`Shared Datasets (${this.sharedDatasets.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'user_node', 'parent':'#', 'text':`Your Datasets (${this.userDatasets.length})`, 'a_attr': {'class':'jstree-ocl'}},
        ];

        // Load datasets into the tree data property
        // NOTE - Datasets can appear in multiple lists, so dataset IDs cannot be used as the node ID
        // otherwise node leaves can turn into "default" type instead of "dataset" type

        $.each(this.domainDatasets, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'domain_node',
                'text': item.text,
                'type': 'dataset',
                'a_attr': {
                    'class': "py-0",
                },
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
           })
        });

        $.each(this.sharedDatasets, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'shared_node',
                'text': item.text,
                'type': 'dataset',
                'a_attr': {
                    'class': "py-0",
                },
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });

        $.each(this.userDatasets, (_i, item) => {
            treeData.push({
                'id': item.value,
                'parent': 'user_node',
                'text': item.text,
                'type': 'dataset',
                'a_attr': {
                    'class': "py-0",
                },
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });
        this.treeData = treeData;
        return this.treeData;
    }

    generateTree() {

        this.generateTreeData();

        // Update existing tree or generate new tree if it doesn't exist
        if (this.tree) {
            this.updateTreeData()
        } else {
            // Instantiate the tree
            $(this.treeDiv).jstree({
                'core':{
                    'data':this.treeData,
                },
                'plugins': ["search", "types", "wholerow"],
                'search': {
                    "show_only_matches": true
                },
                'types': {
                    'default': {
                        'icon': 'fa fa-folder-o'
                    },
                    'dataset': {
                        'icon': 'fa fa-address-card-o',
                        'valid_children':[]
                    }
                }
            })
            this.setTree();
        }

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first()

        this.register_events();
    }

    // Register various DatasetTree events as object properties are updated.
    register_events() {
        const self = this;
        this.register_search();

        // Get layout from the selected node and close dropdown
        $(this.treeDiv).on('select_node.jstree', (_e, data) => {

            // Though you can select multiple nodes in the tree, let's only select the first
            const datasetId = data.selected[0];  // Returns node 'id' property
            if (data.node.type === "default") {
                // Do not toggle if user is navigating a branch node
                // NOTE: If tree is inside a <form>, which cannot be nested inside another <form>, this could toggle closed anyways due to the conflict.
                return;
            }
            const selectedNode = data.instance.get_node(datasetId);
            $(self.dropdownToggleElt).text(selectedNode.text);
            $(self.dropdownToggleElt).val(selectedNode.original.dataset_id);
            $(self.dropdownToggleElt).data("dataset-id", selectedNode.original.dataset_id);
            $(self.dropdownToggleElt).data("organism-id", selectedNode.original.organism_id);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(self.dropdownToggleElt).trigger('change');   // Force the change event to fire, triggering downstream things

        }).jstree(true);
    }

    loadFromDB() {
        //pass
    }

    saveToDB() {
        //pass
    }
}
