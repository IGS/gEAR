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
        , storedValElt
    } = {}) {
        this.treeDiv = treeDiv; // Element to generate the tree structure on
        this.storedValElt = storedValElt;   // Element to store text, vals, and data properties on

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

    addNode(treeData, id, parentID, text, nodeType) {
        let nodeClass ='';

        if (nodeType == 'default') {
            nodeClass = 'jstree-ocl';
        } else if (nodeType == 'genecart') {
            nodeClass = 'py-0'
        }

        treeData.push({
            'id': id,
            'parent': parentID,
            'text': text,
            'type': nodeType,
            'a_attr': {
                'class': nodeClass,
            }
        })
    }

    generateTreeData() {
        // Create JSON tree structure for the data
        const treeKeys = {'domain_node': true, 'user_node': true, 'group_node': true, 'shared_node': true, 'public_node': true};
        const treeData = [
            {'id':'domain_node', 'parent':'#', 'text':`Highlighted gene carts (${this.domainGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'user_node', 'parent':'#', 'text':`Your gene carts (${this.userGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'group_node', 'parent':'#', 'text':`Group gene carts (${this.groupGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'shared_node', 'parent':'#', 'text':`Gene carts shared with you (${this.sharedGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'public_node', 'parent':'#', 'text':`Public carts from other users (${this.publicGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
        ];

        $.each(this.domainGeneCarts, (_i, item) => {
            // TODO: All this parent/grandparent logic should just go into addNode
            // If there's a parent make sure it's added, doesn't currently handle grandparents
            if (item.folder_parent_id && ! treeKeys.hasOwnProperty(item.folder_parent_id)) {
                this.addNode(treeData, item.folder_label, item.folder_parent_id, null, null, 'default');
                treeKeys[item.folder_parent_id] = true;
            }

            // Now do the same for the containing folder itself
            if (item.folder_id) {
                item.folder_id = 'folder-' + item.folder_id;

                if (item.folder_parent_id) {
                    //item.folder_parent_id = 'folder-' + item.folder_parent_id;
                } else {
                    item.folder_parent_id = 'domain_node';
                }

                if (! treeKeys.hasOwnProperty(item.folder_id)) {
                    this.addNode(treeData, item.folder_id, item.folder_parent_id, item.folder_label, 'default');
                    treeKeys[item.folder_id] = true;
                }

                this.addNode(treeData, item.value, item.folder_id, item.text, 'genecart');
            } else {
                // Profile isn't in any kind of folder, so just attach it to the top-level node of this type
                this.addNode(treeData, item.value, 'domain_node', item.text, 'genecart')
            }
        });

        $.each(this.userGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'user_node', item.text, 'genecart')
        });

        $.each(this.groupGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'group_node', item.text, 'genecart')
        });

        $.each(this.sharedGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'shared_node', item.text, 'genecart')
        });

        $.each(this.publicGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'public_node', item.text, 'genecart')
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

        // TODO: Add some configuration to the Tree objects to not rely as heavily on a particular HTML design
        /*
        Example of a tree container structure, which relies on the positioning of the dropdown class
        inspired by https://getbootstrap.com/docs/4.0/components/dropdowns/#examples

        For this example:
         '#datset_tree' is this.treeDiv
         '#dataset_c' is this.dropdownElt
         '#dataset' is this.dropdownToggleElt (and this.storedValElt)

        <div id="dataset_c" class="form-control dropdown" aria-describedby="dataset_help">
          <div id="dataset" value="" data-title="Click to change" class="text-truncate dropdown-toggle" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">Choose dataset</div>
          <div class="dropdown-menu w-100" aria-labelledby="dataset">
            <form class="px-4 py-3" onsubmit="return false;">
              <div class="input-group input-group-sm mb-3">
                <div class="input-group-prepend">
                  <i class="fa fa-search input-group-text" aria-hidden="true"></i>
                </div>
                <input type="text" class="form-control" id="dataset_tree_q" placeholder="Type to search datasets">
              </div>
              <div id="dataset_tree"></div>
            </form>
          </div>
        </div>  <!-- end dataset_c -->
        */

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first();
        // This element will store the text, value, and data properties of the selected node
        this.storedValElt = (this.storedValElt) ? this.storedValElt : this.dropdownToggleElt;
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
                // TODO: If a branch is selected, a max call stack is exceeded
                // Do not toggle if user is navigating a branch node
                // NOTE: If tree is inside a <form>, which cannot be nested inside another <form>, this could toggle closed anyways due to the conflict.
                return;
            }
            const selectedNode = data.instance.get_node(geneCartId);
            $(self.storedValElt).text(selectedNode.text);
            $(self.storedValElt).val(geneCartId);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(self.storedValElt).change(); // Force the change event to fire, triggering downstream things like getting cart members
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
    }={}, domainProfiles, userProfiles, groupProfiles, sharedProfiles, publicProfiles) {
        super(args);
        this.domainProfiles = (domainProfiles) ? domainProfiles : [];
        this.userProfiles = (userProfiles) ? userProfiles : [];
        this.groupProfiles = (groupProfiles) ? groupProfiles : [];
        this.sharedProfiles = (sharedProfiles) ? sharedProfiles : [];
        this.publicProfiles = (publicProfiles) ? publicProfiles : [];
        this.profileIDs = {};
        this.treeKeys = {};
        this.treeData = [];

        // This is needed so we can add folders with labels to the tree
        this.folders = [];
    }

    addNode(item, default_folder) {
        item.tree_id = item.value;
        
        // Modifies the ID to prevent duplicates. We don't use the ID attribute
        //  directly in the link anyway.
        while (this.profileIDs.hasOwnProperty(item.tree_id)) {
            // Add a dash to the exiting ID
            item.tree_id = item.tree_id + '-';
        }

        if (item.folder_id == null) {
            item.folder_id = default_folder;
        }

        let nodeData = {
                'id': item.tree_id,
                'parent': item.folder_id,
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': 'py-0',
                },
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
        };

        this.profileIDs[item.tree_id] = true;
        this.treeData.push(nodeData);
    }

    addFolder(folder) {
        // We skip folder ID 0, since it handled as id:domain_node already
        if (folder.id == 0) {
            return false;
        }
        
        if (folder.parent_id == null) {
            folder.parent_id = '#';
        }
        
        let nodeData = {
                'id': folder.id,
                'parent': folder.parent_id,
                'text': folder.label,
                'type': 'default',
                'a_attr': {
                    'class': 'jstree-ocl',
                },
                'profile_label': null,
                'profile_id': null,
                'profile_share_id': null
        };

        this.treeData.push(nodeData);
    }

    generateTreeData() {
        // Create JSON tree structure for the data
        this.treeKeys = {'domain_node': true, 'user_node': true, 'group_node': true, 'shared_node': true};
        this.treeData = [];

        // Add all the folders first
        $.each(this.folders, (_i, item) => {
            this.addFolder(item);
        });

        // Sort and add profiles into the tree data property
        this.domainProfiles.sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        for (const item of this.domainProfiles) {
            this.addNode(item, 'domain_node');
        };

        this.userProfiles.sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        for (const item of this.userProfiles) {
            this.addNode(item, 'user_node');
        };

        this.groupProfiles.sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        for (const item of this.groupProfiles) {
            this.addNode(item, 'group_node');
        };

        this.sharedProfiles.sort((a, b) => a.text.toLowerCase() > b.text.toLowerCase() ? 1 : -1);
        for (const item of this.sharedProfiles) {
            this.addNode(item, 'shared_node');
        };

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
                        'icon': 'fa fa-th-large',
                        'valid_children':[]
                    }
                }
            })
            this.setTree();
        }

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first();
        // This element will store the text, value, and data properties of the selected node
        this.storedValElt = (this.storedValElt) ? this.storedValElt : this.dropdownToggleElt;
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
            $(self.storedValElt).text(selectedNode.text);
            $(self.storedValElt).val(selectedNode.original.profile_id);
            $(self.storedValElt).data("profile-id", selectedNode.original.profile_id);
            $(self.storedValElt).data("profile-label", selectedNode.original.profile_label);
            $(self.storedValElt).data("profile-share-id", selectedNode.original.profile_share_id);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(self.storedValElt).trigger('change');   // Force the change event to fire, triggering downstream things

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
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first();
        // This element will store the text, value, and data properties of the selected node
        this.storedValElt = (this.storedValElt) ? this.storedValElt : this.dropdownToggleElt;
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
            $(self.storedValElt).text(selectedNode.text);
            $(self.storedValElt).val(selectedNode.original.dataset_id);
            $(self.storedValElt).data("dataset-id", selectedNode.original.dataset_id);
            $(self.storedValElt).data("organism-id", selectedNode.original.organism_id);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(self.storedValElt).trigger('change');   // Force the change event to fire, triggering downstream things

        }).jstree(true);
    }

    loadFromDB() {
        //pass
    }

    saveToDB() {
        //pass
    }
}
