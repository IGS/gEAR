"use strict";

/*
Gene Cart and Profile tree stuff

Just about everything here uses the JSTree library - jstree.com
*/

/** Base class representing a tree structure */
class Tree {
    /**
     * Initialize tree.
     * @param {Object} - Tree data.
     */
    constructor({
        treeDiv
    } = {}) {
        this.treeDiv = treeDiv;

        this.setTree();
        this.setJSON();
    }

    setJSON() {
        this.treeJSON = (this.tree) ? this.tree.get_json() : undefined;
    }

    setTree() {
        this.tree = (this.treeDiv) ? $.jstree.reference(this.treeDiv) : undefined;
        // $(this.treeDiv).jstree(true) returns the same as $.jstree.reference(this.treeDiv)
        // so they may be interchanged in the codebase
    }

    // Update the contents of a JSTree object with new data.
    updateTreeData(newData=null) {
        (this.tree).settings.core.data = newData ? newData : this.treeData;
        (this.tree).refresh();
    }

}

/*
 * Class representing a gene cart selection tree
 * @extends Tree
 */
class GeneCartTree extends Tree {
    /**
     * Initialize GeneCartTree
     * @param {Object} Data - Tree data
     */
    constructor({
        ...args
    }={}, domainGeneCarts, userGeneCarts) {
        super(args);
        this.domainGeneCarts = (domainGeneCarts) ? domainGeneCarts : [];
        this.userGeneCarts = (userGeneCarts) ? userGeneCarts : [] ;

    }

    generateTreeData() {
        // Create JSON tree structure for the data
        let treeData = [
            {'id':'domain_node', 'parent':'#', 'text':"Public Gene Carts"},
            {'id':'user_node', 'parent':'#', 'text':"Your Gene Carts"},
        ];

        $.each(this.domainGeneCarts, function(i, item){
            treeData.push({
                'id': item.value,
                'parent': 'domain_node',  // All carts private for now
                'text':item.text,
                'type': 'gene',
                'a_attr': {
                    'class': "py-0 download-item",
                }
            })
        });

        $.each(this.userGeneCarts, function(i, item){
            treeData.push({
                'id': item.value,
                'parent': 'user_node',
                'text':item.text,
                'type': 'gene',
                'a_attr': {
                    'class': "py-0 download-item",
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
            return;
        } else {
            $(this.treeDiv).jstree({
                'core':{
                    'data':this.treeData,
                    // so that right-click->create works
                    //"check_callback" : true
                },
                'plugins': ["search", "types"],
                //'plugins': ["contextmenu", "dnd", "search", "types"],
                /* Plugins
                    contextmenu - Allows right-click of node for configurable actions
                    dnd - Allows drag-and-drop of nodes to rearrange tree
                    search - search for matching items and expand tree if found
                    types - Allows you to define node types with nesting rules and icons
                */
                'types': {
                    'default': {
                        'icon': 'fa fa-folder-o'
                    },
                    'gene': {
                        'icon': 'fa fa-random',
                        'valid_children':[]
                    }
                }
            })
            // The above JSTree call sets the tree, and this returns it
            this.tree = $(this.treeDiv).jstree(true);
        }

        // Code from "search" section of https://www.jstree.com/plugins/
        // Sets text input to search as tree search box.
        let to = false;
        $(`${this.treeDiv}_q`).keyup(function () {
            if (to) { clearTimeout(to); }
            to = setTimeout(function () {
            let v = $(`${this.treeDiv}_q`).val();
            $(this.treeDiv).jstree(true).search(v);
            }, 250);
        });

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first()
        this.register_events();

    }

    register_events() {
        let self = this;
        // Get genes from the selected gene cart
        $(this.treeDiv).on('select_node.jstree', function(e, data) {
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
     * @param {Object} Data - Tree data
     */
    constructor({
        ...args
    }={}, domainProfiles, userProfiles) {
        super(args);
        this.domainProfiles = (domainProfiles) ? domainProfiles : [];
        this.userProfiles = (userProfiles) ? userProfiles : [];
    }

    generateTreeData() {
        // Create JSON tree structure for the data
        let treeData = [
            {'id':'domain_node', 'parent':'#', 'text':"Public Profiles"},
            {'id':'user_node', 'parent':'#', 'text':"Your Profiles"},
        ];

        // user_profiles/domain_profiles properties - value, text, share_id

        // Load profiles into the tree data property
        $.each(this.domainProfiles, function(i, item){
            treeData.push({
                'id': item.value,
                'parent': 'domain_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': "domain_choice_c py-0 download-item",
                    'data-profile-label': item.text,
                    'data-profile-id': item.value,
                    'data-profile-share-id': item.share_id
                }
            })
        });

        $.each(this.userProfiles, function(i, item){
            treeData.push({
                'id': item.value,
                'parent': 'user_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': "domain_choice_c py-0 download-item",
                    'data-profile-label': item.text,
                    'data-profile-id': item.value,
                    'data-profile-share-id': item.share_id
                }
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
            return;
        } else {
            // Instantiate the tree
            $(this.treeDiv).jstree({
                'core':{
                    'data':this.treeData,
                    // so that right-click->create works
                    //"check_callback" : true
                },
                'plugins': ["search", "types"],
                // 'plugins': ["contextmenu", "dnd", "search", "types"],

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

            // The above JSTree call sets the tree, and this returns it
            this.tree = $(this.treeDiv).jstree(true);

        }

        // Sets text input to search as tree search box.
        // Code from "search" section of https://www.jstree.com/plugins/
        let to = false;
        $(`${this.treeDiv}_q`).keyup(function () {
            if (to) { clearTimeout(to); }
            to = setTimeout(function () {
            let v = $(`${this.treeDiv}_q`).val();
            $(this.treeDiv) .jstree(true).search(v);
            }, 250);
        });

        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first()

        this.register_events();
    }

    // Register various ProfileTree events as object properties are updated.
    register_events() {
        let self = this;

        // Get layout from the selected node and close dropdown
        $(this.treeDiv).on('select_node.jstree', function(e, data) {

            // Though you can select multiple nodes in the tree, let's only select the first
            const layoutId = data.selected[0];  // Returns node 'id' property
            if (data.node.type === "default") {
                // Do not toggle if user is navigating a branch node
                // NOTE: If tree is inside a <form>, which cannot be nested inside another <form>, this could toggle closed anyways due to the conflict.
                return;
            }
            // The dropdown toggle text/val change happens in DatasetCollectionPanel->set_layouts() for the index page
            // but this should be set to assist with other pages.
            const selectedNode = data.instance.get_node(layoutId);
            $(self.dropdownToggleElt).text(selectedNode.text);
            $(self.dropdownToggleElt).val(layoutId);
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown

        }).jstree(true);
    }

    loadFromDB() {
        //pass
    }

    saveToDB() {
        //pass
    }

}