"use strict";

/*
Gene Cart and Profile tree stuff

Just about everything here uses the JSTree library - jstree.com
*/

// TODO: Add options to initialize tree with a selected node.

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
        this.storedValElt = storedValElt || null;   // Element to store text, vals, and data properties on

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

    // Update the contents of a JSTree object with new data.
    updateTreeData(newData=null) {
        (this.tree).settings.core.data = newData || this.treeData;
        (this.tree).refresh();
    }

    // Add a nested node to the tree. Edits "treeData" inplace.
    // Will add a new parent node to the tree if it doesn't exist.
    addNestedNode(treeData, item, parentID, kwargs) {
        /*
        item.folder_id => ID of the folder
        item.folder_parent_id => ID of the parent folder of this folder.
        The item does not have nesting further than the grandparent level.
        */

        // If the item grandparent node does not exist, add the passed-in parent ID
        // Same goes for the if this item is in a folder or not.
        item.folder_parent_id = item.folder_parent_id || parentID;
        item.folder_id = item.folder_id ? `folder-${item.folder_id}` : parentID;

        // Add the grandparent node if it doesn't exist
        // parentID should already be in the tree.
        if (! this.treeKeys.hasOwnProperty(item.folder_parent_id)) {
            this.addNode(treeData, item.folder_parent_id, null, item.folder_label, 'default');
            treeKeys[item.folder_parent_id] = true;
        }

        // If this item is in a folder (branch), add folder to tree.
        if (! this.treeKeys.hasOwnProperty(item.folder_id)) {
            this.addNode(treeData, item.folder_id, item.folder_parent_id, item.folder_label, 'default');
            this.treeKeys[item.folder_id] = true;
        }

        // Add the current item to the tree. Will be a leaf node.
        this.addNode(treeData, item.value, item.folder_id , item.text, this.nodeType, ...kwargs);
    }

    // Add a node to the tree. Edits "treeData" inplace.
    addNode(treeData, id, parentID, text, nodeType, kwargs) {

        // Class is dependent on if node will be a leaf or branch ("default" types)
        const nodeClass = nodeType == 'default' ? 'jstree-ocl' : 'py-0';

        treeData.push({
            'id': id,
            'parent': parentID,
            'text': text,
            'type': nodeType,
            'a_attr': {
                'class': nodeClass,
            },
            ...kwargs
        })
    }

    // Generate the tree structure for the DOM and return the JSTree object.
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
                    // See https://stackoverflow.com/a/11508530 for how to use variables as object keys
                    [this.nodeType]: {
                        'icon': `fa ${this.leafIcon}`,
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
        this.registerEvents();
    }

    // Register various Tree events as object properties are updated.
    registerEvents(dataKeyAsVal="") {
        // datadKeyAsVal - if set, the data-key attribute of the selected node will be used as the value of the dropdown

        this.registerSearch();

        // Get genes from the selected gene cart
        $(this.treeDiv).on('select_node.jstree', (_e, data) => {
            // Though you can select multiple nodes in the tree, let's only select the first
            const selectedID = data.selected[0];  // Returns node 'id' property
            if (data.node.type === "default") {
                // TODO: If a branch is selected, a max call stack is exceeded
                // Do not toggle if user is navigating a branch node
                // NOTE: If tree is inside a <form>, which cannot be nested inside another <form>, this could toggle closed anyways due to the conflict.
                return;
            }
            const selectedNode = data.instance.get_node(selectedID);
            $(this.storedValElt).text(selectedNode.text);

            const val = dataKeyAsVal ? selectedNode.original[dataKeyAsVal] : selectedNode.id;
            $(this.storedValElt).val(val);
            // $(this.storedValElt).val(selectedNode.original.dataset_id);

            // If data attributes were passed into the node, store them in this element for easy retrieval
            for (const key in selectedNode.original) {
                $(this.storedValElt).data(key, selectedNode.original[key]);
            }
            $(this.dropdownToggleElt).dropdown('toggle');  // Close dropdown
            $(this.storedValElt).change(); // Force the change event to fire, triggering downstream things like getting cart members
        }).jstree(true);
    }

    registerSearch() {
        let to = false;
        // Requires searchbox to be named #{treeDiv}_q
        $(`${this.treeDiv}_q`).keyup(() => {
            if (to) { clearTimeout(to); }
            to = setTimeout(() => {
                const v = $(`${this.treeDiv}_q`).val();
                this.tree.search(v);
            }, 250);
        });
    }

}

class ProjectionSourceTree extends Tree {
    /**
     * Initialize projectionSourceTree
     * @constructor
     * @param {Object} Data - Tree data
     */
     constructor({
        ...args
    }={}, weightedDomainGeneCarts, weightedGroupGeneCarts, weightedUserGeneCarts, weightedSharedGeneCarts, weightedPublicGeneCarts,
        unweightedDomainGeneCarts, unweightedGroupGeneCarts, unweightedUserGeneCarts, unweightedSharedGeneCarts, unweightedPublicGeneCarts) {
        super(args);

        this.weighted = {
            domainGeneCarts: (weightedDomainGeneCarts) ? weightedDomainGeneCarts : []
            , groupGeneCarts: (weightedGroupGeneCarts) ? weightedGroupGeneCarts : []
            , userGeneCarts: (weightedUserGeneCarts) ? weightedUserGeneCarts : []
            , sharedGeneCarts: (weightedSharedGeneCarts) ? weightedSharedGeneCarts : []
            , publicGeneCarts: (weightedPublicGeneCarts) ? weightedPublicGeneCarts : []
        };
        this.unweighted = {
            domainGeneCarts: (unweightedDomainGeneCarts) ? unweightedDomainGeneCarts : []
            , groupGeneCarts: (unweightedGroupGeneCarts) ? unweightedGroupGeneCarts : []
            , userGeneCarts: (unweightedUserGeneCarts) ? unweightedUserGeneCarts : []
            , sharedGeneCarts: (unweightedSharedGeneCarts) ? unweightedSharedGeneCarts : []
            , publicGeneCarts: (unweightedPublicGeneCarts) ? unweightedPublicGeneCarts : []
        }
    }

    nodeType = "genecart";
    treeKeys = {'domain_node': true, 'user_node': true, 'group_node': true, 'shared_node': true, 'public_node': true};
    leafIcon = 'fa-shopping-cart';

    addGeneCartTreeData(geneCartTree) {
        // Merge unweighted gene cart properties into this object's properties
        this.unweighted.domainGeneCarts = geneCartTree.domainGeneCarts;
        this.unweighted.groupGeneCarts = geneCartTree.groupGeneCarts;
        this.unweighted.userGeneCarts = geneCartTree.userGeneCarts;
        this.unweighted.sharedGeneCarts = geneCartTree.sharedGeneCarts;
        this.unweighted.publicGeneCarts = geneCartTree.publicGeneCarts;
    }

    getTotalWeightedCarts() {
        // get the total number of weighted gene carts
        return Object.keys(this.weighted).reduce((acc, curr) => acc + this.weighted[curr].length, 0);
    }

    getTotalUnweightedCarts() {
        // get the total number of unweighted gene carts
        return Object.keys(this.unweighted).reduce((acc, curr) => acc + this.unweighted[curr].length, 0);
    }

    generateTreeData() {
        // Create JSON tree structure for the data
        const treeData = [
            {'id':'unweighted_genes_node', 'parent':'#', 'text':`Unweighted Genes(${this.getTotalUnweightedCarts()})`, 'type':'default', 'a_attr':{'class':'jstree-ocl'}},
            {'id':'uw_domain_node', 'parent':'unweighted_genes_node', 'text':`Highlighted gene carts (${this.unweighted.domainGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'uw_user_node', 'parent':'unweighted_genes_node', 'text':`Your gene carts (${this.unweighted.userGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'uw_group_node', 'parent':'unweighted_genes_node', 'text':`Group gene carts (${this.unweighted.groupGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'uw_shared_node', 'parent':'unweighted_genes_node', 'text':`Gene carts shared with you (${this.unweighted.sharedGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'uw_public_node', 'parent':'unweighted_genes_node', 'text':`Public carts from other users (${this.unweighted.publicGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'weighted_genes_node', 'parent':'#', 'text':`Weighted Genes (${this.getTotalWeightedCarts()})`, 'type':'default', 'a_attr':{'class':'jstree-ocl'}},
            {'id':'w_domain_node', 'parent':'weighted_genes_node', 'text':`Highlighted gene carts (${this.weighted.domainGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'w_user_node', 'parent':'weighted_genes_node', 'text':`Your gene carts (${this.weighted.userGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'w_group_node', 'parent':'weighted_genes_node', 'text':`Group gene carts (${this.weighted.groupGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'w_shared_node', 'parent':'weighted_genes_node', 'text':`Gene carts shared with you (${this.weighted.sharedGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
            {'id':'w_public_node', 'parent':'weighted_genes_node', 'text':`Public carts from other users (${this.weighted.publicGeneCarts.length})`, 'a_attr': {'class':'jstree-ocl'}},
        ];

        $.each(this.weighted.domainGeneCarts, (_i, item) => {
            this.addNestedNode(treeData, item, 'w_domain_node');
        });

        $.each(this.weighted.userGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'w_user_node', item.text, this.nodeType)
        });

        $.each(this.weighted.groupGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'w_group_node', item.text, this.nodeType)
        });

        $.each(this.weighted.sharedGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'w_shared_node', item.text, this.nodeType)
        });

        $.each(this.weighted.publicGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'w_public_node', item.text, this.nodeType)
        });

        $.each(this.unweighted.domainGeneCarts, (_i, item) => {
            this.addNestedNode(treeData, item, 'uw_domain_node');
        });

        $.each(this.unweighted.userGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'uw_user_node', item.text, this.nodeType)
        });

        $.each(this.unweighted.groupGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'uw_group_node', item.text, this.nodeType)
        });

        $.each(this.unweighted.sharedGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'uw_shared_node', item.text, this.nodeType)
        });

        $.each(this.unweighted.publicGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'uw_public_node', item.text, this.nodeType)
        });

        this.treeData = treeData;
        return this.treeData;
    }

    generateTree() { super.generateTree(); }

    registerEvents() { super.registerEvents(); }
};

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

    nodeType = 'genecart';
    treeKeys = {'domain_node': true, 'user_node': true, 'group_node': true, 'shared_node': true, 'public_node': true};
    leafIcon = 'fa-shopping-cart';

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
            addNestedNode(treeData, item, "domain_node")
        });

        $.each(this.userGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'user_node', item.text, this.nodeType)
        });

        $.each(this.groupGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'group_node', item.text, this.nodeType)
        });

        $.each(this.sharedGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'shared_node', item.text, this.nodeType)
        });

        $.each(this.publicGeneCarts, (_i, item) => {
            this.addNode(treeData, item.value, 'public_node', item.text, this.nodeType)
        });

        this.treeData = treeData;
        return this.treeData;
    }

    generateTree() { super.generateTree(); }

    registerEvents() { super.registerEvents(); }
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

    nodeType = 'profile';
    treeKeys = {'domain_node': true, 'user_node': true, 'group_node': true, 'shared_node': true};
    leafIcon = "fa-th-large";

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
            this.addNode(treeData, item.value, "domain_node", item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        $.each(this.userProfiles, (_i, item) => {
            this.addNode(treeData, item.value, "user_node", item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        $.each(this.groupProfiles, (_i, item) => {
            this.addNestedNode(treeData, item, "group_node", {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            });

        });

        $.each(this.sharedProfiles, (_i, item) => {
            this.addNode(treeData, item.value, "shared_node", item.text, this.nodeType, {
                'profile_label': item.text,
                'profile_id': item.value,
                'profile_share_id': item.share_id
            })
        });

        this.treeData = treeData;
        return this.treeData;
    }

    generateTree() { super.generateTree(); }

    registerEvents() { super.registerEvents(); }
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

    nodeType = 'dataset';
    treeKeys = {"domain_node": true, "shared_node": true, "user_node": true};
    leafIcon = "fa-address-card-o";

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
            this.addNode(treeData, item.value, 'domain_node', item.text, this.nodeType, {
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });

        $.each(this.sharedDatasets, (_i, item) => {
            this.addNode(treeData, item.value, 'shared_node', item.text, this.nodeType, {
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });

        $.each(this.userDatasets, (_i, item) => {
            this.addNode(treeData, item.value, 'user_node', item.text, this.nodeType, {
                "dataset_id": item.dataset_id,
                'organism_id': item.organism_id
            })
        });
        this.treeData = treeData;
        return this.treeData;
    }

    generateTree() { super.generateTree(); }

    registerEvents() { super.registerEvents("dataset_id"); }
}
