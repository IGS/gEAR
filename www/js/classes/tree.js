"use strict";

/*
Gene Cart and Profile tree stuff
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
        this.tree = (this.treeDiv) ? this.treeDiv.jstree(true) : undefined;

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

        this.selectClass = 'js-genecart-select';

    }

    // Load all saved gene carts for the current user
    // TODO: Change based on gene cart manager page code
    generateGeneCartTree (treeDiv) {

        // Create JSON tree structure for the data
        let treeData = [
            {'id':'domain_node', 'parent':'#', 'text':"Public Gene Carts"},
            {'id':'user_node', 'parent':'#', 'text':"Your Gene Carts"},
        ];

        let self = this;
        $.each(this.domainGeneCarts, function(i, item){
            treeData.push({
                'id': item['id'],
                'parent': 'domain_node',  // Everything private for now
                'text':item['label'],
                'type': 'gene',
                'a_attr': {
                    'class': `py-0 download-item ${self.selectClass}`,
                }
            })
        });

        $.each(this.userGeneCarts, function(i, item){
            treeData.push({
                'id': item['id'],
                'parent': 'user_node',  // Everything private for now
                'text':item['label'],
                'type': 'gene',
                'a_attr': {
                    'class': `py-0 download-item ${self.selectClass}`,

                }
            })
        });

        this.tree = $(treeDiv).jstree({
            'core':{
                'data':treeData,
            },
            /* Plugins
                contextmenu - Allows right-click of node for configurable actions
                dnd - Allows drag-and-drop of nodes to rearrange tree
                search - search for matching items and expand tree if found
                types - Allows you to define node types with nesting rules and icons
            */
            'plugins': ["search", "types"],
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

        // Code from "search" section of https://www.jstree.com/plugins/
        // Sets text input to search as tree search box.
        let to = false;
        $(`#${treeDiv}_q`).keyup(function () {
            if (to) { clearTimeout(to); }
            to = setTimeout(function () {
            let v = $(`#${treeDiv}_q`).val();
            $(treeDiv).jstree(true).search(v);
            }, 250);
        });

        this.treeDiv =treeDiv;
        // NOTE: Using DOM tree traversal to get to the dropdown-toggle feels hacky
        this.dropdownElt = $(this.treeDiv).closest('.dropdown');
        // Get "toggle" for the dropdown tree. Should only be a single element, but "first()" is there for sanity's sake
        this.dropdownToggleElt = $(this.dropdownElt).children('.dropdown-toggle').first()

        this.register_events();

    }

    register_events() {
        let self = this;
        // Get genes from the selected gene cart
        // TODO: Change tree HTML id based on gene cart manager page code
        $(this.treeDiv).on('select_node.jstree', function(e, data) {
            // Though you can select multiple nodes in the tree, let's only select the first
            const geneCartId = data.selected[0];  // Returns node 'id' property
            if (["user_node", "domain_node"].includes(geneCartId)) {
                // Do not toggle if user is navigating a branch node
                return;
            }
            //const selectedNode = data.instance.get_node(geneCartId);
            //$(self.dropdownToggleElt).text(selectedNode.text);  // This already happens via the Bootstrap dropdown code
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown

            const params = { session_id: session_id, gene_cart_id: geneCartId };
            const d = new $.Deferred(); // Causes editable to wait until results are returned

            if (typeof session_id !== 'undefined') {
                // Get the gene cart members and populate the gene symbol search bar
                $.ajax({
                url: './cgi/get_gene_cart_members.cgi',
                type: 'post',
                data: params,
                success: function (data, newValue, oldValue) {
                    if (data.success === 1) {
                        const geneCartSymbols = []
                        // format gene symbols into search string
                        $.each(data.gene_symbols, function (i, item) {
                            geneCartSymbols.push(item.label);
                        });
                        //deduplicate gene cart
                        const geneCartSymbolsSet = [...new Set(geneCartSymbols)]

                        createGeneList(geneCartSymbolsSet);
                    } else {
                        alert("no genes in this cart");
                    }

                    d.resolve();
                }
                });
            } else {
                d.resolve();
            }
            return d.promise();
        }).jstree(true);
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

        this.selectClass = 'js-profile-select';
    }

    generateProfileTree(treeDiv) {

        // Create JSON tree structure for the data
        let treeData = [
            {'id':'domain_node', 'parent':'#', 'text':"Public Profiles"},
            {'id':'user_node', 'parent':'#', 'text':"Your Profiles"},
        ];

        // user_profiles/domain_profiles properties - value, text, share_id

        let self = this;
        // Load profiles into the tree data property
        $.each(this.domainProfiles, function(i, item){
            treeData.push({
                'id': item.value,
                'parent': 'domain_node',
                'text': item.text,
                'type': 'profile',
                'a_attr': {
                    'class': `domain_choice_c py-0 download-item ${self.selectClass}`,
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
                    'class': `domain_choice_c py-0 download-item ${self.selectClass}`,
                    'data-profile-label': item.text,
                    'data-profile-id': item.value,
                    'data-profile-share-id': item.share_id
                }
            })
        });

        // Instantiate the tree
        this.tree = $(treeDiv).jstree({
            'core':{
                'data':treeData,
            },
            'plugins': ["search", "types"],
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

        // Sets text input to search as tree search box.
        // Code from "search" section of https://www.jstree.com/plugins/
        let to = false;
        $(`${treeDiv}_q`).keyup(function () {
            if (to) { clearTimeout(to); }
            to = setTimeout(function () {
            let v = $(`${treeDiv}_q`).val();
            $(treeDiv) .jstree(true).search(v);
            }, 250);
        });

        this.treeDiv =treeDiv;
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
            if (["user_node", "domain_node"].includes(layoutId)) {
                // Do not toggle if user is navigating a branch node
                return;
            }
            //const selectedNode = data.instance.get_node(layoutId);
            //$(self.dropdownToggleElt).text(selectedNode.text);  // This already happens via the Bootstrap dropdown code
            $(self.dropdownToggleElt).dropdown('toggle');  // Close dropdown

        }).jstree(true);
    }

}


// TODO: Change HTML ids based on gene cart manager page code
function createGeneList(genes) {
    const tmpl = $.templates('#genes_tmpl');
    const data = { genes: genes };
    const html = tmpl.render(data);
    $('#genes_container').html(html);
}